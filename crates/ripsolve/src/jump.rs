//! Feasibility jump as Luteberget and Sartor describe it, and as HiGHS runs it.
//!
//! Every column holds one candidate move: the value that minimises its weighted
//! violation with everything else held where it is. The score of that move is kept up
//! to date as other columns move, by the change in each row it shares with them, which
//! is a constant amount of work for each entry of those rows; the value itself is
//! recomputed only when the column moves. Moves are chosen by sampling a few of the
//! columns whose move improves, and at a local minimum the weights of the violated rows
//! go up and a column from one of them moves regardless.
//!
//! `heuristic::feasibility_jump_within` is the older search, which recomputes the
//! value of every column in every row a move touches and always takes the best move.
//! That costs far more for each move and never moves uphill, and on `brazil3`,
//! `CMS750_4`, `rococoC10-001000` and `academictimetablesmall` it finds nothing in 30
//! seconds where HiGHS's search finds a point in under four.

use crate::heuristic::{Incumbent, Limits, is_feasible};
use crate::model::Problem;
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::Instant;

/// Columns sampled when choosing among the improving moves.
const SAMPLE: usize = 25;
/// The chance of taking an improving column at random instead.
const RANDOM_COLUMN: f64 = 0.001;
/// The chance, at a local minimum, of moving a random column of the chosen violated
/// row rather than its best one.
const RANDOM_ENTRY: f64 = 0.01;
/// Work without fewer violated rows or a better point before the search stops, for
/// each nonzero of the model. HiGHS's figure.
const STALL_PER_NONZERO: usize = 256;
/// A better point must improve the objective by this much, relatively, to count.
const MIN_IMPROVEMENT: f64 = 1e-4;
/// Work between readings of the clock and the cancel flag.
const CLOCK_STRIDE: usize = 1 << 18;

/// One side of a row: `lhs >= rhs` when `at_least`, `lhs <= rhs` otherwise.
struct Side {
    row: usize,
    rhs: f64,
    at_least: bool,
    weight: f64,
    lhs: f64,
    /// Position in the violated list, if violated.
    slot: Option<usize>,
}

impl Side {
    /// Zero when satisfied, minus the violation otherwise.
    fn score(&self, lhs: f64) -> f64 {
        if self.at_least {
            -(self.rhs - lhs).max(0.0)
        } else {
            -(lhs - self.rhs).max(0.0)
        }
    }
}

struct Search<'a> {
    lb: Vec<f64>,
    ub: Vec<f64>,
    integer: Vec<bool>,
    obj: Vec<f64>,
    sides: Vec<Side>,
    /// Each column's sides and coefficients.
    cells: Vec<Vec<(u32, f64)>>,
    /// Each side's columns and coefficients, from the row.
    csr: &'a crate::sparse::SparseMatrix,
    x: Vec<f64>,
    objective: f64,
    violated: Vec<u32>,
    move_value: Vec<f64>,
    move_score: Vec<f64>,
    good: Vec<u32>,
    good_slot: Vec<i32>,
    objective_weight: f64,
    increment: f64,
    effort: usize,
    tolerance: f64,
    rng: u64,
    buffer: Vec<(f64, f64)>,
}

impl<'a> Search<'a> {
    fn new(
        problem: &Problem,
        csr: &'a crate::sparse::SparseMatrix,
        start: &[f64],
        tolerance: f64,
        seed: u64,
    ) -> Self {
        let n = problem.n_cols();
        let m = problem.n_rows();
        let integer: Vec<bool> = (0..n).map(|j| problem.is_integer(j)).collect();
        let (lb, ub): (Vec<f64>, Vec<f64>) = (0..n)
            .map(|j| {
                let (lo, hi) = (problem.col_lb[j], problem.col_ub[j]);
                if integer[j] {
                    ((lo - tolerance).ceil(), (hi + tolerance).floor())
                } else {
                    (lo, hi)
                }
            })
            .unzip();
        let mut sides = Vec::new();
        let mut sides_of_row = vec![[u32::MAX; 2]; m];
        for (i, row_sides) in sides_of_row.iter_mut().enumerate() {
            let mut k = 0;
            for (rhs, at_least) in [(problem.row_lb[i], true), (problem.row_ub[i], false)] {
                if rhs.is_finite() && !csr.column(i).0.is_empty() {
                    row_sides[k] = sides.len() as u32;
                    k += 1;
                    sides.push(Side {
                        row: i,
                        rhs,
                        at_least,
                        weight: 1.0,
                        lhs: 0.0,
                        slot: None,
                    });
                }
            }
        }
        let mut cells: Vec<Vec<(u32, f64)>> = vec![Vec::new(); n];
        for (j, column) in cells.iter_mut().enumerate() {
            let (rows, vals) = problem.matrix.column(j);
            for (&i, &a) in rows.iter().zip(vals) {
                for &s in &sides_of_row[i] {
                    if s != u32::MAX {
                        column.push((s, a));
                    }
                }
            }
        }
        let x: Vec<f64> = (0..n)
            .map(|j| {
                let v = if integer[j] { start[j].round() } else { start[j] };
                v.clamp(lb[j], ub[j])
            })
            .collect();
        let mut search = Search {
            lb,
            ub,
            integer,
            obj: problem.obj.clone(),
            sides,
            cells,
            csr,
            x,
            objective: 0.0,
            violated: Vec::new(),
            move_value: vec![0.0; n],
            move_score: vec![0.0; n],
            good: Vec::new(),
            good_slot: vec![-1; n],
            objective_weight: 0.0,
            increment: 1.0,
            effort: 0,
            tolerance,
            rng: seed | 1,
            buffer: Vec::new(),
        };
        search.reset_incumbent();
        search.effort += problem.matrix.nnz();
        for j in 0..n {
            search.reset_move(j);
        }
        search
    }

    fn next(&mut self) -> u64 {
        self.rng ^= self.rng << 13;
        self.rng ^= self.rng >> 7;
        self.rng ^= self.rng << 17;
        self.rng
    }

    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }

    fn chance(&mut self, p: f64) -> bool {
        ((self.next() >> 11) as f64 / (1u64 << 53) as f64) < p
    }

    fn reset_incumbent(&mut self) {
        let csr = self.csr;
        self.objective = self.obj.iter().zip(&self.x).map(|(c, v)| c * v).sum();
        self.violated.clear();
        for s in 0..self.sides.len() {
            let row = self.sides[s].row;
            let (cols, vals) = csr.column(row);
            let lhs: f64 = cols.iter().zip(vals).map(|(&j, &a)| a * self.x[j]).sum();
            self.sides[s].lhs = lhs;
            if self.sides[s].score(lhs) < -self.tolerance {
                self.sides[s].slot = Some(self.violated.len());
                self.violated.push(s as u32);
            } else {
                self.sides[s].slot = None;
            }
        }
    }

    /// The value column `j` would jump to: the minimiser of its weighted violation,
    /// other columns held, and never its current value.
    fn jump_value(&mut self, j: usize) -> f64 {
        self.buffer.clear();
        let current = self.x[j];
        let (lower, upper) = (self.lb[j], self.ub[j]);
        let mut value = lower;
        let mut slope = 0.0;
        for &(s, a) in &self.cells[j] {
            let side = &self.sides[s as usize];
            let residual = side.lhs - a * current;
            let bound = if side.at_least {
                (side.rhs, f64::INFINITY)
            } else {
                (f64::NEG_INFINITY, side.rhs)
            };
            // The range of the column's values that satisfies this side.
            let (p, q) = ((bound.0 - residual) / a, (bound.1 - residual) / a);
            let (mut lo, mut hi) = if a > 0.0 { (p, q) } else { (q, p) };
            if self.integer[j] {
                lo = (lo - self.tolerance).ceil();
                hi = (hi + self.tolerance).floor();
            }
            if lo > hi {
                continue;
            }
            if lo > value {
                slope -= side.weight;
                if lo < upper {
                    self.buffer.push((lo, side.weight));
                }
            }
            if hi <= value {
                slope += side.weight;
            } else if hi < upper {
                self.buffer.push((hi, side.weight));
            }
        }
        if lower.is_finite() {
            self.buffer.push((lower, 0.0));
        }
        if upper.is_finite() {
            self.buffer.push((upper, 0.0));
        }
        self.buffer.sort_by(|a, b| a.0.total_cmp(&b.0));
        let mut score = 0.0;
        value = self.buffer.first().map_or(current, |b| b.0);
        let (mut best_score, mut best_value) = (score, value);
        let same = |a: f64, b: f64| (a - b).abs() < self.tolerance;
        for &(at, step) in &self.buffer {
            score += (at - value) * slope;
            slope += step;
            value = at;
            if same(best_value, current) || (!same(value, current) && score < best_score) {
                best_score = score;
                best_value = value;
            }
            if !same(best_value, current) && slope >= 0.0 {
                break;
            }
        }
        best_value
    }

    fn reset_move(&mut self, j: usize) {
        self.effort += self.cells[j].len();
        let value = self.jump_value(j);
        let delta = value - self.x[j];
        let mut score = self.objective_weight * self.obj[j] * delta;
        for &(s, a) in &self.cells[j] {
            let side = &self.sides[s as usize];
            score += side.weight * (side.score(side.lhs + a * delta) - side.score(side.lhs));
        }
        self.move_value[j] = value;
        self.move_score[j] = score;
        self.update_good(j);
    }

    fn update_good(&mut self, j: usize) {
        let improving = self.move_score[j] > 0.0;
        let slot = self.good_slot[j];
        if improving && slot < 0 {
            self.good_slot[j] = self.good.len() as i32;
            self.good.push(j as u32);
        } else if !improving && slot >= 0 {
            let last = *self.good.last().expect("a slot means the set is not empty");
            self.good.swap_remove(slot as usize);
            if last as usize != j {
                self.good_slot[last as usize] = slot;
            }
            self.good_slot[j] = -1;
        }
    }

    fn select(&mut self) -> Option<usize> {
        let csr = self.csr;
        if !self.good.is_empty() {
            if self.chance(RANDOM_COLUMN) {
                let k = self.below(self.good.len());
                return Some(self.good[k] as usize);
            }
            let samples = SAMPLE.min(self.good.len());
            self.effort += samples;
            let mut best: Option<(usize, f64)> = None;
            for _ in 0..samples {
                let k = self.below(self.good.len());
                let j = self.good[k] as usize;
                if best.is_none_or(|(_, held)| self.move_score[j] > held) {
                    best = Some((j, self.move_score[j]));
                }
            }
            return best.map(|(j, _)| j);
        }
        self.bump_weights();
        if self.violated.is_empty() {
            let k = self.below(self.x.len());
            return Some(k);
        }
        let k = self.below(self.violated.len());
        let row = self.sides[self.violated[k] as usize].row;
        let (cols, _) = csr.column(row);
        if self.chance(RANDOM_ENTRY) {
            let k = self.below(cols.len());
            return Some(cols[k]);
        }
        let mut best: Option<(usize, f64)> = None;
        for &j in cols {
            if best.is_none_or(|(_, held)| self.move_score[j] > held) {
                best = Some((j, self.move_score[j]));
            }
        }
        best.map(|(j, _)| j)
    }

    /// A local minimum: make the violated rows dearer, or, when nothing is violated,
    /// the objective.
    fn bump_weights(&mut self) {
        let csr = self.csr;
        let increment = self.increment;
        if self.violated.is_empty() {
            self.objective_weight += increment;
            self.effort += self.x.len();
            for j in 0..self.x.len() {
                self.move_score[j] += increment * self.obj[j] * (self.move_value[j] - self.x[j]);
                self.update_good(j);
            }
        } else {
            for k in 0..self.violated.len() {
                let s = self.violated[k] as usize;
                self.sides[s].weight += increment;
                let (row, lhs) = (self.sides[s].row, self.sides[s].lhs);
                let (cols, vals) = csr.column(row);
                self.effort += cols.len();
                for (&j, &a) in cols.iter().zip(vals) {
                    let side = &self.sides[s];
                    let after = lhs + a * (self.move_value[j] - self.x[j]);
                    self.move_score[j] += increment * (side.score(after) - side.score(lhs));
                    self.update_good(j);
                }
            }
        }
        // Weights only grow; rescale before they lose the precision to compare.
        let big = self.objective_weight > 1e20 || self.sides.iter().any(|s| s.weight > 1e20);
        if big {
            self.increment *= 1e-20;
            self.objective_weight *= 1e-20;
            for side in &mut self.sides {
                side.weight *= 1e-20;
            }
            for j in 0..self.x.len() {
                self.reset_move(j);
            }
        }
    }

    fn apply(&mut self, j: usize) {
        let csr = self.csr;
        let value = self.move_value[j];
        let delta = value - self.x[j];
        self.x[j] = value;
        self.objective += self.obj[j] * delta;
        for c in 0..self.cells[j].len() {
            let (s, a) = self.cells[j][c];
            let s = s as usize;
            let old = self.sides[s].lhs;
            let new = old + a * delta;
            self.sides[s].lhs = new;
            let violated = self.sides[s].score(new) < -self.tolerance;
            match (violated, self.sides[s].slot) {
                (true, None) => {
                    self.sides[s].slot = Some(self.violated.len());
                    self.violated.push(s as u32);
                }
                (false, Some(slot)) => {
                    let last = *self.violated.last().expect("a slot means the list is not empty");
                    self.violated.swap_remove(slot);
                    if last as usize != s {
                        self.sides[last as usize].slot = Some(slot);
                    }
                    self.sides[s].slot = None;
                }
                _ => {}
            }
            let row = self.sides[s].row;
            let (cols, vals) = csr.column(row);
            self.effort += cols.len();
            for (&k, &b) in cols.iter().zip(vals) {
                if k == j {
                    continue;
                }
                let side = &self.sides[s];
                let shift = b * (self.move_value[k] - self.x[k]);
                let before = side.weight * (side.score(old + shift) - side.score(old));
                let after = side.weight * (side.score(new + shift) - side.score(new));
                self.move_score[k] += after - before;
                self.update_good(k);
            }
        }
        self.reset_move(j);
    }
}

/// Search for a feasible point from `start`, spending at most `budget` work, and
/// return the best one found.
///
/// It keeps going after the first point, as HiGHS's does, with the objective weighted
/// in at each local minimum, and stops when a stretch of work has neither reduced the
/// violated rows below their least nor found a better point.
pub fn jump(
    problem: &Problem,
    start: &[f64],
    limits: &Limits,
    seed: u64,
    budget: &mut usize,
    deadline: Option<Instant>,
    cancel: Option<&AtomicBool>,
) -> Option<Incumbent> {
    let n = problem.n_cols();
    if n == 0 || start.len() != n {
        return None;
    }
    let tolerance = limits.feasibility_tolerance;
    let csr = problem.matrix.to_csr();
    let mut search = Search::new(problem, &csr, start, tolerance, seed);
    if (0..n).any(|j| search.lb[j] > search.ub[j]) {
        return None;
    }
    let stall = STALL_PER_NONZERO.saturating_mul(problem.matrix.nnz().max(1));
    let mut best_violated = usize::MAX;
    let mut best: Option<(f64, Vec<f64>)> = None;
    let mut improved_at = 0usize;
    let mut next_reading = 0usize;
    loop {
        if search.violated.len() < best_violated {
            best_violated = search.violated.len();
            improved_at = search.effort;
        }
        if search.violated.is_empty() {
            let better = best.as_ref().is_none_or(|(held, _)| {
                (held - search.objective) / search.objective.abs().max(1.0) > MIN_IMPROVEMENT
            });
            if better && is_feasible(problem, &search.x, tolerance) {
                improved_at = search.effort;
                best = Some((search.objective, search.x.clone()));
            }
        }
        if search.effort >= *budget || search.effort - improved_at > stall {
            break;
        }
        if search.effort >= next_reading {
            next_reading = search.effort + CLOCK_STRIDE;
            if deadline.is_some_and(|d| Instant::now() >= d)
                || cancel.is_some_and(|c| c.load(Ordering::Relaxed))
            {
                break;
            }
        }
        let Some(j) = search.select() else {
            break;
        };
        search.apply(j);
    }
    *budget = budget.saturating_sub(search.effort);
    best.map(|(_, x)| Incumbent {
        objective: problem.obj.iter().zip(&x).map(|(c, v)| c * v).sum(),
        x,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::generate::{Kind, Spec};
    use lp_parser_rs::problem::LpProblem;

    fn parse(lp: &str) -> Problem {
        Problem::from_lp(&LpProblem::parse(lp).unwrap()).unwrap()
    }

    fn run(p: &Problem, start: &[f64]) -> Option<Incumbent> {
        let mut budget = 10_000_000;
        jump(p, start, &Limits::default(), 7, &mut budget, None, None)
    }

    /// Whatever it returns is a point of the model: checked before it is handed back.
    #[test]
    fn never_returns_an_infeasible_point() {
        let mut found = 0;
        for seed in 0..40u64 {
            for kind in [Kind::Covering, Kind::Signed, Kind::Knapsack] {
                let spec = Spec {
                    kind,
                    n_cols: 30,
                    n_rows: 25,
                    seed,
                };
                let p = parse(&spec.to_lp());
                let start = vec![0.0; p.n_cols()];
                if let Some(point) = run(&p, &start) {
                    found += 1;
                    assert!(is_feasible(&p, &point.x, 1e-6), "{kind:?} seed {seed}");
                }
            }
        }
        assert!(found > 0, "found nothing on any generated model");
    }

    /// A local minimum of the weighted violation is left by moving uphill: from all
    /// zeros, one column on is a local minimum with a row still violated.
    #[test]
    fn leaves_a_local_minimum() {
        let p = parse(
            "Minimize\n obj: x1 + x2 + x3\nSubject To\n r1: x1 + x2 >= 1\n \
             r2: x2 + x3 >= 1\n r3: x1 + x3 >= 1\nBinary\n x1\n x2\n x3\nEnd\n",
        );
        let point = run(&p, &[0.0; 3]).expect("two of the three on satisfies every row");
        assert!(is_feasible(&p, &point.x, 1e-6));
    }

    /// General integers, continuous columns, a negative coefficient and an equality.
    #[test]
    fn moves_general_integers_and_continuous_columns() {
        let p = parse(
            "Minimize\n obj: x + y + z\nSubject To\n r1: 3 x + 2 y >= 7\n \
             r2: y - z <= 1\n r3: x + z >= 4.5\n r4: x - y + w = 0\n\
             Bounds\n 0 <= x <= 10\n 0 <= y <= 10\n 0 <= z <= 5\n -10 <= w <= 10\n\
             General\n x\n y\n w\nEnd\n",
        );
        let point = run(&p, &vec![0.0; p.n_cols()]).expect("x = 2, y = 1, z = 2.5 is one");
        assert!(is_feasible(&p, &point.x, 1e-6));
        for j in p.integer_columns() {
            assert_eq!(point.x[j], point.x[j].round(), "integer column {j} left fractional");
        }
    }

    /// After the first point the objective is weighed in, and what comes back is the
    /// best point seen: here every column on is feasible, and the optimum has one off.
    #[test]
    fn improves_on_its_first_point() {
        let p = parse(
            "Minimize\n obj: x1 + x2 + x3\nSubject To\n r1: x1 + x2 + x3 >= 2\n\
             Binary\n x1\n x2\n x3\nEnd\n",
        );
        let point = run(&p, &[1.0; 3]).expect("the start is already feasible");
        assert_eq!(point.objective, 2.0);
    }
}
