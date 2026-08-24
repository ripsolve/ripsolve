//! A bounded-variable revised primal simplex, with a phase-1 for feasibility.
//!
//! # Computational form
//!
//! The model's range rows are turned into equalities by giving every row a logical
//! (slack) variable:
//!
//! ```text
//!   min  c'z    subject to    [A | -I] z = 0,    l <= z <= u
//! ```
//!
//! where `z = [x; s]` puts the `n` structural variables first and the `m` logicals
//! after. A logical carries its row's bounds, so `<=`, `>=`, `=` and range rows all
//! reduce to a bound pair and the simplex never branches on row sense. Structural
//! bounds start at `[0, 1]`; branching later narrows them to `[0,0]` or `[1,1]`.
//!
//! The starting basis is all-logical, whose matrix is `-I` — trivially invertible,
//! and generally primal infeasible, which is what phase 1 is for.
//!
//! # Phase 1
//!
//! Phase 1 minimizes the sum of bound violations of the basic variables. Its
//! gradient is `-1` for a basic below its lower bound and `+1` for one above its
//! upper bound, and zero elsewhere; nonbasic variables sit on a bound and so never
//! contribute. The ratio test stops at the first breakpoint — the first basic
//! variable to reach a bound, whether it is becoming feasible or losing
//! feasibility. That is the textbook short-step rule: monotone and easy to verify,
//! but it takes more iterations than the long-step piecewise test a mature solver
//! uses. Replacing it is an M1b concern, not a correctness one.

use crate::lp::basis::{Basis, BasisError};
use crate::model::Problem;
use crate::sparse::SparseMatrix;

/// Numerical tolerances. Defaults follow the usual simplex conventions.
#[derive(Clone, Copy, Debug)]
pub struct Tolerances {
    /// A bound violation at or below this is not a violation.
    pub primal_feasibility: f64,
    /// A reduced cost between `-dual_feasibility` and `+dual_feasibility` is zero,
    /// so a column that far from profitable is not worth entering.
    pub dual_feasibility: f64,
    /// Pivots smaller than this are rejected as numerically unsafe.
    pub pivot: f64,
}

impl Default for Tolerances {
    fn default() -> Self {
        Self {
            primal_feasibility: 1e-7,
            dual_feasibility: 1e-7,
            pivot: 1e-9,
        }
    }
}

/// How a simplex solve ended.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LpStatus {
    Optimal,
    /// No point satisfies the constraints and bounds.
    Infeasible,
    /// The objective is unbounded below on the feasible region.
    Unbounded,
    /// The iteration limit was reached before the search concluded.
    IterationLimit,
}

/// The result of solving an LP relaxation.
#[derive(Clone, Debug)]
pub struct LpSolution {
    pub status: LpStatus,
    /// Objective of the internal minimization form, valid when `status` is
    /// [`LpStatus::Optimal`].
    pub objective: f64,
    /// Values of the structural variables only.
    pub x: Vec<f64>,
    pub iterations: usize,
}

/// Where a nonbasic variable is parked.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum At {
    Lower,
    Upper,
    /// A free variable, held at zero.
    Zero,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Status {
    Basic { row: usize },
    NonBasic(At),
}

/// An LP in computational form, ready to solve.
pub struct Lp {
    n_structural: usize,
    m: usize,
    cost: Vec<f64>,
    lower: Vec<f64>,
    upper: Vec<f64>,
    /// The original `A`, column-major, `m x n_structural`.
    matrix: SparseMatrix,
    tol: Tolerances,
}

impl Lp {
    /// Build the LP relaxation of a binary program.
    ///
    /// Column bounds are taken from the problem, so a problem whose bounds have
    /// already been narrowed by presolve or branching relaxes accordingly.
    pub fn relaxation(problem: &Problem) -> Lp {
        let n = problem.n_cols();
        let m = problem.n_rows();
        let mut cost = problem.obj.clone();
        cost.resize(n + m, 0.0);

        let mut lower = problem.col_lb.clone();
        let mut upper = problem.col_ub.clone();
        lower.extend_from_slice(&problem.row_lb);
        upper.extend_from_slice(&problem.row_ub);

        Lp {
            n_structural: n,
            m,
            cost,
            lower,
            upper,
            matrix: problem.matrix.clone(),
            tol: Tolerances::default(),
        }
    }

    pub fn with_tolerances(mut self, tol: Tolerances) -> Self {
        self.tol = tol;
        self
    }

    fn n_total(&self) -> usize {
        self.n_structural + self.m
    }

    /// Scatter column `j` of `[A | -I]` into a dense buffer of length `m`.
    fn column_into(&self, j: usize, out: &mut [f64]) {
        out.fill(0.0);
        if j < self.n_structural {
            let (rows, vals) = self.matrix.column(j);
            for (&i, &v) in rows.iter().zip(vals) {
                out[i] = v;
            }
        } else {
            out[j - self.n_structural] = -1.0;
        }
    }

    /// The bound a nonbasic variable sits on.
    fn value_at(&self, j: usize, at: At) -> f64 {
        match at {
            At::Lower => self.lower[j],
            At::Upper => self.upper[j],
            At::Zero => 0.0,
        }
    }

    /// Solve, returning the optimum of the internal minimization form.
    pub fn solve(&self) -> LpSolution {
        self.solve_with_limit(200_000)
    }

    pub fn solve_with_limit(&self, max_iterations: usize) -> LpSolution {
        Solver::new(self).run(max_iterations)
    }
}

/// Mutable state for one solve.
struct Solver<'a> {
    lp: &'a Lp,
    basis: Basis,
    /// Which variable is basic in each row.
    basic: Vec<usize>,
    status: Vec<Status>,
    /// Current value of every variable.
    z: Vec<f64>,
    /// Scratch buffers, reused across iterations to keep the loop allocation-free.
    col: Vec<f64>,
    alpha: Vec<f64>,
    y: Vec<f64>,
    cost_b: Vec<f64>,
    rhs: Vec<f64>,
}

/// What the ratio test decided.
enum Step {
    /// Pivot: `leaving_row` leaves, entering advances by `step`.
    Pivot {
        leaving_row: usize,
        step: f64,
        to: At,
    },
    /// The entering variable reached its opposite bound first; no basis change.
    BoundFlip { step: f64 },
    /// The objective improves without limit along this ray.
    Unbounded,
}

impl<'a> Solver<'a> {
    fn new(lp: &'a Lp) -> Solver<'a> {
        let n_total = lp.n_total();
        let m = lp.m;

        // Start with every logical basic and every structural nonbasic on the bound
        // nearer zero, which for a binary relaxation is the lower bound.
        let mut status = vec![Status::NonBasic(At::Lower); n_total];
        let mut z = vec![0.0; n_total];
        for j in 0..lp.n_structural {
            let at = if lp.lower[j].is_finite() {
                At::Lower
            } else if lp.upper[j].is_finite() {
                At::Upper
            } else {
                At::Zero
            };
            status[j] = Status::NonBasic(at);
            z[j] = lp.value_at(j, at);
        }
        let mut basic = Vec::with_capacity(m);
        for i in 0..m {
            let j = lp.n_structural + i;
            status[j] = Status::Basic { row: i };
            basic.push(j);
        }

        let mut solver = Solver {
            lp,
            basis: Basis::all_logical(m),
            basic,
            status,
            z,
            col: vec![0.0; m],
            alpha: Vec::new(),
            y: Vec::new(),
            cost_b: vec![0.0; m],
            rhs: vec![0.0; m],
        };
        solver.recompute_basic_values();
        solver
    }

    /// Recompute the basic variables from the nonbasic ones: `z_B = B^-1 (-N z_N)`.
    ///
    /// Called at the start and after every refactorization, so drift from the
    /// incremental per-pivot updates cannot accumulate indefinitely.
    fn recompute_basic_values(&mut self) {
        let lp = self.lp;
        self.rhs.fill(0.0);
        for j in 0..lp.n_total() {
            if let Status::NonBasic(_) = self.status[j] {
                let v = self.z[j];
                if v == 0.0 {
                    continue;
                }
                if j < lp.n_structural {
                    let (rows, vals) = lp.matrix.column(j);
                    for (&i, &a) in rows.iter().zip(vals) {
                        self.rhs[i] -= a * v;
                    }
                } else {
                    self.rhs[j - lp.n_structural] += v;
                }
            }
        }
        let mut out = std::mem::take(&mut self.alpha);
        self.basis.ftran(&self.rhs, &mut out);
        for (i, &v) in out.iter().enumerate() {
            self.z[self.basic[i]] = v;
        }
        self.alpha = out;
    }

    /// Total bound violation across the basic variables.
    fn infeasibility(&self) -> f64 {
        let tol = self.lp.tol.primal_feasibility;
        let mut total = 0.0;
        for (i, &j) in self.basic.iter().enumerate() {
            let _ = i;
            let v = self.z[j];
            if v < self.lp.lower[j] - tol {
                total += self.lp.lower[j] - v;
            } else if v > self.lp.upper[j] + tol {
                total += v - self.lp.upper[j];
            }
        }
        total
    }

    /// The phase-1 gradient of the basic costs, or the true costs in phase 2.
    fn load_basic_costs(&mut self, phase_one: bool) {
        let tol = self.lp.tol.primal_feasibility;
        for i in 0..self.lp.m {
            let j = self.basic[i];
            self.cost_b[i] = if phase_one {
                let v = self.z[j];
                if v < self.lp.lower[j] - tol {
                    -1.0
                } else if v > self.lp.upper[j] + tol {
                    1.0
                } else {
                    0.0
                }
            } else {
                self.lp.cost[j]
            };
        }
    }

    /// Reduced cost of nonbasic `j` given the current duals.
    fn reduced_cost(&mut self, j: usize, phase_one: bool) -> f64 {
        let cj = if phase_one { 0.0 } else { self.lp.cost[j] };
        let dot = if j < self.lp.n_structural {
            let (rows, vals) = self.lp.matrix.column(j);
            rows.iter()
                .zip(vals)
                .map(|(&i, &a)| self.y[i] * a)
                .sum::<f64>()
        } else {
            -self.y[j - self.lp.n_structural]
        };
        cj - dot
    }

    /// Choose an entering variable, returning it with the direction to move.
    ///
    /// Dantzig pricing (steepest reduced cost) normally, switching to Bland's rule
    /// — lowest index that improves — once the solve has stalled. Dantzig is the
    /// faster rule but can cycle on degenerate vertices; Bland's cannot, so falling
    /// back to it guarantees termination at the cost of some speed.
    fn price(&mut self, phase_one: bool, bland: bool) -> Option<(usize, f64)> {
        let tol = self.lp.tol.dual_feasibility;
        let mut best: Option<(usize, f64)> = None;
        let mut best_score = 0.0;

        for j in 0..self.lp.n_total() {
            let Status::NonBasic(at) = self.status[j] else {
                continue;
            };
            // A fixed variable has nowhere to go.
            if self.lp.lower[j] == self.lp.upper[j] {
                continue;
            }
            let d = self.reduced_cost(j, phase_one);
            let sigma = match at {
                At::Lower if d < -tol => 1.0,
                At::Upper if d > tol => -1.0,
                At::Zero if d < -tol => 1.0,
                At::Zero if d > tol => -1.0,
                _ => continue,
            };
            if bland {
                return Some((j, sigma));
            }
            let score = d.abs();
            if score > best_score {
                best_score = score;
                best = Some((j, sigma));
            }
        }
        best
    }

    /// Ratio test: how far can the entering variable move, and what stops it.
    ///
    /// In phase 2 every basic variable must stay inside its bounds. In phase 1 the
    /// step also stops when an *infeasible* basic reaches the bound it is violating,
    /// because the objective's slope changes there.
    ///
    /// Ties are broken on the largest pivot, which is the numerically safest choice.
    /// Under `bland` they are broken on the lowest variable index instead: Bland's
    /// rule only guarantees termination if it governs *both* the entering and the
    /// leaving choice, so pairing a Bland entering rule with a largest-pivot ratio
    /// test does not prevent cycling.
    fn ratio_test(&self, entering: usize, sigma: f64, phase_one: bool, bland: bool) -> Step {
        let lp = self.lp;
        let tol = lp.tol.primal_feasibility;
        let pivot_tol = lp.tol.pivot;

        // The entering variable's own range, if both bounds are finite.
        let mut best = if lp.lower[entering].is_finite() && lp.upper[entering].is_finite() {
            lp.upper[entering] - lp.lower[entering]
        } else {
            f64::INFINITY
        };
        let mut chosen: Option<(usize, At)> = None;
        // Among ties, prefer the largest pivot: it is the numerically safest choice
        // and is the standard cheap defence against degenerate stalling.
        let mut best_pivot = 0.0;

        for i in 0..lp.m {
            let beta = self.alpha[i] * sigma;
            if beta.abs() <= pivot_tol {
                continue;
            }
            let j = self.basic[i];
            let v = self.z[j];
            let (lo, hi) = (lp.lower[j], lp.upper[j]);

            // The bound this basic variable is heading towards, and hence the
            // breakpoint that limits the step. `beta > 0` means the variable falls.
            //
            // The bound it lands on is carried along with the target rather than
            // inferred from the direction of travel: in phase 1 an infeasible
            // variable travels *towards* the bound it violates, so a rising variable
            // can land on its lower bound -- the opposite of the phase-2 case.
            let target = if phase_one {
                if v < lo - tol {
                    // Infeasible below: rising reaches `lo`; falling only worsens it,
                    // which phase 1 permits, so it imposes no limit.
                    if beta < 0.0 {
                        Some((lo, At::Lower))
                    } else {
                        None
                    }
                } else if v > hi + tol {
                    if beta > 0.0 {
                        Some((hi, At::Upper))
                    } else {
                        None
                    }
                } else if beta > 0.0 {
                    Some((lo, At::Lower))
                } else {
                    Some((hi, At::Upper))
                }
            } else if beta > 0.0 {
                Some((lo, At::Lower))
            } else {
                Some((hi, At::Upper))
            };

            let Some((target, to)) = target else { continue };
            if !target.is_finite() {
                continue;
            }

            // Clamp at zero: a basic variable already a hair past its bound would
            // otherwise produce a small negative ratio and a backwards step.
            let ratio = ((v - target) / beta).max(0.0);
            let pivot = self.alpha[i].abs();
            let ties = ratio <= best + 1e-12;
            let wins_tie = if bland {
                chosen.is_none_or(|(prev, _)| j < self.basic[prev])
            } else {
                pivot > best_pivot
            };
            if ratio < best - 1e-12 || (ties && wins_tie) {
                if ratio < best {
                    best = ratio;
                }
                chosen = Some((i, to));
                best_pivot = pivot;
            }
        }

        match chosen {
            Some((row, to)) => Step::Pivot {
                leaving_row: row,
                step: best,
                to,
            },
            None if best.is_finite() => Step::BoundFlip { step: best },
            None => Step::Unbounded,
        }
    }

    /// Apply a step of `sigma * step` on the entering variable.
    fn take_step(&mut self, entering: usize, sigma: f64, step: f64) {
        if step == 0.0 {
            return;
        }
        let delta = sigma * step;
        self.z[entering] += delta;
        for i in 0..self.lp.m {
            let j = self.basic[i];
            self.z[j] -= self.alpha[i] * delta;
        }
    }

    /// Rebuild the basis inverse, repairing singular positions with logicals.
    fn refactorize(&mut self) -> Result<(), LpStatus> {
        let lp = self.lp;
        for _ in 0..lp.m + 1 {
            let columns: Vec<Vec<f64>> = self
                .basic
                .iter()
                .map(|&j| {
                    let mut c = vec![0.0; lp.m];
                    lp.column_into(j, &mut c);
                    c
                })
                .collect();
            match self.basis.refactorize(&columns, lp.tol.pivot) {
                Ok(()) => {
                    self.recompute_basic_values();
                    return Ok(());
                }
                Err(BasisError::Singular { row }) => {
                    // Swap the offending position for its logical, which is always
                    // available and keeps the basis nonsingular by construction.
                    let logical = lp.n_structural + row;
                    if self.basic[row] == logical {
                        return Err(LpStatus::Infeasible);
                    }
                    let displaced = self.basic[row];
                    let at = if lp.lower[displaced].is_finite() {
                        At::Lower
                    } else {
                        At::Upper
                    };
                    self.status[displaced] = Status::NonBasic(at);
                    self.z[displaced] = lp.value_at(displaced, at);
                    self.basic[row] = logical;
                    self.status[logical] = Status::Basic { row };
                }
            }
        }
        Err(LpStatus::Infeasible)
    }

    fn finish(&self, status: LpStatus, iterations: usize) -> LpSolution {
        let x = self.z[..self.lp.n_structural].to_vec();
        let objective = if status == LpStatus::Optimal {
            (0..self.lp.n_structural)
                .map(|j| self.lp.cost[j] * self.z[j])
                .sum()
        } else {
            f64::NAN
        };
        LpSolution {
            status,
            objective,
            x,
            iterations,
        }
    }

    fn run(mut self, max_iterations: usize) -> LpSolution {
        const REFACTOR_EVERY: usize = 50;
        // Degenerate (zero-length) steps in a row before switching to Bland's rule.
        const STALL_LIMIT: usize = 100;

        let mut iterations = 0usize;
        let mut stalled = 0usize;

        if self.refactorize().is_err() {
            return self.finish(LpStatus::Infeasible, 0);
        }

        while iterations < max_iterations {
            let infeasibility = self.infeasibility();
            let phase_one = infeasibility > self.lp.tol.primal_feasibility;

            self.load_basic_costs(phase_one);
            let mut y = std::mem::take(&mut self.y);
            self.basis.btran(&self.cost_b, &mut y);
            self.y = y;

            let bland = stalled > STALL_LIMIT;
            let Some((entering, sigma)) = self.price(phase_one, bland) else {
                // No improving column. In phase 2 that is optimality; in phase 1 it
                // means the infeasibility cannot be reduced further, so the LP has no
                // feasible point at all.
                let status = if phase_one {
                    LpStatus::Infeasible
                } else {
                    LpStatus::Optimal
                };
                return self.finish(status, iterations);
            };

            let mut col = std::mem::take(&mut self.col);
            self.lp.column_into(entering, &mut col);
            let mut alpha = std::mem::take(&mut self.alpha);
            self.basis.ftran(&col, &mut alpha);
            self.col = col;
            self.alpha = alpha;

            match self.ratio_test(entering, sigma, phase_one, bland) {
                Step::Unbounded => {
                    // Phase 1's objective is bounded below by zero, so an unbounded ray
                    // there means the ratio test lost a breakpoint to rounding rather
                    // than that the problem is genuinely unbounded.
                    let status = if phase_one {
                        LpStatus::Infeasible
                    } else {
                        LpStatus::Unbounded
                    };
                    return self.finish(status, iterations);
                }
                Step::BoundFlip { step } => {
                    self.take_step(entering, sigma, step);
                    let at = if sigma > 0.0 { At::Upper } else { At::Lower };
                    self.status[entering] = Status::NonBasic(at);
                    self.z[entering] = self.lp.value_at(entering, at);
                    stalled = if step.abs() <= 1e-12 { stalled + 1 } else { 0 };
                }
                Step::Pivot {
                    leaving_row,
                    step,
                    to,
                } => {
                    self.take_step(entering, sigma, step);
                    let leaving = self.basic[leaving_row];

                    self.basis.update(&self.alpha, leaving_row);
                    self.basic[leaving_row] = entering;
                    self.status[entering] = Status::Basic { row: leaving_row };
                    self.status[leaving] = Status::NonBasic(to);
                    self.z[leaving] = self.lp.value_at(leaving, to);

                    stalled = if step.abs() <= 1e-12 { stalled + 1 } else { 0 };
                }
            }

            iterations += 1;
            if self.basis.updates() >= REFACTOR_EVERY && self.refactorize().is_err() {
                return self.finish(LpStatus::Infeasible, iterations);
            }
        }

        self.finish(LpStatus::IterationLimit, iterations)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{RowSense, Sense};
    use crate::sparse::SparseMatrix;

    /// Build a problem from dense rows, for readable test cases.
    fn problem(obj: &[f64], rows: &[(&[f64], RowSense, f64)]) -> Problem {
        let n = obj.len();
        let m = rows.len();
        let mut triplets = Vec::new();
        let (mut row_lb, mut row_ub) = (Vec::new(), Vec::new());
        for (i, (coeffs, sense, rhs)) in rows.iter().enumerate() {
            for (j, &v) in coeffs.iter().enumerate() {
                if v != 0.0 {
                    triplets.push((i, j, v));
                }
            }
            let (lo, hi) = sense.bounds(*rhs);
            row_lb.push(lo);
            row_ub.push(hi);
        }
        Problem {
            name: "test".into(),
            sense: Sense::Minimize,
            obj: obj.to_vec(),
            obj_offset: 0.0,
            matrix: SparseMatrix::from_triplets(m, n, triplets),
            row_lb,
            row_ub,
            col_lb: vec![0.0; n],
            col_ub: vec![1.0; n],
            col_names: (0..n).map(|j| format!("x{j}")).collect(),
            row_names: (0..m).map(|i| format!("c{i}")).collect(),
        }
    }

    fn solve(obj: &[f64], rows: &[(&[f64], RowSense, f64)]) -> LpSolution {
        Lp::relaxation(&problem(obj, rows)).solve()
    }

    /// Every row's activity lies within its bounds, and every variable within its.
    fn assert_feasible(p: &Problem, x: &[f64]) {
        let csr = p.matrix.to_csr();
        for i in 0..p.n_rows() {
            let (cols, vals) = csr.column(i);
            let activity: f64 = cols.iter().zip(vals).map(|(&j, &a)| a * x[j]).sum();
            assert!(
                activity >= p.row_lb[i] - 1e-7 && activity <= p.row_ub[i] + 1e-7,
                "row {i}: activity {activity} outside [{}, {}]",
                p.row_lb[i],
                p.row_ub[i]
            );
        }
        for (j, &v) in x.iter().enumerate() {
            assert!(
                v >= p.col_lb[j] - 1e-9 && v <= p.col_ub[j] + 1e-9,
                "column {j}: {v} outside [{}, {}]",
                p.col_lb[j],
                p.col_ub[j]
            );
        }
    }

    #[test]
    fn ge_row_forces_a_fractional_value() {
        // Regression: phase 1 parked a variable rising towards its *lower* bound at
        // its upper bound instead, so `>=` rows were silently left violated and the
        // solver returned the infeasible origin as optimal.
        let s = solve(&[1.0], &[(&[1.0], RowSense::Ge, 0.5)]);
        assert_eq!(s.status, LpStatus::Optimal);
        assert!((s.objective - 0.5).abs() < 1e-9, "{}", s.objective);
        assert!((s.x[0] - 0.5).abs() < 1e-9);
    }

    #[test]
    fn picks_the_cheaper_column_to_cover_a_row() {
        let s = solve(&[2.0, 3.0], &[(&[1.0, 1.0], RowSense::Ge, 1.0)]);
        assert_eq!(s.status, LpStatus::Optimal);
        assert!((s.objective - 2.0).abs() < 1e-9, "{}", s.objective);
    }

    #[test]
    fn le_row_with_a_negative_cost_pushes_to_the_upper_bound() {
        let s = solve(&[-1.0], &[(&[1.0], RowSense::Le, 1.0)]);
        assert_eq!(s.status, LpStatus::Optimal);
        assert!((s.objective + 1.0).abs() < 1e-9, "{}", s.objective);
        assert!((s.x[0] - 1.0).abs() < 1e-9);
    }

    #[test]
    fn equality_row_is_satisfied_exactly() {
        let s = solve(&[1.0, 1.0], &[(&[1.0, 1.0], RowSense::Eq, 1.0)]);
        assert_eq!(s.status, LpStatus::Optimal);
        assert!((s.objective - 1.0).abs() < 1e-9, "{}", s.objective);
    }

    #[test]
    fn relaxation_is_fractional_where_the_integer_answer_is_not() {
        // 2a + 3b >= 4 with a, b in [0,1]: the LP takes b = 1 and a = 0.5 for 6,
        // while the integer optimum is a = b = 1 for 11.
        let p = problem(&[10.0, 1.0], &[(&[2.0, 3.0], RowSense::Ge, 4.0)]);
        let s = Lp::relaxation(&p).solve();
        assert_eq!(s.status, LpStatus::Optimal);
        assert!((s.objective - 6.0).abs() < 1e-9, "{}", s.objective);
        assert_feasible(&p, &s.x);
    }

    #[test]
    fn detects_an_infeasible_pair_of_rows() {
        let s = solve(
            &[1.0],
            &[(&[1.0], RowSense::Ge, 1.0), (&[1.0], RowSense::Le, 0.0)],
        );
        assert_eq!(s.status, LpStatus::Infeasible);
    }

    #[test]
    fn respects_bounds_narrowed_by_branching() {
        // Fixing x0 to 1 is how branching will present a subproblem.
        let mut p = problem(&[1.0, 5.0], &[(&[1.0, 1.0], RowSense::Ge, 1.0)]);
        p.col_lb[0] = 1.0;
        let s = Lp::relaxation(&p).solve();
        assert_eq!(s.status, LpStatus::Optimal);
        assert!((s.objective - 1.0).abs() < 1e-9, "{}", s.objective);

        // Fixing it to 0 instead forces the expensive column in.
        let mut p = problem(&[1.0, 5.0], &[(&[1.0, 1.0], RowSense::Ge, 1.0)]);
        p.col_ub[0] = 0.0;
        let s = Lp::relaxation(&p).solve();
        assert_eq!(s.status, LpStatus::Optimal);
        assert!((s.objective - 5.0).abs() < 1e-9, "{}", s.objective);
    }

    #[test]
    fn terminates_on_a_degenerate_problem() {
        // Regression: Bland's rule was applied to the entering variable but not the
        // leaving one, so a degenerate vertex could cycle until the iteration limit.
        // Many identical rows make degeneracy near-certain.
        let rows: Vec<(&[f64], RowSense, f64)> = (0..12)
            .map(|_| (&[1.0, 1.0, 1.0][..], RowSense::Ge, 1.0))
            .collect();
        let p = problem(&[1.0, 1.0, 1.0], &rows);
        let s = Lp::relaxation(&p).solve_with_limit(5_000);
        assert_eq!(
            s.status,
            LpStatus::Optimal,
            "took {} iterations",
            s.iterations
        );
        assert!((s.objective - 1.0).abs() < 1e-9, "{}", s.objective);
        assert_feasible(&p, &s.x);
    }

    #[test]
    fn an_empty_row_set_leaves_the_objective_at_its_best_bound() {
        // No constraints: each variable goes to whichever bound its cost prefers.
        let s = solve(&[1.0, -2.0], &[]);
        assert_eq!(s.status, LpStatus::Optimal);
        assert!((s.objective + 2.0).abs() < 1e-9, "{}", s.objective);
        assert!((s.x[0] - 0.0).abs() < 1e-9 && (s.x[1] - 1.0).abs() < 1e-9);
    }

    #[test]
    fn range_row_is_honoured_on_both_sides() {
        let mut p = problem(&[-1.0], &[(&[1.0], RowSense::Ge, 0.25)]);
        p.row_ub[0] = 0.75;
        let s = Lp::relaxation(&p).solve();
        assert_eq!(s.status, LpStatus::Optimal);
        // Cost is negative, so x rises until the row's upper bound stops it.
        assert!((s.x[0] - 0.75).abs() < 1e-9, "{:?}", s.x);
        assert_feasible(&p, &s.x);
    }
}
