//! Probing in presolve: suppose each binary at zero and at one, propagate fully, and keep
//! what both suppositions agree on.
//!
//! [`crate::presolve`]'s own probing fixes binaries and nothing else, and the node
//! propagator in [`crate::propagate`] recomputes a row's activity every time one of its
//! columns moves and leaves alone any row with an infinite bound in it. Neither is what
//! probing a whole model needs: thousands of suppositions, each undone, over models whose
//! rows are long and whose continuous columns are unbounded. So this keeps every row's
//! activity as it goes, with a count of the infinite terms in place of their sum (a row
//! with one infinite term still bounds that term), and undoes a supposition from a trail
//! rather than from a copy.
//!
//! HiGHS does the same, and on `ex10` it is what takes the model from 15936 columns to 6.

use crate::cuts::{Conflicts, literal};
use crate::model::{Problem, VarType};
use crate::presolve::Link;
use crate::sparse::SparseMatrix;

/// What [`probe_bounds`] did.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct ProbeStats {
    pub probes: usize,
    pub fixed: usize,
    pub tightened: usize,
    pub links: usize,
}

/// Matrix entries one supposition may read before what it has found is taken as it is.
const PROBE_REACH: usize = 1_000_000;

/// Matrix entries a first [`probe_bounds`] may read, as a multiple of the model's.
const PROBE_TOTAL: usize = 300;

/// The least it may read whatever the model's size, and the most: a second or two, and
/// about five on `ex10`, whose 17680 binaries probed in full would take a minute.
const PROBE_TOTAL_RANGE: (usize, usize) = (200_000_000, 400_000_000);

/// Most implications [`probe_bounds`] keeps from its suppositions, counting each way.
const LEARNED_CAPACITY: usize = 20_000_000;

/// Share of the total that may pass in a row without anything found.
const PROBE_PATIENCE: usize = 8;

/// Looser than the model's tolerance, since a supposition refuted by rounding fixes a
/// binary for good: two rows nudging a continuous bound back and forth reach a
/// contradiction of 1e-9 eventually, and one of 1e-6 they do not.
const TOLERANCE: f64 = 1e-6;

/// The share of its range a continuous bound must move by to count, as in
/// `presolve::propagate_rows`; less, and two rows can pass it back and forth for ever.
const CONTINUOUS_STEP: f64 = 1e-3;

/// Probe every binary both ways, with every row's implications on every column.
///
/// - refuted, the binary takes the other value, with everything that value forces;
/// - a column both suppositions fix to one value is that value;
/// - one they fix to two values follows the binary, a [`Link`] for [`crate::reduce`];
/// - an integer column they only bound lies within the wider of the two bounds.
///
/// The bounds proved are written into `problem`: integer bounds, and continuous columns
/// that are fixed. A continuous bound that is only tighter is left out, since the rows
/// imply it anyway and the reductions after this read a column's bounds as given,
/// finding it free only when its rows bound it at least as tightly as they do.
///
/// `None` means the model is infeasible.
///
/// `share` scales the budget, so that a later pass over a model already probed spends
/// less than the first.
pub fn probe_bounds(problem: &mut Problem, share: f64) -> Option<(Vec<Link>, ProbeStats)> {
    let n = problem.n_cols();
    let mut stats = ProbeStats::default();
    let mut links = Vec::new();
    let mut prober = Prober::new(problem);
    // Integer bounds creeping one step at a time through a wide range would run long.
    if prober.propagate(PROBE_REACH.max(10 * problem.matrix.nnz())).is_err() {
        return None;
    }
    prober.commit();

    // Binaries with implications both ways first, as HiGHS orders them: one that
    // forces something only at one value can only be refuted, never linked, and the
    // product finds the ones whose two sides are both busy.
    let mut candidates: Vec<(usize, (usize, usize, u64))> = (0..n)
        .filter(|&j| is_binary(problem, j) && prober.lb[j] < prober.ub[j])
        .map(|j| {
            let up = prober.conflicts.degree(literal(j, true));
            let down = prober.conflicts.degree(literal(j, false));
            let mixed = (j as u64).wrapping_mul(0x9e37_79b9_7f4a_7c15).rotate_left(17);
            (j, ((up * down).min(5000), (up + down).min(100), mixed))
        })
        .collect();
    candidates.sort_by(|a, b| b.1.0.cmp(&a.1.0).then(b.1.1.cmp(&a.1.1)).then(a.1.2.cmp(&b.1.2)));

    let nnz = problem.matrix.nnz();
    let total = ((nnz * PROBE_TOTAL).clamp(PROBE_TOTAL_RANGE.0, PROBE_TOTAL_RANGE.1) as f64 * share) as usize;
    let patience = total / PROBE_PATIENCE;
    let (mut spent, mut found_at) = (0usize, 0usize);
    // Each supposition's bounds on the columns it moved, the first read during the second.
    let mut at_zero: Vec<(f64, f64)> = vec![(f64::NAN, f64::NAN); n];
    let mut moved_zero: Vec<usize> = Vec::new();
    let mut moved_one: Vec<(usize, f64, f64)> = Vec::new();
    let mut proved: Vec<(usize, f64, f64)> = Vec::new();

    for &(y, _) in &candidates {
        if spent > total || spent - found_at > patience {
            break;
        }
        if prober.lb[y] >= prober.ub[y] {
            continue;
        }
        stats.probes += 1;

        let zero = prober.suppose(y, 0.0);
        spent += prober.work;
        moved_zero.clear();
        if zero {
            for k in prober.moved() {
                at_zero[k] = (prober.lb[k], prober.ub[k]);
                moved_zero.push(k);
            }
            prober.learn(y, false, &moved_zero);
        }
        prober.undo();
        let one = prober.suppose(y, 1.0);
        spent += prober.work;
        moved_one.clear();
        if one {
            let moved = prober.moved();
            for &k in &moved {
                moved_one.push((k, prober.lb[k], prober.ub[k]));
            }
            prober.learn(y, true, &moved);
        }
        prober.undo();

        proved.clear();
        match (zero, one) {
            (false, false) => return None,
            (false, true) => proved.push((y, 1.0, 1.0)),
            (true, false) => proved.push((y, 0.0, 0.0)),
            (true, true) => {
                for &(k, l1, u1) in &moved_one {
                    let (l0, u0) = at_zero[k];
                    if k == y || l0.is_nan() {
                        continue;
                    }
                    let fixed_both = u0 - l0 <= TOLERANCE && u1 - l1 <= TOLERANCE;
                    if fixed_both && (l0 - l1).abs() > TOLERANCE {
                        links.push(Link { column: k, switch: y, at_zero: l0, at_one: l1 });
                        stats.links += 1;
                        found_at = spent;
                    } else {
                        proved.push((k, l0.min(l1), u0.max(u1)));
                    }
                }
            }
        }
        for &k in &moved_zero {
            at_zero[k] = (f64::NAN, f64::NAN);
        }

        let mut any = false;
        for &(k, l, u) in &proved {
            let integer = problem.col_type[k] == VarType::Integer;
            let fixing = u - l <= TOLERANCE;
            if !integer && !fixing {
                continue;
            }
            let (l, u) = if integer { (l.round(), u.round()) } else { (l, u) };
            let raise = l > prober.lb[k] + TOLERANCE;
            let lower = u < prober.ub[k] - TOLERANCE;
            if !raise && !lower {
                continue;
            }
            if fixing || (raise && lower && u - l <= TOLERANCE) {
                stats.fixed += 1;
            } else {
                stats.tightened += 1;
            }
            let (l, u) = (l.max(prober.lb[k]), u.min(prober.ub[k]));
            if prober.set_bounds(k, l, u).is_err() {
                return None;
            }
            any = true;
        }
        if any {
            found_at = spent;
            prober.work = 0;
            if prober.propagate(PROBE_REACH).is_err() {
                return None;
            }
            spent += prober.work;
            prober.commit();
        }
    }

    if !enumerate_rows(&mut prober, total / 2, &mut links, &mut stats) {
        return None;
    }

    let Prober { lb, ub, .. } = prober;
    for j in 0..n {
        let integer = problem.col_type[j] == VarType::Integer;
        let (l, u) = (lb[j], ub[j]);
        if integer {
            problem.col_lb[j] = problem.col_lb[j].max(l);
            problem.col_ub[j] = problem.col_ub[j].min(u);
        } else if u - l <= TOLERANCE && (problem.col_lb[j] < l || problem.col_ub[j] > u) {
            let value = 0.5 * (l + u);
            problem.col_lb[j] = value;
            problem.col_ub[j] = value;
        }
    }
    Some((links, stats))
}

/// Longest row [`enumerate_rows`] reads, in free columns.
const ENUMERATE_LENGTH: usize = 32;

/// Most assignments of a row's columns [`enumerate_rows`] collects before giving the
/// row up.
const ENUMERATE_LEAVES: usize = 64;

/// Matrix entries one row's enumeration may read.
const ENUMERATE_REACH: usize = 20_000_000;

/// For each short row over binaries, every assignment of its columns that propagation
/// does not refute, and what holds in all of them.
///
/// Probing one binary splits the model in two; a row splits it as many ways as the row
/// has assignments, and a column can be pinned in every one of them where neither side
/// of any single probe pins it. A set partitioning row of twenty columns has twenty
/// assignments, each a probe that also fixes the other nineteen, and `x <= y + z`
/// has five. HiGHS calls this enumeration, and on `ex10` it is the rule that does the
/// most after probing.
///
/// What holds in every surviving assignment:
/// - a column within the same bounds in all of them is within those bounds;
/// - one fixed to two values, each in step with a column of the row, follows it.
///
/// Returns false when no assignment survives, so the model is infeasible.
fn enumerate_rows(prober: &mut Prober, total: usize, links: &mut Vec<Link>, stats: &mut ProbeStats) -> bool {
    let problem = prober.problem;
    let m = problem.n_rows();
    let mut order: Vec<(usize, usize)> = (0..m)
        .filter_map(|i| {
            let (columns, _) = prober.rows.column(i);
            let free = columns.iter().filter(|&&k| prober.lb[k] < prober.ub[k]).count();
            let binary = columns.iter().all(|&k| is_binary(problem, k));
            (binary && (2..=ENUMERATE_LENGTH).contains(&free)).then_some((free, i))
        })
        .collect();
    order.sort_unstable();
    let (mut spent, mut found_at) = (0usize, 0usize);
    let patience = total / PROBE_PATIENCE * 2;
    let mut leaves: Vec<Vec<(usize, f64, f64)>> = Vec::new();
    // For each column, how many leaves moved it, the union of its bounds in them, and
    // whether every one of them fixed it.
    const NONE: (usize, f64, f64, bool) = (0, f64::INFINITY, f64::NEG_INFINITY, true);
    let mut union: Vec<(usize, f64, f64, bool)> = vec![NONE; problem.n_cols()];
    let mut seen: Vec<usize> = Vec::new();
    // Each moved column's lower bound in each leaf, by its place in `seen`.
    let mut place: Vec<usize> = vec![usize::MAX; problem.n_cols()];
    let mut values: Vec<Vec<f64>> = Vec::new();
    for (_, i) in order {
        if spent > total || spent - found_at > patience {
            break;
        }
        let columns: Vec<usize> =
            prober.rows.column(i).0.iter().copied().filter(|&k| prober.lb[k] < prober.ub[k]).collect();
        if columns.len() < 2 {
            continue;
        }
        leaves.clear();
        prober.work = 0;
        let complete = descend(prober, &columns, 0, &mut leaves);
        spent += prober.work;
        prober.undo();
        if !complete {
            continue;
        }
        if leaves.is_empty() {
            return false;
        }
        seen.clear();
        for (index, leaf) in leaves.iter().enumerate() {
            for &(k, l, u) in leaf {
                let entry = &mut union[k];
                if entry.0 == 0 {
                    place[k] = seen.len();
                    if values.len() <= seen.len() {
                        values.push(Vec::new());
                    }
                    values[seen.len()].clear();
                    values[seen.len()].resize(leaves.len(), f64::NAN);
                    seen.push(k);
                }
                entry.0 += 1;
                entry.1 = entry.1.min(l);
                entry.2 = entry.2.max(u);
                entry.3 &= u - l <= TOLERANCE;
                values[place[k]][index] = l;
            }
        }
        let mut changed = false;
        for &k in &seen {
            let (count, l, u, fixed) = std::mem::replace(&mut union[k], NONE);
            // A leaf that did not move it left it as it was, which proves nothing.
            if count < leaves.len() || prober.lb[k] >= prober.ub[k] {
                continue;
            }
            let integer = problem.col_type[k] == VarType::Integer;
            if !fixed {
                if integer && (l > prober.lb[k] + 0.5 || u < prober.ub[k] - 0.5) {
                    stats.tightened += 1;
                    if prober.set_bounds(k, l.round().max(prober.lb[k]), u.round().min(prober.ub[k])).is_err() {
                        return false;
                    }
                    changed = true;
                }
                continue;
            }
            // Fixed in every leaf: to one value, or to two, perhaps in step with a column
            // of the row.
            let at_k = &values[place[k]];
            let first = at_k[0];
            let other = at_k.iter().copied().find(|v| (v - first).abs() > TOLERANCE);
            match other {
                None => {
                    stats.fixed += 1;
                    if prober.set_bounds(k, first, first).is_err() {
                        return false;
                    }
                    changed = true;
                }
                Some(second) => {
                    let two = at_k.iter().all(|v| (v - first).abs() <= TOLERANCE || (v - second).abs() <= TOLERANCE);
                    let follows = columns.iter().copied().filter(|&y| y != k && two).find_map(|y| {
                        let mut at = [None, None];
                        // `place` outlives the row it was set for; `seen` says whether it is this one's.
                        if seen.get(place[y]) != Some(&y) {
                            return None;
                        }
                        let at_y = &values[place[y]];
                        for (&side, &v) in at_y.iter().zip(at_k) {
                            if side.is_nan() {
                                return None;
                            }
                            let slot = &mut at[usize::from(side > 0.5)];
                            match slot {
                                None => *slot = Some(v),
                                Some(w) if (*w - v).abs() <= TOLERANCE => {}
                                Some(_) => return None,
                            }
                        }
                        Some((y, at[0]?, at[1]?))
                    });
                    if let Some((y, at_zero, at_one)) = follows {
                        links.push(Link { column: k, switch: y, at_zero, at_one });
                        stats.links += 1;
                        found_at = spent;
                    } else if integer && (l > prober.lb[k] + 0.5 || u < prober.ub[k] - 0.5) {
                        // Values in no step with the row: it lies between them.
                        stats.tightened += 1;
                        if prober.set_bounds(k, l.round().max(prober.lb[k]), u.round().min(prober.ub[k])).is_err() {
                            return false;
                        }
                        changed = true;
                    }
                }
            }
        }
        if changed {
            found_at = spent;
            prober.work = 0;
            if prober.propagate(PROBE_REACH).is_err() {
                return false;
            }
            spent += prober.work;
            prober.commit();
        }
    }
    true
}

/// Every assignment of `columns[from..]` propagation leaves standing, from the bounds as
/// they are, into `leaves` as the bounds of every column moved since the commit. False
/// when there are too many to collect or the reading ran out, so that they are not all
/// there.
fn descend(prober: &mut Prober, columns: &[usize], from: usize, leaves: &mut Vec<Vec<(usize, f64, f64)>>) -> bool {
    let Some(at) = (from..columns.len()).find(|&a| prober.lb[columns[a]] < prober.ub[columns[a]]) else {
        if leaves.len() >= ENUMERATE_LEAVES {
            return false;
        }
        let moved = prober.moved();
        leaves.push(moved.into_iter().map(|k| (k, prober.lb[k], prober.ub[k])).collect());
        return true;
    };
    let k = columns[at];
    for value in [1.0, 0.0] {
        let mark = prober.trail.len();
        let survived = prober.set_bounds(k, value, value).is_ok() && prober.propagate(ENUMERATE_REACH).is_ok();
        if prober.work >= ENUMERATE_REACH {
            return false;
        }
        let complete = !survived || descend(prober, columns, at + 1, leaves);
        prober.undo_to(mark);
        if !complete {
            return false;
        }
    }
    true
}

fn is_binary(problem: &Problem, j: usize) -> bool {
    problem.col_type[j] == VarType::Integer && problem.col_lb[j] >= 0.0 && problem.col_ub[j] <= 1.0
}

/// The bounds of every column and the activity range of every row, with a trail to undo
/// to the state last committed.
struct Prober<'a> {
    problem: &'a Problem,
    rows: SparseMatrix,
    conflicts: Conflicts,
    lb: Vec<f64>,
    ub: Vec<f64>,
    /// The finite part of each row's least and greatest activity, and how many terms
    /// are infinite in each.
    min_finite: Vec<f64>,
    max_finite: Vec<f64>,
    min_infinite: Vec<u32>,
    max_infinite: Vec<u32>,
    /// The widest any one term of a row can range over under the committed bounds. A row
    /// whose slack is at least this cannot tighten anything, and is passed over unread.
    widest: Vec<f64>,
    /// Columns' bounds before each change since the commit, oldest first.
    trail: Vec<(u32, f64, f64)>,
    /// Rows' activities before their first change since the commit.
    row_trail: Vec<(u32, f64, f64, u32, u32)>,
    row_saved: Vec<bool>,
    queue: Vec<u32>,
    queued: Vec<bool>,
    /// Implications found by earlier suppositions, as conflicts: for each literal, the
    /// literals it excludes. Supposing `y = 1` fixed `x` to 0 says `y` and `x` exclude
    /// each other, which helps a later supposition `x = 1` that the rows alone would
    /// never take back to `y`.
    learned: Vec<Vec<u32>>,
    learned_count: usize,
    /// The binary supposed, and the columns its own conflicts fixed: those it excludes
    /// already, and learning them again would only make every later read longer.
    seed: usize,
    direct: Vec<bool>,
    /// Binaries fixed since their conflicts were last read.
    fixed: Vec<u32>,
    seen: Vec<bool>,
    /// Matrix entries read since the last supposition began.
    work: usize,
}

impl<'a> Prober<'a> {
    fn new(problem: &'a Problem) -> Self {
        let (n, m) = (problem.n_cols(), problem.n_rows());
        let mut prober = Prober {
            problem,
            rows: problem.matrix.to_csr(),
            conflicts: Conflicts::of(problem),
            lb: problem.col_lb.clone(),
            ub: problem.col_ub.clone(),
            min_finite: vec![0.0; m],
            max_finite: vec![0.0; m],
            min_infinite: vec![0; m],
            max_infinite: vec![0; m],
            widest: vec![0.0; m],
            trail: Vec::new(),
            row_trail: Vec::new(),
            row_saved: vec![false; m],
            queue: (0..m as u32).collect(),
            queued: vec![true; m],
            fixed: Vec::new(),
            learned: vec![Vec::new(); 2 * n],
            learned_count: 0,
            seed: usize::MAX,
            direct: vec![false; n],
            seen: vec![false; n],
            work: 0,
        };
        for i in 0..m {
            prober.recompute(i);
        }
        prober
    }

    /// Row `i`'s activity and widest term, from scratch.
    fn recompute(&mut self, i: usize) {
        let (columns, values) = self.rows.column(i);
        let (mut min_f, mut max_f, mut min_i, mut max_i, mut widest) = (0.0, 0.0, 0, 0, 0.0f64);
        for (&k, &a) in columns.iter().zip(values) {
            let (low, high) = term(a, self.lb[k], self.ub[k]);
            if low.is_finite() {
                min_f += low;
            } else {
                min_i += 1;
            }
            if high.is_finite() {
                max_f += high;
            } else {
                max_i += 1;
            }
            widest = widest.max(high - low);
        }
        self.min_finite[i] = min_f;
        self.max_finite[i] = max_f;
        self.min_infinite[i] = min_i;
        self.max_infinite[i] = max_i;
        self.widest[i] = widest;
    }

    /// Set column `k`'s bounds, which must be no wider than they are, and bring its
    /// rows' activities along.
    fn set_bounds(&mut self, k: usize, l: f64, u: f64) -> Result<(), ()> {
        if l > u + TOLERANCE {
            return Err(());
        }
        let u = u.max(l);
        let (old_l, old_u) = (self.lb[k], self.ub[k]);
        self.trail.push((k as u32, old_l, old_u));
        self.shift(k, l, u);
        if u - l <= TOLERANCE && is_binary(self.problem, k) {
            self.fixed.push(k as u32);
        }
        Ok(())
    }

    /// Move column `k` to `[l, u]`, its rows' activities with it, and queue the rows.
    fn shift(&mut self, k: usize, l: f64, u: f64) {
        let (old_l, old_u) = (self.lb[k], self.ub[k]);
        self.lb[k] = l;
        self.ub[k] = u;
        let (rows, values) = self.problem.matrix.column(k);
        self.work += rows.len();
        for (&i, &a) in rows.iter().zip(values) {
            if !self.row_saved[i] {
                self.row_saved[i] = true;
                self.row_trail.push((
                    i as u32,
                    self.min_finite[i],
                    self.max_finite[i],
                    self.min_infinite[i],
                    self.max_infinite[i],
                ));
            }
            let (old_low, old_high) = term(a, old_l, old_u);
            let (low, high) = term(a, l, u);
            if old_low != low {
                if old_low.is_finite() {
                    self.min_finite[i] -= old_low;
                } else {
                    self.min_infinite[i] -= 1;
                }
                if low.is_finite() {
                    self.min_finite[i] += low;
                } else {
                    self.min_infinite[i] += 1;
                }
            }
            if old_high != high {
                if old_high.is_finite() {
                    self.max_finite[i] -= old_high;
                } else {
                    self.max_infinite[i] -= 1;
                }
                if high.is_finite() {
                    self.max_finite[i] += high;
                } else {
                    self.max_infinite[i] += 1;
                }
            }
            if !self.queued[i] {
                self.queued[i] = true;
                self.queue.push(i as u32);
            }
        }
    }

    /// Drop the changes after the first `mark` of the trail, keeping those before.
    /// The activities come back by subtraction, which [`Prober::undo`] and
    /// [`Prober::commit`] then make exact.
    fn undo_to(&mut self, mark: usize) {
        while self.trail.len() > mark {
            let (k, l, u) = self.trail.pop().unwrap();
            self.direct[k as usize] = false;
            self.shift(k as usize, l, u);
        }
        self.clear_queues();
    }

    /// Suppose binary `y` takes `value` and propagate. Whether that survived.
    fn suppose(&mut self, y: usize, value: f64) -> bool {
        self.work = 0;
        self.seed = y;
        let survived = self.set_bounds(y, value, value).is_ok() && self.propagate(PROBE_REACH).is_ok();
        self.seed = usize::MAX;
        survived
    }

    /// The columns moved since the commit, once each.
    fn moved(&mut self) -> Vec<usize> {
        let mut out = Vec::new();
        for &(k, _, _) in &self.trail {
            let k = k as usize;
            if !self.seen[k] {
                self.seen[k] = true;
                out.push(k);
            }
        }
        for &k in &out {
            self.seen[k] = false;
        }
        out
    }

    /// Keep the binaries the supposition `y = value` fixed as implications of it.
    fn learn(&mut self, y: usize, value: bool, moved: &[usize]) {
        let node = literal(y, value);
        for &k in moved {
            if std::mem::take(&mut self.direct[k]) {
                continue;
            }
            if k == y || !is_binary(self.problem, k) || self.ub[k] - self.lb[k] > TOLERANCE {
                continue;
            }
            if self.learned_count + 2 > LEARNED_CAPACITY {
                return;
            }
            // `y = value` forces `k` to its value: it excludes `k` at the other.
            let excluded = literal(k, self.lb[k] < 0.5);
            self.learned[node as usize].push(excluded);
            self.learned[excluded as usize].push(node);
            self.learned_count += 2;
        }
    }

    /// Drop everything since the commit.
    fn undo(&mut self) {
        while let Some((k, l, u)) = self.trail.pop() {
            self.direct[k as usize] = false;
            self.lb[k as usize] = l;
            self.ub[k as usize] = u;
        }
        for (i, min_f, max_f, min_i, max_i) in self.row_trail.drain(..) {
            let i = i as usize;
            self.min_finite[i] = min_f;
            self.max_finite[i] = max_f;
            self.min_infinite[i] = min_i;
            self.max_infinite[i] = max_i;
            self.row_saved[i] = false;
        }
        self.clear_queues();
    }

    /// Keep everything since the last commit, recomputing the rows it moved so that
    /// rounding in the running sums does not pile up from one probe to the next.
    fn commit(&mut self) {
        for &(k, _, _) in &self.trail {
            self.direct[k as usize] = false;
        }
        self.trail.clear();
        let rows: Vec<u32> = self.row_trail.drain(..).map(|r| r.0).collect();
        for i in rows {
            self.row_saved[i as usize] = false;
            self.recompute(i as usize);
        }
        self.clear_queues();
    }

    fn clear_queues(&mut self) {
        for &i in &self.queue {
            self.queued[i as usize] = false;
        }
        self.queue.clear();
        self.fixed.clear();
    }

    /// Everything the rows and the conflicts imply from the changes queued, until
    /// nothing moves or `budget` entries have been read. `Err` is a contradiction.
    fn propagate(&mut self, budget: usize) -> Result<(), ()> {
        let mut head = 0;
        loop {
            if self.work >= budget {
                return Ok(());
            }
            if let Some(j) = self.fixed.pop() {
                self.conflicts_of(j as usize)?;
                continue;
            }
            if head >= self.queue.len() {
                self.queue.clear();
                return Ok(());
            }
            let i = self.queue[head] as usize;
            head += 1;
            self.queued[i] = false;
            self.row(i)?;
            // Keep the queue from growing without end on a model that propagates long.
            if head > 4096 && head * 2 > self.queue.len() {
                self.queue.drain(..head);
                head = 0;
            }
        }
    }

    /// Fix whatever binary `j`'s value excludes.
    fn conflicts_of(&mut self, j: usize) -> Result<(), ()> {
        let node = literal(j, self.lb[j] > 0.5);
        let mut forced: Vec<(usize, f64)> = Vec::new();
        let (lb, ub) = (&self.lb, &self.ub);
        let mut clash = false;
        let mut read = 0usize;
        self.conflicts.for_each_adjacent(node, |other| {
            read += 1;
            let k = (other / 2) as usize;
            // An even literal is "column k is one"; being excluded forces it to zero.
            let value = if other % 2 == 0 { 0.0 } else { 1.0 };
            if value > ub[k] + TOLERANCE || value < lb[k] - TOLERANCE {
                clash = true;
                return false;
            }
            if ub[k] - lb[k] > TOLERANCE {
                forced.push((k, value));
            }
            true
        });
        if !clash {
            for &other in &self.learned[node as usize] {
                read += 1;
                let k = (other / 2) as usize;
                let value = if other % 2 == 0 { 0.0 } else { 1.0 };
                if value > ub[k] + TOLERANCE || value < lb[k] - TOLERANCE {
                    clash = true;
                    break;
                }
                if ub[k] - lb[k] > TOLERANCE {
                    forced.push((k, value));
                }
            }
        }
        self.work += read;
        if clash {
            return Err(());
        }
        let direct = j == self.seed;
        for (k, value) in forced {
            if direct {
                self.direct[k] = true;
            }
            if self.ub[k] - self.lb[k] > TOLERANCE {
                self.set_bounds(k, value, value)?;
            }
        }
        Ok(())
    }

    /// Tighten the columns of row `i` by what the rest of the row leaves them.
    fn row(&mut self, i: usize) -> Result<(), ()> {
        let (row_lb, row_ub) = (self.problem.row_lb[i], self.problem.row_ub[i]);
        let min_activity = if self.min_infinite[i] == 0 { self.min_finite[i] } else { f64::NEG_INFINITY };
        let max_activity = if self.max_infinite[i] == 0 { self.max_finite[i] } else { f64::INFINITY };
        let scale = 1.0 + row_ub.abs().min(row_lb.abs()).min(1e6);
        if min_activity > row_ub + TOLERANCE * scale || max_activity < row_lb - TOLERANCE * scale {
            return Err(());
        }
        // Read the row only if one side can bound something: it has at most one infinite
        // term, and less slack than its widest term spans.
        let upper = row_ub.is_finite()
            && self.min_infinite[i] <= 1
            && (self.min_infinite[i] == 1 || row_ub - min_activity < self.widest[i]);
        let lower = row_lb.is_finite()
            && self.max_infinite[i] <= 1
            && (self.max_infinite[i] == 1 || max_activity - row_lb < self.widest[i]);
        if !upper && !lower {
            self.work += 1;
            return Ok(());
        }
        let length = self.rows.column(i).0.len();
        self.work += length;
        for idx in 0..length {
            let (k, a) = {
                let (columns, values) = self.rows.column(i);
                (columns[idx], values[idx])
            };
            if a == 0.0 {
                continue;
            }
            let (low, high) = term(a, self.lb[k], self.ub[k]);
            let (mut new_l, mut new_u) = (self.lb[k], self.ub[k]);
            if upper {
                // The least the rest of the row can be.
                let rest = match self.min_infinite[i] {
                    0 => Some(self.min_finite[i] - low),
                    _ if !low.is_finite() => Some(self.min_finite[i]),
                    _ => None,
                };
                if let Some(rest) = rest {
                    let bound = (row_ub - rest) / a;
                    if a > 0.0 {
                        new_u = new_u.min(bound);
                    } else {
                        new_l = new_l.max(bound);
                    }
                }
            }
            if lower {
                let rest = match self.max_infinite[i] {
                    0 => Some(self.max_finite[i] - high),
                    _ if !high.is_finite() => Some(self.max_finite[i]),
                    _ => None,
                };
                if let Some(rest) = rest {
                    let bound = (row_lb - rest) / a;
                    if a > 0.0 {
                        new_l = new_l.max(bound);
                    } else {
                        new_u = new_u.min(bound);
                    }
                }
            }
            let (l, u) = (self.lb[k], self.ub[k]);
            let integer = self.problem.col_type[k] == VarType::Integer;
            let (raise, lower_ub) = if integer {
                new_l = (new_l - TOLERANCE).ceil();
                new_u = (new_u + TOLERANCE).floor();
                (new_l > l + 0.5, new_u < u - 0.5)
            } else {
                // A finite bound where there was none is always worth having; past that
                // only a real share of the range.
                let step = CONTINUOUS_STEP * if (u - l).is_finite() { (u - l).max(1.0) } else { 1.0 };
                (
                    new_l > l + step || (!l.is_finite() && new_l.is_finite()),
                    new_u < u - step || (!u.is_finite() && new_u.is_finite()),
                )
            };
            if !raise && !lower_ub {
                continue;
            }
            let l2 = if raise { new_l } else { l };
            let u2 = if lower_ub { new_u } else { u };
            // Rounded integer bounds that cross are crossed by a whole step; continuous
            // ones may cross by rounding, and are taken as touching.
            if l2 > u2 + if integer { 0.5 } else { TOLERANCE * (1.0 + u2.abs().min(1e6)) } {
                return Err(());
            }
            // A continuous column brought within tolerance of fixed is fixed.
            let (l2, u2) = if l2 > u2 { (u2, u2) } else { (l2, u2) };
            self.set_bounds(k, l2, u2)?;
        }
        Ok(())
    }
}

/// The least and greatest `a x` can be for `x` in `[l, u]`.
fn term(a: f64, l: f64, u: f64) -> (f64, f64) {
    if a >= 0.0 {
        (if l.is_finite() { a * l } else { f64::NEG_INFINITY }, if u.is_finite() { a * u } else { f64::INFINITY })
    } else {
        (if u.is_finite() { a * u } else { f64::NEG_INFINITY }, if l.is_finite() { a * l } else { f64::INFINITY })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn read(text: &str) -> Problem {
        Problem::from_text(text, crate::reader::Format::Lp).unwrap()
    }

    #[test]
    fn a_row_with_one_unbounded_term_still_bounds_it() {
        // z is free, so c1's activity has an infinite term either way; with x fixed it
        // is the only one, and c1 fixes it.
        let mut problem = read(
            "\
Minimize
 obj: z
Subject To
 c1: z - 2 x = 0
Bounds
 z free
Binaries
 x
End
",
        );
        let (links, stats) = probe_bounds(&mut problem, 1.0).unwrap();
        assert_eq!(stats.links, 1, "{stats:?}");
        assert_eq!((links[0].at_zero, links[0].at_one), (0.0, 2.0));
    }

    #[test]
    fn a_column_fixed_to_two_values_follows_the_switch() {
        // y = 1 forces w = 5 through c1; y = 0 forces w = 2 through c2.
        let text = "\
Minimize
 obj: w
Subject To
 c1: w - 3 y >= 2
 c2: w - 3 y <= 2
Bounds
 0 <= w <= 10
Binaries
 y
Generals
 w
End
";
        let mut problem = read(text);
        let (links, stats) = probe_bounds(&mut problem, 1.0).unwrap();
        assert_eq!(stats.links, 1, "{stats:?}");
        let link = links[0];
        assert_eq!((link.at_zero, link.at_one), (2.0, 5.0));
    }

    #[test]
    fn enumerating_a_row_finds_what_no_one_probe_does() {
        // k = x1 + x2, which the partition makes 1 - x3. Supposing x3 = 0 leaves
        // x1 + x2 = 1 and k between them unsettled; only taking the partition's three
        // cases one by one pins k in each.
        let problem = read(
            "\
Minimize
 obj: k + x3
Subject To
 p: x1 + x2 + x3 = 1
 a: k - x1 >= 0
 b: k - x2 >= 0
 c: k - x1 - x2 <= 0
Binaries
 k x1 x2 x3
End
",
        );
        let index = |name: &str| problem.col_names.iter().position(|c| c == name).unwrap();
        let (k, x3) = (index("k"), index("x3"));
        // On a prober that has learned nothing: probing x3 first and then k finds this
        // too, through what the first probe learned.
        let mut prober = Prober::new(&problem);
        let (mut links, mut stats) = (Vec::new(), ProbeStats::default());
        assert!(enumerate_rows(&mut prober, usize::MAX, &mut links, &mut stats));
        let pair = links.iter().find(|l| (l.column, l.switch) == (k, x3) || (l.column, l.switch) == (x3, k));
        let link = pair.unwrap_or_else(|| panic!("{stats:?} {links:?}"));
        assert_eq!((link.at_zero, link.at_one), (1.0, 0.0));
    }

    #[test]
    fn a_refuted_supposition_fixes_the_binary() {
        // y = 1 needs w >= 8, which c2 forbids.
        let text = "\
Minimize
 obj: w
Subject To
 c1: w - 8 y >= 0
 c2: w <= 5
Bounds
 0 <= w <= 10
Binaries
 y
End
";
        let mut problem = read(text);
        let (_, stats) = probe_bounds(&mut problem, 1.0).unwrap();
        let y = problem.col_names.iter().position(|c| c == "y").unwrap();
        assert_eq!(problem.col_ub[y], 0.0, "{stats:?}");
    }
}
