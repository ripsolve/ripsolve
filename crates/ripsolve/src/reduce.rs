//! Reductions that remove columns, and the postsolve that puts them back.
//!
//! `presolve` works in place: it fixes a column by setting `lb == ub` and frees a row
//! by widening its bounds, so nothing is renumbered and no postsolve is needed. That
//! is safe, and it is also why it cannot do the reductions the other open-source
//! solvers lean on hardest. Every one of those eliminates a column *in terms of the
//! others* -- a continuous column alone in an equation is whatever the rest of the
//! equation leaves, one side of a two-column equation is a function of the other --
//! and an eliminated column has no value of its own to be fixed at. On
//! `neos-2987310-joes` HiGHS removes 21762 of 27837 columns that way and solves the
//! root in half a second; we kept all of them and did not finish the root LP in 60.
//!
//! So this pass renumbers, and [`Postsolve`] is the record of what it did. Each
//! reduction pushes one [`Record`] that says how to compute the column it removed from
//! columns still present at that moment. A column present then is either kept or
//! removed by a *later* record, so replaying the records backwards always finds what it
//! reads already known, and no ordering beyond the stack is needed.
//!
//! The map is the whole risk -- one that renumbers wrongly returns a plausible vector
//! for the wrong columns, which still scores an objective -- so, as in `compact` and
//! `aggregate`, the tests solve both models and check the expanded point on the
//! original rather than reading the code and believing it.
//!
//! Only primal values come back. A MIP solver reports no duals, so the reductions here
//! may also drop points that are feasible but never optimal (dual fixing, and turning
//! an inequality into the equation an optimum must satisfy), which is sound only when
//! the rows in the model are all the rows there are. A search with a separator has rows
//! it does not show, and does not run this.

use crate::compact::Infeasible;
use crate::model::{Problem, VarType};
use crate::presolve::Link;
use crate::sparse::SparseMatrix;

/// Below this a coefficient is a cancellation, not an entry.
const ZERO: f64 = 1e-12;

/// The widest ratio of coefficients an elimination may divide by.
///
/// Substituting `x_p = (t - b x_q) / a` multiplies every coefficient `x_p` had by
/// `b / a`, so a small pivot scales the rows it lands in by its inverse. HiGHS refuses
/// a pivot below a hundredth of the row's largest entry for the same reason; this is a
/// little looser because a doubleton has only the one other entry to compare with.
const PIVOT_RATIO: f64 = 1e3;

/// How many passes over the model before giving up on a fixpoint.
const MAX_SWEEPS: usize = 50;

/// Passes of dual bound propagation in [`Work::sweep_dominated`]. Like the primal
/// kind, two constraints can narrow a bound geometrically without end.
const DUAL_PASSES: usize = 10;

/// The longest equation, and the longest column, [`Work::sweep_substitute`] tries.
const SUBSTITUTE_LENGTH: usize = 20;

/// How many more entries one substitution may leave the model with than it had.
const SUBSTITUTE_FILL: usize = 10;

/// The smallest pivot a substitution divides by, as a fraction of the equation's
/// largest entry: HiGHS's choice, for the reason [`PIVOT_RATIO`] gives.
const SUBSTITUTE_PIVOT: f64 = 100.0;

/// Rows longer than this are not searched for a column's dominators.
const DOMINANCE_ROW: usize = 1000;

/// Coefficient comparisons one [`Work::sweep_dominance`] may make.
const DOMINANCE_BUDGET: usize = 50_000_000;

/// How far over the dual box a reduced cost must keep its sign, relative to the cost,
/// before the column it belongs to is fixed. Zero would trust a sign that rounding
/// gave, and a column whose reduced cost can reach zero may be needed off its bound.
const DUAL_TOLERANCE: f64 = 1e-7;

/// How to compute one removed column.
#[derive(Clone, Debug)]
enum Record {
    /// It was fixed at this value.
    Fixed { column: usize, value: f64 },
    /// It was solved for from a row: `lb <= coefficient * x + sum(a x_k) <= ub`, with
    /// `x` within `[col_lb, col_ub]`. For an equation that is one value; for an
    /// inequality whose other side was given away it is an interval, and any point of
    /// it will do, since the column has no cost and no other row.
    Row {
        column: usize,
        coefficient: f64,
        terms: Vec<(usize, f64)>,
        lb: f64,
        ub: f64,
        col_lb: f64,
        col_ub: f64,
        integer: bool,
    },
}

/// How to read a solution of a reduced model as a solution of the model it came from.
#[derive(Clone, Debug)]
pub struct Postsolve {
    /// Original column index of each column of the reduced model.
    columns: Vec<usize>,
    /// What was removed, in the order it was removed.
    records: Vec<Record>,
    n_original: usize,
}

impl Postsolve {
    /// A solution of the reduced model, read as a solution of the original.
    pub fn expand(&self, x: &[f64]) -> Vec<f64> {
        debug_assert_eq!(x.len(), self.columns.len());
        let mut full = vec![0.0; self.n_original];
        for (&j, &value) in self.columns.iter().zip(x) {
            full[j] = value;
        }
        // Backwards: see the module notes on why this order is the right one.
        for record in self.records.iter().rev() {
            match *record {
                Record::Fixed { column, value } => full[column] = value,
                Record::Row {
                    column,
                    coefficient,
                    ref terms,
                    lb,
                    ub,
                    col_lb,
                    col_ub,
                    integer,
                } => {
                    let rest: f64 = terms.iter().map(|&(k, a)| a * full[k]).sum();
                    let (from, to) = ((lb - rest) / coefficient, (ub - rest) / coefficient);
                    let (low, high) = if coefficient > 0.0 { (from, to) } else { (to, from) };
                    let (low, high) = (low.max(col_lb), high.min(col_ub));
                    // Any point of the interval; zero when it is there, since that is
                    // what a slack with nothing to say would sit at, and otherwise the
                    // end nearer it. An equation gives `low == high` up to rounding.
                    let mut value = if low > high {
                        0.5 * (low + high)
                    } else if low > 0.0 {
                        low
                    } else if high < 0.0 {
                        high
                    } else {
                        0.0
                    };
                    if integer {
                        value = value.round();
                    }
                    full[column] = value;
                }
            }
        }
        full
    }

    /// A solution of the original, read as one of the reduced model.
    ///
    /// Used for a caller-supplied starting point. It is only a hint: a column fixed by
    /// its cost may be fixed somewhere the starting point was not.
    pub fn restrict(&self, x: &[f64]) -> Vec<f64> {
        self.columns.iter().map(|&j| x[j]).collect()
    }

    /// Columns of the original that survived, in the reduced model's order.
    pub fn columns(&self) -> &[usize] {
        &self.columns
    }

    /// How many columns were removed.
    pub fn eliminated(&self) -> usize {
        self.records.len()
    }

    /// This map followed by `next`, which reduced the model this one produced.
    ///
    /// `next` names columns of the intermediate model, which are this map's kept
    /// columns by position, so renaming them through `columns` makes it a map from the
    /// original. Its records come after these, so replay undoes it first.
    pub fn then(mut self, next: Postsolve) -> Postsolve {
        debug_assert_eq!(next.n_original, self.columns.len());
        let name = |k: usize| self.columns[k];
        for record in next.records {
            self.records.push(match record {
                Record::Fixed { column, value } => Record::Fixed { column: name(column), value },
                Record::Row { column, coefficient, terms, lb, ub, col_lb, col_ub, integer } => Record::Row {
                    column: name(column),
                    coefficient,
                    terms: terms.into_iter().map(|(k, a)| (name(k), a)).collect(),
                    lb,
                    ub,
                    col_lb,
                    col_ub,
                    integer,
                },
            });
        }
        self.columns = next.columns.iter().map(|&k| self.columns[k]).collect();
        self
    }
}

/// What each rule removed, for reporting.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Stats {
    pub fixed_columns: usize,
    pub row_singletons: usize,
    pub dual_fixed: usize,
    pub column_singletons: usize,
    pub doubletons: usize,
    pub dominated: usize,
    pub dominated_pairs: usize,
    pub switched: usize,
    pub substituted: usize,
    pub linked: usize,
    pub removed_rows: usize,
}

impl Stats {
    fn absorb(&mut self, other: &Stats) {
        self.fixed_columns += other.fixed_columns;
        self.row_singletons += other.row_singletons;
        self.dual_fixed += other.dual_fixed;
        self.column_singletons += other.column_singletons;
        self.doubletons += other.doubletons;
        self.dominated += other.dominated;
        self.dominated_pairs += other.dominated_pairs;
        self.switched += other.switched;
        self.substituted += other.substituted;
        self.linked += other.linked;
        self.removed_rows += other.removed_rows;
    }
}

/// The model while it is being reduced: both orientations, kept in step.
struct Work<'a> {
    problem: &'a Problem,
    rows: Vec<Vec<(usize, f64)>>,
    cols: Vec<Vec<(usize, f64)>>,
    row_lb: Vec<f64>,
    row_ub: Vec<f64>,
    col_lb: Vec<f64>,
    col_ub: Vec<f64>,
    obj: Vec<f64>,
    offset: f64,
    row_alive: Vec<bool>,
    col_alive: Vec<bool>,
    records: Vec<Record>,
    tolerance: f64,
    stats: Stats,
}

/// Lower `bound` to `to` if that is a real step down; true if it moved.
fn lower_to(bound: &mut f64, to: f64) -> bool {
    let step = 1e-6 * to.abs().max(1.0);
    if to < *bound - step {
        *bound = to;
        return true;
    }
    false
}

/// Raise `bound` to `to` if that is a real step up; true if it moved.
fn raise_to(bound: &mut f64, to: f64) -> bool {
    let step = 1e-6 * to.abs().max(1.0);
    if to > *bound + step {
        *bound = to;
        return true;
    }
    false
}

fn take(list: &mut Vec<(usize, f64)>, key: usize) -> Option<f64> {
    let at = list.iter().position(|&(k, _)| k == key)?;
    Some(list.swap_remove(at).1)
}

impl<'a> Work<'a> {
    fn new(problem: &'a Problem, tolerance: f64) -> Self {
        let (n, m) = (problem.n_cols(), problem.n_rows());
        let mut rows = vec![Vec::new(); m];
        let mut cols = vec![Vec::new(); n];
        for (j, col) in cols.iter_mut().enumerate() {
            let (indices, values) = problem.matrix.column(j);
            for (&i, &a) in indices.iter().zip(values) {
                if a != 0.0 {
                    col.push((i, a));
                    rows[i].push((j, a));
                }
            }
        }
        Work {
            problem,
            rows,
            cols,
            row_lb: problem.row_lb.clone(),
            row_ub: problem.row_ub.clone(),
            col_lb: problem.col_lb.clone(),
            col_ub: problem.col_ub.clone(),
            obj: problem.obj.clone(),
            offset: problem.obj_offset,
            row_alive: vec![true; m],
            col_alive: vec![true; n],
            records: Vec::new(),
            tolerance,
            stats: Stats::default(),
        }
    }

    fn integer(&self, j: usize) -> bool {
        self.problem.col_type[j] == VarType::Integer
    }

    fn integral(&self, v: f64) -> bool {
        (v - v.round()).abs() <= self.tolerance
    }

    /// Move a constant `delta` of row `i`'s activity into its bounds.
    fn shift(&mut self, i: usize, delta: f64) {
        if self.row_lb[i].is_finite() {
            self.row_lb[i] -= delta;
        }
        if self.row_ub[i].is_finite() {
            self.row_ub[i] -= delta;
        }
    }

    /// Add `delta` to the coefficient of column `j` in row `i`, both orientations.
    fn add(&mut self, i: usize, j: usize, delta: f64) {
        let row = &mut self.rows[i];
        match row.iter().position(|&(k, _)| k == j) {
            Some(at) => {
                let value = row[at].1 + delta;
                if value.abs() <= ZERO {
                    row.swap_remove(at);
                    take(&mut self.cols[j], i);
                } else {
                    row[at].1 = value;
                    let col = &mut self.cols[j];
                    let at = col.iter().position(|&(k, _)| k == i).unwrap();
                    col[at].1 = value;
                }
            }
            None if delta.abs() > ZERO => {
                row.push((j, delta));
                self.cols[j].push((i, delta));
            }
            None => {}
        }
    }

    fn kill_row(&mut self, i: usize) {
        for (j, _) in std::mem::take(&mut self.rows[i]) {
            take(&mut self.cols[j], i);
        }
        self.row_alive[i] = false;
        self.stats.removed_rows += 1;
    }

    fn kill_col(&mut self, j: usize) {
        for (i, _) in std::mem::take(&mut self.cols[j]) {
            take(&mut self.rows[i], j);
        }
        self.col_alive[j] = false;
    }

    fn fix(&mut self, j: usize, value: f64) {
        self.offset += self.obj[j] * value;
        for k in 0..self.cols[j].len() {
            let (i, a) = self.cols[j][k];
            self.shift(i, a * value);
        }
        self.kill_col(j);
        self.records.push(Record::Fixed { column: j, value });
    }

    /// Narrow column `j` to `[lb, ub]`, rounding inwards for an integer column.
    fn tighten(&mut self, j: usize, mut lb: f64, mut ub: f64) -> Result<(), Infeasible> {
        if self.integer(j) {
            lb = (lb - self.tolerance).ceil();
            ub = (ub + self.tolerance).floor();
        }
        let lb = self.col_lb[j].max(lb);
        let ub = self.col_ub[j].min(ub);
        if lb > ub + self.tolerance {
            return Err(Infeasible);
        }
        // Inside the tolerance the two have met; keep one point rather than an
        // inverted interval nothing downstream knows how to read.
        let ub = ub.max(lb);
        self.col_lb[j] = lb;
        self.col_ub[j] = ub;
        Ok(())
    }

    /// The range of row `i`'s activity without column `skip`, over the bounds.
    fn activity_without(&self, i: usize, skip: usize) -> (f64, f64) {
        let (mut low, mut high) = (0.0, 0.0);
        for &(k, a) in &self.rows[i] {
            if k == skip {
                continue;
            }
            let (l, u) = (self.col_lb[k], self.col_ub[k]);
            if a > 0.0 {
                low += a * l;
                high += a * u;
            } else {
                low += a * u;
                high += a * l;
            }
        }
        (low, high)
    }

    /// Rows that are empty, free or a single entry. False when nothing changed.
    fn sweep_rows(&mut self) -> Result<bool, Infeasible> {
        let mut changed = false;
        for i in 0..self.rows.len() {
            if !self.row_alive[i] {
                continue;
            }
            let (lo, hi) = (self.row_lb[i], self.row_ub[i]);
            if !lo.is_finite() && !hi.is_finite() {
                self.kill_row(i);
                changed = true;
                continue;
            }
            match self.rows[i].len() {
                0 => {
                    if lo > self.tolerance || hi < -self.tolerance {
                        return Err(Infeasible);
                    }
                    self.kill_row(i);
                    changed = true;
                }
                1 => {
                    let (j, a) = self.rows[i][0];
                    let (from, to) = (lo / a, hi / a);
                    let (l, u) = if a > 0.0 { (from, to) } else { (to, from) };
                    self.tighten(j, l, u)?;
                    self.kill_row(i);
                    self.stats.row_singletons += 1;
                    changed = true;
                }
                _ => {}
            }
        }
        Ok(changed)
    }

    /// Columns that are fixed, empty, or that no row stops from moving to the end
    /// their cost prefers.
    fn sweep_columns(&mut self) -> Result<bool, Infeasible> {
        let mut changed = false;
        for j in 0..self.cols.len() {
            if !self.col_alive[j] {
                continue;
            }
            let (l, u) = (self.col_lb[j], self.col_ub[j]);
            if u - l <= self.tolerance {
                let value = if self.integer(j) { l.round() } else { l };
                if self.integer(j) && !self.integral(l) {
                    return Err(Infeasible);
                }
                self.fix(j, value);
                self.stats.fixed_columns += 1;
                changed = true;
                continue;
            }
            // Dual fixing: if lowering the column never takes any row out of its
            // bounds and the cost does not reward raising it, some optimum has it at
            // its lower bound, and the same the other way round. An empty column is
            // the case with no rows to ask.
            let c = self.obj[j];
            let down = self.cols[j].iter().all(|&(i, a)| {
                if a > 0.0 { self.row_lb[i] == f64::NEG_INFINITY } else { self.row_ub[i] == f64::INFINITY }
            });
            let up = self.cols[j].iter().all(|&(i, a)| {
                if a > 0.0 { self.row_ub[i] == f64::INFINITY } else { self.row_lb[i] == f64::NEG_INFINITY }
            });
            let value = if c >= 0.0 && down && l.is_finite() {
                Some(l)
            } else if c <= 0.0 && up && u.is_finite() {
                Some(u)
            } else if c == 0.0 && down && up {
                // Unconstrained and free of charge: anywhere will do.
                Some(0.0f64.clamp(l, u))
            } else {
                None
            };
            if let Some(value) = value {
                self.fix(j, value);
                self.stats.dual_fixed += 1;
                changed = true;
            }
        }
        Ok(changed)
    }

    /// A continuous column in exactly one row is whatever that row leaves for it.
    ///
    /// Three cases. In an equation it is determined, and the row becomes a range on
    /// the others that keeps it inside its own bounds; its cost moves onto them. In an
    /// inequality with no cost it is a slack, and the row loses as much as the slack
    /// could absorb. In an inequality with a cost it is pushed against the side its
    /// cost prefers, and when its bounds can never stop it getting there (it is
    /// implied free) that side holds at every optimum, so the row is that equation and
    /// the first case applies.
    fn sweep_singletons(&mut self) -> bool {
        let mut changed = false;
        for j in 0..self.cols.len() {
            if !self.col_alive[j] || self.integer(j) || self.cols[j].len() != 1 {
                continue;
            }
            let (i, a) = self.cols[j][0];
            if self.rows[i].len() < 2 {
                continue;
            }
            let largest = self.rows[i].iter().fold(0.0f64, |m, &(_, v)| m.max(v.abs()));
            if a.abs() * PIVOT_RATIO < largest {
                continue;
            }
            let c = self.obj[j];
            let (l, u) = (self.col_lb[j], self.col_ub[j]);
            let equation = self.row_ub[i] - self.row_lb[i] <= self.tolerance;
            if !equation && c != 0.0 {
                let target = if c / a > 0.0 { self.row_lb[i] } else { self.row_ub[i] };
                if !target.is_finite() {
                    continue;
                }
                let (low, high) = self.activity_without(i, j);
                let (from, to) = ((target - high) / a, (target - low) / a);
                let (need_l, need_u) = if a > 0.0 { (from, to) } else { (to, from) };
                if need_l < l - self.tolerance || need_u > u + self.tolerance {
                    continue;
                }
                self.row_lb[i] = target;
                self.row_ub[i] = target;
            }
            let (lo, hi) = (self.row_lb[i], self.row_ub[i]);
            let terms: Vec<(usize, f64)> =
                self.rows[i].iter().copied().filter(|&(k, _)| k != j).collect();
            if c != 0.0 {
                // Equation `a x + r = t`: c x = (c / a)(t - r).
                let ratio = c / a;
                self.offset += ratio * lo;
                for &(k, ak) in &terms {
                    self.obj[k] -= ratio * ak;
                }
            }
            self.records.push(Record::Row {
                column: j,
                coefficient: a,
                terms,
                lb: lo,
                ub: hi,
                col_lb: l,
                col_ub: u,
                integer: false,
            });
            // What `a x` can take up, and so how far the rest may move.
            let (ends_l, ends_u) = (a * l, a * u);
            let (least, most) = if a > 0.0 { (ends_l, ends_u) } else { (ends_u, ends_l) };
            self.row_lb[i] = if lo.is_finite() && most.is_finite() { lo - most } else { f64::NEG_INFINITY };
            self.row_ub[i] = if hi.is_finite() && least.is_finite() { hi - least } else { f64::INFINITY };
            self.kill_col(j);
            self.stats.column_singletons += 1;
            changed = true;
        }
        changed
    }

    /// An equation on two columns makes one an affine function of the other.
    ///
    /// The one removed is a continuous one when there is a choice, since a continuous
    /// column can be any affine function of anything. An integer column can be removed
    /// only into another integer column and only by integral multiples and offsets, or
    /// it would come back fractional at a point the search thought was integral.
    fn sweep_doubletons(&mut self) -> Result<bool, Infeasible> {
        let mut changed = false;
        for i in 0..self.rows.len() {
            if !self.row_alive[i] || self.rows[i].len() != 2 {
                continue;
            }
            let (lo, hi) = (self.row_lb[i], self.row_ub[i]);
            if !lo.is_finite() || hi - lo > self.tolerance {
                continue;
            }
            let t = 0.5 * (lo + hi);
            let (p, a) = self.rows[i][0];
            let (q, b) = self.rows[i][1];
            // The order to try them in: continuous first, then the shorter column.
            let key = |j: usize| (self.integer(j), self.cols[j].len(), j);
            let options = if key(p) <= key(q) { [(p, a, q, b), (q, b, p, a)] } else { [(q, b, p, a), (p, a, q, b)] };
            for (drop, ad, keep, ak) in options {
                // x_drop = alpha + beta x_keep
                let alpha = t / ad;
                let beta = -ak / ad;
                if beta.abs() > PIVOT_RATIO || beta.abs() * PIVOT_RATIO < 1.0 {
                    continue;
                }
                if self.integer(drop)
                    && !(self.integer(keep) && self.integral(alpha) && self.integral(beta))
                {
                    continue;
                }
                let (alpha, beta) = if self.integer(drop) { (alpha.round(), beta.round()) } else { (alpha, beta) };
                let (from, to) = ((self.col_lb[drop] - alpha) / beta, (self.col_ub[drop] - alpha) / beta);
                let (l, u) = if beta > 0.0 { (from, to) } else { (to, from) };
                self.records.push(Record::Row {
                    column: drop,
                    coefficient: ad,
                    terms: vec![(keep, ak)],
                    lb: t,
                    ub: t,
                    col_lb: self.col_lb[drop],
                    col_ub: self.col_ub[drop],
                    integer: self.integer(drop),
                });
                self.tighten(keep, l, u)?;
                self.kill_row(i);
                for k in 0..self.cols[drop].len() {
                    let (r, ar) = self.cols[drop][k];
                    self.shift(r, ar * alpha);
                }
                let entries = self.cols[drop].clone();
                for (r, ar) in entries {
                    self.add(r, keep, ar * beta);
                }
                self.obj[keep] += self.obj[drop] * beta;
                self.offset += self.obj[drop] * alpha;
                self.kill_col(drop);
                self.stats.doubletons += 1;
                changed = true;
                break;
            }
        }
        Ok(changed)
    }

    /// Are column `j`'s lower and upper bounds infinite, or implied by one of its rows
    /// from the bounds of the others?
    fn implied_free(&self, j: usize) -> (bool, bool) {
        let (l, u) = (self.col_lb[j], self.col_ub[j]);
        let (mut lower, mut upper) = (!l.is_finite(), !u.is_finite());
        for &(i, a) in &self.cols[j] {
            if lower && upper {
                break;
            }
            let (low, high) = self.activity_without(i, j);
            // a x in [lb - high, ub - low].
            let (from, to) = ((self.row_lb[i] - high) / a, (self.row_ub[i] - low) / a);
            let (least, most) = if a > 0.0 { (from, to) } else { (to, from) };
            lower |= least.is_finite() && least >= l - self.tolerance;
            upper |= most.is_finite() && most <= u + self.tolerance;
        }
        (lower, upper)
    }

    /// A column whose bounds hold anyway is whatever an equation it is in leaves.
    ///
    /// The doubleton rule's elimination for longer equations: `a x_j + sum a_k x_k = t`
    /// gives `x_j = (t - sum a_k x_k) / a`, which goes into every other row of `x_j`
    /// and into the objective, and the equation goes. That loses `x_j`'s bounds, which
    /// is why they have to be implied (by some row, from the bounds of the rest) or
    /// absent, and fills in the rows `x_j` was in, which is why the rows and columns
    /// it is tried on are short and the fill it may cause is limited. HiGHS calls this
    /// the aggregator; on `app1-1` it is most of the difference between its 1059
    /// columns and the 2477 we kept.
    ///
    /// An integer column goes only when the rest of the equation makes it integral:
    /// integer columns with coefficients and a right-hand side that are integral
    /// multiples of its own.
    fn sweep_substitute(&mut self) -> Result<bool, Infeasible> {
        let mut changed = false;
        for r in 0..self.rows.len() {
            let len_r = self.rows[r].len();
            if !self.row_alive[r] || len_r < 3 || len_r > SUBSTITUTE_LENGTH {
                continue;
            }
            let (lo, hi) = (self.row_lb[r], self.row_ub[r]);
            if !lo.is_finite() || hi - lo > self.tolerance {
                continue;
            }
            let t = 0.5 * (lo + hi);
            let largest = self.rows[r].iter().fold(0.0f64, |m, &(_, v)| m.max(v.abs()));
            let mut order: Vec<(usize, f64)> = self.rows[r].clone();
            order.sort_by_key(|&(j, _)| (self.integer(j), self.cols[j].len(), j));
            let pick = order.into_iter().find(|&(j, a)| {
                let len_j = self.cols[j].len();
                if len_j > SUBSTITUTE_LENGTH || a.abs() * SUBSTITUTE_PIVOT < largest {
                    return false;
                }
                let added = (len_r - 1) * (len_j - 1);
                if added > len_r + len_j - 1 + SUBSTITUTE_FILL {
                    return false;
                }
                if self.integer(j) {
                    let whole = self.integral(t / a)
                        && self.rows[r]
                            .iter()
                            .all(|&(k, ak)| k == j || (self.integer(k) && self.integral(ak / a)));
                    if !whole {
                        return false;
                    }
                }
                self.implied_free(j) == (true, true)
            });
            let Some((j, a)) = pick else {
                continue;
            };
            let terms: Vec<(usize, f64)> = self.rows[r].iter().copied().filter(|&(k, _)| k != j).collect();
            self.records.push(Record::Row {
                column: j,
                coefficient: a,
                terms: terms.clone(),
                lb: t,
                ub: t,
                col_lb: self.col_lb[j],
                col_ub: self.col_ub[j],
                integer: self.integer(j),
            });
            self.kill_row(r);
            // x_j = t / a - sum (a_k / a) x_k.
            let entries = self.cols[j].clone();
            for (s, aj) in entries {
                self.shift(s, aj * t / a);
                for &(k, ak) in &terms {
                    self.add(s, k, -aj * ak / a);
                }
            }
            let c = self.obj[j];
            if c != 0.0 {
                self.offset += c * t / a;
                for &(k, ak) in &terms {
                    self.obj[k] -= c * ak / a;
                }
            }
            self.kill_col(j);
            self.stats.substituted += 1;
            changed = true;
        }
        Ok(changed)
    }

    /// Continuous columns whose reduced cost has one sign at every dual solution.
    ///
    /// With the integer columns at their values in an optimum, what is left is an LP,
    /// and the box below holds every dual solution of it. A row only ever bounded
    /// below has a dual of at least zero, one bounded above at most zero. A continuous
    /// column with no upper bound needs a reduced cost `c - a y` of at least zero, and
    /// with no lower bound at most zero, which bounds the duals further, the same way
    /// a row bounds the columns in it. A continuous column whose reduced cost is then
    /// positive over the whole box sits at its lower bound in every optimum of that LP,
    /// by complementary slackness, and so in every optimum of the model. One negative
    /// over the box sits at its upper bound.
    ///
    /// A bound some row already implies can be left out of the model without changing
    /// it, and is, for this: on `neos-2987310-joes` almost no column has an infinite
    /// bound of its own. But two bounds can each imply the other, and leaving both out
    /// changes the model, so a row a bound is taken from keeps every other bound in it.
    ///
    /// Integer columns take no part as columns. An integer column cannot move by the
    /// small amount the argument rests on, so its reduced cost proves nothing about it,
    /// and it states no dual constraint, being fixed in the LP the argument is about.
    fn sweep_dominated(&mut self) -> Result<bool, Infeasible> {
        let (m, n) = (self.rows.len(), self.cols.len());
        // Row activity ranges over the bounds, with the infinite contributions
        // counted apart so that one entry's share can be taken out exactly.
        let mut low = vec![(0.0f64, 0usize); m];
        let mut high = vec![(0.0f64, 0usize); m];
        for i in (0..m).filter(|&i| self.row_alive[i]) {
            for &(k, a) in &self.rows[i] {
                let (l, u) = if a > 0.0 { (self.col_lb[k], self.col_ub[k]) } else { (self.col_ub[k], self.col_lb[k]) };
                if l.is_finite() { low[i].0 += a * l } else { low[i].1 += 1 }
                if u.is_finite() { high[i].0 += a * u } else { high[i].1 += 1 }
            }
        }
        let without = |(sum, infinite): (f64, usize), share: f64| match (infinite, share.is_finite()) {
            (0, _) => sum - share,
            (1, false) => sum,
            _ => f64::NAN,
        };
        // Which bounds of the continuous columns can be left out.
        let mut free_lb = vec![false; n];
        let mut free_ub = vec![false; n];
        let mut locked = vec![false; n];
        // The column whose bound a row's columns have lost, or `MANY` past one; a row
        // can still imply a bound of that one column, from the others' bounds.
        const NONE: usize = usize::MAX;
        const MANY: usize = usize::MAX - 1;
        let mut lost = vec![NONE; m];
        let mut row_locked = vec![false; m];
        for j in 0..n {
            if !self.col_alive[j] || self.integer(j) {
                continue;
            }
            free_lb[j] = !self.col_lb[j].is_finite();
            free_ub[j] = !self.col_ub[j].is_finite();
            for upper in [false, true] {
                if (if upper { free_ub[j] } else { free_lb[j] }) || locked[j] {
                    continue;
                }
                let bound = if upper { self.col_ub[j] } else { self.col_lb[j] };
                let source = self.cols[j].iter().copied().find(|&(i, a)| {
                    if lost[i] != NONE && lost[i] != j {
                        return false;
                    }
                    // a x <= ub - (least of the rest), a x >= lb - (most of the rest).
                    let (lo_a, hi_a) = if a > 0.0 { (self.col_lb[j] * a, self.col_ub[j] * a) } else { (self.col_ub[j] * a, self.col_lb[j] * a) };
                    let want_hi = (a > 0.0) == upper;
                    let implied = if want_hi {
                        let rest = without(low[i], lo_a);
                        (self.row_ub[i] - rest) / a
                    } else {
                        let rest = without(high[i], hi_a);
                        (self.row_lb[i] - rest) / a
                    };
                    let implied = if implied.is_nan() { return false } else { implied };
                    if upper { implied <= bound + self.tolerance } else { implied >= bound - self.tolerance }
                });
                if let Some((i, _)) = source {
                    if upper { free_ub[j] = true } else { free_lb[j] = true }
                    for &(r, _) in &self.cols[j] {
                        lost[r] = if lost[r] == NONE || lost[r] == j { j } else { MANY };
                    }
                    if !row_locked[i] {
                        row_locked[i] = true;
                        for &(k, _) in &self.rows[i] {
                            locked[k] = true;
                        }
                    }
                    // Its own bound stays out; only the others' must stay in.
                    locked[j] = false;
                }
            }
        }
        // The box on the duals: signs from the rows, then the columns' constraints.
        let mut y_lb = vec![f64::NEG_INFINITY; m];
        let mut y_ub = vec![f64::INFINITY; m];
        for i in (0..m).filter(|&i| self.row_alive[i]) {
            match (self.row_lb[i].is_finite(), self.row_ub[i].is_finite()) {
                (true, false) => y_lb[i] = 0.0,
                (false, true) => y_ub[i] = 0.0,
                (false, false) => (y_lb[i], y_ub[i]) = (0.0, 0.0),
                (true, true) => {}
            }
        }
        // The dual activity `a y` of column j over the box, as (least, most), each a
        // finite sum and a count of infinite terms.
        let range = |y_lb: &[f64], y_ub: &[f64], j: usize| {
            let (mut least, mut most) = ((0.0, 0), (0.0, 0));
            for &(i, a) in &self.cols[j] {
                let (l, u) = if a > 0.0 { (a * y_lb[i], a * y_ub[i]) } else { (a * y_ub[i], a * y_lb[i]) };
                if l.is_finite() { least.0 += l } else { least.1 += 1 }
                if u.is_finite() { most.0 += u } else { most.1 += 1 }
            }
            (least, most)
        };
        for _ in 0..DUAL_PASSES {
            let mut moved = false;
            for j in 0..n {
                if !self.col_alive[j] || self.integer(j) || !(free_lb[j] || free_ub[j]) {
                    continue;
                }
                let c = self.obj[j];
                // Kept up to date as the loop moves the duals of j's rows, rather than
                // summed afresh for every entry, which is quadratic in a long column.
                let (mut all_least, mut all_most) = range(&y_lb, &y_ub, j);
                let share = |y_lb: &[f64], y_ub: &[f64], i: usize, a: f64| {
                    if a > 0.0 { (a * y_lb[i], a * y_ub[i]) } else { (a * y_ub[i], a * y_lb[i]) }
                };
                let swap = |sum: &mut (f64, usize), old: f64, new: f64| {
                    if old.is_finite() { sum.0 -= old } else { sum.1 -= 1 }
                    if new.is_finite() { sum.0 += new } else { sum.1 += 1 }
                };
                for k in 0..self.cols[j].len() {
                    let (i, a) = self.cols[j][k];
                    let (l, u) = share(&y_lb, &y_ub, i, a);
                    let (least, most) = (without(all_least, l), without(all_most, u));
                    // No upper bound: a y <= c, so a y_i <= c - least of the rest.
                    if free_ub[j] && least.is_finite() {
                        let limit = (c - least) / a;
                        moved |= if a > 0.0 { lower_to(&mut y_ub[i], limit) } else { raise_to(&mut y_lb[i], limit) };
                    }
                    // No lower bound: a y >= c.
                    if free_lb[j] && most.is_finite() {
                        let limit = (c - most) / a;
                        moved |= if a > 0.0 { raise_to(&mut y_lb[i], limit) } else { lower_to(&mut y_ub[i], limit) };
                    }
                    let (new_l, new_u) = share(&y_lb, &y_ub, i, a);
                    swap(&mut all_least, l, new_l);
                    swap(&mut all_most, u, new_u);
                    if y_lb[i] > y_ub[i] + DUAL_TOLERANCE {
                        // No dual solution: the LP is infeasible or unbounded, and
                        // which is not this rule's to say.
                        return Ok(false);
                    }
                }
            }
            if !moved {
                break;
            }
        }
        // Decided together on the one box and applied after: every optimum has all of
        // them at once.
        let mut fixes = Vec::new();
        for j in 0..n {
            if !self.col_alive[j] || self.integer(j) || self.cols[j].is_empty() {
                continue;
            }
            let c = self.obj[j];
            let (least, most) = range(&y_lb, &y_ub, j);
            let total = |(sum, infinite): (f64, usize), sign: f64| if infinite > 0 { sign * f64::INFINITY } else { sum };
            let (least, most) = (total(least, -1.0), total(most, 1.0));
            let margin = DUAL_TOLERANCE * c.abs().max(1.0);
            let value = if c - most > margin && self.col_lb[j].is_finite() {
                self.col_lb[j]
            } else if c - least < -margin && self.col_ub[j].is_finite() {
                self.col_ub[j]
            } else {
                continue;
            };
            fixes.push((j, value));
        }
        for &(j, value) in &fixes {
            self.fix(j, value);
        }
        self.stats.dominated += fixes.len();
        Ok(!fixes.is_empty())
    }

    /// A column another does better than in every row, at no more cost.
    ///
    /// Column `j` dominates `k` when `c_j <= c_k` and moving an amount from `k` onto
    /// `j` keeps every row as satisfied as it was: `a_ij >= a_ik` in a row bounded only
    /// below, `<=` in one bounded only above, equal in one bounded on both sides. Then
    /// in any solution with `k` above its lower bound the excess can go onto `j`, as
    /// long as `j` has room, and the cost does not rise. With no upper bound on `j`, or
    /// one some row implies anyway, it always has room, and `k` can be fixed at its
    /// lower bound. With no lower bound on `k`, the same move the other way fixes `j`
    /// at its upper bound. The argument moves values, not duals, so it holds for
    /// integer columns too, provided an integral amount lands on an integer column:
    /// a continuous `k` goes only onto a continuous `j`.
    ///
    /// `neos-2987310-joes` has its continuous columns in identical pairs, each in a
    /// row `sum x <= 1` that makes its own bound of 1 redundant, so one of each pair
    /// goes: this is most of the 12420 columns HiGHS's base rules remove from it.
    ///
    /// Candidates for dominating `k` are the columns of its shortest row, which keeps
    /// the pairs compared to those that share a row.
    fn sweep_dominance(&mut self) -> bool {
        let (m, n) = (self.rows.len(), self.cols.len());
        let mut low = vec![(0.0f64, 0usize); m];
        let mut high = vec![(0.0f64, 0usize); m];
        for i in (0..m).filter(|&i| self.row_alive[i]) {
            for &(k, a) in &self.rows[i] {
                let (l, u) = if a > 0.0 { (self.col_lb[k], self.col_ub[k]) } else { (self.col_ub[k], self.col_lb[k]) };
                if l.is_finite() { low[i].0 += a * l } else { low[i].1 += 1 }
                if u.is_finite() { high[i].0 += a * u } else { high[i].1 += 1 }
            }
        }
        let without = |(sum, infinite): (f64, usize), share: f64| match (infinite, share.is_finite()) {
            (0, _) => sum - share,
            (1, false) => sum,
            _ => f64::NAN,
        };
        // Is column j's upper (lower) bound infinite or implied by one of its rows,
        // from the bounds of the others? Fixing columns at values inside their bounds
        // keeps an implication true, so these hold for the whole sweep.
        let implied = |j: usize, upper: bool| {
            let bound = if upper { self.col_ub[j] } else { self.col_lb[j] };
            if !bound.is_finite() {
                return true;
            }
            self.cols[j].iter().any(|&(i, a)| {
                let (lo_a, hi_a) = if a > 0.0 { (self.col_lb[j] * a, self.col_ub[j] * a) } else { (self.col_ub[j] * a, self.col_lb[j] * a) };
                let value = if (a > 0.0) == upper {
                    (self.row_ub[i] - without(low[i], lo_a)) / a
                } else {
                    (self.row_lb[i] - without(high[i], hi_a)) / a
                };
                if upper { value <= bound + self.tolerance } else { value >= bound - self.tolerance }
            })
        };
        let roomy_up: Vec<bool> = (0..n).map(|j| self.col_alive[j] && implied(j, true)).collect();
        let roomy_down: Vec<bool> = (0..n).map(|j| self.col_alive[j] && implied(j, false)).collect();
        for col in &mut self.cols {
            col.sort_unstable_by_key(|&(i, _)| i);
        }
        // How a row lets its coefficients differ: +1 bounded below only (more is
        // better), -1 above only, 0 both (they must match), None free.
        // Fixing a column moves a finite bound and never makes one infinite, so these
        // stay right while the sweep fixes columns.
        let kinds: Vec<Option<i8>> = (0..m)
            .map(|i| match (self.row_lb[i].is_finite(), self.row_ub[i].is_finite()) {
                (true, false) => Some(1),
                (false, true) => Some(-1),
                (true, true) => Some(0),
                (false, false) => None,
            })
            .collect();
        let dominates = |cols: &[Vec<(usize, f64)>], j: usize, k: usize| {
            let (cj, ck) = (&cols[j], &cols[k]);
            let (mut p, mut q) = (0, 0);
            while p < cj.len() || q < ck.len() {
                let (i, aj, ak) = match (cj.get(p), ck.get(q)) {
                    (Some(&(ri, a)), Some(&(rk, b))) if ri == rk => { p += 1; q += 1; (ri, a, b) }
                    (Some(&(ri, a)), Some(&(rk, _))) if ri < rk => { p += 1; (ri, a, 0.0) }
                    (Some(&(ri, a)), None) => { p += 1; (ri, a, 0.0) }
                    (_, Some(&(rk, b))) => { q += 1; (rk, 0.0, b) }
                    (None, None) => unreachable!(),
                };
                let fine = match kinds[i] {
                    Some(1) => aj >= ak,
                    Some(-1) => aj <= ak,
                    Some(_) => aj == ak,
                    None => true,
                };
                if !fine {
                    return false;
                }
            }
            true
        };
        let mut budget = DOMINANCE_BUDGET;
        let mut changed = false;
        for k in 0..n {
            if !self.col_alive[k] || self.cols[k].is_empty() || budget == 0 {
                continue;
            }
            let Some(&(shortest, _)) = self.cols[k].iter().min_by_key(|&&(i, _)| self.rows[i].len()) else {
                continue;
            };
            if self.rows[shortest].len() > DOMINANCE_ROW {
                continue;
            }
            let candidates: Vec<usize> = self.rows[shortest].iter().map(|&(j, _)| j).collect();
            for j in candidates {
                if j == k || !self.col_alive[j] || !self.col_alive[k] {
                    continue;
                }
                if self.obj[j] > self.obj[k] {
                    continue;
                }
                // Which of the two moves would fix anything. Each moves an amount the
                // column it empties decides, so that amount is fractional when that
                // column is continuous and must not land on an integer one.
                let lower_k = roomy_up[j]
                    && self.col_lb[k].is_finite()
                    && (self.integer(k) || !self.integer(j));
                let raise_j = roomy_down[k]
                    && self.col_ub[j].is_finite()
                    && (self.integer(j) || !self.integer(k));
                if !lower_k && !raise_j {
                    continue;
                }
                budget = budget.saturating_sub(self.cols[j].len() + self.cols[k].len());
                if !dominates(&self.cols, j, k) {
                    continue;
                }
                if lower_k {
                    self.fix(k, self.col_lb[k]);
                } else {
                    self.fix(j, self.col_ub[j]);
                }
                self.stats.dominated_pairs += 1;
                changed = true;
                if !self.col_alive[k] {
                    break;
                }
            }
        }
        changed
    }

    /// A row that pins every other column in it at both values of an integer column
    /// with only two.
    ///
    /// With `y` in `{l, l + 1}` the row is two rows on the rest, one for each value.
    /// When each is forcing, holding the rest at the ends of their bounds that its
    /// one feasible side needs, every other column takes one value at `y = l` and one
    /// at `y = l + 1`, so it is the affine function of `y` through those two, at every
    /// integer point. Each is substituted out, and what is left of the row then holds
    /// at both values and goes too. This is a reduction for the MIP, not for its
    /// relaxation, which allows `y` between the two and the rest anything summing to
    /// the row: it tightens the relaxation as well as shrinking it.
    ///
    /// `neos-1445765` has rows `x_1 + ... + x_k = k y` with every `x` in `[0, 1]`:
    /// `y = 0` holds them all at 0 and `y = 1` at 1, so each is `y`. That removes all
    /// 18467 of its continuous columns, and HiGHS's base rules remove the same.
    fn sweep_switched(&mut self) -> Result<bool, Infeasible> {
        let mut changed = false;
        for i in 0..self.rows.len() {
            if !self.row_alive[i] || self.rows[i].len() < 2 {
                continue;
            }
            let Some((y, values)) = self.switch_of(i) else {
                continue;
            };
            let l = self.col_lb[y];
            for (k, (at_l, at_u)) in values {
                self.follow(k, y, l, at_l, at_u);
                self.stats.switched += 1;
            }
            // Holds at both values of y by construction; what rounding left of it is
            // not a constraint.
            self.kill_row(i);
            changed = true;
        }
        Ok(changed)
    }

    /// Substitute column `k` out as the affine function of `y` that is `at_l` at
    /// `y = l` and `at_u` at `y = l + 1`.
    fn follow(&mut self, k: usize, y: usize, l: f64, at_l: f64, at_u: f64) {
        // x_k = at_l + d (y - l).
        let d = at_u - at_l;
        let alpha = at_l - d * l;
        self.records.push(Record::Row {
            column: k,
            coefficient: 1.0,
            terms: vec![(y, -d)],
            lb: alpha,
            ub: alpha,
            col_lb: self.col_lb[k],
            col_ub: self.col_ub[k],
            integer: self.integer(k),
        });
        let entries = self.cols[k].clone();
        for &(r, a) in &entries {
            self.shift(r, a * alpha);
        }
        for (r, a) in entries {
            self.add(r, y, a * d);
        }
        self.obj[y] += self.obj[k] * d;
        self.offset += self.obj[k] * alpha;
        self.kill_col(k);
    }

    /// Substitute out the columns probing found to follow a binary. A link is
    /// skipped when either end has gone, or when the binary no longer has both
    /// values, since then propagation has already fixed what it would decide.
    fn apply_links(&mut self, links: &[Link]) -> bool {
        let mut changed = false;
        for link in links {
            let (k, y) = (link.column, link.switch);
            if k == y || !self.col_alive[k] || !self.col_alive[y] {
                continue;
            }
            if !self.integer(y) || self.col_lb[y] != 0.0 || self.col_ub[y] != 1.0 {
                continue;
            }
            if self.col_lb[k] >= self.col_ub[k] {
                continue;
            }
            self.follow(k, y, 0.0, link.at_zero, link.at_one);
            self.stats.linked += 1;
            changed = true;
        }
        changed
    }

    /// For [`Work::sweep_switched`]: the two-valued integer column of row `i` at
    /// both of whose values the row is forcing, and each other column's value at
    /// `y = l` and at `y = l + 1`.
    fn switch_of(&self, i: usize) -> Option<(usize, Vec<(usize, (f64, f64))>)> {
        let row = &self.rows[i];
        if row.iter().any(|&(k, _)| !self.col_lb[k].is_finite() || !self.col_ub[k].is_finite()) {
            return None;
        }
        // Activity range of the whole row, to take one column's share out of.
        let (mut low, mut high) = (0.0, 0.0);
        for &(k, a) in row {
            let (p, q) = (a * self.col_lb[k], a * self.col_ub[k]);
            low += p.min(q);
            high += p.max(q);
        }
        for &(y, b) in row {
            if !self.integer(y) || self.col_ub[y] - self.col_lb[y] != 1.0 {
                continue;
            }
            let (p, q) = (b * self.col_lb[y], b * self.col_ub[y]);
            let (rest_low, rest_high) = (low - p.min(q), high - p.max(q));
            if rest_high - rest_low <= self.tolerance {
                continue;
            }
            // At each value of y: true when the rest is held at its least activity,
            // false at its most, None when the row does not force it.
            let forced = |v: f64| {
                let (lo, hi) = (self.row_lb[i] - b * v, self.row_ub[i] - b * v);
                let near = |a: f64, t: f64| (a - t).abs() <= self.tolerance * t.abs().max(1.0);
                if hi.is_finite() && near(rest_low, hi) {
                    Some(true)
                } else if lo.is_finite() && near(rest_high, lo) {
                    Some(false)
                } else {
                    None
                }
            };
            let (Some(at_l), Some(at_u)) = (forced(self.col_lb[y]), forced(self.col_ub[y])) else {
                continue;
            };
            let end = |k: usize, a: f64, least: bool| if (a > 0.0) == least { self.col_lb[k] } else { self.col_ub[k] };
            let values = row
                .iter()
                .filter(|&&(k, _)| k != y)
                .map(|&(k, a)| (k, (end(k, a, at_l), end(k, a, at_u))))
                .collect();
            return Some((y, values));
        }
        None
    }

    fn finish(self) -> (Problem, Postsolve) {
        let problem = self.problem;
        let kept_columns: Vec<usize> = (0..self.cols.len()).filter(|&j| self.col_alive[j]).collect();
        let kept_rows: Vec<usize> = (0..self.rows.len()).filter(|&i| self.row_alive[i]).collect();
        let mut row_of = vec![usize::MAX; self.rows.len()];
        for (new, &old) in kept_rows.iter().enumerate() {
            row_of[old] = new;
        }
        let triplets = kept_columns.iter().enumerate().flat_map(|(new_j, &j)| {
            let row_of = &row_of;
            self.cols[j].iter().map(move |&(i, a)| (row_of[i], new_j, a))
        });
        let matrix = SparseMatrix::from_triplets(kept_rows.len(), kept_columns.len(), triplets);
        let reduced = Problem {
            name: problem.name.clone(),
            sense: problem.sense,
            obj: kept_columns.iter().map(|&j| self.obj[j]).collect(),
            obj_offset: self.offset,
            matrix,
            row_lb: kept_rows.iter().map(|&i| self.row_lb[i]).collect(),
            row_ub: kept_rows.iter().map(|&i| self.row_ub[i]).collect(),
            col_lb: kept_columns.iter().map(|&j| self.col_lb[j]).collect(),
            col_ub: kept_columns.iter().map(|&j| self.col_ub[j]).collect(),
            col_type: kept_columns.iter().map(|&j| problem.col_type[j]).collect(),
            col_names: kept_columns.iter().map(|&j| problem.col_names[j].clone()).collect(),
            row_names: kept_rows.iter().map(|&i| problem.row_names[i].clone()).collect(),
        };
        let postsolve = Postsolve {
            columns: kept_columns,
            records: self.records,
            n_original: problem.n_cols(),
        };
        (reduced, postsolve)
    }
}

/// Alternate the cheap half of `presolve` with [`reduce`] until neither finds more.
///
/// The cheap reductions settle the rows an elimination leaves behind (a forcing row,
/// a redundant one, a bound that now propagates) and those can expose the next
/// elimination, so a few rounds go further than either alone: `binkar10_1` keeps 91%
/// of its columns after presolve and 66% after this.
pub fn reduce_with_presolve(
    problem: &Problem,
    tolerance: f64,
) -> Result<Option<(Problem, Postsolve, Stats)>, Infeasible> {
    let mut model = problem.clone();
    let mut map: Option<Postsolve> = None;
    let mut total = Stats::default();
    let mut probes = 0;
    for _ in 0..ROUNDS + PROBE_ROUNDS + 1 {
        if crate::presolve::presolve_cheap(&mut model, 20) == crate::presolve::Outcome::Infeasible {
            return Err(Infeasible);
        }
        let mut reduced = reduce(&model, tolerance)?;
        if reduced.is_none() && probes < PROBE_ROUNDS {
            // The cheap rules have run dry: probe, for the columns that follow a
            // binary, which only a renumbering pass can take out, and for the fixings
            // and bounds that set the cheap rules going again.
            let mut links = Vec::new();
            if probes == 0 {
                let (outcome, found) = crate::presolve::presolve_links(&mut model, 20, None);
                if outcome == crate::presolve::Outcome::Infeasible {
                    return Err(Infeasible);
                }
                links = found;
            }
            let Some((more, found)) = crate::probe::probe_bounds(&mut model, 0.5f64.powi(probes as i32))
            else {
                return Err(Infeasible);
            };
            probes += 1;
            if found.fixed + found.tightened + found.links == 0 {
                // Another pass would read the same model the same way.
                probes = PROBE_ROUNDS;
            }
            links.extend(more);
            reduced = reduce_linked(&model, tolerance, &links)?;
        }
        let Some((smaller, next, stats)) = reduced else {
            break;
        };
        total.absorb(&stats);
        model = smaller;
        map = Some(match map {
            None => next,
            Some(map) => map.then(next),
        });
    }
    Ok(map.map(|map| (model, map, total)))
}

/// Rounds of [`reduce_with_presolve`]; the third rarely finds anything.
const ROUNDS: usize = 4;

/// Passes of probing [`reduce_with_presolve`] may make, each with half the budget of
/// the one before, while each finds something. On `ex10` a pass's fixings let the rules
/// remove a fifth of the model, which leaves the next pass less to read.
const PROBE_ROUNDS: usize = 3;

/// Reduce `problem` to a smaller model and the map back.
///
/// Returns `None` when nothing was removed, so the caller keeps the model it has
/// rather than paying for a copy of it.
pub fn reduce(
    problem: &Problem,
    tolerance: f64,
) -> Result<Option<(Problem, Postsolve, Stats)>, Infeasible> {
    reduce_linked(problem, tolerance, &[])
}

/// [`reduce`], first substituting out the columns `links` says follow a binary.
pub fn reduce_linked(
    problem: &Problem,
    tolerance: f64,
    links: &[Link],
) -> Result<Option<(Problem, Postsolve, Stats)>, Infeasible> {
    let mut work = Work::new(problem, tolerance);
    let mut any = work.apply_links(links);
    for _ in 0..MAX_SWEEPS {
        let mut changed = work.sweep_rows()?;
        changed |= work.sweep_columns()?;
        changed |= work.sweep_singletons();
        changed |= work.sweep_doubletons()?;
        changed |= work.sweep_switched()?;
        if !changed {
            changed |= work.sweep_substitute()?;
        }
        if !changed {
            // The most expensive, and only once the cheap rules have run dry.
            changed |= work.sweep_dominated()?;
            changed |= work.sweep_dominance();
        }
        if !changed {
            break;
        }
        any = true;
    }
    if !any {
        return Ok(None);
    }
    let stats = work.stats;
    let (reduced, postsolve) = work.finish();
    Ok(Some((reduced, postsolve, stats)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{Builder, RowSense, Sense};
    use crate::search::{self, Options, Status};

    /// A random mixed model, seeded with the structure each rule looks for.
    ///
    /// Continuous columns alone in a row (with and without a cost, in equations and
    /// inequalities), two-column equations of every type pairing and with non-unit
    /// ratios, columns no row restrains, and fixed columns: everything here should be
    /// reached by some seed, and the test checks that it is.
    fn model(seed: u64) -> Problem {
        let mut state = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        let mut b = Builder::new(if next() % 2 == 0 { Sense::Minimize } else { Sense::Maximize });
        let n = 5 + (next() % 5) as usize;
        let mut columns = Vec::new();
        let mut bounds = Vec::new();
        for j in 0..n {
            let (lo, hi) = match next() % 5 {
                0 => (0.0, 1.0),
                1 => (-2.0, 4.0),
                2 => (0.0, 5.0),
                3 => (-3.0, 3.0),
                // No upper bound, which is what gives the duals a bound.
                _ => (0.0, f64::INFINITY),
            };
            let column = if lo == -2.0 || hi == 1.0 {
                b.integer(&format!("x{j}"), lo, hi)
            } else {
                b.continuous(&format!("x{j}"), lo, hi)
            };
            columns.push(column);
            bounds.push((lo, hi));
        }
        let mut objective: Vec<(usize, f64)> =
            columns.iter().map(|&j| (j, (next() % 9) as f64 - 4.0)).collect();
        // Sometimes a column with bounds too wide to matter, defined by an equation on
        // a few finite columns and used in other rows: one to substitute out.
        let wide = (next() % 2 == 0).then(|| {
            let w = b.continuous("wide", -100.0, 100.0);
            objective.push((w, (next() % 5) as f64 - 2.0));
            w
        });
        // Sometimes a pair of binaries two inequalities make complements, which no one
        // row says and only probing finds: `p + q <= 1`, `p + q >= 1`.
        let pair = (next() % 2 == 0).then(|| {
            let (p, q) = (b.binary("p"), b.binary("q"));
            objective.push((p, (next() % 5) as f64 - 2.0));
            objective.push((q, (next() % 5) as f64 - 2.0));
            (p, q)
        });
        // Sometimes a near copy of the first column, for one to dominate the other.
        let twin = (next() % 2 == 0).then(|| {
            let t = if next() % 2 == 0 { b.continuous("twin", 0.0, 5.0) } else { b.integer("twin", 0.0, 3.0) };
            objective.push((t, objective[0].1 + (next() % 3) as f64 - 1.0));
            t
        });
        let sense = |k: u64| match k % 3 {
            0 => RowSense::Le,
            1 => RowSense::Ge,
            _ => RowSense::Eq,
        };
        for _ in 0..2 + (next() % 4) as usize {
            let mut terms: Vec<(usize, f64)> = Vec::new();
            for &j in &columns {
                if next() % 2 == 0 {
                    terms.push((j, (next() % 7) as f64 - 3.0));
                }
            }
            terms.retain(|&(_, a)| a != 0.0);
            if terms.is_empty() {
                continue;
            }
            if let Some((p, q)) = pair
                && next() % 3 == 0
            {
                terms.push((if next() % 2 == 0 { p } else { q }, (next() % 3) as f64 + 1.0));
            }
            if let Some(w) = wide
                && next() % 3 == 0
            {
                terms.push((w, if next() % 2 == 0 { 1.0 } else { -2.0 }));
            }
            if let (Some(t), Some(&(_, a))) = (twin, terms.iter().find(|&&(j, _)| j == columns[0])) {
                let nudge = [0.0, 0.0, 1.0, -1.0][(next() % 4) as usize];
                if a + nudge != 0.0 {
                    terms.push((t, a + nudge));
                }
            }
            // A right-hand side the origin nearly satisfies, so most models are feasible.
            let rhs = (next() % 7) as f64 - 3.0;
            let s = sense(next());
            // Sometimes a private continuous column in it: a slack, or a cost-carrying
            // singleton that may be implied free.
            if next() % 2 == 0 {
                let (lo, hi) = if next() % 2 == 0 { (0.0, 20.0) } else { (-1.0, 1.5) };
                let k = b.continuous(&format!("s{}", columns.len() + terms.len()), lo, hi);
                let c = if next() % 2 == 0 { 0.0 } else { (next() % 5) as f64 - 2.0 };
                objective.push((k, c));
                terms.push((k, if next() % 2 == 0 { 1.0 } else { -2.5 }));
            }
            b.row(&terms, s, rhs);
        }
        for _ in 0..(next() % 4) {
            let p = columns[(next() % n as u64) as usize];
            let q = columns[(next() % n as u64) as usize];
            if p == q {
                continue;
            }
            let ratio = [1.0, -1.0, 2.0, -3.0, 0.5][(next() % 5) as usize];
            let rhs = (next() % 5) as f64 - 2.0;
            b.row(&[(p, 1.0), (q, ratio)], RowSense::Eq, rhs);
        }
        if next() % 3 == 0 {
            // A column fixed by its bounds, and one in no row at all.
            let f = b.continuous("fixed", 1.5, 1.5);
            let e = b.integer("empty", -1.0, 3.0);
            objective.push((f, 2.0));
            objective.push((e, (next() % 3) as f64 - 1.0));
            let j = columns[0];
            b.row(&[(f, 1.0), (j, 1.0)], RowSense::Le, 4.0);
        }
        if let Some((p, q)) = pair {
            b.row(&[(p, 1.0), (q, 1.0)], RowSense::Le, 1.0);
            b.row(&[(p, 1.0), (q, 1.0)], RowSense::Ge, 1.0);
        }
        if let Some(w) = wide {
            let mut terms = vec![(w, if next() % 2 == 0 { 1.0 } else { -1.5 })];
            for (&j, &(_, hi)) in columns.iter().zip(&bounds) {
                if hi.is_finite() && next() % 2 == 0 {
                    terms.push((j, (next() % 5) as f64 - 2.0));
                }
            }
            terms.retain(|&(_, a)| a != 0.0);
            b.row(&terms, RowSense::Eq, (next() % 5) as f64 - 2.0);
        }
        if next() % 2 == 0 {
            // A row a binary switches between holding the rest at their least and at
            // their most: `a x + (least - most) y = least`.
            let y = b.binary("switch");
            objective.push((y, (next() % 5) as f64 - 2.0));
            let mut terms = Vec::new();
            let (mut least, mut most) = (0.0, 0.0);
            for (&j, &(lo, hi)) in columns.iter().zip(&bounds) {
                if hi.is_finite() && next() % 3 == 0 {
                    let a = [1.0, -1.0, 2.0][(next() % 3) as usize];
                    let (p, q) = (a * lo, a * hi);
                    least += f64::min(p, q);
                    most += f64::max(p, q);
                    terms.push((j, a));
                }
            }
            if !terms.is_empty() && most > least {
                terms.push((y, least - most));
                b.row(&terms, RowSense::Eq, least);
            }
        }
        b.objective(&objective);
        b.build()
    }

    fn plain() -> Options {
        Options { aggregate: false, ..Options::default() }
    }

    #[test]
    fn a_reduced_model_has_the_same_optimum_and_the_map_recovers_it() {
        let mut reduced_any = 0;
        let mut totals = Stats::default();
        for seed in 0..600u64 {
            let problem = model(seed);
            // A column with no upper bound can make the model unbounded, which neither
            // search reports as such and which leaves nothing to compare.
            if crate::lp::Lp::relaxation(&problem).solve().status == crate::lp::LpStatus::Unbounded {
                continue;
            }
            let direct = search::solve(&problem, plain());
            let reduction = reduce(&problem, 1e-9);
            let (smaller, map, stats) = match reduction {
                Ok(Some(r)) => r,
                Ok(None) => continue,
                Err(Infeasible) => {
                    assert_eq!(direct.status, Status::Infeasible, "seed {seed}: reduce said infeasible");
                    continue;
                }
            };
            reduced_any += 1;
            totals.fixed_columns += stats.fixed_columns;
            totals.row_singletons += stats.row_singletons;
            totals.dual_fixed += stats.dual_fixed;
            totals.column_singletons += stats.column_singletons;
            totals.doubletons += stats.doubletons;
            totals.dominated += stats.dominated;
            totals.dominated_pairs += stats.dominated_pairs;
            totals.switched += stats.switched;
            totals.substituted += stats.substituted;
            let small = search::solve(&smaller, plain());
            assert_eq!(small.status, direct.status, "seed {seed}: status, {stats:?}");
            match (small.objective, direct.objective) {
                (Some(a), Some(b)) => {
                    assert!(
                        (a - b).abs() <= 1e-6 * b.abs().max(1.0),
                        "seed {seed}: reduced {a}, original {b}, {stats:?}"
                    );
                    let expanded = map.expand(&small.x);
                    assert_eq!(expanded.len(), problem.n_cols(), "seed {seed}: width");
                    let scored: f64 = problem.obj.iter().zip(&expanded).map(|(c, v)| c * v).sum();
                    let scored = problem.objective_value(scored + problem.obj_offset);
                    assert!(
                        (scored - a).abs() <= 1e-6 * a.abs().max(1.0),
                        "seed {seed}: the expanded point scores {scored} where the reduced \
                         search claimed {a}, {stats:?}"
                    );
                    assert!(
                        crate::heuristic::is_feasible(&problem, &expanded, 1e-6),
                        "seed {seed}: the expanded point is not feasible for the original, \
                         {stats:?}: {expanded:?}"
                    );
                    for j in problem.integer_columns() {
                        assert!(
                            (expanded[j] - expanded[j].round()).abs() <= 1e-6,
                            "seed {seed}: column {j} came back fractional at {}",
                            expanded[j]
                        );
                    }
                }
                (None, None) => {}
                (a, b) => panic!("seed {seed}: reduced {a:?}, original {b:?}, {stats:?}"),
            }
        }
        assert!(reduced_any > 300, "only {reduced_any} models reduced");
        // Every rule has to have been reached, or this checked less than it says.
        for (rule, count) in [
            ("fixed", totals.fixed_columns),
            ("row singleton", totals.row_singletons),
            ("dual fixing", totals.dual_fixed),
            ("column singleton", totals.column_singletons),
            ("doubleton", totals.doubletons),
            ("dominated", totals.dominated),
            ("dominated pair", totals.dominated_pairs),
            ("switched", totals.switched),
            ("substituted", totals.substituted),
        ] {
            assert!(count > 20, "{rule} was reached only {count} times");
        }
    }

    /// The pass the solver runs: presolve, probing and `reduce` in rounds, with the
    /// columns probing links to a binary substituted out between them.
    #[test]
    fn the_full_pass_has_the_same_optimum_and_the_map_recovers_it() {
        let mut totals = Stats::default();
        for seed in 0..600u64 {
            let problem = model(seed);
            if crate::lp::Lp::relaxation(&problem).solve().status == crate::lp::LpStatus::Unbounded {
                continue;
            }
            let direct = search::solve(&problem, plain());
            let (smaller, map, stats) = match reduce_with_presolve(&problem, 1e-9) {
                Ok(Some(r)) => r,
                Ok(None) => continue,
                Err(Infeasible) => {
                    assert_eq!(direct.status, Status::Infeasible, "seed {seed}: said infeasible");
                    continue;
                }
            };
            totals.absorb(&stats);
            let small = search::solve(&smaller, plain());
            assert_eq!(small.status, direct.status, "seed {seed}: status, {stats:?}");
            if let (Some(a), Some(b)) = (small.objective, direct.objective) {
                assert!((a - b).abs() <= 1e-6 * b.abs().max(1.0), "seed {seed}: {a} against {b}, {stats:?}");
                let expanded = map.expand(&small.x);
                let scored: f64 = problem.obj.iter().zip(&expanded).map(|(c, v)| c * v).sum();
                let scored = problem.objective_value(scored + problem.obj_offset);
                assert!((scored - a).abs() <= 1e-6 * a.abs().max(1.0), "seed {seed}: scores {scored}, {stats:?}");
                assert!(
                    crate::heuristic::is_feasible(&problem, &expanded, 1e-6),
                    "seed {seed}: infeasible for the original, {stats:?}: {expanded:?}"
                );
                for j in problem.integer_columns() {
                    assert!((expanded[j] - expanded[j].round()).abs() <= 1e-6, "seed {seed}: column {j} fractional");
                }
            }
        }
        assert!(totals.linked > 20, "probing linked only {} columns", totals.linked);
    }

    #[test]
    fn a_doubleton_with_a_fractional_ratio_keeps_the_integer_column() {
        // 2x - y = 0 on integers: x = y / 2 is not integral for odd y, so x cannot go;
        // y = 2x can, and the optimum has to survive the trip.
        let mut b = Builder::new(Sense::Maximize);
        let x = b.column("x", 0.0, 10.0, VarType::Integer);
        let y = b.column("y", 0.0, 7.0, VarType::Integer);
        b.objective(&[(x, 1.0), (y, 1.0)]);
        b.row(&[(x, 2.0), (y, -1.0)], RowSense::Eq, 0.0);
        let problem = b.build();
        let (smaller, map, stats) = reduce(&problem, 1e-9).unwrap().unwrap();
        assert_eq!(stats.doubletons, 1);
        // With y gone x is in no row, and its cost puts it at its bound, so nothing is
        // left; what matters is that x was the one kept, and that its bound came from y.
        assert!(!map.columns().contains(&y));
        let solved = search::solve(&smaller, plain());
        let expanded = map.expand(&solved.x);
        assert_eq!(expanded, vec![3.0, 6.0], "y <= 7 leaves x at most 3");
    }

    #[test]
    fn an_implied_free_singleton_turns_its_row_into_an_equation() {
        // min s with s >= 3 - x - y, s in [-10, 10], x and y binary: s is implied free
        // (3 - x - y is in [1, 3]), so the row holds with equality and s goes.
        let mut b = Builder::new(Sense::Minimize);
        let x = b.binary("x");
        let y = b.binary("y");
        let s = b.continuous("s", -10.0, 10.0);
        b.objective(&[(s, 1.0), (x, 0.5), (y, 0.25)]);
        b.row(&[(s, 1.0), (x, 1.0), (y, 1.0)], RowSense::Ge, 3.0);
        let problem = b.build();
        let (smaller, map, stats) = reduce(&problem, 1e-9).unwrap().unwrap();
        assert_eq!(stats.column_singletons, 1);
        let solved = search::solve(&smaller, plain());
        let expanded = map.expand(&solved.x);
        // s = 3 - x - y, cost 3 - 0.5x - 0.75y: both binaries at one.
        assert_eq!(expanded, vec![1.0, 1.0, 1.0]);
        assert_eq!(solved.objective, Some(1.75));
    }
}
