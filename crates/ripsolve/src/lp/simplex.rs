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

/// Factorizations kept per LP.
///
/// Swept over 1, 4, 8, 32, 64, 128 and 256: solve time falls all the way up, but so
/// does memory, and each worker in a parallel search holds its own cache. 64 keeps
/// most of the gain — v064c1000n100 11.2s to 9.2s, v128c1000n100 5.1s to 4.1s — at
/// 50MB on one thread and 200MB across sixteen. Larger is faster if memory is free.
const FACTOR_CACHE: usize = 64;

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
    /// Product-form updates to accumulate before refactorizing.
    ///
    /// Trades two costs against each other: every solve replays the whole eta file,
    /// so a long interval makes solves expensive, while refactorizing is a fixed
    /// cost amortized over the interval.
    pub refactor_interval: usize,
}

impl Default for Tolerances {
    fn default() -> Self {
        Self {
            primal_feasibility: 1e-7,
            dual_feasibility: 1e-7,
            pivot: 1e-9,
            refactor_interval: 50,
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
    /// The dual objective rose past the caller's cutoff, so the true optimum is
    /// worse than a solution the caller already holds. Only the dual simplex can
    /// report this, since only there does the objective climb monotonically.
    CutOff,
}

/// A basis, saved so a related LP can start from it instead of from scratch.
///
/// Branch-and-bound lives on this: a child differs from its parent only in one
/// column's bounds, which leaves the parent's basis dual feasible. Re-solving with
/// the dual simplex from here typically takes a handful of pivots, where a cold
/// start would repeat the whole phase-1.
#[derive(Clone, Debug)]
pub struct BasisState {
    basic: Vec<usize>,
    status: Vec<Status>,
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
    /// The final basis, for warm-starting a related solve.
    pub basis: BasisState,
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

/// What a warm-started solve should do first.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Entry {
    /// Repair primal infeasibility with the dual simplex, keeping dual feasibility.
    Dual,
    /// Repair primal infeasibility with phase 1, then optimize.
    Primal,
}

/// An LP in computational form, ready to solve.
///
/// Cloneable so that each thread of a parallel search can hold its own, since
/// solving a node mutates the column bounds.
pub struct Lp {
    n_structural: usize,
    m: usize,
    cost: Vec<f64>,
    lower: Vec<f64>,
    upper: Vec<f64>,
    /// The original `A`, column-major, `m x n_structural`.
    matrix: SparseMatrix,
    /// Which structural columns must take integer values. Logicals are continuous.
    integer: Vec<bool>,
    tol: Tolerances,
    /// Recently used factorizations, most recent first, keyed by basis columns.
    ///
    /// Branch and bound re-solves the same basis constantly: a child differs from
    /// its parent only in a column's *bounds*, which do not enter the basis matrix
    /// at all, so the parent's factors are valid for the child verbatim. Without
    /// this, every warm start refactorized from scratch — measured at 9.1ms of a
    /// 10.9ms node on a 1000-row model, against 1.8ms of actual pivoting.
    ///
    /// Several entries rather than one, because best-bound node selection does not
    /// visit the tree in an order that keeps a single entry warm: consecutive nodes
    /// are usually unrelated, and a one-entry cache measured an 8-20% hit rate. What
    /// does recur is siblings, which share a parent's basis and tend to be reached
    /// near each other.
    factors: Vec<(Vec<usize>, Basis)>,
}

impl Clone for Lp {
    /// The cache is deliberately not cloned: it is a scratch optimization, and a
    /// fresh clone has no basis history worth carrying.
    fn clone(&self) -> Self {
        Self {
            n_structural: self.n_structural,
            m: self.m,
            cost: self.cost.clone(),
            lower: self.lower.clone(),
            upper: self.upper.clone(),
            matrix: self.matrix.clone(),
            integer: self.integer.clone(),
            tol: self.tol,
            factors: Vec::new(),
        }
    }
}

impl Lp {
    /// Build the LP relaxation of a binary program.
    ///
    /// Column bounds are taken from the problem, so a problem whose bounds have
    /// already been narrowed by presolve or branching relaxes accordingly.
    pub fn relaxation(problem: &Problem) -> Lp {
        let integer: Vec<bool> = (0..problem.n_cols())
            .map(|j| problem.is_integer(j))
            .collect();
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
            integer,
            tol: Tolerances::default(),
            factors: Vec::new(),
        }
    }

    pub fn with_tolerances(mut self, tol: Tolerances) -> Self {
        self.tol = tol;
        self
    }

    /// Narrow a structural column's bounds, as branching does.
    pub fn set_column_bounds(&mut self, j: usize, lo: f64, hi: f64) {
        debug_assert!(j < self.n_structural);
        self.lower[j] = lo;
        self.upper[j] = hi;
    }

    /// Discard any cached factorization.
    ///
    /// Only needed if the constraint matrix changes; bounds and costs do not affect
    /// the basis matrix, which is the whole reason the cache is worth having.
    pub fn invalidate_factors(&mut self) {
        self.factors.clear();
    }

    /// Replace the structural objective coefficients.
    ///
    /// The feasibility pump repeatedly re-optimizes the same constraint set under a
    /// distance objective, so it needs to swap the costs without rebuilding the LP.
    /// Logical columns keep their zero cost.
    pub fn set_costs(&mut self, costs: &[f64]) {
        debug_assert_eq!(costs.len(), self.n_structural);
        self.cost[..self.n_structural].copy_from_slice(costs);
    }

    /// The structural objective coefficients, for saving and restoring.
    pub fn costs(&self) -> &[f64] {
        &self.cost[..self.n_structural]
    }

    /// The current bounds of a structural column.
    pub fn column_bounds(&self, j: usize) -> (f64, f64) {
        (self.lower[j], self.upper[j])
    }

    pub fn n_columns(&self) -> usize {
        self.n_structural
    }

    fn n_total(&self) -> usize {
        self.n_structural + self.m
    }

    /// Column `j` of `[A | -I]` as `(rows, values)`.
    ///
    /// Refactorization wants the sparsity, not a dense scatter — the whole point of
    /// the LU is to touch only the nonzeros.
    fn column_sparse(&self, j: usize) -> (Vec<usize>, Vec<f64>) {
        if j < self.n_structural {
            let (rows, vals) = self.matrix.column(j);
            (rows.to_vec(), vals.to_vec())
        } else {
            (vec![j - self.n_structural], vec![-1.0])
        }
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

    /// Solve from a cold start, returning the optimum of the internal minimization
    /// form.
    pub fn solve(&mut self) -> LpSolution {
        self.solve_with_limit(200_000)
    }

    pub fn solve_with_limit(&mut self, max_iterations: usize) -> LpSolution {
        self.run_solver(None, None, max_iterations)
    }

    /// Generate Gomory mixed-integer cuts from the tableau at `basis`.
    ///
    /// Returned in terms of the *structural* columns, as `(coefficients, lower
    /// bound)` pairs meaning `sum a_j x_j >= lb`.
    ///
    /// # Derivation
    ///
    /// Each tableau row reads `z_B(i) + sum_j a_j w_j = beta`, where `w_j >= 0` is
    /// the nonbasic's distance from the bound it sits on — `z_j - l_j` at a lower
    /// bound, `u_j - z_j` at an upper one. That shift is what lets a
    /// bounded-variable simplex use the textbook formula, which assumes nonbasics at
    /// zero, and it is also where the sign of `a_j` flips for an at-upper column.
    ///
    /// With `f0` the fractional part of `beta`, the cut is `sum_j c_j w_j >= 1` for
    /// coefficients `c_j` that differ between integer and continuous `w_j`. Every
    /// structural here is binary and so integer; the logicals are continuous.
    ///
    /// Substituting `w_j` back gives a row over `z = [x; s]`, and the logical part
    /// is eliminated with `s = A x`, leaving a cut over the structural columns only.
    pub fn gomory_cuts(
        &self,
        basis: &BasisState,
        max_cuts: usize,
    ) -> Vec<(Vec<(usize, f64)>, f64)> {
        let mut solver = Solver::warm(self, basis, None);
        if solver.refactorize().is_err() {
            return Vec::new();
        }
        solver.gomory_cuts(max_cuts)
    }

    /// Re-solve starting from a saved basis.
    ///
    /// `cutoff` lets the solve abandon a node early: the dual simplex's objective
    /// only ever climbs, so once it passes the cutoff the node's true bound must be
    /// worse still and the result is [`LpStatus::CutOff`].
    pub fn solve_warm(
        &mut self,
        start: &BasisState,
        cutoff: Option<f64>,
        max_iterations: usize,
    ) -> LpSolution {
        self.run_solver(Some(start), cutoff, max_iterations)
    }

    /// Shared body of the solve entry points, wrapping the factorization cache.
    ///
    /// The cache hits when the requested basis is the one the last solve *ended*
    /// on, which in a tree search is the common case: both children of a node ask
    /// for their parent's final basis, and a dive or a probe asks for the basis it
    /// just produced.
    fn run_solver(
        &mut self,
        start: Option<&BasisState>,
        cutoff: Option<f64>,
        max_iterations: usize,
    ) -> LpSolution {
        let factors = start.and_then(|start| {
            let found = self
                .factors
                .iter()
                .position(|(columns, _)| *columns == start.basic)?;
            Some(self.factors.remove(found).1)
        });

        let (solution, basic, basis) = {
            let lp = &*self;
            match start {
                Some(start) => Solver::warm(lp, start, factors).run(max_iterations, cutoff),
                None => Solver::new(lp).run(max_iterations, None),
            }
        };
        // Most recent first, oldest evicted. The key comparison is a slice equality
        // over the basis columns, which is negligible against the cost of the
        // factorization it avoids.
        self.factors.insert(0, (basic, basis));
        self.factors.truncate(FACTOR_CACHE);
        solution
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
    /// Which method to try first; see [`Solver::run`].
    entry_hint: Entry,
    /// True when `basis` already factorizes `basic`, so entry can skip the rebuild.
    factorized: bool,
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
            entry_hint: Entry::Primal,
            factorized: false,
        };
        solver.recompute_basic_values();
        solver
    }

    /// Rebuild a solver around a saved basis.
    ///
    /// Nonbasic variables are re-parked on their (possibly just-changed) bounds. A
    /// variable that has become fixed is pulled onto the fixed value regardless of
    /// which side it was parked on, so the recomputed basic values reflect the new
    /// bounds rather than the parent's.
    fn warm(lp: &'a Lp, start: &BasisState, factors: Option<Basis>) -> Solver<'a> {
        let m = lp.m;
        let mut z = vec![0.0; lp.n_total()];
        let mut status = start.status.clone();

        for (j, st) in status.iter_mut().enumerate() {
            if let Status::NonBasic(at) = st {
                // Re-park anything sitting on a bound that is now infinite, which can
                // happen only if the caller widened bounds rather than narrowing them.
                let mut here = *at;
                if !lp.value_at(j, here).is_finite() {
                    here = match here {
                        At::Lower if lp.upper[j].is_finite() => At::Upper,
                        At::Upper if lp.lower[j].is_finite() => At::Lower,
                        _ => At::Zero,
                    };
                    *at = here;
                }
                z[j] = lp.value_at(j, here);
            }
        }

        // Reusing the caller's factors is the point of the cache; without them the
        // solver falls back to factorizing on entry, as it always did.
        let factorized = factors.is_some();
        let mut solver = Solver {
            lp,
            basis: factors.unwrap_or_else(|| Basis::all_logical(m)),
            factorized,
            basic: start.basic.clone(),
            status,
            z,
            col: vec![0.0; m],
            alpha: Vec::new(),
            y: Vec::new(),
            cost_b: vec![0.0; m],
            rhs: vec![0.0; m],
            entry_hint: Entry::Dual,
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
            let columns: Vec<(Vec<usize>, Vec<f64>)> =
                self.basic.iter().map(|&j| lp.column_sparse(j)).collect();
            match self.basis.refactorize(&columns, lp.tol.pivot) {
                Ok(()) => {
                    self.factorized = true;
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

    /// Current value of the internal minimization objective.
    fn objective(&self) -> f64 {
        (0..self.lp.n_total())
            .map(|j| self.lp.cost[j] * self.z[j])
            .sum()
    }

    /// Is every nonbasic reduced cost on the correct side of zero?
    ///
    /// A basis inherited from a parent node is dual feasible, because changing a
    /// column's *bounds* leaves reduced costs untouched. This checks that rather
    /// than assuming it, so a caller passing an unrelated basis falls back to the
    /// primal method instead of getting a wrong answer from the dual one.
    fn is_dual_feasible(&mut self) -> bool {
        let tol = self.lp.tol.dual_feasibility;
        self.load_basic_costs(false);
        let mut y = std::mem::take(&mut self.y);
        self.basis.btran(&self.cost_b, &mut y);
        self.y = y;

        for j in 0..self.lp.n_total() {
            let Status::NonBasic(at) = self.status[j] else {
                continue;
            };
            // A fixed variable cannot move, so its reduced cost is unconstrained.
            if self.lp.lower[j] == self.lp.upper[j] {
                continue;
            }
            let d = self.reduced_cost(j, false);
            let ok = match at {
                At::Lower => d >= -tol,
                At::Upper => d <= tol,
                At::Zero => d.abs() <= tol,
            };
            if !ok {
                return false;
            }
        }
        true
    }

    /// The most primal-infeasible basic variable, with the bound it violates.
    ///
    /// Under `bland` it returns the *first* infeasible row by basic variable index
    /// instead of the worst, which is the leaving half of Bland's rule.
    fn most_infeasible_row(&self, bland: bool) -> Option<(usize, At)> {
        let tol = self.lp.tol.primal_feasibility;
        let mut worst = tol;
        let mut found: Option<(usize, At)> = None;
        for i in 0..self.lp.m {
            let j = self.basic[i];
            let v = self.z[j];
            let (below, above) = (self.lp.lower[j] - v, v - self.lp.upper[j]);
            let violation = below.max(above);
            if violation <= tol {
                continue;
            }
            let side = if below > above { At::Lower } else { At::Upper };
            if bland {
                // Lowest basic variable index wins, regardless of how bad it is.
                if found.is_none_or(|(prev, _)| j < self.basic[prev]) {
                    found = Some((i, side));
                }
            } else if violation > worst {
                worst = violation;
                found = Some((i, side));
            }
        }
        found
    }

    /// The dual simplex: repair primal infeasibility while holding dual feasibility.
    ///
    /// Each iteration takes the most primal-infeasible basic variable out of the
    /// basis at the bound it violates, and chooses the entering column by the dual
    /// ratio test — the smallest `|d_j / alpha_rj|` among columns that can move in
    /// the required direction. No such column means no assignment can repair the
    /// violation, which is exactly primal infeasibility.
    ///
    /// The objective climbs monotonically here, which is what makes `cutoff` sound:
    /// once it passes, the node's true bound can only be worse.
    fn run_dual(
        &mut self,
        max_iterations: usize,
        cutoff: Option<f64>,
        iterations: &mut usize,
    ) -> Option<LpStatus> {
        let tol = self.lp.tol.dual_feasibility;
        let pivot_tol = self.lp.tol.pivot;
        let mut rho: Vec<f64> = Vec::new();
        // Degenerate steps in a row before switching to Bland's rule. The primal
        // method has had this since it first cycled; the dual method had nothing,
        // and since every warm start enters through the dual method, that is nearly
        // every node in the search. MIPLIB's pk1 -- 86 columns -- burned 100,000
        // iterations across four nodes before the caller gave up.
        const STALL_LIMIT: usize = 100;
        let mut stalled = 0usize;

        while *iterations < max_iterations {
            let bland = stalled > STALL_LIMIT;
            if cutoff.is_some_and(|limit| self.objective() > limit) {
                return Some(LpStatus::CutOff);
            }

            // Refresh the duals every iteration. Each pivot changes the basis and so
            // changes every reduced cost; pricing the ratio test against duals left
            // over from a previous iteration silently picks the wrong entering column
            // and the solve converges to a point that is not the optimum.
            self.load_basic_costs(false);
            let mut y = std::mem::take(&mut self.y);
            self.basis.btran(&self.cost_b, &mut y);
            self.y = y;

            let Some((r, violated)) = self.most_infeasible_row(bland) else {
                // Primal feasible and still dual feasible, so this is the optimum.
                return Some(LpStatus::Optimal);
            };

            self.basis.btran_unit(r, &mut rho);
            // The leaving variable must rise to its lower bound, or fall to its upper.
            let rising = violated == At::Lower;

            let mut best: Option<(usize, f64)> = None;
            let mut best_ratio = f64::INFINITY;
            let mut best_pivot = 0.0;
            for j in 0..self.lp.n_total() {
                let Status::NonBasic(at) = self.status[j] else {
                    continue;
                };
                if self.lp.lower[j] == self.lp.upper[j] {
                    continue;
                }
                let arj = if j < self.lp.n_structural {
                    let (rows, vals) = self.lp.matrix.column(j);
                    rows.iter()
                        .zip(vals)
                        .map(|(&i, &a)| rho[i] * a)
                        .sum::<f64>()
                } else {
                    -rho[j - self.lp.n_structural]
                };
                if arj.abs() <= pivot_tol {
                    continue;
                }
                // Moving z_j by t moves the leaving variable by -arj * t. A column at
                // its lower bound can only increase, one at its upper only decrease;
                // keep those that push the leaving variable the way it needs to go.
                let usable = match at {
                    At::Lower => (arj < 0.0) == rising,
                    At::Upper => (arj > 0.0) == rising,
                    At::Zero => true,
                };
                if !usable {
                    continue;
                }
                let d = self.reduced_cost(j, false);
                let ratio = (d.abs() / arj.abs()).max(0.0);
                let pivot = arj.abs();
                // Ties go to the largest pivot, which is numerically safest. Under
                // Bland's they go to the lowest index instead: the rule only
                // guarantees termination when it governs both choices, which is the
                // same mistake that made strong branching cycle earlier.
                let wins_tie = if bland {
                    best.is_none_or(|(prev, _)| j < prev)
                } else {
                    pivot > best_pivot
                };
                if ratio < best_ratio - tol || (ratio <= best_ratio + tol && wins_tie) {
                    if ratio < best_ratio {
                        best_ratio = ratio;
                    }
                    best = Some((j, arj));
                    best_pivot = pivot;
                }
            }

            let Some((entering, _)) = best else {
                // The dual is unbounded, so the primal has no feasible point.
                return Some(LpStatus::Infeasible);
            };

            // Pivot the chosen column into row r.
            let mut col = std::mem::take(&mut self.col);
            self.lp.column_into(entering, &mut col);
            let mut alpha = std::mem::take(&mut self.alpha);
            self.basis.ftran(&col, &mut alpha);
            self.col = col;
            self.alpha = alpha;

            if self.alpha[r].abs() <= pivot_tol {
                // The pivot row and column disagree about this element's magnitude,
                // which means the inverse has drifted. Rebuild and try again.
                if self.refactorize().is_err() {
                    return Some(LpStatus::Infeasible);
                }
                *iterations += 1;
                continue;
            }

            let leaving = self.basic[r];
            let target = self.lp.value_at(leaving, violated);
            // Move the entering variable exactly far enough to place the leaving one
            // on the bound it was violating.
            let step = (target - self.z[leaving]) / -self.alpha[r];
            stalled = if step.abs() <= 1e-12 { stalled + 1 } else { 0 };
            self.z[entering] += step;
            for i in 0..self.lp.m {
                let bj = self.basic[i];
                self.z[bj] -= self.alpha[i] * step;
            }
            self.z[leaving] = target;

            self.basis.update(&self.alpha, r);
            self.basic[r] = entering;
            self.status[entering] = Status::Basic { row: r };
            self.status[leaving] = Status::NonBasic(violated);

            *iterations += 1;
            if self.basis.updates() >= self.lp.tol.refactor_interval && self.refactorize().is_err()
            {
                return Some(LpStatus::Infeasible);
            }
        }
        None
    }

    /// See [`Lp::gomory_cuts`].
    fn gomory_cuts(&mut self, max_cuts: usize) -> Vec<(Vec<(usize, f64)>, f64)> {
        // A fractionality this close to an integer makes the cut's coefficients blow
        // up: they carry 1/f0 and 1/(1 - f0) factors.
        const MIN_FRACTIONALITY: f64 = 0.01;
        // Cuts whose coefficients span more than this are numerically worthless and
        // actively harmful to the LPs that must then carry them.
        const MAX_DYNAMISM: f64 = 1e6;
        // How many columns a cut may touch: a fraction of the model, but never
        // fewer than a floor.
        //
        // GMI cuts come out far denser than the rows they are derived from -- on
        // MIPLIB's pk1 the model's rows average 23% of the columns while its cuts
        // reach 81%. A dense row destroys the sparsity the LU factorization depends
        // on, so every later solve pays for it: pk1 went from 171k nodes without
        // cuts to 24 with them, each node burning tens of thousands of iterations.
        //
        // The floor matters as much as the fraction. On a sixteen-column model no
        // useful GMI cut is under 40% dense, and a pure fraction rejects every one
        // of them -- while the LPs there are so cheap that a dense row costs nothing
        // worth guarding against.
        const MAX_DENSITY: f64 = 0.4;
        const MIN_SUPPORT: usize = 30;
        const TINY: f64 = 1e-11;

        let lp = self.lp;
        let mut rho: Vec<f64> = Vec::new();
        let mut cuts = Vec::new();
        // Row-major view of A, built once. Substituting the logicals out needs the
        // rows of A, and the matrix is stored by column.
        let rows_of_a = lp.matrix.to_csr();

        for i in 0..lp.m {
            if cuts.len() >= max_cuts {
                break;
            }
            // Only an integer column carries a requirement to exploit.
            let basic = self.basic[i];
            if basic >= lp.n_structural || !lp.integer[basic] {
                continue;
            }
            let beta = self.z[basic];
            let f0 = beta - beta.floor();
            if !(MIN_FRACTIONALITY..=1.0 - MIN_FRACTIONALITY).contains(&f0) {
                continue;
            }

            self.basis.btran_unit(i, &mut rho);

            // Accumulate the cut over z = [x; s], then fold the logicals into x.
            let mut z_coeff = vec![0.0f64; lp.n_total()];
            let mut rhs = 1.0f64;
            let mut usable = true;

            for j in 0..lp.n_total() {
                let Status::NonBasic(at) = self.status[j] else {
                    continue;
                };
                if lp.lower[j] == lp.upper[j] {
                    // Fixed: contributes a constant, and w_j is identically zero.
                    continue;
                }
                let alpha = if j < lp.n_structural {
                    let (rows, vals) = lp.matrix.column(j);
                    rows.iter()
                        .zip(vals)
                        .map(|(&r, &v)| rho[r] * v)
                        .sum::<f64>()
                } else {
                    -rho[j - lp.n_structural]
                };
                // At an upper bound the shift w_j = u_j - z_j reverses the sign.
                let a = match at {
                    At::Upper => -alpha,
                    _ => alpha,
                };
                if a.abs() <= TINY {
                    continue;
                }

                // The two coefficient formulas differ, and picking by position
                // rather than by declared integrality would silently produce an
                // invalid cut on any model with continuous columns.
                let c = if j < lp.n_structural && lp.integer[j] {
                    // Integer column: the piecewise-linear integer coefficient.
                    let f = a - a.floor();
                    if f <= f0 {
                        f / f0
                    } else {
                        (1.0 - f) / (1.0 - f0)
                    }
                } else {
                    // Continuous column.
                    if a >= 0.0 { a / f0 } else { -a / (1.0 - f0) }
                };
                if c.abs() <= TINY {
                    continue;
                }

                // Undo the shift: w_j = z_j - l_j, or u_j - z_j at an upper bound.
                match at {
                    At::Upper => {
                        let u = lp.upper[j];
                        if !u.is_finite() {
                            usable = false;
                            break;
                        }
                        z_coeff[j] -= c;
                        rhs -= c * u;
                    }
                    _ => {
                        let l = lp.lower[j];
                        if !l.is_finite() {
                            usable = false;
                            break;
                        }
                        z_coeff[j] += c;
                        rhs += c * l;
                    }
                }
            }
            if !usable {
                continue;
            }

            // Eliminate the logicals with s_k = sum_j A[k][j] x_j, walking row k of A
            // directly. Searching every column for row k instead turns this into
            // O(columns * nonzeros-per-column) per logical, which on a dense 256x256
            // model was most of the solve time.
            let mut x_coeff = z_coeff[..lp.n_structural].to_vec();
            for k in 0..lp.m {
                let h = z_coeff[lp.n_structural + k];
                if h == 0.0 {
                    continue;
                }
                let (columns, values) = rows_of_a.column(k);
                for (&j, &v) in columns.iter().zip(values) {
                    x_coeff[j] += h * v;
                }
            }

            let coefficients: Vec<(usize, f64)> = x_coeff
                .iter()
                .enumerate()
                .filter(|&(_, &v)| v.abs() > 1e-9)
                .map(|(j, &v)| (j, v))
                .collect();
            if coefficients.is_empty() {
                continue;
            }
            let largest = coefficients
                .iter()
                .map(|&(_, v)| v.abs())
                .fold(0.0f64, f64::max);
            let smallest = coefficients
                .iter()
                .map(|&(_, v)| v.abs())
                .fold(f64::INFINITY, f64::min);
            let allowed = MIN_SUPPORT.max((MAX_DENSITY * lp.n_structural as f64) as usize);
            if largest / smallest > MAX_DYNAMISM || !rhs.is_finite() || coefficients.len() > allowed
            {
                continue;
            }

            cuts.push((coefficients, rhs));
        }
        cuts
    }

    /// Finish, handing the final basis and its factorization back for caching.
    fn done(self, status: LpStatus, iterations: usize) -> (LpSolution, Vec<usize>, Basis) {
        let solution = self.finish(status, iterations);
        (solution, self.basic, self.basis)
    }

    fn finish(&self, status: LpStatus, iterations: usize) -> LpSolution {
        let x = self.z[..self.lp.n_structural].to_vec();
        let objective = if status == LpStatus::Optimal {
            self.objective()
        } else {
            f64::NAN
        };
        LpSolution {
            status,
            objective,
            x,
            iterations,
            basis: BasisState {
                basic: self.basic.clone(),
                status: self.status.clone(),
            },
        }
    }

    fn run(
        mut self,
        max_iterations: usize,
        cutoff: Option<f64>,
    ) -> (LpSolution, Vec<usize>, Basis) {
        // Degenerate (zero-length) steps in a row before switching to Bland's rule.
        const STALL_LIMIT: usize = 100;

        let mut iterations = 0usize;
        let mut stalled = 0usize;

        if !self.factorized && self.refactorize().is_err() {
            return self.done(LpStatus::Infeasible, 0);
        }

        // A warm start inherits dual feasibility from its parent, so the dual method
        // repairs the one bound the branch just changed in a few pivots. When that
        // does not hold -- a cold start, or a basis from an unrelated problem -- fall
        // through to the primal method, which needs no assumptions.
        let entry = if self.entry_hint == Entry::Dual && self.is_dual_feasible() {
            Entry::Dual
        } else {
            Entry::Primal
        };
        if entry == Entry::Dual {
            match self.run_dual(max_iterations, cutoff, &mut iterations) {
                Some(status) => return self.done(status, iterations),
                // Ran out of iterations inside the dual method.
                None => return self.done(LpStatus::IterationLimit, iterations),
            }
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
                return self.done(status, iterations);
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
                    return self.done(status, iterations);
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
            if self.basis.updates() >= self.lp.tol.refactor_interval && self.refactorize().is_err()
            {
                return self.done(LpStatus::Infeasible, iterations);
            }
        }

        self.done(LpStatus::IterationLimit, iterations)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{RowSense, Sense};
    use crate::sparse::SparseMatrix;

    /// Build a problem from dense rows, for readable test cases.
    pub(super) fn problem(obj: &[f64], rows: &[(&[f64], RowSense, f64)]) -> Problem {
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
            col_type: vec![crate::model::VarType::Integer; n],
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
    fn warm_start_agrees_with_a_cold_solve_on_every_branch() {
        // The differential test the dual simplex actually needs: for each column,
        // fix it to 0 and to 1 in turn, then check that re-solving from the parent
        // basis reaches the same objective as solving that subproblem from scratch.
        // This is exactly the operation branch-and-bound performs at every node.
        let p = problem(
            &[3.0, 5.0, 2.0, 7.0],
            &[
                (&[2.0, 3.0, 1.0, 4.0], RowSense::Ge, 5.0),
                (&[1.0, 1.0, 1.0, 1.0], RowSense::Le, 3.0),
                (&[4.0, 1.0, 2.0, 1.0], RowSense::Ge, 3.0),
            ],
        );
        let mut lp = Lp::relaxation(&p);
        let root = lp.solve();
        assert_eq!(root.status, LpStatus::Optimal);

        for j in 0..p.n_cols() {
            for fix in [0.0, 1.0] {
                let saved = lp.column_bounds(j);
                lp.set_column_bounds(j, fix, fix);

                let warm = lp.solve_warm(&root.basis, None, 10_000);
                let cold = lp.solve();
                assert_eq!(warm.status, cold.status, "column {j} fixed to {fix}");
                if cold.status == LpStatus::Optimal {
                    assert!(
                        (warm.objective - cold.objective).abs() < 1e-7,
                        "column {j} fixed to {fix}: warm {} vs cold {}",
                        warm.objective,
                        cold.objective
                    );
                    assert_feasible(&p, &warm.x);
                }

                lp.set_column_bounds(j, saved.0, saved.1);
            }
        }
    }

    #[test]
    fn warm_start_detects_an_infeasible_branch() {
        // Fixing both columns to 0 cannot satisfy the row, and the dual simplex must
        // report that rather than returning a bogus optimum.
        let p = problem(&[1.0, 1.0], &[(&[1.0, 1.0], RowSense::Ge, 1.0)]);
        let mut lp = Lp::relaxation(&p);
        let root = lp.solve();
        lp.set_column_bounds(0, 0.0, 0.0);
        lp.set_column_bounds(1, 0.0, 0.0);
        let warm = lp.solve_warm(&root.basis, None, 10_000);
        assert_eq!(warm.status, LpStatus::Infeasible);
    }

    #[test]
    fn cutoff_abandons_a_node_whose_bound_is_already_too_weak() {
        // Forcing the expensive column in pushes the bound to 5; a cutoff below that
        // must stop the solve rather than complete it.
        let p = problem(&[1.0, 5.0], &[(&[1.0, 1.0], RowSense::Ge, 1.0)]);
        let mut lp = Lp::relaxation(&p);
        let root = lp.solve();
        lp.set_column_bounds(0, 0.0, 0.0);

        let cut = lp.solve_warm(&root.basis, Some(2.0), 10_000);
        assert_eq!(cut.status, LpStatus::CutOff);

        // With a cutoff above the true bound the same solve completes normally.
        let full = lp.solve_warm(&root.basis, Some(10.0), 10_000);
        assert_eq!(full.status, LpStatus::Optimal);
        assert!((full.objective - 5.0).abs() < 1e-9, "{}", full.objective);
    }

    #[test]
    fn warm_start_costs_fewer_iterations_than_a_cold_solve() {
        // The whole point of warm starting. Not a tight bound -- just that inheriting
        // the parent basis avoids repeating the work.
        let rows: Vec<(&[f64], RowSense, f64)> = vec![
            (&[2.0, 3.0, 1.0, 4.0, 1.0, 2.0][..], RowSense::Ge, 6.0),
            (&[1.0, 2.0, 3.0, 1.0, 2.0, 1.0][..], RowSense::Ge, 5.0),
            (&[3.0, 1.0, 2.0, 2.0, 1.0, 3.0][..], RowSense::Ge, 7.0),
            (&[1.0, 1.0, 1.0, 1.0, 1.0, 1.0][..], RowSense::Le, 4.0),
        ];
        let p = problem(&[4.0, 6.0, 3.0, 8.0, 2.0, 5.0], &rows);
        let mut lp = Lp::relaxation(&p);
        let root = lp.solve();

        lp.set_column_bounds(3, 1.0, 1.0);
        let warm = lp.solve_warm(&root.basis, None, 10_000);
        let cold = lp.solve();
        assert_eq!(warm.status, LpStatus::Optimal);
        assert!(
            warm.iterations < cold.iterations,
            "warm {} iterations, cold {}",
            warm.iterations,
            cold.iterations
        );
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

#[cfg(test)]
mod cache_tests {
    use super::tests::problem;
    use super::*;
    use crate::model::RowSense;

    #[test]
    fn a_cached_factorization_does_not_change_the_answer() {
        // Re-solving the same bases repeatedly must give the same results whether
        // the factors came from the cache or were rebuilt. Alternating between two
        // unrelated bases and back exercises eviction and reuse together.
        let p = problem(
            &[3.0, 5.0, 2.0, 7.0, 1.0],
            &[
                (&[2.0, 3.0, 1.0, 4.0, 1.0], RowSense::Ge, 5.0),
                (&[1.0, 1.0, 1.0, 1.0, 1.0], RowSense::Le, 3.0),
                (&[4.0, 1.0, 2.0, 1.0, 3.0], RowSense::Ge, 3.0),
            ],
        );
        let mut lp = Lp::relaxation(&p);
        let root = lp.solve();

        // Two different subproblems, each solved twice, interleaved.
        let mut results = Vec::new();
        for round in 0..3 {
            for column in [0usize, 3] {
                let saved = lp.column_bounds(column);
                lp.set_column_bounds(column, 1.0, 1.0);
                let solved = lp.solve_warm(&root.basis, None, 10_000);
                lp.set_column_bounds(column, saved.0, saved.1);
                if round == 0 {
                    results.push((solved.status, solved.objective));
                } else {
                    let (status, objective) = results[usize::from(column == 3)];
                    assert_eq!(solved.status, status, "column {column}, round {round}");
                    assert!(
                        (solved.objective - objective).abs() < 1e-9,
                        "column {column}, round {round}: {} vs {objective}",
                        solved.objective
                    );
                }
            }
        }
    }

    #[test]
    fn a_cloned_lp_starts_with_no_cache() {
        // Each worker in a parallel search clones the LP; carrying a cache across
        // would share nothing useful and only cost memory.
        let p = problem(&[1.0, 1.0], &[(&[1.0, 1.0], RowSense::Ge, 1.0)]);
        let mut lp = Lp::relaxation(&p);
        let _ = lp.solve();
        assert!(!lp.factors.is_empty(), "solving should populate the cache");
        assert!(
            lp.clone().factors.is_empty(),
            "a clone carried the cache over"
        );
    }

    #[test]
    fn invalidation_empties_the_cache() {
        let p = problem(&[1.0, 1.0], &[(&[1.0, 1.0], RowSense::Ge, 1.0)]);
        let mut lp = Lp::relaxation(&p);
        let _ = lp.solve();
        lp.invalidate_factors();
        assert!(lp.factors.is_empty());
        // And it still solves correctly afterwards.
        assert_eq!(lp.solve().status, LpStatus::Optimal);
    }
}
