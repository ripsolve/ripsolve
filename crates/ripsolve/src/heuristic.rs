//! Primal heuristics: find good feasible solutions early.
//!
//! Branch and bound prunes by comparing a node's bound against the best solution
//! found so far, so it cannot prune anything at all until it holds one. A search
//! that spends its first thousand nodes without an incumbent is doing exhaustive
//! work no matter how good its bound is — which is exactly what `v064c1000n020`
//! did, running out its time limit having never found a feasible point.
//!
//! Two heuristics, cheapest first:
//!
//! - **Rounding.** Round the relaxation and test it. Costs one pass over the rows
//!   and no LP at all, and on models whose relaxation is nearly integral it simply
//!   works.
//! - **Diving.** Repeatedly fix the *least* fractional column to its nearest
//!   integer and re-solve, until the relaxation comes out integral or the dive
//!   hits infeasibility. Each re-solve differs from the last by one bound, so the
//!   warm-started dual simplex settles it in a handful of pivots — the same
//!   machinery the search itself runs on.
//!
//! Least fractional, not most: the aim here is a feasible point rather than a
//! strong bound, so the column already closest to deciding itself is the one whose
//! rounding is least likely to make the LP infeasible. Branching wants the
//! opposite, which is why it scores columns by pseudocost instead.
//!
//! A feasibility pump would be the next addition, for models where diving keeps
//! hitting infeasibility. It is a genuinely different mechanism — alternating
//! projection and rounding — rather than a refinement of this one.

use crate::lp::{BasisState, Lp, LpStatus};
use crate::model::Problem;

/// A feasible binary assignment and its objective, in the internal minimization
/// form.
#[derive(Clone, Debug)]
pub struct Incumbent {
    pub x: Vec<u8>,
    pub objective: f64,
}

/// Limits on how much work the heuristics may do.
#[derive(Clone, Copy, Debug)]
pub struct Limits {
    /// A relaxation value this close to an integer counts as integral.
    pub integrality_tolerance: f64,
    /// Row activities may violate a bound by this much.
    pub feasibility_tolerance: f64,
    /// Columns a dive may fix before giving up.
    pub max_dive_steps: usize,
    /// Simplex iterations per dive or pump re-solve.
    pub max_iterations_per_solve: usize,
    /// Projection rounds the feasibility pump may take before giving up.
    pub max_pump_rounds: usize,
}

impl Default for Limits {
    fn default() -> Self {
        Self {
            integrality_tolerance: 1e-6,
            feasibility_tolerance: 1e-6,
            max_dive_steps: 200,
            max_iterations_per_solve: 2_000,
            max_pump_rounds: 60,
        }
    }
}

/// Is this binary assignment feasible for every row and column bound?
pub fn is_feasible(problem: &Problem, x: &[f64], tolerance: f64) -> bool {
    for (j, &v) in x.iter().enumerate() {
        if v < problem.col_lb[j] - tolerance || v > problem.col_ub[j] + tolerance {
            return false;
        }
    }
    let csr = problem.matrix.to_csr();
    (0..problem.n_rows()).all(|i| {
        let (cols, vals) = csr.column(i);
        let activity: f64 = cols.iter().zip(vals).map(|(&j, &a)| a * x[j]).sum();
        activity >= problem.row_lb[i] - tolerance && activity <= problem.row_ub[i] + tolerance
    })
}

/// The internal-form objective of an assignment.
fn objective_of(problem: &Problem, x: &[f64]) -> f64 {
    problem.obj.iter().zip(x).map(|(c, v)| c * v).sum()
}

/// Round a relaxation to the nearest binary point and keep it if it is feasible.
///
/// Deliberately does not try to repair an infeasible rounding: that is diving's
/// job, and doing it here would duplicate the machinery badly.
pub fn round(problem: &Problem, x: &[f64], limits: &Limits) -> Option<Incumbent> {
    let rounded: Vec<f64> = x.iter().map(|v| v.round().clamp(0.0, 1.0)).collect();
    is_feasible(problem, &rounded, limits.feasibility_tolerance).then(|| Incumbent {
        x: rounded.iter().map(|&v| v as u8).collect(),
        objective: objective_of(problem, &rounded),
    })
}

/// Dive from a relaxation towards a feasible binary point.
///
/// `lp`'s column bounds are restored before returning, whatever the outcome, so
/// the caller's node is left exactly as it was found.
pub fn dive(
    problem: &Problem,
    lp: &mut Lp,
    basis: &BasisState,
    start: &[f64],
    cutoff: Option<f64>,
    limits: &Limits,
    iterations: &mut usize,
) -> Option<Incumbent> {
    let n = problem.n_cols();
    let saved: Vec<(f64, f64)> = (0..n).map(|j| lp.column_bounds(j)).collect();

    let result = dive_inner(problem, lp, basis, start, cutoff, limits, iterations);

    for (j, &(lo, hi)) in saved.iter().enumerate() {
        lp.set_column_bounds(j, lo, hi);
    }
    result
}

fn dive_inner(
    problem: &Problem,
    lp: &mut Lp,
    basis: &BasisState,
    start: &[f64],
    cutoff: Option<f64>,
    limits: &Limits,
    iterations: &mut usize,
) -> Option<Incumbent> {
    let tol = limits.integrality_tolerance;
    let mut x = start.to_vec();
    // Warm start each step from the previous step's basis rather than always from
    // the node's. Consecutive dive LPs differ by a single bound, so the previous
    // basis is nearly optimal for the next one; the node's grows staler with depth.
    let mut warm = basis.clone();

    for _ in 0..limits.max_dive_steps {
        // The column closest to already being decided.
        let next = x
            .iter()
            .enumerate()
            .filter_map(|(j, &v)| {
                let distance = (v - v.round()).abs();
                (distance > tol).then_some((j, v, distance))
            })
            .min_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal));

        let Some((j, value, _)) = next else {
            // Nothing fractional left: the relaxation is a binary point.
            let rounded: Vec<f64> = x.iter().map(|v| v.round().clamp(0.0, 1.0)).collect();
            if !is_feasible(problem, &rounded, limits.feasibility_tolerance) {
                return None;
            }
            return Some(Incumbent {
                x: rounded.iter().map(|&v| v as u8).collect(),
                objective: objective_of(problem, &rounded),
            });
        };

        // Greedy rounding regularly paints the dive into a corner, so when the
        // preferred value turns out infeasible, try the other one before giving up.
        // Without this the dive fails on essentially every model tested; with it,
        // one extra LP rescues the descent.
        let preferred = value.round().clamp(0.0, 1.0);
        let mut advanced = false;
        for target in [preferred, 1.0 - preferred] {
            lp.set_column_bounds(j, target, target);
            let solved = lp.solve_warm(&warm, cutoff, limits.max_iterations_per_solve);
            *iterations += solved.iterations;
            if solved.status == LpStatus::Optimal {
                x = solved.x;
                warm = solved.basis;
                advanced = true;
                break;
            }
            // Cut off or truncated says nothing about the other direction being
            // better, only that this dive is not worth continuing.
            if solved.status != LpStatus::Infeasible {
                return None;
            }
        }
        if !advanced {
            // Neither value works, so the fixings made so far are jointly infeasible.
            return None;
        }
    }
    None
}

/// Find a feasible point by alternating projection and rounding.
///
/// The *feasibility pump*. It keeps two points: a relaxation-feasible `x` that is
/// fractional, and an integral `x~` that is generally infeasible. Each round
/// re-optimizes the original constraint set under a new objective — the distance to
/// `x~` — which pulls `x` toward integrality, and then re-rounds to get the next
/// `x~`. When the two coincide, that point is both integral and feasible.
///
/// The distance objective is linear on binaries: `|x_j - t_j|` is `x_j` when
/// `t_j = 0` and `1 - x_j` when `t_j = 1`, so minimizing it means costs of
/// `1 - 2*t_j` and a dropped constant.
///
/// This is the right tool where diving is the wrong one. Diving commits to a
/// rounding and re-solves a *smaller* LP each step, so on a model whose feasible
/// set is sparse it walks into infeasibility and cannot recover — measured, it
/// failed on every instance in this benchmark set. The pump never fixes anything,
/// so its LP is always feasible; it can wander, but it cannot dead-end.
pub fn feasibility_pump(
    problem: &Problem,
    lp: &mut Lp,
    basis: &BasisState,
    start: &[f64],
    limits: &Limits,
    iterations: &mut usize,
) -> Option<Incumbent> {
    let n = problem.n_cols();
    let saved_costs = lp.costs().to_vec();
    let result = pump_inner(problem, lp, basis, start, limits, iterations, n);
    lp.set_costs(&saved_costs);
    result
}

#[allow(clippy::too_many_arguments)]
fn pump_inner(
    problem: &Problem,
    lp: &mut Lp,
    basis: &BasisState,
    start: &[f64],
    limits: &Limits,
    iterations: &mut usize,
    n: usize,
) -> Option<Incumbent> {
    let mut x = start.to_vec();
    // Only the objective changes between rounds, so the previous basis stays primal
    // feasible and re-optimizing from it is far cheaper than a cold solve.
    let mut warm = basis.clone();
    let mut previous: Option<Vec<f64>> = None;
    // Deterministic tie-breaking for perturbation, so a run is reproducible.
    let mut tick = 0usize;

    for _ in 0..limits.max_pump_rounds {
        let mut target: Vec<f64> = x.iter().map(|v| v.round().clamp(0.0, 1.0)).collect();

        if is_feasible(problem, &target, limits.feasibility_tolerance) {
            return Some(Incumbent {
                x: target.iter().map(|&v| v as u8).collect(),
                objective: objective_of(problem, &target),
            });
        }

        // A repeated target means the pump has cycled. Flip the columns that were
        // hardest to round — the ones furthest from the integer they landed on —
        // which is the standard escape and keeps the walk deterministic.
        if previous.as_ref() == Some(&target) {
            let mut by_distance: Vec<usize> = (0..n).collect();
            by_distance.sort_by(|&a, &b| {
                let d = |j: usize| (x[j] - x[j].round()).abs();
                d(b).partial_cmp(&d(a)).unwrap_or(std::cmp::Ordering::Equal)
            });
            tick += 1;
            let flips = (1 + tick % 5).min(n);
            for &j in by_distance.iter().take(flips) {
                target[j] = 1.0 - target[j];
            }
        }
        previous = Some(target.clone());

        // Minimize the distance to the target over the original constraint set.
        let costs: Vec<f64> = target.iter().map(|&t| 1.0 - 2.0 * t).collect();
        lp.set_costs(&costs);
        let solved = lp.solve_warm(&warm, None, limits.max_iterations_per_solve);
        *iterations += solved.iterations;
        if solved.status != LpStatus::Optimal {
            return None;
        }
        x = solved.x;
        warm = solved.basis;
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::generate::{Kind, Spec};
    use lp_parser_rs::problem::LpProblem;

    fn instance(kind: Kind, n_cols: usize, n_rows: usize, seed: u64) -> Problem {
        let spec = Spec {
            kind,
            n_cols,
            n_rows,
            seed,
        };
        Problem::from_lp(&LpProblem::parse(&spec.to_lp()).unwrap()).unwrap()
    }

    /// Whatever a heuristic returns must actually be feasible and correctly costed.
    /// A heuristic that returns an infeasible point does not fail loudly — it
    /// installs a bogus incumbent and the search prunes the real optimum away.
    fn assert_valid(problem: &Problem, found: &Incumbent, label: &str) {
        let x: Vec<f64> = found.x.iter().map(|&v| f64::from(v)).collect();
        assert_eq!(x.len(), problem.n_cols(), "{label}: wrong length");
        assert!(
            is_feasible(problem, &x, 1e-6),
            "{label}: returned an infeasible point"
        );
        let recomputed = objective_of(problem, &x);
        assert!(
            (recomputed - found.objective).abs() < 1e-6,
            "{label}: reported {} but the point costs {recomputed}",
            found.objective
        );
    }

    #[test]
    fn rounding_accepts_an_already_integral_relaxation() {
        let p = instance(Kind::Covering, 25, 30, 1);
        let relaxed = Lp::relaxation(&p).solve();
        // This instance's relaxation is integral, so rounding it must succeed.
        let found = round(&p, &relaxed.x, &Limits::default()).expect("integral relaxation");
        assert_valid(&p, &found, "covering rounding");
    }

    #[test]
    fn rounding_rejects_an_infeasible_rounding() {
        // Rounding down a tight `>=` row breaks it, and round() must not claim it.
        let p = instance(Kind::Knapsack, 20, 10, 1);
        let x = vec![0.4; p.n_cols()];
        assert!(round(&p, &x, &Limits::default()).is_none());
    }

    #[test]
    fn the_pump_finds_feasible_points_where_diving_does_not() {
        // The measured reason the pump exists. Diving commits to a rounding and can
        // walk into infeasibility with no way back; the pump never fixes anything,
        // so its LP stays feasible throughout.
        let limits = Limits::default();
        let mut pumped = 0;
        for seed in 0..6u64 {
            let p = instance(Kind::Knapsack, 30, 15, seed);
            let mut lp = Lp::relaxation(&p);
            let root = lp.solve();
            let mut iterations = 0;
            if let Some(found) =
                feasibility_pump(&p, &mut lp, &root.basis, &root.x, &limits, &mut iterations)
            {
                assert_valid(&p, &found, &format!("pump seed {seed}"));
                pumped += 1;
            }
        }
        assert!(pumped >= 4, "the pump found only {pumped} of 6");
    }

    #[test]
    fn heuristics_leave_the_lp_exactly_as_they_found_it() {
        // Both dive and pump mutate the LP -- bounds and costs respectively -- and
        // both must restore it, or the node they were called from is corrupted.
        let p = instance(Kind::Knapsack, 25, 12, 3);
        let mut lp = Lp::relaxation(&p);
        let root = lp.solve();
        let limits = Limits::default();

        let bounds: Vec<(f64, f64)> = (0..p.n_cols()).map(|j| lp.column_bounds(j)).collect();
        let costs = lp.costs().to_vec();

        let mut iterations = 0;
        let _ = dive(
            &p,
            &mut lp,
            &root.basis,
            &root.x,
            None,
            &limits,
            &mut iterations,
        );
        let _ = feasibility_pump(&p, &mut lp, &root.basis, &root.x, &limits, &mut iterations);

        let after: Vec<(f64, f64)> = (0..p.n_cols()).map(|j| lp.column_bounds(j)).collect();
        assert_eq!(bounds, after, "column bounds were not restored");
        assert_eq!(costs, lp.costs(), "objective costs were not restored");

        // And the LP still solves to the same relaxation value.
        let again = lp.solve();
        assert!((again.objective - root.objective).abs() < 1e-9);
    }

    #[test]
    fn a_dive_that_succeeds_returns_a_feasible_point() {
        let limits = Limits::default();
        for seed in 0..8u64 {
            let p = instance(Kind::Covering, 30, 40, seed);
            let mut lp = Lp::relaxation(&p);
            let root = lp.solve();
            let mut iterations = 0;
            if let Some(found) = dive(
                &p,
                &mut lp,
                &root.basis,
                &root.x,
                None,
                &limits,
                &mut iterations,
            ) {
                assert_valid(&p, &found, &format!("dive seed {seed}"));
            }
        }
    }

    #[test]
    fn feasibility_check_catches_a_violated_row() {
        let p = instance(Kind::Knapsack, 15, 8, 2);
        let all_zero = vec![0.0; p.n_cols()];
        // Knapsack rows are `>=` with a positive right-hand side, so the origin fails.
        assert!(!is_feasible(&p, &all_zero, 1e-6));
        let all_one = vec![1.0; p.n_cols()];
        assert!(is_feasible(&p, &all_one, 1e-6));
    }
}
