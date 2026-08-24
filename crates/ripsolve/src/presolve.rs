//! Presolve: tighten the model before the search ever sees it.
//!
//! Every reduction here is applied **in place**. A fixed column becomes
//! `lb == ub` rather than being deleted, and a redundant row has its bounds freed
//! to `(-inf, +inf)` rather than being dropped. Nothing is renumbered, so the
//! solution vector still maps one-to-one onto the original columns and there is no
//! postsolve step at all.
//!
//! That trade is deliberate. Physically removing rows and columns would make the
//! LPs smaller, but it requires an index map and a postsolve pass to rebuild the
//! original solution — historically the most bug-prone part of a presolver. The
//! simplex already skips fixed columns during pricing and a freed row's logical is
//! basic and never binds, so most of the benefit arrives anyway. Physical removal
//! is a later optimization, worth doing when instance sizes make the wasted rows
//! matter.
//!
//! The reductions, applied to a fixpoint:
//!
//! - **Activity bounds.** A row whose minimum activity already exceeds its upper
//!   bound (or maximum falls short of its lower) proves the model infeasible; one
//!   whose whole activity range sits inside its bounds is redundant.
//! - **Forcing rows.** When a row can only be satisfied at the extreme of its
//!   activity range, every column in it is pinned to the value that achieves it.
//! - **Bound propagation.** The residual activity of the other columns implies a
//!   bound on each column, which for a binary rounds to a fixing.
//! - **Coefficient tightening.** A coefficient larger than the row can actually
//!   use is reduced. This never changes the binary feasible set, but it tightens
//!   the LP relaxation, which is the whole point.

use crate::model::Problem;
use crate::sparse::SparseMatrix;

/// What presolve did.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Stats {
    pub rounds: usize,
    pub fixed_columns: usize,
    pub redundant_rows: usize,
    pub tightened_coefficients: usize,
}

/// The result of presolving.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Outcome {
    /// The model was reduced (possibly not at all) and remains possibly feasible.
    Reduced(Stats),
    /// Presolve proved the model has no feasible binary assignment.
    Infeasible,
}

const TOL: f64 = 1e-9;

/// A row's coefficients, as `(column, value)` pairs.
type Row = Vec<(usize, f64)>;

/// Reduce `problem` in place.
pub fn presolve(problem: &mut Problem, max_rounds: usize) -> Outcome {
    let n = problem.n_cols();
    let m = problem.n_rows();
    let mut stats = Stats::default();

    // Work row-wise: every reduction here reasons about one row at a time.
    let csr = problem.matrix.to_csr();
    let mut rows: Vec<Row> = (0..m)
        .map(|i| {
            let (cols, vals) = csr.column(i);
            cols.iter().copied().zip(vals.iter().copied()).collect()
        })
        .collect();

    for round in 0..max_rounds {
        stats.rounds = round + 1;
        let mut changed = false;

        for (i, row) in rows.iter().enumerate() {
            if is_free(problem.row_lb[i], problem.row_ub[i]) {
                continue;
            }
            let (min_activity, max_activity) = activity(row, &problem.col_lb, &problem.col_ub);

            if min_activity > problem.row_ub[i] + TOL || max_activity < problem.row_lb[i] - TOL {
                return Outcome::Infeasible;
            }

            if min_activity >= problem.row_lb[i] - TOL && max_activity <= problem.row_ub[i] + TOL {
                free_row(problem, i);
                stats.redundant_rows += 1;
                changed = true;
                continue;
            }

            // Forcing rows: the row is satisfiable only at one end of its activity
            // range, which pins every column appearing in it.
            if (max_activity - problem.row_lb[i]).abs() <= TOL {
                if !pin_row(problem, row, Extreme::Max, &mut stats) {
                    return Outcome::Infeasible;
                }
                free_row(problem, i);
                stats.redundant_rows += 1;
                changed = true;
                continue;
            }
            if (problem.row_ub[i] - min_activity).abs() <= TOL {
                if !pin_row(problem, row, Extreme::Min, &mut stats) {
                    return Outcome::Infeasible;
                }
                free_row(problem, i);
                stats.redundant_rows += 1;
                changed = true;
                continue;
            }

            match propagate(problem, row, i, min_activity, max_activity, &mut stats) {
                Propagation::Infeasible => return Outcome::Infeasible,
                Propagation::Changed => changed = true,
                Propagation::Unchanged => {}
            }
        }

        for (i, row) in rows.iter_mut().enumerate() {
            if tighten_coefficients(problem, row, i, &mut stats) {
                changed = true;
            }
        }

        if fix_columns_absent_from_every_row(problem, &rows, &mut stats) {
            changed = true;
        }

        if !changed {
            break;
        }
    }

    // Rebuild the column-major matrix from the (possibly retightened) rows.
    let triplets = rows
        .iter()
        .enumerate()
        .flat_map(|(i, row)| row.iter().map(move |&(j, v)| (i, j, v)));
    problem.matrix = SparseMatrix::from_triplets(m, n, triplets);

    debug_assert_eq!(problem.n_cols(), n);
    Outcome::Reduced(stats)
}

fn is_free(lb: f64, ub: f64) -> bool {
    lb == f64::NEG_INFINITY && ub == f64::INFINITY
}

fn free_row(problem: &mut Problem, i: usize) {
    problem.row_lb[i] = f64::NEG_INFINITY;
    problem.row_ub[i] = f64::INFINITY;
}

/// The minimum and maximum a row's activity can reach within the column bounds.
fn activity(row: &Row, col_lb: &[f64], col_ub: &[f64]) -> (f64, f64) {
    let mut min = 0.0;
    let mut max = 0.0;
    for &(j, a) in row {
        let (lo, hi) = (a * col_lb[j], a * col_ub[j]);
        min += lo.min(hi);
        max += lo.max(hi);
    }
    (min, max)
}

/// Which end of its activity range a row is being pinned to.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Extreme {
    Min,
    Max,
}

/// Fix every column of a forcing row to the value that achieves the extreme.
///
/// Returns false if that demands two different values for one column.
fn pin_row(problem: &mut Problem, row: &Row, extreme: Extreme, stats: &mut Stats) -> bool {
    for &(j, a) in row {
        let (lo, hi) = (problem.col_lb[j], problem.col_ub[j]);
        // A positive coefficient reaches its minimum contribution at the column's
        // lower bound; a negative one at its upper.
        let target = match (extreme, a > 0.0) {
            (Extreme::Min, true) | (Extreme::Max, false) => lo,
            (Extreme::Min, false) | (Extreme::Max, true) => hi,
        };
        if !fix_column(problem, j, target, stats) {
            return false;
        }
    }
    true
}

/// Narrow a column to a single value, reporting false if that is impossible.
fn fix_column(problem: &mut Problem, j: usize, value: f64, stats: &mut Stats) -> bool {
    if problem.col_lb[j] > value + TOL || problem.col_ub[j] < value - TOL {
        return false;
    }
    let already_fixed = problem.col_lb[j] == problem.col_ub[j];
    problem.col_lb[j] = value;
    problem.col_ub[j] = value;
    if !already_fixed {
        stats.fixed_columns += 1;
    }
    true
}

enum Propagation {
    Unchanged,
    Changed,
    Infeasible,
}

/// Derive per-column bounds from the residual activity of the rest of the row.
///
/// With `row_lb <= rest + a*x_j <= row_ub`, the term `a*x_j` is confined to
/// `[row_lb - rest_max, row_ub - rest_min]`. Dividing by `a` turns that into a
/// bound on the column itself, which for a binary rounds inwards to an integer and
/// so often fixes it outright.
fn propagate(
    problem: &mut Problem,
    row: &Row,
    i: usize,
    min_activity: f64,
    max_activity: f64,
    stats: &mut Stats,
) -> Propagation {
    let (row_lb, row_ub) = (problem.row_lb[i], problem.row_ub[i]);
    let mut changed = false;

    for &(j, a) in row {
        if a == 0.0 {
            continue;
        }
        let (lo, hi) = (a * problem.col_lb[j], a * problem.col_ub[j]);
        let (own_min, own_max) = (lo.min(hi), lo.max(hi));
        // What the *other* columns of this row can contribute.
        let rest_min = min_activity - own_min;
        let rest_max = max_activity - own_max;

        // An infinite row bound imposes nothing on that side; `inf - finite` is
        // still infinite, so the arithmetic carries that through correctly.
        let contribution_lo = row_lb - rest_max;
        let contribution_hi = row_ub - rest_min;

        let (implied_lo, implied_hi) = if a > 0.0 {
            (contribution_lo / a, contribution_hi / a)
        } else {
            (contribution_hi / a, contribution_lo / a)
        };

        // An integer column's implied bound rounds inwards, which is what turns a
        // weak implication into a fixing. A continuous column must not be rounded:
        // doing so cuts off values it is entitled to take, and the search then
        // proves a worse solution optimal.
        let (implied_lo, implied_hi) = if problem.is_integer(j) {
            (
                ceil_with_tolerance(implied_lo),
                floor_with_tolerance(implied_hi),
            )
        } else {
            (implied_lo, implied_hi)
        };
        let new_lb = problem.col_lb[j].max(implied_lo);
        let new_ub = problem.col_ub[j].min(implied_hi);
        if new_lb > new_ub + TOL {
            return Propagation::Infeasible;
        }
        if new_lb > problem.col_lb[j] + TOL || new_ub < problem.col_ub[j] - TOL {
            let was_fixed = problem.col_lb[j] == problem.col_ub[j];
            problem.col_lb[j] = new_lb;
            problem.col_ub[j] = new_ub;
            if !was_fixed && new_lb == new_ub {
                stats.fixed_columns += 1;
            }
            changed = true;
        }
    }

    if changed {
        Propagation::Changed
    } else {
        Propagation::Unchanged
    }
}

fn ceil_with_tolerance(v: f64) -> f64 {
    if v.is_finite() { (v - TOL).ceil() } else { v }
}

fn floor_with_tolerance(v: f64) -> f64 {
    if v.is_finite() { (v + TOL).floor() } else { v }
}

/// A column appearing in no still-active row is decided entirely by its cost.
fn fix_columns_absent_from_every_row(
    problem: &mut Problem,
    rows: &[Row],
    stats: &mut Stats,
) -> bool {
    let n = problem.n_cols();
    let mut appears = vec![false; n];
    for (i, row) in rows.iter().enumerate() {
        if is_free(problem.row_lb[i], problem.row_ub[i]) {
            continue;
        }
        for &(j, a) in row {
            if a != 0.0 {
                appears[j] = true;
            }
        }
    }

    let mut changed = false;
    for (j, &seen) in appears.iter().enumerate() {
        if seen || problem.col_lb[j] == problem.col_ub[j] {
            continue;
        }
        // Minimization: a non-negative cost prefers the lower bound.
        let target = if problem.obj[j] > 0.0 {
            problem.col_lb[j]
        } else {
            problem.col_ub[j]
        };
        if fix_column(problem, j, target, stats) {
            changed = true;
        }
    }
    changed
}

/// Reduce coefficients that are larger than their row can use.
fn tighten_coefficients(problem: &mut Problem, row: &mut Row, i: usize, stats: &mut Stats) -> bool {
    // Only single-sided rows. Tightening one side of a range row invalidates the
    // reasoning for the other, and an equality row has no slack to exploit.
    let lb = problem.row_lb[i];
    let ub = problem.row_ub[i];
    let ge = lb.is_finite() && ub == f64::INFINITY;
    let le = ub.is_finite() && lb == f64::NEG_INFINITY;
    if !ge && !le {
        return false;
    }

    let mut changed = false;
    for k in 0..row.len() {
        let (j, a) = row[k];
        // Only free binaries. The reduction reasons about a column being switched
        // fully on or fully off, which is meaningless for a continuous column and
        // wrong for a general integer that can land in between.
        if a == 0.0 || !problem.is_binary(j) || problem.col_lb[j] != 0.0 || problem.col_ub[j] != 1.0
        {
            continue;
        }
        let (min_activity, max_activity) = activity(row, &problem.col_lb, &problem.col_ub);

        if ge {
            let slack = problem.row_lb[i] - min_activity;
            if slack <= TOL {
                continue;
            }
            if a > 0.0 && a > slack + TOL {
                // `x_j = 1` already satisfies the row on its own; the excess above
                // `slack` can never be needed, so spending it is pure LP looseness.
                row[k].1 = slack;
                changed = true;
                stats.tightened_coefficients += 1;
            } else if a < 0.0 && -a > slack + TOL {
                // Mirror image: `x_j = 0` is the non-binding case.
                row[k].1 = -slack;
                problem.row_lb[i] = min_activity - a;
                changed = true;
                stats.tightened_coefficients += 1;
            }
        } else {
            let slack = max_activity - problem.row_ub[i];
            if slack <= TOL {
                continue;
            }
            if a > 0.0 && a > slack + TOL {
                row[k].1 = slack;
                problem.row_ub[i] = max_activity - a;
                changed = true;
                stats.tightened_coefficients += 1;
            } else if a < 0.0 && -a > slack + TOL {
                row[k].1 = -slack;
                changed = true;
                stats.tightened_coefficients += 1;
            }
        }
    }
    changed
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::generate::{Kind, Spec};
    use crate::model::{RowSense, Sense};
    use lp_parser_rs::problem::LpProblem;

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
            col_type: vec![crate::model::VarType::Integer; n],
            col_names: (0..n).map(|j| format!("x{j}")).collect(),
            row_names: (0..m).map(|i| format!("c{i}")).collect(),
        }
    }

    /// Every binary assignment the model admits, by exhaustive enumeration.
    ///
    /// Only usable for small `n`, but it is ground truth: presolve is allowed to
    /// change coefficients and bounds however it likes, provided this set — and the
    /// objective value of each member — comes out identical.
    fn feasible_set(p: &Problem) -> Vec<(u32, f64)> {
        let n = p.n_cols();
        assert!(n <= 20, "exhaustive enumeration needs a small model");
        let csr = p.matrix.to_csr();
        let mut out = Vec::new();

        for mask in 0u32..(1u32 << n) {
            let x: Vec<f64> = (0..n).map(|j| f64::from((mask >> j) & 1)).collect();
            let within_bounds =
                (0..n).all(|j| x[j] >= p.col_lb[j] - 1e-9 && x[j] <= p.col_ub[j] + 1e-9);
            if !within_bounds {
                continue;
            }
            let satisfied = (0..p.n_rows()).all(|i| {
                let (cols, vals) = csr.column(i);
                let activity: f64 = cols.iter().zip(vals).map(|(&j, &a)| a * x[j]).sum();
                activity >= p.row_lb[i] - 1e-9 && activity <= p.row_ub[i] + 1e-9
            });
            if satisfied {
                let obj: f64 = p.obj.iter().zip(&x).map(|(c, v)| c * v).sum();
                out.push((mask, obj));
            }
        }
        out
    }

    /// Presolve's actual contract, checked exhaustively.
    ///
    /// It is *not* that the feasible set is preserved — fixing a column that appears
    /// in no constraint to its cheaper value legitimately discards feasible points.
    /// The contract is two-sided:
    ///
    /// 1. **Sound**: every point the reduced model admits was admitted by the
    ///    original, at the same cost. Presolve may never invent solutions.
    /// 2. **Optimality-preserving**: the best achievable objective is unchanged, so
    ///    no optimal solution was reduced away.
    fn assert_presolve_is_sound(before: &Problem, label: &str) -> Stats {
        let expected = feasible_set(before);
        let mut after = before.clone();

        let stats = match presolve(&mut after, 20) {
            Outcome::Infeasible => {
                assert!(
                    expected.is_empty(),
                    "{label}: presolve rejected a model with {} feasible points",
                    expected.len()
                );
                return Stats::default();
            }
            Outcome::Reduced(stats) => stats,
        };

        after
            .validate()
            .unwrap_or_else(|e| panic!("{label}: invalid after presolve: {e}"));
        let got = feasible_set(&after);

        for (mask, obj) in &got {
            let original = expected.iter().find(|(m, _)| m == mask);
            match original {
                None => panic!("{label}: presolve admits {mask:0b}, which the original rejects"),
                Some((_, want)) => assert!(
                    (want - obj).abs() < 1e-9,
                    "{label}: point {mask:0b} costs {obj} after presolve, {want} before"
                ),
            }
        }

        let best = |set: &[(u32, f64)]| set.iter().map(|(_, o)| *o).fold(f64::INFINITY, f64::min);
        let (before_best, after_best) = (best(&expected), best(&got));
        assert!(
            (before_best - after_best).abs() < 1e-9
                || (before_best.is_infinite() && after_best.is_infinite()),
            "{label}: optimum moved from {before_best} to {after_best}"
        );
        stats
    }

    /// Coefficient tightening on its own must preserve the feasible set *exactly*.
    ///
    /// Unlike the pipeline as a whole, this reduction is only ever allowed to change
    /// the relaxation, never the integer solutions — so it gets the strict check the
    /// full presolve cannot be held to.
    fn assert_tightening_preserves_feasible_set(before: &Problem, label: &str) {
        let expected = feasible_set(before);
        let mut after = before.clone();

        let csr = after.matrix.to_csr();
        let mut rows: Vec<Row> = (0..after.n_rows())
            .map(|i| {
                let (cols, vals) = csr.column(i);
                cols.iter().copied().zip(vals.iter().copied()).collect()
            })
            .collect();
        let mut stats = Stats::default();
        for (i, row) in rows.iter_mut().enumerate() {
            tighten_coefficients(&mut after, row, i, &mut stats);
        }
        let triplets = rows
            .iter()
            .enumerate()
            .flat_map(|(i, row)| row.iter().map(move |&(j, v)| (i, j, v)));
        after.matrix = SparseMatrix::from_triplets(after.n_rows(), after.n_cols(), triplets);

        assert_eq!(
            expected,
            feasible_set(&after),
            "{label}: coefficient tightening changed the feasible set"
        );
    }

    #[test]
    fn is_sound_on_generated_instances() {
        // The broad net: presolve is allowed to do anything to these models except
        // change which binary points satisfy them.
        for kind in [Kind::Knapsack, Kind::Covering, Kind::Signed] {
            for seed in 0..12u64 {
                let spec = Spec {
                    kind,
                    n_cols: 12,
                    n_rows: 6,
                    seed,
                };
                let parsed = LpProblem::parse(&spec.to_lp()).unwrap();
                let p = Problem::from_lp(&parsed).unwrap();
                assert_presolve_is_sound(&p, &spec.name());
            }
        }
    }

    /// The coefficient of column `j` in row `i`, read back out of the matrix.
    fn coefficient(p: &Problem, i: usize, j: usize) -> f64 {
        let csr = p.matrix.to_csr();
        let (cols, vals) = csr.column(i);
        cols.iter()
            .zip(vals)
            .find(|(c, _)| **c == j)
            .map(|(_, v)| *v)
            .unwrap_or(0.0)
    }

    #[test]
    fn tightens_an_oversized_coefficient_on_a_ge_row() {
        // `5x0 + 2x1 + 2x2 + 2x3 >= 3`: x0 alone already satisfies the row, so 5 is
        // more than the row can use and reduces to 3. The binary solutions are
        // unchanged; the relaxation is not, since x0 = 0.6 was feasible and is not now.
        //
        // The spare columns are load-bearing: with fewer of them, bound propagation
        // fixes x0 = 1 outright before tightening gets a chance -- a stronger
        // reduction, but not the one under test.
        let p = problem(
            &[1.0, 1.0, 1.0, 1.0],
            &[(&[5.0, 2.0, 2.0, 2.0], RowSense::Ge, 3.0)],
        );
        let stats = assert_presolve_is_sound(&p, "ge tightening");
        assert_eq!(stats.tightened_coefficients, 1);
        assert_eq!(
            stats.fixed_columns, 0,
            "propagation should not have decided anything"
        );

        let mut after = p.clone();
        presolve(&mut after, 20);
        assert!(
            (coefficient(&after, 0, 0) - 3.0).abs() < 1e-9,
            "{}",
            coefficient(&after, 0, 0)
        );
        assert!(
            (after.row_lb[0] - 3.0).abs() < 1e-9,
            "row bound moved: {}",
            after.row_lb[0]
        );
        assert_tightening_preserves_feasible_set(&p, "ge tightening");
    }

    #[test]
    fn tightens_an_oversized_coefficient_on_a_le_row() {
        // `5x0 + 2x1 + 2x2 + 2x3 <= 8`: the row can only ever be violated by 3 above
        // its bound, so x0's coefficient reduces to 3 and the bound follows to 6.
        let p = problem(
            &[1.0, 1.0, 1.0, 1.0],
            &[(&[5.0, 2.0, 2.0, 2.0], RowSense::Le, 8.0)],
        );
        let stats = assert_presolve_is_sound(&p, "le tightening");
        assert_eq!(stats.tightened_coefficients, 1);

        let mut after = p.clone();
        presolve(&mut after, 20);
        assert!(
            (coefficient(&after, 0, 0) - 3.0).abs() < 1e-9,
            "{}",
            coefficient(&after, 0, 0)
        );
        assert!(
            (after.row_ub[0] - 6.0).abs() < 1e-9,
            "row bound is {}",
            after.row_ub[0]
        );
        assert_tightening_preserves_feasible_set(&p, "le tightening");
    }

    #[test]
    fn tightening_strengthens_the_relaxation() {
        // The point of the reduction: the same integer solutions, a higher LP bound.
        use crate::lp::{Lp, LpStatus};
        let p = problem(
            &[1.0, 1.0, 1.0, 1.0],
            &[(&[5.0, 2.0, 2.0, 2.0], RowSense::Ge, 3.0)],
        );
        let before = Lp::relaxation(&p).solve();
        let mut after = p.clone();
        presolve(&mut after, 20);
        let reduced = Lp::relaxation(&after).solve();

        assert_eq!(before.status, LpStatus::Optimal);
        assert_eq!(reduced.status, LpStatus::Optimal);
        assert!(
            reduced.objective > before.objective + 1e-9,
            "relaxation did not improve: {} -> {}",
            before.objective,
            reduced.objective
        );
    }

    #[test]
    fn tightening_handles_negative_coefficients() {
        for sense in [RowSense::Ge, RowSense::Le] {
            for rhs in [-4.0, -1.0, 0.0, 1.0, 3.0] {
                let p = problem(&[1.0, 2.0, 1.0], &[(&[-6.0, 3.0, -2.0], sense, rhs)]);
                assert_presolve_is_sound(&p, &format!("{sense:?} rhs {rhs}"));
            }
        }
    }

    #[test]
    fn detects_an_infeasible_row() {
        // Maximum activity is 3, so the row can never reach 5.
        let p = problem(&[1.0, 1.0], &[(&[2.0, 1.0], RowSense::Ge, 5.0)]);
        let mut after = p.clone();
        assert_eq!(presolve(&mut after, 20), Outcome::Infeasible);
    }

    #[test]
    fn removes_a_redundant_row() {
        // Maximum activity is 3, so `<= 10` can never bind.
        let p = problem(&[1.0, 1.0], &[(&[2.0, 1.0], RowSense::Le, 10.0)]);
        let stats = assert_presolve_is_sound(&p, "redundant");
        assert_eq!(stats.redundant_rows, 1);
    }

    #[test]
    fn a_forcing_row_pins_every_column() {
        // Only x0 = x1 = 1 reaches 3, so both are forced.
        let p = problem(&[1.0, 1.0], &[(&[2.0, 1.0], RowSense::Ge, 3.0)]);
        let mut after = p.clone();
        assert!(matches!(presolve(&mut after, 20), Outcome::Reduced(_)));
        assert_eq!((after.col_lb[0], after.col_ub[0]), (1.0, 1.0));
        assert_eq!((after.col_lb[1], after.col_ub[1]), (1.0, 1.0));
        assert_presolve_is_sound(&p, "forcing");
    }

    #[test]
    fn a_singleton_row_fixes_its_column() {
        let p = problem(&[1.0, 1.0], &[(&[3.0, 0.0], RowSense::Ge, 2.0)]);
        let mut after = p.clone();
        assert!(matches!(presolve(&mut after, 20), Outcome::Reduced(_)));
        assert_eq!((after.col_lb[0], after.col_ub[0]), (1.0, 1.0));
        assert_presolve_is_sound(&p, "singleton");
    }

    #[test]
    fn a_column_in_no_row_follows_its_cost() {
        // x1 appears nowhere, so cost decides: positive goes to 0, negative to 1.
        let p = problem(&[1.0, 4.0], &[(&[1.0, 0.0], RowSense::Ge, 1.0)]);
        let mut after = p.clone();
        presolve(&mut after, 20);
        assert_eq!((after.col_lb[1], after.col_ub[1]), (0.0, 0.0));

        let p = problem(&[1.0, -4.0], &[(&[1.0, 0.0], RowSense::Ge, 1.0)]);
        let mut after = p.clone();
        presolve(&mut after, 20);
        assert_eq!((after.col_lb[1], after.col_ub[1]), (1.0, 1.0));
    }

    #[test]
    fn presolve_is_idempotent() {
        // A second pass over an already-reduced model must find nothing new.
        let spec = Spec {
            kind: Kind::Knapsack,
            n_cols: 14,
            n_rows: 7,
            seed: 4,
        };
        let mut p = Problem::from_lp(&LpProblem::parse(&spec.to_lp()).unwrap()).unwrap();
        presolve(&mut p, 20);
        let once = p.clone();
        let stats = presolve(&mut p, 20);
        if let Outcome::Reduced(stats) = stats {
            assert_eq!(stats.fixed_columns, 0);
            assert_eq!(stats.tightened_coefficients, 0);
            assert_eq!(stats.redundant_rows, 0);
        }
        assert_eq!(once.col_lb, p.col_lb);
        assert_eq!(once.col_ub, p.col_ub);
        assert_eq!(once.matrix, p.matrix);
    }
}
