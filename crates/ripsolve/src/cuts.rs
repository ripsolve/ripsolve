//! Cutting planes: inequalities valid for every binary solution, but violated by
//! the current fractional relaxation.
//!
//! Presolve reduces the *model*; cuts strengthen the *relaxation*. On the dense
//! random instances this solver is aimed at, presolve finds nothing at all — and
//! neither does a leading commercial solver's — so cuts are where the remaining
//! bound gap lives.
//!
//! # Knapsack covers
//!
//! Take a row in knapsack form, `sum w_j y_j <= c` with every `w_j > 0` and every
//! `y_j` binary. A *cover* is a set `C` whose weights overrun the capacity,
//! `sum_{j in C} w_j > c`. No binary solution can switch on all of `C` at once, so
//!
//! ```text
//!     sum_{j in C} y_j <= |C| - 1
//! ```
//!
//! holds for every feasible point. It is useful when the relaxation violates it,
//! which happens exactly when the fractional values across `C` sum to more than
//! `|C| - 1`.
//!
//! # Getting rows into that form
//!
//! Real rows are neither one-sided nor sign-uniform, so each is normalized first:
//!
//! - A `>=` side becomes `<=` by negating the row.
//! - A negative coefficient is made positive by complementing its column,
//!   `x_j = 1 - y_j`, which moves `a_j` into the capacity.
//!
//! A range or equality row yields a knapsack from each of its two finite sides.
//! Cuts are found in the complemented `y` space and mapped back to `x` at the end,
//! where a complemented column contributes `-1` and shifts the right-hand side.

use crate::model::Problem;

/// A generated inequality, in range form so it drops straight into a [`Problem`].
#[derive(Clone, Debug, PartialEq)]
pub struct Cut {
    /// `(column, coefficient)` pairs, sorted by column.
    pub coefficients: Vec<(usize, f64)>,
    pub lb: f64,
    pub ub: f64,
}

impl Cut {
    /// The row's activity at a point.
    pub fn activity(&self, x: &[f64]) -> f64 {
        self.coefficients.iter().map(|&(j, a)| a * x[j]).sum()
    }

    /// How far `x` violates this cut; zero or less when satisfied.
    pub fn violation(&self, x: &[f64]) -> f64 {
        let activity = self.activity(x);
        let over = if self.ub.is_finite() {
            activity - self.ub
        } else {
            f64::NEG_INFINITY
        };
        let under = if self.lb.is_finite() {
            self.lb - activity
        } else {
            f64::NEG_INFINITY
        };
        over.max(under)
    }
}

/// A row rewritten as `sum w_j y_j <= capacity`, with all weights positive.
struct Knapsack {
    /// `(column, weight, complemented)`.
    terms: Vec<(usize, f64, bool)>,
    capacity: f64,
}

const TOL: f64 = 1e-9;
/// A cut must beat this violation to be worth adding; smaller ones are noise.
const MIN_VIOLATION: f64 = 1e-4;

/// Rewrite one side of a row as a knapsack over its *binary* columns, or `None` if
/// it cannot be one.
///
/// `negate` selects the `>=` side, which becomes `<=` under negation.
///
/// A cover argument counts how many columns can be switched on at once, which only
/// means anything for 0/1 columns. A general integer or continuous column in the
/// row is therefore folded into the capacity at its most favourable value, so the
/// resulting cut stays valid whatever that column does. A row with no binary
/// columns yields nothing.
fn to_knapsack(
    problem: &Problem,
    coefficients: &[(usize, f64)],
    bound: f64,
    negate: bool,
) -> Option<Knapsack> {
    let (col_lb, col_ub) = (&problem.col_lb, &problem.col_ub);
    if !bound.is_finite() {
        return None;
    }
    let sign = if negate { -1.0 } else { 1.0 };
    let mut capacity = sign * bound;
    let mut terms = Vec::with_capacity(coefficients.len());

    for &(j, a) in coefficients {
        let a = sign * a;
        if a == 0.0 {
            continue;
        }
        // A fixed column is a constant, not a decision: fold it into the capacity.
        if col_lb[j] == col_ub[j] {
            capacity -= a * col_lb[j];
            continue;
        }
        // A non-binary column is not something a cover can reason about. Charge the
        // capacity for the least it could contribute, so the cut holds for every
        // value it might take. An infinite worst case means no usable bound.
        if !problem.is_binary(j) {
            let least = if a > 0.0 {
                a * col_lb[j]
            } else {
                a * col_ub[j]
            };
            if !least.is_finite() {
                return None;
            }
            capacity -= least;
            continue;
        }
        if a > 0.0 {
            terms.push((j, a, false));
        } else {
            // x_j = 1 - y_j, so a_j*x_j = a_j - a_j*y_j; the constant moves right.
            capacity -= a;
            terms.push((j, -a, true));
        }
    }

    // A non-positive capacity means no subset fits, which is a feasibility question
    // for presolve rather than a cut to separate.
    if terms.is_empty() || capacity <= TOL {
        return None;
    }
    Some(Knapsack { terms, capacity })
}

/// The value of a knapsack term in the complemented `y` space.
fn y_value(knapsack: &Knapsack, i: usize, x: &[f64]) -> f64 {
    let (col, _, complemented) = knapsack.terms[i];
    if complemented { 1.0 - x[col] } else { x[col] }
}

/// Find a violated minimal cover, optionally forced to contain `seed`.
///
/// The separation problem — find the most violated cover — is itself a knapsack
/// problem, so this uses the standard greedy: take terms in increasing order of
/// `(1 - y_j) / w_j`, since a term near 1 adds most violation per unit of capacity
/// consumed, and stop once the weights overrun.
///
/// Seeding matters more than it looks. One greedy pass per row yields one cover;
/// forcing each fractional term into the cover in turn yields a family of them, and
/// several usually separate where the unseeded one does not.
fn find_cover(knapsack: &Knapsack, x: &[f64], seed: Option<usize>) -> Option<Vec<usize>> {
    let cost = |i: usize| (1.0 - y_value(knapsack, i, x)) / knapsack.terms[i].1;
    let mut order: Vec<usize> = (0..knapsack.terms.len()).collect();
    order.sort_by(|&a, &b| {
        cost(a)
            .partial_cmp(&cost(b))
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let mut chosen = Vec::new();
    let mut weight = 0.0;
    if let Some(seed) = seed {
        chosen.push(seed);
        weight += knapsack.terms[seed].1;
    }
    if weight <= knapsack.capacity + TOL {
        for &i in &order {
            if Some(i) == seed {
                continue;
            }
            chosen.push(i);
            weight += knapsack.terms[i].1;
            if weight > knapsack.capacity + TOL {
                break;
            }
        }
    }
    if weight <= knapsack.capacity + TOL {
        return None;
    }

    // Minimalize: drop the heaviest members the cover can spare, keeping any seed.
    // A minimal cover dominates every cover containing it.
    let mut by_weight: Vec<usize> = chosen.clone();
    by_weight.sort_by(|&a, &b| {
        knapsack.terms[b]
            .1
            .partial_cmp(&knapsack.terms[a].1)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    for i in by_weight {
        if Some(i) == seed || chosen.len() <= 2 {
            continue;
        }
        let w = knapsack.terms[i].1;
        if weight - w > knapsack.capacity + TOL {
            weight -= w;
            chosen.retain(|&k| k != i);
        }
    }

    let total: f64 = chosen.iter().map(|&i| y_value(knapsack, i, x)).sum();
    if total > chosen.len() as f64 - 1.0 + MIN_VIOLATION {
        chosen.sort_unstable();
        Some(chosen)
    } else {
        None
    }
}

/// Largest capacity a lifting DP will be built over.
const MAX_LIFT_CAPACITY: usize = 20_000;
/// Ceiling on the total work one lifting pass may do, as `terms * capacity` summed
/// over the columns lifted.
///
/// The DP is `O(terms * capacity)` *per column lifted*, and both grow with the
/// model, so on a dense 256-column knapsack the unbounded version reached billions
/// of operations per cover — separation cost 1.3 seconds a round to save a handful
/// of nodes. Lifting is an improvement to a cut that is already valid, so running
/// out of budget just means a weaker cut, never a wrong one.
const MAX_LIFT_WORK: usize = 400_000;
/// Fractional columns used to seed cover searches in a single row.
///
/// Every fractional column was a seed originally, which is quadratic in the row's
/// support and bought very little: the covers found from the most fractional few
/// are the ones that separate.
const MAX_SEEDS_PER_ROW: usize = 4;

/// Sequential up-lifting of the columns outside the cover.
///
/// The plain cover inequality `sum_{j in C} y_j <= |C| - 1` says nothing about
/// columns outside `C`, which makes it weak. Lifting gives each outside column the
/// largest coefficient that keeps the inequality valid:
///
/// ```text
///     a_j = rhs - max { sum a_k y_k : sum w_k y_k <= capacity - w_j }
/// ```
///
/// taken over the columns already in the inequality. That inner problem is a 0/1
/// knapsack, solved here exactly by dynamic programming — which needs integer
/// weights, so lifting is skipped when a row has fractional coefficients. Lifting
/// *sequentially* (each column against the inequality built so far, not against the
/// original cover) is what keeps the result valid.
fn lift(knapsack: &Knapsack, cover: &[usize]) -> Option<(Vec<(usize, f64)>, f64)> {
    let integral = |v: f64| (v - v.round()).abs() <= 1e-9 && v.round() >= 0.0;
    if !knapsack.terms.iter().all(|&(_, w, _)| integral(w)) || !integral(knapsack.capacity) {
        return None;
    }
    let capacity = knapsack.capacity.round() as usize;
    if capacity > MAX_LIFT_CAPACITY {
        return None;
    }

    let rhs = cover.len() as f64 - 1.0;
    // (term index, weight, coefficient) for everything currently in the inequality.
    let mut inequality: Vec<(usize, usize, i64)> = cover
        .iter()
        .map(|&i| (i, knapsack.terms[i].1.round() as usize, 1i64))
        .collect();
    let mut work = 0usize;
    // Reused across columns; reallocating a table per column was most of the cost.
    let mut dp: Vec<i64> = Vec::new();

    // Heaviest first is the conventional lifting order and tends to give the
    // strongest coefficients to the columns that matter most.
    let mut outside: Vec<usize> = (0..knapsack.terms.len())
        .filter(|i| !cover.contains(i))
        .collect();
    outside.sort_by(|&a, &b| {
        knapsack.terms[b]
            .1
            .partial_cmp(&knapsack.terms[a].1)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let mut lifted: Vec<(usize, f64)> = Vec::new();
    for j in outside {
        let wj = knapsack.terms[j].1.round() as usize;
        let coefficient = if wj > capacity {
            // The column cannot be switched on at all, so it may carry the whole
            // right-hand side.
            rhs as i64
        } else {
            let budget = capacity - wj;
            work = work.saturating_add(inequality.len().saturating_mul(budget));
            if work > MAX_LIFT_WORK {
                // Out of budget: stop lifting and keep what has been lifted so far,
                // which is still a valid inequality.
                break;
            }
            dp.clear();
            dp.resize(budget + 1, 0);
            for &(_, w, a) in &inequality {
                if w > budget {
                    continue;
                }
                for cap in (w..=budget).rev() {
                    dp[cap] = dp[cap].max(dp[cap - w] + a);
                }
            }
            rhs as i64 - dp[budget]
        };
        if coefficient > 0 {
            inequality.push((j, wj, coefficient));
            lifted.push((j, coefficient as f64));
        }
    }

    if lifted.is_empty() {
        return None;
    }
    Some((lifted, rhs))
}

/// Turn a cover over a knapsack's `y` space back into a cut over the original `x`.
fn cover_to_cut(knapsack: &Knapsack, cover: &[usize]) -> Cut {
    // sum_{j in C} y_j <= |C| - 1, plus any lifted terms, where y_j is x_j or
    // (1 - x_j) depending on whether the column was complemented.
    let mut rhs = cover.len() as f64 - 1.0;
    let mut terms: Vec<(usize, f64)> = cover.iter().map(|&i| (i, 1.0)).collect();
    if let Some((lifted, _)) = lift(knapsack, cover) {
        terms.extend(lifted);
    }

    let mut coefficients = Vec::with_capacity(terms.len());
    for (i, coefficient) in terms {
        let (col, _, complemented) = knapsack.terms[i];
        if complemented {
            // a*(1 - x_j) contributes -a*x_j and moves a to the right-hand side.
            coefficients.push((col, -coefficient));
            rhs -= coefficient;
        } else {
            coefficients.push((col, coefficient));
        }
    }
    coefficients.sort_by_key(|&(j, _)| j);
    Cut {
        coefficients,
        lb: f64::NEG_INFINITY,
        ub: rhs,
    }
}

/// Separate Gomory mixed-integer cuts from the tableau at `basis`.
///
/// Covers are combinatorial and need a row that reads as a knapsack; GMI cuts come
/// from the simplex tableau and need no structure at all, so they reach models the
/// cover separator cannot see. The derivation lives in [`crate::lp::Lp::gomory_cuts`],
/// which has the tableau; this wraps the result and applies the same
/// worth-adding test as every other family.
pub fn separate_gomory(
    lp: &crate::lp::Lp,
    basis: &crate::lp::BasisState,
    x: &[f64],
    limit: usize,
) -> Vec<Cut> {
    let mut found: Vec<Cut> = lp
        .gomory_cuts(basis, limit * 4)
        .into_iter()
        .map(|(coefficients, lb)| Cut {
            coefficients,
            lb,
            ub: f64::INFINITY,
        })
        .filter(|cut| cut.violation(x) > MIN_VIOLATION)
        .collect();

    found.sort_by(|a, b| {
        b.violation(x)
            .partial_cmp(&a.violation(x))
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    found.truncate(limit);
    found
}

/// Separate knapsack cover cuts violated by `x`.
///
/// Returns at most `limit` cuts, strongest violation first.
pub fn separate(problem: &Problem, x: &[f64], limit: usize) -> Vec<Cut> {
    separate_until(problem, x, limit, None)
}

/// As [`separate`], stopping early if `deadline` passes.
pub fn separate_until(
    problem: &Problem,
    x: &[f64],
    limit: usize,
    deadline: Option<std::time::Instant>,
) -> Vec<Cut> {
    let csr = problem.matrix.to_csr();
    let mut found: Vec<Cut> = Vec::new();

    for i in 0..problem.n_rows() {
        // Separation over a large model is itself a long operation; stop rather
        // than run past a budget the caller set.
        if i.is_multiple_of(256) && deadline.is_some_and(|d| std::time::Instant::now() >= d) {
            break;
        }
        let (cols, vals) = csr.column(i);
        if cols.is_empty() {
            continue;
        }
        let coefficients: Vec<(usize, f64)> =
            cols.iter().copied().zip(vals.iter().copied()).collect();

        // Each finite side of the row is its own knapsack.
        for (bound, negate) in [(problem.row_ub[i], false), (problem.row_lb[i], true)] {
            let Some(knapsack) = to_knapsack(problem, &coefficients, bound, negate) else {
                continue;
            };
            // One unseeded pass, then one seeded by each fractional term. Distinct
            // seeds usually produce distinct covers, and only some of them separate.
            // Seeded by the most fractional columns only. Seeding from every
            // fractional column is quadratic in the row's support and finds little
            // the top few do not.
            let mut fractional: Vec<usize> = (0..knapsack.terms.len())
                .filter(|&t| {
                    let v = y_value(&knapsack, t, x);
                    v > 1e-6 && v < 1.0 - 1e-6
                })
                .collect();
            fractional.sort_by(|&a, &b| {
                let d = |t: usize| (y_value(&knapsack, t, x) - 0.5).abs();
                d(a).partial_cmp(&d(b)).unwrap_or(std::cmp::Ordering::Equal)
            });
            fractional.truncate(MAX_SEEDS_PER_ROW);
            let seeds = std::iter::once(None).chain(fractional.into_iter().map(Some));
            for seed in seeds {
                if let Some(cover) = find_cover(&knapsack, x, seed) {
                    let cut = cover_to_cut(&knapsack, &cover);
                    if cut.violation(x) > MIN_VIOLATION {
                        found.push(cut);
                    }
                }
            }
        }
    }

    found.sort_by(|a, b| {
        b.violation(x)
            .partial_cmp(&a.violation(x))
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    found.dedup_by(|a, b| a.coefficients == b.coefficients && a.ub == b.ub);
    found.truncate(limit);
    found
}

impl Problem {
    /// Append cuts as new rows.
    ///
    /// Cuts add rows but never columns, so the solution vector's meaning is
    /// unchanged and nothing downstream needs remapping.
    pub fn add_cuts(&mut self, cuts: &[Cut]) {
        if cuts.is_empty() {
            return;
        }
        use crate::sparse::SparseMatrix;

        let n = self.n_cols();
        let old_rows = self.n_rows();
        let csr = self.matrix.to_csr();

        let mut triplets: Vec<(usize, usize, f64)> = Vec::with_capacity(self.matrix.nnz());
        for i in 0..old_rows {
            let (cols, vals) = csr.column(i);
            triplets.extend(cols.iter().zip(vals).map(|(&j, &v)| (i, j, v)));
        }
        for (k, cut) in cuts.iter().enumerate() {
            let row = old_rows + k;
            triplets.extend(cut.coefficients.iter().map(|&(j, a)| (row, j, a)));
            self.row_lb.push(cut.lb);
            self.row_ub.push(cut.ub);
            self.row_names.push(format!("cut{k}"));
        }

        self.matrix = SparseMatrix::from_triplets(old_rows + cuts.len(), n, triplets);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::generate::{Kind, Spec};
    use crate::lp::{Lp, LpStatus};
    use crate::model::{RowSense, Sense};
    use crate::sparse::SparseMatrix;
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

    /// Every binary point the model admits, by exhaustive enumeration.
    fn feasible_points(p: &Problem) -> Vec<Vec<f64>> {
        let n = p.n_cols();
        assert!(n <= 18);
        let csr = p.matrix.to_csr();
        (0u32..(1u32 << n))
            .map(|mask| {
                (0..n)
                    .map(|j| f64::from((mask >> j) & 1))
                    .collect::<Vec<f64>>()
            })
            .filter(|x| {
                (0..n).all(|j| x[j] >= p.col_lb[j] - 1e-9 && x[j] <= p.col_ub[j] + 1e-9)
                    && (0..p.n_rows()).all(|i| {
                        let (cols, vals) = csr.column(i);
                        let a: f64 = cols.iter().zip(vals).map(|(&j, &v)| v * x[j]).sum();
                        a >= p.row_lb[i] - 1e-9 && a <= p.row_ub[i] + 1e-9
                    })
            })
            .collect()
    }

    /// The property a cut generator must never break: no feasible binary point may
    /// be cut off. A violated cut is not merely suboptimal, it silently produces the
    /// wrong answer, so this is checked exhaustively rather than by spot-check.
    fn assert_cuts_are_valid(p: &Problem, label: &str) -> usize {
        let relaxed = Lp::relaxation(p).solve();
        if relaxed.status != LpStatus::Optimal {
            return 0;
        }
        let cuts = separate(p, &relaxed.x, 64);
        let points = feasible_points(p);

        for cut in &cuts {
            for point in &points {
                assert!(
                    cut.violation(point) <= 1e-9,
                    "{label}: cut {:?} <= {} removes feasible point {point:?} (activity {})",
                    cut.coefficients,
                    cut.ub,
                    cut.activity(point)
                );
            }
            // A cut that does not cut is wasted work.
            assert!(
                cut.violation(&relaxed.x) > MIN_VIOLATION,
                "{label}: cut does not separate the relaxation"
            );
        }
        cuts.len()
    }

    #[test]
    fn lifting_strengthens_a_cover_without_breaking_it() {
        // Lifting is the step most likely to produce a subtly invalid cut, so it gets
        // its own exhaustive check across a spread of capacities.
        for capacity in [4.0, 5.0, 7.0, 9.0, 11.0] {
            let p = problem(
                &[-3.0, -2.0, -2.0, -1.0],
                &[(&[4.0, 3.0, 3.0, 2.0], RowSense::Le, capacity)],
            );
            assert_cuts_are_valid(&p, &format!("lifting cap {capacity}"));
        }
    }

    #[test]
    fn lifting_is_skipped_for_fractional_weights() {
        // The lifting DP indexes by integer weight, so a row it cannot represent must
        // fall back to the unlifted cover rather than round and produce a bad cut.
        let p = problem(
            &[-1.0, -1.0, -1.0],
            &[(&[2.5, 2.5, 2.5], RowSense::Le, 4.2)],
        );
        assert_cuts_are_valid(&p, "fractional weights");
    }

    /// GMI cuts get the same exhaustive treatment as covers: no feasible binary
    /// point may be cut off. They are derived from floating-point tableau entries
    /// rather than from combinatorial reasoning, so they are the family most likely
    /// to produce a subtly invalid inequality.
    fn assert_gomory_is_valid(p: &Problem, label: &str) -> usize {
        let mut lp = Lp::relaxation(p);
        let relaxed = lp.solve();
        if relaxed.status != LpStatus::Optimal {
            return 0;
        }
        let cuts = separate_gomory(&lp, &relaxed.basis, &relaxed.x, 32);
        let points = feasible_points(p);
        for cut in &cuts {
            for point in &points {
                assert!(
                    cut.violation(point) <= 1e-6,
                    "{label}: gomory cut {:?} >= {} removes feasible point {point:?}",
                    cut.coefficients,
                    cut.lb
                );
            }
            assert!(
                cut.violation(&relaxed.x) > MIN_VIOLATION,
                "{label}: gomory cut does not separate the relaxation"
            );
        }
        cuts.len()
    }

    #[test]
    fn gomory_cuts_never_remove_a_feasible_point() {
        let mut total = 0;
        for kind in [Kind::Knapsack, Kind::Covering, Kind::Signed] {
            for seed in 0..14u64 {
                let spec = Spec {
                    kind,
                    n_cols: 13,
                    n_rows: 7,
                    seed,
                };
                let p = Problem::from_lp(&LpProblem::parse(&spec.to_lp()).unwrap()).unwrap();
                total += assert_gomory_is_valid(&p, &spec.name());
            }
        }
        assert!(
            total > 0,
            "no gomory cuts generated anywhere, so nothing was tested"
        );
    }

    #[test]
    fn gomory_cuts_tighten_the_relaxation() {
        let spec = Spec {
            kind: Kind::Knapsack,
            n_cols: 16,
            n_rows: 8,
            seed: 5,
        };
        let mut p = Problem::from_lp(&LpProblem::parse(&spec.to_lp()).unwrap()).unwrap();
        let mut lp = Lp::relaxation(&p);
        let before = lp.solve();
        let cuts = separate_gomory(&lp, &before.basis, &before.x, 16);
        assert!(
            !cuts.is_empty(),
            "expected gomory cuts on a fractional relaxation"
        );

        p.add_cuts(&cuts);
        let after = Lp::relaxation(&p).solve();
        assert_eq!(after.status, LpStatus::Optimal);
        assert!(
            after.objective > before.objective + 1e-9,
            "bound did not improve: {} -> {}",
            before.objective,
            after.objective
        );
    }

    #[test]
    fn cuts_never_remove_a_feasible_point() {
        // The broad net, across all three instance families and both row senses.
        let mut total = 0;
        for kind in [Kind::Knapsack, Kind::Covering, Kind::Signed] {
            for seed in 0..14u64 {
                let spec = Spec {
                    kind,
                    n_cols: 14,
                    n_rows: 7,
                    seed,
                };
                let p = Problem::from_lp(&LpProblem::parse(&spec.to_lp()).unwrap()).unwrap();
                total += assert_cuts_are_valid(&p, &spec.name());
            }
        }
        assert!(
            total > 0,
            "no cuts were generated anywhere, so nothing was tested"
        );
    }

    #[test]
    fn separates_a_textbook_cover() {
        // 3x0 + 3x1 + 3x2 <= 5, maximizing the count. Any two columns overrun the
        // capacity, so x_a + x_b <= 1 is valid for each pair, while the relaxation
        // can afford 5/3 columns in total and so violates it.
        //
        // The weights must not divide the capacity evenly: with 6/5/5 <= 10 the
        // relaxation lands on the integral point [0, 1, 1] and there is correctly
        // nothing to separate.
        let p = problem(
            &[-1.0, -1.0, -1.0],
            &[(&[3.0, 3.0, 3.0], RowSense::Le, 5.0)],
        );
        let relaxed = Lp::relaxation(&p).solve();
        assert_eq!(relaxed.status, LpStatus::Optimal);
        let cuts = separate(&p, &relaxed.x, 8);
        assert!(!cuts.is_empty(), "no cover found for {:?}", relaxed.x);

        // The cover is {a, b} giving `x_a + x_b <= 1`, and lifting then brings the
        // third column in at coefficient 1 as well: no two of the three fit, so
        // `x0 + x1 + x2 <= 1` is valid and strictly stronger than the cover alone.
        let cut = &cuts[0];
        assert_eq!(cut.coefficients.len(), 3, "{:?}", cut.coefficients);
        assert!(
            cut.coefficients.iter().all(|&(_, a)| a == 1.0),
            "{:?}",
            cut.coefficients
        );
        assert!((cut.ub - 1.0).abs() < 1e-9, "rhs is {}", cut.ub);
        assert_cuts_are_valid(&p, "textbook cover");
    }

    #[test]
    fn handles_ge_rows_by_negation() {
        // A `>=` row is a knapsack after negating, which also flips every sign.
        let p = problem(
            &[1.0, 1.0, 1.0],
            &[(&[-6.0, -5.0, -5.0], RowSense::Ge, -10.0)],
        );
        assert_cuts_are_valid(&p, "ge row");
    }

    #[test]
    fn handles_negative_coefficients_by_complementing() {
        let p = problem(
            &[1.0, -1.0, 1.0, -1.0],
            &[(&[6.0, -5.0, 5.0, -4.0], RowSense::Le, 3.0)],
        );
        assert_cuts_are_valid(&p, "complemented");
    }

    #[test]
    fn a_cut_tightens_the_relaxation() {
        // The point of the exercise: adding cuts must raise the bound.
        let spec = Spec {
            kind: Kind::Knapsack,
            n_cols: 16,
            n_rows: 8,
            seed: 5,
        };
        let mut p = Problem::from_lp(&LpProblem::parse(&spec.to_lp()).unwrap()).unwrap();
        let before = Lp::relaxation(&p).solve();
        assert_eq!(before.status, LpStatus::Optimal);

        let cuts = separate(&p, &before.x, 32);
        assert!(!cuts.is_empty());
        p.add_cuts(&cuts);
        let after = Lp::relaxation(&p).solve();

        assert_eq!(after.status, LpStatus::Optimal);
        assert!(
            after.objective > before.objective + 1e-9,
            "bound did not improve: {} -> {}",
            before.objective,
            after.objective
        );
    }

    #[test]
    fn adding_cuts_preserves_every_feasible_point() {
        // add_cuts rebuilds the matrix; this checks the rebuild keeps the model's
        // meaning as well as the cuts keeping their validity.
        let spec = Spec {
            kind: Kind::Knapsack,
            n_cols: 14,
            n_rows: 6,
            seed: 2,
        };
        let p = Problem::from_lp(&LpProblem::parse(&spec.to_lp()).unwrap()).unwrap();
        let before = feasible_points(&p);

        let relaxed = Lp::relaxation(&p).solve();
        let cuts = separate(&p, &relaxed.x, 32);
        let mut with_cuts = p.clone();
        with_cuts.add_cuts(&cuts);
        with_cuts.validate().unwrap();

        assert_eq!(
            before,
            feasible_points(&with_cuts),
            "cuts changed the feasible set"
        );
        assert_eq!(with_cuts.n_rows(), p.n_rows() + cuts.len());
        assert_eq!(with_cuts.n_cols(), p.n_cols());
    }

    #[test]
    fn no_cuts_from_an_integral_relaxation() {
        // Nothing to separate when the relaxation is already integral.
        let p = problem(&[1.0, 1.0], &[(&[1.0, 1.0], RowSense::Ge, 1.0)]);
        let relaxed = Lp::relaxation(&p).solve();
        let cuts = separate(&p, &relaxed.x, 8);
        for cut in &cuts {
            assert!(cut.violation(&relaxed.x) > MIN_VIOLATION);
        }
    }
}
