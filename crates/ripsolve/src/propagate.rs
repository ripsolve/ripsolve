//! Domain propagation at a node.
//!
//! Until this existed, a node's only inference was its LP relaxation: `Worker::process`
//! applied the branch's bound changes and solved. Everything the model *says* about those
//! bounds -- that a row cannot reach its right-hand side any more, that a column set to
//! one forbids every column it conflicts with -- was left for the simplex to rediscover
//! as arithmetic, or not discovered at all.
//!
//! That is the gap the other solvers' own statistics point at. On MIPLIB's
//! `chromaticindex32-8` a reference solver reduces the model by nothing, separates no
//! cuts, and closes it in 0.71 seconds where this search spends sixty and 3959 nodes with
//! the optimum already in hand. What it carries that this did not is a clique table it
//! propagates at every node -- and the conflict graph here is the same size, 12642 edges
//! against its 8438 entries. The information was already built and only ever used to
//! separate cuts, seed heuristics and probe at the root.
//!
//! Two rules, run together to a fixpoint:
//!
//! - **Conflicts.** A binary fixed to a value excludes every literal the conflict graph
//!   joins it to, which forces the columns behind those literals to their other value.
//!   One clique of a set partitioning row decides its whole row from a single branch.
//! - **Row activity.** A row's least and greatest attainable values, given the bounds in
//!   force, bound every column in it; on an integer column the implied bound rounds
//!   inwards, which is what turns a fractional implication into a fixing.
//!
//! Both are ordinary and neither is interesting on its own. What matters is that they run
//! *per node* against bounds the branch has just changed, and that they are bounded: this
//! runs on every node of the tree, so it is given a work budget in matrix entries read and
//! stops when it is spent rather than when it is finished. Stopping early is sound -- the
//! bounds it has proved are still implied, it has simply not proved all of them.

use crate::cuts::{Conflicts, literal};
use crate::model::{Problem, VarType};
use crate::sparse::SparseMatrix;

/// The node's bounds cannot all hold at once, so nothing below it can be feasible.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Infeasible;

/// Scratch reused across nodes, so propagating costs no allocation.
///
/// A worker keeps one of these for the life of the search. The queue and the flag vector
/// are cleared per call rather than rebuilt.
pub struct Scratch {
    /// Columns whose bounds have moved and whose consequences are not drawn yet.
    queue: Vec<u32>,
    /// Whether a column is already waiting, so it is queued once however often it moves.
    queued: Vec<bool>,
}

impl Scratch {
    pub fn new(n_cols: usize) -> Self {
        Scratch {
            queue: Vec::new(),
            queued: vec![false; n_cols],
        }
    }
}

/// Is `j` a binary column under the bounds it started with?
fn is_binary(problem: &Problem, j: usize) -> bool {
    problem.col_type[j] == VarType::Integer
        && problem.col_lb[j] >= -0.5
        && problem.col_ub[j] <= 1.5
}

/// Tighten `lb` and `ub` to a fixpoint from the columns in `seeds`.
///
/// `lb` and `ub` are this node's bounds, already carrying the branch. They are tightened
/// in place. Returns how many columns became fixed, which is what the caller reports;
/// `Err(Infeasible)` means the node is empty and needs no LP at all.
///
/// `budget` is in matrix entries read. It bounds the whole call, not one round.
pub fn propagate(
    problem: &Problem,
    rows_of_a: &SparseMatrix,
    conflicts: &Conflicts,
    scratch: &mut Scratch,
    lb: &mut [f64],
    ub: &mut [f64],
    seeds: &[u32],
    budget: usize,
    tolerance: f64,
) -> Result<usize, Infeasible> {
    scratch.queue.clear();
    for &j in seeds {
        if !scratch.queued[j as usize] {
            scratch.queued[j as usize] = true;
            scratch.queue.push(j);
        }
    }
    let mut fixed = 0usize;
    let mut work = 0usize;
    let mut head = 0usize;

    // A local closure cannot hold `scratch` and the bounds at once, so pushing is written
    // out where it is needed rather than shared.
    while head < scratch.queue.len() {
        if work >= budget {
            break;
        }
        let j = scratch.queue[head] as usize;
        head += 1;
        scratch.queued[j] = false;

        if lb[j] > ub[j] + tolerance {
            return Err(Infeasible);
        }

        // Conflicts, when this column is fixed and binary.
        if is_binary(problem, j) && (ub[j] - lb[j]).abs() <= tolerance {
            let value = lb[j] > 0.5;
            let node = literal(j, value);
            let mut contradiction = false;
            conflicts.for_each_adjacent(node, |excluded| {
                let k = (excluded / 2) as usize;
                // An even literal is "column k is one"; being excluded forces it to zero.
                let forced = if excluded.is_multiple_of(2) { 0.0 } else { 1.0 };
                if forced > ub[k] + tolerance || forced < lb[k] - tolerance {
                    contradiction = true;
                    return false;
                }
                if (ub[k] - lb[k]).abs() > tolerance {
                    lb[k] = forced;
                    ub[k] = forced;
                    fixed += 1;
                    if !scratch.queued[k] {
                        scratch.queued[k] = true;
                        scratch.queue.push(k as u32);
                    }
                }
                true
            });
            if contradiction {
                return Err(Infeasible);
            }
            work += conflicts.degree(node);
        }

        // Row activity, over the rows this column appears in.
        let (rows, _) = problem.matrix.column(j);
        for &i in rows {
            if work >= budget {
                break;
            }
            let (columns, values) = rows_of_a.column(i);
            work += columns.len();
            let (row_lb, row_ub) = (problem.row_lb[i], problem.row_ub[i]);
            if !row_lb.is_finite() && !row_ub.is_finite() {
                continue;
            }

            // Least and greatest the row can reach under the bounds in force. A column
            // with an infinite bound makes one of them infinite, and an infinite
            // activity implies nothing about anything, so those rows are left alone.
            let (mut min_activity, mut max_activity) = (0.0f64, 0.0f64);
            let mut finite = true;
            for (&k, &a) in columns.iter().zip(values) {
                let (low, high) = if a >= 0.0 {
                    (a * lb[k], a * ub[k])
                } else {
                    (a * ub[k], a * lb[k])
                };
                if !low.is_finite() || !high.is_finite() {
                    finite = false;
                    break;
                }
                min_activity += low;
                max_activity += high;
            }
            if !finite {
                continue;
            }
            if min_activity > row_ub + tolerance || max_activity < row_lb - tolerance {
                return Err(Infeasible);
            }

            for (&k, &a) in columns.iter().zip(values) {
                if a == 0.0 {
                    continue;
                }
                let (low, high) = if a >= 0.0 {
                    (a * lb[k], a * ub[k])
                } else {
                    (a * ub[k], a * lb[k])
                };
                // What the rest of the row leaves for this term.
                let head_room = row_ub - (min_activity - low);
                let floor_room = row_lb - (max_activity - high);
                let (mut new_lb, mut new_ub) = (lb[k], ub[k]);
                if a > 0.0 {
                    if row_ub.is_finite() {
                        new_ub = new_ub.min(head_room / a);
                    }
                    if row_lb.is_finite() {
                        new_lb = new_lb.max(floor_room / a);
                    }
                } else {
                    if row_ub.is_finite() {
                        new_lb = new_lb.max(head_room / a);
                    }
                    if row_lb.is_finite() {
                        new_ub = new_ub.min(floor_room / a);
                    }
                }
                if problem.col_type[k] == VarType::Integer {
                    new_lb = (new_lb - tolerance).ceil();
                    new_ub = (new_ub + tolerance).floor();
                }
                if new_lb > new_ub + tolerance {
                    return Err(Infeasible);
                }
                // Only a real move is worth queueing; a bound that repeats itself would
                // otherwise cycle the worklist forever on a degenerate row.
                let moved = new_lb > lb[k] + tolerance || new_ub < ub[k] - tolerance;
                if moved {
                    let was_free = (ub[k] - lb[k]).abs() > tolerance;
                    lb[k] = lb[k].max(new_lb);
                    ub[k] = ub[k].min(new_ub);
                    if was_free && (ub[k] - lb[k]).abs() <= tolerance {
                        fixed += 1;
                    }
                    if !scratch.queued[k] {
                        scratch.queued[k] = true;
                        scratch.queue.push(k as u32);
                    }
                }
            }
        }
    }

    // Whatever is still waiting is simply not drawn; clear the flags for the next call.
    for &j in &scratch.queue[head..] {
        scratch.queued[j as usize] = false;
    }
    Ok(fixed)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{Builder, RowSense, Sense};

    /// Propagate from a set of branch fixings, returning the bounds it proves.
    fn run(problem: &Problem, fixings: &[(usize, f64)]) -> Result<(Vec<f64>, Vec<f64>), Infeasible> {
        let csr = problem.matrix.to_csr();
        let conflicts = Conflicts::of(problem);
        let mut scratch = Scratch::new(problem.n_cols());
        let mut lb = problem.col_lb.clone();
        let mut ub = problem.col_ub.clone();
        let mut seeds = Vec::new();
        for &(j, v) in fixings {
            lb[j] = v;
            ub[j] = v;
            seeds.push(j as u32);
        }
        propagate(
            problem, &csr, &conflicts, &mut scratch, &mut lb, &mut ub, &seeds, usize::MAX, 1e-9,
        )?;
        Ok((lb, ub))
    }

    #[test]
    fn a_set_partitioning_row_is_decided_by_one_branch() {
        // x + y + z = 1: setting x to one leaves nothing for y and z, which is the
        // reduction a clique carries and the one a colouring model is made of.
        let mut b = Builder::new(Sense::Maximize);
        let x = b.binary("x");
        let y = b.binary("y");
        let z = b.binary("z");
        b.objective(&[(x, 1.0), (y, 1.0), (z, 1.0)]);
        b.row(&[(x, 1.0), (y, 1.0), (z, 1.0)], RowSense::Eq, 1.0);
        let problem = b.build();
        let (lb, ub) = run(&problem, &[(x, 1.0)]).unwrap();
        assert_eq!((lb[y], ub[y]), (0.0, 0.0), "y must go to zero");
        assert_eq!((lb[z], ub[z]), (0.0, 0.0), "z must go to zero");
    }

    #[test]
    fn a_chain_of_implications_runs_to_a_fixpoint() {
        // x=1 forces y=0 through the first row, and y=0 forces z=1 through the second.
        // One round would stop at y; the fixpoint is the point of the worklist.
        let mut b = Builder::new(Sense::Maximize);
        let x = b.binary("x");
        let y = b.binary("y");
        let z = b.binary("z");
        b.objective(&[(x, 1.0), (y, 1.0), (z, 1.0)]);
        b.row(&[(x, 1.0), (y, 1.0)], RowSense::Le, 1.0);
        b.row(&[(y, 1.0), (z, 1.0)], RowSense::Ge, 1.0);
        let problem = b.build();
        let (lb, ub) = run(&problem, &[(x, 1.0)]).unwrap();
        assert_eq!((lb[y], ub[y]), (0.0, 0.0), "y follows from x");
        assert_eq!((lb[z], ub[z]), (1.0, 1.0), "z follows from y, which needs a second round");
    }

    #[test]
    fn a_row_that_cannot_be_satisfied_is_reported_infeasible() {
        // x + y >= 2 with x driven to zero leaves y alone to carry 2, which it cannot.
        let mut b = Builder::new(Sense::Minimize);
        let x = b.binary("x");
        let y = b.binary("y");
        b.objective(&[(x, 1.0), (y, 1.0)]);
        b.row(&[(x, 1.0), (y, 1.0)], RowSense::Ge, 2.0);
        let problem = b.build();
        assert_eq!(run(&problem, &[(x, 0.0)]), Err(Infeasible));
    }

    #[test]
    fn an_integer_bound_rounds_inwards() {
        // 2x <= 3 over an integer x in [0,5] implies x <= 1.5, and therefore x <= 1.
        // Without the rounding this proves nothing a branch would act on.
        let mut b = Builder::new(Sense::Maximize);
        let x = b.column("x", 0.0, 5.0, VarType::Integer);
        let y = b.column("y", 0.0, 5.0, VarType::Integer);
        b.objective(&[(x, 1.0), (y, 1.0)]);
        b.row(&[(x, 2.0), (y, 1.0)], RowSense::Le, 3.0);
        let problem = b.build();
        let (_, ub) = run(&problem, &[(y, 0.0)]).unwrap();
        assert_eq!(ub[x], 1.0, "2x <= 3 makes x at most 1, not 1.5");
    }

    #[test]
    fn the_budget_stops_it_without_making_it_wrong() {
        // A budget of zero draws no consequence at all, which has to be sound: the node
        // is simply not tightened. Anything it *had* proved would still hold.
        let mut b = Builder::new(Sense::Maximize);
        let x = b.binary("x");
        let y = b.binary("y");
        b.objective(&[(x, 1.0), (y, 1.0)]);
        b.row(&[(x, 1.0), (y, 1.0)], RowSense::Le, 1.0);
        let problem = b.build();
        let csr = problem.matrix.to_csr();
        let conflicts = Conflicts::of(&problem);
        let mut scratch = Scratch::new(problem.n_cols());
        let (mut lb, mut ub) = (problem.col_lb.clone(), problem.col_ub.clone());
        lb[x] = 1.0;
        ub[x] = 1.0;
        let fixed = propagate(
            &problem, &csr, &conflicts, &mut scratch, &mut lb, &mut ub, &[x as u32], 0, 1e-9,
        )
        .unwrap();
        assert_eq!(fixed, 0);
        assert_eq!((lb[y], ub[y]), (0.0, 1.0), "y is untouched, not wrongly fixed");
    }

    /// Every point the original admits must survive propagation from a partial assignment
    /// that is consistent with it. This is the property that matters: propagation may
    /// prove bounds, and may prove infeasibility, but it may never remove a feasible
    /// completion of the branch it was given.
    #[test]
    fn propagation_never_removes_a_feasible_point() {
        for seed in 0..400u64 {
            let mut state = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
            let mut next = move || {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                state
            };
            let n = 4 + (next() % 5) as usize;
            let mut b = Builder::new(Sense::Maximize);
            let columns: Vec<usize> = (0..n).map(|j| b.binary(&format!("x{j}"))).collect();
            b.objective(&columns.iter().map(|&j| (j, 1.0)).collect::<Vec<_>>());
            for _ in 0..2 + (next() % 4) {
                let mut terms = Vec::new();
                for &j in &columns {
                    if next() % 2 == 0 {
                        terms.push((j, (next() % 5) as f64 - 2.0));
                    }
                }
                if terms.is_empty() {
                    continue;
                }
                let span: f64 = terms.iter().map(|&(_, a)| a.abs()).sum();
                let rhs = (next() % (2 * span as u64 + 1)) as f64 - span;
                let sense = match next() % 3 {
                    0 => RowSense::Le,
                    1 => RowSense::Ge,
                    _ => RowSense::Eq,
                };
                b.row(&terms, sense, rhs);
            }
            let problem = b.build();

            // Enumerate every binary point, keep the feasible ones.
            let feasible: Vec<Vec<f64>> = (0..1u32 << n)
                .map(|mask| {
                    (0..n)
                        .map(|j| f64::from((mask >> j) & 1))
                        .collect::<Vec<f64>>()
                })
                .filter(|x| crate::heuristic::is_feasible(&problem, x, 1e-9))
                .collect();

            // Branch on one column at each value, and check propagation keeps every
            // feasible point that agrees with the branch.
            for j in 0..n {
                for value in [0.0, 1.0] {
                    let agreeing: Vec<&Vec<f64>> =
                        feasible.iter().filter(|x| x[j] == value).collect();
                    match run(&problem, &[(j, value)]) {
                        Err(Infeasible) => assert!(
                            agreeing.is_empty(),
                            "seed {seed}: propagation called x{j}={value} infeasible, but \
                             {} feasible points agree with it",
                            agreeing.len()
                        ),
                        Ok((lb, ub)) => {
                            for point in agreeing {
                                for k in 0..n {
                                    assert!(
                                        point[k] >= lb[k] - 1e-9 && point[k] <= ub[k] + 1e-9,
                                        "seed {seed}: branch x{j}={value} propagated x{k} to \
                                         [{}, {}], which excludes a feasible point at {}",
                                        lb[k], ub[k], point[k]
                                    );
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
