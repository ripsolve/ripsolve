//! Branch and bound over the LP relaxation.
//!
//! Depth-first with warm starts. Each node fixes one column to 0 or 1, which leaves
//! the parent's basis dual feasible, so the child re-solves with the dual simplex in
//! a handful of pivots rather than repeating a cold phase-1. That is what makes an
//! LP per node affordable, and it is the whole reason this design beats enumeration:
//! the relaxation's bound prunes subtrees that implicit enumeration would have to
//! walk.
//!
//! What is deliberately not here yet: presolve, cutting planes, pseudocost
//! branching, and primal heuristics. Those are the milestones that take this from
//! "correct and far better than enumeration" to "competitive", and each slots in
//! without disturbing this loop.

use std::time::{Duration, Instant};

use crate::lp::{BasisState, Lp, LpStatus};
use crate::model::Problem;
use crate::presolve::{self, Outcome};

/// How a search ended.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Status {
    /// The incumbent is proven optimal.
    Optimal,
    /// No binary assignment satisfies the constraints.
    Infeasible,
    /// Stopped at the node limit; any incumbent is best-so-far, not proven.
    NodeLimit,
    /// Stopped at the time limit; any incumbent is best-so-far, not proven.
    TimeLimit,
}

/// Limits and tolerances for a search.
#[derive(Clone, Copy, Debug)]
pub struct Options {
    pub max_nodes: usize,
    pub time_limit: Option<Duration>,
    /// A relaxation value within this of an integer counts as integral.
    pub integrality_tolerance: f64,
    /// Stop once `(incumbent - bound) / |incumbent|` falls to this.
    pub gap_tolerance: f64,
    /// Simplex iteration limit for a single node.
    pub max_iterations_per_node: usize,
    /// Tighten the model before searching.
    pub presolve: bool,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            max_nodes: usize::MAX,
            time_limit: None,
            integrality_tolerance: 1e-6,
            gap_tolerance: 0.0,
            max_iterations_per_node: 100_000,
            presolve: true,
        }
    }
}

/// The outcome of a search.
#[derive(Clone, Debug)]
pub struct Solution {
    pub status: Status,
    /// Objective in the user's original sense; `None` if no solution was found.
    pub objective: Option<f64>,
    /// The incumbent assignment, in column order.
    pub x: Vec<u8>,
    /// Best proven bound, in the user's original sense.
    pub bound: f64,
    pub nodes: usize,
    pub simplex_iterations: usize,
    /// What presolve achieved, if it ran.
    pub presolve: Option<presolve::Stats>,
}

impl Solution {
    /// Remaining optimality gap, relative to the incumbent. Zero when proven.
    pub fn gap(&self) -> f64 {
        match self.objective {
            Some(obj) if obj.is_finite() => {
                let denom = obj.abs().max(1e-10);
                ((obj - self.bound).abs() / denom).max(0.0)
            }
            _ => f64::INFINITY,
        }
    }
}

/// One open subproblem: the columns fixed so far, and where to resume from.
struct Node {
    /// Every column fixed on the path from the root, as `(column, value)`.
    ///
    /// Held in full rather than as a delta against a parent, which costs a little
    /// memory and removes a whole class of undo bugs from the search loop.
    fixings: Vec<(u32, u8)>,
    /// The parent's relaxation value, used to prune before solving anything.
    bound: f64,
    /// The parent's final basis, to warm start from.
    basis: BasisState,
}

/// Solve a binary integer program to proven optimality.
pub fn solve(problem: &Problem, options: Options) -> Solution {
    let started = Instant::now();

    // Presolve reduces in place and introduces no renumbering, so the reduced
    // model's solution vector is directly the original's -- there is no postsolve.
    // It is sound (it never admits a point the original rejects) and preserves the
    // optimum, so searching the reduced model answers the original question.
    let mut reduced;
    let (problem, presolve_stats) = if options.presolve {
        reduced = problem.clone();
        match presolve::presolve(&mut reduced, 20) {
            Outcome::Infeasible => {
                return Solution {
                    status: Status::Infeasible,
                    objective: None,
                    x: Vec::new(),
                    bound: f64::NAN,
                    nodes: 0,
                    simplex_iterations: 0,
                    presolve: None,
                };
            }
            Outcome::Reduced(stats) => (&reduced, Some(stats)),
        }
    } else {
        (problem, None)
    };

    let n = problem.n_cols();
    let mut lp = Lp::relaxation(problem);

    let mut nodes = 0usize;
    let mut iterations = 0usize;
    // Everything below is in the internal minimization form; conversion to the
    // user's sense happens once, on the way out.
    let mut incumbent = f64::INFINITY;
    let mut incumbent_x: Option<Vec<u8>> = None;

    let root = lp.solve_with_limit(options.max_iterations_per_node);
    iterations += root.iterations;
    nodes += 1;

    if root.status != LpStatus::Optimal {
        let status = if root.status == LpStatus::Infeasible {
            Status::Infeasible
        } else {
            Status::NodeLimit
        };
        return Solution {
            status,
            objective: None,
            x: Vec::new(),
            bound: f64::NAN,
            nodes,
            simplex_iterations: iterations,
            presolve: presolve_stats,
        };
    }

    let root_bound = root.objective;
    let mut stack: Vec<Node> = vec![Node {
        fixings: Vec::new(),
        bound: root_bound,
        basis: root.basis.clone(),
    }];

    // Consider the root itself first, so an integral relaxation is picked up without
    // a branching step.
    if let Some(x) = integral_solution(&root.x, options.integrality_tolerance) {
        incumbent = root.objective;
        incumbent_x = Some(x);
        stack.clear();
    }

    let mut status = Status::Optimal;

    while let Some(node) = stack.pop() {
        if nodes >= options.max_nodes {
            status = Status::NodeLimit;
            stack.push(node);
            break;
        }
        if options
            .time_limit
            .is_some_and(|limit| started.elapsed() >= limit)
        {
            status = Status::TimeLimit;
            stack.push(node);
            break;
        }
        // The parent's bound may have been overtaken by an incumbent found since this
        // node was pushed, in which case it needs no solve at all.
        if !improves(node.bound, incumbent, options.gap_tolerance) {
            continue;
        }

        // Rebuild this node's bounds from the root. Resetting first costs O(n) and
        // makes the node independent of whatever the previous one left behind.
        for j in 0..n {
            lp.set_column_bounds(j, problem.col_lb[j], problem.col_ub[j]);
        }
        for &(j, v) in &node.fixings {
            lp.set_column_bounds(j as usize, f64::from(v), f64::from(v));
        }

        let cutoff = if incumbent.is_finite() {
            Some(incumbent)
        } else {
            None
        };
        let solved = lp.solve_warm(&node.basis, cutoff, options.max_iterations_per_node);
        iterations += solved.iterations;
        nodes += 1;

        match solved.status {
            LpStatus::Infeasible | LpStatus::CutOff => continue,
            LpStatus::Unbounded => {
                // Impossible for a bounded binary relaxation; treat defensively.
                continue;
            }
            LpStatus::IterationLimit => {
                status = Status::NodeLimit;
                break;
            }
            LpStatus::Optimal => {}
        }

        if !improves(solved.objective, incumbent, options.gap_tolerance) {
            continue;
        }

        match integral_solution(&solved.x, options.integrality_tolerance) {
            Some(x) => {
                incumbent = solved.objective;
                incumbent_x = Some(x);
            }
            None => {
                let Some(branch) = most_fractional(&solved.x, options.integrality_tolerance) else {
                    continue;
                };
                let value = solved.x[branch];

                // Push both children; the stack pops last-in first, so the branch
                // nearer the relaxation's own value is pushed last and explored first.
                // Diving towards it tends to reach a feasible incumbent sooner, which
                // in turn makes the bound prune harder.
                let (first, second) = if value > 0.5 { (0u8, 1u8) } else { (1u8, 0u8) };
                for v in [first, second] {
                    let mut fixings = node.fixings.clone();
                    fixings.push((branch as u32, v));
                    stack.push(Node {
                        fixings,
                        bound: solved.objective,
                        basis: solved.basis.clone(),
                    });
                }
            }
        }
    }

    // The best bound is the weakest still-open node; with none left, the incumbent
    // is proven and the bound is the incumbent itself.
    let open_bound = stack
        .iter()
        .map(|node| node.bound)
        .fold(f64::INFINITY, f64::min);
    let internal_bound = if stack.is_empty() {
        if incumbent.is_finite() {
            incumbent
        } else {
            root_bound
        }
    } else {
        open_bound.min(incumbent)
    };

    let status = match (status, &incumbent_x) {
        (Status::Optimal, None) => Status::Infeasible,
        (other, _) => other,
    };

    Solution {
        status,
        objective: incumbent_x
            .as_ref()
            .map(|_| problem.objective_value(incumbent)),
        x: incumbent_x.unwrap_or_default(),
        bound: problem.objective_value(internal_bound),
        nodes,
        simplex_iterations: iterations,
        presolve: presolve_stats,
    }
}

/// Would a node with this bound improve on the incumbent by enough to be worth
/// exploring?
fn improves(bound: f64, incumbent: f64, gap_tolerance: f64) -> bool {
    if !incumbent.is_finite() {
        return true;
    }
    // An absolute floor as well as the relative gap: without it, an incumbent of
    // zero would make the relative test vacuous.
    let slack = (gap_tolerance * incumbent.abs()).max(1e-9);
    bound < incumbent - slack
}

/// Round a relaxation to a binary assignment, or `None` if it is fractional.
fn integral_solution(x: &[f64], tolerance: f64) -> Option<Vec<u8>> {
    let mut out = Vec::with_capacity(x.len());
    for &v in x {
        let rounded = v.round();
        if (v - rounded).abs() > tolerance {
            return None;
        }
        out.push(rounded as u8);
    }
    Some(out)
}

/// The column whose value sits closest to one half.
fn most_fractional(x: &[f64], tolerance: f64) -> Option<usize> {
    let mut best = None;
    let mut best_distance = 0.5;
    for (j, &v) in x.iter().enumerate() {
        let distance = (v - 0.5).abs();
        if (v - v.round()).abs() > tolerance && distance < best_distance {
            best_distance = distance;
            best = Some(j);
        }
    }
    best
}
