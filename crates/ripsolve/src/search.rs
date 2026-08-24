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

use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::time::{Duration, Instant};

use crate::branch::{self, Pseudocosts};
use crate::cuts;
use crate::heuristic::{self, Limits};
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
    /// Rounds of cut separation at the root. Zero disables cuts.
    pub cut_rounds: usize,
    /// Most cuts to add in a single round.
    pub cuts_per_round: usize,
    /// Consecutive depth-first steps before the search jumps to the best-bound
    /// open node. Zero makes the search pure best-bound.
    ///
    /// Defaults to zero, which measured best on four of five instances and turned
    /// `v064c1000n100` from a 77% gap into a 14-second solve. The textbook argument
    /// for plunging is that depth-first reaches incumbents sooner, but the primal
    /// heuristics already supply those, so the plunge buys little and costs bound
    /// progress. Raise it if memory becomes the binding constraint: a plunging
    /// search keeps its open set to roughly the tree depth, while best-bound holds
    /// every unexplored node.
    pub plunge_limit: usize,
    /// Run primal heuristics at the root, and every `heuristic_frequency` nodes.
    /// Zero disables them.
    pub heuristic_frequency: usize,
    /// Limits on the heuristics themselves.
    pub heuristic_limits: Limits,
    /// Total strong-branching probes allowed across the whole search.
    ///
    /// Defaults to zero — strong branching is implemented but **off**, because on
    /// every instance measured it made the search substantially worse, not better:
    /// v081c162n018 goes from 3828 nodes to 23686, and v128c256n100 from 10 to 917.
    /// Pure pseudocost scoring is the better default here. See the module docs on
    /// why the probes are suspected to be uninformative on these instances.
    pub strong_branching_budget: usize,
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
            // Separation converges after two or three rounds on every instance
            // measured; a larger budget finds no more cuts and only costs time.
            cut_rounds: 3,
            cuts_per_round: 32,
            plunge_limit: 0,
            heuristic_frequency: 100,
            heuristic_limits: Limits::default(),
            strong_branching_budget: 0,
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
    /// Cuts added at the root.
    pub cuts_added: usize,
    /// Incumbents found by a primal heuristic rather than by the search itself.
    pub heuristic_solutions: usize,
    /// The root relaxation before any cuts, in the user's original sense.
    pub root_bound: f64,
    /// The root relaxation after cutting, in the user's original sense.
    pub root_bound_after_cuts: f64,
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

/// A node's key in the best-bound pool: weakest bound first, deeper nodes winning
/// ties so the pool does not grow with shallow work.
///
/// `f64` is not `Ord`, and the bounds here are always finite, so `total_cmp` gives
/// the total order a heap needs without a newtype over the bit pattern.
#[derive(PartialEq)]
struct Key {
    bound: f64,
    depth: usize,
}

impl Eq for Key {}

impl Ord for Key {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.bound
            .total_cmp(&other.bound)
            .then(other.depth.cmp(&self.depth))
    }
}

impl PartialOrd for Key {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

/// The open node set: a depth-first plunge stack plus a best-bound pool.
///
/// Neither rule alone is good enough. Pure depth-first finds incumbents quickly —
/// a child differs from its parent by one bound, so its LP re-solves in a few
/// pivots — but it will happily spend an entire run inside one subtree while the
/// global bound never moves. Pure best-bound moves the bound as fast as possible
/// but wanders across the tree, warm-starting badly and finding incumbents late.
///
/// So: plunge for a bounded number of steps, then flush what is left of the dive
/// into the pool and take the weakest-bounded node in it. That is what makes the
/// reported gap close rather than just the incumbent improve.
struct OpenNodes {
    plunge: Vec<Node>,
    pool: BinaryHeap<(Reverse<Key>, usize)>,
    /// Nodes held out of the heap, indexed by the id stored alongside the key.
    store: Vec<Option<Node>>,
    steps: usize,
    limit: usize,
}

impl OpenNodes {
    fn new(limit: usize) -> Self {
        Self {
            plunge: Vec::new(),
            pool: BinaryHeap::new(),
            store: Vec::new(),
            steps: 0,
            limit,
        }
    }

    fn push(&mut self, node: Node) {
        self.plunge.push(node);
    }

    /// Move the whole dive into the pool, ending the plunge.
    fn flush(&mut self) {
        while let Some(node) = self.plunge.pop() {
            let key = Key {
                bound: node.bound,
                depth: node.fixings.len(),
            };
            let id = self.store.len();
            self.store.push(Some(node));
            self.pool.push((Reverse(key), id));
        }
        self.steps = 0;
    }

    fn pop(&mut self) -> Option<Node> {
        if self.steps >= self.limit {
            self.flush();
        }
        if let Some(node) = self.plunge.pop() {
            self.steps += 1;
            return Some(node);
        }
        self.steps = 0;
        while let Some((_, id)) = self.pool.pop() {
            if let Some(node) = self.store[id].take() {
                return Some(node);
            }
        }
        None
    }

    fn is_empty(&self) -> bool {
        self.plunge.is_empty() && self.pool.iter().all(|&(_, id)| self.store[id].is_none())
    }

    /// The weakest bound among all open nodes, which is the search's global bound.
    fn best_bound(&self) -> f64 {
        let from_plunge = self
            .plunge
            .iter()
            .map(|n| n.bound)
            .fold(f64::INFINITY, f64::min);
        let from_pool = self
            .pool
            .iter()
            .filter(|&&(_, id)| self.store[id].is_some())
            .map(|(Reverse(key), _)| key.bound)
            .fold(f64::INFINITY, f64::min);
        from_plunge.min(from_pool)
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
    /// The branch that created this node: `(column, went_up, parent objective,
    /// parent fractional value)`. Recorded so that solving this node measures the
    /// real cost of that branching decision and feeds it back.
    origin: Option<(usize, bool, f64, f64)>,
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
                    cuts_added: 0,
                    heuristic_solutions: 0,
                    root_bound: f64::NAN,
                    root_bound_after_cuts: f64::NAN,
                };
            }
            Outcome::Reduced(stats) => (&reduced, Some(stats)),
        }
    } else {
        (problem, None)
    };

    let n = problem.n_cols();
    let mut lp = Lp::relaxation(problem);
    // Cuts add rows but never columns, so `n` stays valid across the cut loop.

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
            cuts_added: 0,
            heuristic_solutions: 0,
            root_bound: f64::NAN,
            root_bound_after_cuts: f64::NAN,
        };
    }

    let first_bound = root.objective;

    // Cut at the root only. Separating deeper in the tree ("local cuts") is
    // stronger still, but it needs the cut pool to track which nodes each cut is
    // valid for; these are globally valid, so adding them once to the model is
    // both simpler and correct everywhere.
    let mut cuts_added = 0usize;
    let mut root = root;
    let mut with_cuts;
    let mut problem = problem;
    for _ in 0..options.cut_rounds {
        if root.status != LpStatus::Optimal
            || integral_solution(&root.x, options.integrality_tolerance).is_some()
        {
            break;
        }
        // Two families with different reach: covers need a row that reads as a
        // knapsack, while GMI comes off the tableau and applies to any fractional
        // basic column. On dense random rows the second is usually the only one that
        // finds anything.
        let mut found = cuts::separate(problem, &root.x, options.cuts_per_round);
        found.extend(cuts::separate_gomory(
            &lp,
            &root.basis,
            &root.x,
            options.cuts_per_round,
        ));
        if found.is_empty() {
            break;
        }
        with_cuts = problem.clone();
        with_cuts.add_cuts(&found);
        cuts_added += found.len();

        let candidate = Lp::relaxation(&with_cuts);
        let resolved = candidate.solve_with_limit(options.max_iterations_per_node);
        iterations += resolved.iterations;
        if resolved.status != LpStatus::Optimal {
            // Keep the model that is known to solve rather than pressing on with one
            // that does not; the bound already gained is still sound.
            break;
        }
        reduced = with_cuts;
        problem = &reduced;
        lp = candidate;
        root = resolved;
    }

    let root_bound = root.objective;
    let mut open = OpenNodes::new(options.plunge_limit);
    open.push(Node {
        fixings: Vec::new(),
        bound: root_bound,
        basis: root.basis.clone(),
        origin: None,
    });
    let mut pseudocosts = Pseudocosts::new(n);
    let mut strong_budget = options.strong_branching_budget;
    let mut heuristic_solutions = 0usize;

    // An incumbent before the first branch is worth more than one found later: the
    // search cannot prune anything until it holds one.
    if options.heuristic_frequency > 0 && root.status == LpStatus::Optimal {
        // Cheapest first. Rounding costs no LP at all; diving costs a short chain of
        // them; the pump costs the most but is the only one that reliably finds
        // anything on models whose feasible set is sparse.
        let found = heuristic::round(problem, &root.x, &options.heuristic_limits)
            .or_else(|| {
                heuristic::dive(
                    problem,
                    &mut lp,
                    &root.basis,
                    &root.x,
                    None,
                    &options.heuristic_limits,
                    &mut iterations,
                )
            })
            .or_else(|| {
                heuristic::feasibility_pump(
                    problem,
                    &mut lp,
                    &root.basis,
                    &root.x,
                    &options.heuristic_limits,
                    &mut iterations,
                )
            });
        if let Some(found) = found
            && found.objective < incumbent
        {
            incumbent = found.objective;
            incumbent_x = Some(found.x);
            heuristic_solutions += 1;
        }
    }

    // Consider the root itself first, so an integral relaxation is picked up without
    // a branching step.
    if let Some(x) = integral_solution(&root.x, options.integrality_tolerance) {
        incumbent = root.objective;
        incumbent_x = Some(x);
        open = OpenNodes::new(options.plunge_limit);
    }

    let mut status = Status::Optimal;

    while let Some(node) = open.pop() {
        if nodes >= options.max_nodes {
            status = Status::NodeLimit;
            open.push(node);
            break;
        }
        if options
            .time_limit
            .is_some_and(|limit| started.elapsed() >= limit)
        {
            status = Status::TimeLimit;
            open.push(node);
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

        // Whatever this node's LP turned out to be, it is the measured consequence
        // of the branch that created it; feed that back before doing anything else.
        if let Some((column, up, parent_objective, fraction)) = node.origin
            && solved.status == LpStatus::Optimal
        {
            let step = if up { 1.0 - fraction } else { fraction };
            let degradation = (solved.objective - parent_objective).max(0.0);
            pseudocosts.record(column, up, degradation, step);
        }

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

        // Dive from this node every so often. Nodes deep in the tree already have
        // many columns fixed, so a dive from one is short and often lands somewhere
        // the search would take a long time to reach.
        if options.heuristic_frequency > 0
            && nodes.is_multiple_of(options.heuristic_frequency)
            && integral_solution(&solved.x, options.integrality_tolerance).is_none()
        {
            let cutoff = incumbent.is_finite().then_some(incumbent);
            let found = heuristic::dive(
                problem,
                &mut lp,
                &solved.basis,
                &solved.x,
                cutoff,
                &options.heuristic_limits,
                &mut iterations,
            );
            if let Some(found) = found
                && found.objective < incumbent
            {
                incumbent = found.objective;
                incumbent_x = Some(found.x);
                heuristic_solutions += 1;
            }
        }

        match integral_solution(&solved.x, options.integrality_tolerance) {
            Some(x) => {
                incumbent = solved.objective;
                incumbent_x = Some(x);
            }
            None => {
                let decision = branch::select(
                    &mut lp,
                    &solved.basis,
                    &solved.x,
                    solved.objective,
                    &mut pseudocosts,
                    options.integrality_tolerance,
                    &mut strong_budget,
                    &mut iterations,
                );
                let Some(decision) = decision else { continue };
                // Strong branching can prove a node has no feasible completion at
                // all, which prunes it outright.
                if decision.dead {
                    continue;
                }

                // A probe that proved one side infeasible decides the column outright,
                // so descend into the single surviving child instead of branching.
                let children: Vec<u8> = match decision.forced {
                    Some(value) => vec![value],
                    None => {
                        // The stack pops last-in first, so the side nearer the
                        // relaxation's own value is pushed last and explored first.
                        // Diving towards it reaches a feasible incumbent sooner,
                        // which makes the bound prune harder.
                        if decision.fraction > 0.5 {
                            vec![0, 1]
                        } else {
                            vec![1, 0]
                        }
                    }
                };

                for v in children {
                    let mut fixings = node.fixings.clone();
                    fixings.push((decision.column as u32, v));
                    open.push(Node {
                        fixings,
                        bound: solved.objective,
                        basis: solved.basis.clone(),
                        origin: Some((
                            decision.column,
                            v == 1,
                            solved.objective,
                            decision.fraction,
                        )),
                    });
                }
            }
        }
    }

    // The best bound is the weakest still-open node; with none left, the incumbent
    // is proven and the bound is the incumbent itself.
    let internal_bound = if open.is_empty() {
        if incumbent.is_finite() {
            incumbent
        } else {
            root_bound
        }
    } else {
        open.best_bound().min(incumbent)
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
        cuts_added,
        heuristic_solutions,
        root_bound: problem.objective_value(first_bound),
        root_bound_after_cuts: problem.objective_value(root_bound),
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lp::Lp;
    use crate::model::{RowSense, Sense};
    use crate::sparse::SparseMatrix;

    fn node(bound: f64, depth: usize) -> Node {
        // A node's basis is irrelevant to the ordering under test.
        let lp = Lp::relaxation(&trivial());
        Node {
            fixings: (0..depth).map(|j| (j as u32, 0u8)).collect(),
            bound,
            basis: lp.solve().basis,
            origin: None,
        }
    }

    fn trivial() -> Problem {
        Problem {
            name: "t".into(),
            sense: Sense::Minimize,
            obj: vec![1.0],
            obj_offset: 0.0,
            matrix: SparseMatrix::from_triplets(1, 1, [(0, 0, 1.0)]),
            row_lb: vec![RowSense::Ge.bounds(0.0).0],
            row_ub: vec![RowSense::Ge.bounds(0.0).1],
            col_lb: vec![0.0],
            col_ub: vec![1.0],
            col_names: vec!["x".into()],
            row_names: vec!["c".into()],
        }
    }

    #[test]
    fn a_zero_plunge_limit_is_pure_best_bound() {
        let mut open = OpenNodes::new(0);
        for bound in [5.0, 1.0, 3.0] {
            open.push(node(bound, 0));
        }
        // Weakest bound first, regardless of push order.
        let order: Vec<f64> = std::iter::from_fn(|| open.pop().map(|n| n.bound)).collect();
        assert_eq!(order, vec![1.0, 3.0, 5.0]);
    }

    #[test]
    fn plunging_takes_the_most_recent_node_first() {
        let mut open = OpenNodes::new(10);
        open.push(node(5.0, 0));
        open.push(node(1.0, 0));
        // Depth-first ignores the bound while the plunge lasts.
        assert_eq!(open.pop().map(|n| n.bound), Some(1.0));
        assert_eq!(open.pop().map(|n| n.bound), Some(5.0));
    }

    #[test]
    fn the_plunge_ends_after_its_limit() {
        let mut open = OpenNodes::new(2);
        for bound in [9.0, 8.0, 1.0] {
            open.push(node(bound, 0));
        }
        // Two depth-first steps, then a jump to the weakest remaining bound.
        assert_eq!(open.pop().map(|n| n.bound), Some(1.0));
        assert_eq!(open.pop().map(|n| n.bound), Some(8.0));
        assert_eq!(open.pop().map(|n| n.bound), Some(9.0));
    }

    #[test]
    fn ties_prefer_the_deeper_node() {
        // Otherwise the pool fills with shallow work that never gets finished.
        let mut open = OpenNodes::new(0);
        open.push(node(2.0, 1));
        open.push(node(2.0, 7));
        assert_eq!(open.pop().map(|n| n.fixings.len()), Some(7));
    }

    #[test]
    fn the_global_bound_is_the_weakest_open_node() {
        let mut open = OpenNodes::new(0);
        open.push(node(4.0, 0));
        open.push(node(2.5, 0));
        assert_eq!(open.best_bound(), 2.5);
        // Once that node is taken, the bound rises to what is left.
        open.pop();
        assert_eq!(open.best_bound(), 4.0);
        open.pop();
        assert!(open.is_empty());
        assert_eq!(open.best_bound(), f64::INFINITY);
    }

    #[test]
    fn emptiness_accounts_for_nodes_already_taken() {
        // Popping leaves a spent entry in the heap; it must not read as still open.
        let mut open = OpenNodes::new(0);
        open.push(node(1.0, 0));
        open.pop();
        assert!(open.is_empty());
        assert!(open.pop().is_none());
    }
}
