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
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::sync::{Condvar, Mutex};
use std::time::{Duration, Instant};

use crate::branch::{self, Pseudocosts};
use crate::cuts;
use crate::heuristic::{self, Limits, Schedule};
use crate::lp::{BasisState, Lp, LpSolution, LpStatus};
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
    ///
    /// Defaults to `1e-4`, which is what the established solvers use, and the reason
    /// to match them is not convention but arithmetic: proving the last ten-thousandth
    /// of a percent can cost more than everything before it. On MIPLIB's `app2-1` the
    /// search reaches a 0.0006% gap quickly and then cannot close it, so at zero it
    /// runs out the clock with the answer already in hand; at `1e-4` it finishes in
    /// 1.1s, which is faster than HiGHS on the same model.
    ///
    /// Set it to zero to demand a proof. `Status::Optimal` then means exactly that;
    /// otherwise it means optimal to within this gap, as it does elsewhere.
    pub gap_tolerance: f64,
    /// Simplex iteration limit for a single node.
    pub max_iterations_per_node: usize,
    /// Tighten the model before searching.
    pub presolve: bool,
    /// Separate cuts at one node in every `local_cut_frequency`, not only at the
    /// root. Zero disables node-local cutting.
    ///
    /// Root cuts turn out to be a shallow-depth phenomenon: measured over six models
    /// they bind at 33-50% of rows at depth one and 0-4% by depth ten, while the tree
    /// carries them through every node. Cuts derived at a node use that node's bounds,
    /// so they bind where they were made, but they are valid only for that subtree,
    /// which is why they never enter the shared model.
    ///
    /// Ten is from a sweep of 0, 1, 3, 10, 50 and 200 over eight models. It is the only
    /// setting that beats not cutting at all, 22.98s against 24.75s, while still taking
    /// real chunks out of the tree, `v064c1000n100` 1106 nodes to 786, `mkp_200` 72150
    /// to 63896. Separating at *every* node shrinks trees far harder (`v256c256n100`
    /// 288 nodes to 86, `v064c200` 2690 to 1136) but costs 36.50s: worth reaching for
    /// on a model where nodes are the bottleneck, not as a default.
    pub local_cut_frequency: usize,
    /// Most cuts to derive at one node.
    pub local_cuts_per_node: usize,
    /// Rounds of cut separation at the root. Zero disables cuts.
    ///
    /// Defaults to zero, which is not where a branch-and-*cut* solver expects to
    /// end up. Measured over eleven models spanning this solver's target range,
    /// cutting was slower on every one of them:
    ///
    /// ```text
    ///                  no cuts      3 rounds of 32
    ///     mkp_200      17.1s          47.8s
    ///     v064c200      1.7s           2.4s
    ///     v081c162n009  0.8s           1.2s
    ///     v048c048     0.01s          0.04s
    /// ```
    ///
    /// The cuts are not useless, they raise the root bound substantially, taking
    /// v064c200 from 72.1 to 95.4 against an optimum of 225. They simply do not pay
    /// for themselves here: separation is expensive, the cuts come out dense enough
    /// to slow every subsequent LP, and best-bound node selection had already
    /// captured most of what a better bound was worth. On `mkp_200` cutting even
    /// *raised* the node count, 72150 to 91346, and on MIPLIB's markshare_4_0 it was
    /// the difference between proving optimality in 21 seconds and not proving it at
    /// all.
    ///
    /// Turning this on is worthwhile when node count matters more than wall clock,
    /// and would become worthwhile generally with cheaper separation and proper cut
    /// selection by efficacy and orthogonality. Neither is implemented.
    pub cut_rounds: usize,
    /// Most cuts to add in a single round.
    pub cuts_per_round: usize,
    /// Product-form updates before refactorizing the basis; see
    /// [`crate::lp::Tolerances::refactor_interval`].
    pub refactor_interval: usize,
    /// Worker threads for the tree search. One (or zero) runs the serial driver.
    ///
    /// Only the tree is parallel: presolve, cut generation and the root heuristics
    /// all run once, before any thread is spawned.
    pub threads: usize,
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
    /// Node interval between improvement searches. Zero disables them.
    ///
    /// Once an incumbent exists, the cheapest place to look for a better one is the
    /// neighbourhood where the incumbent and the relaxation already agree. Fixing the
    /// integer columns they agree on leaves a much smaller model whose optimum is a
    /// better solution to the original, and it is reached by the same search rather
    /// than by a separate mechanism.
    ///
    /// This is what the other heuristics here do not do: rounding, diving and the pump
    /// all *find* a solution and none of them improves one. On MIPLIB's
    /// graphdraw-gemcutter the search reaches 13176 against an optimum of 7118 and then
    /// sits there, not for want of nodes but for want of anything that looks near a
    /// good solution rather than near the relaxation.
    pub improvement_frequency: usize,
    /// Nodes an improvement search may spend before giving the budget back.
    pub improvement_nodes: usize,
    /// Base node interval between in-tree heuristic attempts. Zero disables them.
    ///
    /// This is a starting point, not a fixed cadence: the interval doubles after
    /// each attempt that finds nothing and snaps back here after one that does.
    pub heuristic_frequency: usize,
    /// Limits on the heuristics themselves.
    pub heuristic_limits: Limits,
    /// Total strong-branching probes allowed across the whole search.
    ///
    /// Defaults to zero. Under best-bound selection strong branching does reduce
    /// node counts (by 10% to 32% on seven of eight instances) but two extra LPs
    /// per candidate cost more time than those nodes save. Worth raising on large
    /// models, where it pays: v128c1000n100 goes from 13.3s to 9.5s at a budget of
    /// 100. See `branch.rs` for why it was catastrophic under depth-first search
    /// and is merely unprofitable now.
    pub strong_branching_budget: usize,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            max_nodes: usize::MAX,
            time_limit: None,
            integrality_tolerance: 1e-6,
            gap_tolerance: 1e-4,
            max_iterations_per_node: 100_000,
            presolve: true,
            local_cut_frequency: 10,
            local_cuts_per_node: 8,
            // Off by default. See `cut_rounds`.
            cut_rounds: 0,
            cuts_per_round: 32,
            refactor_interval: 200,
            threads: 1,
            plunge_limit: 0,
            improvement_frequency: 500,
            improvement_nodes: 2_000,
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
    pub x: Vec<f64>,
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
/// Neither rule alone is good enough. Pure depth-first finds incumbents quickly,
/// a child differs from its parent by one bound, so its LP re-solves in a few
/// pivots, but it will happily spend an entire run inside one subtree while the
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

/// What processing one node produced.
struct NodeOutcome {
    /// A feasible assignment found here, by the LP landing integral or by a
    /// heuristic.
    incumbent: Option<(f64, Vec<f64>)>,
    children: Vec<Node>,
    /// The node's LP hit its iteration limit, so the search cannot claim to have
    /// examined this subtree and must not report optimality.
    exhausted: bool,
    /// Incumbents that came from a heuristic rather than from the relaxation.
    heuristic_hits: usize,
    /// This node's relaxation, kept for the improvement search to compare against.
    relaxation: Option<Vec<f64>>,
}

/// Everything needed to process nodes: an LP to re-solve them in, and the
/// branching history built up while doing so.
///
/// Node processing lives here rather than inline in the search loop so that the
/// serial and parallel drivers run *the same* code. The logic is subtle, cutoff
/// handling, pseudocost feedback, forced fixings, and two copies of it would
/// diverge.
struct Worker<'a> {
    problem: &'a Problem,
    lp: Lp,
    pseudocosts: Pseudocosts,
    strong_budget: usize,
    options: Options,
    iterations: usize,
    dives: Schedule,
}

impl<'a> Worker<'a> {
    fn new(problem: &'a Problem, lp: Lp, options: Options) -> Self {
        let n = problem.n_cols();
        Self {
            problem,
            lp,
            pseudocosts: Pseudocosts::new(n),
            strong_budget: options.strong_branching_budget,
            options,
            iterations: 0,
            dives: Schedule::new(options.heuristic_frequency),
        }
    }

    /// Solve one node and decide what becomes of it.
    ///
    /// `incumbent` is the best objective known *now*; in a parallel search that may
    /// be better than when the node was created, which only ever prunes more.
    /// Tighten one node's bound with cuts derived from its own relaxation.
    ///
    /// Returns the improved bound, or `None` if nothing was gained. Cuts separated
    /// here are read off a tableau built under this node's branching bounds, so they
    /// are valid for this subtree and not for the tree at large. They are therefore
    /// never added to the shared model: the augmented LP lives for exactly one solve,
    /// and the only thing that outlives it is the bound, which *is* valid everywhere
    /// below this node and so is safe to prune and order children with.
    ///
    /// The node's own basis and solution are left untouched, so branching still reads
    /// the same relaxation it would have without this. That is deliberate, the
    /// augmented vertex is arguably the better branching point, but mixing it with a
    /// basis taken from the un-augmented model is a correctness trap not worth
    /// setting for a first cut at this.
    fn separate_locally(
        &mut self,
        solved: &LpSolution,
        options: &Options,
        cutoff: Option<f64>,
    ) -> Option<f64> {
        let limit = options.local_cuts_per_node;
        let found = cuts::separate_gomory(&self.lp, &solved.basis, &solved.x, limit);
        let found = cuts::select(found, &solved.x, limit);
        if found.is_empty() {
            return None;
        }

        let rows: Vec<crate::lp::RangeRow> = found
            .iter()
            .map(|c| (c.coefficients.clone(), c.lb, c.ub))
            .collect();
        let resolved = self.lp.solve_with_rows(
            &solved.basis,
            &rows,
            cutoff,
            options.max_iterations_per_node,
        );
        self.iterations += resolved.iterations;

        match resolved.status {
            // The tightened relaxation is already worse than the incumbent, so
            // everything below this node is too.
            LpStatus::CutOff => Some(f64::INFINITY),
            LpStatus::Optimal if resolved.objective > solved.objective => Some(resolved.objective),
            _ => None,
        }
    }

    fn process(&mut self, node: &Node, incumbent: f64, index: usize) -> NodeOutcome {
        let problem = self.problem;
        let options = self.options;
        let mut out = NodeOutcome {
            incumbent: None,
            relaxation: None,
            children: Vec::new(),
            exhausted: false,
            heuristic_hits: 0,
        };

        // Rebuild this node's bounds from the root. Resetting first costs O(n) and
        // makes the node independent of whatever the previous one left behind.
        for j in 0..problem.n_cols() {
            self.lp
                .set_column_bounds(j, problem.col_lb[j], problem.col_ub[j]);
        }
        for &(j, lo, hi) in &node.fixings {
            self.lp.set_column_bounds(j as usize, lo, hi);
        }

        let cutoff = incumbent.is_finite().then_some(incumbent);
        let solved = self
            .lp
            .solve_warm(&node.basis, cutoff, options.max_iterations_per_node);
        self.iterations += solved.iterations;

        // Whatever this node's LP turned out to be, it is the measured consequence
        // of the branch that created it; feed that back before doing anything else.
        if let Some((column, up, parent_objective, fraction)) = node.origin
            && solved.status == LpStatus::Optimal
        {
            let step = if up { 1.0 - fraction } else { fraction };
            let degradation = (solved.objective - parent_objective).max(0.0);
            self.pseudocosts.record(column, up, degradation, step);
        }

        match solved.status {
            // Unbounded is impossible for a bounded binary relaxation; treated as a
            // dead node defensively rather than trusted.
            LpStatus::Infeasible | LpStatus::CutOff | LpStatus::Unbounded => return out,
            LpStatus::IterationLimit => {
                out.exhausted = true;
                return out;
            }
            LpStatus::Optimal => {}
        }

        if !improves(solved.objective, incumbent, options.gap_tolerance) {
            return out;
        }
        // Kept only when an improvement search might use it, since it is a copy of the
        // whole column vector at every node otherwise.
        if options.improvement_frequency > 0 && index.is_multiple_of(options.improvement_frequency)
        {
            out.relaxation = Some(solved.x.clone());
        }

        // Cuts derived from this node's own relaxation, if it is one of the nodes
        // chosen for it. Only the bound escapes; see `separate_locally`.
        let mut bound = solved.objective;
        if options.local_cut_frequency > 0
            && index.is_multiple_of(options.local_cut_frequency)
            && integral_solution(problem, &solved.x, options.integrality_tolerance).is_none()
            && let Some(tightened) = self.separate_locally(&solved, &options, cutoff)
        {
            bound = tightened;
            if !improves(bound, incumbent, options.gap_tolerance) {
                return out;
            }
        }

        // Dive from this node every so often. Nodes deep in the tree already have
        // many columns fixed, so a dive from one is short and often lands somewhere
        // the search would take a long time to reach.
        if self.dives.due(index)
            && integral_solution(problem, &solved.x, options.integrality_tolerance).is_none()
        {
            let found = heuristic::dive(
                problem,
                &mut self.lp,
                &solved.basis,
                &solved.x,
                cutoff,
                &options.heuristic_limits,
                &mut self.iterations,
            );
            let found = found.map(|found| {
                let polished = heuristic::polish(
                    problem,
                    &mut self.lp,
                    &solved.basis,
                    &found.x,
                    &options.heuristic_limits,
                    &mut self.iterations,
                );
                match polished {
                    Some(better) if better.objective < found.objective => better,
                    _ => found,
                }
            });
            let improved = match found {
                Some(found) if found.objective < incumbent => {
                    out.incumbent = Some((found.objective, found.x));
                    out.heuristic_hits += 1;
                    true
                }
                _ => false,
            };
            self.dives.record(index, improved);
        }

        match integral_solution(problem, &solved.x, options.integrality_tolerance) {
            Some(x) => {
                let better = out
                    .incumbent
                    .as_ref()
                    .is_none_or(|(o, _)| solved.objective < *o);
                if better {
                    out.incumbent = Some((solved.objective, x));
                }
            }
            None => {
                let decision = branch::select(
                    problem,
                    &mut self.lp,
                    &solved.basis,
                    &solved.x,
                    solved.objective,
                    &mut self.pseudocosts,
                    options.integrality_tolerance,
                    &mut self.strong_budget,
                    &mut self.iterations,
                );
                let Some(decision) = decision else { return out };
                // Strong branching can prove a node has no feasible completion at
                // all, which prunes it outright.
                if decision.dead {
                    return out;
                }

                let column = decision.column;
                let value = solved.x[column];
                let (lo, hi) = self.lp.column_bounds(column);
                // Split the column's range at the fractional value: `x <= floor(v)`
                // against `x >= ceil(v)`. Together these cover every integer the
                // column could take and overlap in none, so no solution is lost and
                // none is counted twice. On a binary column this is exactly fixing
                // to 0 and to 1.
                let down = (lo, value.floor());
                let up = (value.ceil(), hi);

                let mut sides = match decision.forced {
                    Some(0) => vec![down],
                    Some(_) => vec![up],
                    // Explored last-pushed first under a plunge, so the side nearer
                    // the relaxation's own value goes last. Immaterial under
                    // best-bound selection.
                    None if value - value.floor() > 0.5 => vec![down, up],
                    None => vec![up, down],
                };
                sides.retain(|&(lo, hi)| lo <= hi);

                for (lo, hi) in sides {
                    let mut fixings = node.fixings.clone();
                    fixings.push((column as u32, lo, hi));
                    // `up` is the side that raised the lower bound.
                    let went_up = lo > value.floor();
                    out.children.push(Node {
                        fixings,
                        bound,
                        basis: solved.basis.clone(),
                        origin: Some((column, went_up, solved.objective, decision.fraction)),
                    });
                }
            }
        }
        out
    }
}

/// What a tree search produced, whether it ran on one thread or several.
struct TreeResult {
    status: Status,
    incumbent: f64,
    incumbent_x: Option<Vec<f64>>,
    nodes: usize,
    iterations: usize,
    heuristic_hits: usize,
    /// The weakest bound left open, or infinity if the tree was exhausted.
    open_bound: f64,
}

/// The open node set and the count of workers currently holding a node.
///
/// `active` is what makes termination decidable: an empty pool does not mean the
/// search is finished, only that every remaining node is currently being expanded
/// by some worker and may yet produce children. The search is over when the pool is
/// empty *and* no worker is active.
struct SharedPool {
    open: OpenNodes,
    active: usize,
    finished: bool,
}

/// State shared by every worker in a parallel search.
struct Shared {
    pool: Mutex<SharedPool>,
    /// Signalled whenever a worker pushes children or goes idle.
    wake: Condvar,
    /// The incumbent, guarded for writing.
    best: Mutex<(f64, Option<Vec<f64>>)>,
    /// The incumbent objective again, as raw bits, for lock-free reads on the hot
    /// path. Always written under `best`, so it can only lag, never lead, and a
    /// stale-but-worse cutoff prunes less, never wrongly.
    best_bits: AtomicU64,
    nodes: AtomicUsize,
    iterations: AtomicUsize,
    heuristic_hits: AtomicUsize,
    /// Set when a limit stops the search; read by every worker to wind down.
    stopped: AtomicUsize,
    /// Nodes whose LP never resolved. Any of them could hold the optimum.
    unresolved: AtomicUsize,
}

/// Reasons a parallel search stops early, encoded for the atomic.
/// Consecutive slack resolves after which a cut is dropped from the model.
///
/// Swept over 0 (never drop) through 3 with the root budget at three rounds. Two is
/// the best of them, though not by much: it takes `v128c1000n100` from 740 nodes to
/// 610 and `v064c200` from 2932 to 2716, and is neutral on everything else measured.
/// Dropping at the first slack resolve is too eager, a cut can go slack for one
/// round and bind again once the next round's cuts move the vertex.
const CUT_MAX_AGE: u32 = 2;

/// Row activity within this of a cut's bound counts as binding. This is a tolerance
/// on constraint activity, which is why it is not `integrality_tolerance`.
const CUT_SLACK_TOLERANCE: f64 = 1e-7;

const STOP_NONE: usize = 0;
const STOP_NODES: usize = 1;
const STOP_TIME: usize = 2;

impl Shared {
    fn incumbent(&self) -> f64 {
        f64::from_bits(self.best_bits.load(Ordering::Relaxed))
    }

    /// Install a better incumbent, if it really is better.
    ///
    /// Two workers can find improvements concurrently, so the comparison is redone
    /// under the lock rather than trusted from the atomic read that motivated it.
    fn offer(&self, objective: f64, x: Vec<f64>) -> bool {
        let mut best = self.best.lock().expect("incumbent lock");
        if objective < best.0 {
            best.0 = objective;
            best.1 = Some(x);
            self.best_bits.store(objective.to_bits(), Ordering::Relaxed);
            return true;
        }
        false
    }

    /// Take the next node, or `None` once the search is over.
    fn take(&self) -> Option<Node> {
        let mut pool = self.pool.lock().expect("pool lock");
        loop {
            if pool.finished {
                return None;
            }
            if self.stopped.load(Ordering::Relaxed) != STOP_NONE {
                pool.finished = true;
                self.wake.notify_all();
                return None;
            }
            if let Some(node) = pool.open.pop() {
                pool.active += 1;
                return Some(node);
            }
            if pool.active == 0 {
                // Nothing open and nobody working: the tree is exhausted.
                pool.finished = true;
                self.wake.notify_all();
                return None;
            }
            pool = self.wake.wait(pool).expect("pool wait");
        }
    }

    /// Return a node's children and release the worker.
    fn give_back(&self, children: Vec<Node>) {
        let mut pool = self.pool.lock().expect("pool lock");
        pool.active -= 1;
        for child in children {
            pool.open.push(child);
        }
        self.wake.notify_all();
    }
}

/// Run the tree search across `threads` workers.
///
/// Node counts and iteration counts vary between runs, because which node a worker
/// takes depends on timing. The *answer* does not: every worker prunes against a
/// shared incumbent and every cut and bound is globally valid, so the proven
/// optimum is the same however the work is divided.
#[allow(clippy::too_many_arguments)]
fn run_parallel(
    problem: &Problem,
    lp: &Lp,
    open: OpenNodes,
    options: Options,
    started: Instant,
    incumbent: f64,
    incumbent_x: Option<Vec<f64>>,
    threads: usize,
) -> TreeResult {
    let shared = Shared {
        pool: Mutex::new(SharedPool {
            open,
            active: 0,
            finished: false,
        }),
        wake: Condvar::new(),
        best: Mutex::new((incumbent, incumbent_x)),
        best_bits: AtomicU64::new(incumbent.to_bits()),
        nodes: AtomicUsize::new(0),
        iterations: AtomicUsize::new(0),
        heuristic_hits: AtomicUsize::new(0),
        stopped: AtomicUsize::new(STOP_NONE),
        unresolved: AtomicUsize::new(0),
    };

    std::thread::scope(|scope| {
        for _ in 0..threads {
            // Each worker gets its own LP, since solving a node mutates the column
            // bounds, and its own branching history. Sharing the pseudocosts would
            // pool more evidence but put a lock on the hot path; per-worker history
            // is the cheaper trade at these thread counts.
            let mut worker = Worker::new(problem, lp.clone(), options);
            let shared = &shared;
            scope.spawn(move || {
                while let Some(node) = shared.take() {
                    let index = shared.nodes.fetch_add(1, Ordering::Relaxed) + 1;

                    if index > options.max_nodes {
                        shared.stopped.store(STOP_NODES, Ordering::Relaxed);
                        shared.give_back(vec![node]);
                        continue;
                    }
                    if options
                        .time_limit
                        .is_some_and(|limit| started.elapsed() >= limit)
                    {
                        shared.stopped.store(STOP_TIME, Ordering::Relaxed);
                        shared.give_back(vec![node]);
                        continue;
                    }

                    let best = shared.incumbent();
                    if !improves(node.bound, best, options.gap_tolerance) {
                        shared.give_back(Vec::new());
                        continue;
                    }

                    let outcome = worker.process(&node, best, index);
                    if outcome.exhausted {
                        // Skip the node and keep going; see the serial driver. The
                        // count is what stops the search claiming optimality.
                        shared.unresolved.fetch_add(1, Ordering::Relaxed);
                    }
                    if let Some((objective, x)) = outcome.incumbent
                        && shared.offer(objective, x)
                        && outcome.heuristic_hits > 0
                    {
                        shared.heuristic_hits.fetch_add(1, Ordering::Relaxed);
                    }
                    shared.give_back(outcome.children);
                }
                shared
                    .iterations
                    .fetch_add(worker.iterations, Ordering::Relaxed);
            });
        }
    });

    let pool = shared.pool.into_inner().expect("pool lock");
    let (incumbent, incumbent_x) = shared.best.into_inner().expect("incumbent lock");
    let status = match shared.stopped.load(Ordering::Relaxed) {
        STOP_NODES => Status::NodeLimit,
        STOP_TIME => Status::TimeLimit,
        // A tree that was exhausted apart from a skipped node proves nothing.
        _ if shared.unresolved.load(Ordering::Relaxed) > 0 => Status::NodeLimit,
        _ => Status::Optimal,
    };

    TreeResult {
        status,
        incumbent,
        incumbent_x,
        nodes: shared.nodes.load(Ordering::Relaxed).min(options.max_nodes),
        iterations: shared.iterations.load(Ordering::Relaxed),
        heuristic_hits: shared.heuristic_hits.load(Ordering::Relaxed),
        open_bound: if pool.open.is_empty() {
            f64::INFINITY
        } else {
            pool.open.best_bound()
        },
    }
}

/// One open subproblem: the columns fixed so far, and where to resume from.
struct Node {
    /// Every bound change on the path from the root, as `(column, lower, upper)`.
    ///
    /// Bound changes rather than fixings, because branching on a general integer
    /// splits a range (`x <= floor(v)` against `x >= ceil(v)`) instead of pinning
    /// a value. On a binary column that split *is* a fixing, so the binary case
    /// needs no separate handling.
    ///
    /// Held in full rather than as a delta against a parent, which costs a little
    /// memory and removes a whole class of undo bugs from the search loop.
    fixings: Vec<(u32, f64, f64)>,
    /// The parent's relaxation value, used to prune before solving anything.
    bound: f64,
    /// The parent's final basis, to warm start from.
    basis: BasisState,
    /// The branch that created this node: `(column, went_up, parent objective,
    /// parent fractional value)`. Recorded so that solving this node measures the
    /// real cost of that branching decision and feeds it back.
    origin: Option<(usize, bool, f64, f64)>,
}

/// Solve a mixed-integer program to proven optimality.
pub fn solve(problem: &Problem, options: Options) -> Solution {
    let started = Instant::now();

    // Presolve reduces in place and introduces no renumbering, so the reduced
    // model's solution vector is directly the original's, there is no postsolve.
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

    let _ = problem.n_cols();
    let mut lp = Lp::relaxation(problem).with_tolerances(crate::lp::Tolerances {
        refactor_interval: options.refactor_interval,
        ..Default::default()
    });
    // The budget belongs to the LP as well as the node loop. Checking only between
    // nodes is no limit at all once a single solve outlives the whole budget.
    let deadline = options.time_limit.map(|limit| started + limit);
    lp.set_deadline(deadline);
    // Cuts add rows but never columns, so `n` stays valid across the cut loop.

    let mut nodes = 0usize;
    let mut iterations = 0usize;
    // Everything below is in the internal minimization form; conversion to the
    // user's sense happens once, on the way out.
    let mut incumbent = f64::INFINITY;
    let mut incumbent_x: Option<Vec<f64>> = None;

    let root = lp.solve_with_limit(options.max_iterations_per_node);
    iterations += root.iterations;
    nodes += 1;

    if root.status != LpStatus::Optimal {
        // The LP reports only that it gave up, not why. If the budget has run out,
        // that is what stopped it, and saying NodeLimit would send the caller
        // looking for a node limit they never set.
        let status = match root.status {
            LpStatus::Infeasible => Status::Infeasible,
            _ if deadline.is_some_and(|d| Instant::now() >= d) => Status::TimeLimit,
            _ => Status::NodeLimit,
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

    // The model without any cut rows. Each round rebuilds from here rather than
    // appending to the previous round's model, which is what lets a cut leave again.
    let base = problem.clone();
    // Cuts currently carried in the model, each with a count of consecutive resolves
    // it has sat slack through.
    let mut active: Vec<(cuts::Cut, u32)> = Vec::new();

    for _ in 0..options.cut_rounds {
        if deadline.is_some_and(|d| Instant::now() >= d) {
            break;
        }
        if root.status != LpStatus::Optimal
            || integral_solution(problem, &root.x, options.integrality_tolerance).is_some()
        {
            break;
        }
        // Two families with different reach: covers need a row that reads as a
        // knapsack, while GMI comes off the tableau and applies to any fractional
        // basic column. On dense random rows the second is usually the only one that
        // finds anything.
        let mut found = cuts::separate_until(problem, &root.x, options.cuts_per_round, deadline);
        found.extend(cuts::separate_gomory(
            &lp,
            &root.basis,
            &root.x,
            options.cuts_per_round,
        ));
        // Separation is deliberately generous; this is where the model's row count is
        // actually decided. Ranking by efficacy and dropping near-parallel duplicates
        // keeps the rows that move the bound and discards the ones that only cost an
        // LP column scan at every node below.
        let found = cuts::select(found, &root.x, options.cuts_per_round);
        if found.is_empty() {
            break;
        }
        cuts_added += found.len();
        active.extend(found.into_iter().map(|c| (c, 0)));

        with_cuts = base.clone();
        let rows: Vec<cuts::Cut> = active.iter().map(|(c, _)| c.clone()).collect();
        with_cuts.add_cuts(&rows);

        let mut candidate = Lp::relaxation(&with_cuts);
        candidate.set_deadline(deadline);
        let resolved = candidate.solve_with_limit(options.max_iterations_per_node);
        iterations += resolved.iterations;
        if resolved.status != LpStatus::Optimal {
            // Keep the model that is known to solve rather than pressing on with one
            // that does not; the bound already gained is still sound.
            break;
        }

        // A cut that has gone slack is no longer shaping the relaxation, but it still
        // costs work in every LP below it. Later rounds routinely make earlier rounds'
        // cuts redundant, so age them out rather than accumulating.
        for (cut, age) in &mut active {
            if cut.is_tight(&resolved.x, CUT_SLACK_TOLERANCE) {
                *age = 0;
            } else {
                *age += 1;
            }
        }
        active.retain(|&(_, age)| age < CUT_MAX_AGE);

        reduced = with_cuts;
        problem = &reduced;
        lp = candidate;
        root = resolved;
    }

    // One last purge before the tree opens. Every remaining slack cut would otherwise
    // be carried through every node LP for the rest of the search, and a row that is
    // inactive at the optimum of a convex program cannot be holding the bound up:
    // dropping it leaves the same point optimal. The re-solve is a guard against that
    // reasoning meeting a degenerate basis, not an expectation of change.
    if root.status == LpStatus::Optimal
        && active
            .iter()
            .any(|(c, _)| !c.is_tight(&root.x, CUT_SLACK_TOLERANCE))
    {
        let rows: Vec<cuts::Cut> = active
            .iter()
            .filter(|(c, _)| c.is_tight(&root.x, CUT_SLACK_TOLERANCE))
            .map(|(c, _)| c.clone())
            .collect();
        let mut trimmed = base.clone();
        trimmed.add_cuts(&rows);
        let mut candidate = Lp::relaxation(&trimmed);
        candidate.set_deadline(deadline);
        let resolved = candidate.solve_with_limit(options.max_iterations_per_node);
        iterations += resolved.iterations;
        // Keep the trimmed model only if it really did hold the bound.
        if resolved.status == LpStatus::Optimal
            && resolved.objective >= root.objective - 1e-9 * root.objective.abs().max(1.0)
        {
            reduced = trimmed;
            problem = &reduced;
            lp = candidate;
            root = resolved;
        }
    }

    let root_bound = root.objective;
    let mut open = OpenNodes::new(options.plunge_limit);
    open.push(Node {
        fixings: Vec::new(),
        bound: root_bound,
        basis: root.basis.clone(),
        origin: None,
    });
    let mut heuristic_solutions = 0usize;
    // Nodes whose LP never resolved. Any of them could hold the optimum.
    let mut unresolved = 0usize;

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
        // Whatever produced the point, its continuous columns are still sitting where
        // the relaxation left them. Re-optimizing them is one LP and is often worth
        // far more than the choice of heuristic that found the integers.
        let found = found.map(|found| {
            let polished = heuristic::polish(
                problem,
                &mut lp,
                &root.basis,
                &found.x,
                &options.heuristic_limits,
                &mut iterations,
            );
            match polished {
                Some(better) if better.objective < found.objective => better,
                _ => found,
            }
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
    if let Some(x) = integral_solution(problem, &root.x, options.integrality_tolerance) {
        // Scored from the rounded point, not the relaxation value; see
        // `objective_at`.
        incumbent = objective_at(problem, &x);
        incumbent_x = Some(x);
        open = OpenNodes::new(options.plunge_limit);
    }

    let mut status = Status::Optimal;

    let threads = options.threads.max(1);
    if threads > 1 {
        let result = run_parallel(
            problem,
            &lp,
            open,
            options,
            started,
            incumbent,
            incumbent_x,
            threads,
        );
        nodes += result.nodes;
        iterations += result.iterations;
        heuristic_solutions += result.heuristic_hits;
        let internal_bound = if result.open_bound.is_finite() {
            result.open_bound.min(result.incumbent)
        } else if result.incumbent.is_finite() {
            result.incumbent
        } else {
            root_bound
        };
        let status = match (result.status, &result.incumbent_x) {
            (Status::Optimal, None) => Status::Infeasible,
            (other, _) => other,
        };
        return Solution {
            status,
            objective: result
                .incumbent_x
                .as_ref()
                .map(|_| problem.objective_value(result.incumbent)),
            x: result.incumbent_x.unwrap_or_default(),
            bound: problem.objective_value(internal_bound),
            nodes,
            simplex_iterations: iterations,
            presolve: presolve_stats,
            cuts_added,
            heuristic_solutions,
            root_bound: problem.objective_value(first_bound),
            root_bound_after_cuts: problem.objective_value(root_bound),
        };
    }

    let mut worker = Worker::new(problem, lp, options);

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
        // The parent's bound may have been overtaken by an incumbent found since
        // this node was pushed, in which case it needs no solve at all.
        if !improves(node.bound, incumbent, options.gap_tolerance) {
            continue;
        }

        nodes += 1;
        let outcome = worker.process(&node, incumbent, nodes);
        if outcome.exhausted {
            // This node's LP ran out of iterations, so its subtree was never
            // examined and nothing can be concluded about what is in it. That
            // forfeits the *optimality claim*, but abandoning the whole search over
            // it is far worse: on MIPLIB's pk1 the run stopped after four nodes with
            // an incumbent of 35, where continuing reaches 21.
            unresolved += 1;
            continue;
        }
        heuristic_solutions += outcome.heuristic_hits;
        if let Some((objective, x)) = outcome.incumbent
            && objective < incumbent
        {
            incumbent = objective;
            incumbent_x = Some(x);
        }

        // Improvement runs only where it can: it needs an incumbent to improve on and
        // a relaxation to compare it against, and this node's own relaxation is the
        // freshest one available.
        if options.improvement_frequency > 0
            && nodes.is_multiple_of(options.improvement_frequency)
            && let Some(current) = &incumbent_x
            && let Some(relaxation) = &outcome.relaxation
            && let Some((objective, x)) = improve(problem, current, incumbent, relaxation, &options)
        {
            incumbent = objective;
            incumbent_x = Some(x);
            heuristic_solutions += 1;
        }
        for child in outcome.children {
            open.push(child);
        }
    }
    iterations += worker.iterations;

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
        // A search that exhausted its tree but skipped a node has not proven
        // anything, and must not say it has.
        (Status::Optimal, _) if unresolved > 0 => Status::NodeLimit,
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

/// The internal-form objective of an assignment.
/// Look for a better incumbent in the neighbourhood the incumbent and the relaxation
/// agree on.
///
/// Every integer column where the two already agree is fixed, and the search is turned
/// loose on what is left with a node budget and the incumbent as a cutoff. The two
/// points agreeing on a column is weak evidence that a good solution has it there, and
/// weak evidence over hundreds of columns leaves a model small enough to search
/// properly. This is Danna, Rothberg and Le Pape's RINS.
///
/// Returns a strictly better point, in internal minimization form, or `None`.
///
/// The sub-search runs with improvement off, which is what stops this recursing: a
/// neighbourhood of a neighbourhood is the same idea applied to less and less, and the
/// budget is better spent on the original.
fn improve(
    problem: &Problem,
    incumbent: &[f64],
    incumbent_objective: f64,
    relaxation: &[f64],
    options: &Options,
) -> Option<(f64, Vec<f64>)> {
    let mut neighbourhood = problem.clone();
    let mut fixed = 0usize;
    for j in problem.integer_columns() {
        if (incumbent[j] - relaxation[j]).abs() <= options.integrality_tolerance {
            let value = incumbent[j];
            neighbourhood.col_lb[j] = value;
            neighbourhood.col_ub[j] = value;
            fixed += 1;
        }
    }
    // Nothing agreed, or everything did. Neither leaves a model worth searching: the
    // first is the original problem again, the second is the incumbent again.
    let integers = problem.integer_columns().count();
    if fixed == 0 || fixed == integers {
        return None;
    }

    let found = solve(
        &neighbourhood,
        Options {
            improvement_frequency: 0,
            max_nodes: options.improvement_nodes,
            threads: 1,
            // The sub-search inherits the deadline, so a budget spent here is a budget
            // taken from the search that asked for it, not added to the run.
            time_limit: options.time_limit,
            ..*options
        },
    );

    let x = found.x;
    if x.len() != problem.n_cols() {
        return None;
    }
    let objective = objective_at(problem, &x);
    (objective < incumbent_objective).then_some((objective, x))
}

fn objective_at(problem: &Problem, x: &[f64]) -> f64 {
    problem.obj.iter().zip(x).map(|(c, &v)| c * v).sum()
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

/// Snap a relaxation to an integral assignment, or `None` if some integer column is
/// still fractional.
///
/// Continuous columns pass through untouched: a MIP solution is integral in its
/// integer columns only, and rounding a continuous one would leave the point
/// infeasible.
fn integral_solution(problem: &Problem, x: &[f64], tolerance: f64) -> Option<Vec<f64>> {
    let mut out = Vec::with_capacity(x.len());
    for (j, &v) in x.iter().enumerate() {
        if problem.is_integer(j) {
            let rounded = v.round();
            if (v - rounded).abs() > tolerance {
                return None;
            }
            out.push(rounded);
        } else {
            out.push(v);
        }
    }
    Some(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lp::Lp;
    use crate::model::{RowSense, Sense};
    use crate::sparse::SparseMatrix;

    /// Maximize `x0 + 2 x1 + 3 x2` subject to `x0 + x1 + x2 <= 2`, all binary.
    ///
    /// The optimum takes `x1` and `x2`, for 5. An incumbent that took `x0` and `x1`,
    /// for 3, is the sort of thing a rounding heuristic produces.
    fn knapsackish() -> Problem {
        let mut model = crate::model::Builder::new(Sense::Maximize);
        let x0 = model.binary("x0");
        let x1 = model.binary("x1");
        let x2 = model.binary("x2");
        model.objective(&[(x0, 1.0), (x1, 2.0), (x2, 3.0)]);
        model.row(&[(x0, 1.0), (x1, 1.0), (x2, 1.0)], RowSense::Le, 2.0);
        model.build()
    }

    #[test]
    fn improvement_searches_where_the_incumbent_and_relaxation_agree() {
        let p = knapsackish();
        let options = Options::default();
        // Internal form is minimization, so the incumbent's 3 is -3 and better is more
        // negative.
        let incumbent = vec![1.0, 1.0, 0.0];
        let value = objective_at(&p, &incumbent);
        assert_eq!(value, -3.0);

        // The relaxation agrees on `x0` and disagrees elsewhere, so `x0` is fixed and
        // the rest is searched. With `x0` at one the best remaining is `x2`, for 4.
        let relaxation = vec![1.0, 0.5, 1.0];
        let (improved, x) = improve(&p, &incumbent, value, &relaxation, &options)
            .expect("a better point exists in that neighbourhood");
        assert_eq!(improved, -4.0);
        assert_eq!(x[0], 1.0, "the agreed column stays where both put it");
        assert!(improved < value);
    }

    #[test]
    fn improvement_declines_when_the_neighbourhood_is_the_whole_problem_or_none_of_it() {
        let p = knapsackish();
        let options = Options::default();
        let incumbent = vec![1.0, 1.0, 0.0];
        let value = objective_at(&p, &incumbent);

        // Agreeing on nothing leaves the original problem, which the search is already
        // doing.
        let disagrees = vec![0.0, 0.0, 1.0];
        assert!(improve(&p, &incumbent, value, &disagrees, &options).is_none());

        // Agreeing on everything leaves the incumbent, which cannot be improved on.
        assert!(improve(&p, &incumbent, value, &incumbent, &options).is_none());
    }

    fn node(bound: f64, depth: usize) -> Node {
        // A node's basis is irrelevant to the ordering under test.
        let mut lp = Lp::relaxation(&trivial());
        Node {
            fixings: (0..depth).map(|j| (j as u32, 0.0, 0.0)).collect(),
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
            col_type: vec![crate::model::VarType::Integer],
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
