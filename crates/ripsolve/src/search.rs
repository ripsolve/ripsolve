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
    /// Restart a stalled root relaxation from a perturbed one.
    pub perturb_stalled_root: bool,
    /// Narrow columns the root's bound and the incumbent together already decide.
    ///
    /// A toggle because the claim it makes is a proof about the whole tree, and a proof
    /// is worth being able to turn off and compare against: the test that keeps it
    /// honest solves every sample both ways and requires the same answer.
    pub fix_by_reduced_cost: bool,
    /// What percentage of a pass's free integer columns must become decided before the
    /// solve starts over on the model those fixings leave behind. Zero never restarts.
    pub restart_share: usize,
    /// Flips the LP-free feasibility search may make per column, per attempt, before
    /// giving up.
    ///
    /// Per column because a flip is a column: an absolute count is a different budget
    /// on every model, and on the large ones it is no budget at all. What bounds a run
    /// that is going nowhere is the stall cutoff inside the search, and what bounds the
    /// whole of it is the deadline.
    pub jump_moves: usize,
    /// Rounds of cut separation at the root. Zero disables cuts.
    ///
    /// This defaulted to zero for a long time, on the strength of eleven generated
    /// models over which cutting was slower on every one. That measurement still
    /// holds, and it is still the wrong default, because the models it was taken over
    /// are not the ones the answer turns on.
    ///
    /// On the generated knapsack and covering instances cutting costs about twice the
    /// wall clock and buys nothing, since both configurations solve them either way:
    ///
    /// ```text
    ///                          no cuts    50 rounds of 64
    ///     mkp_200               12.7s         27.5s
    ///     covering_c60_r80_s3    8.2ms        19.2ms
    ///     mkp_500               0.183% gap    0.194% gap
    /// ```
    ///
    /// On MIPLIB the same comparison decides whether an instance is solved at all:
    ///
    /// ```text
    ///                          no cuts    50 rounds of 64
    ///     nexp-50-20-1-1        47.3% gap     optimal in 1.3s
    ///     neos-911970           72.8% gap      6.2% gap
    ///     beavma                36.3% gap      8.2% gap
    ///     n13-3                 47.0% gap     26.2% gap
    ///     decomp1                1.9s         20.8s
    /// ```
    ///
    /// Twice the time on a model that solves regardless is a smaller loss than never
    /// solving one at all, so the default follows the second table. `decomp1` is the
    /// price, and it is a real one.
    ///
    /// What changed in between is that the cuts got better rather than the reasoning:
    /// mixed-integer rounding with bound substitution reaches structure the earlier
    /// families could not, and selection by efficacy and orthogonality keeps the dense
    /// ones out of the model. The old note said this would become worthwhile with
    /// proper cut selection, and it did.
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
            jump_moves: 200,
            perturb_stalled_root: true,
            fix_by_reduced_cost: true,
            restart_share: 25,
            cut_rounds: 50,
            cuts_per_round: 64,
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
    /// A result carrying no solution: the search ended before it held one.
    ///
    /// Written once because the fields that must be absent are easy to get subtly
    /// wrong when spelled out at each exit, and a bound left over from a previous
    /// phase would be read as a proof of something.
    fn without_solution(
        status: Status,
        nodes: usize,
        simplex_iterations: usize,
        presolve: Option<presolve::Stats>,
    ) -> Self {
        Self {
            status,
            objective: None,
            x: Vec::new(),
            bound: f64::NAN,
            nodes,
            simplex_iterations,
            presolve,
            cuts_added: 0,
            heuristic_solutions: 0,
            root_bound: f64::NAN,
            root_bound_after_cuts: f64::NAN,
        }
    }

    /// Remaining optimality gap, relative to the incumbent. Zero when proven.
    ///
    /// A run that never got a bound reports no gap rather than a gap of zero. The
    /// bound is deliberately `NaN` where nothing was proven, so that a leftover from an
    /// earlier phase cannot be read as a proof, and `f64::max` returns its non-`NaN`
    /// argument: clamping the gap at zero turned that `NaN` into exactly the "proven
    /// optimal" reading the `NaN` was there to prevent. A run of `neos-954925` whose
    /// root relaxation never finished reported `gap 0.0000%`.
    pub fn gap(&self) -> f64 {
        match self.objective {
            Some(obj) if obj.is_finite() && self.bound.is_finite() => {
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
    /// Bounds the root's reduced costs imply once the incumbent is good enough, and how
    /// many of them the last incumbent this worker saw had earned. Recomputed when that
    /// incumbent moves, which is rarely, rather than per node.
    lurking: &'a Lurking,
    in_force: usize,
    priced_at: f64,
    /// The spacing of the objective values this model can take; see
    /// [`objective_granularity`].
    granularity: Option<f64>,
    /// Integer columns this pass arrived with still free, which is what a restart's
    /// worth is measured against.
    free_at_entry: usize,
}

impl<'a> Worker<'a> {
    /// Has the incumbent decided enough of the model to be worth starting over on?
    ///
    /// The lurking table already knows: the entries in force at a cutoff are a prefix,
    /// so the count is a binary search and costs nothing to ask on every node. What it
    /// measures is what a fresh pass would find *already fixed* before it did anything,
    /// which is exactly what a restart is buying.
    ///
    /// Against the columns this pass started with, not against the original model. A
    /// pass that arrives on a model three quarters fixed and fixes another tenth has
    /// found a tenth, and it is that tenth which decides whether presolving and cutting
    /// again is worth the time it takes.
    fn wants_restart(&self, incumbent: f64, options: &Options) -> bool {
        if self.lurking.is_empty() || !incumbent.is_finite() {
            return false;
        }
        let cutoff = match self.granularity {
            Some(g) => lift_to_granularity(incumbent, self.granularity) - g,
            None => incumbent,
        };
        let free = self.free_at_entry;
        if free == 0 {
            return false;
        }
        let entries = &self.lurking.entries[..self.lurking.in_force(cutoff)];
        let decided = entries
            .iter()
            .filter(|entry| {
                let j = entry.column as usize;
                entry.lo >= entry.hi && self.problem.col_lb[j] < self.problem.col_ub[j]
            })
            .count();
        decided * 100 >= free * options.restart_share
    }

    fn new(problem: &'a Problem, lp: Lp, options: Options, lurking: &'a Lurking) -> Self {
        let n = problem.n_cols();
        Self {
            problem,
            lp,
            pseudocosts: Pseudocosts::new(n),
            strong_budget: options.strong_branching_budget,
            options,
            iterations: 0,
            dives: Schedule::new(options.heuristic_frequency),
            lurking,
            in_force: 0,
            priced_at: f64::INFINITY,
            granularity: objective_granularity(problem),
            free_at_entry: problem
                .integer_columns()
                .filter(|&j| problem.col_lb[j] < problem.col_ub[j])
                .count(),
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
        // Then whatever the incumbent has earned since the root, which is a statement
        // about the whole tree and so belongs before this node's own branching.
        if !self.lurking.is_empty() && incumbent != self.priced_at {
            let cutoff = match self.granularity {
                Some(g) => lift_to_granularity(incumbent, self.granularity) - g,
                None => incumbent,
            };
            self.in_force = self.lurking.in_force(cutoff);
            self.priced_at = incumbent;
        }
        for entry in &self.lurking.entries[..self.in_force] {
            self.lp
                .set_column_bounds(entry.column as usize, entry.lo, entry.hi);
        }
        // The branch last, and intersected rather than assigned: a branch that asks for
        // a value the incumbent has already ruled out leaves an empty range, and an
        // empty range is a node with nothing in it rather than a bound to be honoured.
        for &(j, lo, hi) in &node.fixings {
            let (held_lo, held_hi) = self.lp.column_bounds(j as usize);
            let (lo, hi) = (lo.max(held_lo), hi.min(held_hi));
            if lo > hi {
                return out;
            }
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

        // Lifted to the next objective the model can take before it is compared with
        // anything. A relaxation worth 52000.3 on a model whose objectives step by 200
        // is a proof that nothing here beats 52200.
        let relaxed = lift_to_granularity(solved.objective, self.granularity);
        if !improves(relaxed, incumbent, options.gap_tolerance) {
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
        let mut bound = relaxed;
        if options.local_cut_frequency > 0
            && index.is_multiple_of(options.local_cut_frequency)
            && integral_solution(problem, &solved.x, options.integrality_tolerance).is_none()
            && let Some(tightened) = self.separate_locally(&solved, &options, cutoff)
        {
            bound = lift_to_granularity(tightened, self.granularity);
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
    /// The search stopped because starting over is worth more than continuing.
    wants_restart: bool,
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
    /// The weakest bound among the nodes skipped for running out of iterations, which
    /// is what decides whether skipping them forfeited the optimality claim.
    unresolved_bound: AtomicU64,
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
/// The incumbent has decided enough of the model to be worth starting over on.
const STOP_RESTART: usize = 3;

impl Shared {
    fn incumbent(&self) -> f64 {
        f64::from_bits(self.best_bits.load(Ordering::Relaxed))
    }

    /// The incumbent and the assignment that achieves it.
    ///
    /// Takes the lock, unlike [`Shared::incumbent`], because the improvement search
    /// needs the point and not just its value. Called once every few hundred nodes,
    /// so the copy is not on any hot path.
    fn incumbent_solution(&self) -> (f64, Option<Vec<f64>>) {
        let best = self.best.lock().expect("incumbent lock");
        (best.0, best.1.clone())
    }

    /// Remember the weakest bound skipped, so the claim can be settled at the end.
    fn weaken_unresolved(&self, bound: f64) {
        let mut held = self.unresolved_bound.load(Ordering::Relaxed);
        while bound < f64::from_bits(held) {
            match self.unresolved_bound.compare_exchange_weak(
                held,
                bound.to_bits(),
                Ordering::Relaxed,
                Ordering::Relaxed,
            ) {
                Ok(_) => return,
                Err(current) => held = current,
            }
        }
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
    lurking: &Lurking,
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
        unresolved_bound: AtomicU64::new(f64::INFINITY.to_bits()),
    };

    std::thread::scope(|scope| {
        for _ in 0..threads {
            // Each worker gets its own LP, since solving a node mutates the column
            // bounds, and its own branching history. Sharing the pseudocosts would
            // pool more evidence but put a lock on the hot path; per-worker history
            // is the cheaper trade at these thread counts.
            let mut worker = Worker::new(problem, lp.clone(), options, lurking);
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
                    // Whatever the incumbent has proved since this pass began is worth
                    // more spent on a smaller model than on this one; see `solve`.
                    if worker.wants_restart(shared.incumbent(), &options) {
                        shared.stopped.store(STOP_RESTART, Ordering::Relaxed);
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
                        // Skip the node and keep going; see the serial driver. What is
                        // remembered is the node's own bound, which is what decides
                        // later whether skipping it cost anything.
                        shared.unresolved.fetch_add(1, Ordering::Relaxed);
                        shared.weaken_unresolved(node.bound);
                    }
                    if let Some((objective, x)) = outcome.incumbent
                        && shared.offer(objective, x)
                        && outcome.heuristic_hits > 0
                    {
                        shared.heuristic_hits.fetch_add(1, Ordering::Relaxed);
                    }

                    // Improve the incumbent from this node's relaxation, which the
                    // serial driver has always done here and this one never did. The
                    // omission left the parallel search with no way to repair a poor
                    // first incumbent: on MIPLIB's cap6000 at two threads it sat on
                    // -2355245 for the whole minute where the serial search reaches
                    // -2451274 in three seconds and stops.
                    //
                    // `index` is drawn from the shared counter, so this fires once per
                    // `improvement_frequency` nodes across the search as a whole,
                    // exactly as it does serially, rather than once per worker.
                    if let Some(relaxation) = &outcome.relaxation {
                        let (value, current) = shared.incumbent_solution();
                        if let Some(current) = current
                            && let Some((objective, x)) = improve(
                                problem,
                                &current,
                                value,
                                relaxation,
                                &options,
                                remaining_of(options.time_limit, started),
                            )
                            && shared.offer(objective, x)
                        {
                            shared.heuristic_hits.fetch_add(1, Ordering::Relaxed);
                        }
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
        // Not an answer, and not reported as one: the caller starts over and the status
        // that matters is whatever the pass after this one arrives at.
        STOP_RESTART => Status::NodeLimit,
        // A tree exhausted apart from a skipped node proves nothing, unless the
        // skipped node could not have held anything better anyway; see `improves`.
        _ if shared.unresolved.load(Ordering::Relaxed) > 0
            && improves(
                f64::from_bits(shared.unresolved_bound.load(Ordering::Relaxed)),
                incumbent,
                options.gap_tolerance,
            ) =>
        {
            Status::NodeLimit
        }
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
        wants_restart: shared.stopped.load(Ordering::Relaxed) == STOP_RESTART,
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
/// What separating at the root produced.
struct RootCuts {
    /// The model carrying the cuts that survived, when any were kept. `None` means the
    /// caller's model is still the one to search.
    model: Option<Problem>,
    lp: Lp,
    root: LpSolution,
    added: usize,
    iterations: usize,
}

/// Separate cuts at the root, re-solving after each round.
///
/// Cut at the root only. Separating deeper in the tree is stronger still and is what
/// [`Options::local_cut_frequency`] does; those cuts are valid for one subtree, while
/// these are valid everywhere and so are added to the model once.
///
/// Returns the model to search rather than rebinding the caller's, which is what keeps
/// "the presolved model" and "the presolved model plus cuts" from being the same
/// variable.
fn separate_at_root(
    problem: &Problem,
    mut lp: Lp,
    mut root: LpSolution,
    options: &Options,
    deadline: Option<Instant>,
    // When to stop cutting and let the search have the rest, distinct from the run's
    // own deadline; see the caller.
    stop_cutting: Option<Instant>,
) -> RootCuts {
    // Nothing to do, and in particular no reason to clone the model, which is the
    // common case since root cutting is off by default.
    if options.cut_rounds == 0 {
        return RootCuts {
            model: None,
            lp,
            root,
            added: 0,
            iterations: 0,
        };
    }

    let mut added = 0usize;
    let mut iterations = 0usize;
    let mut problem = problem;
    // Cut at the root only. Separating deeper in the tree ("local cuts") is
    // stronger still, but it needs the cut pool to track which nodes each cut is
    // valid for; these are globally valid, so adding them once to the model is
    // both simpler and correct everywhere.

    // The model without any cut rows. Each round rebuilds from here rather than
    // appending to the previous round's model, which is what lets a cut leave again.
    let base = problem.clone();
    // Which columns exclude one another, read once from the model that arrived. Cut
    // rows are added below but never create conflicts: they are implied by the model
    // the graph was built from, so anything they would say is already in it.
    let conflicts = cuts::Conflicts::of(problem);
    let mut model: Option<Problem> = None;
    let mut with_cuts;
    // Cuts currently carried in the model, each with a count of consecutive resolves
    // it has sat slack through.
    let mut active: Vec<(cuts::Cut, u32)> = Vec::new();

    for _ in 0..options.cut_rounds {
        if deadline.is_some_and(|d| Instant::now() >= d)
            || stop_cutting.is_some_and(|d| Instant::now() >= d)
        {
            break;
        }
        if root.status != LpStatus::Optimal
            || integral_solution(problem, &root.x, options.integrality_tolerance).is_some()
        {
            break;
        }
        // Five families with different reach: covers need a row that reads as a
        // knapsack, MIR needs a row mixing integer and continuous columns, cliques need
        // columns that exclude one another, GMI comes off the tableau and applies to any
        // fractional basic column, and mod-2 combines whole rows. On dense random rows the last is usually the
        // only one that finds anything; on mixed models from MIPLIB the second carries
        // most of the bound; and on the pure binary set, where there is no continuous
        // structure for MIR to reach, the third is the one with anything to say.
        let mut found = cuts::separate_until(problem, &root.x, options.cuts_per_round, deadline);
        found.extend(cuts::separate_mir(problem, &root.x, options.cuts_per_round));
        found.extend(cuts::separate_cliques(
            problem,
            &conflicts,
            &root.x,
            options.cuts_per_round,
        ));
        found.extend(cuts::separate_gomory(
            &lp,
            &root.basis,
            &root.x,
            options.cuts_per_round,
        ));
        // Sixth family, and the only one that combines rows. The four above each read a
        // single row or a single tableau row, and on a model whose optimal face is large
        // that removes one vertex of it and the relaxation steps to the next one worth
        // the same. See `separate_mod2`.
        //
        // On a smaller allowance than the rest; see `MOD2_SHARE`.
        found.extend(cuts::separate_mod2(
            problem,
            &root.x,
            options.cuts_per_round / MOD2_SHARE,
        ));
        // Separation is deliberately generous; this is where the model's row count is
        // actually decided. Ranking by efficacy and dropping near-parallel duplicates
        // keeps the rows that move the bound and discards the ones that only cost an
        // LP column scan at every node below.
        let mut found = cuts::select(found, &root.x, options.cuts_per_round);

        // A row aggregated through the variable upper bounds the model implies, and it
        // is selected on its own rather than against the families above. Efficacy
        // divides violation by the cut's norm, which prices an aggregated row of six
        // hundred terms against a two-term clique and loses every time; on
        // `neos-787933` the two-term cliques win that ranking, move the bound from 8.07
        // to 8.16, and crowd out the 77 aggregated rows that move it to 30, which is
        // the optimum. The ranking is not wrong in general -- a short cut really is
        // cheaper to carry -- so this family gets its own budget instead of a thumb on
        // the scale. See `separate_implied_aggregations`.
        found.extend(cuts::select(
            cuts::separate_implied_aggregations(
                problem,
                &conflicts,
                &root.x,
                options.cuts_per_round,
            ),
            &root.x,
            options.cuts_per_round / AGGREGATION_SHARE,
        ));
        if found.is_empty() {
            break;
        }
        added += found.len();
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

        model = Some(with_cuts);
        problem = model.as_ref().expect("just set");
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
            model = Some(trimmed);
            lp = candidate;
            root = resolved;
        }
    }
    RootCuts {
        model,
        lp,
        root,
        added,
        iterations,
    }
}

/// What one pass at a model produced, and whether it is worth starting over.
struct Attempt {
    solution: Solution,
    /// The model to restart on: the one just searched, with everything the incumbent
    /// has since proved about it written in. `None` when nothing was proved, or when
    /// there is no time to spend on a fresh start.
    restart: Option<Problem>,
}

/// Solve, restarting on the model the search has narrowed while it ran.
///
/// A pass fixes columns as the incumbent improves, through reduced cost fixing and the
/// lurking table, but it fixes them *in the tree it is already walking*: the model that
/// was presolved, cut and branched on is the one that arrived, and none of the work that
/// shaped it is redone on the smaller model those fixings leave behind. Starting over is
/// how that work gets redone, and on a model where three quarters of the columns have
/// gone it is a different problem.
pub fn solve(problem: &Problem, options: Options) -> Solution {
    solve_from(problem, options, None)
}

/// Solve, starting from a point already known to be feasible.
///
/// The usual way to hand a solver an answer it should not have to rediscover: a solution
/// from a previous run, from a related model, or from somewhere outside the solver
/// entirely. The point is checked before it is believed, and a point that does not
/// satisfy the model is ignored rather than trusted, because an incumbent is a claim the
/// search prunes against and a wrong one removes the answer.
///
/// It is also how the primal side is measured. Handing the search a known optimum and
/// asking whether it can then *prove* optimality separates an instance that is short of a
/// point from one that is short of a bound, and those want completely different work.
pub fn solve_from(problem: &Problem, options: Options, start: Option<&[f64]>) -> Solution {
    let started = Instant::now();
    let mut model: Option<Problem> = None;
    let mut seed: Option<(f64, Vec<f64>)> = start
        .filter(|x| {
            x.len() == problem.n_cols()
                && heuristic::is_feasible(
                    problem,
                    x,
                    options.heuristic_limits.feasibility_tolerance,
                )
                && problem
                    .integer_columns()
                    .all(|j| (x[j] - x[j].round()).abs() <= options.integrality_tolerance)
        })
        .map(|x| (objective_at(problem, x), x.to_vec()));
    let (mut nodes, mut iterations, mut heuristics) = (0usize, 0usize, 0usize);
    for _ in 0..=MAX_RESTARTS {
        let current = model.as_ref().unwrap_or(problem);
        let attempt = solve_once(current, options, started, seed.take());
        // The counters belong to the whole solve, not to the last pass of it.
        nodes += attempt.solution.nodes;
        iterations += attempt.solution.simplex_iterations;
        heuristics += attempt.solution.heuristic_solutions;
        let mut solution = attempt.solution;
        solution.nodes = nodes;
        solution.simplex_iterations = iterations;
        solution.heuristic_solutions = heuristics;
        let out_of_time = options
            .time_limit
            .is_some_and(|limit| started.elapsed() >= limit);
        let Some(narrowed) = attempt.restart.filter(|_| !out_of_time) else {
            return solution;
        };
        // Nothing carries over but the point, which is still feasible: the narrowed
        // model only ever has tighter bounds, and every fixing written into it was
        // proved against this incumbent or a worse one.
        seed = (!solution.x.is_empty()).then(|| {
            let internal = objective_at(problem, &solution.x);
            (internal, solution.x.clone())
        });
        model = Some(narrowed);
    }
    // Out of restarts. One more pass, and whatever it says is the answer.
    let attempt = solve_once(
        model.as_ref().unwrap_or(problem),
        options,
        started,
        seed,
    );
    let mut solution = attempt.solution;
    solution.nodes += nodes;
    solution.simplex_iterations += iterations;
    solution.heuristic_solutions += heuristics;
    solution
}

fn solve_once(
    problem: &Problem,
    options: Options,
    started: Instant,
    seed: Option<(f64, Vec<f64>)>,
) -> Attempt {

    // Presolve reduces in place and introduces no renumbering, so the reduced
    // model's solution vector is directly the original's, there is no postsolve.
    // It is sound (it never admits a point the original rejects) and preserves the
    // optimum, so searching the reduced model answers the original question.
    let mut reduced;
    let (problem, presolve_stats) = if options.presolve {
        reduced = problem.clone();
        // A backstop, not the budget: what bounds probing is the work it may do per
        // column. This catches a pathological model, not the ordinary case.
        let presolve_deadline = options
            .time_limit
            .map(|limit| started + limit.mul_f64(PRESOLVE_SHARE));
        match presolve::presolve_until(&mut reduced, 20, presolve_deadline) {
            Outcome::Infeasible => {
                return Attempt {
                    solution: Solution::without_solution(Status::Infeasible, 0, 0, None),
                    restart: None,
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
    // A restart begins with the point the pass before it ended on. That point is still
    // feasible for this model, whose bounds are only tighter, and having it from the
    // first node is most of what a restart is for: reduced cost fixing and the lurking
    // table both read the incumbent, and a fresh pass that had to find one again would
    // spend the run re-earning what it already knew.
    let (mut incumbent, mut incumbent_x) = match seed {
        Some((objective, x)) => (objective, Some(x)),
        None => (f64::INFINITY, None),
    };

    // Run only where it is needed, which is the two places nothing else reaches: when
    // the relaxation does not finish, and when it does and every heuristic built on it
    // comes back empty.
    //
    // Not before the relaxation unconditionally, which is what it did first. The models
    // it cannot crack take the whole budget it is given and return nothing, and most
    // models are ones it is never needed on: paying for it up front took `f2gap401600`
    // from 0.27 seconds to 11.5 and `mod010` from 0.78 to 12.9, on instances whose own
    // heuristics find a point in under a second. Asked only after those have failed, it
    // costs nothing on either.
    let jump = |from: Option<&[f64]>, elapsed_share: f64| -> Option<heuristic::Incumbent> {
        if options.heuristic_frequency == 0 {
            return None;
        }
        let start: Vec<f64> = (0..problem.n_cols())
            .map(|j| match from {
                Some(x) => x[j].clamp(problem.col_lb[j], problem.col_ub[j]),
                // No relaxation to round, so start from whichever bound is nearer zero
                // and let the weights do the rest.
                None if problem.col_lb[j] > 0.0 => problem.col_lb[j],
                None => 0.0f64.clamp(problem.col_lb[j], problem.col_ub[j]),
            })
            .collect();
        // From now rather than from the start of the run. This is asked for late, after
        // the relaxation and the cutting and the whole chain above it have had their
        // turn, and a budget measured from the beginning of a minute has already been
        // spent by the time it is reached: every instance this heuristic wins came back
        // empty when the deadline was anchored that way, having been given no time at
        // all rather than too little.
        let jump_deadline = options
            .time_limit
            .map(|limit| Instant::now() + limit.mul_f64(elapsed_share))
            .or(deadline);
        // Never past the run's own end.
        let jump_deadline = match (jump_deadline, deadline) {
            (Some(own), Some(hard)) => Some(own.min(hard)),
            (own, hard) => own.or(hard),
        };
        // Restarted from randomised corners when the first attempt fails, which is what
        // the search stopping is for: a run that has settled has settled, and the
        // weights that got it there are the reason it will not move again. Where it
        // begins decides which local minimum it has to climb out of, and it is cheaper
        // to begin somewhere else than to keep climbing.
        //
        // `neos-3226448-wkra` is the case: nothing from the model's own bounds however
        // long it is given, and a feasible point on the seventeenth random start. The
        // budget is the same either way, so a model this does not help pays only the
        // scheduling.
        let mut state = 0x9E37_79B9_7F4A_7C15u64;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        for attempt in 0..JUMP_RESTARTS {
            let from: Vec<f64> = if attempt == 0 {
                start.clone()
            } else {
                (0..problem.n_cols())
                    .map(|j| {
                        let (lo, hi) = (problem.col_lb[j], problem.col_ub[j]);
                        let pick = if next() & 1 == 0 { lo } else { hi };
                        if pick.is_finite() { pick } else { start[j] }
                    })
                    .collect()
            };
            let found = heuristic::feasibility_jump(
                problem,
                &from,
                &options.heuristic_limits,
                // Per attempt, and per *column*: an absolute count of flips is a
                // budget of two and a half flips a column on a model of ten thousand
                // and two hundred and fifty on a model of a hundred, which is the same
                // mistake this file complains about elsewhere. What actually bounds a
                // run that is going nowhere is the stall cutoff inside the search, and
                // what bounds the whole of it is the deadline; this is the third guard
                // and it should not be the binding one. `neos-3226448-wkra` closes at
                // two hundred flips a column and does not at two and a half.
                options.jump_moves.saturating_mul(problem.n_cols().max(1)),
                jump_deadline,
            );
            if found.is_some() {
                return found;
            }
            if jump_deadline.is_some_and(|d| Instant::now() >= d) {
                break;
            }
        }
        None
    };

    // A model with nothing to optimise is asking for a feasible point, and the whole
    // apparatus above the feasibility search exists to find a *good* point. So it goes
    // first here and gets the run rather than a twentieth of it, because there is
    // nothing else for the run to be spent on: the bound is a constant that needs no
    // relaxation to establish, and any point at all attains it.
    //
    // `neos-3226448-wkra` is the case. Its objective is empty, HiGHS answers it in two
    // tenths of a second without solving a single LP, and this solver spent sixty
    // seconds on a relaxation whose value it could have written down.
    if objective_is_constant(problem)
        && let Some(found) = jump(None, CONSTANT_OBJECTIVE_SHARE)
    {
        let value = problem.objective_value(constant_objective(problem));
        let mut solution =
            Solution::without_solution(Status::Optimal, nodes, iterations, presolve_stats);
        solution.objective = Some(value);
        solution.bound = value;
        solution.root_bound = value;
        solution.root_bound_after_cuts = value;
        solution.x = found.x;
        solution.heuristic_solutions = 1;
        return Attempt {
            solution,
            restart: None,
        };
    }

    // The relaxation gets a first attempt on a share of the run. Most models finish it
    // and pay nothing more.
    if let Some(limit) = options.time_limit {
        let first = started + limit.mul_f64(ROOT_LP_FIRST_SHARE);
        lp.set_deadline(Some(match deadline {
            Some(hard) => hard.min(first),
            None => first,
        }));
    }
    let mut root = lp.solve_with_limit(options.max_iterations_per_node);
    lp.set_deadline(deadline);
    iterations += root.iterations;

    // A relaxation that did not finish is usually not slow, it is stuck. Two rescues,
    // each bounded by what the first attempt spent: the same model at the same size is
    // worth about as much again, and no more. Held to the caller's deadline instead,
    // the first of them took everything that was left, and on all eight of the models
    // that need a rescue here it came back at IterationLimit with the whole run gone.
    // The second rescue was unreachable for that reason alone.
    let budget = root.iterations.clamp(1, options.max_iterations_per_node);

    // The dual method first, because it is cheap where it works: `air04`'s relaxation
    // goes from not finishing in 130 seconds to finishing in 1.5, and `tanglegram6`'s
    // in 0.5. Where it does not work it is bounded like everything else here.
    //
    // Asked here and nowhere else, and that is the whole of why it is safe. The primal
    // method reaches feasibility through phase 1, whose cost vector scores a basic
    // variable by whether it currently violates a bound and so is blind to the kink at
    // one sitting exactly on one; no ratio test repairs a column choice made that way,
    // and the dual method has no phase 1 to get stuck in. Made the default it loses
    // anyway: the two methods end on *different* optimal vertices, and everything
    // downstream reads the root's vertex rather than its objective, since Gomory cuts
    // come off that tableau and branching reads its fractional values. On this set the
    // vertex it reaches is usually the worse one to start from. Asked only where the
    // primal method produced no vertex at all, there is nothing to be worse than.
    if root.status != LpStatus::Optimal && deadline.is_none_or(|d| Instant::now() < d) {
        let rescued = lp.solve_cold_dual(budget);
        iterations += rescued.iterations;
        if rescued.status == LpStatus::Optimal {
            root = rescued;
        }
    }

    // Then perturbation. The models this catches are set partitioning problems whose
    // every coefficient is one, where a crowd of bases describes one point and the
    // steps between them have length zero: `ex9` spends 8871 consecutive iterations
    // taking steps of length zero without the worst violation moving off 1.0.
    //
    // Perturbing the bounds by a random amount too small to matter breaks that: no two
    // variables sit on the same bound any more, so the steps have somewhere to go. What
    // comes back is the wrong problem's answer and the right problem's *basis*, which is
    // the part worth having. Restoring the true bounds and re-solving from there gives
    // the true optimum, and gives it quickly, because that basis is optimal for a
    // problem next to this one rather than stuck inside it.
    //
    // Measured on the relaxations that stall: `neos-1324574` from 212471 iterations and
    // 214 seconds to 5, and `tanglegram6` from not finishing at all to 172. This is the
    // same idea as an earlier attempt that was reverted, and differs in the one way that
    // matters: that one perturbed and warm started from the *stalled* basis, which puts
    // the search back in the degeneracy it was trying to leave.
    if root.status != LpStatus::Optimal
        && deadline.is_none_or(|d| Instant::now() < d)
        && options.perturb_stalled_root
        && let Some(rescued) = perturbed_root(problem, &lp, deadline, budget)
    {
        iterations += rescued.iterations;
        root = rescued;
    }
    // Both rescues bounded and both failed, with the caller's clock still running. The
    // budget above is what makes a *second* attempt worth trying, not a ceiling on the
    // relaxation itself: what the caller asked for is the best answer available by their
    // deadline, so the last resort is the first attempt continued to it. Without this a
    // model whose rescues are cheap and hopeless returns at 41 seconds of a 60 second
    // budget, having given up with a third of the run unspent.
    if root.status != LpStatus::Optimal && deadline.is_none_or(|d| Instant::now() < d) {
        let again = lp.solve_with_limit(options.max_iterations_per_node);
        iterations += again.iterations;
        if again.status == LpStatus::Optimal {
            root = again;
        }
    }
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
        // A point found before the relaxation survives the relaxation failing. It is
        // not proven optimal and never claims to be, but reporting nothing when
        // something feasible is in hand throws away the only result there was.
        let rescued = if status == Status::Infeasible {
            None
        } else {
            jump(None, JUMP_SHARE)
        };
        let solution = match rescued {
            Some(found) if status != Status::Infeasible => {
                let mut solution =
                    Solution::without_solution(status, nodes, iterations, presolve_stats);
                solution.objective = Some(problem.objective_value(found.objective));
                solution.x = found.x;
                solution.heuristic_solutions = 1;
                solution
            }
            _ => Solution::without_solution(status, nodes, iterations, presolve_stats),
        };
        return Attempt {
            solution,
            restart: None,
        };
    }

    let first_bound = lift_to_granularity(root.objective, objective_granularity(problem));

    // Cutting gets a share of the run rather than the run. Each round re-solves the
    // whole model, so on a large one the loop can spend everything it is given: on
    // MIPLIB's mitre the search reaches one node with fifty rounds of cutting and 5126
    // with none, the difference being that the rounds consume the minute before the
    // tree is ever entered. A bound that is never spent against is worth less than no
    // bound and a search.
    //
    // A share rather than a fixed count because the cost of a round is the cost of an
    // LP solve, which is what varies between the models where fifty rounds are free and
    // the models where five are not.
    let stop_cutting = options
        .time_limit
        .map(|limit| started + limit.mul_f64(ROOT_CUT_SHARE));
    let cut = separate_at_root(problem, lp, root, &options, deadline, stop_cutting);
    iterations += cut.iterations;
    let cuts_added = cut.added;
    let mut lp = cut.lp;
    let root = cut.root;
    // Held here so the borrow below outlives the branch that produced it.
    let with_cuts_model = cut.model;
    let problem = with_cuts_model.as_ref().unwrap_or(problem);

    let granularity = objective_granularity(problem);
    let root_bound = lift_to_granularity(root.objective, granularity);
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
    let mut unresolved_bound = f64::INFINITY;

    // An incumbent before the first branch is worth more than one found later: the
    // search cannot prune anything until it holds one.
    if options.heuristic_frequency > 0 && root.status == LpStatus::Optimal {
        // Rounding costs no LP at all, diving a short chain of them, the pump the
        // most. Fixing with propagation costs no LP either and still comes last: it is
        // ordered by what its points are worth rather than by what they cost.
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
            })
            .or_else(|| {
                // Last rather than first despite being the cheapest of them. It finds
                // points the others cannot, and finds worse ones than they do where
                // both succeed: put ahead of diving it took `cap6000` from 500 nodes to
                // 1500, having installed a weaker first incumbent and left the search
                // less to prune with. Where the chain above it works there is nothing
                // here worth having, and where it does not this is the whole of what
                // there is.
                let conflicts = cuts::Conflicts::of(problem);
                heuristic::fix_and_propagate(
                    problem,
                    &conflicts,
                    &root.x,
                    &options.heuristic_limits,
                )
            })
            // The point found before the relaxation, for the same reason and with the
            // same evidence: installed ahead of the chain rather than behind it, it
            // took `eil33-2` from solved in 96 seconds to unsolved in 150, having been
            // good enough that nothing above replaced it and worse than what the chain
            // would have reached on its own. Every instance it wins is one where
            // everything above returns nothing, so nothing is given up by asking it
            // last.
            // The point found before the relaxation, kept only when the relaxation
            // agrees it is a good one. Where this heuristic wins it wins outright, its
            // point landing on the optimum with the bound already there to prove it:
            // `acc-tight2`, `disctom` and `neos-913984` all close in two nodes from a
            // point at zero gap.
            //
            // Its point used to be thrown away when it sat more than a tenth off the
            // root bound, because a poor one had cost `eil33-2` its solve, 96 seconds
            // to unsolved in 150. That is no longer true of any instance in the set,
            // and a point is worth more than it was when the threshold was measured:
            // reduced cost fixing reads the incumbent, so a poor point that prunes
            // nothing by itself still decides columns. Every one of these arrives where
            // the search would otherwise report no incumbent at all, and accepting
            // them all costs none of the 31 its solve while `neos-820879` goes from no
            // bound to a gap of 1.27%.
            .or_else(|| {
                // From the bounds rather than from the relaxation, which is measured
                // rather than assumed: rounding the relaxation looks like the better
                // start and is not, losing `acc-tight2`, `disctom` and `neos-913984`
                // outright. The weights are what find the point, and where they start
                // from decides which local minimum they have to climb out of first.
                jump(None, JUMP_SHARE)
            })
            // Last of all, and cheapest of all. A corner of the box is usually a poor
            // point and is not offered ahead of anything: put in front of diving, a
            // cheap heuristic whose points are poor takes the search's incumbent away
            // from a better one, which is what a corner point mostly is. Reached only
            // when everything above has failed, it is the difference between reporting
            // a feasible solution and reporting none, and no quality bar is applied
            // because at that point there is nothing to compare it against.
            .or_else(|| heuristic::corners(problem, &options.heuristic_limits));
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
        //
        // Compared rather than assigned. This used to take the relaxation's point
        // outright, which was safe while every solve began with no incumbent at all and
        // is not now: a restart begins holding the point the pass before it found, and
        // an integral relaxation on the narrowed model can be worse than it. On
        // `samples/simple.lp` it is — the restarted root is integral at 7 against a
        // seed of 10 — and assigning threw the answer away.
        let objective = objective_at(problem, &x);
        if objective < incumbent {
            incumbent = objective;
            incumbent_x = Some(x);
            open = OpenNodes::new(options.plunge_limit);
        }
    }

    // With a proven bound and a point in hand, some columns are already decided.
    //
    // A nonbasic column's reduced cost is the rate at which the objective rises as it
    // leaves the bound it is parked on, so if travelling to its other bound would push
    // the root's own bound past the incumbent, no solution better than the one already
    // held has it there. That is a proof about the whole tree, made once, from numbers
    // the root LP has already computed.
    //
    // The bounds go into the model rather than into a node, because every node rebuilds
    // its bounds from the model before solving; one pass here reaches all of them.
    // What a *better* solution has to reach, which on a model whose objectives step by
    // 200 is 200 below the incumbent rather than a whisker below it. That is the whole
    // of the room reduced cost fixing has to work in, so a step's worth of it matters:
    // on `n2seq36f` the room goes from 200 to nothing at all, and a room of nothing
    // pins every column the relaxation left on a bound.
    let cutoff = |incumbent: f64| match granularity {
        Some(g) => lift_to_granularity(incumbent, granularity) - g,
        None => incumbent,
    };
    let tightened =
        (options.fix_by_reduced_cost && root.status == LpStatus::Optimal && incumbent.is_finite())
            .then(|| fix_by_reduced_cost(problem, &lp, &root, cutoff(incumbent), &options))
            .flatten();
    let problem = tightened.as_ref().unwrap_or(problem);

    // The same reasoning carried forward. The one pass above spends the incumbent the
    // root heuristics happened to find, which on this set is far weaker than the one the
    // search ends with -- `n2seq36f` is at a 39.7% gap here and 0.38% at the end -- so
    // most of what these reduced costs prove is not provable yet. Working out in advance
    // what each of them will prove, and at which incumbent, costs one pass and lets the
    // search collect the rest as it earns it.
    let lurking = if options.fix_by_reduced_cost && root.status == LpStatus::Optimal {
        lp.reduced_costs(&root.basis)
            .map(|costs| {
                Lurking::build(
                    problem,
                    &costs,
                    root.objective,
                    options.integrality_tolerance,
                )
            })
            .unwrap_or(Lurking {
                entries: Vec::new(),
            })
    } else {
        Lurking {
            entries: Vec::new(),
        }
    };

    // With the table in hand, the cheapest thing to try next is what it would prove if
    // the answer were better than anything found so far; see
    // `reduced_cost_neighbourhood`. This is where the instances whose bound is already
    // exact are waiting: they need a point, not a bound.
    // Not when the incumbent already matches the bound. A search asked to beat a point
    // the bound has already matched cannot succeed and does not stop trying, and the
    // instances in that state are the ones that need what is left of the run to
    // certify what they are holding: `acc-tight2`, `disctom` and `neos-913984` each
    // hold their optimum after one node, and an unguarded improvement search at the
    // root takes all three from optimal to a timeout reporting a gap of zero. This is
    // the same trap recorded against the earlier attempt at root improvement, and it
    // was walked into again by assuming the situation had changed.
    if options.improvement_frequency > 0
        && improves(root_bound, incumbent, options.gap_tolerance)
        && let Some((objective, x)) = reduced_cost_neighbourhood(
            problem,
            &lurking,
            root_bound,
            &options,
            remaining_of(options.time_limit, started),
        )
        && objective < incumbent
    {
        incumbent = objective;
        incumbent_x = Some(x);
        heuristic_solutions += 1;
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
            &lurking,
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
        let solution = Solution {
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
        let restart = result
            .wants_restart
            .then(|| narrow_for_restart(problem, &lurking, result.incumbent))
            .flatten();
        return Attempt { solution, restart };
    }

    let mut worker = Worker::new(problem, lp, options, &lurking);
    let mut restart_wanted = false;

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
        if worker.wants_restart(incumbent, &options) {
            restart_wanted = true;
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
            //
            // Its own bound is kept, because it decides later whether the claim was
            // really forfeited: the bound this node inherited holds over its whole
            // subtree, so if the incumbent has since overtaken it there was nothing in
            // there worth having and skipping it cost nothing.
            unresolved += 1;
            unresolved_bound = unresolved_bound.min(node.bound);
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
            && let Some((objective, x)) = improve(
                problem,
                current,
                incumbent,
                relaxation,
                &options,
                remaining_of(options.time_limit, started),
            )
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
        // A search that exhausted its tree but skipped a node has not proven anything,
        // and must not say it has -- unless the incumbent has since overtaken every
        // bound it skipped, which is the same test the search applies before opening a
        // node at all.
        (Status::Optimal, _)
            if unresolved > 0 && improves(unresolved_bound, incumbent, options.gap_tolerance) =>
        {
            Status::NodeLimit
        }
        (Status::Optimal, None) => Status::Infeasible,
        (other, _) => other,
    };

    let solution = Solution {
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
    };
    let restart = restart_wanted
        .then(|| narrow_for_restart(problem, &lurking, incumbent))
        .flatten();
    Attempt { solution, restart }
}

/// The model a restart should begin from: this one, with everything the incumbent has
/// proved written into its bounds.
///
/// Only bounds the lurking table has *earned* at this incumbent, so every one of them is
/// a proof rather than a supposition, and the narrowed model has the same optimum as the
/// one handed in wherever that optimum beats the incumbent. Returns `None` when nothing
/// was earned, which makes starting over pointless.
fn narrow_for_restart(problem: &Problem, lurking: &Lurking, incumbent: f64) -> Option<Problem> {
    if !incumbent.is_finite() {
        return None;
    }
    let granularity = objective_granularity(problem);
    let cutoff = match granularity {
        Some(g) => lift_to_granularity(incumbent, granularity) - g,
        None => incumbent,
    };
    let mut narrowed = problem.clone();
    let mut changed = 0usize;
    for entry in &lurking.entries[..lurking.in_force(cutoff)] {
        let j = entry.column as usize;
        let (lo, hi) = (narrowed.col_lb[j], narrowed.col_ub[j]);
        narrowed.col_lb[j] = entry.lo.max(lo);
        narrowed.col_ub[j] = entry.hi.min(hi);
        if narrowed.col_lb[j] != lo || narrowed.col_ub[j] != hi {
            changed += 1;
        }
    }
    (changed > 0).then_some(narrowed)
}

/// The share of a run the root cut loop may spend before the search takes over.
///
/// Cutting earns its place on most models and cannot be allowed to earn it on all of
/// them: a round costs an LP solve, and on a large model fifty of those is the whole
/// budget. A third leaves the bound most of what cutting was going to give it while
/// guaranteeing the search two thirds of the run to spend it in.
const ROOT_CUT_SHARE: f64 = 0.33;

/// How much of the run the relaxation gets before it is treated as stuck.
///
/// Generous: a relaxation that is merely slow should be allowed to finish, and only one
/// that is going nowhere is worth restarting from a different point.
const ROOT_LP_FIRST_SHARE: f64 = 0.40;

/// The reciprocal of the share of a round's cut allowance mod-2 separation may take.
///
/// An eighth. These are not comparable to the other families: every one is violated by
/// exactly one half, so ranking them alongside families whose violation varies admits
/// all of them or none, and a handful is all they need to be, since eighteen close
/// `n2seq36f` in two rounds. Taken at a quarter they crowd the node LPs on a model with
/// few, wide rows, and `irp`, 39 rows against 20315 columns, stops closing at all.
const MOD2_SHARE: usize = 8;

/// The share of a round's cuts reserved for rows aggregated through implied bounds.
///
/// A quarter. They are the whole answer where they apply, so the share only has to be
/// wide enough to carry a model's covering rows over a few rounds -- `neos-787933` wants
/// 77 of them and has 133 in all -- and narrow enough that a model where they are merely
/// valid does not spend its round on them.
const AGGREGATION_SHARE: usize = 4;

/// A backstop on presolve's expensive half, not its budget.
///
/// A quarter until probing was given a budget it could finish on: `ex9` wants 22 of its
/// 60 seconds in presolve and then closes in three more, so a quarter was refusing it
/// the reduction rather than capping a runaway. Probing's own guards are what stop the
/// models that should be stopped, and they leave every model the solver closes under a
/// second here, so this only has to be wide enough not to cut off a model that is still
/// paying its way.
const PRESOLVE_SHARE: f64 = 0.4;

/// How far the bounds are moved when a relaxation has stalled.
///
/// Small enough that the perturbed problem is next to the real one, large enough to
/// separate variables sitting on a shared bound. Swept over 1e-7 to 1e-4; 1e-6 finishes
/// the stalled relaxations and leaves the clean-up solve little to undo.
const PERTURBATION_REACH: f64 = 1e-6;

/// Solve a stalled relaxation by way of a perturbed one.
///
/// Returns the true model's solution, warm started from the perturbed model's basis, or
/// `None` if either half failed to finish. Nothing here is approximate: the bound handed
/// back is the true relaxation's, proved on the true bounds.
fn perturbed_root(
    problem: &Problem,
    exact: &Lp,
    deadline: Option<Instant>,
    max_iterations: usize,
) -> Option<LpSolution> {
    // SplitMix64, so a run is reproducible.
    let mut state = 0x5DEE_CE66Du64;
    let mut noise = move || {
        state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        ((z ^ (z >> 31)) >> 11) as f64 / (1u64 << 53) as f64
    };
    let mut widen = |lo: &mut f64, hi: &mut f64| {
        let reach = PERTURBATION_REACH * (1.0 + lo.abs().max(hi.abs()).min(1e6)) * (1.0 + noise());
        if lo.is_finite() {
            *lo -= reach;
        }
        if hi.is_finite() {
            *hi += reach;
        }
    };

    let mut loosened = problem.clone();
    for j in 0..loosened.n_cols() {
        let (mut lo, mut hi) = (loosened.col_lb[j], loosened.col_ub[j]);
        widen(&mut lo, &mut hi);
        loosened.col_lb[j] = lo;
        loosened.col_ub[j] = hi;
    }
    for i in 0..loosened.n_rows() {
        let (mut lo, mut hi) = (loosened.row_lb[i], loosened.row_ub[i]);
        widen(&mut lo, &mut hi);
        loosened.row_lb[i] = lo;
        loosened.row_ub[i] = hi;
    }

    let mut loose = Lp::relaxation(&loosened);
    loose.set_deadline(deadline);
    let solved = loose.solve_with_limit(max_iterations);
    if solved.status != LpStatus::Optimal {
        return None;
    }
    // The perturbed objective is a bound on a weaker problem and is not used. Only the
    // basis crosses back.
    let cleaned = exact.solve_with_rows(&solved.basis, &[], None, max_iterations);
    let mut cleaned = cleaned;
    cleaned.iterations += solved.iterations;
    (cleaned.status == LpStatus::Optimal).then_some(cleaned)
}

/// The spacing of the objective values this model can actually take.
///
/// When every column carrying an objective coefficient is an integer column, and every
/// one of those coefficients is a whole multiple of `g`, no feasible point scores
/// anything but a multiple of `g`. A relaxation bound of 52000.3 is then really a bound
/// of 52200, and a node whose relaxation lands a hair above 52200 is not a node whose
/// subtree might hold 52200.1, it is a node whose subtree holds nothing better than
/// 52400.
///
/// `n2seq36f` is the case that makes the point: every objective coefficient is a
/// multiple of 200, its relaxation is worth exactly 52000 and its optimum is 52200, one
/// step apart. Four thousand cuts do not move that bound and do not need to.
///
/// Returns `None` where the reasoning does not apply, which is a fractional coefficient
/// or one on a column free to take fractional values.
fn objective_granularity(problem: &Problem) -> Option<f64> {
    // Below this a coefficient is read as zero rather than as a very fine spacing, and
    // above it a whole multiple has to look like one.
    const WHOLE: f64 = 1e-9;
    let mut granularity: u64 = 0;
    for j in 0..problem.n_cols() {
        let c = problem.obj[j];
        if c.abs() <= WHOLE {
            continue;
        }
        // A continuous column moves the objective continuously, whatever its
        // coefficient, unless the model has already pinned it.
        if !problem.is_integer(j) && problem.col_lb[j] < problem.col_ub[j] {
            return None;
        }
        let rounded = c.round();
        if (c - rounded).abs() > WHOLE * c.abs().max(1.0) {
            return None;
        }
        let magnitude = rounded.abs();
        // Beyond this the coefficient no longer converts to an integer exactly, and a
        // greatest common divisor of approximations is not a granularity.
        if magnitude > 2.0f64.powi(53) {
            return None;
        }
        granularity = gcd(granularity, magnitude as u64);
        // No early exit once the divisor reaches one, tempting as it is. Whether this
        // reasoning applies at all is decided by *every* column, and returning as soon
        // as the spacing stops improving skips the ones not yet looked at. A fuzz
        // instance whose first two integer coefficients were 8 and 3 returned a spacing
        // of one before reaching the continuous column with a coefficient of 6.5, and
        // the search then proved an optimum of 34.375 where the answer is 34.25.
    }
    (granularity > 0).then_some(granularity as f64)
}

fn gcd(a: u64, b: u64) -> u64 {
    if b == 0 { a } else { gcd(b, a % b) }
}

/// Does every feasible point of this model score the same?
///
/// A model whose objective is empty, or whose only objective coefficients sit on columns
/// the model has already pinned, is asking for a feasible point rather than a good one.
/// Two of the instances this solver misses are exactly that, `neos-3226448-wkra` and
/// `supportcase4`, and on both of them it spends the whole run on a relaxation whose
/// value it could have written down: the bound is that constant, no search can improve
/// on it, and the first feasible point is optimal.
fn objective_is_constant(problem: &Problem) -> bool {
    (0..problem.n_cols()).all(|j| {
        problem.obj[j].abs() <= 1e-12 || problem.col_lb[j] >= problem.col_ub[j] - 1e-12
    })
}

/// What every feasible point scores, where they all score the same.
fn constant_objective(problem: &Problem) -> f64 {
    (0..problem.n_cols())
        .filter(|&j| problem.obj[j].abs() > 1e-12)
        .map(|j| problem.obj[j] * problem.col_lb[j])
        .sum()
}

/// Lift a bound to the next objective value the model can actually take.
///
/// Sound because every feasible point of every subtree scores a multiple of `g`, so the
/// least multiple at or above a valid bound is also a valid bound. The tolerance is
/// what keeps it sound in floating point: a bound sitting a whisker above a multiple
/// because of rounding must round *to* that multiple and not past it, which would claim
/// a whole step more than was proven.
fn lift_to_granularity(bound: f64, granularity: Option<f64>) -> f64 {
    let Some(g) = granularity else { return bound };
    if !bound.is_finite() {
        return bound;
    }
    let steps = bound / g;
    let tolerance = 1e-9 * steps.abs().max(1.0);
    (steps - tolerance).ceil() * g
}

/// A bound the root's reduced costs will imply, once the incumbent is good enough.
///
/// Everything reduced cost fixing needs is fixed at the root except one number. The
/// reduced costs come from the root basis and do not change; the root's value does not
/// change; only the room `u - root` shrinks as the search finds better points. So the
/// bound a column will eventually take, and the incumbent at which it takes it, can both
/// be worked out once and read off later.
///
/// Which is what makes this affordable in a parallel search. Re-deriving the fixing
/// whenever the incumbent improved would mean narrowing a model that every worker holds
/// immutably; a table computed once and read against whatever incumbent a worker
/// currently sees needs no coordination at all, and a worker that is behind on the
/// incumbent simply applies fewer of them.
struct Lurking {
    /// Sorted by `threshold` descending, so the entries in force are always a prefix.
    /// Later entries for the same column are the tighter ones, since a tighter bound
    /// needs a better incumbent, so applying a prefix in order leaves each column on
    /// the best bound it has earned.
    entries: Vec<LurkingBound>,
}

struct LurkingBound {
    /// The incumbent at or below which this bound holds.
    threshold: f64,
    column: u32,
    lo: f64,
    hi: f64,
}

/// Lurking bounds recorded per column for a general integer.
///
/// A binary needs one: the only bound worth recording is the one that fixes it. A column
/// with a wide range could have one per value it might be narrowed to, which is a table
/// the size of the model's total range, so it is sampled instead.
const LURKING_STEPS: usize = 8;

impl Lurking {
    /// Work out, for every nonbasic column, what the root's reduced cost will prove
    /// about it and when.
    ///
    /// A column parked at its lower bound with reduced cost `d > 0` can be capped at
    /// `lo + k` once travelling `k + 1` would pass the incumbent, which is
    /// `u < root + d (k + 1)`. Recording that for each reachable `k` costs nothing now
    /// and saves re-deriving it from a factorization later.
    fn build(problem: &Problem, costs: &[Option<(f64, bool)>], root: f64, tolerance: f64) -> Self {
        let mut entries: Vec<LurkingBound> = Vec::new();
        for (j, entry) in costs.iter().enumerate() {
            let Some((d, at_upper)) = *entry else {
                continue;
            };
            let (lo, hi) = (problem.col_lb[j], problem.col_ub[j]);
            if lo >= hi || !lo.is_finite() || !hi.is_finite() || !problem.is_integer(j) {
                continue;
            }
            let reach = d.abs();
            if reach <= tolerance {
                continue;
            }
            // The range is sampled rather than enumerated when it is wide, which for a
            // binary is one step and the whole story.
            let span = (hi - lo).round().max(1.0);
            let step = (span / LURKING_STEPS as f64).ceil().max(1.0);
            let mut travel = 0.0f64;
            while travel < span {
                // Slightly inside what the arithmetic proves, so that a threshold met
                // by a hair does not fix a column. Every margin here is on the side
                // that fixes later.
                let threshold = (travel + 1.0).mul_add(reach, root) - tolerance * reach;
                if at_upper {
                    entries.push(LurkingBound {
                        threshold,
                        column: j as u32,
                        lo: hi - travel,
                        hi,
                    });
                } else {
                    entries.push(LurkingBound {
                        threshold,
                        column: j as u32,
                        lo,
                        hi: lo + travel,
                    });
                }
                travel += step;
            }
        }
        entries.sort_by(|a, b| b.threshold.total_cmp(&a.threshold));
        Self { entries }
    }

    /// How many entries are in force at this cutoff. They are the leading ones, because
    /// the table is ordered by the cutoff each entry needs.
    fn in_force(&self, cutoff: f64) -> usize {
        self.entries
            .partition_point(|entry| entry.threshold >= cutoff)
    }

    fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}

/// Narrow the columns the root's bound and the incumbent together already decide.
///
/// At an optimal basis every nonbasic column sits on a bound with a reduced cost `d`
/// whose sign says the objective cannot fall by leaving it. Moving distance `t` off that
/// bound therefore raises the objective by at least `|d| t`, and since the root's value
/// is a bound on every solution in the tree, anything better than the incumbent `u`
/// satisfies `root + |d| t < u`. That caps `t` at `(u - root) / |d|`, and for an integer
/// column the cap rounds inwards, which on a binary usually means fixing it outright.
///
/// Returns the narrowed model, or `None` when nothing moved and the caller should keep
/// the one it has rather than pay for a copy.
fn fix_by_reduced_cost(
    problem: &Problem,
    lp: &Lp,
    root: &LpSolution,
    cutoff: f64,
    options: &Options,
) -> Option<Problem> {
    let costs = lp.reduced_costs(&root.basis)?;
    // The room a better solution has to move in. Any margin on this belongs on the
    // generous side: a *smaller* room caps the travel harder and so fixes more, which
    // is the direction that can cut off a solution nobody has seen yet. The tolerances
    // below are all widenings for the same reason.
    let room = cutoff - root.objective;
    if room < 0.0 || !room.is_finite() {
        return None;
    }

    let mut narrowed = problem.clone();
    let mut changed = 0usize;
    for (j, entry) in costs.iter().enumerate() {
        let Some((d, at_upper)) = *entry else {
            continue;
        };
        let (lo, hi) = (narrowed.col_lb[j], narrowed.col_ub[j]);
        if lo >= hi || !lo.is_finite() || !hi.is_finite() {
            continue;
        }
        let reach = d.abs();
        if reach <= 1e-9 {
            continue;
        }
        // Widened by a relative and an absolute slack, because `room` and `reach` are
        // both floating point results of a long solve and the claim being made is that
        // no better solution exists past this point.
        let travel = (room / reach).mul_add(1.0 + 1e-9, 1e-9);
        if travel >= hi - lo {
            continue;
        }
        // Integer columns round the cap inwards; a fractional cap of 0.4 on a binary
        // parked at zero means it can only be zero.
        let travel = if narrowed.is_integer(j) {
            (travel + options.integrality_tolerance).floor()
        } else {
            travel
        };
        if at_upper {
            narrowed.col_lb[j] = (hi - travel).min(hi).max(lo);
        } else {
            narrowed.col_ub[j] = (lo + travel).max(lo).min(hi);
        }
        if narrowed.col_lb[j] != lo || narrowed.col_ub[j] != hi {
            changed += 1;
        }
    }
    (changed > 0).then_some(narrowed)
}

/// How much of the run the LP-free feasibility search may have.
///
/// It runs before the relaxation, so whatever it takes is taken from everything else.
/// A share of the limit is the wrong yardstick and this is deliberately a small one.
/// The models this cracks it cracks in about a second: `acc-tight2`, `disctom` and
/// `neos-913984` each cost under a second and a half of jumping to close outright. The
/// models it cannot crack will take whatever they are given, and at a tenth of a two
/// minute limit that came to twelve seconds spent on instances that solve in under one,
/// taking `f2gap401600` from 0.27 seconds to 11.5. The stall cutoff inside the search
/// is what actually bounds it; this is the backstop.
const JUMP_SHARE: f64 = 0.05;

/// Randomised starts the feasibility search may take before giving up.
///
/// It stops when it has settled, and where it settles is decided by where it began, so
/// beginning somewhere else is the cheapest thing to try next. The count is generous
/// because the deadline is what actually stops this: `neos-3226448-wkra` needs
/// seventeen starts and nine seconds, and a model that is never going to yield spends
/// exactly the same budget either way.
const JUMP_RESTARTS: usize = 200;

/// How many times a solve may start over on the model its own search narrowed.
///
/// Each one repeats presolve, the root relaxation and the cut loop, so it is only worth
/// taking when what has been fixed since the last start is a large part of the model.
/// The count is a backstop; what actually ends the sequence is running out of columns
/// to fix or time to spend.
const MAX_RESTARTS: usize = 3;

/// How much of the model the reduced cost neighbourhood may decide, as a percentage.
///
/// Supposing a better answer is a supposition, and the further it is pushed the less
/// likely the neighbourhood is to hold anything at all. With half the model decided,
/// `air05` and `neos-820879` come back with no point at all; with three tenths, `air05`
/// returns 26827 against an optimum of 26374 and `neos-820879` 26348 against 25468.
const NEIGHBOURHOOD_CEILING: usize = 30;

/// The share of what is left that the reduced cost neighbourhood may spend.
///
/// Three tenths, and the window is narrow at both ends. At a twentieth
/// `neos-3045796-mogo` comes back with an incumbent of 1380; at a fifth it reaches its
/// optimum of -175 in two runs of four; at three tenths in three of three; at two fifths
/// in none of three, the search having been left too little to certify what it holds.
const NEIGHBOURHOOD_SHARE: f64 = 0.3;

/// The share of the run a model with nothing to optimise gives its feasibility search.
///
/// Most of it, because there is nothing else to spend it on. What is left over goes to
/// the ordinary path, which can still prove the model infeasible where this cannot.
const CONSTANT_OBJECTIVE_SHARE: f64 = 0.75;

/// What is left of a time limit that started at `started`, or `None` for no limit.
fn remaining_of(limit: Option<Duration>, started: Instant) -> Option<Duration> {
    limit.map(|limit| limit.saturating_sub(started.elapsed()))
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
/// Search the neighbourhood the root's reduced costs would fix if the answer were
/// better than anything yet found.
///
/// The lurking table already says, for every column, the bound it takes and the
/// incumbent at which it takes it. Read in the other direction it is a construction:
/// *suppose* the answer beats the weakest threshold in the table, and the entries at or
/// above it are all true. Applying them leaves a much smaller model, and its optimum is
/// either a better solution to the original or a proof that the supposition was wrong.
/// Nothing about the sub-model is unsound -- its bounds are only tighter -- so any point
/// it returns is a genuine solution to the original.
///
/// This is what the ordinary improvement search cannot do here. That one fixes columns
/// where the incumbent and the relaxation agree, and where the incumbent came from a
/// feasibility search rather than from the relaxation, they agree about nothing useful:
/// on `neos-3045796-mogo` it walks an incumbent of 1180 down to 300 against an optimum
/// of -175. This construction does not read the incumbent at all.
fn reduced_cost_neighbourhood(
    problem: &Problem,
    lurking: &Lurking,
    root_bound: f64,
    options: &Options,
    remaining: Option<Duration>,
) -> Option<(f64, Vec<f64>)> {
    if lurking.is_empty() || remaining.is_some_and(|left| left.is_zero()) {
        return None;
    }
    let integers = problem.integer_columns().count();
    if integers == 0 {
        return None;
    }
    let mut narrowed = problem.clone();
    let mut fixed = 0usize;
    // Ordered by the incumbent each entry needs, so this walks from the entries a
    // barely-better answer would justify towards the ones only a much better answer
    // would, and stops as soon as enough of the model is decided.
    for entry in &lurking.entries {
        if entry.threshold <= root_bound {
            break;
        }
        let j = entry.column as usize;
        let (lo, hi) = (narrowed.col_lb[j], narrowed.col_ub[j]);
        if lo >= hi {
            continue;
        }
        narrowed.col_lb[j] = entry.lo.max(lo);
        narrowed.col_ub[j] = entry.hi.min(hi);
        if narrowed.col_lb[j] >= narrowed.col_ub[j] {
            fixed += 1;
            if fixed * 100 >= integers * NEIGHBOURHOOD_CEILING {
                break;
            }
        }
    }
    // Too little of the model decided and this is the original problem again, which the
    // caller is already searching.
    // A tenth, not the third HiGHS uses for the same construction. Its neighbourhood is
    // grown by propagation, which turns each fixed column into several; nothing here
    // propagates, so the same share of the model needs far more entries to reach and
    // the cases that pay do not reach it. `neos-3045796-mogo` fixes 17% and is worth
    // searching: its incumbent goes from 1180 to -155 against an optimum of -175.
    if fixed * 10 < integers {
        return None;
    }

    let found = solve(
        &narrowed,
        Options {
            improvement_frequency: 0,
            max_nodes: options.improvement_nodes,
            threads: 1,
            // With its own feasibility search, unlike the improvement search next door.
            // That one starts from an incumbent and is looking for a better point near
            // it, so re-establishing feasibility is wasted; this one starts from no
            // point at all, and finding one is the whole job. Inheriting the rule from
            // the neighbour left `neos-3045796-mogo` at 150 where it reaches -155.
            jump_moves: options.jump_moves,
            // A share of what is left, not all of it. This runs before the tree opens,
            // so whatever it spends the search never gets.
            time_limit: remaining.map(|left| left.mul_f64(NEIGHBOURHOOD_SHARE)),
            ..*options
        },
    );
    let x = found.x;
    if x.len() != problem.n_cols() {
        return None;
    }
    // Feasible for the original by construction, but checked rather than assumed: the
    // sub-search reports a point from whatever model it actually ran on.
    let tolerance = options.heuristic_limits.feasibility_tolerance;
    heuristic::is_feasible(problem, &x, tolerance).then(|| (objective_at(problem, &x), x))
}

fn improve(
    problem: &Problem,
    incumbent: &[f64],
    incumbent_objective: f64,
    relaxation: &[f64],
    options: &Options,
    remaining: Option<Duration>,
) -> Option<(f64, Vec<f64>)> {
    // Out of time is out of time, and a sub-search with nothing left to spend would
    // still pay for the model copy below.
    if remaining.is_some_and(|left| left.is_zero()) {
        return None;
    }
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
            // The sub-search starts from an incumbent and is looking for a better point
            // near it, which is not what a feasibility search is for. Inherited, it
            // runs again on every neighbourhood: on `eil33-2` that came to 51 seconds
            // of a 96 second solve, spent re-establishing feasibility that the
            // neighbourhood already had.
            jump_moves: 0,
            // What is left of the caller's budget, not a fresh copy of it. A time
            // limit is a duration and the sub-search starts its own clock, so passing
            // the original handed it the whole limit again: an improvement beginning
            // at fifty-five seconds of a sixty second run could spend another sixty.
            // Measured, n13-3 returned at 71.6s against a 60s limit.
            time_limit: remaining,
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

    #[test]
    fn a_fractional_coefficient_anywhere_rules_out_an_objective_spacing() {
        // The column that rules the reasoning out can be anywhere, including after the
        // spacing has already fallen to one. A version that stopped looking at that
        // point read this model as integral and proved an optimum of 34.375 where the
        // answer is 34.25.
        let mut model = crate::model::Builder::new(Sense::Minimize);
        let a = model.integer("a", 0.0, 4.0);
        let b = model.integer("b", 0.0, 4.0);
        let c = model.continuous("c", 0.0, 4.0);
        // 8 and 3 are coprime, so the spacing is already 1 before `c` is reached.
        model.objective(&[(a, 8.0), (b, 3.0), (c, 6.5)]);
        model.row(&[(a, 1.0), (b, 1.0), (c, 1.0)], RowSense::Ge, 1.0);
        assert_eq!(objective_granularity(&model.build()), None);
    }

    #[test]
    fn a_whole_objective_spacing_lifts_a_bound_to_the_next_value() {
        let mut model = crate::model::Builder::new(Sense::Minimize);
        let a = model.integer("a", 0.0, 4.0);
        let b = model.integer("b", 0.0, 4.0);
        model.objective(&[(a, 200.0), (b, 600.0)]);
        model.row(&[(a, 1.0), (b, 1.0)], RowSense::Ge, 1.0);
        let problem = model.build();
        let granularity = objective_granularity(&problem);
        assert_eq!(granularity, Some(200.0));
        // Between two attainable values, so the bound is really the upper one.
        assert_eq!(lift_to_granularity(52000.3, granularity), 52200.0);
        // Already attainable, and must not be pushed a whole step past what was proven,
        // including when rounding leaves it a whisker above.
        assert_eq!(lift_to_granularity(52200.0, granularity), 52200.0);
        assert_eq!(lift_to_granularity(52200.000000001, granularity), 52200.0);
        assert_eq!(lift_to_granularity(-400.5, granularity), -400.0);
    }

    #[test]
    fn a_point_without_a_bound_reports_no_gap_rather_than_none_left() {
        // The state a time limit at an unfinished root leaves behind: a heuristic point
        // in hand and nothing proven about it. Reported as a gap of zero, that reads as
        // a proof of optimality, which is the one thing it is not.
        let mut solution =
            Solution::without_solution(Status::TimeLimit, 1, 0, None);
        solution.objective = Some(0.0);
        assert!(!solution.bound.is_finite());
        assert_eq!(solution.gap(), f64::INFINITY);

        // A bound that was proven still reports the gap it proves.
        solution.bound = -1.0;
        assert!(solution.gap().is_finite());
    }

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
        let (improved, x) = improve(&p, &incumbent, value, &relaxation, &options, None)
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
        assert!(improve(&p, &incumbent, value, &disagrees, &options, None).is_none());

        // Agreeing on everything leaves the incumbent, which cannot be improved on.
        assert!(improve(&p, &incumbent, value, &incumbent, &options, None).is_none());
    }

    /// A time limit is a duration and the sub-search starts its own clock, so it has
    /// to be handed what is left of the caller's budget rather than a fresh copy of
    /// it. Handing over the original let an improvement beginning near the end of a
    /// run spend the whole limit again: n13-3 returned at 71.6s against a 60s limit,
    /// and 60.3s once this was passed through.
    ///
    /// What this pins is that `remaining` is honoured at all. The early return and the
    /// limit handed to the sub-search each produce the right answer on their own, so
    /// neither fails this alone; removing both does. The early return is there to skip
    /// a whole-model clone once the budget is gone, which is work rather than answer
    /// and so is not visible from here.
    #[test]
    fn the_improvement_search_gets_what_is_left_of_the_budget_and_no_more() {
        let p = knapsackish();
        // A caller with a generous limit overall, which is the value the sub-search
        // must *not* be handed once the run is nearly over.
        let options = Options {
            time_limit: Some(Duration::from_secs(30)),
            ..Options::default()
        };
        let incumbent = vec![1.0, 1.0, 0.0];
        let value = objective_at(&p, &incumbent);
        let relaxation = vec![1.0, 0.5, 1.0];

        // Nothing left means nothing spent, and no answer either.
        assert!(
            improve(
                &p,
                &incumbent,
                value,
                &relaxation,
                &options,
                Some(Duration::ZERO)
            )
            .is_none(),
            "an exhausted budget still ran a sub-search"
        );

        // With budget it does the work it is there to do.
        assert!(
            improve(
                &p,
                &incumbent,
                value,
                &relaxation,
                &options,
                Some(Duration::from_secs(30))
            )
            .is_some(),
            "a funded improvement search found nothing"
        );
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
