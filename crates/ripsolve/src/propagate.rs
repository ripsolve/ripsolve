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
///
/// Carries where the contradiction was found, which is where [`explain`] starts from.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Infeasible {
    source: Source,
}

/// Scratch reused across nodes, so propagating costs no allocation.
///
/// A worker keeps one of these for the life of the search. The queue and the flag vector
/// are cleared per call rather than rebuilt.
pub struct Scratch {
    /// Columns whose bounds have moved and whose consequences are not drawn yet.
    queue: Vec<u32>,
    /// Whether a column is already waiting, so it is queued once however often it moves.
    queued: Vec<bool>,
    /// Why each column's bounds are where they are, for [`explain`]. Only meaningful for
    /// a column this call actually moved, or seeded.
    reason: Vec<Option<Reason>>,
    /// Columns whose `reason` this call wrote, so the next one can clear just those.
    touched: Vec<u32>,
    /// Working stack for [`explain`], kept here so explaining allocates nothing.
    stack: Vec<u32>,
    seen: Vec<bool>,
    seen_list: Vec<u32>,
}

impl Scratch {
    pub fn new(n_cols: usize) -> Self {
        Scratch {
            queue: Vec::new(),
            queued: vec![false; n_cols],
            reason: vec![None; n_cols],
            touched: Vec::new(),
            stack: Vec::new(),
            seen: vec![false; n_cols],
            seen_list: Vec::new(),
        }
    }

    fn clear_reasons(&mut self) {
        for &j in &self.touched {
            self.reason[j as usize] = None;
        }
        self.touched.clear();
    }

    fn set_reason(&mut self, j: usize, why: Reason) {
        if self.reason[j].is_none() {
            self.touched.push(j as u32);
        }
        self.reason[j] = Some(why);
    }
}

/// Why a column's bounds are where they are.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Reason {
    /// A bound this node's branching chose. A clause may rest on it: the branch is a
    /// statement about where in the tree we are, not about which points are feasible,
    /// so "these branchings cannot all hold" follows from the model alone.
    Decision,
    /// A bound the incumbent earned, from reduced cost fixing. **A clause may not rest
    /// on one.** It is only valid while that incumbent stands, so a conflict explained
    /// through it is not a property of the model and must not be learned.
    Assumption,
    /// Forced because the named column was fixed and the two conflict.
    Conflict(u32),
    /// Forced by the named row's activity under the bounds then in force.
    Row(u32),
    /// Forced by the named learned clause, every other literal of which had failed.
    Clause(u32),
    /// Moved, but not to a single value, so what justifies its present bounds is not one
    /// implication. Nothing is learned through one of these; on a binary model they do
    /// not arise, because any bound change on a binary fixes it.
    Opaque,
}

/// What a seeded bound change is: see [`Reason`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Seed {
    Decision,
    Assumption,
}

/// Where propagation ran out of room, which is where an explanation starts.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Source {
    /// A row whose activity cannot meet its bounds.
    Row(u32),
    /// Two columns the conflict graph says cannot both hold as they are fixed.
    Pair(u32, u32),
    /// A learned clause with every literal failed.
    Clause(u32),
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
    clauses: &Clauses,
    scratch: &mut Scratch,
    lb: &mut [f64],
    ub: &mut [f64],
    seeds: &[(u32, Seed)],
    budget: usize,
    tolerance: f64,
) -> Result<usize, Infeasible> {
    scratch.queue.clear();
    scratch.clear_reasons();
    for &(j, seed) in seeds {
        scratch.set_reason(
            j as usize,
            match seed {
                Seed::Decision => Reason::Decision,
                Seed::Assumption => Reason::Assumption,
            },
        );
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
            return Err(Infeasible { source: Source::Pair(j as u32, j as u32) });
        }

        // Conflicts, when this column is fixed and binary.
        if is_binary(problem, j) && (ub[j] - lb[j]).abs() <= tolerance {
            let value = lb[j] > 0.5;
            let node = literal(j, value);
            let mut contradiction = None;
            let mut forced_columns: Vec<(usize, f64)> = Vec::new();
            conflicts.for_each_adjacent(node, |excluded| {
                let k = (excluded / 2) as usize;
                // An even literal is "column k is one"; being excluded forces it to zero.
                let forced = if excluded.is_multiple_of(2) { 0.0 } else { 1.0 };
                if forced > ub[k] + tolerance || forced < lb[k] - tolerance {
                    contradiction = Some(k);
                    return false;
                }
                if (ub[k] - lb[k]).abs() > tolerance {
                    lb[k] = forced;
                    ub[k] = forced;
                    fixed += 1;
                    forced_columns.push((k, forced));
                }
                true
            });
            if let Some(k) = contradiction {
                return Err(Infeasible { source: Source::Pair(j as u32, k as u32) });
            }
            for (k, _) in forced_columns {
                scratch.set_reason(k, Reason::Conflict(j as u32));
                if !scratch.queued[k] {
                    scratch.queued[k] = true;
                    scratch.queue.push(k as u32);
                }
            }
            work += conflicts.degree(node);
        }

        // The clauses this column appears in. A clause with every literal but one
        // failed forces that one; with every literal failed it is the conflict itself.
        if !clauses.is_empty() {
            for &id in &clauses.mentions[j] {
                let clause = &clauses.clauses[id as usize];
                work += clause.len();
                let mut unassigned = None;
                let mut satisfied = false;
                for &(k, want) in clause {
                    let k = k as usize;
                    let fixed_here = (ub[k] - lb[k]).abs() <= tolerance;
                    if fixed_here {
                        if (lb[k] > 0.5) == want {
                            satisfied = true;
                            break;
                        }
                    } else if unassigned.is_some() {
                        unassigned = Some(usize::MAX); // more than one: nothing to force
                        break;
                    } else {
                        unassigned = Some(k);
                    }
                }
                if satisfied {
                    continue;
                }
                match unassigned {
                    None => return Err(Infeasible { source: Source::Clause(id) }),
                    Some(usize::MAX) => {}
                    Some(k) => {
                        let want = clause
                            .iter()
                            .find(|&&(c, _)| c as usize == k)
                            .map(|&(_, w)| w)
                            .unwrap_or(false);
                        let forced = if want { 1.0 } else { 0.0 };
                        if forced > ub[k] + tolerance || forced < lb[k] - tolerance {
                            return Err(Infeasible { source: Source::Clause(id) });
                        }
                        lb[k] = forced;
                        ub[k] = forced;
                        fixed += 1;
                        scratch.set_reason(k, Reason::Clause(id));
                        if !scratch.queued[k] {
                            scratch.queued[k] = true;
                            scratch.queue.push(k as u32);
                        }
                    }
                }
            }
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
                return Err(Infeasible { source: Source::Row(i as u32) });
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
                    return Err(Infeasible { source: Source::Row(i as u32) });
                }
                // Only a real move is worth queueing; a bound that repeats itself would
                // otherwise cycle the worklist forever on a degenerate row.
                let moved = new_lb > lb[k] + tolerance || new_ub < ub[k] - tolerance;
                if moved {
                    let was_free = (ub[k] - lb[k]).abs() > tolerance;
                    lb[k] = lb[k].max(new_lb);
                    ub[k] = ub[k].min(new_ub);
                    let now_fixed = (ub[k] - lb[k]).abs() <= tolerance;
                    if was_free && now_fixed {
                        fixed += 1;
                    }
                    // Only a fixing has one implication behind it; a partial tightening
                    // does not, and nothing is learned through it. See `Reason::Opaque`.
                    scratch.set_reason(k, if now_fixed { Reason::Row(i as u32) } else { Reason::Opaque });
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

/// A learned clause: the conflict holds only if *every* one of these fails, so at least
/// one must hold in any feasible point. `(column, value)` reads "column takes value".
pub type Clause = Vec<(u32, bool)>;

/// Clauses learned from conflicts, and which columns each mentions.
///
/// Kept per worker rather than shared. A shared pool would let one thread's conflict
/// prune another's subtree, which is worth having, but it puts a lock on the propagation
/// path that every node runs; per-worker costs nothing and still pays, because a worker
/// takes nodes from all over the tree. Sharing is the obvious next step and is not this.
pub struct Clauses {
    clauses: Vec<Clause>,
    /// For each column, the clauses that mention it.
    mentions: Vec<Vec<u32>>,
    /// Most clauses to keep, and the longest one worth keeping. A long clause forbids
    /// little and costs a scan at every node that touches any of its columns.
    capacity: usize,
    max_length: usize,
}

impl Clauses {
    pub fn new(n_cols: usize, capacity: usize, max_length: usize) -> Self {
        Clauses {
            clauses: Vec::new(),
            mentions: vec![Vec::new(); n_cols],
            capacity,
            max_length,
        }
    }

    /// Keep `clause`, unless it is too long or the pool is full. Returns whether it was.
    pub fn add(&mut self, clause: Clause) -> bool {
        if clause.is_empty() || clause.len() > self.max_length || self.clauses.len() >= self.capacity
        {
            return false;
        }
        let id = self.clauses.len() as u32;
        for &(j, _) in &clause {
            self.mentions[j as usize].push(id);
        }
        self.clauses.push(clause);
        true
    }

    pub fn len(&self) -> usize {
        self.clauses.len()
    }

    pub fn is_empty(&self) -> bool {
        self.clauses.is_empty()
    }
}

/// A pool holding nothing, for callers that do not learn.
impl Default for Clauses {
    fn default() -> Self {
        Clauses::new(0, 0, 0)
    }
}

/// Why the node was empty, expressed over the branchings that made it so.
///
/// Walks the implication graph back from where propagation broke, through the rows and
/// conflicts that forced each bound, until nothing is left but the branch's own choices.
/// Those choices are then a *sufficient* reason for the contradiction, and since every
/// step of the walk is an implication of the model itself, no feasible point of the model
/// satisfies all of them. Negating them is therefore a globally valid clause.
///
/// Returns `None` when nothing may be learned, which is the common case and must stay
/// cheap. That happens when the reason reaches a bound the incumbent earned rather than
/// one the branch chose ([`Reason::Assumption`]), because such a bound is not a property
/// of the model and a clause resting on it would be wrong the moment the incumbent moved;
/// when it reaches a tightening that is not a fixing ([`Reason::Opaque`]); or when a
/// column's bounds differ from the root's with nothing recorded to say why.
pub fn explain(
    problem: &Problem,
    rows_of_a: &SparseMatrix,
    clauses: &Clauses,
    scratch: &mut Scratch,
    conflict: Infeasible,
    lb: &[f64],
    ub: &[f64],
    budget: usize,
    tolerance: f64,
) -> Option<Clause> {
    // A column whose bounds are the root's took no part in anything.
    let moved = |k: usize| {
        (lb[k] - problem.col_lb[k]).abs() > tolerance || (ub[k] - problem.col_ub[k]).abs() > tolerance
    };

    scratch.stack.clear();
    for &j in &scratch.seen_list {
        scratch.seen[j as usize] = false;
    }
    scratch.seen_list.clear();

    let mut push = |scratch: &mut Scratch, k: usize| {
        if !scratch.seen[k] {
            scratch.seen[k] = true;
            scratch.seen_list.push(k as u32);
            scratch.stack.push(k as u32);
        }
    };

    match conflict.source {
        Source::Pair(a, b) => {
            push(scratch, a as usize);
            push(scratch, b as usize);
        }
        Source::Row(i) => {
            let (columns, _) = rows_of_a.column(i as usize);
            for &k in columns {
                if moved(k) {
                    push(scratch, k);
                }
            }
        }
        Source::Clause(id) => {
            for &(k, _) in &clauses.clauses[id as usize] {
                if moved(k as usize) {
                    push(scratch, k as usize);
                }
            }
        }
    }

    let mut clause: Clause = Vec::new();
    let mut work = 0usize;
    while let Some(k) = scratch.stack.pop() {
        work += 1;
        if work > budget {
            return None;
        }
        let k = k as usize;
        match scratch.reason[k] {
            Some(Reason::Decision) => {
                // The branch put this column where it is. Any feasible point that agrees
                // with the rest of the branch must differ here, so the clause asks for
                // the other value.
                if !is_binary(problem, k) {
                    return None;
                }
                let held = lb[k] > 0.5;
                clause.push((k as u32, !held));
            }
            // Not a property of the model; see this function's note.
            Some(Reason::Assumption) | Some(Reason::Opaque) | None => return None,
            Some(Reason::Conflict(from)) => push(scratch, from as usize),
            Some(Reason::Row(i)) => {
                let (columns, _) = rows_of_a.column(i as usize);
                for &other in columns {
                    if other != k && moved(other) {
                        push(scratch, other);
                    }
                }
            }
            Some(Reason::Clause(id)) => {
                for &(other, _) in &clauses.clauses[id as usize] {
                    if other as usize != k && moved(other as usize) {
                        push(scratch, other as usize);
                    }
                }
            }
        }
    }

    // A conflict explained by nothing is a model that is infeasible outright, which is
    // not this function's business to report and not a clause.
    if clause.is_empty() {
        return None;
    }
    clause.sort_unstable();
    clause.dedup();
    Some(clause)
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
            seeds.push((j as u32, Seed::Decision));
        }
        propagate(
            problem,
            &csr,
            &conflicts,
            &Clauses::default(),
            &mut scratch,
            &mut lb,
            &mut ub,
            &seeds,
            usize::MAX,
            1e-9,
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
        assert!(run(&problem, &[(x, 0.0)]).is_err());
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
            &problem,
            &csr,
            &conflicts,
            &Clauses::default(),
            &mut scratch,
            &mut lb,
            &mut ub,
            &[(x as u32, Seed::Decision)],
            0,
            1e-9,
        )
        .unwrap();
        assert_eq!(fixed, 0);
        assert_eq!((lb[y], ub[y]), (0.0, 1.0), "y is untouched, not wrongly fixed");
    }

    /// Propagate from a branch and, if it conflicts, explain it.
    fn learn(problem: &Problem, fixings: &[(usize, f64)]) -> Option<Clause> {
        let csr = problem.matrix.to_csr();
        let conflicts = Conflicts::of(problem);
        let mut scratch = Scratch::new(problem.n_cols());
        let mut lb = problem.col_lb.clone();
        let mut ub = problem.col_ub.clone();
        let mut seeds = Vec::new();
        for &(j, v) in fixings {
            lb[j] = v;
            ub[j] = v;
            seeds.push((j as u32, Seed::Decision));
        }
        let conflict = propagate(
            problem,
            &csr,
            &conflicts,
            &Clauses::default(),
            &mut scratch,
            &mut lb,
            &mut ub,
            &seeds,
            usize::MAX,
            1e-9,
        )
        .err()?;
        explain(
            problem,
            &csr,
            &Clauses::default(),
            &mut scratch,
            conflict,
            &lb,
            &ub,
            usize::MAX,
            1e-9,
        )
    }

    #[test]
    fn a_clause_names_the_branchings_that_caused_the_conflict() {
        // x + y >= 1 with both branched to zero is a conflict, and the clause has to say
        // so: at least one of them must be one.
        let mut b = Builder::new(Sense::Minimize);
        let x = b.binary("x");
        let y = b.binary("y");
        let z = b.binary("z");
        b.objective(&[(x, 1.0), (y, 1.0), (z, 1.0)]);
        b.row(&[(x, 1.0), (y, 1.0)], RowSense::Ge, 1.0);
        let problem = b.build();
        let clause = learn(&problem, &[(x, 0.0), (y, 0.0)]).expect("this conflicts");
        assert_eq!(clause, vec![(x as u32, true), (y as u32, true)]);
        // And z, which had nothing to do with it, is not in the clause.
        assert!(clause.iter().all(|&(j, _)| j != z as u32));
    }

    #[test]
    fn nothing_is_learned_through_a_bound_the_incumbent_earned() {
        // The same conflict, but y is fixed by reduced cost rather than by branching.
        // That bound holds only while the incumbent does, so a clause rested on it would
        // be wrong the moment a better incumbent arrived. It must refuse to learn.
        let mut b = Builder::new(Sense::Minimize);
        let x = b.binary("x");
        let y = b.binary("y");
        b.objective(&[(x, 1.0), (y, 1.0)]);
        b.row(&[(x, 1.0), (y, 1.0)], RowSense::Ge, 1.0);
        let problem = b.build();
        let csr = problem.matrix.to_csr();
        let conflicts = Conflicts::of(&problem);
        let mut scratch = Scratch::new(problem.n_cols());
        let (mut lb, mut ub) = (problem.col_lb.clone(), problem.col_ub.clone());
        lb[x] = 0.0;
        ub[x] = 0.0;
        lb[y] = 0.0;
        ub[y] = 0.0;
        let conflict = propagate(
            &problem,
            &csr,
            &conflicts,
            &Clauses::default(),
            &mut scratch,
            &mut lb,
            &mut ub,
            &[(x as u32, Seed::Decision), (y as u32, Seed::Assumption)],
            usize::MAX,
            1e-9,
        )
        .expect_err("this conflicts");
        assert_eq!(
            explain(
                &problem,
                &csr,
                &Clauses::default(),
                &mut scratch,
                conflict,
                &lb,
                &ub,
                usize::MAX,
                1e-9
            ),
            None,
            "a clause resting on an incumbent-earned bound must not be learned"
        );
    }

    /// **The property the whole mechanism rests on.** A learned clause is added to the
    /// model and every node below obeys it, so a clause that excludes even one feasible
    /// point silently deletes the optimum. Enumerate every binary point of a random
    /// model, learn from every branch that conflicts, and check each clause against every
    /// feasible point.
    #[test]
    fn a_learned_clause_never_excludes_a_feasible_point() {
        let mut learned_any = 0;
        for seed in 0..500u64 {
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
            let feasible: Vec<Vec<f64>> = (0..1u32 << n)
                .map(|mask| (0..n).map(|j| f64::from((mask >> j) & 1)).collect::<Vec<f64>>())
                .filter(|x| crate::heuristic::is_feasible(&problem, x, 1e-9))
                .collect();

            // Branch on every pair, which is where conflicts actually arise.
            for a in 0..n {
                for va in [0.0, 1.0] {
                    for c in (a + 1)..n {
                        for vc in [0.0, 1.0] {
                            let Some(clause) = learn(&problem, &[(a, va), (c, vc)]) else {
                                continue;
                            };
                            learned_any += 1;
                            for point in &feasible {
                                let satisfied = clause
                                    .iter()
                                    .any(|&(j, want)| (point[j as usize] > 0.5) == want);
                                assert!(
                                    satisfied,
                                    "seed {seed}: clause {clause:?} learned from \
                                     x{a}={va}, x{c}={vc} excludes the feasible point \
                                     {point:?}"
                                );
                            }
                        }
                    }
                }
            }
        }
        assert!(
            learned_any > 200,
            "only {learned_any} clauses learned, so little was checked"
        );
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
                        Err(_) => assert!(
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
