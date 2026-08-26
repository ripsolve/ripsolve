# Design notes

A record of what each decision in `ripsolve` was worth when it was made, including
the ones that were measured and rejected. It is kept separately from the README so
that the README can stay about using the solver.

Figures here are historical. Each was taken against whatever the instance set and
the rest of the solver looked like at the time, so a number below will not always
match the current benchmark table. `mkp_200` in particular names a harder generated
instance now than it did early in the project. Where a later measurement overturned
an earlier one, both are kept, because the reversal is usually the useful part.

## Why an LP relaxation rather than enumeration

The project began as an implicit-enumeration solver in the style of Balas, which
plateaus around 64 variables. Enumeration prunes with logical tests on partial
assignments; branch and cut prunes with a dual bound from a solved LP relaxation.
The difference is not incremental. On a 60-variable instance the enumerating solver
expanded roughly 4.7 million times as many nodes as the relaxation-based one. Every
later decision in this solver assumes a bound worth branching on, so the relaxation
is the foundation rather than a component.

## Presolve

Presolve reduces in place. A fixed column becomes `lb == ub`, and a redundant row
has its bounds freed. Nothing is renumbered, so there is no postsolve pass.

It is worth a great deal on structured models (`v006c*` is solved outright, and
`bin_10var_5con` drops from 24 nodes to 1) and nothing at all on the dense random
families in the benchmark set. That is not a shortfall in the implementation. A
commercial solver's presolve also reduces those instances to exactly their original
dimensions.

One bug is worth recording because it was invisible without a differential test.
Implied bounds were being rounded for every column, including continuous ones, which
cut off values those columns were entitled to take. A randomized comparison caught
an objective of 28.0 where the true optimum was 24.4.

## Cutting planes

Two families are implemented. Lifted knapsack covers are combinatorial and need a row
that reads as a knapsack. Gomory mixed-integer cuts come off the simplex tableau and
need no structure at all, which is what reaches the dense random rows that covers
cannot see.

### Cuts at the root do not pay

Early measurements were encouraging: adding GMI cuts to covers moved
`v128c1000n100` from a 63% gap after 60s to proven optimal in 6s, took
`v256c256n100` from 2290 nodes to 250, and `v128c256n100` from 368 to 1.

Every one of those was taken under most-fractional branching and depth-first search.
Re-measured after pseudocost branching and best-bound node selection were in place,
root cutting was slower on all eleven models in the target range:

| | no cuts | 3 rounds of 32 |
|---|---:|---:|
| `mkp_200` | 16.6s | 48.2s |
| `v064c1000n100` | 7.9s | 10.5s |
| `v064c200` | 1.5s | 2.3s |
| total, nine models | 32.1s | 67.8s |

The earlier wins were a good bound rescuing a bad search order. The search order is
no longer bad, so the bound has less to rescue. On `mkp_200` cutting even raised the
node count, from 72150 to 91346, and on MIPLIB's `markshare_4_0` it was the
difference between proving optimality in 21s and not proving it within 60s at all.

A separate bug was hiding underneath. The GMI density cap was absolute, at 40% of
columns with a floor of 30, so on models that are 99.5% dense it rejected every GMI
cut ever generated. Making it relative to average row support brought GMI back and
lifted the `v064c200` root bound from 82.7 to 95.4 against an optimum of 225.

### Managing the pool is most of the cost

Every separated cut was being added, permanently. Near-duplicate Gomory rows off
neighbouring tableau entries each took a row, and a cut that stopped binding after a
later round still rode along in every one of `mkp_200`'s seventy thousand node LPs.

Three changes address that. Candidates are ranked by efficacy, which is the distance
from the relaxation optimum to the cut's hyperplane. Unlike raw violation, efficacy
does not change when a row is rescaled, so it is comparable between cuts.
Near-parallel candidates are skipped, because the second of two cuts removing the
same region costs a row at every node and buys nothing. Cuts are aged out after two
consecutive slack resolves, and everything still slack is purged before the tree
opens.

The purge is free by construction. A row inactive at the optimum of a convex program
is not holding the bound up, so dropping it leaves the same point optimal. Node
counts come out identical across the purge, to the node, which is the check that the
reasoning is right.

Together these cost `mkp_200` four cuts instead of dozens and took it from 48.2s to
10.1s, faster than not cutting at all. Across the suite root cutting went from a
2.1x net loss to about 1.15x. Still a loss, so `cut_rounds` defaults to zero, but
`--cut-rounds N` is worth setting on knapsack-structured models.

### Two rejected follow-ups

**Capping the cut count.** If cuts pay when they are cheap to carry, a hard cap
should make them pay generally. Swept at 2, 3, 5, 8, 16 and 32 cuts per round, no cap
beats not cutting (38.3s against 44.7s at the best cap), and the existing default of
8 is already best on the one instance where cutting wins. The reasoning was wrong
about the mechanism, not just the number: a cap does not hold the cuts fixed and take
fewer of them. Keeping fewer rows in round one changes the vertex round two separates
from, so it changes which cuts exist at all. Capping at 2 yields five cuts on
`mkp_200` where capping at 8 yields four, a different and worse set.

**Orthogonality across rounds.** Selection compared candidates only against others
from the same round, so a round-two cut nearly parallel to one already held passed the
filter built to catch it. Since each round separates from the vertex the previous
round produced, near-copies are the common case. Checking against the held cuts too
works mechanically, dropping `v064c200` from 22 rows to 15 and `v064c1000n100` from
19 to 12, and is still 7% slower over the suite: three regressions against two gains,
with `v128c1000n100` going from 4.0s to 5.9s. At a 0.1 orthogonality bar, the
direction two cuts do not share carries more bound than the duplication costs.

### Cutting at nodes instead

Counting how many root cuts are still binding as the tree descends explains the root
result and points at the fix:

| instance | depth 1 | depth 3 | depth 10 |
|---|---:|---:|---:|
| `v064c200` | 36% | 12% | 1% |
| `v081c162n009` | 50% | 25% | 0% |
| `mkp_200` | 50% | 17% | 0% |
| `v256c256n100` | 50% | 25% | 16% |

Binding roughly halves per level. Weighted by where the nodes actually are, these
rows are carried through about 99% of the tree and bind in about 2% of it. Fifteen
of `mkp_200`'s 73436 nodes are at depth three or less. A root cut is a shallow-depth
object, and the shallow depths are a rounding error in the node count. That is the
whole explanation for a better root bound buying no smaller tree: on `v064c200`,
cutting lifts the root from 72.1 to 92.4 and the tree comes out at 2716 nodes against
2690 without.

Cuts read off a node's own tableau bind where they were made. They are valid for that
node's subtree only, so they never enter the shared model. The grown LP lives for one
solve, and the only thing that outlives it is the bound, which is valid everywhere
below that node and so is safe to prune and order children with. The node's own basis
and solution are left untouched, so branching still reads the relaxation it would
have read anyway.

| nodes | no local cuts | every 10th node | every node |
|---|---:|---:|---:|
| `v064c200` | 2690 | 2116 | 1136 |
| `v256c256n100` | 288 | 214 | 86 |
| `v064c1000n100` | 1106 | 786 | 458 |
| `mkp_200` | 72150 | 63896 | 40200 |
| total time | 24.75s | 22.98s | 36.50s |

Separating at every node halves trees but does not pay for itself. Every tenth node
is a 7% win outright and is the default, swept over 0, 1, 3, 10, 50 and 200.

### Growing a basis without refactorizing it

Adding a row to an LP would normally mean refactorizing, which at every separating
node is most of the cost. It does not have to be. A basis grown by `k` cuts is block
lower triangular against the one already factorized:

```text
    B' = [ B    0 ]        S = -I, the appended rows' own logicals
         [ R_B  S ]
```

so `B'^-1` is the existing `B^-1` plus a sparse rank-`k` correction. FTRAN
substitutes forward through the block, and BTRAN transposes it to upper triangular
and mirrors that.

The one trap is that the correction needs `B^-1` and not `LU^-1`, so the extension
wraps the whole base operator including its eta file. Pivots taken after the
extension cannot join the etas beneath it and need their own layer above. Collapsing
those two layers is silent until the first post-extension pivot, so there is a test
for exactly that case.

Measured against the same binary with the reuse forced to miss, the saving tracks row
count, which is the signature to expect when what has been removed is an `O(m*fill)`
factorization:

| rows | instance | reuse | refactorize | |
|---|---|---:|---:|---:|
| 30 | `mkp_200` | 22.21s | 23.84s | 1.07x |
| 200 | `v064c200` | 1.51s | 1.60s | 1.07x |
| 1000 | `v128c1000n100` | 4.07s | 5.42s | 1.33x |
| 1000 | `v064c1000n100` | 6.64s | 9.83s | 1.48x |
| | total, eight models | 35.90s | 42.29s | 1.18x |

## Node selection

The open node set is a depth-first plunge stack plus a best-bound pool. By default
the plunge length is zero, so the search is pure best-bound.

That default was measured rather than assumed. The textbook argument for plunging is
that depth-first reaches incumbents sooner and a child re-solves in a few pivots from
its parent's basis. The primal heuristics already supply incumbents, so the plunge
buys little and costs bound progress:

| instance | depth-first | best-bound |
|---|---:|---:|
| `v081c162n018` | 13570 nodes, 7.9s | 302 nodes, 0.5s |
| `v081c162n009` | 20584 nodes, 12.7s | 1286 nodes, 1.6s |
| `v064c200` | 17472 nodes, 14.8s | 2846 nodes, 2.9s |
| `v064c1000n100` | 77% gap after 60s | solved in 14s |

The trade is real in two places. Best-bound finds incumbents later, so a run that
hits its time limit reports a worse one, with `mkp_500` ending at a 3% gap rather
than 1%. It also holds every unexplored node in memory, where plunging keeps the open
set to roughly the tree depth. `plunge_limit` exists to raise when that binds.

### Diving, and why the default is still zero

MIPLIB's `graphdraw-gemcutter` is 166 columns and 474 rows, and the reference solver
closes it in 6.5 seconds where this one does not close it at all. It is small enough
that the failure is worth understanding.

The relaxation is not the problem: both solvers get 4310.0 against an optimum of 7118.5.
The failure is primal. Pure best-bound finishes 85% above the optimum, at 13176, and the
bound has crawled to 4926 after 107902 nodes.

Diving fixes the primal side and nothing else:

| dive length | incumbent |
|---|---:|
| none (pure best-bound) | 13176.5 |
| 10 | 9369.5 |
| 50 | 7893.5 |

and it costs the dense binary models 2 to 4x: `v081c162n009` goes from 0.62s to 2.40s,
`v064c064` from 0.03s to 0.13s. The two families want opposite things.

An adaptive dive length was implemented to have both: dive while dives are finding
incumbents, back off when they are not, with a periodic full-length probe so a model
that needs long dives can still discover it. It protects the dense models almost
completely, and it does not work well enough to keep. Part of the gain is lost, 8336
against 7893 at the same ceiling, because the models that need diving need it
*sustained*, and a rule that withdraws it after an unproductive dive is withdrawing it
from exactly the models it was meant for. Diving also never closes `graphdraw-gemcutter`
at any setting, so the cost to the target class bought no additional solved instance.

`plunge_limit` therefore stays at zero and stays a fixed length when set, so a caller
who asks for diving gets what they asked for. The finding this leaves behind is that the
gap on models like this one is a primal-heuristic gap, not a search-order gap: diving
helps only because it is acting as a heuristic, and the honest fix is a better heuristic.

## Branching

Pseudocost branching scores a column by the objective degradation it has actually
caused rather than by how fractional it looks. It is a large win on the hardest
instances and roughly neutral elsewhere.

Strong branching was tested twice, once under depth-first and again under best-bound.
It reduces node counts by 10 to 32% and does not pay for the probes either time, so
`strong_branching_budget` defaults to zero.

## The basis factorization

The basis is a sparse LU with Markowitz pivoting and threshold partial pivoting, plus
a product-form eta file for per-pivot updates. It replaced a dense explicit inverse
costing `O(m^2)` per solve and `O(m^3)` to rebuild, which at `m = 1000` was 0.2
seconds per node.

On the 1000-row models, in a fixed 100-second budget:

| | dense inverse | sparse LU |
|---|---:|---:|
| `v064c1000n100` nodes | 1,370 | 13,451 |
| `v064c1000n100` simplex iterations | 6,594 | 49,245 |
| `v064c1000n100` remaining gap | 93.7% | 51.6% |

That is roughly 7.5x the simplex throughput. Smaller models gain between 1.2x and
2.4x in wall clock.

Where the remaining time went was rebuilding factors that were already correct. A
warm re-solve performing zero pivots cost 9.1ms on a 1000-row model against 10.9ms
for a real node doing eleven, so 84% of node time was setup. A child's basis is
identical to its parent's, because bounds do not enter the basis matrix, so the
factors are reusable verbatim. Each LP now keeps its recent factorizations:

| | before | after |
|---|---:|---:|
| `v128c1000n100` | 9.2s | 4.2s |
| `v064c1000n100` | 10.9s | 9.2s |
| `v256c256n100` | 0.48s | 0.33s |

The cache holds several entries rather than one, because best-bound selection does
not visit the tree in an order that keeps a single entry warm. One entry measured an
8 to 20% hit rate. What recurs is siblings, which share a parent's basis.

Two things were measured and not adopted. Forrest-Tomlin updates would shorten the
eta file, and the `Basis` interface was built so the swap needs no change to the
simplex, but sweeping the refactorization interval over an 80x range moves solve time
by only 3 to 7%. Neither replaying the eta file nor rebuilding the factors is where
the time goes. Equilibration scaling also gained nothing, because threshold partial
pivoting is already scale-invariant per column.

## Primal heuristics

Branch and bound cannot prune anything until it holds a feasible solution, so finding
one early matters independently of the bound. Three are tried, cheapest first:
rounding the relaxation, diving, and a feasibility pump.

Diving is the wrong tool for this instance family and fails on every model in it. It
commits to a rounding and re-solves a smaller LP each step, so where the feasible set
is sparse it walks into infeasibility with no way back. A one-level backtrack was not
enough to save it. The feasibility pump never fixes anything, alternating instead
between rounding and re-optimizing the original constraint set under a distance
objective, so its LP stays feasible throughout. It finds solutions on five of the six
models where diving finds none.

The pump was given restarts, and they did not do what they were meant to. Three
defects turned up on the way and are worth separating from the result.

It only noticed a cycle of period one, by comparing against the previous target, so a
longer cycle ran out the round limit doing nothing. It now fingerprints every target
since the last restart.

Small flips escape a short cycle but not a basin: once the walk keeps returning to the
same neighbourhood, nudging a handful of columns returns it there again. After three
cycles it now restarts, moving a tenth of the integer columns and forgetting where it
has been. The perturbation is a seeded SplitMix64 so a run stays reproducible.

The third was a plain bug. A re-solve that hit the per-solve iteration cap ended the
whole pump, and the distance LP is a full re-optimization rather than the handful of
pivots a dive step takes. On `nursesched-sprint02` and `piperout-27` the first round hit
the cap and the heuristic returned having done nothing, which from outside is
indistinguishable from a pump that ran and failed. It now restarts on an unfinished
re-solve and gives up only after three in a row.

None of that finds a first solution on the models it was aimed at. `neos-555001` runs
10000 rounds and 81000 pivots without reaching feasibility; the other two run properly
now and still fail. The pump is bounded to eight re-solves' worth of work in total for
that reason, because a heuristic that fails should not first spend 26 seconds of a 45
second budget doing it.

What this leaves is that these models need feasibility machinery of a different kind,
propagation-based rather than rounding-based, and that the pump's failures here are
genuine rather than an artefact of it being throttled.

Those solutions are poor, 2791 against an optimum of 137 on `v064c064`, but an
incumbent of any quality switches pruning on. `v064c200` drops from 9916 nodes and
5.7s to 3388 nodes and 2.25s.

None of those three improves a solution; they only find one. That is the gap
`graphdraw-gemcutter` exposed: the search reached 13176 against an optimum of 7118 and
sat there, not for want of nodes but for want of anything looking near a good solution
rather than near the relaxation.

The improvement search fixes every integer column where the incumbent and the current
relaxation already agree, and turns the search loose on what is left with a node budget
and the incumbent as a cutoff. Two points agreeing on a column is weak evidence that a
good solution has it there, and weak evidence over hundreds of columns leaves a model
small enough to search properly. It is Danna, Rothberg and Le Pape's RINS, and it reuses
the search rather than adding a mechanism.

| | before | after |
|---|---:|---:|
| `graphdraw-gemcutter` (optimum 7118.5) | 13176.5 | 8150.5 |
| `australia-abs-cta` (optimum 106.9) | 10865 | 2332.6 |

The target class is unaffected, because there the heuristics already supply a good
incumbent and the neighbourhood searches finish almost immediately: `v064c200` 1.19s,
`mkp_200` 14.5s, both unchanged.

What it does not do is find a first solution. On the instances where nothing is found at
all, `neos-555001`, `nursesched-sprint02`, `piperout-27`, `hypothyroid-k1`, it has no
incumbent to improve and changes nothing. Improvement and discovery are separate
problems, and only one of them is addressed here.

In-tree attempts are scheduled adaptively rather than on a fixed cadence. The
interval doubles after an attempt that finds nothing and snaps back after one that
succeeds. Because diving fails on whole instance families rather than the occasional
node, the wasted attempts were measurable: running unconditionally cost
`v064c1000n100` its incumbent quality. Backing off leaves the tree identical and
removes the overhead, taking `v128c1000n100` from 13.3s to 9.8s and `v081c162n009`
from 1.7s to 1.4s at unchanged node counts.

## Parallelism

The tree search runs across worker threads sharing one node pool and one incumbent.
Presolve, cut generation and the root heuristics run once, before any thread is
spawned, so only the tree is parallel.

Node counts vary between runs, because which node a worker takes depends on timing.
The answer does not vary. Every bound and cut is globally valid and every worker
prunes against the shared incumbent, so the proven optimum is the same however the
work is divided. That is asserted across every sample at 2, 4 and 8 threads.

| | 1 thread | 4 threads | 16 threads |
|---|---:|---:|---:|
| `mkp_200` | 8.3s | 3.3s | 2.6s |
| `mkp_500` gap after 120s | 3.07% | 3.03% | 0.30% |
| `mkp_500` simplex iterations | 834k | 3.27M | 6.70M |

Throughput scales about 8x on 16 threads. The shortfall against linear is inherent
rather than incidental: workers expand nodes that a serial search would have pruned
against an incumbent it had already found, so total node count rises with thread
count even as wall-clock time falls. `mkp_200` goes from 19522 nodes to 49432.

## Correctness practices

Several bugs in this project were invisible to unit tests and were caught only by
differential testing against another solver. Two patterns recur.

Tests that cannot fail. An early set of LU scale tests gave every basis a strong
diagonal, which made pivot choice irrelevant. They passed with the bad pivot search
deliberately reintroduced. Any test asserting that an optimization is correct should
be checked against the unoptimized path, and any threshold should be checked by
disabling it and confirming the test then fails.

Circular verification. A check written to confirm two pivot orderings agreed compared
a shortlist against a stable sort of the same ordering it came from, so it could not
detect the difference it existed to find. Ground truth has to come from somewhere the
code under test did not produce it.
