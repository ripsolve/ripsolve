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

### Where the remaining throughput goes

The MIPLIB instances that find no solution are not short of heuristics, they are short of
nodes: 1, 17, 50 and 350 of them in twenty-five seconds. Profiling `neos-555001`, 3474
rows, puts the time in the triangular solves and in pricing, roughly a fifth each, with
factorization and node setup behind them.

The factors are the surprise. Measured at the point of refactorization, they carry
essentially no fill:

| instance | rows | factor nonzeros | basis nonzeros | fill |
|---|---:|---:|---:|---:|
| `v064c1000n100` | 1000 | 2986 | 2986 | 1.00 |
| `neos-555001` | 3474 | 3861 | 3856 | 1.00 |
| `piperout-27` | 18442 | 50837 | 50441 | 1.01 |

So the solves are not doing arithmetic, they are walking. Each `ftran` makes four passes
of length `m`: permute in, solve `L`, solve `U`, permute out. The numeric work along the
way touches about `m` nonzeros in total, and a single entering column has a handful. The
iteration rate follows the row count and not the nonzero count, dropping 14x between
3474 rows and 18442 while the model grows 5x.

The obvious fix is a hyper-sparse triangular solve: find the positions reachable from
the right-hand side's nonzeros and visit only those, rather than scanning `0..m` four
times. It was implemented, and it is slower.

The reasoning that motivated it was wrong, and the mistake is worth keeping. Near-zero
fill does not imply that the reachable set from a sparse right-hand side is small. The
factors being sparse says each pivot touches few others; it says nothing about how far a
nonzero travels along the chain. Measured, `B^-1 a_q` comes out at 32.6% of the rows on
`piperout-27` and 12.4% on `neos-555001`, averaged over every solve of the relaxation.
At those densities the reachability walk and the sort that orders it cost more than the
scan they replace:

| | dense | sparse |
|---|---:|---:|
| `neos-555001` relaxation | 0.48s | 0.56s |
| `piperout-27` relaxation | 27.6s | 30.5s |

Iteration counts were identical, so the sparse path was numerically exact and simply not
worth its bookkeeping. It is reverted.

That leaves the dual computation, a dense BTRAN of the basic costs on every iteration,
and pricing, which is `O(n)`. Reduced costs updated from the pivot row remove the first
and shrink the second, so that was implemented too, and it is also slower.

The premise checked out this time. A pivot row needs `B^-T e_r`, and a *unit* vector's
solution really is hyper-sparse where a column's is not: 0.4% of the rows on
`piperout-27` and 2.7% on `neos-555001`, against 32.6% and 12.4% for `B^-1 a_q`. Reading
the pivot row from the row-wise model, restricted to those rows, touches a few hundred
entries where pricing every column touches every nonzero. What went wrong is elsewhere:

| | pricing as it is | reduced costs maintained |
|---|---:|---:|
| `mkp_200` | 11.2s | 13.4s |
| `v064c200` | 1.69s | 1.75s |
| `neos-555001`, nodes in 25s | 1289 | 443 |

Iteration counts were identical on the first two, so the pivots were the same and the
difference is bookkeeping.

Three things are worth taking from it. Three quarters of `piperout-27`'s iterations are
phase one, whose costs are the gradient of the bound violations and change whenever a
step changes which basic variables are infeasible, so they cannot be carried across a
pivot at all: the technique only ever applies to the other quarter. The saving is
smaller than it looks, because this solver's models are extremely sparse, 1.77 nonzeros
per column on `piperout-27`, so a reduced cost is a two-element dot product and reading
one from an array is barely cheaper than computing it. And against that, maintaining
them adds a BTRAN, a row-wise product and pattern bookkeeping to every pivot.

### Devex pricing

Devex attacks the other axis: not the cost of an iteration but how many are needed. It
scores a column by `d^2 / w`, where the reference weight `w` approximates how far the
objective actually moves per unit of the step, so that pricing stops preferring columns
that look steep only because they are badly scaled and then move almost nowhere.

It was implemented, with the pivot row read the cheap way described above and a
reference framework that restarts when the weights drift above the true column norm.
The result is a coin flip with a bad tail:

| | Dantzig | devex |
|---|---:|---:|
| `nursesched-sprint02` | 20167 pivots, 15.5s | 10000 pivots, 7.9s |
| `hypothyroid-k1` | over 200s | 92.8s |
| `v064c200` | 59 pivots | 46 pivots |
| `neos-555001` | 3955 pivots, 0.47s | 4376 pivots, 0.73s |
| `piperout-27` | 15679 pivots, 27.5s | 16107 pivots, 38.6s |
| `neos-3075395-nile` | 30827 pivots, 24s | over 200s |

Halving the pivots on one model and turning a 24 second solve into a timeout on another
is not something to ship, even behind an option, while the reason for the bad case is
unknown. What is known is that the first attempt was far worse still, 3955 pivots to
24299, because the weights only ever grow and nothing restarted the framework; and that
the drift check has to be periodic, since a pass over `alpha` on every pivot costs more
at 27756 rows than the rule saves. Neither of those explains `nile`, which stayed slow
with both fixed.

Taken with the sparse-solve result, the picture is that the LP here is near a local
optimum for the models it sees. The classical sparse-simplex techniques all trade
per-column arithmetic for bookkeeping, and on a matrix this sparse there is not enough
per-column arithmetic left to trade. Getting further would mean changing what is
computed rather than how, and both attempts in that direction have now failed for
different reasons: partial pricing returned a wrong answer, and devex is unpredictable.
What they share is that both give up Dantzig's rule, which is doing more work here than
choosing a column: it is also the reason the chosen column is numerically safe, and the
reason the iteration count is stable across models. A replacement needs to supply both,
not just a better score.

### Re-measuring pricing once the ratio test was sound

The conclusion above, that Dantzig's rule was doing double duty as a numerical
safeguard, turned out to be half right in a way worth correcting. The safeguard was
real, but it was compensating for a defect in the ratio test rather than being a
property of the pricing rule.

The ratio test accepted any pivot above the absolute tolerance, including an entry a
billionth the size of the largest in its own transformed column. On MIPLIB's
neos-850681 one such pivot took the basis inverse from entries of 1e4 to 1e15 in two
iterations, after which the solve reported a feasible relaxation infeasible. Dantzig's
rule avoided the worst of that only by tending to pick columns that did not provoke it.
So both earlier pricing rejections rested on evidence gathered over a broken ratio
test: partial pricing's wrong answer, and devex's unexplained bad case on `nile`, are
both the signature of that bug.

With the guard in place devex was implemented again and re-measured. It is no longer
erratic, and the objective it reaches now agrees with Dantzig's everywhere, but it is
not worth shipping for a plainer reason: the pivot row costs a BTRAN and a pass over
the nonbasic columns, which roughly doubles what an iteration spends on pricing, and
the iteration count does not fall nearly enough to pay for it.

| LP relaxation | Dantzig | devex |
|---|---:|---:|
| `misc07` | 464 pivots, 25.8ms | 201 pivots, 10.7ms |
| `neos-3610173-itata` | 2496 pivots, 238ms | 1280 pivots, 74ms |
| `decomp1` | 4210 pivots, 1.08s | 3817 pivots, 1.40s |
| `neos-1445532` | 4562 pivots, 794ms | 4089 pivots, 1.06s |
| `neos-1582420` | 3730 pivots, 2.97s | 4605 pivots, 4.49s |
| `n13-3` | 1811 pivots, 111ms | 2185 pivots, 198ms |
| `neos-850681` | iteration limit | iteration limit |
| `s55` | iteration limit | iteration limit |

Two models gain outright and the rest lose wall time to the pivot row. Crucially it
does nothing for the three relaxations that do not finish at all, which is what it was
re-tried for.

### What the iteration counts are actually being compared against

A correction to the numbers quoted elsewhere here. Comparing our relaxation solves
against a reference solver's default settings compares against a *presolved* model,
which is not the same input at all. With the reference solver's presolve turned off:

| LP relaxation | ripsolve | reference, presolve on | reference, presolve off |
|---|---:|---:|---:|
| `p200x1188c` | 454 | 46 | 525 |
| `dsbmip` | 1733 | 1623 | 1719 |
| `decomp1` | 4210 | 2230 | 3268 |
| `neos-1445532` | 4562 | 1197 | 3320 |
| `n13-3` | 1811 | 327 | 1051 |
| `neos-1582420` | 3730 | 1178 | 1165 |
| `neos-595904` | 5756 | 660 | 742 |

Like for like the simplex is within 1.0 to 1.4x on several of these and ahead on one,
rather than the 2 to 9x the default comparison suggests. The gap that comparison was
really measuring is presolve, and the relaxation path here does not presolve at all.

That does not excuse the models that never finish, because presolve is not what rescues
those either: without it the reference solver still takes 4502 iterations on
`neos-850681`, 10422 on `s55` and 23666 on `gasprod1-2`, all of which it finishes.

### Phase 1 is where the real gap is

Splitting the iterations by phase locates it. `s55` spends every one of 60000
iterations in phase 1 and never reaches a feasible basis at all, on a model the
reference solver solves outright in 10422; `neos-850681` spends 54%. Phase 1 also
dominates several relaxations that do finish: 92% of `decomp1`, 81% of `n13-3`, 73% of
`misc07`.

The reason is structural. Phase 1 shares the phase-2 ratio test, which stops at the
first breakpoint. That fixes at most one infeasibility per pivot, and at a degenerate
vertex, where many basic variables sit exactly on their bounds and the first breakpoint
is at zero, it fixes none and moves nothing. The textbook answer is the long-step rule:
the phase-1 objective is piecewise linear and convex in the step, falling at the
entering column's reduced cost and flattening by `|beta_i|` each time a basic variable
crosses a bound, so the minimum is the first breakpoint at which the slope reaches
zero, which may be many breakpoints along.

It was implemented and reverted. The rule itself works: on `lu_pivot_regression` the
sum of violations falls from 1.35e4 to 81.25 within 2000 iterations. It then cycles,
and the model goes from 622 iterations to not finishing in 500000.

What the cycle exposes is worth keeping even though the change is not. The anti-cycling
here counts steps of length zero, which is a proxy for the objective not moving, and
the short-step rule makes the two coincide. The long-step rule breaks that equivalence:
it takes a full step of 0.157 every iteration, trading one bound for another, and
leaves the total violation at exactly 8.1254231825e1 for as long as it is allowed to
run. Because the step is not zero, Bland's rule is never reached. Measuring objective
progress rather than step length makes the detector fire, and was tried, but by then
the long steps have walked the basis somewhere the short-step rule cannot recover from
either, so it still does not finish.

The rule was then finished properly on the `longstep-phase1` branch, which is kept
unmerged. Two real defects in it turned up, both worth knowing about:

- A crossing that computes to a small *negative* step was being discarded rather than
  clamped to zero. A basic variable a hair outside a bound, but not far enough outside
  to count as violating it, is pushed further out the moment the step begins, so its
  kink belongs at zero; dropping it loses that slope entirely. The short-step rule
  already clamps, and says why in a comment. Without the clamp the walk starts uphill
  while believing it starts downhill: at one vertex the true slope is +0.49 where the
  reduced cost alone reports -0.56.
- The slope is a reduced cost plus a running sum of `|beta|` terms of the same size, so
  at a degenerate vertex it reaches zero by exact cancellation. Tested against a hard
  zero it misses by an ulp and the walk steps straight past the minimum. On `s55` a
  slope of -1 met a first crossing of weight 1 at a step of zero, and the walk went on
  to a step of 1.6e-3 that improved nothing, indefinitely.

With both fixed the rule is correct and does what it claims, and it still is not worth
shipping:

| LP relaxation | short step | long step |
|---|---:|---:|
| `misc07` | 464 | 409 |
| `lu_pivot_regression` | 622 | 694 |
| `decomp1` | 4210 | 4307 |
| `n13-3` | 1811 | 2106 |
| `s55` | never finishes | never finishes |

Breakpoints per pivot were not the bottleneck. It also needs the stall detector to
measure the objective rather than the step length, since it takes full-length steps
that improve nothing, and making that change causes Bland's rule to engage earlier and
leaves `s55` stuck at a worse point than before, a violation of 9.875 against 5.225.

### Where phase 1 actually gets stuck, and what it points at

What `s55` does is now clear. Phase 1 removes 99.99% of the infeasibility, 9.7e4 down
to about 5, within 30000 iterations, and then stops dead at a degenerate vertex with
three violating rows out of 9892. Every column pricing likes there has a true
directional derivative of exactly zero, so nothing enters that can improve the
objective, and Bland's rule walks between bases at the same point taking steps of
length zero.

The reason pricing likes them is worth stating, because it is not a tuning matter. The
phase-1 cost vector scores a basic variable by whether it currently violates a bound,
so it is blind to the kink at a variable sitting exactly *on* one, and blindest of all
at a fixed variable, where the two bounds coincide and any movement at all is a
violation. On `s55` the entering column's reduced cost claims a slope of -1 where the
true one-sided slope is 0, and the whole of that difference is one fixed basic variable
whose bounds are `[0, 0]`. A smooth reduced cost cannot see a kink, so no ratio test,
short, long or Harris, can repair a column choice made this way.

That is the ceiling on the primal method here, and it is why the reference solver does
not hit it: its default is the dual simplex, which has no primal phase 1 to get stuck
in. It solves `s55` in 10422 iterations and `neos-850681` in 4502 with its presolve
turned off. A dual simplex already exists here and is used for warm starts; what it
lacks for a cold start is a dual-feasible basis to begin from, which for a bounded
model is mostly a matter of parking each nonbasic column on the bound that gives its
reduced cost the right sign. That, rather than another primal ratio rule, is where the
next attempt should go.

### Entering cold solves through the dual method

Since phase 1 is where the primal method gets stuck, and the dual method has no phase 1
to get stuck in, the obvious move is to start there. It was built on the
`dual-cold-start` branch, which is kept unmerged.

The starting basis costs nothing to work out. With every logical basic the basis is the
identity and their costs are zero, so the duals are zero and each structural's reduced
cost is exactly its objective coefficient; dual feasibility is then a sign condition per
column, no factorization needed. A column that costs something to increase belongs at
its lower bound, one that pays to increase at its upper, and a column whose cost points
at a bound it does not have rules the method out for that model. Row selection is dual
steepest edge, with the chosen row's weight taken exactly from the pivot row already in
hand, `rho . rho`, rather than carried forward: the carried value took
`drayage-100-23` 487230 iterations where the exact one takes 2374.

Three defects turned up, and the first of them is now fixed on the main line because it
was never about cold starts at all.

- The dual method counted stalls by primal step length. Degeneracy there is a zero
  *ratio*, meaning the entering column was already priced at zero and the dual objective
  cannot move whatever the primal step is. On `drayage-100-23` the dual objective sat at
  240.4124 for half a million iterations with the stall counter never leaving zero.
- Steepest edge weights could grow until the score underflowed, and a selection rule
  requiring a strict improvement over zero then skipped every violating row. That
  reported a basis with a violation of 1.6e5 as optimal, returning -1.0 for a relaxation
  whose optimum is 2087. A violating row must never be passed over, whatever it scores.
- Primal feasibility means optimality only if the basis is also dual feasible, which the
  ratio test preserves in exact arithmetic and drifts out of in this one. Assuming it
  returned wrong optima on `gasprod1-2` and `s55`. Checked instead, and where the
  invariant has lapsed the dual method now hands its basis to the primal loop rather
  than ruling on it, which makes the two a sequence rather than a choice.

The last two were caught by checking every relaxation objective against a reference
solver rather than reading iteration counts, which is the only reason they did not
survive. Eighteen of the twenty now agree; `gasprod1-2` solves for the first time, and
`s55` and `neos-850681` no longer answer wrongly, they just do not finish.

On relaxations the result is genuinely mixed: `neos-1582420` 3730 pivots to 937,
`neos-595904` 5756 to 1473, `n13-3` 1811 to 1257, against `decomp1` 4210 to 8724,
`neos-1445532` 4562 to 22197 and `cap6000` 810 to 6216.

On the search it is worse, and that is what decides it. The same four instances solve,
`mik-250-20-75-3` improves from 2.95% to 1.16% and two instances that had found no
incumbent at all now find one, but most gaps widen: `decomp1` from 20.7s to 53.0s,
`n13-3` from 21.8% to 37.1%, `neos-911970` from 48.1% to 77.6%.

The tempting explanation, that node LPs are warm starts priced on steepest edge weights
of one that are exact only at the identity basis, is wrong: restricting the rule to cold
starts changes nothing. What is left is that the root relaxation is solved once and
everything downstream reads its basis. Gomory cuts come off that tableau and branching
reads its fractional values, so ending on a different optimal vertex redirects the whole
search. The dual method reaches a different vertex, and on this set the vertex it
reaches is usually the worse one to start from.

That is a real finding about where the leverage is, and it is not in the relaxation.
Getting the root LP to a better *vertex* rather than a better *objective* is the
question the branch leaves open.

### The Harris ratio test

Those three relaxations stall rather than mislead: thousands of consecutive
zero-length steps, with the objective on `neos-850681` creeping from 2087.475 towards
2087.0 over hundreds of thousands of pivots. That is degeneracy, and the textbook
answer is Harris's two-pass ratio test. Pass one finds the furthest the entering
variable could move if every basic variable were allowed to overshoot its bound by the
feasibility tolerance; pass two takes the largest pivot within that band. It should
both stabilize the basis and let a step be nonzero where the exact rule can only stall.

It was implemented and reverted. It failed `a_feasible_relaxation_is_not_reported_infeasible`,
the regression test kept for exactly this class of bug, and it left `neos-850681` at the
iteration limit and slower than before. The mechanism is inherent to the plain form:
the chosen step is the selected row's breakpoint rather than the shortest one, so every
pivot may leave a basic variable up to a tolerance outside its bounds, and those
violations accumulate until phase 1 gives up on a feasible model. Controlling that is
what bound shifting is for, and Harris without shifting is not a partial version of the
technique but a broken one.

So the degeneracy on those three relaxations remains open, and the next thing to try is
shifting itself rather than another pricing or ratio rule.

One smaller thing was measured and rejected on the way. Both solves allocated an
`m`-length buffer per call, which on `neos-555001` is 128000 allocations of 27KB and
showed as 5.8% in `malloc`. Reusing a thread-local buffer removed the allocations and
changed the iteration rate by nothing measurable, so it was not kept: the allocator
returns the same size class each time and the cost was never really there.

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

### Why the pump does not rescue the starved instances

Eight of the twenty instances in the tractable set finish with no incumbent at all, so
no bound is worth anything on them. Diving is what fails there: given room to run it
dead-ends rather than running out of steps, on `pigeon-08` with two columns left to
decide and on `neos-1445532` after 1287 successful fixings. Backtracking the dive over
its own trail was implemented and did not change the outcome.

That leaves the feasibility pump, which is the tool for exactly this case because it
fixes nothing and so cannot dead-end. It ran once, at the root, and only if rounding
and diving had both already failed there. Running it again from node relaxations, so
that it starts somewhere other than the root, was implemented and measured and produced
no incumbent on any of the six.

The reason is that the pump is not near-missing, it is not converging. Instrumenting how
close its rounded point ever comes to feasible:

| relaxation | closest snapped point, 60 rounds | at 600 rounds |
|---|---:|---:|
| `neos-1582420` | violates a row by 0.5 | 0.2 |
| `neos-1445532` | 27.0 | 18.0 |
| `neos-595904` | 187.0 | 119.0 |

Ten times the rounds barely moves it, so the round limit is not what binds. Three
explanations were tested and none holds: it is not the iteration budget, it is not the
starting point, and it is not the linearization, which is only valid for binary columns
but which every one of these models is, `neos-1582420` aside with a hundred general
integers of span thirteen. The walk reaches a basin and the restart mechanism does not
leave it.

What the literature offers from here is a pump that carries a decaying fraction of the
real objective, and a third stage that takes the closest point and searches a small
model around it for the remaining violation. Both are real work rather than tuning, and
neither is a small change. Until one of them exists, these instances have no route to a
first solution, and their bounds cannot be spent.

### Presolve cannot pay until the model is compacted

Comparing our presolve against a reference solver's on the tractable set says we leave a
great deal on the table: it keeps 9% of the columns of `neos-1445532` and 33% of the
rows of `n13-3`, where we remove seventy-six columns from the first and nothing at all
from the second. Counting structure says where that comes from. Doubleton equalities,
which let a variable be substituted out, number 530 on `n13-3`, 127 on `neos-1445532`
and 90 on `neos-850681`; duplicate rows number 1586 of `decomp1`'s 8357.

Parallel row merging was implemented for the second of those, since it needs nothing
that is not already here: one row a scalar multiple of another, so the survivor takes
the intersection of both bound sets and the other is freed. It works, folding away 1565
rows of `decomp1` and 59 of `dsbmip`, and it made the whole set slightly and uniformly
slower. `decomp1` itself went from 20.8s to 25.6s having shed a fifth of its rows.

The reason is that freeing a row does not remove it. `Lp::relaxation` takes `m` from
`problem.n_rows()` and clones the matrix whole, so a row whose bounds have been widened
to infinity still becomes a logical variable, still occupies a row of the basis and
still enters every factorization. The reduction costs presolve time and returns nothing.
That applies equally to the redundant-row and forcing-row reductions already here: they
have been cosmetic at the LP level all along.

So the order of work looked the other way round from how it appears: compacting the
model, dropping freed rows and fixed columns behind a map back to the original indices,
looked like the prerequisite any row reduction needs before it can pay.

It was built, measured over all 140 pure binary instances, and reverted. It works. A
fixed column's contribution moves into the row bounds and the objective offset, a row
left empty and violated is reported rather than discarded, and the expansion puts the
solution back; `cap6000` sheds 62% of its nonzeros, five dense cardinality rows across
all six thousand columns, and `mitre` loses 392 rows. What it does not do is help:

| | closed | total time on those |
|---|---:|---:|
| without compaction | 24 of 140 | 151.2s |
| with compaction | 24 of 140 | 172.5s |

Not one instance changed hands, and the models that do close take 14% longer overall.
`decomp1` improves by 31% and `cap6000` by 12% against `irp` at 71% worse and `decomp2`
at 34%, and the `decomp1` and `decomp2` directions reproduce single threaded, so this is
the change and not the parallel search's noise. Removing rows alters the basis ordering,
which alters cut generation and branching, and the search that follows is simply a
different one; some models come out ahead and more come out behind.

The reason the systematic gain never arrives is that there is nothing to compact.
`air04`, `air05`, `cod105`, `neos-820879` and `eil33-2`, all of them models the
benchmark says we lose, come out of presolve at exactly their original size: no row, no
column, no nonzero removed. Compaction can only remove what presolve has already
finished with, and on the models that matter presolve has finished with nothing.

That leaves the dependency intact but the order reversed. Compaction is still what
doubleton aggregation needs, since substituting a variable out is a column removal, and
the work is in the history to bring back with it. What it is not is a throughput fix in
its own right, and shipping it on that expectation would have cost 14% for nothing.

### Why the pure binary instances are not lost for one reason

The pure binary benchmark is 140 instances, every column of every one confined to
{0, 1}, and this solver closes 24 of them against 46 for HiGHS and 45 for SCIP. On the
24 it does close it is usually the quickest of the open-source field, median rank one,
so the machinery is not the difficulty; coverage is.

Diagnosing the 33 it loses that an open-source solver closes, thirteen of fourteen
sampled end with no feasible point at all. That looked like one problem with one answer,
and it is not. Counting nodes explored in a minute separates it into three:

- `ex9` and `ex10` reach **one node**, with or without cuts. The root relaxation alone
  spends the whole budget; at 517112 and 1162000 nonzeros they are simply larger than
  this LP is quick on.
- `mitre` reaches **one node with root cutting on and 5126 with it off**. Nothing is
  wrong with its LP; the root cut loop, fifty rounds each re-solving the model, eats the
  minute before the search starts.
- The `acc-tight` family reaches twenty-odd nodes either way and finds nothing in them.
  That is the primal failure the count of starved instances suggested.

Three causes, and each needs its own answer. Four attempts aimed at the wrong one:

| attempt | what it did | result |
|---|---|---|
| clique cuts | 1296 to 10640 conflict edges found | shipped, no instance closed |
| implication closure | neos18 182 triangles to 55245 | reverted, bound *worse* |
| RENS | held 96.9% of air05's columns | reverted, one instance gained an incumbent |
| probing with compaction | ex9 1867 columns and 7208 rows removed, LP 22% smaller | reverted |

The last is the sharpest lesson. Probing found real reductions where nothing had before,
compaction carried them into the LP, and the node count did not move: `ex9` explores one
node either way, 10970 simplex iterations against 11036. A model 22% smaller does not
help a search that never gets past its first node. Two instances that do close, `decomp2`
and `cap6000`, came out 45% slower for it.

The order to work in follows from the split rather than from the totals. `mitre` says the
root cut loop needs a budget, since spending an entire run on it is worse than not cutting
at all. `ex9` and `ex10` say the LP needs to be quicker on large models before anything
built on top of it can matter. Only the `acc-tight` family is the primal problem that all
four attempts above were aimed at.

### What stops ex9 and ex10, and why nothing above the LP can help them

These two reach one node in a minute whether cutting is on or off, and the reason is
that their root relaxation never finishes. It does not finish in 900 seconds either,
against 82s for a reference solver on ex9 and 238s on ex10 with its presolve turned off.

Tracing the solve says where it goes. Neither model gets out of phase 1:

```text
iter 0      phase1 true  worst 1.0000e0  stalled 0
iter 20000  phase1 true  worst 1.0000e0  stalled 8871
```

The worst bound violation does not move from 1.0 across twenty thousand iterations, and
8871 of them in an unbroken run take a step of length zero. This is the phase-1
degeneracy recorded above, in its plainest form: a crowd of bases describing one point,
with every step between them going nowhere.

Two things were tried against it and neither worked. Entering through the dual method,
which has no phase 1 to get stuck in, does not finish ex9 in 400 seconds either. Bound
perturbation, triggered on a degenerate run of 2000 where healthy relaxations peak at
296, fires where it should and then the loosened solve does not finish either.

Worth noting against the earlier perturbation experiment, which took `neos-850681` from
not finishing to 11 seconds: that one perturbed the model and solved it cold, while this
one perturbs and warm starts from the stalled basis. The difference has not been chased
and may be the whole of it.

The conclusion that does hold is about ordering. Anything built above the relaxation,
cuts, heuristics, presolve reductions, branching, can only act on instances whose
relaxation is solved, and for this pair it never is. That is what made the probing and
compaction work on ex9 pointless in advance: a model 22% smaller still has to get
through phase 1, and 22% off a solve that does not finish is a solve that does not
finish.

Also worth recording: a reference solver never solves this relaxation either. It closes
ex9 and ex10 as MIPs in zero simplex iterations, its presolve collapsing the models
before any LP is attempted. Matching it here is a presolve problem wearing an LP
problem's clothes, and the presolve it would take is well beyond what probing reached.

### Reading conflicts out of long rows

The conflict graph was built from two-column rows alone, and the cost of that was not
visible until it was measured. On `air04`, 823 set partitioning rows over 8904 columns,
it held four edges and no triangle at all. A clique needs three literals, so the clique
separator had nothing to work with and emitted nothing. The same held on `air05`,
`cod105`, `mitre`, `chromaticindex32-8`, `neos-913984` and `graph20-20-1rand`: an empty
graph, a shipped separator, and no cuts from it on any of them.

A row `sum x_j = 1` is the richest conflict a binary model has, every pair in it
excluded. Skipping it was a deliberate guard against the quadratic blowup of writing
those pairs out, one row here spanning six thousand columns, and the guard was aimed at
the right hazard with the wrong remedy. Holding the clique whole rather than as edges
costs the row's length instead of its square, and reading a row for its longest
overshooting prefix costs a sort.

The graph fills as expected: `air04` from 4 edges to 4560125 across 821 cliques, `air05`
to 4673004, `cod105` to 1576960, `mitre` to 144774.

What that bought is much less than the numbers suggest, and the reason is worth keeping.
On `air04`, `cod105` and `chromaticindex32-8` the filled graph still yields **no cuts at
all**. A set partitioning row is satisfied exactly by the relaxation, `sum x_j = 1`, so
a clique drawn inside one row is never violated. Only a clique spanning several rows can
cut, and those are rare. The empty graph was never why those instances fail, which
retires the theory that filling it would close them.

What it did buy, on the pure binary set at one thread: `eil33-2` from not finishing in
120 seconds to 92, the root bound on `ab71-20-100` from unmoved by its 25 cuts to 1.5%
tighter with 339, and `air05` from no cuts to a root bound of 25877.6 -> 25957.1. Across
the 33 instances some open source solver closes and this one does not, it converts none
of them. Across the 23 it already closes, nothing regresses.

One thing tried alongside it and reverted. Rounds that add cuts without moving the bound
look like pure waste, and on `acc-tight4`, whose objective is constant so that no cut can
ever tighten it, they are: 255 cuts and three times the root time to move a bound from
zero to zero. Stopping after three flat rounds fixes that and costs more elsewhere,
`decomp1` from 20.9 seconds to 31.9 and `decomp2` from 32.4 to 58.5. The bound does not
climb steadily, it jumps, and a run of flat rounds is not evidence the next one is flat.

A caution for anything built on this graph next. The extraction takes each row's longest
overshooting prefix and stops, which is exactly right for set packing and partitioning
rows and leaves cliques behind on general knapsack rows. It also refuses any row holding
a non-binary column outright, where the binaries in such a row can still conflict.

### Where the binary losses actually are

Two numbers reframe the problem. Of the 140 pure binary instances, this solver closes 24
against 46 for the strongest open source solver. But of the 116 it misses, **83 are
missed by every open source solver too**. The set that separates this solver from them
is 33 instances, not 116, and passing the best of them means converting 23 of those 33.

Measured across those 116 at 30 seconds and 16 threads, 86% finish with **no incumbent
at all**, and 85% of the addressable 33 do. Read alone that says feasibility is the
blocker and heuristics are the answer. Read with the node counts it says something else,
because the instances with no incumbent divide into two groups that need opposite work.

Twenty of the addressable 33 reach one to three nodes, and timing the root alone shows
why: every one of them spends the entire 30 second budget there and never gets out.

```text
air04         30.1s root,  38313 simplex iterations, 1 node
supportcase4  30.0s root,  20166 simplex iterations, 1 node
tanglegram6   30.1s root,   1324 simplex iterations, 1 node
bley_xl1      31.4s root,   1863 simplex iterations, 1 node
```

No heuristic can help these. A diving or pumping heuristic needs LP solves, and there is
no budget left to give it. The blocker is root throughput, the relaxation and the
separation on top of it, and it is the largest group by some way.

Eight more search thousands of nodes and still find nothing: `acc-tight2` at 43 nodes,
`air05` at 289, `mine-166-5` at 317, `neos-820879` at 3408, `disctom` at 3141, `mitre`
at 6371, `neos-3045796-mogo` at 6752, `graph20-20-1rand` at 133. These are the genuine
feasibility failures, the ones where rounding, diving and the pump all run and all come
back empty, and they are where fix and propagate would earn its place.

The remaining five find an incumbent and lose on the bound: `chromaticindex32-8`,
`n2seq36f`, `neos-1516309`, `neos-1599274`, `neos18`.

So the ordering that suggested itself from the incumbent count alone, heuristics first,
is wrong. Twenty instances need the root to finish before anything above it can run at
all, eight need feasibility, five need the bound. That ordering is what the numbers
support, and it is worth re-deriving before the next build rather than inheriting.

### What the root actually spends its time on

The twenty instances that reach one to three nodes were assumed to be losing to the
relaxation, and the assumption is wrong for most of them. Timing the phases apart, with
`relax` for the relaxation alone and `--cut-rounds 0` for the root without separation,
splits them cleanly.

Seven are genuinely relaxation bound: `air04`, `bley_xl1`, `cod105`, `neos-1324574`,
`neos-3226448-wkra`, `supportcase4` and `tanglegram6` do not finish their relaxation in
120 seconds. Nothing above the LP reaches them, which is the same wall `ex9` and `ex10`
sit behind.

**Read this with "The rescue that spent the whole run" below.** It was measured with
`ripsolve relax`, which solves the relaxation and stops. The solver perturbs and tries
again when the first attempt fails, and did so already when this was written, so what the
seven are actually blocked on is that second attempt rather than the first.

The other thirteen solve their relaxation quickly and then spend the root on heuristics
that find nothing:

```text
                relax it   pre-cut root it   heuristics   incumbents
neos-787933         2163             26179          92%            0
mine-166-5          2278             16914          87%            0
ab71-20-100         4366             30457          86%            0
neos-957143         5366             29688          82%            0
mitre               1825              7596          76%            0
```

Across the group 71% of pre-cut root iterations go to the chain of rounding, diving and
pumping, and on all thirteen it returns nothing. Each solve inside the chain is bounded
and the sequence is not, so on a model whose relaxation is cheap it runs a very long way.

Two separate hypotheses died on the way here, both from reading rather than measuring.
Gomory separation looked like the cost, being a BTRAN and a pass over every column per
candidate row, on the order of 1e11 operations a round for `bley_xl1`. Timing each
separator showed the cut loop never runs at all on those models: it breaks the moment
the root LP is not optimal. The relaxation then looked like the cost, and it is, for
seven of twenty and not the rest.

Bounding the chain was tried twice and reverted both times. As a share of the time limit
it does not bind where it needs to: the models whose chain runs longest are the ones
whose relaxation is cheap, so at a 60 second limit a share of 0.43 left `ab71-20-100`
byte identical at 30457 iterations. Measured against the relaxation's own cost instead,
at twice what the relaxation took, it binds the same way at every time limit, 30457
iterations down to 13700. It converts nothing: no instance newly solved of the 33, no
incumbent gained or lost, `cap6000` 1.44 times slower, and one unreproducible run
reporting `supportcase14` infeasible where it is optimal. A change that buys nothing and
carries an anomaly nobody can explain is not worth its constant.

What the split does say is that the thirteen have root time to reclaim and the seven do
not, and that reclaiming it is worth nothing until something exists that can use it.
Heuristics that find nothing and cuts that move no bound are the same problem seen twice.

### Fixing with propagation, for the models diving cannot reach

Thirteen of the instances this solver loses on solve their relaxation quickly and then
spend the root on heuristics that find nothing, 71% of the pre-cut root iterations and no
incumbent on any of them. What they have in common is a sparse feasible set expressed as
set partitioning rows, and what diving does with those is pay an LP solve to be told
about one more fractional vertex.

The alternative is to stop asking the relaxation. A fixing propagates two ways to a
common fixed point. Through the conflict graph, a literal that holds excludes every
literal it conflicts with, and excluding `x_k = 1` is fixing `x_k = 0`; this is where a
partitioning row pays, since fixing one of its columns to one settles every other column
in the row at once without an LP and without looking at the row again. Through row
activities, a row whose remaining slack cannot absorb a coefficient forces that column
to the end that fits, which catches the rows that exclude no pair but still leave one
value open once enough of their columns are pinned.

This is what the long-row clique extraction was actually for. Reading conflicts out of
partitioning rows bought two root bounds and no solves as a cut source; as a propagator
it is the whole mechanism.

The relaxation still picks the order, most nearly decided first, and the value to try.
What it does not do is get asked again after every fixing.

Safety here is structural rather than argued. Propagation may be wrong about where the
feasible points are, and cannot be wrong about what it returns, because the completed
assignment is checked against the model before it is handed back. A heuristic that
returns an infeasible point does not fail loudly, it installs a bogus incumbent and the
search prunes the real optimum away.

Measured on the pure binary set, four instances that had never found a feasible point now
find one: `ab71-20-100`, `neos-787933`, `neos-953928` and `neos-957143`, each confirmed
at one thread where the result is reproducible. None of them closes. Nothing regresses
among the 23 already closed.

It runs last in the chain despite being the cheapest link in it. Ordered first, ahead of
diving, it took `cap6000` from 500 nodes to 1500: it found a point, short circuited the
chain, and left the search a weaker incumbent to prune with than diving would have given
it. Cheapest first is the wrong rule for a heuristic that is cheap and whose points are
poor. Where the chain above it works there is nothing here worth having, and where that
chain fails this is the whole of what there is.

A caution on reading its results. The two instances that appeared to lose an incumbent
when this was measured at 16 threads, `mine-166-5` and `neos-820879`, find none in either
version at one thread. Incumbent counts on the parallel search are not stable enough to
read a regression from.

### A feasibility search that does not use the relaxation

Every heuristic here asks the relaxation where to look. On a model whose feasible set is
sparse the relaxation does not know, and on seven of the instances this solver loses on
the relaxation does not finish at all, so nothing built on it ever runs.

Feasibility Jump, from Luteberget and Sandvik, asks the constraints instead. Every row
carries a weight, the thing being minimized is the weighted sum of how far the rows are
from satisfied, and each step flips whichever column reduces that sum the most. Arriving
somewhere no single flip improves is not the end: the weights of the rows still violated
go up, reshaping the surface until a flip helps again. What is here is its feasibility
half, with no objective term, because the measured blocker is finding any point at all.

It closes three instances outright, `acc-tight2`, `disctom` and `neos-913984`, each in a
single node, and gives `cod105` the first feasible point it has ever had here. The
reason those three close is worth keeping: the point it returns *is* the optimum, and the
root bound already matches it, so the search only has to agree.

That is also what decides when to keep its answer. Where this heuristic wins it wins
outright, at a gap of zero; where its point is poor it is worse than none, installed
early enough to steer the search and too loose to prune with. `eil33-2` at 37% off the
bound went from solved in 96 seconds to unsolved in 150. So the point is kept only when
the relaxation agrees it is a good one.

Where it runs matters more than anything inside it, and this took four attempts to get
right. Run before the relaxation unconditionally, which is what a reference solver's log
appears to show, it is paid for on every model including the great majority that never
need it: `f2gap401600` went from 0.27 seconds to 11.5 and `mod010` from 0.78 to 12.9,
both of which find a point by themselves in under a second. It is now asked only in the
two places nothing else reaches, when the relaxation fails and when it succeeds and every
heuristic built on it comes back empty, and it costs those models nothing.

Four bugs, each found by measurement rather than by reading:

- A move budget does not bound time. A flip costs the columns of the rows it touches, so
  `mitre` spent 34 seconds against a limit of 20 and never reached the relaxation.
- Its deadline was anchored to the start of the run. Once it moved later in the pipeline
  that budget was already spent, and every instance it wins came back empty, given no
  time rather than too little.
- The improvement search re-enters the solver and inherited the jump budget, re-running
  the whole feasibility search on every neighbourhood: 51 seconds of a 96 second solve.
- Its candidate queue superseded entries rather than removing them, so each step popped
  through the accumulated staleness: 6.9 seconds on a model of 3200 nonzeros.

A stall cutoff was tried and abandoned. Cutting off a run that has stopped reducing
violation sounds right and is not: `acc-tight2` sits at a violation of 19 for thousands
of flips before escaping to zero, which is exactly what the weights are for, and any
cutoff tight enough to catch a hopeless run also catches that one.

### Searching the neighbourhood at the root, which did not work

The improvement search is reached once every few hundred nodes, and every instance this
solver loses sits at one or two nodes when its time runs out, so it has never run on any
of them. Moving it to the root, immediately after the first point, follows directly and
is what a reference solver does: on `mitre` its first point is 140535 and the one its
sub-search returns half a second later is 115155, which is optimal.

Here it found nothing, on any instance, and cost `cap6000` 1.83 times its solve. Worse
before it was guarded: asked to beat a point the bound already matched, it could not
succeed and would not stop trying, and four rounds of that spent the budget `disctom` and
`neos-913984` needed to prove the optimum they already held. Guarding it against that
left a heuristic that costs and does not pay, so it is gone.

The gap it was aimed at is real and remains. What the reference solver does at that point
is not what this did: its neighbourhoods come from several constructions, and the one
here fixes columns where the incumbent agrees with the relaxation, which is undefined on
precisely the instances that most need it, the ones whose relaxation never finishes.

### Restarting a stalled relaxation from a perturbed one

Seven of the instances this solver loses never finish their root relaxation, and nothing
built above the LP can reach them: three separate neighbourhood searches were written
for that group over one sitting and every one of them failed for the same reason, that
it asks the relaxation a question and there is no relaxation to ask.

Those models are not slow so much as stuck. Their matrices are all ones, set
partitioning and set covering, where a crowd of bases describes one point and every step
between them has length zero. `ex9` spends 8871 consecutive iterations taking steps of
length zero without its worst violation moving off 1.0.

Moving every bound by a random amount too small to matter breaks the tie. No two
variables sit on the same bound any more, so the steps have somewhere to go.

The part that took two attempts to understand is what the perturbed solve is *for*. Its
objective is worthless: it is a valid bound on a weaker problem and can be far from the
real one, `neos-1324574` reporting -0.0001 where the truth is 4.5, and `tanglegram6`
reporting a number that scales linearly with the size of the perturbation. What is worth
having is the **basis**. A basis optimal for a problem next to this one is a warm start
this one can use, and re-solving the true model from there recovers the true optimum
almost immediately.

```text
                    plain                     perturbed, then cleaned
neos-1324574    212471 iterations, 214.2s     419 + 4244 iterations, 5.2s
tanglegram6     does not finish               10754 iterations, 172.3s
```

Both cleaned solves return the true optimum, 4.5 and 0.0, proved on the true bounds.

An earlier attempt at this was measured, reverted, and recorded as a failure. It differs
in one respect: it perturbed and warm started from the **stalled** basis, which puts the
search back inside the degeneracy it was trying to leave. The note left behind at the
time said the difference was unexplained and might be the whole of it. It was the whole
of it.

In the solver it is asked for only when the relaxation has failed on its first attempt,
so models that finish theirs pay nothing, and the regression set is unchanged to within
noise across all 23. It converts no instance at a minute. What it does is take
`neos-1324574` from one node and an incumbent of 36 to nine nodes and an incumbent of 14,
which is the difference between a search that cannot start and one that can.

### The rescue that spent the whole run, and the correction it forces

"Seven are genuinely relaxation bound" above was measured with `ripsolve relax`, which
solves the relaxation and nothing else. The solver does not do that. It solves the
relaxation, and when that fails it perturbs the bounds and tries again, and the second
attempt was written after the measurement that the seven come from. So the claim has been
stale ever since, and reading it as "these seven never get a root" is reading a
measurement of a code path the solver does not take.

What it actually takes, traced on all eight models that need a rescue here:

```text
air04              first attempt IterationLimit at 24s, perturbed rescue ends at 60.0s
bley_xl1           ...                                  ends at 60.2s
cod105                                                  ends at 61.0s
supportcase4                                            ends at 60.0s
neos-3226448-wkra                                       ends at 60.0s
ex9                                                     ends at 60.0s
ex10                                                    ends at 60.0s
tanglegram6                                             ends at 60.0s
```

The rescue is reached on every one of them and is handed the caller's entire remaining
clock, because it is given the run's deadline and nothing else. Every one of them spends
it and returns IterationLimit. This is the first lesson of the handoff made one level
further down: not a share of the time limit this time, but no budget at all, which is the
same mistake with the share set to one.

The perturbed solve is the same model at the same size as the attempt it is repeating, so
what that attempt spent is what it is worth spending again, and it is now bounded by
exactly that. Nothing that works is anywhere near the bound: `neos-1324574` needs five
iterations after the perturbation and `tanglegram6` 172, against first attempts of tens of
thousands.

### A second rescue, and why the dual method is safe here and not in general

Bounding the first rescue makes room for a second. Phase 1 is where the primal method gets
stuck, and the reason is structural rather than tunable: the phase-1 cost vector scores a
basic variable by whether it currently *violates* a bound, so it is blind to the kink at
one sitting exactly on one, and no ratio test — short, long or Harris — repairs a column
choice made that way. The dual method has no phase 1 to be stuck in.

That is why `dual-cold-start` was built, and it is not why it was shelved. It was shelved
because the two methods end on *different* optimal vertices and everything downstream
reads the root's vertex rather than its objective: Gomory cuts come off that tableau and
branching reads its fractional values. Made the default, it reaches the worse vertex on
most of this set and widens most gaps even where the relaxation improves sharply.

Reached only where the primal method produced no vertex at all, there is nothing for its
vertex to be worse than. That is the whole of the safety argument, and it is why the entry
is gated rather than switched: warm starts, which are every node of the search, keep the
row selection they had, and only a cold entry pays the extra FTRAN per iteration that
steepest edge row pricing costs.

```text
                 primal          dual cold
air04            >130s           1.5s,  3644 iterations
tanglegram6      >130s           0.5s,   207 iterations
bley_xl1         >130s           >130s
cod105           >130s           >130s
supportcase4     >130s           >130s
ex9              >130s           >130s
ex10             >130s           >130s
neos-3226448     >130s           >130s
```

Two of eight. In the search `air04` goes from one node to 380 and `tanglegram6` from no
incumbent to one.

Three details had to come across with it, and each of them was a defect somewhere:

- Steepest edge row selection is most of why it works, and the chosen row's weight has to
  be read exactly from the pivot row already in hand, `rho . rho`, rather than carried
  forward through the recurrence: carried, `drayage-100-23` takes 487230 iterations where
  the exact value takes 2374.
- A violating row must never be passed over. Scores underflow to zero when a weight is
  enormous, and a rule accepting only a strict improvement over zero then selects nothing,
  which is read as primal feasibility and ends the solve at a basis violating by 1.6e5.
- Primal feasible is the optimum only if the basis is still dual feasible, which the ratio
  test holds in exact arithmetic and drifts out of in this one. Checked rather than
  assumed, and where it has lapsed the basis goes to the primal loop instead of being
  ruled on.

A cold entry also needs its own way out of a degenerate patch. Warm starts here hand the
node to the primal loop when the dual objective freezes, which is cheaper than Bland's
rule and right for them. For a cold entry the primal loop is the thing it was reached in
place of, and handing back means handing back to something stuck: `tanglegram6` finishes
in 207 pivots when the dual method escalates to Bland's rule instead, and does not finish
at all when it gives up.

### What this converts, which is nothing, and what that says

No instance is newly closed. The reason is visible in the trace and is the same lesson a
third time: at a sixty second limit the first attempt is given 40% of the run before
anything else may be tried, so `air04`'s root is answered at twenty-five seconds and the
search gets what is left. `ROOT_LP_FIRST_SHARE` is a share of the caller's time limit, and
the models it hurts are exactly the ones where the primal method is not slow but stuck,
which is a property of the model and knowable without the clock.

So the next question is not another root method. It is why a method that is going to fail
is given nearly half the run to fail in, when whether it will fail is decidable from
whether phase 1 is making progress. The shape that suggests itself is the one that worked
for probing: a small budget first, escalating only where the small one paid, with the
primal method warm started across the escalation so it pays nothing for being interrupted.

One thing was gained on the way that is worth keeping either way. Both rescues now cost
what they are worth rather than what is left, and a model whose rescues are cheap and
hopeless no longer returns at 41 seconds of a 60 second budget with a third of the run
unspent.

### A test that could not fail, again

The dual entry is checked by comparing its relaxation values against a reference solver's,
which is the only reason the branch's two wrong answers were ever found. The obvious way
to write that check is to skip whatever comes back non-optimal, since a model that cannot
be parked dual feasibly declines to run at all. That version passes with the dual ratio
test deliberately inverted. A broken dual method does not return a wrong optimum, it stops
reaching an optimum, so every model gets skipped and the test reports success — slowly,
which is the only visible symptom.

Every sample is now required to be either answered or declined, and a declined one to have
spent no iterations at all, which is what distinguishes "cannot be parked" from "tried and
failed". With that, the inverted ratio test fails on the first sample.

### The addressable set, re-derived

The breakdown in "Where the binary losses actually are" is from an older build and has
been quietly wrong for a while: `mitre` has since moved out of it, and the seven called
relaxation bound were measured through a code path the solver does not take. This is the
current one, taken by running all 31 at 60 seconds and 16 threads and reading nodes, cuts
and whether any feasible point was found. The 31 are the instances some open source
solver closes and this one does not.

```text
instance                   nodes   cuts  point  gap
bley_xl1                       1      0  no     
cod105                         1      0  yes    -
ex10                           1      0  no     
ex9                            1      0  no     
neos-3226448-wkra              1      0  no     
neos-780889                    1      0  no     
supportcase4                   1      0  no     
neos-1324574                   2      0  yes    67.8571
neos-5129192-manaia            2      0  no     
neos-633273                    2    109  no     
neos-953928                    2      0  yes    104.3261
neos-957143                    2    128  yes    68.9214
neos-960392                    2     27  yes    238000000000000.6250
tanglegram6                    2      0  yes    100.0000
acc-tight4                    12    128  no     
acc-tight5                    52    192  no     
neos-787933                  117    192  yes    89.1029
neos18                       134    190  yes    60.1449
ab71-20-100                  135    246  yes    579.2901
graph20-20-1rand             290      0  no     
air04                        364      0  no     
mine-166-5                   877      0  yes    35.5003
air05                       1064     48  no     
nw04                        1574     80  yes    0.3992
chromaticindex32-8          7653      0  yes    25.0000
neos-3045796-mogo           7819   2688  no     
neos-820879                 8720      4  no     
irp                        16058    340  yes    -
neos-1599274               33419    432  yes    3.8533
n2seq36f                  158581    978  yes    0.3831
neos-1516309              468406      7  yes    1.4254
```

Three groups, and they want opposite work.

**Fourteen stop at one or two nodes.** Eleven of them add no cuts at all, which is the
tell that the root relaxation never finished: the cut loop breaks the moment the root LP
is not optimal. Those eleven are the group the two root rescues above are aimed at, and
the rescues answer two of them. Three more — `neos-633273`, `neos-957143`,
`neos-960392` — do finish their relaxation and separate cuts, and then stop anyway, so
they belong with the group below rather than here.

**Seven search and find nothing:** `acc-tight4`, `acc-tight5`, `graph20-20-1rand`,
`air04`, `air05`, `neos-3045796-mogo`, `neos-820879`, from twelve nodes to 8720. These
are the genuine feasibility failures, where rounding, diving, the pump, fixing with
propagation and the LP-free jump all run and all come back empty.

**Ten search, hold a point, and lose on the bound.** Four of them are close enough that
the bound is the only thing between them and a proof:

```text
n2seq36f          0.38%   158581 nodes    HiGHS closes it in 3.6s
nw04              0.40%     1574 nodes    CBC closes it in 15.4s
neos-1516309      1.43%   468406 nodes    HiGHS closes it in 0.3s
neos-1599274      3.85%    33419 nodes    HiGHS closes it in 0.6s
```

`irp` is a fifth, and it is not really in this set: it closes about half the time at 57
to 60 seconds, so the headline figure reads 26 or 27 depending on it. Anything measured
against that figure has to say which.

That is where the next work belongs, and it is not where the last three sessions have
been looking. The 31 divide 14 / 7 / 10, and the ten holding a point and losing on the
bound are the only group with instances within a percent of closing. Reduced cost fixing
is the obvious untried thing for them and does not exist here: with an incumbent and the
root's duals in hand, a nonbasic binary whose reduced cost exceeds the remaining gap
cannot move in any better solution and can be fixed for good. At a gap of 0.4% that is
most of the model. Nothing in this solver reads a reduced cost outside the simplex.

### Scaling, which does not answer this

Row equilibration was written for the same group and does not fit it. Five of the seven
have matrices of all ones: `air04`, `cod105`, `neos-1324574`, `neos-3226448-wkra`,
`tanglegram6`, along with `ex9`. There is nothing to scale, and `air04` bears that out,
95784 iterations before and 94302 after.

Where coefficients do vary it buys throughput and no solves: `supportcase4` 13% more LP
iterations in the same budget, `mitre` 6%, and `bley_xl1` unchanged and still unsolved
despite spanning 0.5 to 1e8. Held back rather than kept, since what it treats is
conditioning and what stops these models is degeneracy.

### Probing, which reduces the models and cannot be afforded

Presolve here fixes columns by reasoning forwards from the bounds it has. Probing
reasons backwards: suppose a column takes a value, work out what follows, and if that
ends in a contradiction the value is refuted and the column is fixed to the other one.
Nothing is guessed, and a value refuted this way appears in no feasible solution.

Reasoning it with the search's own propagator, conflict graph and row activities
together rather than activities alone, is what makes it strong here. The conflicts of
these models live in long set partitioning rows that pairwise reasoning never sees, and
`ex9` carries 1.9 million of them.

```text
                 presolve alone    with probing     a reference solver
ex9              1844 columns      6907 of 10404    all of them
mitre               0              3677 of 10724    4693
air04               0              1375 of  8904    1400
```

On `air04` that is the reference solver's own reduction, matched. It also solves `mitre`,
which nothing else here has.

It is reverted anyway, because a probe costs a propagation sweep per column per value
and no budget was found that keeps the reduction without paying for it somewhere worse:

- Unbounded, presolve runs past the caller's entire time limit. `air03` solves in 1.2
  seconds and was still probing at 200; five instances of the regression set were lost
  outright.
- Bounded by a share of the time limit, the reduction collapses, `ex9` to 1920 columns
  and `air04` to 114, while a 0.14 second instance still pays six seconds. A share of
  the limit prices work against how long the caller happens to be willing to wait, which
  is unrelated to what the work is worth. This is the third time that mistake appears in
  these notes.
- Bounded by a run of unproductive probes, which is the right shape of budget, three
  instances come back to par and `air03` goes to 27 seconds instead, with the reduction
  down again.

What is missing is a way to tell, before probing, which models it will pay on. Ordering
the candidates was tried next and found the answer without being able to use it.

The cost of a probe is set by how far its consequences reach, and reach is the size of
the cliques the column sits in. That separates the two groups cleanly, and it is a
property of the model rather than of the clock:

```text
   pays                        ruins
   mitre       28 literals     air03     3861
   ex9        196              eil33-2   2568
   air04      368              irp       6146
   decomp2     16              nw04     42032
```

Which makes "most constrained first" exactly backwards. The widest columns are where the
refutations are and also the ones that cannot be afforded; ordering by reach and taking
the widest first probes the most expensive columns in the model before any others.

Skipping the columns above a reach of a thousand and taking the rest widest first
restores most of the reduction, `mitre` to its full 3677 and `air04` to 787, and still
costs `air03` 3.7 times its solve, because the cap gates the column probing starts from
and not the ones its propagation walks into. Counting each literal visited as a unit of
work rather than each column fixed, which is the honest way to bound it, moves the cost
somewhere else again: `mod010` from 0.63 seconds to 11.

Six budgets were tried. Every one of them trades the reduction against the cost, and the
coupling is not incidental: probing pays on models whose conflicts are dense, and dense
conflicts are what make it expensive. Something has to break that coupling before this
is worth having, and a budget is not it.

One thing found on the way is worth keeping in mind whatever happens to probing: reading
a literal's exclusions through `neighbours` alone reports zero for every column of a set
partitioning model, whose conflicts are all held as cliques. Any future caller wanting a
column's true degree has to count clique membership too, or it will conclude that the
most constrained columns in the set are the least.

The measurement stands whatever happens to the code: the reduction is available, it
matches a reference solver on one instance and gets most of the way on two more, and it
is the only thing tried here that closes `mitre`.

**Superseded.** Everything above about the cost is true of the implementation and false
of the idea; the next three sections are what happened when the sweep stopped rebuilding
itself. The reasoning about reach, and about what "most constrained first" means on a set
partitioning model, still holds.

### The corners of the box

The cheapest heuristic here, and the last one added. Putting every column on one of its
bounds costs a single pass over the matrix, and is feasible more often than it deserves
to be: on the pure binary set the all-zero point satisfies `mine-166-5`, `neos-953928`,
`neos-957143` and `neos-960392`, and the all-upper point satisfies `neos-787933`. A
reference solver reports one of these as its first incumbent on eight of the thirty
instances this solver loses to it.

`mine-166-5` and `neos-957143` went from reporting no solution at all to reporting one.
Neither closes, and the points are poor, which is the entire reason this runs last rather
than first: a corner of the box is exactly the kind of cheap, weak point that takes the
search's incumbent away from a better one. Put ahead of diving it would repeat what
putting the propagating heuristic ahead of diving did to `cap6000`, 500 nodes to 1500.
Reached only when everything above has failed, it is the difference between a feasible
answer and none, and nothing is compared against it because there is nothing to compare.

That ordering mistake has now been made three times in this file, with two of them
already written down when the third was made. Cheapest first is the wrong rule whenever
a heuristic is cheap because its answers are bad.

### What probing actually cost, which was not probing

This supersedes "Probing, which reduces the models and cannot be afforded" above.

Probing was parked twice for cost, both times after trying another budget. The cost was
never in the idea. It was in what a sweep rebuilt before it did any reasoning, and a
probe does two sweeps per column.

Three things were being rebuilt. `adjacent` allocated a vector, cloned a literal's edges
into it, appended each of its cliques, sorted the result and deduplicated it, once per
column fixed; the widest literal of `air03` excludes 3861 others, so that is a
3861-element sort performed to spare the propagator a comparison it does not mind making.
A trial assignment was a copy of the model's bounds, so a probe paid a pass over the
columns before looking at anything, and a model of ten thousand columns was copied twenty
thousand times over a pass. And a row's activity range was recomputed from its
coefficients whenever any column in it moved, which on a set partitioning model is most
of the matrix per fixing.

None of it is needed. Literals are visited in place. The assignment records the value it
overwrote and undoes by walking that record backwards, which is what "does not rebuild
its state per probe" meant when this file asked for it. Row activities are carried on the
same record, restored to the value they held rather than by re-adding what was subtracted,
because an activity walked forwards and backwards through a million floating point
additions does not come back to where it started and probing fixes columns for good on
what these numbers say. And a row whose slack at both ends exceeds what any single column
could contribute to it can force nothing whatever is fixed in it, so it is skipped on two
comparisons against numbers it already carries, instead of being read.

Measured on presolve alone, with the search taken out of the picture:

```text
                 before              after
mitre            3677 in 0.29s       3677 in 0.06s
mod010            168 in 8.97s        168 in 0.10s
air04             787 in 10.19s      1343 in 2.33s
air03               0 in 3.14s          0 in 0.51s
```

`air04` is the interesting column. The old budget was spending ten seconds to get 787
columns because most of the ten seconds was overhead; the same budget spent honestly gets
1343, against a reference solver's 1400.

### A budget for probing that is not another share of something

With the overhead gone the budgets could be reconsidered, and the six that failed have
one thing in common: they counted probes. A probe's cost is how far its consequences
reach, and reach spans three orders of magnitude within one model, so counting probes
prices a column in a clique of four the same as one in a clique of four thousand. A
model whose probes are all expensive and all fruitless then spends its whole pass earning
the right to stop.

The unit is matrix entries read. Two bounds on it, because they catch different failures:

- **Patience.** What may be read without proving anything, reset on every column proved.
  A model that keeps proving things keeps its budget; one that has stopped proving them
  loses it quickly. This is what a run of unproductive probes was trying to be.
- **Total.** What may be read in all. Patience alone does not bound a model with nine
  thousand proofs to find, each individually earned and affordable, because every one of
  them resets the count: `ex10` goes on finding them for two minutes.

Both are absolute counts rather than multiples of the matrix. Multiples of the matrix
were tried and are the wrong shape: a model with eight million nonzeros then gets a
hundred times the work a model with eighty thousand gets for the same reduction, and
`rail02` and `opm2-z12-s8` spend twenty seconds proving nothing on that account. What the
absolute figure encodes is a judgement about how much reasoning about a model is worth
doing before solving it, which is a property of neither the caller's clock nor the
matrix's size.

A third bound, on one probe, is separate from both rather than a share of them. It says
how far a single consequence is worth chasing; shrinking it because the model has already
proved everything it is going to prove would abandon every remaining probe at once, which
is what happened when it was derived from patience.

The distinction that makes patience work is between a sweep that finished and found
nothing and a sweep that ran out of work. The first says the column has nothing to give.
The second says only that it was not allowed to look. Treating them alike stops the pass
on models that were about to pay: with the two conflated, `ex9` stopped after 512 probes
having proved nothing, because all 512 had been abandoned rather than completed.

### What the reduction is worth, and what it is not

At the figures chosen, every instance this solver already closes pays under half a second
for probing and most pay under a tenth, `mitre` gets the whole of its 3677 column
reduction, and `mitre` closes at 18.2 seconds. Nothing else tried here has closed it.

Ten times the budget buys several thousand more fixings on `ex9`, `ex10`, `air04`,
`neos-4754521-awarau` and `rail01`, and closes none of them. Those first four all stop at
one node with the root relaxation unfinished, which is the wall already recorded for them
and is not a presolve problem: `air04` fixed 1343 columns and still spent sixty seconds
on 37994 simplex iterations without leaving the root. The same ten times costs `irp`,
which closes with two seconds to spare, more than two seconds.

So the honest reading is that the reduction is now affordable and only one instance was
waiting for it. The rest of the twelve that lose a fifth of their columns are waiting for
something else, and it is the same something the twenty root bound instances are waiting
for.

### Two corrections found by rewriting it

Row propagation offered `0` and `1` to any column in the row, taking whichever was the
only one that fitted. On a continuous column bounded in `[0, 5]` that "forces" it to 0
whenever 1 does not fit, which is harmless in the heuristic, where the point it produces
is checked for feasibility afterwards, and unsound in presolve, where it is a proof that
fixes a column for good. Values are now offered only to binaries.

The clock backstop was read at every 64th candidate but placed after the skip for
candidates with no conflicts, so a model whose first candidates were all skipped ran up
to 63 probes past its deadline before noticing.

### Fixing the columns the bound and the incumbent already decide

Every reduction before this one reasons about feasibility: presolve asks what the rows
allow, probing asks what a value implies. This one reasons about *optimality*, and it is
the first thing here that does.

At an optimal basis a nonbasic column sits on a bound whose reduced cost `d` says the
objective cannot fall by leaving it. Moving a distance `t` off that bound therefore raises
the objective by at least `|d| t`. The root's value is a bound on every solution in the
tree, so any solution better than the incumbent `u` satisfies `root + |d| t < u`, which
caps `t` at `(u - root) / |d|`. For an integer column the cap rounds inwards, and on a
binary a cap below one decides the column outright.

Nothing here read a reduced cost outside the simplex before this. `LpSolution` carried the
status, the objective, the primal values and the basis, and no duals at all.

```text
                fixed at the root
nw04            25471 of 87482
neos-1516309      350 of 4500
n2seq36f           41 of 8100
irp                 3 of 20315
```

`nw04` closes, in 32 to 38 seconds, where it had sat at a 0.40% gap for the whole minute
and had never closed in any run of any version. Five consecutive runs give four optimal
and one timeout, which is the honest form of that claim; a single benchmark run caught the
one and reported no change at all. `neos-1599274` goes from 3.85% to 2.62%. `irp` closes
in 46 to 56 seconds where it used to take 57 to 60, and is still a coin flip, two runs in
five.

Which is worth stating separately from the change, because it is a fact about the
measurement rather than about the solver: **two of the twenty-eight instances this solver
can close sit on the sixty second line**, so the headline count reads 26, 27 or 28
depending on where two flips land, and a single run of the benchmark cannot distinguish a
real gain of one from noise of one. The number to quote is the set that closes every time,
which went from 25 to 26 this session, with `nw04`'s four-in-five stated beside it.

The narrowed bounds go into the model rather than into the root node, because a node
rebuilds its bounds from the model before every solve. One pass reaches the whole tree,
and it costs one copy of the model.

**What it does not reach, and why that is most of it.** The bound this uses is the root's,
and the incumbent it uses is whatever the root heuristics found, which is far weaker than
what the search eventually holds: `n2seq36f` is at a 39.7% gap at the root and 0.38% at
the end. The room `(u - root)` is therefore near its widest exactly when the fixing runs,
and near its narrowest when it would pay most. Re-running it as the incumbent improves is
where the rest of this is, and the awkward part is that the search is parallel and the
model is shared immutably across its threads.

### A margin on the wrong side is not a small mistake

The room a better solution has to move in is `u - root`, and it is tempting to subtract a
tolerance from it "to be safe". That is backwards, and the direction is worth stating
plainly because both readings feel conservative:

- A **smaller** room caps the travel harder, so **more** columns are fixed. That is the
  direction that can remove a solution nobody has seen yet.
- A **larger** room fixes fewer columns and can only lose reduction.

A first version subtracted the feasibility tolerance and fixed 1566 columns of `irp` where
the sound version fixes 3. It passed the tests. Every widening here — a relative slack and
an absolute one on the travel cap, the integrality tolerance added before the floor — is
on the generous side for the same reason.

### The test, and what it does and does not catch

Solving every sample twice, with the fixing on and off, and requiring the same answer.
Checked against the search with the claim disabled rather than against a recorded value,
because a recorded value would also have to be trusted.

It was then checked for whether it can fail, which is the practice this file keeps
recommending and this session kept needing:

- Cap set to zero, pinning every nonbasic column where it sits: **caught**, on
  `v064c200`, which answers 1039 against 225. But only once the generated families were
  added; the bundled samples alone do not catch even that.
- Cap halved, which is plainly unsound: **not caught** by either set.

So it is a guard against a proof that is wrong in kind, not one that is wrong at the
margin. That is worth having and worth not overstating, and the test says so in a comment
rather than leaving the next reader to assume more of it.

### Reduced cost fixing beyond the root, which is a table and not a mutation

The root fixing spends whatever incumbent the root heuristics happened to find, and on
this set that is far weaker than the one the search ends with: `n2seq36f` is at a 39.7%
gap at the root and 0.38% at the end. So most of what the root's reduced costs prove is
not provable when the only chance to use them goes by.

This was written up as blocked, on the grounds that re-deriving the fixing whenever the
incumbent improved means narrowing a model every worker holds immutably. It is not
blocked, and HiGHS's `HighsRedcostFixing` says why: nothing in the derivation changes
except one number. The reduced costs come from the root basis, the root's value is fixed,
and only the room `u - root` shrinks. The bound each column will eventually take, and the
incumbent at which it takes it, can both be worked out at the root and read off later.

That makes it a table rather than a mutation, which is what makes it free in a parallel
search. Ordered by the incumbent each entry needs, the entries in force are always a
prefix: a worker recomputes the prefix length when the incumbent it sees moves, which is
rarely, and applies it per node for the cost of the bound resets it was already doing. No
coordination, and a worker behind on the incumbent applies fewer entries, which is wrong
only in the direction that fixes less.

```text
n2seq36f     5910 of 6642 lurking bounds in force by the end, 158581 nodes a minute -> 314790
neos-1599274  500 of 1950
nw04         five closes in five at 31 to 33s, against four in five at 32 to 38
irp          five closes in five at 34 to 47s, against two in five at 46 to 56
```

`nw04` and `irp` were both coin flips against the sixty second limit and are neither any
more, which is worth two instances on the standing and is the first time this file has
been able to say a headline figure is not sitting on a boundary.

Nothing that was far from closing closes, and `n2seq36f` shows why in one line: **its bound never leaves 52000,
the root's own value, however many columns are fixed underneath it.** Fixing columns makes
the tree cheaper to walk. It does not make the bound move.

### Three things that do not move that bound, measured rather than assumed

**A restart does not.** The obvious amplifier, and what HiGHS does at 55.8% inactive
columns. Simulated before building it, with `cargo run --example restartsim`, which solves
for one budget, applies the fixing the incumbent earned, and solves what is left for
another. `n2seq36f` narrows to 6348 of 8100 columns decided and the restarted search
returns a bound of **52000**, to the unit. `neos-1516309` comes back *worse*, 35487
against 35525, having thrown away the cuts the first pass accumulated. This is a large
refactor that the measurement says not to do.

**Objective granularity does not, here.** Every coefficient of `n2seq36f` is a multiple of
200, so 52000 and 52200 are adjacent attainable values and the bound is already on one of
them. It is worth having for other reasons and is in; see below.

**Accumulating more cuts does not.** The cut loop ages a cut out after two rounds sitting
slack, and on a degenerate model that looked like the reason nothing accumulates. Holding
every cut instead, 1019 of them, moves the bound by nothing at all. So it is not cut
management:

```text
cuts held    2 rounds     10 rounds    forever
bound        52000        52000        52000
```

### What HiGHS does on the same model, which is the target

```text
 J    0   0.00%   -inf          123000    Large       0 cuts     0 inlp    0.3s
      0   0.00%   52000         123000    57.72%      0          0         0.3s
 C    0   0.00%   52000          63200    17.72%    200         33         0.4s
 L    0   0.00%   52200          57000     8.42%   1120        101         2.0s
 L    0   0.00%   52200          52400     0.38%   1121        108         2.8s
55.8% inactive integer columns, restarting
      ... 254 rows, 1427 cols after restart ...
 L    0   0.00%   52200          52200     0.00%    975        102         3.6s
```

It starts from the same bound this solver has, 52000, and **cuts move it to 52200 at two
seconds**, closing the whole of the gap from the bound side. 1120 cuts generated, 101 held
in the LP. This solver generates 1019 in the same place and moves the bound zero.

So the difference is not the number of cuts, nor how long they are kept. It is what they
cut. `n2seq36f` has 285 rows and 8100 columns, and its optimal face is enormous: every cut
here removes one vertex of it and the LP steps to another vertex worth the same 52000.
Cutting a face rather than a vertex is what single-row separation cannot do, and it is
what this solver's four families all are. HiGHS's extra reach on this model comes from
aggregation, `HighsPathSeparator` and `HighsLpAggregator` building MIR cuts over
*combinations* of rows rather than one row at a time.

That is the next thing to build for the bound group, and it is the first time this file
has had a specific target for it rather than "stronger cuts".

### Bounds lifted to the objective values a model can take

Not aimed at `n2seq36f`, though it is what found it. When every column carrying an
objective coefficient is an integer column and every such coefficient is a whole multiple
of `g`, no feasible point scores anything but a multiple of `g`, so a bound of 52000.3 is
really a bound of 52200. Twenty-four of the thirty-one addressable instances have an
integral objective and four have a spacing above one: `n2seq36f` 200, `neos-780889` and
`bley_xl1` 5, `nw04` 2.

It also tightens reduced cost fixing, which is where most of what it does shows up. A
solution better than the incumbent is better by a whole step, so the cutoff is `u - g` and
the room is a step smaller. On `n2seq36f` the room goes from 200 to nothing, and a room of
nothing pins every column the relaxation left on a bound.

Small, measured at one thread where the search is deterministic: `nw04` optimal in 29.2s
against 30.9, `neos-1516309` a better bound over more nodes, `air05` 127 nodes against 85,
nothing worse. No instance newly closed.

The defect in it is the part worth keeping. The divisor loop returned as soon as the
spacing reached one, which is where it stops improving and so looks like the right place
to stop. It is not: whether the reasoning applies at all is decided by *every* column, and
returning early skips the ones not yet looked at. A fuzz instance whose first two integer
coefficients were 8 and 3 was read as integral before the loop reached the continuous
column with a coefficient of 6.5, and the search then proved an optimum of 34.375 where
the answer is 34.25. The unit tests did not catch it and the differential fuzz did, on the
first run after the change. An early exit on the *value* being computed is safe; an early
exit that also skips a *validity* check is not, and the two are easy to write as one line.

### The cut filter that rejected everything, silently

The root cut loop reports one line: cuts added, bound before, bound after. On `n2seq36f`
that reads "1019 added, 52000 -> 52000", which says the loop is not working and does not
say which part of it. `cargo run --example cutfamilies` runs each family on its own
against the same relaxation, adds what it finds and re-solves, so the families cannot
hide behind one another:

```text
n2seq36f      root 52000
cover      20 cuts, worst violation 0.80, bound 52000
mir        60 cuts, worst violation 0.80, bound 52000
clique      0 cuts,                       bound 52000
gomory      0 cuts,                       bound 52000
```

Gomory reporting zero on a relaxation with twenty fractional columns is not a result, it
is a symptom. It was generating twenty-six cuts and discarding all twenty-six, on a
density filter, and the largest missed by eight terms out of 647.

The filter compares a cut against three times the model's average row. The reasoning
behind "relative rather than absolute" is sound and is written up above: an absolute cap
rejected every GMI cut on a 99.5% dense model, where dense cuts are what there is. What
was not noticed is that the same rule is far too strict at the other end. On a sparse
model a cut denser than the rows is still an ordinary cut, and GMI cuts always come out
denser than the rows they derive from.

`neos-1516309` is the clean case: rows of 62 columns in a model of 4500, an allowance of
186, and 77 of the 80 cuts the root produced thrown away. The three that survived were
worth 1300 of bound between them, and the ones thrown away another 640.

```text
allowance     cuts kept    root bound      (optimum 35954)
186 (3x row)      3        34739
450                10        35187
1125               28        35379
```

### Loosening it needs a bound at both ends

Either bound alone breaks a model that the other one holds:

- **Relative alone.** At any factor loose enough to help, the allowance passes the column
  count on a model with few, wide rows: `irp` has 39 rows and 20315 columns, an average
  row of 2519, and at twelve times that the filter is not a filter. It stops closing
  entirely, nought runs in four, having closed five in five.
- **A share of the columns alone.** On `decomp2`, whose rows average five columns of
  fourteen thousand, a quarter of the columns allows a cut seven hundred times denser
  than anything in the model, which is exactly what the original reasoning warns about.
  It went from closing in 25 seconds to closing once in three attempts.

Twenty-four times the average row, capped at a third of the columns, holds both and
closes two more. `neos-1516309` goes from a 1.4% gap after 450000 nodes to **optimal at
one node in 0.22 seconds**, and `neos-1599274` from 2.6% to optimal in 20 seconds. Three
runs of each of the four instances involved agree, and every instance that closed before
still closes.

The general shape is worth keeping in mind for the next filter: a threshold with one free
end is a threshold that will be wrong on half the models, and the half it is wrong on is
the half nobody is looking at.

### What this does not fix, and why that is the more useful half

`n2seq36f` still does not close, and the diagnostic now says why in a way that four
thousand cuts did not. With the filter out of the way it produces twenty GMI cuts, each
violated by a full unit, and the bound stays at 52000 to the digit:

```text
allowance         cuts    bound
639 (as shipped)     0    52000
810                  1    52000
2025                 4    52000
8100 (all)          20    52000
```

Twenty violated cuts that move nothing is a statement about the shape of the problem. The
optimal face of `n2seq36f` is large enough that each single-row cut removes one vertex of
it and the relaxation steps to another vertex worth the same 52000. Cutting a face rather
than a vertex is what a single-row separator cannot do, and all four families here are
single-row or single-tableau-row. That is the case for aggregation, and it now rests on a
measurement rather than on the observation that HiGHS has one.

### Mod-2 cuts, and a hypothesis corrected by reading rather than reasoning

The section above concluded that `n2seq36f` needs cuts that combine rows, and named
aggregation as the thing to build, on the strength of HiGHS having `HighsPathSeparator`
and `HighsLpAggregator`. That was reasoning from a name. Reading the file says
`HighsPathSeparator` aggregates along **continuous** columns -- every loop in it is over
`continuous_cols` -- and `n2seq36f` has none, so on that model it produces nothing at
all. The right answer was two files away in `HighsModkSeparator`.

The construction is simple and does combine rows. Take tight rows `a_i x <= b_i` with
integer coefficients, add a subset `S`, halve, and round the coefficients down, which is
valid because every column is a nonnegative integer:

```text
    floor(A/2) . x  <=  floor(B/2),      A = sum_S a_i,  B = sum_S b_i
```

The rounding is worth something exactly when it loses nothing on the left and a half on
the right: every `A_j` even and `B` odd. Because the rows were tight, the cut is then
violated by exactly one half.

The step that makes it cheap is that a column sitting *on* a bound costs nothing when its
coefficient is odd, since it contributes zero to the activity either way. So the parity
condition is needed only on the **fractional** columns. `n2seq36f` has twenty-one of them
against 8100 columns, and choosing the subset becomes a linear system over GF(2) with
twenty-one equations and one per tight row -- wide, cheap, and with a large null space
that is where the cuts are. Complementing the binaries sitting at one is what puts every
variable at zero or fractional so that this holds.

A Python prototype against the same relaxation settled it in an afternoon rather than a
build:

```text
round  0: +12 cuts, bound 52000
round  1: + 6 cuts, bound 52200      <- the optimum
...
round 20: no violated cut, stopping
```

Two rounds and eighteen cuts, where four thousand single-row cuts moved nothing. In the
solver `n2seq36f` closes in 8 to 15 seconds, three runs in three.

### An eighth of the allowance, not the whole of it

These cuts are not comparable to the other four families, and giving them the same
allowance is wrong in a way worth naming. Every mod-2 cut is violated by *exactly* one
half, by construction. Ranking a set of identically-violated cuts against families whose
violation varies admits all of them or none of them, so the allowance is doing the
selecting rather than the ranking, and a handful is all they need to be: eighteen closed
`n2seq36f`.

Taken at a quarter of the allowance they crowd the node LPs instead. `irp` is 39 rows
against 20315 columns, so cuts multiply its row count rather than adding to it, and it
went from closing five runs in five to none in three. An eighth holds both.

What does *not* separate the two cases is worth recording, because it was the first
guess: cut density relative to the model's rows. `irp`'s mod-2 cuts run 645 to 2962 terms
against an average row of 2519, so they are *sparser* relative to the model than
`n2seq36f`'s at 936 to 1584 against 213. Density was the right rule for the GMI filter
and is the wrong rule here.

### A skipped node does not always forfeit the claim

A node whose LP runs out of simplex iterations leaves its subtree unexamined, and the
search has always downgraded `Optimal` to `NodeLimit` on account of it. That is right
when the subtree could have held something, and `neos-1599274` is the case where it is
not: with mod-2 cuts it reaches a gap of zero, holds the optimum, and reports NodeLimit
against a node it skipped early on.

The bound that node inherited holds over its whole subtree. If the incumbent has since
overtaken it, the subtree contained nothing worth having, and skipping it cost nothing --
which is the same test the search applies before opening any node at all, so nothing new
is being assumed. Keeping the weakest such bound and checking it at the end turns
`neos-1599274` into four closes in four.

Relaxing a rule about when optimality may be claimed is only safe if it never claims it
wrongly, so the test squeezes the per-node iteration limit until nodes are actually
skipped, and requires every run that says `Optimal` to match the reference optimum.
Without the guard, `v064c064` claims 165 against a true 137.

### Four budgets between a feasibility search and the point it could reach

`neos-3226448-wkra` closes in 8.7 seconds and no simplex iterations at all, having
previously spent sixty seconds on a relaxation. HiGHS answers it in two tenths the same
way. Nothing about the search changed; four separate guards had each stopped meaning what
they were written to mean, and all four had to go before the point was reachable.

**Its objective is empty**, and so is `supportcase4`'s. Every feasible point scores the
same, the bound is a constant that no relaxation is needed to establish, and the first
point found is optimal. Such a model is asking for a feasible point, and the whole
apparatus above the feasibility search exists to find a *good* one. It now goes first and
gets the run rather than a twentieth of it.

**The move budget was an absolute count.** Twenty-five thousand flips is two and a half
per column on a model of ten thousand columns and two hundred and fifty on a model of a
hundred. A flip is a column, so the budget is per column now. What bounds a run going
nowhere is the stall cutoff, and what bounds the whole of it is the deadline; this was
the third guard and should never have been the binding one.

**It never restarted.** The search stops when it has settled, and the weights that got it
there are why it will not move again, so where it begins decides where it ends. From this
model's own bounds it yields nothing however long it is given, and yields a point on the
seventeenth random start.

**Its point would have been discarded anyway.** A jumped point more than a tenth off the
root bound was thrown away, because a poor one had cost `eil33-2` its solve. Re-measured,
that costs nothing now on any instance in the set, and a point is worth more than it was
when the threshold was set: reduced cost fixing reads the incumbent, so a point that
prunes nothing by itself still decides columns. Every point the filter rejected arrived
where the search reports no incumbent at all. `neos-820879` goes from no bound to a gap
of 1.03%, `neos-3045796-mogo` to a bound already equal to its optimum.

### Where the remaining instances actually are, primal against dual

Worth doing before building anything, and it takes one run of a reference solver per
instance. Comparing this solver's bound and incumbent against the true optimum splits the
set cleanly, and the two halves want opposite work:

```text
                    our bound      optimum    waiting on
neos-3045796-mogo       -175          -175    a point
neos-953928           -99.920       -99.904   a point
air05                   26297         26374   a point, and 0.3% of bound
chromaticindex32-8          3             4   a bound; the incumbent is already exact
neos-820879             25342         25468   a little of both
```

Three instances hold an exact or near-exact bound and are waiting entirely on a good
point, and in HiGHS's log all three are answered by a sub-MIP.

### The sub-MIP construction that does pay, and the trap next to it

HiGHS's `rootReducedCost` builds a neighbourhood out of the lurking bounds this solver
already computes: *suppose* the answer beats the weakest threshold in the table, apply
every entry at or above it, and search what is left. The table read backwards is a
construction rather than a filter, and it was built here because the table was in hand.

`neos-3045796-mogo` closes. Its incumbent goes from 930 to its optimum of -175, where the
ordinary improvement search walks it to 930 and stops: that one fixes columns where the
incumbent and the relaxation agree, and where the incumbent came from a feasibility search
rather than from the relaxation they agree about nothing useful. This construction does
not read the incumbent at all.

The share of the remaining budget it may spend has a narrow window at both ends, which is
the whole of the tuning: at a twentieth `mogo` comes back with 1380, at a fifth it reaches
-175 in two runs of four, at three tenths in three of three, and at two fifths in none of
three, the search having been left too little time to certify what it holds.

**This was nearly discarded on a confounded measurement, which is worth more than the
result.** It was built at the same time as an unrelated change, root improvement, and when
`irp` stopped closing the two were separated by an environment switch that turned this one
off. `irp` still failed, so this was written up as costing `irp` its solve and reverted.
The switch was being read by a binary that had not been rebuilt; the culprit was the other
change, and reverting the pair discarded a working one. **One change at a time, and a
switch that disables something is a claim to be verified like any other.**

The trap beside it is worth recording more than the construction is. The ordinary
improvement search was also moved to the root, on the reasoning that "Searching the
neighbourhood at the root, which did not work" above had been measured when these
instances had no incumbent to improve, and that keeping the jumped points changed the
situation. The situation had changed; the failure had not. `acc-tight2`, `disctom` and
`neos-913984` each hold their optimum after one node and need what is left of the run to
*certify* it, and an improvement search asked to beat a point the bound has already
matched cannot succeed and does not stop trying. All three went from optimal to a timeout
reporting a gap of zero -- which is exactly what the earlier note says happened, in the
same words.

The lesson is not "do not retry reverted work". Retrying it is how `JUMP_QUALITY` was
found to be costing five instances. It is that a note saying *why* something failed
deserves to be read as a prediction and tested against, rather than treated as spent once
the circumstances look different: the guard the earlier attempt needed was described in
that note, and adding it afterwards cost a round of measurement that reading it would
have saved.

### What `neos-820879` needs, which is two things and neither of them small

The closest instance in the set at 0.5%, and worth writing down in full because the chain
that closes it is visible in a reference solver's log and only one link of it is missing
here.

```text
HiGHS       root LP 24874, same as this solver
 J   0.5s   a point, 39392
 L   4.7s   sub-MIP finds 25505, bound moves to 25140
     ...    74.0% inactive integer columns, restarting
            9522 columns -> 2476
            bound 25150 -> 25229 -> 25236, restarting twice more
```

Everything after the first line depends on the second. A point of 25505 against a bound
of 24874 leaves a room of 631, and at that room reduced cost fixing decides three quarters
of the model; restarting on what is left re-presolves and re-cuts a model a quarter the
size, and *that* is what walks the bound from 25140 to the optimum of 25468.

This solver has the fixing and the lurking table and gets no benefit from either here,
because it has no point until the run is nearly over. Its first incumbent arrives after 45
seconds of a 60 second limit; the ordering that produces that is deliberate and is
explained above, the chain being ordered by what its points are worth rather than by what
they cost. On this model nothing above the feasibility search finds anything at all, so
the chain runs to the end before reaching the one thing that would.

Running the feasibility search early was tried and is not the answer either. It reaches a
point in 2.5 seconds, but the point is 38973, which leaves a room of 4000 and fixes
**nothing**. What matters is not having a point early but having a *good* one, which is
what the sub-MIP gives HiGHS and what nothing here produces in time.

It costs, too, and the cost falls where this file has recorded it falling before:
`f2gap401600` 0.25 seconds to 1.8 and `mod010` 0.53 to 4.5. Bounding it by flips per
column rather than by a share of the clock, which is the correction this file applies
everywhere else, changes nothing: the cost *is* the flips, and 200 per column on a model
of 2655 is four seconds. Reverted.

So the two missing pieces are a sub-MIP good enough to reach 25500 within a few seconds,
and restarts to spend what it earns. The second was measured earlier against `n2seq36f`
and `neos-1516309` and found not to move their bounds, which is still true and is not the
same claim: here the bound demonstrably does move, four times, and the difference is that
those two had nothing to fix and this one has three quarters of its columns to fix. A
restart is worth what the fixing before it was worth.

### The sub-MIP for `neos-820879`, and four ways it does not work

The instance needs a point near 25500 within a few seconds, which a reference solver gets
from a sub-MIP at 4.7 seconds. Four constructions were built and measured, and the
negative results are worth more than the attempt.

**RENS, fixing the relaxation's integral columns.** 9260 of `neos-820879`'s 9522 integer
columns sit exactly on an integer in the relaxation. Fixing all of them leaves a model
whose **LP is feasible and whose MIP is not**, and the same holds at every rate tried from
30% to 90%. Fixing fewer does not help, because the infeasibility is not caused by the
last few fixings; it is the shape of the neighbourhood.

**RENS with propagation.** Fixing them one at a time, propagating each through the
conflict graph and row activities and skipping any whose consequences contradict what is
known — which is the expensive, careful version — reaches 3333 of 9522 decided and the
sub-MIP is *still* infeasible. Propagation sees what one row can prove about another and
that is not what is wrong here.

**RENS with the LP as the oracle.** Bisecting on the number fixed and asking the LP after
each attempt finds that the largest feasible fixing is *all of them*: the LP stays optimal
with 9260 columns pinned. The neighbourhood is LP-feasible and integer-infeasible, which
is exactly the case an LP cannot detect and the reason the careful version above did not
help either.

**The dive's leftovers.** A dive that dead-ends has fixed columns that each survived an LP
solve, which is a stronger statement than propagation makes, and it throws them away when
it reports failure. Keeping them as a neighbourhood is nearly free to build. It reaches
200 columns of 9522, because that is the dive's step budget, and two per cent of a model
is not a neighbourhood.

What is left is RINS, which this solver has and runs every 500 nodes. Giving its
sub-search a feasibility search of its own — it is disabled there, on the reasoning that a
search starting from an incumbent is not looking for feasibility — does help: `neos-820879`
goes from no incumbent at all to 25860 and `air05` from 47030 to 45318, and `eil33-2`, the
instance that reasoning came from, is unchanged. It also costs `neos-3045796-mogo` three
closes in three, dropping it to one in four, because the feasibility search inside every
neighbourhood spends the budget `mogo` needs. Reverted on that.

So the thing the reference solver does in four seconds is not any of: fixing what the
relaxation is sure of, fixing it carefully, fixing as much as the LP allows, or keeping
what a dive has already paid for. That is a real narrowing of the search space and the
next attempt should start from it rather than from the same four ideas.

### The primal side is not the gap, measured by removing it

Before building a suite of primal heuristics, the thing to establish is what a perfect
one would be worth. `search::solve_from` takes a point the caller already has, and
`cargo run --release --example primalceiling` hands the search a reference solver's
*optimum* and asks whether it can then prove optimality. An instance short of a point and
one short of a bound both time out with a gap and want opposite work; this separates them.

```text
                     seeded with the true optimum      gap left
neos-3045796-mogo    Optimal                            0%
neos-953928          TimeLimit                          0.0156%
neos-820879          TimeLimit                          0.69%
air05                TimeLimit                          0.77%
neos18               TimeLimit                          50%
mine-166-5           TimeLimit                          22%
```

**Handed the answer, five of the six still cannot prove it.** So a primal heuristic suite
converts none of them, and the classification above -- which called three of these primal
bound -- was wrong. It compared this solver's bound against the optimum, and this solver's
bound was measured in runs where the incumbent was poor; with a perfect incumbent the
bound improves and still falls short. The right question was never "is our bound below the
optimum" but "would a better point let us prove it", and only the second one can be
answered by asking.

`neos-953928` is the near miss and is worth its own line: with the optimum in hand it is
0.0156% short of proving it, against a default gap tolerance of 0.01%. Its objective is
not integral, so no rounding of the bound helps.

### A permuted vector looks exactly like an infeasible point

The first run of this measurement said `neos18` was seeded with a point scoring 40 against
a true 16, violating 536 rows, and `neos-953928` with one violating 92. Both were read as
infeasible reference solutions.

They were not. **Two readers of the same MPS file need not agree on column order, and
these two do not**: on `neos18` this solver's first column is `r_0` and the reference
solver's is `x_1_0`. Matching by position makes the seed a permutation of the answer,
which scores wrongly and violates rows exactly as an infeasible point would. Neither
reader is wrong; column order is an internal matter, and the file names its columns for
precisely this reason.

Two instances happened to agree on order, which is what made the failure look like a
property of the other two rather than of the method. Anything crossing a solver boundary
has to be keyed by name, and a check that the reference point is feasible before it is
believed is what turns this from a wrong conclusion into a caught mistake.

### `air05` and `neos-953928`: the filter is right and the face is large

`air05` produces no Gomory cut at all from a relaxation with 224 fractional columns,
which is the symptom that turned out last time to be a filter discarding them silently.
It is that again, and this time the filter is correct.

```text
air05        367 basic integer rows, 217 usably fractional, 150 too near an integer
             all 217 cuts dropped for density, every one of them exactly 6557 terms
             against an allowance of 2398 and a model of 7195 columns
```

Ninety-one per cent dense, and uniformly so. One such cut carries 6557 nonzeros into a
matrix that has 52121 altogether, so it is an eighth of the model per row; the two hundred
of them are twenty-seven times the matrix, and adding them does not finish an LP solve in
two hundred seconds. Taken a handful at a time they do move the bound -- three cuts are
worth 10, ten are worth 15, thirty are worth 18 -- and `air05`'s *clique* cuts are worth
29 for a small fraction of the cost and are already getting in. So the filter is keeping
out cuts that are both ruinous and worse than what is already there, which is what it is
for.

None of it is nearly enough either way. `air05`'s root is 25877.6 against an optimum of
26374, and seeded with that optimum the search still reaches only 26170. The root is short
by five hundred and the cut families here argue about thirty.

`neos-953928` answers even more plainly. Its Gomory cuts are dropped the same way, 348 of
them at sizes from 1519 to 21475 against an allowance of 312, and allowing all of them
changes the bound by **nothing at all**: -99.9200 with five, with thirty, and with a
hundred and nine across every family. That is the large-optimal-face case again, and what
answered it on `n2seq36f` does not apply here, because mod-2 separation finds no cut on
this model at all.

So neither instance is short of a filter setting. `air05` wants a bound five hundred
better than any family here produces, and `neos-953928` wants a family that can cut a face
it has no cut for.

### `neos-820879`: the cut families are exhausted, and restarting does not refill them

The one instance left where more of what already works looked like it would be enough. It
is not, and the reason is worth having because it closes the line of enquiry.

The root cut loop **converges in three rounds** and then no family finds a violated cut at
all:

```text
round 1   70 offered, 31 kept    bound 24874.27
round 2   66 offered, 20 kept          24969.21
round 3   18 offered,  7 kept          25067.83
round 4    0 offered                   25113.87
```

Raising the round limit from 10 to 200 changes nothing: 58 cuts and 25114 every time. The
five families here are exhausted at 25114 and the bound needs 25468.

Handed its own optimum, the search reaches 25258 in a minute from that root, which is 144
of the 354 it needs. So neither the bound nor the point is the whole story; the root is
short and the search closes less than half of what remains.

**Two defects were found chasing this and neither is worth its cost.** Recorded so the
next attempt does not rediscover them and does not reapply them without a reason:

- The restart trigger measures columns decided *since the pass began*, and a pass begins
  after the root's own reduced cost fixing. On `neos-820879` handed its optimum that
  fixing decides 6871 columns of 9522 — 72%, against the 74% a reference solver reports at
  exactly this point before restarting — and the trigger saw 2651 free and nothing newly
  decided. Firing on the root's own fixing instead makes the restart happen.
- Every root budget that is a share of the run is a share measured from the *solve's*
  start, so on any pass after the first they are all already spent. With the restart
  firing, the second pass separated not one cut, because the cut loop's third of the run
  had ended twenty seconds earlier. A restart that cannot redo the root work is a restart
  for nothing.

Both were built, and with both in place the restarted pass runs its cut loop on a model
72% fixed and offers **zero cuts at the same bound of 25113.87**. The families are
exhausted whatever the size of the model, so the restart has nothing to refill them with.
The pair also cost `irp` four closes in five, dropping it to one, and were reverted on
that: two fixes to machinery that does not pay, for an instance they do not convert.

What a reference solver does differently at this point is not only restart but *physically
shrink* — 9522 columns to 2476 — where this solver keeps every fixed column in the matrix,
by the deliberate choice recorded under "Presolve" not to renumber. A restart is the
moment that choice costs something, and testing whether a compacted model yields cuts the
uncompacted one does not is the only remaining lead on this instance.

### Compaction, and the question it finally answers

"Physical removal is a later optimization, worth doing when instance sizes make the wasted
rows matter" has stood under "Presolve" since the beginning, and a reference solver does
exactly that on restart, taking `neos-820879` from 9522 columns to 2476 where this solver
keeps every fixed column in the matrix. That is the last untested explanation for why its
cut families run dry.

`compact` is built, and the answer is no.

Simulated round by round -- cut, fix against the incumbent on the cut model's own reduced
costs, compact, repeat -- which is the sequence a restart performs:

```text
uncompacted   9522 cols, 82 cuts, bound 24958.02
round 0       3934 cols, 86 cuts, 24874.27 -> 24968.23
round 1       3404 cols, 110 cuts,          -> 25032.53
round 2       3058 cols, 87 cuts,           -> 25108.20
round 3       2665 cols, 0 cuts             -> 25108.20
round 4       2647 cols, 0 cuts             -> 25108.20
```

The model shrinks to within a couple of hundred columns of what the reference solver
restarts on, and the bound converges to 25108, which is where the *uncompacted* loop also
stops. From round three no family finds a single violated cut. So the exhaustion is a
property of the families, not of the representation, and a reference solver reaching 25468
does it with cuts this one does not have rather than with a model shape this one cannot
build.

That is worth knowing beyond this instance: three separate explanations for `neos-820879`
-- restarts, sub-MIPs, and now compaction -- have each been built and each turned out to
be downstream of the same thing.

**The module is kept and nothing calls it.** It is the tested answer to a question this
file has carried from the start, it is the prerequisite for any restart that does pay, and
the general version was measured before and found "correct, and 14% slower for nothing",
so wiring it in without a reason would repeat that. Its own cost is a module and its tests.

The tests are the part worth reading. A postsolve that renumbers wrongly returns a
plausible vector for the wrong columns, which scores as a valid objective and is not one,
so correctness is checked by solving both models and comparing rather than by reading the
map: the compacted optimum against the original's, *and* the expanded point scored on the
original against what the compacted search claimed for it. Only the second catches a
permutation, since a permuted vector still scores something. Both were checked against
deliberate breakage -- expanding by position instead of through the map, and dropping the
fixed columns' contribution to the objective -- and each fails on the first seed that
compacts.

### The losses re-derived from node counts, which points at presolve and then away again

Every cut section above works on the three instances closest on the bound. Counting nodes
across all eighteen that HiGHS closes and this solver does not says those three are not
the largest group. Seven of the eighteen reach **one or two nodes in a minute**: their
root relaxation, or the loop around it, spends the whole budget, and nothing built above
the LP can reach them.

Comparing presolved sizes says where that comes from, and it is not subtle. HiGHS's
presolve leaves `ex9` at **0 rows and 0 columns**, `ex10` at 14 x 6, `bley_xl1` at 9.7%
of its rows, `neos-780889` at 33%, `neos-633273` at 25%, `neos18` at 23% of its columns.
The instances it *cannot* reduce -- `neos-820879`, `chromaticindex32-8`, `ab71-20-100`
and `neos-953928`, all at 97% or more -- are exactly the ones already answered here as
bound-limited. The two groups do not overlap at all.

That is a strong enough signal to check before building anything, and checking it mostly
refuted it. Dumping HiGHS's presolved models and handing them to this solver:

```text
ex10           optimal in 1ms          (presolved to 14 x 6)
neos18         531 nodes, gap 52%      (127 nodes and 66% on the original)
supportcase4   2 nodes                 (4370 x 1112)
neos-633273    2 nodes                 (5520 x 5541)
bley_xl1       1 node                  (17039 x 751)
```

Four of the five still reach one or two nodes on a model a fifth the size, so for them
presolve was never the obstacle. `bley_xl1` is the sharpest: 17039 rows, **751 columns**,
57327 nonzeros, and its root relaxation alone takes 337 seconds over 83581 iterations at
247 a second. That is an LP defect, not a presolve one, and it is the same shape the
hyper-sparse experiment above was rejected on -- but rejected against models with many
columns, where `B^-1 a_q` is 12-33% dense. A basis of 17039 rows holding at most 751
structural columns is nearly all logical, which is the case that measurement did not
cover. Entering through the dual instead does not finish it in ten minutes either.

So the presolve lead survives on `ex9` and `ex10` alone, and those two were the pair the
note above retired as needing "presolve well beyond what probing reached".

### Which rule collapses `ex9`, asked of the solver that does it

HiGHS takes a `presolve_rule_off` bitmask, so the question can be put to it directly
rather than guessed at from names -- the mistake the mod-2 section records. Turning each
rule off in turn changes nothing: `ex9` still collapses to nothing, so no rule is
necessary. Turning every rule off *but* one is the question worth asking:

```text
only (none)                     33846 x 8560
only Probing                      394 x 17
only Enumeration                  548 x 87
every other rule alone          33846 x 8560
```

Two things fall out. Probing alone does essentially the whole of it. And "no rules at
all" leaves 33846 x 8560, which is **exactly what this solver's own cheap presolve
already reaches** -- 7116 rows and 1848 columns removed from 40962 x 10404. The gap on
`ex9` was never the reductions this solver is missing. It was probing, which it has.

### What probing costs, which is one of its two directions

Instrumenting the two sweeps of each probe separately gives the shape of it on `ex9`:

```text
zero sweeps   1612 fixings for      99,496 work
one  sweeps   2,714,582 fixings for 915,722,433 work
```

Supposing `x_j = 0` costs 62 entries and derives *nothing at all* -- one fixing per
probe, which is the probed column itself. Supposing `x_j = 1` costs 568,000 and fixes
1684 columns, a fifth of the model. In a set packing row a one excludes everything beside
it and a zero excludes nothing, so all the cost and all the information sit on one side.

The information is then thrown away. Unless the *other* sweep is refuted, the loop
rewinds and moves on, so 1684 implications per probe are re-derived and discarded, and
only the 158 refutations of 1612 probes are kept.

### A column both suppositions force is forced

The model takes one of the two values, so a column both sweeps fix to the same value is
fixed outright, whichever way the probed column goes. Both sweeps have been paid for
already, and the second one's fixings are already listed on the trail, so the rule costs
one lookup table and a pass over what the sweeps touched.

It is not more *powerful* than probing: a column it decides could also have been decided
by probing that column directly. What it has is reach. Probing starts only from a binary
the conflict graph has an edge for, and decides the column it started from, so a column
with no edge of its own is never the subject of a probe. On `ex9` the rule fires 28 times
and decides 2920 columns, and the difference it makes is **8244 columns fixed against
10076** -- which is the difference between timing out and closing.

It cannot reach a general integer, because row propagation here only forces binaries, and
that is worth knowing before reaching for it on a mixed model.

### What was actually stopping `ex9`, which was not what the budget note said

Given budgets it can finish on, probing fixes 10076 of `ex9`'s 10404 columns in 22
seconds, and the whole model then closes at 25 seconds against a 60 second limit,
objective 81, in two nodes. `ex10` reaches 17272 of 17680 but wants 86 seconds, which no
60 second limit affords.

Three separate caps were in the way, and only one of them was the one the notes blamed.

`PROBE_TOTAL`, the total work cap, was **not** it: raising it alone changes nothing
anywhere. `PROBE_PATIENCE`, the work allowed without proving anything, was. `ex9`'s
proofs are front-loaded -- four in its first eighteen probes -- and it then goes 57.8
million entries between two proofs before finding another 2600. Thirty million cut it off
inside that gap, at 63 probes of 8560, and what has been recorded here twice as a model
probing could not reduce was a model probing was not allowed to finish reducing.

Raising patience is what the old note warned costs `irp` two seconds, and it does: at 200
times, `irp`, `nw04`, `air03` and `eil33-2` spend 9, 23, 10 and 9 seconds each and prove
**nothing whatever**. But the reason is visible once the sweeps are counted rather than
the work, and it separates the set exactly:

| | sweeps abandoned | proofs |
|---|---:|---:|
| `irp`, `nw04`, `air03`, `eil33-2` | **every one** | 0 |
| `ex9` | 0 of 7548 | 2636 |
| `ex10` | 45 of 15170 | 3888 |
| every other model in the set | 0 | -- |

Those four abandon every sweep they start: each runs into `PROBE_REACH` and gives up. A
sweep that did not finish has not said the column has nothing to give, only that it was
not allowed to look, and a run of sixteen of them says the model is out of probing's
reach. That is `ABANDON_RUN`, and it is what `irp` needed all along -- a guard on the
shape of the failure rather than on the size of the bill. It costs the models that pay
nothing at all, because they abandon nothing.

With patience at 120 million, `ABANDON_RUN` at sixteen, and `PRESOLVE_SHARE` widened from
a quarter to two fifths so a 22 second presolve fits inside a 60 second run, `ex9` closes.

### `ex10` is `ex9` with the bill four times over, and the second half is not the model

`ex9` and `ex10` are the same instance family and answer differently. Given patience of
1200 million rather than 120 -- `ex10`'s dry stretch is between 600 million and that,
where `ex9`'s is 57.8 -- `ex10` closes: objective 100, two nodes, **161 seconds**. So it
is winnable, and not at sixty.

The 161 splits into 78 seconds of probing and 83 of everything after it, and the second
half is the surprise. After probing, `ex10` has 408 free columns of 17680 and still
carries all 69608 rows and every fixed column into the basis, which is exactly the case
compaction was built for and has never been tried on. Compacting it reaches **16760 rows
by 408 columns**, a fortieth of the nonzeros, and:

```text
compacted      Optimal in 119.0s, 2 nodes
as presolved   Optimal in 121.0s, 2 nodes
```

Two seconds of 120. The cost was never the rows being carried; it is one root relaxation
that takes two minutes whatever shape the model is written in. That is the second time
compaction has been measured against a case that looked made for it and returned nothing,
and the two together are a fair summary of it: it is correct, it is cheap, and the models
where the size looks like the problem are models where the size is not the problem.

Parallel probing would take the 78 seconds to perhaps ten. It would not touch the 119.

### `supportcase4`, and why better presolve is not its answer either

Its objective is zero, so the whole model is a feasibility question and the
constant-objective path already gives the feasibility search three quarters of the run.
Feasibility Jump takes the violation from 1558 to about 300 and stops there, restart after
restart. Nothing else finds a point either: at a **600 second** limit it reaches one node
and no incumbent, because on this solver's own presolve the root relaxation does not
finish in ten minutes.

HiGHS's presolve reaches 4370 x 1112 where this one reaches 6418 x 2136, and asking which
rule does it says **Parallel rows and columns** -- alone, it is the whole reduction, as
are Probing and Enumeration alone. Parallel row merging is in this history, built and
reverted because freeing a row did not remove it, with the note that compaction was the
prerequisite it needed. Compaction now exists, so that is a coherent piece of work.

It is also not worth doing for this instance, and the experiment that says so was already
run: handed HiGHS's presolved `supportcase4` directly, this solver still reaches **two
nodes**. The root relaxation on 4370 x 1112 takes 10.8 seconds and the run then goes
nowhere. Presolve is not what `supportcase4` is short of.

### Where the remaining conversions actually are, which is one place

Putting the sixty-second losses together after all of the above, the instances that are
neither bound-limited nor already answered come down to one cause:

```text
ex10           root relaxation 119s on 16760 x 408 after everything
bley_xl1       root relaxation 337s on 17039 x 751, 83581 iterations at 247/s
supportcase4   root relaxation does not finish in 600s on 6418 x 2136
neos-633273    two nodes in sixty seconds
```

All four are many rows against few columns, and all four are LP throughput and nothing
else. Presolve has been tested against them and converts one instance, `ex9`; compaction
has been tested against them and converts none; the cut families do not apply, because a
bound needs a relaxation to bound. The next conversion at this limit is a faster
relaxation on a row-heavy model, and the rest of the machinery is waiting on it.

One caution before starting there. The hyper-sparse triangular solve was implemented and
reverted, on the measured grounds that `B^-1 a_q` runs 12 to 33% dense; a basis of 17039
rows holding at most 751 structural columns is nearly all logicals and is the shape that
measurement did not cover. That is a hypothesis and not a finding: an attempt to measure
it here counted every `ftran`, including the dense `B^-1 b`, and came out at 95% on all
three models tried, which says only that the measurement has to isolate the entering
column's solve to mean anything.

### Asking the other two solvers what they know that the first one does not

Every comparison here has been against HiGHS. Ten of the 140 are closed by SCIP or CBC
and missed by HiGHS, which is a different question and a better-posed one: whatever does
those is a technique HiGHS lacks. Four of the ten this solver already closes -- `disctom`,
`decomp2`, `nw04`, `neos-4754521-awarau`. SCIP takes parameters for its own families, so
the remaining six can be put to it directly rather than guessed at:

```text
cod105             baseline 57.7s   no symmetry -> timeout, gap 1.03
graph20-20-1rand   baseline  6.1s   no symmetry -> timeout, gap 4.61
neos-787933        baseline  1.2s   no presolve -> timeout, gap 23.1
tanglegram6        baseline 19.0s   nothing changes it
```

Two answers. `cod105` and `graph20-20-1rand` are **symmetry**, and neither has a single
duplicate column, so it is a real automorphism group rather than interchangeable columns
-- nauty-class machinery, and nothing smaller will do. `tanglegram6` is nothing in
particular; every knob off, still nineteen seconds.

Worth noting for both of the symmetric ones: turning SCIP's *cuts* off makes them
**faster**, `cod105` 57.7s to 11.5 and `graph20-20-1rand` 6.1 to 2.3. Cuts are not what
that group wants.

### A relaxation that is weak by construction, and the family for it

`neos-787933` is presolve, and narrowing it further is worth the detail. Disabling SCIP's
constraint-handler presolving leaves 63708 columns -- which is exactly where HiGHS's
presolve lands, and exactly where this solver's own lands. All three agree; SCIP then goes
to 1764 columns and 439 rows, and handed that model this solver closes it in 22 seconds.

The model is 1764 big-M linking rows, `sum_{j in S_k} x_j - 133 y_k <= 0` with
`|S_k| = 133`, minimising `sum_k y_k`, over covering rows wanting three of a set. Its LP
relaxation is **3.0** against an integer optimum of **30**. That is not a relaxation that
happens to be weak; `y_k` is free to sit at a hundred-and-thirty-third of what any one
member of its group holds, and no amount of separating from single rows repairs it.

Each such row says `x_j <= y_k` for all 133 members. The conflict graph took each row's
longest overshooting prefix and stopped, which for this row is `{not y_k, x_1}` -- **one
implication of the 133**. The caution recorded under "Reading conflicts out of long rows"
said exactly this ("leaves cliques behind on general knapsack rows") without anything
depending on it. Sorted by weight, a literal's conflicting partners are a prefix of the
list, so all of them come out in one walk: the graph on `neos-787933` goes from about 1764
edges to 234612.

Then the family. For a row `sum a_j x_j >= b` with `a_j > 0`, every term obeys
`a_j x_j <= a_j y_j`, so

```text
    sum_j a_j y_j  >=  sum_j a_j x_j  >=  b
```

is valid, and columns sharing a bounding variable collect onto it. Measured before any
Rust: 133 such rows take the LP from 3.0 to **30.0**, the integer optimum, on the nose.
Built, it is 77 cuts and the root bound goes 9 -> 30.

One thing had to change for it to reach the LP at all. `select` ranks by efficacy, which
divides violation by the cut's norm and so prices a six-hundred-term aggregated row
against a two-term clique; on this model the two-term cliques win, move the bound from
8.07 to 8.16, and crowd out the rows that move it to 30. The ranking is not wrong in
general -- a short cut really is cheaper to carry at every node below -- so the family
gets a reserved quarter of the round rather than a thumb on the scale.

With the bound exact the incumbent still stalled at 64 after 600 seconds, so what was a
bound failure became purely a primal one, and the answer was the other half of the same
observation. Where `x_j <= y` holds, the model rewritten with every member replaced by its
gate is a *restriction* of the original: every point of it maps to a point of the model,
so it can be searched and the point it returns checked rather than trusted. `neos-787933`
rewrites to 1764 columns and 133 rows, which this solver closes in **a fifth of a second**.

One thing had to be chained on. Substitution leaves every column presolve has already
decided sitting in the model at its fixed value -- 172668 of the 174432 that come out of
it -- so the rewritten model *looks* almost as big as the original and the guard against
searching the same question again rejects it. Compaction is what turns that back into the
1764 the model is about, and this is the first thing in the solver to call it. It was
committed as a tested capability nothing used, on the argument that it was the
prerequisite for work that had not been done yet; this is that work.

`neos-787933` closes: **objective 30, optimal, 39 seconds**, where HiGHS does not close it
at all.

### What a marginal instance costs a count, and `neos-3045796-mogo`

The run that gained `neos-787933` lost `neos-3045796-mogo`, so the count stayed at 34.
That is not what it looks like. `mogo` has closed at 38 to 41 seconds of a 60 second
budget in every baseline since it was gained, and it is now a coin flip: three runs of the
current build close one, and four runs with the aggregated search *removed* close two. It
is marginal on its own account and was already sitting on the line.

The claim in `handoff.md` that "the 34 all close on every run, with no instance left
sitting on the sixty second line" is therefore wrong, and was wrong when written. The
discipline it states is right and was not applied to `mogo`: quote the set that closes
every time. That set is 34, and `neos-787933` -- three runs, three optimal -- is in it
while `mogo` is not.

### Three of three that was one of six

`neos-3045796-mogo` is the cheapest-looking instance left: its root bound is already the
optimum before a single cut, so it is purely primal, and it closes at 38 to 41 seconds of
a 60 second budget. Making it reliable is the 35th instance.

Tracing the root cut loop on it says where the time goes. It adds 64 cuts a round for
forty rounds, 2240 of them, and the bound does not move by a digit. They are not
duplicates -- 2179 of the 2240 are distinct and 1837 are still carried -- so the
separator keeps finding genuinely new violated cuts on a large optimal face, and the model
the search then runs on has grown by **82% of its own row count**. `decomp1` reaches 10%
and `decomp2` 5%, which looked like a clean discriminator and a principled rule: a pool
that size is replacing the relaxation rather than strengthening it.

Three settings, three runs each, said cut-rounds 5 closed **three of three** at 23 seconds
against a default of two of four at 37. Six runs each says otherwise:

```text
default              closed 2 of 6
cut pool <= rows/2   closed 3 of 6
cut-rounds 5         closed 1 of 6
```

The three-of-three was luck, and reverses. None of the settings is distinguishable from
the default or from each other; `mogo` is a coin flip at roughly two in five whatever the
cut budget is, because what decides it is whether the parallel search stumbles onto -175,
and the workers do not visit the tree in the same order twice.

This history already says "repeated runs at sixteen threads spread by several points of
final gap" and to quote the set that closes every time. What it did not have was a case
where a three-of-three sweep reversed outright, and that is the useful part: three runs is
not a measurement of anything on this set, and a change that looks like it converts an
instance on three runs has not been measured at all. The experiment is reverted.

Also worth keeping, because it retires a plausible-looking rule: flat root rounds are
**not** waste. `decomp1` and `decomp2` never move their root bound during cutting either,
and stopping early cost them 50 to 80% when it was tried. Cuts that do nothing at the root
still bind at nodes below, where more columns are fixed. "The bound is not moving" is
measuring the wrong thing, which is why the earlier attempt at this failed and why the
cut-pool version was worth trying instead -- and it fails for a different reason, which is
that there is nothing there to measure.

### One seed in a thousand generated no model at all

Sweeping the model generator for a shape that reaches the agreement rule hung it.
`Kind::Signed` draws a binary witness `x*` and then redraws each row's coefficients until
`a'x*` is nonzero, so an all-zero witness makes that sum zero for every draw and the
redraw never ends. It is one witness in `2^n_cols`: at ten columns, `signed_c10_r8_s56`
is the first seed to hit it, which is rare enough to have sat unnoticed and common enough
for any sweep over a couple of hundred seeds to find. Drawing the witness again when it
comes up empty is the whole fix, and it changes only the seeds that produced nothing.

Worth recording for the method rather than the bug: this took an hour to find because the
symptom appeared inside presolve's probe loop, three layers from the cause, and two
rounds of instrumentation went into presolve before the trace showed the last probe
completing and nothing after it. The print that would have found it immediately was the
one naming the *next* model, and it was sitting after the parse rather than before it.

## Measuring the search

The parallel search is not deterministic, and single runs of it are not evidence.
Workers expand nodes in whatever order they finish, so the incumbent found first
differs between runs of the same binary on the same model, and everything downstream
follows it. Repeated runs at sixteen threads spread by several points of final gap:
`neos-911970` over 19.01%, 21.59% and 22.61%, `beavma` over 9.99%, 13.17% and 11.54%.

That is wider than most of the changes worth measuring, and it has already produced a
false result recorded in this history. A single sweep showed `cap6000` improving from a
4.08% gap to 0.58% and the improvement was attributed to a change in the dual method.
Rerunning it three times on that build and three more on the build without the change
gives 4.0816% every time. The 0.58% was one lucky ordering.

Single threaded, the same solve is exactly reproducible: `beavma` returns 8.7117% three
times running. So an A/B of an algorithmic change belongs at one thread, and the
sixteen-thread numbers belong only in headline results where the comparison is against
another solver rather than against another build of this one.

The trap this falls into is measuring the loudest number rather than the one the change
acts on. A related case earlier: a change to the dual method's stall detection was read
as harmless because final gaps barely moved, when it had cut the node count on
`drayage-25-23` from 98 to 4. The gaps could not move, because that instance finds no
incumbent either way, so the only number that could show the damage was the one not
being watched.

## Correctness practices

Several bugs in this project were invisible to unit tests and were caught only by
differential testing against another solver. Two patterns recur.

Tests that cannot fail. An early set of LU scale tests gave every basis a strong
diagonal, which made pivot choice irrelevant. They passed with the bad pivot search
deliberately reintroduced. Any test asserting that an optimization is correct should
be checked against the unoptimized path, and any threshold should be checked by
disabling it and confirming the test then fails.

A guard that quietly undoes itself. A search that proves nothing sets its bound to
`NaN`, deliberately, so that a bound left over from an earlier phase cannot be read as a
proof. The gap was then computed with a `.max(0.0)` on the end, and `f64::max` returns
its non-`NaN` argument, so the guard came out the other side as `gap 0.0000%` — the
reading it existed to prevent, and the most reassuring one available. `neos-954925`
reported it on every run whose root relaxation did not finish, which is every run of it.
A sentinel is only a guard if every path that consumes it knows it is a sentinel.

Circular verification. A check written to confirm two pivot orderings agreed compared
a shortlist against a stable sort of the same ordering it came from, so it could not
detect the difference it existed to find. Ground truth has to come from somewhere the
code under test did not produce it.
