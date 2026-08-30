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

Circular verification. A check written to confirm two pivot orderings agreed compared
a shortlist against a stable sort of the same ordering it came from, so it could not
detect the difference it existed to find. Ground truth has to come from somewhere the
code under test did not produce it.
