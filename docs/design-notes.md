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
