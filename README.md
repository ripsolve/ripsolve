# ripsolve

A branch-and-cut solver for **mixed-integer programs**, in pure Rust.

`ripsolve` is an LP-relaxation-based branch-and-cut solver: every node solves a
bounded LP relaxation with the simplex method, and the resulting dual bound —
strengthened by presolve and cutting planes — is what prunes the search.

Columns may be binary, general integer, or continuous. A binary variable is just an
integer bounded to `[0, 1]`; the solver has no separate notion of one, because
branching splits a range (`x <= floor(v)` against `x >= ceil(v)`) and that
degenerates to fixing at 0 or 1 on its own.

## Scope

`ripsolve` aims to be **fast on small and medium models** — up to roughly a
thousand rows — and it is: on models in that range it is quicker than HiGHS, SCIP
and CBC on several, ties a leading commercial solver on some, and beats it on a few.

It is **not** a large-model solver, and the distance is not small. On a seeded
random sample of 25 models from the MIPLIB 2017 benchmark set, whose median size is
11,697 columns, it proved optimality on one and failed to finish even the root node
on eleven. Closing that would take an efficient sparse simplex with partial pricing
and steepest-edge, Forrest-Tomlin updates, a presolve that does real work at scale,
and numerics to match. None of that is here.

So: reach for it on models of a few hundred to a few thousand rows. Reach for HiGHS
or SCIP above that.

What is implemented: LP and MPS readers, presolve, lifted knapsack cover cuts,
Gomory mixed-integer cuts, pseudocost branching, primal heuristics, a sparse LU
factorization, best-bound search, parallel tree search, and a `gurobipy`-shaped
Python API. Quadratic terms, SOS constraints, callbacks and lazy constraints are
not.

Correctness is checked against a leading commercial solver as an oracle —
relaxation values, integer optima, and 400 randomized mixed-integer models — and
against HiGHS, SCIP and CBC, which must all agree wherever more than one proves
optimality. Nothing links against any of them and the test suite runs without them.

## Why not pure enumeration

The predecessor to this project implemented Balas' additive algorithm — implicit
enumeration with a feasibility-based bound. Measured on a dense 60-variable
instance, that approach explored 10.3 billion nodes in 54.7s across 16 threads.
A leading commercial solver did the same instance in 125 nodes and 0.15s.

The enumeration engine was not slow; at 188M nodes/sec it processed nodes roughly
1000x faster than the commercial solver does. It lost because it visited 82 million times more of
them. No amount of SIMD or threading closes a gap of that shape — only a stronger
bound does, and that means solving an LP relaxation at each node. That measurement
is the reason this project exists and why the simplex is its centrepiece.

On that same 60-variable instance, `ripsolve` now proves optimality in **89ms on a
single thread** — against the enumeration solver's 49.7s across sixteen. On a family
of dense random instances with 30 rows, enumeration on sixteen threads against the
other two on one:

| variables | enumeration (16t) | ripsolve (1t) | commercial (1t) |
|---:|---:|---:|---:|
| 60 | 49.7s | **0.09s** | 0.33s |
| 64 | >120s | **0.16s** | 0.17s |
| 80 | >120s | **0.24s** | 0.25s |
| 200 | >120s | 6.8s | 1.6s |
| 500 | — | 1.7% gap at 300s | 294s |

## Against the open-source solvers

HiGHS, SCIP and CBC, every solver on one thread with the same 60-second limit,
timing only the solve. `bench/open_compare.py` reproduces it, and fails if any two
solvers that both claim optimality disagree on the optimum — it is a differential
test as much as a benchmark.

| instance | ripsolve | HiGHS | SCIP | CBC |
|---|---:|---:|---:|---:|
| `v032c032` | **0.04s** | 0.04s | 0.04s | 0.04s |
| `v048c048` | **0.04s** | 0.1s | 0.3s | 0.04s |
| `v048c128` | **0.1s** | 0.1s | 0.2s | 0.1s |
| `v064c064` | **0.1s** | 0.3s | 0.8s | 0.2s |
| `v064c200` | 2.8s | 2.8s | **2.0s** | 4.1s |
| `v081c162n009` | **1.1s** | 1.6s | 2.5s | 2.5s |
| `v081c162n018` | **0.4s** | 0.7s | 1.6s | 0.7s |
| `v128c256n100` | **0.1s** | 0.2s | 0.2s | none |
| `v256c256n100` | **0.5s** | 1.0s | 1.6s | none |
| `v128c1000n100` | 6.1s | **2.9s** | 2.9s | none |

Fastest on seven of ten, and every proven optimum agrees. The two it loses are the
largest, which is the same boundary the scope section describes.

## Against a commercial solver

Both solvers given the same thread count and the same 60-second limit, timing only
the solve. Thirteen of the fourteen models are solved to proven optimality, matching
the commercial solver's objective exactly every time.

| | single-threaded | sixteen threads |
|---|---|---|
| within 1x | `v032c032`, `v048c048`, `v048c128`, `v064c064`, `v081c162n018`, `v128c256n100` | those, plus `v064c200`, `v064c1000n100`, `v128c1000n100`, `v081c162n009` |
| 2–4x | `v081c162n009`, `v256c256n100`, `v064c200`, `v128c1000n100`, `v064c1000n100` | `v048c128`, `v256c256n100` |
| worst | `mkp_200` at 20x | `mkp_200` at 6x |

`v081c162n009` is *faster* than the commercial solver at sixteen threads — 0.4s against 1.0s.

Two models are not solved by either. `v064c1000n020` yields no feasible point to
anything, the commercial solver included. `mkp_500` it solves in 294s single-threaded and
by `ripsolve` not at all; at the 60-second limit it sits at a 5% gap.

`mkp_200` is the remaining outlier and worth naming rather than averaging away: 91k
nodes against the commercial solver's 18k, and a per-node cost that parallelism improves but does
not close.

### How these were measured

The commercial solver is deliberately unnamed. Its licence here is an academic one,
and rather than read its terms as permission to publish benchmarks under its name,
the comparison omits it. Nothing else is withheld: every number is real and
reproducible by anyone holding a licence, and `bench/compare.py` names the solver it
calls.

One machine, 16 cores, `--release`. The commercial solver is pinned to the same thread count as
`ripsolve`; left to itself it takes every core, which would compare deployments
rather than algorithms. Only the solve is timed on both sides — gurobipy's
`Runtime` excludes parsing, and so does `ripsolve`'s reported elapsed.

`bench/compare.py <seconds> <threads>` reproduces the table. The `v*` models come
from `samples/` and from a separate generator; the `mkp_N` models are generated on
demand by `ripsolve gen --kind knapsack --cols N --rows 30 --seed 42`.

The sections below record what each design decision was worth *at the time it was
made*, measured against whatever the instance set was then. They are kept as
reasoning rather than restated as current numbers, so a figure there will not always
match the table above — `mkp_200` in particular now names a harder generated
instance than it did earlier in the project.

Presolve reduces in place — a fixed column becomes `lb == ub` and a redundant row
has its bounds freed — so nothing is renumbered and there is no postsolve pass. It
is worth a great deal on structured models (`v006c*` is solved outright;
`bin_10var_5con` drops from 24 nodes to 1) and nothing at all on the dense random
families above. That is not a shortfall in the implementation: the commercial solver's presolve
also reduces those instances to exactly their original dimensions.

Two cut families run at the root. Lifted knapsack covers are combinatorial and need
a row that reads as a knapsack; Gomory mixed-integer cuts come off the simplex
tableau and need no structure at all, which is what reaches the dense random rows
covers cannot see. Together they raise the root bound substantially where presolve
finds nothing — on `v064c064` from 72.47 to 84.83 against an optimum of 137, and on
`v064c200` from 72.13 to 82.70 against 225. Converting that bound into a smaller
tree is another matter: with most-fractional branching the node counts do not track
the bound at all. On `v064c200`, successive cut budgets give 9304, then 16918, then
6230 nodes while the bound only improves. Branching choice, not the bound,
dominated the tree's shape.

Adding GMI cuts to covers moved `v128c1000n100` from a 63% gap after 60s to proven
optimal in 6s, and cut `v256c256n100` from 2290 nodes to 250 and `v128c256n100`
from 368 to 1. Two instances regressed instead — the node counts still do not track
the bound, which improves monotonically with every cut added.

**Cuts are nonetheless off by default**, which is not where a branch-and-*cut* solver
expects to land. Every measurement above was taken under most-fractional branching and
depth-first search. Re-measured after pseudocost branching and best-bound node selection
were in place, cutting was slower on all eleven models in the target range:

| | no cuts | 3 rounds of 32 |
|---|---|---|
| `mkp_200` | **16.6s** | 48.2s |
| `v064c1000n100` | **7.9s** | 10.5s |
| `v064c200` | **1.5s** | 2.3s |
| `v081c162n009` | **0.8s** | 1.2s |
| total, nine models | **32.1s** | 67.8s |

Most of that cost turned out to be management rather than cutting. Every separated cut
was being added, permanently: near-duplicate Gomory rows off neighbouring tableau
entries each took a row, and a cut that stopped binding after a later round still rode
along in every one of `mkp_200`'s seventy thousand node LPs. Three changes address it —
candidates are ranked by **efficacy**, the distance from the relaxation optimum to the
cut's hyperplane, which unlike raw violation does not change when a row is rescaled;
near-parallel candidates are skipped, since the second of two cuts removing the same
region costs a row and buys nothing; and cuts are **aged** out once they have sat slack
for two resolves, with a final purge of everything still slack before the tree opens.
That last one is free by construction: a row inactive at the optimum of a convex program
is not holding the bound up, so dropping it leaves the same point optimal. The node
counts confirm it, coming out identical across the purge to the last node.

Together they cost `mkp_200` four cuts instead of dozens and take it from 48.2s to 10.1s,
which is now *faster* than not cutting at all. Across the suite cutting went from a 2.1x
net loss to roughly 1.15x, and the earlier failure on `markshare_4_0` is gone. It is
still a loss, so the default does not change — but it is now close enough that the
remaining gap is worth naming precisely.

That gap is not the bound. Cutting takes `v064c200`'s root from 72.1 to 92.4 against an
optimum of 225, closing 28% of what a commercial solver closes 62% of — and the tree
comes out at 2716 nodes against 2690 without cuts. A materially better root bound buys
*no* smaller tree. The most likely reading is that these cuts are tight at the root
vertex and go slack as soon as branching fixes variables, so they stop binding exactly
where the pruning would have to happen. If that is right, the way forward is separation
at nodes rather than more of it at the root, and it needs the cut pool to track which
subtree each cut was derived under.

The cuts still work: they take the `v064c200` root bound from 72.1 to 95.4. They just no
longer pay for themselves. Separation is expensive, the cuts come out dense enough to slow
every subsequent LP, and best-bound selection had already collected most of what a better
bound was worth — the earlier wins were largely a better bound rescuing a bad search order,
and the search order is no longer bad. On `mkp_200` cutting even *raises* the node count,
72150 to 91346. On MIPLIB's `markshare_4_0` it is the difference between proving optimality
in 21s and not proving it within 60s at all.

`--cut-rounds N` turns them back on, and they are worth turning on for knapsack-structured
models like `mkp_200`, where they now win outright.

### Cutting at nodes instead

The measurement that explains the root-cut result also says what to do about it. Counting
how many root cuts are still binding as the tree descends:

| instance | depth 1 | depth 3 | depth 10 |
|---|---|---|---|
| `v064c200` | 36% | 12% | 1% |
| `v081c162n009` | 50% | 25% | 0% |
| `mkp_200` | 50% | 17% | 0% |
| `v256c256n100` | 50% | 25% | 16% |

Binding roughly halves per level. Weighted by where the nodes actually are, these rows are
carried through about 99% of the tree and bind in about 2% of it -- 15 of `mkp_200`'s 73436
nodes are at depth three or less. A root cut is a shallow-depth object, and the shallow
depths are a rounding error in the node count.

Cuts read off a *node's* tableau bind where they were made. They are valid only for that
node's subtree, so they never enter the shared model: the grown LP lives for one solve and
only the bound escapes, which is valid everywhere below that node and so is safe to prune
and order children with. That does what root cutting could not:

| | nodes, no local cuts | every 10th node | every node |
|---|---|---|---|
| `v064c200` | 2690 | 2116 | **1136** |
| `v256c256n100` | 288 | 214 | **86** |
| `v064c1000n100` | 1106 | 786 | **458** |
| `mkp_200` | 72150 | 63896 | **40200** |
| total time | 24.75s | **22.98s** | 36.50s |

Separating at every node halves trees but does not pay for itself; every tenth node is a
7% win outright, and is the default. `--local-cut-frequency N` sets it.

Growing an LP by a row would normally mean refactorizing it, which at every separating node
is most of the cost. It does not have to: a basis grown by *k* cuts is block lower
triangular against the one already factorized,

```text
    B' = [ B    0 ]        S = -I, the appended rows' own logicals
         [ R_B  S ]
```

so `B'^-1` is the existing `B^-1` plus a sparse rank-*k* correction -- FTRAN substitutes
forward through the block, BTRAN transposes it to upper triangular and mirrors. The one
trap is that the correction needs `B^-1` and not `LU^-1`, so the extension wraps the whole
base operator including its etas, and pivots taken afterwards need their own layer above
it; collapsing those two layers is silent until the first post-extension pivot. Measured
against the same binary with the reuse forced to miss, the saving tracks row count, which
is the signature to expect when what you have removed is an `O(m * fill)` factorization:

| rows | instance | reuse | refactorize | |
|---|---|---|---|---|
| 30 | `mkp_200` | 22.21s | 23.84s | 1.07x |
| 200 | `v064c200` | 1.51s | 1.60s | 1.07x |
| 1000 | `v128c1000n100` | 4.07s | 5.42s | 1.33x |
| 1000 | `v064c1000n100` | 6.64s | 9.83s | **1.48x** |
| | total, eight models | 35.90s | 42.29s | 1.18x |

Two follow-ups were measured and rejected, both of which sounded better than they were.

**Capping the cut count.** Since `mkp_200` wins with three cuts and everything else loses
with twenty, an obvious reading is that cuts pay when they are cheap to carry, so a hard
cap should make them pay generally. It does not: swept at 2, 3, 5, 8, 16 and 32 cuts per
round, no cap beats not cutting at all (38.3s against 44.7s at the best cap), and the
current default of 8 is already the best setting on the instance cutting wins. The reading
was wrong about the mechanism, not just the number. A cap is not a cost knob holding the
cuts fixed: keeping fewer rows in the first round changes the vertex the second round
separates from, so it changes which cuts get derived. Capping at 2 yields *five* cuts on
`mkp_200` where capping at 8 yields four -- a different and worse set, not a cheaper subset.

**Orthogonality across rounds.** Selection compares candidates only against others from
the same round, so a second-round cut nearly parallel to one already in the model passes a
filter built to catch exactly that -- and since each round separates from the vertex the
previous round produced, near-copies are the common case. Checking candidates against the
held cuts too does work mechanically, dropping `v064c200` from 22 rows to 15 and
`v064c1000n100` from 19 to 12. It is still 7% slower over the suite, with three instances
regressing against two gains: `v128c1000n100` goes 4.0s to 5.9s and `v256c256n100` 0.31s
to 0.43s. Two cuts sharing most of a direction still differ in the part they do not share,
and at a 0.1 orthogonality bar that remainder is carrying more bound than the duplication
is costing.

Pseudocost branching addresses that, scoring a column by the objective degradation
it has actually caused rather than by how fractional it looks. It is a large win on
the hardest instances and roughly neutral elsewhere — `v081c162n018` drops from
14602 nodes to 3828, `v081c162n009` from 87138 to 74570, while `v048c048` and
`v064c200` move slightly the wrong way.

Strong branching is implemented alongside it and is **off by default**, though the
reason changed once best-bound node selection landed.

Under depth-first selection it was catastrophic — `v128c256n100` went from 10 nodes
to 917. The suspicion at the time was dual degeneracy making the probes
uninformative. That was measured and is false: instrumenting the probe outcomes
shows every probe returning a finite, strictly positive degradation, none of them
degenerate, truncated, or infeasible.

The real cause was the interaction with depth-first search. Under a plunge a
branching decision compounds — the search descends into whichever subtree the rule
picked and stays there — so any change to the rule swings the tree enormously.
Best-bound selection is self-correcting, since the next node comes from the global
pool regardless of the last decision. Re-measured under best-bound, strong branching
*reduces* node counts on seven of eight instances, by 10% to 32%.

It stays off because those savings do not pay for the probes: two extra LPs per
candidate cost more than 10-30% fewer nodes saves, on five of eight instances. It is
worth enabling on large models — `v128c1000n100` goes from 13.3s to 9.5s at a budget
of 100 — which is what `strong_branching_budget` is for.

## Parallelism

The tree search runs across worker threads sharing one node pool and one incumbent.
Presolve, cut generation and the root heuristics all run once, before any thread is
spawned — only the tree is parallel. The CLI uses the machine's parallelism by
default; the library defaults to one thread so its behaviour is predictable.

Node counts vary between runs, because which node a worker takes depends on timing.
The *answer* does not: every bound and cut is globally valid and every worker prunes
against the shared incumbent, so the proven optimum is the same however the work is
divided. That is asserted across every sample at 2, 4 and 8 threads.

| | 1 thread | 4 threads | 16 threads |
|---|---:|---:|---:|
| `mkp_200` | 8.3s | 3.3s | 2.6s |
| `mkp_500` gap after 120s | 3.07% | 3.03% | 0.30% |
| `mkp_500` simplex iterations | 834k | 3.27M | 6.70M |

Throughput scales about 8x on 16 threads. The shortfall against linear is inherent
rather than incidental: workers expand nodes that a serial search would have pruned
against an incumbent it had already found, so total node count rises with thread
count even as wall-clock time falls. `mkp_200` goes from 19522 nodes to 49432.

## Node selection

The open node set is a depth-first plunge stack plus a best-bound pool, and by
default the plunge length is zero — the search is pure best-bound.

That default was measured, not assumed. The textbook argument for plunging is that
depth-first reaches incumbents sooner, and a child's LP re-solves in a few pivots
from its parent's basis. But the primal heuristics already supply incumbents, so
the plunge buys little and costs bound progress. Switching from pure depth-first:

| instance | depth-first | best-bound |
|---|---:|---:|
| `v081c162n018` | 13570 nodes, 7.9s | 302 nodes, 0.5s |
| `v081c162n009` | 20584 nodes, 12.7s | 1286 nodes, 1.6s |
| `v064c200` | 17472 nodes, 14.8s | 2846 nodes, 2.9s |
| `v064c1000n100` | 77% gap after 60s | solved in 14s |

The trade is real, in two places. Best-bound finds incumbents later, so a run that
hits its time limit reports a worse one — `mkp_500` ends at a 3% gap rather than 1%.
And it holds every unexplored node in memory where plunging keeps the open set to
roughly the tree depth, so `plunge_limit` is there to raise when that binds.

## The basis factorization

The basis is held as a sparse LU with Markowitz pivoting, plus a product-form eta
file for the per-pivot updates. It replaced a dense explicit inverse, which cost
`O(m^2)` per solve and `O(m^3)` to rebuild — at `m = 1000` that was 0.2 seconds per
branch-and-bound node.

On the 1000-row models, in a fixed 100-second budget:

| | dense inverse | sparse LU |
|---|---:|---:|
| `v064c1000n100` nodes | 1,370 | 13,451 |
| `v064c1000n100` simplex iterations | 6,594 | 49,245 |
| `v064c1000n100` remaining gap | 93.7% | 51.6% |
| `v064c1000n020` simplex iterations | 10,486 | 59,926 |

That is roughly 7.5x the simplex throughput, which narrows the per-iteration gap
against the commercial solver on these models from about 53x to about 9x. Smaller models gain
between 1.2x and 2.4x in wall clock.

Forrest-Tomlin updates would shorten the eta file further, and the `Basis` interface
was built so that swap needs no change to the simplex — but measurement says not to
bother. Sweeping the refactorization interval over an 80x range moves solve time by
3-7%, so neither replaying the eta file nor rebuilding the factors is where the time
goes.

Where it went was rebuilding factors that were already correct. A warm re-solve
performing *zero* pivots cost 9.1ms on a 1000-row model, against 10.9ms for a real
node doing eleven — 84% of node time was setup. A child's basis is identical to its
parent's, since bounds do not enter the basis matrix, so the factors are reusable
verbatim. Each LP now keeps its recent factorizations:

| | before | after |
|---|---:|---:|
| `v128c1000n100` | 9.2s | 4.2s |
| `v064c1000n100` | 10.9s | 9.2s |
| `v256c256n100` | 0.48s | 0.33s |

The cache holds many entries rather than one because best-bound node selection does
not visit the tree in an order that keeps a single entry warm — measured at an 8-20%
hit rate for one entry. What recurs is siblings, which share a parent's basis.

## Primal heuristics

Branch and bound cannot prune anything until it holds a feasible solution, so
finding one early matters independently of the bound. Three are tried, cheapest
first: rounding the relaxation, diving, and a feasibility pump.

Diving turned out to be the wrong tool for this instance family and **fails on
every model in it**. It commits to a rounding and re-solves a smaller LP each step,
so where the feasible set is sparse it walks into infeasibility with no way back —
adding a one-level backtrack was not enough to save it. The feasibility pump never
fixes anything, alternating instead between rounding and re-optimizing the original
constraint set under a distance objective, so its LP stays feasible throughout. It
finds solutions on five of the six models where diving finds none.

The solutions are poor — 2791 against an optimum of 137 on `v064c064` — but an
incumbent of any quality switches pruning on. `v064c200` drops from 9916 nodes and
5.7s to 3388 nodes and 2.25s. `v064c1000n020` remains unsolved, with no feasible
point found by anything, the commercial solver included.

In-tree attempts are scheduled adaptively rather than on a fixed cadence: the
interval doubles after each attempt that finds nothing and snaps back to the base
after one that succeeds. A fixed cadence is wrong in both directions, and diving
fails on whole instance families rather than the occasional node, so the wasted
attempts were measurable — running unconditionally cost `v064c1000n100` its
incumbent quality. Backing off leaves the search tree identical and removes the
overhead: `v128c1000n100` 13.3s to 9.8s, `v081c162n009` 1.7s to 1.4s, at unchanged
node counts.

## Layout

| Path | Contents |
|---|---|
| `crates/ripsolve` | The solver library |
| `crates/ripsolve-cli` | The `ripsolve` command-line application |
| `samples/` | Example models in LP and MPS format |
| `crates/ripsolve-py` | The Python extension module |
| `python/` | Python build script and test suite |
| `bench/` | Benchmark and fixture-refresh tooling |

Inside the library, `lp` holds the simplex: `lp::basis` is the basis inverse and
the two solves against it, `lp::simplex` the bounded-variable primal method that
drives them.

## Building

```sh
cargo build --release
cargo test
```

Builds on stable Rust. The test suite is self-contained — no external solver is
needed to run it.

## Usage

```sh
ripsolve info  samples/v064c064.lp
ripsolve relax samples/v064c064.lp
ripsolve gen --kind knapsack --cols 60 --rows 30 --seed 42 -o hard.lp
```

## Python

A `gurobipy`-shaped interface, so a gurobipy script that stays inside the supported
feature set runs unchanged after swapping the import:

```python
import ripsolve as gp
from ripsolve import GRB

m = gp.Model("knapsack")
x = m.addVars(n, vtype=GRB.BINARY, name="x")
m.setObjective(gp.quicksum(value[j] * x[j] for j in range(n)), GRB.MAXIMIZE)
m.addConstr(gp.quicksum(weight[j] * x[j] for j in range(n)) <= capacity)
m.optimize()
```

Build it with `python/build.sh` and see `python/README.md` for what is and is not
supported. The test suite runs the same model through both solvers and requires the
answers to match.

## Testing against a reference solver

Correctness is checked against a leading commercial solver, used strictly as a **test oracle**: nothing
in `ripsolve` links against it, and `cargo test` does not require it.
`bench/refresh_fixtures.py` records reference LP-relaxation and MIP-optimal values
into `crates/ripsolve/tests/fixtures/reference.json`, and the tests read only that
file. Each entry carries a digest of its instance, so a change to the generator
fails the suite as a stale fixture rather than silently pairing new instances with
old reference values.

Refresh them (needs `gurobipy`) with:

```sh
cargo build --release && python3 bench/refresh_fixtures.py
```

## Licence

MIT — see [LICENSE](LICENSE).
