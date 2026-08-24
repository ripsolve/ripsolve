# ripsolve

A branch-and-cut solver for **mixed-integer programs**, in pure Rust.

`ripsolve` is an LP-relaxation-based branch-and-cut solver: every node solves a
bounded LP relaxation with the simplex method, and the resulting dual bound —
strengthened by presolve and cutting planes — is what prunes the search.

Columns may be binary, general integer, or continuous. A binary variable is just an
integer bounded to `[0, 1]`; the solver has no separate notion of one, because
branching splits a range (`x <= floor(v)` against `x >= ceil(v)`) and that
degenerates to fixing at 0 or 1 on its own.

Status: **early, but solving**. LP and MPS readers, presolve, lifted knapsack cover
cuts, Gomory mixed-integer cuts, pseudocost branching, primal heuristics, a sparse
LU factorization, best-bound search, parallel tree search, and a `gurobipy`-shaped
Python API are all in place. Quadratic terms, SOS constraints, callbacks and lazy
constraints are not.

Correctness is checked against a leading commercial solver as an oracle — relaxation values, integer
optima, and 400 randomized mixed-integer models — but nothing links against it and
the test suite runs without it.

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
