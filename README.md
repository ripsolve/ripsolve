# ripsolve

A branch-and-cut solver for mixed-integer programs, in pure Rust.

[![crates.io](https://img.shields.io/crates/v/ripsolve.svg)](https://crates.io/crates/ripsolve)
[![docs.rs](https://docs.rs/ripsolve/badge.svg)](https://docs.rs/ripsolve)

Columns may be binary, general integer, or continuous. Every node solves a bounded LP
relaxation with the simplex method, and the resulting dual bound, strengthened by
presolve and cutting planes, is what prunes the search.

There is a command-line application, a Rust library, and a Python module shaped like
`gurobipy`. Nothing links against another solver.

## Is it the right tool?

`ripsolve` is built to be fast on small and medium models, up to roughly a thousand
rows. In that range it is quicker than HiGHS, SCIP and CBC on most of the benchmark
set, and competitive with Gurobi.

Outside that shape it is not competitive, and the distance is not small. On a seeded
random sample of MIPLIB 2017's `easy` list, sixteen threads and a sixty-second limit for
every solver, it closed **three of the sixteen instances Gurobi closes**. Six of the
thirteen it missed are ones Gurobi finishes in under three seconds. Reach for HiGHS or
SCIP there.

`bench/miplib_sample.py` reproduces that. It runs the reference first and skips any
instance the reference cannot close within the same budget, because a timeout on both
sides measures the instance rather than this solver, and it caches the reference's
answers so a re-run costs only ripsolve's time.

Implemented: LP and MPS readers, presolve, knapsack cover cuts, Gomory mixed-integer
cuts, node-local cut separation, pseudocost branching, primal heuristics, a sparse LU
factorization, best-bound search, and parallel tree search.

Not implemented: quadratic terms, SOS constraints, callbacks, lazy constraints,
multi-objective.

## Install

Command-line application:

```sh
cargo install ripsolve-cli
```

Rust library:

```sh
cargo add ripsolve
```

Python module, which needs a Python with development headers:

```sh
git clone https://github.com/ripsolve/ripsolve
cd ripsolve/python && ./build.sh
```

That writes `ripsolve.so` next to the script. Put it on your `PYTHONPATH`.

## Command line

```sh
ripsolve solve model.lp                 # solve to proven optimality
ripsolve solve model.mps --time-limit 60 --gap 0.01
ripsolve info  model.lp                 # dimensions and column types
ripsolve relax model.lp                 # LP relaxation bound only
ripsolve gen --kind knapsack --cols 60 --rows 30 --seed 42 -o hard.lp
```

Output names the objective, the status, and what the search cost:

```
objective: 225
status:    optimal
presolve:  0 columns fixed, 0 rows removed, 0 coefficients tightened
heuristic: 1 incumbents
2116 nodes, 9284 simplex iterations, 1.21s
```

Useful flags for `solve`:

| Flag | Effect |
|---|---|
| `--time-limit <s>` | Stop after this many seconds and report the gap |
| `--gap <g>` | Stop once the relative gap reaches `g` |
| `--threads <n>` | Worker threads. Defaults to the machine's parallelism |
| `--local-cut-frequency <n>` | Separate cuts at one node in every `n`. Default 10, `0` disables |
| `--cut-rounds <n>` | Rounds of root cut separation. Default 0 |
| `--no-presolve` | Skip presolve |
| `--values <none\|nonzero\|all>` | Print column values, one per line |
| `--solution <PATH>` | Write the solution to a file as `name value` lines |
| `-v` | Shorthand for `--values nonzero` |

A run that hits its time limit reports the best solution found and the remaining gap
rather than failing.

## Python

The interface is shaped like `gurobipy`, so a script that stays inside the supported
feature set runs unchanged after swapping the import.

```python
import ripsolve as gp
from ripsolve import GRB

value  = [12, 9, 7, 5, 3]
weight = [ 6, 5, 4, 3, 2]

m = gp.Model("knapsack")
x = m.addVars(len(value), vtype=GRB.BINARY, name="x")
m.setObjective(gp.quicksum(value[j] * x[j] for j in range(len(value))), GRB.MAXIMIZE)
m.addConstr(gp.quicksum(weight[j] * x[j] for j in range(len(value))) <= 10)
m.optimize()

print(m.ObjVal)                              # 19.0
print([j for j in range(len(value)) if x[j].X > 0.5])
```

Mixed-integer models use the same `vtype` values as gurobipy:

```python
m = gp.Model()
b = m.addVar(vtype=GRB.BINARY, name="b")
n = m.addVar(vtype=GRB.INTEGER, lb=0, ub=10, name="n")
c = m.addVar(vtype=GRB.CONTINUOUS, lb=0.0, name="c")

m.addConstr(2 * b + n + 0.5 * c <= 12)
m.setObjective(3 * b + 2 * n + c, GRB.MAXIMIZE)
m.setParam("TimeLimit", 30)
m.setParam("MIPGap", 0.01)
m.optimize()

if m.Status == GRB.OPTIMAL:
    print(m.ObjVal, m.NodeCount, m.Runtime)
```

Reading a model from a file:

```python
m = gp.read("model.mps")
m.optimize()
```

### Supported surface

| | |
|---|---|
| Model | `Model(name)`, `optimize`, `getVars`, `write`, `setParam`, `read` |
| Variables | `addVar`, `addVars`, `.X`, `.VarName`, `.VType`, `.Obj`, `.LB`, `.UB` |
| Expressions | `+ - *`, unary `-`, `<= >= ==`, `quicksum` |
| Constraints | `addConstr`, `addConstrs`, `.ConstrName` |
| Objective | `setObjective(expr, GRB.MINIMIZE / GRB.MAXIMIZE)` |
| Attributes | `ObjVal`, `ObjBound`, `Status`, `SolCount`, `NodeCount`, `Runtime`, `MIPGap`, `NumVars`, `NumConstrs`, `ModelName`, `ModelSense` |
| Parameters | `TimeLimit`, `Threads`, `MIPGap`, `OutputFlag` |

Two deliberate differences from gurobipy's behaviour. `vtype` defaults to continuous,
as gurobipy's does, even though binary would suit this solver's history: a ported
script calling `addVar()` should not silently get a different model. Unknown parameter
names raise `KeyError` rather than being ignored, so a misspelling fails loudly.

`optimize()` releases the GIL, so a solve does not block other Python threads.

See [`python/README.md`](python/README.md) for the full list.

## Rust library

```rust
use ripsolve::{Problem, search};
use std::path::Path;

let problem = Problem::from_file(Path::new("model.lp"))?;
problem.validate()?;

let solution = search::solve(&problem, search::Options::default());
println!("{:?} {:?}", solution.status, solution.objective);
```

`Builder` assembles a model a column and a row at a time. Objective coefficients are
written in the sense you ask for. Maximizing `3b + 2n` subject to `2b + n <= 12`:

```rust
use ripsolve::model::{Builder, RowSense, Sense};
use ripsolve::search;

let mut model = Builder::new(Sense::Maximize).named("example");
let b = model.binary("b");
let n = model.integer("n", 0.0, 10.0);
model.objective(&[(b, 3.0), (n, 2.0)]);
model.row(&[(b, 2.0), (n, 1.0)], RowSense::Le, 12.0);

let problem = model.build();
problem.validate()?;

let solution = search::solve(&problem, search::Options::default());
assert_eq!(solution.objective, Some(23.0));
```

`continuous` adds a column with no integrality requirement, and `range` adds a row
bounded on both sides. A binary column is an integer column bounded to `[0, 1]`;
nothing treats that as a distinct case, because branching splits a range and
degenerates to fixing at 0 or 1.

`Problem` is also a plain struct with public fields, so it can be filled in directly
when translating a model from somewhere else.

`search::Options` carries the tuning knobs, including `threads`, `time_limit`,
`gap_tolerance`, `local_cut_frequency` and `cut_rounds`. The library defaults to one
thread so its behaviour is predictable; the CLI defaults to the machine's
parallelism.

Full API documentation is on [docs.rs](https://docs.rs/ripsolve).

## Benchmarks

Every solver on one thread with the same 60-second limit, timing only the solve.
`bench/open_compare.py` reproduces this, and fails if any two solvers that both claim
optimality disagree on the optimum, so it is a differential test as much as a
benchmark.

| instance | rows | ripsolve | HiGHS | SCIP | CBC |
|---|---:|---:|---:|---:|---:|
| `v032c032` | 32 | **0.0s** | 0.0s | 0.0s | 0.0s |
| `v048c048` | 48 | **0.0s** | 0.1s | 0.2s | 0.0s |
| `v048c128` | 128 | **0.0s** | 0.1s | 0.1s | 0.1s |
| `v064c064` | 64 | **0.0s** | 0.2s | 0.5s | 0.1s |
| `v081c162n018` | 162 | **0.2s** | 0.5s | 0.9s | 0.5s |
| `v128c256n100` | 256 | **0.0s** | 0.1s | 0.1s | none |
| `v256c256n100` | 256 | **0.3s** | 0.7s | 1.0s | none |
| `v081c162n009` | 162 | **1.0s** | 1.2s | 1.7s | 1.6s |
| `mkp_200` | 30 | **11.5s** | 12.6s | 22.8s | 20.3s |
| `v064c200` | 200 | 1.7s | 1.8s | **1.4s** | 2.6s |
| `v128c1000n100` | 1000 | 5.6s | 1.8s | **1.8s** | none |
| `v064c1000n100` | 1000 | 14.0s | 3.2s | **2.2s** | 6.6s |

Fastest or tied on nine of twelve, and every proven optimum agrees. The models it loses
are the widest: both thousand-row instances, and `v064c200` by a fifth of a second. CBC
finds nothing at all on three.

Against Gurobi on the same models and the same limit, matching thread counts so the
comparison is of algorithms rather than of core usage:

| | single-threaded | sixteen threads |
|---|---|---|
| within 1x | six of twelve | eight of twelve |
| 2x to 3x | `v081c162n009`, `v256c256n100`, `v064c200` | `v048c128`, `v128c1000n100`, `mkp_200`, `v064c1000n100` |
| worst | `v064c1000n100` and `mkp_200` at 7x | `v064c1000n100` at 3x |

Every objective it proves matches. Two models sit outside: `v064c1000n020` yields no
feasible point to either solver, and `mkp_500` reaches the same 6191 Gurobi does at one
thread, but does not close it at sixteen, where the extra workers find more incumbents
and prove less.

All three tables come from one sitting on one machine, which matters more than it
sounds: absolute times here drifted by a factor of two across a day of benchmarking, so
figures are comparable within a table and not against an older copy of this file. Every
comparison runs both solvers in the same session for that reason.

Gurobi is run under an academic licence. It is named here rather than anonymised, and
the figures for it are reported the same way as everyone else's: same machine, same
limit, same thread count, one sitting. `bench/compare.py` and `bench/binary_bench.py`
reproduce these for anyone holding a licence. Nothing is withheld.

### MIPLIB, against the open-source field

The tables above are generated models. This one is MIPLIB: the 140 instances of its
`easy` set whose every column is binary, at a 60-second limit and sixteen threads.
Membership is a property of the model rather than of what any solver manages on it, so
the set cannot flatter the solver being measured. `bench/binary_bench.py` reproduces it.

Of the 140, **57 are closed by at least one open-source solver**, and those are the rows
below; the other 83 defeat all four. Counting over the whole 140:

| | closes |
|---|---:|
| Gurobi | 94 |
| HiGHS | 46 |
| SCIP | 45 |
| **ripsolve** | **34** |
| CBC | 30 |

So this solver is ahead of CBC and twelve behind HiGHS, and the honest way to read that
is coverage rather than speed: on the instances it does close it is usually quick, and
the gap is instances it does not close at all. Six of its 34 are ones HiGHS does not
close -- `decomp2`, `disctom`, `eil33-2`, `neos-4754521-awarau`, `neos-787933`, `nw04` --
and two are ones SCIP does not, so the sets are not nested.

Bold marks the fastest of the four open-source solvers on that row. `--` is a solver that
did not close the instance inside the limit; `inf` marks an instance proven infeasible,
which is a complete answer and counts as closed.

| instance | rows | cols | ripsolve | HiGHS | SCIP | CBC | Gurobi |
|---|---:|---:|---:|---:|---:|---:|---:|
| `supportcase14` | 234 | 304 | **0.0s** | 0.1s | 0.2s | 0.2s | 0.0s |
| `supportcase16` | 130 | 319 | **0.0s** | 0.1s | 0.2s | 0.2s | 0.0s |
| `app2-2` | 335 | 1226 | **0.1s** | 0.6s | 1.5s | -- | 0.0s |
| `f2gap40400` | 40 | 400 | 0.1s | **0.0s** | 0.1s | 0.4s | 0.0s |
| `stein9inf` | 14 | 9 | 0.2s inf | 0.0s inf | 0.0s inf | **0.0s inf** | 0.0s inf |
| `f2gap401600` | 40 | 1600 | 0.3s | 0.2s | **0.1s** | 0.4s | 0.1s |
| `neos-1516309` | 489 | 4500 | 0.5s | **0.3s** | 3.2s | 2.7s | 0.1s |
| `stein15inf` | 37 | 15 | 0.5s inf | 0.1s inf | **0.0s inf** | 0.1s inf | 0.0s inf |
| `f2gap201600` | 20 | 1600 | 0.6s | 1.4s | **0.2s** | 0.3s | 0.1s |
| `mod010` | 146 | 2655 | 0.6s | 0.6s | 0.6s | **0.4s** | 0.2s |
| `neos-3437289-erdre` | 7358 | 12162 | 0.6s | 1.1s | **0.0s** | -- | 0.0s |
| `p0201` | 133 | 201 | **0.6s** | 1.2s | 0.9s | 2.3s | 0.4s |
| `app2-1` | 1038 | 3283 | 1.2s | 3.1s | 1.1s | **1.1s** | 0.1s |
| `air03` | 124 | 10757 | 1.6s | 2.6s | 1.6s | **1.1s** | 0.3s |
| `neos-3148108-pahi` | 154774 | 12241 | 1.6s | 1.6s | **0.2s** | 7.7s | 0.0s |
| `stein45inf` | 332 | 45 | 3.2s inf | 1.3s inf | **1.0s inf** | 1.5s inf | 0.2s inf |
| `mod008inf` | 7 | 319 | 3.9s inf | 0.1s inf | 0.1s inf | **0.1s inf** | 0.0s inf |
| `neos-4754521-awarau` | 255386 | 41680 | **6.9s** | -- | 52.6s | -- | 0.1s |
| `f2gap801600` | 80 | 1600 | 7.0s | **0.1s** | 0.2s | 0.7s | 0.1s |
| `neos-3226448-wkra` | 2208 | 10012 | 8.8s | **0.2s** | -- | -- | 25.5s |
| `disctom` | 399 | 10000 | 15.9s | -- | 6.2s | **1.5s** | 0.2s |
| `n2seq36f` | 285 | 8100 | 18.4s | **3.6s** | 10.1s | -- | 0.2s |
| `ex9` | 40962 | 10404 | 25.0s | **3.5s** | 23.1s | -- | 2.7s |
| `mitre` | 2054 | 10724 | 27.4s | 1.4s | 10.2s | **1.2s** | 0.1s |
| `cap6000` | 2176 | 6000 | 29.1s | 2.8s | 5.5s | **1.5s** | 0.3s |
| `eil33-2` | 32 | 4516 | **31.4s** | -- | -- | -- | 2.1s |
| `decomp1` | 8357 | 10735 | 32.2s | 18.1s | **1.6s** | 20.0s | 0.3s |
| `neos-1599274` | 1237 | 4500 | 33.5s | **0.6s** | 1.5s | 1.1s | 0.2s |
| `acc-tight2` | 2520 | 1620 | 37.2s | **0.1s** | 16.6s | 3.2s | 0.6s |
| `neos-913984` | 1076 | 76000 | 37.3s | **4.3s** | 4.8s | 5.1s | 0.2s |
| `neos-787933` | 1897 | 236376 | 39.3s | -- | **1.7s** | -- | 0.4s |
| `decomp2` | 10765 | 14387 | 42.3s | -- | **1.8s** | 18.9s | 0.3s |
| `nw04` | 36 | 87482 | 49.2s | -- | 40.1s | **15.4s** | 2.2s |
| `irp` | 39 | 20315 | 53.7s | 23.7s | 18.8s | **10.5s** | 1.0s |
| `ab71-20-100` | 6380 | 6689 | -- | **41.6s** | 58.0s | -- | 1.9s |
| `acc-tight4` | 3285 | 1620 | -- | **24.2s** | -- | -- | 5.2s |
| `acc-tight5` | 3052 | 1339 | -- | 49.7s | **19.5s** | -- | 6.3s |
| `air04` | 823 | 8904 | -- | **55.0s** | -- | -- | 8.5s |
| `air05` | 426 | 7195 | -- | 34.6s | 36.6s | **13.1s** | 3.3s |
| `bley_xl1` | 175620 | 5831 | -- | **17.7s** | -- | -- | 2.9s |
| `chromaticindex32-8` | 2111 | 2304 | -- | 26.9s | **1.0s** | -- | 0.1s |
| `cod105` | 1024 | 1024 | -- | -- | **57.2s** | -- | 42.7s |
| `ex10` | 69608 | 17680 | -- | **33.3s** | -- | -- | 18.2s |
| `graph20-20-1rand` | 5587 | 2183 | -- | -- | **6.2s** | -- | -- |
| `mine-166-5` | 8429 | 830 | -- | **52.6s** | -- | -- | 1.2s |
| `neos-1324574` | 5904 | 5256 | -- | -- | **43.5s** | -- | 12.0s |
| `neos-3045796-mogo` | 2245 | 11016 | -- | **1.7s** | 2.4s | 2.2s | 0.1s |
| `neos-5129192-manaia` | 531569 | 164229 | -- | -- | **50.9s** | -- | 14.5s |
| `neos-633273` | 21781 | 11154 | -- | 5.7s | **0.5s** | 2.6s | 0.2s |
| `neos-780889` | 73910 | 182700 | -- | **7.2s** | -- | 33.1s | 1.4s |
| `neos-820879` | 361 | 9522 | -- | **26.8s** | -- | -- | 4.4s |
| `neos-953928` | 12498 | 23305 | -- | **15.0s** | -- | -- | 1.1s |
| `neos-957143` | 2767 | 57756 | -- | **42.4s** | -- | -- | 2.2s |
| `neos-960392` | 4744 | 59376 | -- | **21.9s** | -- | -- | 4.7s |
| `neos18` | 11402 | 3312 | -- | 43.5s | **14.4s** | -- | 1.3s |
| `supportcase4` | 9492 | 3162 | -- | 23.1s | **2.1s** | -- | -- |
| `tanglegram6` | 17712 | 9182 | -- | -- | **19.2s** | -- | 0.1s |

Two cautions about reading this table. It is one run: the parallel search is not
deterministic, and an instance finishing near the limit may not finish on the next run --
`neos-3045796-mogo` closes at 38 to 41 seconds and is a coin flip, which is why it is
absent here and present in earlier baselines. And Gurobi, at the same limit and thread
count, closes 55 of these 57 and 94 of the 140 -- which is a more useful thing to know
about the state of the art than any ranking among the rest.

## Building and testing

```sh
cargo build --release
cargo test
```

Builds on stable Rust. The test suite is self-contained and needs no external solver.

Correctness is additionally checked against other solvers used strictly as oracles.
`bench/refresh_fixtures.py` records reference relaxation and optimal values into
`crates/ripsolve/tests/fixtures/reference.json`, and the tests read only that file.
Each entry carries a digest of its instance, so changing the generator fails the suite
as a stale fixture rather than silently pairing new instances with old values.
`bench/mip_fuzz.py` compares randomized mixed-integer models against a reference
solver.

## Layout

| Path | Contents |
|---|---|
| `crates/ripsolve` | The solver library |
| `crates/ripsolve-cli` | The `ripsolve` command-line application |
| `crates/ripsolve-py` | The Python extension module |
| `python/` | Python build script and test suite |
| `samples/` | Example models in LP and MPS format |
| `bench/` | Benchmark and fixture-refresh tooling |
| `docs/design-notes.md` | Why the solver is built the way it is |

## Licence

MIT, see [LICENSE](LICENSE).
