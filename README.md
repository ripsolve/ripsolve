# bipper

A branch-and-cut solver for **binary integer programs**, in pure Rust.

`bipper` is an LP-relaxation-based branch-and-cut solver: every node solves a
bounded LP relaxation with the simplex method, and the resulting dual bound —
strengthened by presolve and cutting planes — is what prunes the search. Every
variable is binary; general integer and continuous variables are out of scope.

Status: **early, but solving**. The model layer, LP/MPS readers, instance
generator, test oracle, LP relaxation solver, presolve, and a depth-first
branch-and-bound search are in place. Cutting planes, better branching, and primal
heuristics are not — those are what would make it competitive rather than merely
correct.

Both the relaxation and the integer optimum are checked against an independent
solver on every bundled sample and generated instance.

## Why not pure enumeration

The predecessor to this project implemented Balas' additive algorithm — implicit
enumeration with a feasibility-based bound. Measured on a dense 60-variable
instance, that approach explored 10.3 billion nodes in 54.7s across 16 threads.
Gurobi solved the same instance in 125 nodes and 0.15s.

The enumeration engine was not slow; at 188M nodes/sec it processed nodes roughly
1000x faster than Gurobi does. It lost because it visited 82 million times more of
them. No amount of SIMD or threading closes a gap of that shape — only a stronger
bound does, and that means solving an LP relaxation at each node. That measurement
is the reason this project exists and why the simplex is its centrepiece.

On that same 60-variable instance, `bipper` now proves optimality in **2204 nodes
and 0.10s** — 4.7 million times fewer nodes than the enumeration approach, and 543x
faster in wall clock. On a family of dense random instances with 30 rows:

| variables | enumeration | bipper | Gurobi |
|---:|---:|---:|---:|
| 60 | 54.7s | 0.10s | 0.15s |
| 64 | >120s | 0.20s | 0.24s |
| 80 | — | 0.52s | 0.23s |
| 200 | — | 19.0s | 0.87s |
| 500 | — | >400s | 82s |

The gap that remains against Gurobi is roughly 20x in nodes, and it is exactly the
work that has not been done yet: cuts and pseudocost branching.

Presolve reduces in place — a fixed column becomes `lb == ub` and a redundant row
has its bounds freed — so nothing is renumbered and there is no postsolve pass. It
is worth a great deal on structured models (`v006c*` is solved outright;
`bin_10var_5con` drops from 24 nodes to 1) and nothing at all on the dense random
families above. That is not a shortfall in the implementation: Gurobi's presolve
also reduces those instances to exactly their original dimensions. On that family
the entire advantage is cutting planes.

## Layout

| Path | Contents |
|---|---|
| `crates/bipper` | The solver library |
| `crates/bipper-cli` | The `bipper` command-line application |
| `samples/` | Example models in LP and MPS format |
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
bipper info  samples/v064c064.lp
bipper relax samples/v064c064.lp
bipper gen --kind knapsack --cols 60 --rows 30 --seed 42 -o hard.lp
```

## Testing against a reference solver

Correctness is checked against Gurobi, used strictly as a **test oracle**: nothing
in `bipper` links against it, and `cargo test` does not require it.
`bench/refresh_fixtures.py` records reference LP-relaxation and MIP-optimal values
into `crates/bipper/tests/fixtures/reference.json`, and the tests read only that
file. Each entry carries a digest of its instance, so a change to the generator
fails the suite as a stale fixture rather than silently pairing new instances with
old reference values.

Refresh them (needs `gurobipy`) with:

```sh
cargo build --release && python3 bench/refresh_fixtures.py
```

## Licence

MIT.
