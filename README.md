# bipper

A branch-and-cut solver for **binary integer programs**, in pure Rust.

`bipper` is an LP-relaxation-based branch-and-cut solver: every node solves a
bounded LP relaxation with the simplex method, and the resulting dual bound —
strengthened by presolve and cutting planes — is what prunes the search. Every
variable is binary; general integer and continuous variables are out of scope.

Status: **early**. The model layer, LP/MPS readers, instance generator, and test
oracle are in place. The simplex and the search are not yet written.

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

## Layout

| Path | Contents |
|---|---|
| `crates/bipper` | The solver library |
| `crates/bipper-cli` | The `bipper` command-line application |
| `samples/` | Example models in LP and MPS format |
| `bench/` | Benchmark and fixture-refresh tooling |

## Building

```sh
cargo build --release
cargo test
```

Builds on stable Rust. The test suite is self-contained — no external solver is
needed to run it.

## Usage

```sh
bipper info samples/v064c064.lp
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
