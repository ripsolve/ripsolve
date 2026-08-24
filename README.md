# ripsolve

A branch-and-cut solver for **binary integer programs**, in pure Rust.

`ripsolve` is an LP-relaxation-based branch-and-cut solver: every node solves a
bounded LP relaxation with the simplex method, and the resulting dual bound —
strengthened by presolve and cutting planes — is what prunes the search. Every
variable is binary; general integer and continuous variables are out of scope.

Status: **early, but solving**. The model layer, LP/MPS readers, instance
generator, test oracle, LP relaxation solver, presolve, lifted knapsack cover
cuts, pseudocost branching, and a depth-first branch-and-bound search are in place.
Primal heuristics are not, and neither is a sparse LU factorization — which the
measurements now say is the binding constraint on larger models.

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

On that same 60-variable instance, `ripsolve` now proves optimality in **2204 nodes
and 0.10s** — 4.7 million times fewer nodes than the enumeration approach, and 543x
faster in wall clock. On a family of dense random instances with 30 rows:

| variables | enumeration | ripsolve | Gurobi |
|---:|---:|---:|---:|
| 60 | 54.7s | 0.10s | 0.15s |
| 64 | >120s | 0.20s | 0.24s |
| 80 | — | 0.52s | 0.23s |
| 200 | — | 19.0s | 0.87s |
| 500 | — | >400s | 82s |

The gap that remains against Gurobi is roughly 20x in nodes, and the work that
would close it is better branching and primal heuristics.

Presolve reduces in place — a fixed column becomes `lb == ub` and a redundant row
has its bounds freed — so nothing is renumbered and there is no postsolve pass. It
is worth a great deal on structured models (`v006c*` is solved outright;
`bin_10var_5con` drops from 24 nodes to 1) and nothing at all on the dense random
families above. That is not a shortfall in the implementation: Gurobi's presolve
also reduces those instances to exactly their original dimensions.

Lifted knapsack cover cuts raise the root bound substantially where presolve finds
nothing — on `v064c064` from 72.47 to 84.83 against an optimum of 137, and on
`v064c200` from 72.13 to 82.70 against 225. Converting that bound into a smaller
tree is another matter: with most-fractional branching the node counts do not track
the bound at all. On `v064c200`, successive cut budgets give 9304, then 16918, then
6230 nodes while the bound only improves. Branching choice, not the bound,
dominated the tree's shape.

Pseudocost branching addresses that, scoring a column by the objective degradation
it has actually caused rather than by how fractional it looks. It is a large win on
the hardest instances and roughly neutral elsewhere — `v081c162n018` drops from
14602 nodes to 3828, `v081c162n009` from 87138 to 74570, while `v048c048` and
`v064c200` move slightly the wrong way.

Strong branching is implemented alongside it and is **off by default**, because
measurement said so: probing made every instance tried markedly worse, taking
`v128c256n100` from 10 nodes to 917. Two explanations were tested and ruled out.
See `branch.rs` for what is and is not yet known about why.

## Layout

| Path | Contents |
|---|---|
| `crates/ripsolve` | The solver library |
| `crates/ripsolve-cli` | The `ripsolve` command-line application |
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
ripsolve info  samples/v064c064.lp
ripsolve relax samples/v064c064.lp
ripsolve gen --kind knapsack --cols 60 --rows 30 --seed 42 -o hard.lp
```

## Testing against a reference solver

Correctness is checked against Gurobi, used strictly as a **test oracle**: nothing
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

MIT.
