# ripsolve-cli

The command-line application for [`ripsolve`](https://crates.io/crates/ripsolve), a
branch-and-cut solver for mixed-integer programs in pure Rust.

```sh
cargo install ripsolve-cli
```

## Usage

```sh
ripsolve solve model.lp                 # solve to proven optimality
ripsolve solve model.mps --time-limit 60 --gap 0.01
ripsolve info  model.lp                 # dimensions and column types
ripsolve relax model.lp                 # LP relaxation bound only
ripsolve gen --kind knapsack --cols 60 --rows 30 --seed 42 -o hard.lp
```

Both LP and MPS files are read, chosen by extension. Output names the objective, the
status, and what the search cost:

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
| `-v` | Print the value of every column |

A run that hits its time limit reports the best solution found and the remaining gap
rather than failing.

## Scope

Built to be fast on small and medium models, up to roughly a thousand rows. Above that
range it is not competitive; reach for HiGHS or SCIP there.

See the [repository](https://github.com/ripsolve/ripsolve) for benchmarks, the Rust
library, and a `gurobipy`-shaped Python module.

## Licence

MIT.
