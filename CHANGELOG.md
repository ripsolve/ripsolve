# Changelog

## Unreleased

### Fixed

- **Time limits held on models with very dense rows.** Three loops were quadratic in
  row length and none of them checked the clock: coefficient tightening in presolve
  recomputed a row's activity for every entry, and the conflict-driven search and the
  fixing heuristics reread every row whose slack had not changed. `pb-fit2d` at a
  60-second limit returned after 285s and now returns after 63s; `nw04` at a 10-second
  limit on one thread ran to 14s and now stops on time. The conflict-driven search also
  checks the deadline on every conflict, not every 64th.
- Separating implied-bound cuts reads each clique once rather than once per column in
  it. The cuts are the same.

### Benchmarks

- The pure binary MIPLIB set, 140 instances at sixty seconds, rerun with every answer
  checked against MIPLIB's solution file: 51 closed at sixteen threads (HiGHS 46, SCIP
  45, CBC 30, Gurobi 94) and 47 at one thread with the 0.4.0 build (HiGHS 39, SCIP 41,
  CBC 31, Gurobi 81). Optima agree with MIPLIB's to within the default 1e-4 gap. The
  results are in `docs/baselines/binary-2026-10-23.json` and the README.
- `bench/benchmark_set.py` runs MIPLIB 2017's 240 benchmark instances and checks each
  answer against MIPLIB.

## 0.4.0

Over this release the pure binary MIPLIB set used to track progress, 140 instances at
sixty seconds and sixteen threads, went from 22 closed to 36. For comparison on the
same run: HiGHS 46, SCIP 45, CBC 30.

### Changed behaviour you will notice

- **Root cut separation is on**, at 50 rounds of 64 cuts (`cut_rounds`,
  `cuts_per_round`). It doubles the time on generated knapsack and covering models
  that solve either way, and decides whether a MIPLIB instance is solved at all:
  `nexp-50-20-1-1` from a 47.3% gap to optimal in 1.3s.
- **A conflict-driven search runs beside the tree** on the pure binary part of a model
  (`cdcl`, `--no-cdcl`), on its own thread. Its points prune the tree and its proof
  ends the solve: the `easter_monster` sudoku from 2.4s to 14ms.
- **Doubleton aggregation is on** (`aggregate`, `--no-aggregate`), with a postsolve.
- **Node propagation, reduced cost fixing, restarts and the LP-free feasibility jump**
  are on by default; see the new options below.
- An LP file's variable with no `Bounds` entry is read as `[0, inf)`, as the format
  says, not as free; a bound at `1e30` or beyond reads as infinite.
- `search::Options` has new public fields, so code that spells every field out needs
  `..Options::default()`.

### Added

- `search::solve_with` and the `Separator` trait: rows supplied while the search runs,
  asked for at the root, at every node, and of every point before it may become the
  incumbent. It serves as a cut callback and as lazy constraints at once. A clique
  partitioning model whose 52000 transitivity rows are left to the separator solves in
  one call rather than in an outer loop of whole MIP solves.
- `Lp::fit`, `BasisState::without_rows` and `lu::Scratch`.
- Python: `Model.optimize(callback)` with gurobipy's callback API for lazy constraints
  and cuts: `GRB.Callback.MIPNODE` and `MIPSOL`, `cbGetSolution`, `cbGetNodeRel`,
  `cbLazy`, `cbCut` and `cbGet`, the `LazyConstraints` parameter, `model.Params`, and
  attributes set on the model for the callback to read (`model._vars`). The same
  lazy-constraint script returns Gurobi's answer under either module.
- `presolve::presolve_with`, which can leave out the reductions that fix a column by
  its cost.
- `Problem::to_lp`, the LP writer, in the core crate rather than private to the
  Python module, and `Problem::from_text`.
- Cut families: mixed-integer rounding with variable bound substitution, mod-2 cuts
  that combine rows, clique cuts from a conflict graph over literals, and rows
  aggregated through the bounds the model implies.
- Heuristics: fixing with propagation, a feasibility jump that needs no relaxation,
  the corners of the box, searches over the neighbourhood the root's reduced costs
  would fix and over the model rewritten through its implied bounds.
- Reduced cost fixing at the root and as the incumbent earns more, restarts on the
  model those fixings leave, and a stalled root relaxation restarted from a perturbed
  one.
- Optional, off by default and documented with what they were measured to be worth:
  symmetry handling by orbital branching (`symmetry`), root probing, rapid learning,
  and conflict analysis.
- CLI flags for each of them: `--no-node-propagation`, `--conflict-analysis`,
  `--root-probing`, `--rapid-learning`, `--no-aggregate`, `--no-cdcl`,
  `--jump-moves`, `--symmetry` and `--min-orbit`; `relax --dual` enters the
  relaxation through the dual method.
- Modules `aggregate`, `cdcl`, `compact`, `propagate` and `symmetry`.

### Python

- Packaged for PyPI as `ripsolve`, so `pip install ripsolve` works once a release is
  tagged: one abi3 wheel per platform for every CPython from 3.9, built, tested and
  published by CI on a version tag.
- `addLConstr`, `GRB.LESS_EQUAL`, `GREATER_EQUAL` and `EQUAL`; attribute names ignore
  case as gurobipy's do; and the model reads variable and constraint attributes as
  lists, `model.X`.

### Performance

- Refactorization is triggered by what the eta file has cost against what the last
  factorization cost, not by the eta file's size against the factors'. With the
  factorization's working memory kept between calls, root relaxations over the cached
  MIPLIB set are 25% faster by geometric mean on those taking half a second or more,
  identical in value on all 128 both finish, and five more finish inside 20 seconds.
- Dropping a slack cut at the root keeps the basis instead of re-solving cold, in the
  ageing and in the final purge.
- Warm-started dual solves price their leaving row by steepest edge
  (`Options::warm_steepest_edge`, on by default). Over the cached MIPLIB set at thirty
  seconds it closes `n9-3` and `ex9` and is 12% slower by geometric mean on the 32
  instances both settings close.
- The factorization no longer searches a column for rows it already holds, which was
  quadratic on a basis that fills in: `cod105`'s relaxation from 57.8s to 9.0s.
- The cut loop grows its model and warm-starts instead of rebuilding and solving cold,
  and its final purge is bounded.
- A stalling relaxation switches itself to devex pricing, and Bland's rule is left
  after a short spell.

### Correctness

- Fixing with propagation decides only binaries, and returned general integer columns
  at their fractional relaxation values; nothing downstream checked integrality.
  `neos-3381206-awhea` was reported optimal at 452 against an optimum of 453, and
  `neos-3024952-loue` optimal at a different wrong value on each run.
- An LP built for a round of root cutting did not inherit the flag the search beside
  the tree raises when it proves its point, so cutting ran on after the answer was
  known: `bnatt400` took 28 seconds where it takes 1.5.
- An integral node's incumbent is scored from its rounded point, as the root's is,
  rather than from the relaxation's objective.
- Failing to factorize a basis no longer proves a model infeasible, and neither does
  phase 1's arithmetic noise: `supportcase4` and `gasprod1-2` had been reported
  infeasible.
- The ratio test refuses a pivot negligible within its column, which on `neos-850681`
  took the basis inverse to entries of 1e15.
- A column no row constrains is no longer fixed at an infinite bound.
- A search that proved no bound reports no gap, where it reported 0.0000%, and nodes it
  skipped count against the bound it reports.
- The incumbent survives a root relaxation that does not finish.
- Mod-2 separation no longer lists a column twice when rows cancel it part way through
  a combination, which produced a cut that removed the optimum.

## 0.3.0

Neither 0.2.0 nor 0.1.x is affected by the correctness fixes below in the sense that
matters: they were wrong, and this release is where they stop being wrong.

### Changed behaviour you will notice

- **The default optimality gap is now `1e-4`, not zero.** Every established solver
  stops there, and proving the last ten-thousandth of a percent can cost more than
  everything before it: MIPLIB's `app2-1` reaches a 0.0006% gap almost at once and then
  cannot close it, running out a sixty second clock with the answer already in hand. At
  `1e-4` it finishes in about a second. Pass `--gap 0`, or set `gap_tolerance` to zero,
  to demand a proof.
- **Node-local cut separation is on**, at one node in ten. It shrinks trees where root
  cutting did not: `v064c200` 2690 nodes to 2116, `mkp_200` 72150 to 63896.
- **Root cut separation stays off** (`--cut-rounds N` enables it). It is worth turning
  on for knapsack-structured models.
- **MPS files with sections that change the model are refused** rather than read with
  the section dropped: `INDICATORS`, `SOS`, `QUADOBJ`, `QMATRIX`, `QCMATRIX`,
  `LAZYCONS`. Reading one before produced a confident answer to a different model.

### Correctness

- An MPS `INDICATORS` section was silently skipped, leaving conditional constraints
  standing unconditionally. MIPLIB's `cvrpsimple2i`, optimum 34, was reported infeasible.
- A failed factorization was reported as `Infeasible`. Failing to factorize a basis says
  nothing about whether a feasible point exists; `hypothyroid-k1` came back "Infeasible"
  where it is solvable.
- A free column displaced during singular-basis repair was parked at its upper bound,
  which for a free column is infinity.
- Time limits are now honoured during factorization and basis repair. Both are single
  operations that increment no iteration counter, so a deadline expiring inside one went
  unnoticed: a five second limit ran for 193 seconds on a 255386-row model.

### Performance

- Presolve pins singleton columns whose objective and single row agree on a direction.
  `neos-3048764-nadi` goes from a sixty second timeout to optimal in 0.33s at one node.
- Factorizing the all-logical starting basis was quadratic in the row count: 65 seconds
  for a 255386-row diagonal, now 50 milliseconds.
- Bland's rule is now left after a short spell and re-entered with escalation, rather
  than running an entire solve. `neos-555001`'s relaxation goes from 116953 pivots to
  4841.
- Cuts are selected by efficacy and orthogonality, aged out, and purged when slack.
  `mkp_200` with cutting on goes from 48.2s to about 10s.

### Added

- `model::Builder` assembles a model a column and a row at a time, in the sense you ask
  for, rather than filling `Problem`'s fields and negating a maximization by hand.
- `--values none|nonzero|all` prints solution values one per line, and `--solution PATH`
  writes them to a file.
- `--cuts-per-round` and `--local-cut-frequency`.

### Breaking

- `lu::Singular` is now `lu::FactorError`, with an `OutOfTime` variant.
- `BasisError` has an `OutOfTime` variant.
- `Lu::factor` and `Basis::refactorize` take a deadline.
- The CLI's `--no-cuts` became `--cut-rounds N` (0.2.0).

## 0.2.0

Never published.

## 0.1.1

Allocation reductions in the LU factorization.

## 0.1.0

First release: LP and MPS readers, presolve, cover and Gomory cuts, pseudocost
branching, primal heuristics, sparse LU, best-bound search, parallel tree search, and a
`gurobipy`-shaped Python module.
