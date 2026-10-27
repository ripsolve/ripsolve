# Changelog

## Unreleased

### Changed behaviour you will notice

- **`--time-limit` counts from when the command starts**, so reading the model is part
  of it. On a large file that is much of the wait: `neos-3402454-bohle`'s 346MB take
  about twenty seconds to read, and a 60-second limit returned after 106s. It now
  returns after 62s. The library's `Options::time_limit` is unchanged and still covers
  the solve alone.
- **Heuristics run on work budgets, not shares of the time limit.** The feasibility
  jump is given work in proportion to the model's nonzeros, and the RINS and reduced
  cost neighbourhood searches a number of simplex iterations scaled from the root
  relaxation, with ten cut rounds at their roots. The time shares they had remain as
  backstops. What they find no longer depends on the limit or on how busy the machine
  is wherever the run has time to spare: `mc11` at 60 seconds on one thread, run
  twice under load, now finds the same points in the same order. Where the root's own
  cut loop runs into its share of the limit, which on short runs is common, results
  still vary, so at 30 seconds two runs agree about as often as before. Over 54
  instances at 30 seconds `drayage-100-23` goes from 506846 to 178247, `n5-3` from
  9331 to 8938 and `neos-911970` from 56.87 to 55.33; none is worse beyond what
  repeated runs of one build show.
  `Options` gains `deadline`, a moment to stop at that nothing takes a share of,
  `max_iterations` and `improvement_iterations`.

### Fixed

- **A bound proven before a restart is kept.** A later pass that stopped before proving
  a bound of its own reported none: `h80x6320d` at 30 seconds on one thread said "no
  bound proven" after its first pass had proven 6344, and now reports a gap of 1.55%.
- **The feasibility jump reads the clock by work done, not by moves.** A move on a
  dense model rescores thousands of columns, and 512 moves between readings took
  `ns1952667` to 64 seconds against a 60-second limit. It now stops at 60.1.
- **Time limits held on models with very dense rows.** Three loops were quadratic in
  row length and none of them checked the clock: coefficient tightening in presolve
  recomputed a row's activity for every entry, and the conflict-driven search and the
  fixing heuristics reread every row whose slack had not changed. `pb-fit2d` at a
  60-second limit returned after 285s and now returns after 63s; `nw04` at a 10-second
  limit on one thread ran to 14s and now stops on time. The conflict-driven search also
  checks the deadline on every conflict, not every 64th.
- **A search stopped by the clock no longer says NodeLimit.** A node whose LP runs into
  the deadline is skipped, and when the skipped nodes were the last ones the tree
  emptied and reported a node limit nobody had set. Five of MIPLIB's benchmark
  instances did this at 60 seconds, among them `tbfp-network` and `eilA101-2`; they
  now report TimeLimit.
- **A pass that stopped to restart no longer claims to be optimal.** The serial search
  left its status as it stood when it broke off to restart with nodes still open. The
  restart loop discarded it, but a pass that settled its model is now the answer rather
  than a reason to restart, and that exposed it: `h80x6320d` reported an optimum of 6386
  against MIPLIB's 6382. The search also asks whether to restart only of a node it would
  solve, so a root whose cuts had closed the gap ends the search instead.
- **Feasibility jump no longer panics on crossed bounds.** A column whose bounds arrive
  crossed by rounding, `[25, 24.999999999999986]` on `h80x6320d`, made `f64::clamp`
  panic.
- Separating implied-bound cuts reads each clique once rather than once per column in
  it. The cuts are the same.
- **The simplex reads the clock often enough on large LPs.** It read it every 256
  pivots whatever the model's size, and a pivot on `neos-3402454-bohle` takes about a
  tenth of a second, so half a minute could pass between readings. The interval now
  shrinks with the LP's size, down to every pivot.
- **A phase 1 that never ended on wide-ranging columns.** The primal ratio test leaves
  out rows whose entry is too small to pivot on, and left them out of limiting the step
  as well, so a step could push a feasible variable in such a row out of bounds. On
  `leo1`, whose columns run from 1 to 9e7, phase 1 flipped one binary back and forth
  for a million iterations; its relaxation now solves in 1103, and `leo2`'s in 1956.
  Where the step was already feasible nothing changes.
- A root relaxation that ran out of iterations before the deadline reported
  `NodeLimit`, a limit nobody had set. It now keeps solving until the deadline, and
  reports `TimeLimit` if that comes first.

### Performance

General MIP, which the pure binary set does not exercise. Measured on instances from
MIPLIB 2017's benchmark set that HiGHS and SCIP both close in seconds, at thirty
seconds on one thread.

- **Feasibility jump moves general integer and continuous columns.** It was binary
  only, declining any model with a general integer and leaving continuous columns where
  they started; now each such column jumps to the value in its range that minimises the
  weighted violation of its rows. It finds a first point in under three seconds on
  `neos-3381206-awhea`, `neos-1582420`, `neos-860300`, `app1-1` and `neos8`, none of
  which had one. Binary columns flip as before.
- **RINS at the root**, as soon as there is an incumbent, and not only every 500 nodes:
  the models with a poor first point are often the ones whose nodes are slow enough
  never to reach 500. Its neighbourhood is searched with the fixed columns removed from
  the model rather than pinned in it. `beasleyC3` goes from 1783 to 807 against an
  optimum of 754, and `neos-860300` from no point to 4161 against 3201.
- **MIR cuts from equality rows, from rows combined to cancel continuous columns, and
  from rows whose columns sit away from zero.** Equalities were skipped, which left a
  network model's node balances uncut; combining rows is the aggregation of Marchand and
  Wolsey; and a column with a nonzero lower bound, every time variable of a scheduling
  model, is now measured from its nearer bound instead of keeping its row out.
  `sp150x300d` closes at the root, 35 to 69, the optimum; `exp-1-500-5-5` closes where
  it ended at a 48.6% gap; `neos-3381206-awhea` closes at 453.
- **MIR puts idle fixed-charge arcs on their binaries.** A continuous column with a
  variable upper bound was substituted only when strictly nearer that bound than zero,
  and an arc carrying nothing with its binary off is at both at once, so every idle arc
  stayed as flow and the cuts read `y + x'/2 + ... >= 1`, which the next relaxation
  escapes by moving the flow. Ties now go by the coefficient's sign. `p200x1188c`'s root
  goes from 7954 to 15078, the optimum, and the model closes in 6 seconds at 60 on 12
  threads; `cbs-cta` closes in 29; `h80x6320d` ends at a 1% gap where it ended at 17%.
- **Factorizing a basis of mostly logicals is linear rather than quadratic.** The
  column singletons a basis starts with are pivoted in one sweep before the Markowitz
  elimination, which took them one at a time and rebuilt every column each one's row
  touched. On `neos8`, 46324 rows, a factorization took 1.5 seconds; the root relaxation
  goes from 52 seconds to 3.7, and the model now closes.
- **A warm dual solve stalled on zero ratios perturbs its costs** before handing the
  node to the primal method. After ten pivots in a row at a zero ratio, each nonbasic
  column's cost is shifted slightly onto the side its bound needs, the dual method runs
  on for a bounded spell, and the true costs come back at the end. Only on models where
  at least a tenth of the columns carry a cost; on nearly flat objectives it invented
  one and spent its spell optimizing it. Over the warm child solves of a fixed probe,
  pivots fall to 0.84 and time to 0.80 of what they were: `nursesched-sprint02` from
  736 pivots a child to 50, `n2seq36q` from 424 to 123, `neos-1171448` from 217 to 63.
  Timed searches are mixed.
- **Reduced costs are carried across pivots** instead of recomputed. Both simplex
  methods solved for the duals and priced every column afresh each pivot, a BTRAN and
  a pass over the whole matrix; they now update them from the pivot row, which is
  formed from the rows of the matrix it touches, and recompute them only at a
  refactorization. Optimality and infeasibility are still decided on duals computed
  from scratch. Devex reads the same pivot row instead of forming its own.
- **The LU factors are kept in both orientations**, so that every triangular solve
  skips the zeros of its right-hand side; FTRAN through `U` and BTRAN through `L'`
  read every entry of the factor whatever they were given.
- **Pivots go into the factors as Forrest-Tomlin updates**, not into an eta file
  every BTRAN replays whole. An update that would eliminate with a multiplier above a
  hundred is refused and goes into the eta file as before. The pivot element is now
  compared between the entering column and the pivot row on every pivot, and the
  basis refactorized when they disagree past 1e-7.

  Together, on the root relaxations of the 41 benchmark instances HiGHS closes within
  a minute and this solver does not, the 33 that finish take 0.58 of the time they
  did (geometric mean), 333 seconds in all against 588: `hypothyroid-k1` 103s to 39,
  `qap10` 35 to 17, `nursesched-sprint02` 13.5 to 5.4, `dano3_3` 75 to 40. Every
  optimum is unchanged. Three, `ex10`, `physiciansched6-2` and `unitcal_7`, still do
  not finish in 150 seconds.
- **Pricing and the ratio tests read state carried per column and per row.** Which
  side of its bounds each nonbasic sits on, and how far each basic variable is outside
  its bounds, are now updated as pivots change them instead of being worked out from
  the status and bounds at every scan. The primal method's pricing scan had been a
  third of its own time on `neos-827175` and `dano3_3`. On the 29 of the 41 above whose
  root relaxations finish in five seconds, root time falls to 0.90 with the same
  pivots; `neos-827175` 54s to 44, `dano3_3` 45 to 38. Child solves are unchanged.
- **A warm start whose factors are not cached no longer factorizes the all-logical
  basis** only to discard it at the refactorization that follows. Factors from a
  refactorization that failed are no longer put in the cache under the basis they do
  not belong to.

  Over the 240 benchmark instances at 80 seconds on 12 threads, 24 close against 23,
  `neos-1582420` the new one. Root cuts depend on which optimal vertex the relaxation
  ends at, and the new engine ends at different ones: `neos8` now takes 27 seconds
  against 7, its first cut round leaving the bound where it was instead of closing the
  gap, while `neos-911970` and `neos-860300` end with better bounds and
  `neos-4738912-atrato`, `swath3` and `neos17` with worse.

### Benchmarks

- The pure binary MIPLIB set, 140 instances at sixty seconds, rerun with every answer
  checked against MIPLIB's solution file: 51 closed at sixteen threads (HiGHS 46, SCIP
  45, CBC 30, Gurobi 94) and 47 at one thread with the 0.4.0 build (HiGHS 39, SCIP 41,
  CBC 31, Gurobi 81). Optima agree with MIPLIB's to within the default 1e-4 gap. The
  results are in `docs/baselines/binary-2026-10-23.json` and the README.
- `bench/benchmark_set.py` runs MIPLIB 2017's 240 benchmark instances and checks each
  answer against MIPLIB.
- MIPLIB 2017's 240 benchmark instances at 80 seconds on twelve threads, the 60
  seconds of Mittelmann's runs scaled to this machine, which runs HiGHS 1.33 times
  slower than his: 23 closed, none wrong, and no run past 82.2 seconds. The check now
  also counts a bound that passes MIPLIB's optimum as wrong. `ex10`, which closes in
  about 26 seconds, ran out of time in two of seven runs with its root relaxation
  unfinished, for reasons not yet found. The results are in
  `docs/baselines/benchmark-2026-10-25.json`.

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
