# Changelog

## Unreleased

### Added

- `search::solve_with` and the `Separator` trait: rows supplied while the search runs,
  asked for at the root, at every node, and of every point before it may become the
  incumbent. It serves as a cut callback and as lazy constraints at once. A clique
  partitioning model whose 52000 transitivity rows are left to the separator solves in
  one call rather than in an outer loop of whole MIP solves.
- `Lp::fit`, `BasisState::without_rows` and `lu::Scratch`.
- `presolve::presolve_with`, which can leave out the reductions that fix a column by
  its cost.

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
