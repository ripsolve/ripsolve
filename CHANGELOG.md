# Changelog

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
