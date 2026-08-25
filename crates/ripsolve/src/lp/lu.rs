//! Sparse LU factorization of a basis matrix, with Markowitz pivoting.
//!
//! The dense explicit inverse this replaces costs `O(m^2)` per solve and `O(m^3)`
//! to rebuild. At `m = 1000` the rebuild alone is ~10^9 operations every 50 pivots,
//! which measured out at 0.2 seconds per branch-and-bound node — 53x slower per
//! simplex iteration than a leading commercial solver on the same model.
//!
//! # Pivot choice
//!
//! Pivots are chosen by the Markowitz criterion: among the remaining entries,
//! minimize `(r_i - 1) * (c_j - 1)`, the number of positions the elimination could
//! turn nonzero. That single rule subsumes the special cases — a column singleton
//! scores zero and is taken immediately, so an LP basis's large triangular part
//! (every logical variable contributes one) is peeled off without any separate
//! triangularization pass.
//!
//! Sparsity alone is not safe, so a candidate must also satisfy a relative
//! threshold, `|a_ij| >= tau * max|a_.j|`. Pure Markowitz would happily pivot on an
//! arbitrarily small entry and produce enormous multipliers.
//!
//! The search is bounded: it examines columns in increasing order of count and
//! stops once a pivot cannot be beaten. Searching every remaining entry at every
//! step would itself be quadratic.

/// The factorization `P B Q = L U`, with `L` unit lower triangular.
pub struct Lu {
    m: usize,
    /// Original row index of pivot `k`.
    prow: Vec<usize>,
    /// Basis position of pivot `k`.
    pcol: Vec<usize>,
    /// `L` by pivot step, in permuted index space: `(k', multiplier)` with `k' > k`.
    l_cols: Vec<Vec<(usize, f64)>>,
    /// `U` by pivot step, in permuted index space: `(k', value)` with `k' > k`.
    u_rows: Vec<Vec<(usize, f64)>>,
    /// Diagonal of `U`.
    diag: Vec<f64>,
}

/// The basis could not be factored.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Singular {
    /// A basis position with no acceptable pivot; the caller should replace it.
    pub position: usize,
}

/// Relative threshold a pivot must meet against the largest entry in its column.
const PIVOT_THRESHOLD: f64 = 0.01;
/// Columns examined per pivot search before settling for the best seen.
const MAX_COLUMN_SEARCH: usize = 4;

impl Lu {
    pub fn dimension(&self) -> usize {
        self.m
    }

    /// Factor the basis whose columns are given as `(rows, values)` pairs.
    pub fn factor(
        m: usize,
        columns: &[(Vec<usize>, Vec<f64>)],
        zero_tol: f64,
    ) -> Result<Lu, Singular> {
        debug_assert_eq!(columns.len(), m);

        // Active submatrix, column-major by basis position.
        let mut cols: Vec<Vec<(usize, f64)>> = columns
            .iter()
            .map(|(rows, vals)| {
                rows.iter()
                    .copied()
                    .zip(vals.iter().copied())
                    .filter(|&(_, v)| v != 0.0)
                    .collect()
            })
            .collect();
        // Which columns touch each row. Entries may go stale; readers filter.
        let mut row_cols: Vec<Vec<usize>> = vec![Vec::new(); m];
        let mut row_nz = vec![0usize; m];
        for (j, col) in cols.iter().enumerate() {
            for &(i, _) in col {
                row_cols[i].push(j);
                row_nz[i] += 1;
            }
        }

        let mut row_done = vec![false; m];
        let mut col_done = vec![false; m];
        let mut prow = Vec::with_capacity(m);
        let mut pcol = Vec::with_capacity(m);
        // Collected in original index space; permuted once the order is known.
        let mut l_raw: Vec<Vec<(usize, f64)>> = Vec::with_capacity(m);
        let mut u_raw: Vec<Vec<(usize, f64)>> = Vec::with_capacity(m);
        let mut diag = Vec::with_capacity(m);

        // Reused scatter buffers, so column updates never allocate.
        let mut scatter = vec![0.0f64; m];
        let mut present = vec![false; m];

        // The still-active columns, compacted as they are eliminated.
        let mut order: Vec<usize> = (0..m).collect();
        // The few shortest of them, refilled each step.
        let mut shortlist: Vec<(usize, usize)> = Vec::with_capacity(MAX_COLUMN_SEARCH);

        for _step in 0..m {
            order.retain(|&j| !col_done[j]);

            // Select the shortest few columns rather than sorting all of them.
            //
            // This was a full sort at every one of the m elimination steps, which is
            // O(m * n log n) per factorization. It did not show on small models and
            // dominated everything else on real ones: profiling a 3904-row instance
            // put 60% of the entire solve in that one sort. Only MAX_COLUMN_SEARCH
            // columns are ever examined, so sorting the rest was pure waste.
            shortlist.clear();
            for &j in &order {
                let len = cols[j].len();
                if shortlist.len() == MAX_COLUMN_SEARCH && len >= shortlist[MAX_COLUMN_SEARCH - 1].1
                {
                    continue;
                }
                let at = shortlist
                    .iter()
                    .position(|&(_, l)| len < l)
                    .unwrap_or(shortlist.len());
                shortlist.insert(at, (j, len));
                shortlist.truncate(MAX_COLUMN_SEARCH);
            }

            let mut best: Option<(usize, usize, f64, usize)> = None; // (row, col, value, cost)
            for (examined, &(j, _)) in shortlist.iter().enumerate() {
                if cols[j].is_empty() {
                    return Err(Singular { position: j });
                }
                let max_abs = cols[j].iter().map(|&(_, v)| v.abs()).fold(0.0f64, f64::max);
                if max_abs <= zero_tol {
                    return Err(Singular { position: j });
                }
                for &(i, v) in &cols[j] {
                    // Stability first: a sparser pivot is worthless if it is tiny.
                    if v.abs() < PIVOT_THRESHOLD * max_abs {
                        continue;
                    }
                    let cost = (row_nz[i].saturating_sub(1)) * (cols[j].len().saturating_sub(1));
                    if best.is_none_or(|(_, _, _, c)| cost < c) {
                        best = Some((i, j, v, cost));
                    }
                }
                // A zero-cost pivot cannot be improved on, and the search is bounded
                // regardless so that finding a pivot stays sub-quadratic.
                if best.is_some_and(|(_, _, _, c)| c == 0) || examined + 1 >= MAX_COLUMN_SEARCH {
                    break;
                }
            }

            let Some((pi, pj, pivot, _)) = best else {
                let position = order.first().copied().unwrap_or(0);
                return Err(Singular { position });
            };

            row_done[pi] = true;
            col_done[pj] = true;
            prow.push(pi);
            pcol.push(pj);
            diag.push(pivot);

            // Scatter the pivot column, dropping its own pivot row.
            let pivot_col: Vec<(usize, f64)> =
                cols[pj].iter().copied().filter(|&(i, _)| i != pi).collect();
            for &(i, v) in &pivot_col {
                scatter[i] = v;
                present[i] = true;
            }

            // L takes the multipliers below the pivot.
            let l_entries: Vec<(usize, f64)> =
                pivot_col.iter().map(|&(i, v)| (i, v / pivot)).collect();

            // Eliminate the pivot row from every other active column.
            let mut u_entries: Vec<(usize, f64)> = Vec::new();
            let touching: Vec<usize> = row_cols[pi]
                .iter()
                .copied()
                .filter(|&j| !col_done[j])
                .collect();
            for j in touching {
                let Some(pos) = cols[j].iter().position(|&(i, _)| i == pi) else {
                    continue;
                };
                let alpha = cols[j][pos].1;
                if alpha == 0.0 {
                    continue;
                }
                u_entries.push((j, alpha));
                let factor = alpha / pivot;

                // col_j <- col_j - factor * pivot_col, with the pivot row removed.
                let mut updated: Vec<(usize, f64)> = Vec::with_capacity(cols[j].len());
                for &(i, v) in &cols[j] {
                    if i == pi {
                        continue;
                    }
                    let value = if present[i] {
                        v - factor * scatter[i]
                    } else {
                        v
                    };
                    if value.abs() > zero_tol {
                        updated.push((i, value));
                    } else {
                        row_nz[i] -= 1;
                    }
                }
                // Fill-in: pivot-column rows this column did not already have.
                for &(i, v) in &pivot_col {
                    if !cols[j].iter().any(|&(r, _)| r == i) {
                        let value = -factor * v;
                        if value.abs() > zero_tol {
                            updated.push((i, value));
                            row_cols[i].push(j);
                            row_nz[i] += 1;
                        }
                    }
                }
                row_nz[pi] -= 1;
                cols[j] = updated;
            }

            for &(i, _) in &pivot_col {
                scatter[i] = 0.0;
                present[i] = false;
                row_nz[i] -= 1;
            }
            row_nz[pi] = 0;
            cols[pj].clear();

            l_raw.push(l_entries);
            u_raw.push(u_entries);
        }

        // Translate to permuted index space, now that the pivot order is known.
        let mut row_pos = vec![usize::MAX; m];
        let mut col_pos = vec![usize::MAX; m];
        for k in 0..m {
            row_pos[prow[k]] = k;
            col_pos[pcol[k]] = k;
        }
        let l_cols: Vec<Vec<(usize, f64)>> = l_raw
            .into_iter()
            .map(|entries| entries.into_iter().map(|(i, v)| (row_pos[i], v)).collect())
            .collect();
        let u_rows: Vec<Vec<(usize, f64)>> = u_raw
            .into_iter()
            .map(|entries| entries.into_iter().map(|(j, v)| (col_pos[j], v)).collect())
            .collect();

        Ok(Lu {
            m,
            prow,
            pcol,
            l_cols,
            u_rows,
            diag,
        })
    }

    /// Solve `B d = a`, writing `d` (indexed by basis position) into `out`.
    pub fn ftran(&self, a: &[f64], out: &mut Vec<f64>) {
        let m = self.m;
        debug_assert_eq!(a.len(), m);
        // Permute the right-hand side into pivot order, then solve L z = b.
        let mut work: Vec<f64> = (0..m).map(|k| a[self.prow[k]]).collect();
        for k in 0..m {
            let z = work[k];
            if z != 0.0 {
                for &(k2, mult) in &self.l_cols[k] {
                    work[k2] -= mult * z;
                }
            }
        }
        // Back-substitute U w = z.
        for k in (0..m).rev() {
            let mut acc = work[k];
            for &(k2, v) in &self.u_rows[k] {
                acc -= v * work[k2];
            }
            work[k] = acc / self.diag[k];
        }
        out.clear();
        out.resize(m, 0.0);
        for k in 0..m {
            out[self.pcol[k]] = work[k];
        }
    }

    /// Solve `B' y = c`, writing `y` (indexed by row) into `out`.
    pub fn btran(&self, c: &[f64], out: &mut Vec<f64>) {
        let m = self.m;
        debug_assert_eq!(c.len(), m);
        // Permute by column order, then solve U' v = c (U' is lower triangular).
        let mut work: Vec<f64> = (0..m).map(|k| c[self.pcol[k]]).collect();
        for k in 0..m {
            let v = work[k] / self.diag[k];
            work[k] = v;
            if v != 0.0 {
                for &(k2, u) in &self.u_rows[k] {
                    work[k2] -= u * v;
                }
            }
        }
        // Then L' u = v, which is upper triangular in this ordering.
        for k in (0..m).rev() {
            let mut acc = work[k];
            for &(k2, mult) in &self.l_cols[k] {
                acc -= mult * work[k2];
            }
            work[k] = acc;
        }
        out.clear();
        out.resize(m, 0.0);
        for k in 0..m {
            out[self.prow[k]] = work[k];
        }
    }

    /// Stored nonzeros in both factors, for diagnostics and refactorization policy.
    pub fn nnz(&self) -> usize {
        self.m
            + self.l_cols.iter().map(Vec::len).sum::<usize>()
            + self.u_rows.iter().map(Vec::len).sum::<usize>()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Deterministic PRNG, so a failing case is reproducible from its seed.
    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = self.0;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            z ^ (z >> 31)
        }
        fn below(&mut self, n: usize) -> usize {
            (self.next() % n as u64) as usize
        }
        fn value(&mut self) -> f64 {
            (self.next() % 19) as f64 - 9.0
        }
    }

    /// Dense columns from sparse ones, for the reference solves.
    fn densify(m: usize, columns: &[(Vec<usize>, Vec<f64>)]) -> Vec<Vec<f64>> {
        columns
            .iter()
            .map(|(rows, vals)| {
                let mut col = vec![0.0; m];
                for (&i, &v) in rows.iter().zip(vals) {
                    col[i] = v;
                }
                col
            })
            .collect()
    }

    /// Check `B (B^-1 a) == a` and `(B' y) == c`, which pins both solves without
    /// needing a reference inverse.
    fn assert_solves(m: usize, columns: &[(Vec<usize>, Vec<f64>)], lu: &Lu, tol: f64) {
        let dense = densify(m, columns);
        let mut rng = Rng(7);
        for _ in 0..4 {
            let a: Vec<f64> = (0..m).map(|_| rng.value()).collect();

            let mut d = Vec::new();
            lu.ftran(&a, &mut d);
            // B d must reproduce a.
            for i in 0..m {
                let got: f64 = (0..m).map(|p| dense[p][i] * d[p]).sum();
                assert!((got - a[i]).abs() < tol, "ftran row {i}: {got} vs {}", a[i]);
            }

            let c: Vec<f64> = (0..m).map(|_| rng.value()).collect();
            let mut y = Vec::new();
            lu.btran(&c, &mut y);
            // y' B must reproduce c.
            for p in 0..m {
                let got: f64 = (0..m).map(|i| dense[p][i] * y[i]).sum();
                assert!(
                    (got - c[p]).abs() < tol,
                    "btran column {p}: {got} vs {}",
                    c[p]
                );
            }
        }
    }

    fn sparse(m: usize, dense: &[Vec<f64>]) -> Vec<(Vec<usize>, Vec<f64>)> {
        dense
            .iter()
            .map(|col| {
                let mut rows = Vec::new();
                let mut vals = Vec::new();
                for (i, &v) in col.iter().enumerate().take(m) {
                    if v != 0.0 {
                        rows.push(i);
                        vals.push(v);
                    }
                }
                (rows, vals)
            })
            .collect()
    }

    #[test]
    fn factors_the_logical_basis() {
        // The all-logical starting basis is -I, the case every solve begins from.
        let m = 5;
        let columns: Vec<(Vec<usize>, Vec<f64>)> = (0..m).map(|i| (vec![i], vec![-1.0])).collect();
        let lu = Lu::factor(m, &columns, 1e-12).unwrap();
        assert_solves(m, &columns, &lu, 1e-9);
        // Perfectly triangular: no fill beyond the diagonal.
        assert_eq!(lu.nnz(), m);
    }

    #[test]
    fn factors_a_permuted_triangular_basis() {
        // Markowitz should peel this off with zero fill, choosing singletons first.
        let m = 4;
        let dense = vec![
            vec![0.0, 0.0, 2.0, 0.0],
            vec![1.0, 0.0, 0.0, 0.0],
            vec![0.0, 3.0, 1.0, 0.0],
            vec![0.0, 0.0, 0.0, 4.0],
        ];
        let columns = sparse(m, &dense);
        let lu = Lu::factor(m, &columns, 1e-12).unwrap();
        assert_solves(m, &columns, &lu, 1e-9);
    }

    #[test]
    fn factors_a_dense_basis() {
        let m = 6;
        let mut rng = Rng(11);
        loop {
            let dense: Vec<Vec<f64>> = (0..m)
                .map(|_| (0..m).map(|_| rng.value()).collect())
                .collect();
            let columns = sparse(m, &dense);
            if let Ok(lu) = Lu::factor(m, &columns, 1e-12) {
                assert_solves(m, &columns, &lu, 1e-6);
                return;
            }
        }
    }

    #[test]
    fn factors_random_sparse_bases() {
        // The broad net: many shapes, each verified by round-tripping both solves.
        let mut rng = Rng(2024);
        let mut factored = 0;
        for _ in 0..60 {
            let m = 8 + rng.below(20);
            // Start from the identity so the basis is generally nonsingular, then
            // scatter extra entries to force fill-in and real pivot choices.
            let mut dense: Vec<Vec<f64>> = (0..m)
                .map(|j| {
                    let mut col = vec![0.0; m];
                    col[j] = 1.0 + rng.below(4) as f64;
                    col
                })
                .collect();
            for _ in 0..(m * 3) {
                let (j, i) = (rng.below(m), rng.below(m));
                let v = rng.value();
                if v != 0.0 {
                    dense[j][i] = v;
                }
            }
            let columns = sparse(m, &dense);
            if let Ok(lu) = Lu::factor(m, &columns, 1e-12) {
                assert_solves(m, &columns, &lu, 1e-5);
                factored += 1;
            }
        }
        assert!(factored > 40, "only {factored} of 60 bases factored");
    }

    #[test]
    fn reports_a_singular_basis() {
        // Third column is the sum of the first two.
        let m = 3;
        let dense = vec![
            vec![1.0, 0.0, 1.0],
            vec![0.0, 1.0, 1.0],
            vec![1.0, 1.0, 2.0],
        ];
        assert!(Lu::factor(m, &sparse(m, &dense), 1e-12).is_err());
    }

    #[test]
    fn reports_an_empty_column_as_singular() {
        let m = 3;
        let columns = vec![
            (vec![0], vec![1.0]),
            (Vec::new(), Vec::new()),
            (vec![2], vec![1.0]),
        ];
        let err = Lu::factor(m, &columns, 1e-12)
            .err()
            .expect("empty column is singular");
        assert_eq!(err, Singular { position: 1 });
    }

    #[test]
    fn markowitz_keeps_an_arrow_matrix_sparse() {
        // An arrow with the dense row and column last: pivoting on the diagonal
        // spokes first costs no fill, while pivoting on the hub first would make the
        // whole remaining matrix dense.
        let m = 30;
        let mut dense: Vec<Vec<f64>> = (0..m)
            .map(|j| {
                let mut col = vec![0.0; m];
                col[j] = 2.0;
                col
            })
            .collect();
        for col in dense.iter_mut().take(m - 1) {
            col[m - 1] = 1.0;
        }
        let hub = &mut dense[m - 1];
        for entry in hub.iter_mut().take(m - 1) {
            *entry = 1.0;
        }
        let columns = sparse(m, &dense);
        let lu = Lu::factor(m, &columns, 1e-12).unwrap();
        assert_solves(m, &columns, &lu, 1e-7);
        // Fill-free would be 3m - 2 nonzeros; allow slack but not densification.
        assert!(
            lu.nnz() < 4 * m,
            "fill blew up: {} nonzeros for m = {m}",
            lu.nnz()
        );
    }
}
