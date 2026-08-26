//! Sparse LU factorization of a basis matrix, with Markowitz pivoting.
//!
//! The dense explicit inverse this replaces costs `O(m^2)` per solve and `O(m^3)`
//! to rebuild. At `m = 1000` the rebuild alone is ~10^9 operations every 50 pivots,
//! which measured out at 0.2 seconds per branch-and-bound node, 53x slower per
//! simplex iteration than a leading commercial solver on the same model.
//!
//! # Pivot choice
//!
//! Pivots are chosen by the Markowitz criterion: among the remaining entries,
//! minimize `(r_i - 1) * (c_j - 1)`, the number of positions the elimination could
//! turn nonzero. That single rule subsumes the special cases, a column singleton
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
#[derive(Clone)]
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
pub enum FactorError {
    /// A basis position with no acceptable pivot; the caller should replace it.
    Singular { position: usize },
    /// The caller's deadline passed mid-factorization.
    ///
    /// Factorizing is one call however long it takes, and on a large basis that is
    /// long: 255386 rows of MIPLIB's neos-4754521-awarau took 166 seconds, during
    /// which the simplex had not begun and so had nothing to check a clock between.
    /// A five second limit ran for 193.
    OutOfTime,
}

/// Relative threshold a pivot must meet against the largest entry in its column.
///
/// Relative to the column maximum, which makes the test exactly invariant to column
/// scaling, both sides scale together. That is worth knowing before reaching for
/// equilibration: scaling the basis before factorizing was implemented and
/// measured, and changed nothing. It reduced the coefficient range of a badly
/// scaled basis from 2e15 to 8e1 and left the residual where it was (1.4e-7 against
/// 3.3e-7, i.e. marginally worse), because the pivoting was never scale-naive.
///
/// Two things it does not address, either:
///
/// - A basis of independently random wide-range entries cannot be equilibrated at
///   all. Diagonal scaling has 2m degrees of freedom against m^2 entries, so the
///   residuals measured in `scale_tests` are ill-conditioning rather than bad
///   scaling, and no scaling step will fix them.
/// - Scaling the whole *model* rather than the basis is a different thing and
///   remains untried. It is what production solvers do, and the benefit is that
///   feasibility and optimality tolerances then mean the same thing in every row,
///   which this, operating inside the factorization, never touches.
const PIVOT_THRESHOLD: f64 = 0.01;

/// Columns examined per pivot search before settling for the best seen.
const MAX_COLUMN_SEARCH: usize = 4;

/// Elimination steps between clock reads during a factorization.
const FACTOR_CLOCK_INTERVAL: usize = 1024;

impl Lu {
    pub fn dimension(&self) -> usize {
        self.m
    }

    /// Factor the basis whose columns are given as `(rows, values)` pairs.
    pub fn factor(
        m: usize,
        columns: &[(Vec<usize>, Vec<f64>)],
        zero_tol: f64,
        deadline: Option<std::time::Instant>,
    ) -> Result<Lu, FactorError> {
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
        // Reused for every column update. Allocating one of these per touched
        // column per elimination step put ~18% of a solve in the allocator.
        let mut updated: Vec<(usize, f64)> = Vec::new();

        // Columns bucketed by nonzero count, so the search reaches the sparsest few
        // without ordering the rest.
        //
        // Sorting every remaining column at every step, which is what this replaced,
        // costs `O(m^2 log m)` over a factorization and exists only to read off the
        // first `MAX_COLUMN_SEARCH` entries. On a 3186-row basis that was a third of
        // the solve. Entries are left in place when a column's length changes and
        // skipped on the way past, so an update costs a push rather than a search.
        let mut buckets: Vec<Vec<usize>> = vec![Vec::new(); m + 1];
        let mut col_len: Vec<usize> = (0..m).map(|j| cols[j].len()).collect();
        for j in 0..m {
            buckets[col_len[j].min(m)].push(j);
        }
        // Reused across steps so the shortlist does not allocate per pivot.
        let mut shortlist: Vec<usize> = Vec::with_capacity(MAX_COLUMN_SEARCH);
        // Lower bound on any live column's length. Scanning from zero each step would
        // walk `m` empty buckets per pivot and cost exactly what the sort did.
        // Elimination can only shorten a column by removing the pivot row, so this
        // moves down by one at a time and is corrected on the way past.
        let mut min_len = 0usize;

        for step in 0..m {
            // One factorization is a single call to the caller, so the deadline has to
            // be read here or not at all until it finishes.
            if step.is_multiple_of(FACTOR_CLOCK_INTERVAL)
                && deadline.is_some_and(|d| std::time::Instant::now() >= d)
            {
                return Err(FactorError::OutOfTime);
            }
            shortlist.clear();
            let mut shortest_live = None;
            // Candidates are *popped* rather than scanned. Walking a bucket to compact
            // it costs its whole length at every step, and with `-I` every column sits
            // in the length-one bucket, so that was quadratic: factorizing the
            // all-logical basis of MIPLIB's neos-4754521-awarau, 255386 rows of
            // nothing but a diagonal, took 65 seconds. Popping stops at the few
            // candidates actually wanted, and stale entries are discarded as they
            // surface.
            // Indexing rather than iterating: the loop pops from the bucket it reads,
            // which an iterator would hold borrowed.
            #[allow(clippy::needless_range_loop)]
            'buckets: for len in min_len..=m {
                while shortlist.len() < MAX_COLUMN_SEARCH {
                    let Some(j) = buckets[len].pop() else {
                        break;
                    };
                    // The column has since been pivoted out, or moved bucket.
                    if col_done[j] || col_len[j] != len {
                        continue;
                    }
                    if shortest_live.is_none() {
                        shortest_live = Some(j);
                    }
                    shortlist.push(j);
                }
                if !shortlist.is_empty() {
                    // Everything below this length is exhausted, so later steps start
                    // here rather than walking the empty buckets again.
                    min_len = len;
                    break 'buckets;
                }
            }

            let mut best: Option<(usize, usize, f64, usize)> = None; // (row, col, value, cost)
            for (examined, &j) in shortlist.iter().enumerate() {
                if cols[j].is_empty() {
                    return Err(FactorError::Singular { position: j });
                }
                let max_abs = cols[j].iter().map(|&(_, v)| v.abs()).fold(0.0f64, f64::max);
                if max_abs <= zero_tol {
                    return Err(FactorError::Singular { position: j });
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
                // The sparsest column still active is the one to report, matching what
                // the caller repairs by swapping in a logical.
                let position = shortest_live.unwrap_or(0);
                return Err(FactorError::Singular { position });
            };

            // Candidates that were looked at but not pivoted are still live and were
            // taken out of their bucket by the pop, so they go back.
            for &j in &shortlist {
                if j != pj {
                    buckets[col_len[j].min(m)].push(j);
                }
            }

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
                updated.clear();
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
                // Swap rather than assign: `cols[j]`'s allocation becomes the next
                // iteration's buffer instead of being freed.
                std::mem::swap(&mut cols[j], &mut updated);
                if cols[j].len() != col_len[j] {
                    col_len[j] = cols[j].len();
                    let bucket = col_len[j].min(m);
                    buckets[bucket].push(j);
                    min_len = min_len.min(bucket);
                }
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
        // Renumber in place. Collecting into fresh vectors here allocated 2m of
        // them per factorization to change nothing but the indices.
        let mut l_cols = l_raw;
        for entries in &mut l_cols {
            for (i, _) in entries.iter_mut() {
                *i = row_pos[*i];
            }
        }
        let mut u_rows = u_raw;
        for entries in &mut u_rows {
            for (j, _) in entries.iter_mut() {
                *j = col_pos[*j];
            }
        }

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
        let lu = Lu::factor(m, &columns, 1e-12, None).unwrap();
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
        let lu = Lu::factor(m, &columns, 1e-12, None).unwrap();
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
            if let Ok(lu) = Lu::factor(m, &columns, 1e-12, None) {
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
            if let Ok(lu) = Lu::factor(m, &columns, 1e-12, None) {
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
        assert!(Lu::factor(m, &sparse(m, &dense), 1e-12, None).is_err());
    }

    #[test]
    fn reports_an_empty_column_as_singular() {
        let m = 3;
        let columns = vec![
            (vec![0], vec![1.0]),
            (Vec::new(), Vec::new()),
            (vec![2], vec![1.0]),
        ];
        let err = Lu::factor(m, &columns, 1e-12, None)
            .err()
            .expect("empty column is singular");
        assert_eq!(err, FactorError::Singular { position: 1 });
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
        let lu = Lu::factor(m, &columns, 1e-12, None).unwrap();
        assert_solves(m, &columns, &lu, 1e-7);
        // Fill-free would be 3m - 2 nonzeros; allow slack but not densification.
        assert!(
            lu.nnz() < 4 * m,
            "fill blew up: {} nonzeros for m = {m}",
            lu.nnz()
        );
    }
}

#[cfg(test)]
mod scale_tests {
    use super::*;

    /// SplitMix64, so a failure is reproducible from its seed.
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
            let v = (self.next() % 2001) as f64 / 100.0 - 10.0;
            if v.abs() < 0.5 { 1.0 } else { v }
        }
    }

    /// A basis shaped like one the simplex actually holds: mostly unit columns from
    /// logical variables, with a minority of sparse structural columns mixed in.
    ///
    /// Structure is the point. A dense random matrix makes every pivot choice
    /// equivalent, so it cannot distinguish a good Markowitz search from a bad one.
    fn realistic_basis(
        m: usize,
        structural: usize,
        per_column: usize,
        seed: u64,
    ) -> Vec<(Vec<usize>, Vec<f64>)> {
        let mut rng = Rng(seed);
        let mut columns: Vec<(Vec<usize>, Vec<f64>)> =
            (0..m).map(|i| (vec![i], vec![-1.0])).collect();
        for _ in 0..structural {
            let target = rng.below(m);
            let mut rows: Vec<usize> = Vec::new();
            // Keep the diagonal entry so the basis stays comfortably nonsingular.
            rows.push(target);
            for _ in 0..per_column.saturating_sub(1) {
                let r = rng.below(m);
                if !rows.contains(&r) {
                    rows.push(r);
                }
            }
            let mut pairs: Vec<(usize, f64)> = rows.into_iter().map(|r| (r, rng.value())).collect();
            pairs.sort_unstable_by_key(|&(r, _)| r);
            // A strong diagonal keeps the basis well conditioned, so a residual this
            // test rejects is the factorization's fault and not the matrix's.
            for (r, v) in &mut pairs {
                if *r == target {
                    *v = 25.0;
                }
            }
            columns[target] = (
                pairs.iter().map(|&(r, _)| r).collect(),
                pairs.iter().map(|&(_, v)| v).collect(),
            );
        }
        columns
    }

    fn dense_column(m: usize, column: &(Vec<usize>, Vec<f64>)) -> Vec<f64> {
        let mut out = vec![0.0; m];
        for (&r, &v) in column.0.iter().zip(&column.1) {
            out[r] = v;
        }
        out
    }

    /// The largest relative residual of `B (B^-1 a) = a` over random right-hand
    /// sides, how much accuracy the factorization actually delivers.
    ///
    /// Checking that a solve *runs* says nothing. A factorization built on bad
    /// pivots still produces numbers; they are simply wrong enough, far enough
    /// downstream, to make the simplex conclude something false.
    fn worst_residual(m: usize, columns: &[(Vec<usize>, Vec<f64>)], lu: &Lu) -> f64 {
        let dense: Vec<Vec<f64>> = columns.iter().map(|c| dense_column(m, c)).collect();
        let mut rng = Rng(99);
        let mut worst: f64 = 0.0;
        let mut d = Vec::new();
        for _ in 0..3 {
            let a: Vec<f64> = (0..m).map(|_| rng.value()).collect();
            lu.ftran(&a, &mut d);
            let scale = a.iter().fold(0.0f64, |acc, v| acc.max(v.abs())).max(1.0);
            for i in 0..m {
                let got: f64 = (0..m).map(|p| dense[p][i] * d[p]).sum();
                worst = worst.max((got - a[i]).abs() / scale);
            }
        }
        worst
    }

    /// Coefficients spanning `10^spread` either side of one, which is what makes a
    /// basis hard to factorize accurately. A matrix of ones and twos does not
    /// distinguish a good factorization from a bad one.
    fn spread_basis(
        m: usize,
        structural: usize,
        per_column: usize,
        spread: f64,
        seed: u64,
    ) -> Vec<(Vec<usize>, Vec<f64>)> {
        let mut rng = Rng(seed);
        let mut columns: Vec<(Vec<usize>, Vec<f64>)> =
            (0..m).map(|i| (vec![i], vec![-1.0])).collect();
        for _ in 0..structural {
            let target = rng.below(m);
            let mut rows = vec![target];
            for _ in 0..per_column.saturating_sub(1) {
                let r = rng.below(m);
                if !rows.contains(&r) {
                    rows.push(r);
                }
            }
            let mut pairs: Vec<(usize, f64)> = rows
                .into_iter()
                .map(|r| {
                    let exponent = ((rng.next() % 2001) as f64 / 1000.0 - 1.0) * spread;
                    let magnitude = 10f64.powf(exponent);
                    (
                        r,
                        if rng.next().is_multiple_of(2) {
                            magnitude
                        } else {
                            -magnitude
                        },
                    )
                })
                .collect();
            pairs.sort_unstable_by_key(|&(r, _)| r);
            columns[target] = (
                pairs.iter().map(|&(r, _)| r).collect(),
                pairs.iter().map(|&(_, v)| v).collect(),
            );
        }
        columns
    }

    #[test]
    fn stays_accurate_on_a_basis_large_enough_for_pivoting_to_matter() {
        // The gap this closes: every other LU test here uses a handful of rows,
        // where any pivot order works and a Markowitz search cannot be wrong. A
        // change to pivot selection once passed the whole suite and then reported a
        // feasible LP infeasible on a 402-row model.
        //
        // Checked as a *residual*, not as "the solve ran". A factorization built on
        // bad pivots still returns numbers; they are simply wrong enough, far enough
        // downstream, for the simplex to conclude something false.
        for (m, structural) in [(200, 60), (500, 150), (900, 250)] {
            let columns = spread_basis(m, structural, 8, 0.5, 7 + m as u64);
            let lu =
                Lu::factor(m, &columns, 1e-12, None).unwrap_or_else(|e| panic!("m = {m}: {e:?}"));
            let residual = worst_residual(m, &columns, &lu);
            assert!(
                residual < 1e-7,
                "m = {m}: residual {residual:.3e}, the factorization is losing accuracy"
            );
        }
    }

    /// What the factorization does as the data gets worse conditioned.
    ///
    /// Not an assertion of quality, a record of a known limitation, measured. The
    /// factorization has no scaling step, so accuracy tracks the spread of the
    /// coefficients directly:
    ///
    /// | coefficient range | residual |
    /// |---|---|
    /// | 1x        | 6e-14 |
    /// | 10x       | 6e-10 |
    /// | 100x      | 7e-8  |
    /// | 1000x     | 2e-3  |
    /// | 10000x    | 9e0   |
    ///
    /// By a range of 1e4 the result is meaningless, and real models routinely carry
    /// 1e4 to 1e6. Equilibrating the matrix before factorizing is the standard
    /// remedy and is not implemented. This test asserts only the part that holds
    /// today, so the boundary is recorded rather than assumed away.
    #[test]
    fn accuracy_tracks_how_badly_scaled_the_data_is() {
        let m = 400;
        let well = worst_residual(
            m,
            &spread_basis(m, 120, 8, 0.0, 3),
            &Lu::factor(m, &spread_basis(m, 120, 8, 0.0, 3), 1e-12, None).unwrap(),
        );
        assert!(
            well < 1e-10,
            "uniform coefficients should factorize cleanly, got {well:.3e}"
        );

        // And the documented failure, so a future scaling step has something to beat.
        let bad_columns = spread_basis(m, 120, 8, 2.0, 3);
        let bad = worst_residual(
            m,
            &bad_columns,
            &Lu::factor(m, &bad_columns, 1e-12, None).unwrap(),
        );
        assert!(
            bad > well,
            "a 1e4 coefficient range is expected to be worse than a uniform one; \
             if this fails, scaling has been added and the bound above should be tightened"
        );
    }

    #[test]
    fn markowitz_keeps_fill_bounded_at_scale() {
        // What the pivot search is *for*. A search that picks poorly still produces
        // a usable factorization on a small matrix, and buries the solver in fill on
        // a real one, so this is checked where the difference shows.
        for (m, structural) in [(300, 90), (800, 240)] {
            let columns = spread_basis(m, structural, 8, 0.5, 11 + m as u64);
            let original: usize = columns.iter().map(|c| c.0.len()).sum();
            let lu = Lu::factor(m, &columns, 1e-12, None).unwrap();
            assert!(
                lu.nnz() < 6 * original,
                "m = {m}: {} nonzeros from {original}, fill has blown up",
                lu.nnz()
            );
        }
    }

    #[test]
    fn accuracy_survives_a_long_run_of_updates() {
        // The product-form file grows between refactorizations, and its error grows
        // with it. This is the regime the simplex spends its life in.
        let m = 400;
        let mut columns = realistic_basis(m, 120, 8, 31);
        let mut basis = crate::lp::basis::Basis::all_logical(m);
        basis.refactorize(&columns, 1e-9, None).unwrap();

        let mut rng = Rng(5);
        for step in 0..60 {
            let position = rng.below(m);
            let replacement = realistic_basis(m, 1, 6, 100 + step)[0].clone();
            let dense = dense_column(m, &replacement);
            let mut d = Vec::new();
            basis.ftran(&dense, &mut d);
            // Skip a pivot too small to be safe, exactly as the simplex would.
            if d[position].abs() < 1e-7 {
                continue;
            }
            basis.update(&d, position);
            columns[position] = replacement;
        }

        // B^-1 B must still be the identity after all of that.
        let dense: Vec<Vec<f64>> = columns.iter().map(|c| dense_column(m, c)).collect();
        let mut out = Vec::new();
        let mut worst: f64 = 0.0;
        for (j, column) in dense.iter().enumerate() {
            basis.ftran(column, &mut out);
            for (i, &got) in out.iter().enumerate() {
                let expected = if i == j { 1.0 } else { 0.0 };
                worst = worst.max((got - expected).abs());
            }
        }
        assert!(
            worst < 1e-6,
            "after 60 updates the inverse is off by {worst:.3e}"
        );
    }
}
