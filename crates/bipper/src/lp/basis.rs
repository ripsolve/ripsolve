//! The basis inverse, and the two solves the simplex needs against it.
//!
//! Representation: a dense, explicitly maintained `B^-1`, updated in place by
//! elementary row operations after each pivot and periodically recomputed from
//! scratch to stop error accumulating.
//!
//! That is deliberately the simplest thing that is *correct*. A production simplex
//! keeps a sparse LU factorization with Forrest-Tomlin updates instead, which is
//! both faster and numerically stronger. The point of this module's narrow surface
//! — [`Basis::ftran`], [`Basis::btran`], [`Basis::update`], [`Basis::refactorize`]
//! — is that the swap can happen later without the simplex driver noticing. Until
//! the benchmarks demand it, an explicit inverse costs O(m^2) per iteration and is
//! far easier to verify.

/// Why a refactorization could not produce a usable basis.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BasisError {
    /// The basis matrix is singular even after pivoting; the caller must repair the
    /// basis (typically by swapping the offending column for a logical) and retry.
    Singular { row: usize },
}

/// A dense inverse of the current basis matrix.
pub struct Basis {
    m: usize,
    /// `B^-1`, row-major, `m x m`.
    inv: Vec<f64>,
    /// Pivots applied since the last refactorization.
    updates: usize,
}

impl Basis {
    /// The inverse of the all-logical starting basis.
    ///
    /// Logical `i` enters the computational matrix as `-e_i`, so that basis is `-I`
    /// and is its own inverse.
    pub fn all_logical(m: usize) -> Self {
        let mut inv = vec![0.0; m * m];
        for i in 0..m {
            inv[i * m + i] = -1.0;
        }
        Self { m, inv, updates: 0 }
    }

    pub fn dimension(&self) -> usize {
        self.m
    }

    /// Pivots applied since the last refactorization.
    pub fn updates(&self) -> usize {
        self.updates
    }

    fn row(&self, i: usize) -> &[f64] {
        &self.inv[i * self.m..(i + 1) * self.m]
    }

    /// FTRAN: solve `B d = a`, returning `d = B^-1 a`.
    pub fn ftran(&self, a: &[f64], out: &mut Vec<f64>) {
        debug_assert_eq!(a.len(), self.m);
        out.clear();
        out.resize(self.m, 0.0);
        for (i, slot) in out.iter_mut().enumerate() {
            let row = &self.inv[i * self.m..(i + 1) * self.m];
            // Skipping zeros pays off because entering columns are typically sparse.
            *slot = a
                .iter()
                .zip(row)
                .filter(|(v, _)| **v != 0.0)
                .map(|(v, r)| r * v)
                .sum();
        }
    }

    /// BTRAN: solve `B' y = c`, returning `y' = c' B^-1`.
    pub fn btran(&self, c: &[f64], out: &mut Vec<f64>) {
        debug_assert_eq!(c.len(), self.m);
        out.clear();
        out.resize(self.m, 0.0);
        // y_k = sum_i c_i * inv[i][k]; skipping zero c_i matters because in phase 1
        // the cost vector is mostly zeros.
        for (i, &ci) in c.iter().enumerate() {
            if ci == 0.0 {
                continue;
            }
            let row = self.row(i);
            for k in 0..self.m {
                out[k] += ci * row[k];
            }
        }
    }

    /// Apply the rank-one update for a pivot on row `r`, where `d = B^-1 a_q` was
    /// computed for the entering column `a_q`.
    ///
    /// `d[r]` is the pivot element; the caller is responsible for having rejected a
    /// pivot too small to be safe.
    pub fn update(&mut self, d: &[f64], r: usize) {
        debug_assert_eq!(d.len(), self.m);
        let m = self.m;
        if m == 0 {
            return;
        }
        let pivot = d[r];
        debug_assert!(pivot != 0.0, "pivot on a zero element");

        let scale = 1.0 / pivot;
        for v in &mut self.inv[r * m..(r + 1) * m] {
            *v *= scale;
        }
        // Copied once, outside the loop: it aliases the mutable row borrow below, and
        // hoisting it turns m allocations per pivot into one.
        let pivot_row: Vec<f64> = self.row(r).to_vec();
        for (i, (&factor, row)) in d.iter().zip(self.inv.chunks_exact_mut(m)).enumerate() {
            if i == r || factor == 0.0 {
                continue;
            }
            for (target, p) in row.iter_mut().zip(&pivot_row) {
                *target -= factor * p;
            }
        }
        self.updates += 1;
    }

    /// Recompute `B^-1` from scratch by Gauss-Jordan elimination with partial
    /// pivoting, given the basis columns as dense vectors.
    ///
    /// Returns the row whose pivot was unusable if the basis is singular, so the
    /// caller can repair that position and try again.
    pub fn refactorize(&mut self, columns: &[Vec<f64>], pivot_tol: f64) -> Result<(), BasisError> {
        let m = self.m;
        debug_assert_eq!(columns.len(), m);
        // A model with no rows has an empty basis, which is already its own inverse.
        if m == 0 {
            self.updates = 0;
            return Ok(());
        }

        // Work on [B | I] and reduce the left half to the identity.
        let mut work = vec![0.0f64; m * 2 * m];
        for (j, col) in columns.iter().enumerate() {
            debug_assert_eq!(col.len(), m);
            for i in 0..m {
                work[i * 2 * m + j] = col[i];
            }
        }
        for i in 0..m {
            work[i * 2 * m + m + i] = 1.0;
        }

        for c in 0..m {
            let (mut best, mut best_val) = (c, work[c * 2 * m + c].abs());
            for r in c + 1..m {
                let v = work[r * 2 * m + c].abs();
                if v > best_val {
                    best = r;
                    best_val = v;
                }
            }
            if best_val <= pivot_tol {
                return Err(BasisError::Singular { row: c });
            }
            if best != c {
                for k in 0..2 * m {
                    work.swap(c * 2 * m + k, best * 2 * m + k);
                }
            }
            let pivot = work[c * 2 * m + c];
            let scale = 1.0 / pivot;
            for k in 0..2 * m {
                work[c * 2 * m + k] *= scale;
            }
            for r in 0..m {
                if r == c {
                    continue;
                }
                let factor = work[r * 2 * m + c];
                if factor == 0.0 {
                    continue;
                }
                for k in 0..2 * m {
                    work[r * 2 * m + k] -= factor * work[c * 2 * m + k];
                }
            }
        }

        // Copy the reduced right half of [B | I], which is now B^-1.
        for (i, row) in self.inv.chunks_exact_mut(m).enumerate() {
            let start = i * 2 * m + m;
            row.copy_from_slice(&work[start..start + m]);
        }
        self.updates = 0;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dense_identity_check(basis: &Basis, columns: &[Vec<f64>]) {
        // B^-1 B must be the identity.
        let m = basis.dimension();
        let mut out = Vec::new();
        for (j, col) in columns.iter().enumerate() {
            basis.ftran(col, &mut out);
            assert_eq!(out.len(), m);
            for (i, &got) in out.iter().enumerate() {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (got - expected).abs() < 1e-9,
                    "B^-1 B [{i}][{j}] = {got}, expected {expected}"
                );
            }
        }
    }

    #[test]
    fn logical_basis_is_negative_identity() {
        let b = Basis::all_logical(3);
        let mut out = Vec::new();
        b.ftran(&[1.0, 2.0, 3.0], &mut out);
        assert_eq!(out, vec![-1.0, -2.0, -3.0]);
        b.btran(&[1.0, 2.0, 3.0], &mut out);
        assert_eq!(out, vec![-1.0, -2.0, -3.0]);
    }

    #[test]
    fn refactorize_inverts() {
        let columns = vec![
            vec![2.0, 1.0, 1.0],
            vec![1.0, 3.0, 2.0],
            vec![1.0, 0.0, 4.0],
        ];
        let mut b = Basis::all_logical(3);
        b.refactorize(&columns, 1e-9).unwrap();
        dense_identity_check(&b, &columns);
        assert_eq!(b.updates(), 0);
    }

    #[test]
    fn refactorize_reports_a_singular_basis() {
        // Third column is the sum of the first two.
        let columns = vec![
            vec![1.0, 0.0, 1.0],
            vec![0.0, 1.0, 1.0],
            vec![1.0, 1.0, 2.0],
        ];
        let mut b = Basis::all_logical(3);
        assert!(matches!(
            b.refactorize(&columns, 1e-9),
            Err(BasisError::Singular { .. })
        ));
    }

    #[test]
    fn refactorize_pivots_around_a_zero_leading_entry() {
        // Needs a row swap: the (0,0) entry is zero but the matrix is nonsingular.
        let columns = vec![vec![0.0, 1.0], vec![1.0, 0.0]];
        let mut b = Basis::all_logical(2);
        b.refactorize(&columns, 1e-9).unwrap();
        dense_identity_check(&b, &columns);
    }

    #[test]
    fn update_matches_a_refactorization_of_the_same_basis() {
        // Replacing basis column 1 by `entering` must give the same inverse whether
        // reached by a rank-one update or by inverting the new basis outright.
        let mut columns = vec![
            vec![2.0, 1.0, 1.0],
            vec![1.0, 3.0, 2.0],
            vec![1.0, 0.0, 4.0],
        ];
        let entering = vec![3.0, -1.0, 2.0];

        let mut updated = Basis::all_logical(3);
        updated.refactorize(&columns, 1e-9).unwrap();
        let mut d = Vec::new();
        updated.ftran(&entering, &mut d);
        updated.update(&d, 1);

        columns[1] = entering;
        let mut fresh = Basis::all_logical(3);
        fresh.refactorize(&columns, 1e-9).unwrap();

        for (a, b) in updated.inv.iter().zip(&fresh.inv) {
            assert!((a - b).abs() < 1e-9, "{a} vs {b}");
        }
        dense_identity_check(&updated, &columns);
        assert_eq!(updated.updates(), 1);
    }

    #[test]
    fn btran_solves_the_transposed_system() {
        let columns = vec![vec![2.0, 1.0], vec![1.0, 3.0]];
        let mut b = Basis::all_logical(2);
        b.refactorize(&columns, 1e-9).unwrap();

        let c = [5.0, 7.0];
        let mut y = Vec::new();
        b.btran(&c, &mut y);
        // y' B == c'
        for (j, col) in columns.iter().enumerate() {
            let dot: f64 = y.iter().zip(col).map(|(a, b)| a * b).sum();
            assert!((dot - c[j]).abs() < 1e-9, "column {j}: {dot} vs {}", c[j]);
        }
    }
}
