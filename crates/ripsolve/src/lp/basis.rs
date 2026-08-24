//! The basis inverse, and the solves the simplex needs against it.
//!
//! Representation: a sparse LU factorization (see [`crate::lp::lu`]) plus a
//! *product form* update — an "eta file" of rank-one corrections, one per pivot,
//! replayed on top of the factors. Periodically the etas are discarded and the
//! basis refactorized, which both bounds the replay cost and stops error
//! accumulating.
//!
//! # Why product form rather than Forrest-Tomlin
//!
//! Forrest-Tomlin updates the `U` factor in place and keeps the eta file shorter,
//! which is what a mature solver does. Product form is a few lines by comparison
//! and is what makes the correctness argument easy to see: after a pivot on row
//! `r`, the new basis is `B_new = B_old E`, where `E` is the identity with column
//! `r` replaced by `d = B_old^-1 a_q`. So `B_new^-1 = E^-1 B_old^-1`, and after `k`
//! pivots `B_k^-1 = E_k^-1 ... E_1^-1 B_0^-1`. FTRAN therefore applies the factors
//! first and the etas in order; BTRAN applies the transposed etas in reverse and
//! the factors last.
//!
//! The interface — [`Basis::ftran`], [`Basis::btran`], [`Basis::update`],
//! [`Basis::refactorize`] — is unchanged from the dense explicit inverse this
//! replaced, so Forrest-Tomlin can supersede product form later without the simplex
//! driver noticing, exactly as this change did.

use crate::lp::lu::{Lu, Singular};

/// Why a refactorization could not produce a usable basis.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BasisError {
    /// The basis is singular at this position; the caller must repair it (typically
    /// by swapping in a logical) and retry.
    Singular { row: usize },
}

/// One product-form update: the transformed entering column and its pivot row.
struct Eta {
    row: usize,
    pivot: f64,
    /// The column's other nonzeros, `(row, value)` with `row != self.row`.
    column: Vec<(usize, f64)>,
}

/// A factorized basis.
pub struct Basis {
    m: usize,
    lu: Lu,
    etas: Vec<Eta>,
}

/// Entries smaller than this are treated as structurally absent.
const ZERO_TOL: f64 = 1e-12;

impl Basis {
    /// The all-logical starting basis, whose matrix is `-I`.
    pub fn all_logical(m: usize) -> Self {
        let columns: Vec<(Vec<usize>, Vec<f64>)> = (0..m).map(|i| (vec![i], vec![-1.0])).collect();
        let lu = Lu::factor(m, &columns, ZERO_TOL).expect("-I is nonsingular");
        Self {
            m,
            lu,
            etas: Vec::new(),
        }
    }

    pub fn dimension(&self) -> usize {
        self.m
    }

    /// Pivots applied since the last refactorization.
    pub fn updates(&self) -> usize {
        self.etas.len()
    }

    /// Nonzeros held in the factors, for diagnostics.
    pub fn nnz(&self) -> usize {
        self.lu.nnz() + self.etas.iter().map(|e| e.column.len() + 1).sum::<usize>()
    }

    /// FTRAN: solve `B d = a`, returning `d = B^-1 a`.
    pub fn ftran(&self, a: &[f64], out: &mut Vec<f64>) {
        self.lu.ftran(a, out);
        for eta in &self.etas {
            // v_r <- v_r / pivot, then v_i <- v_i - d_i * v_r.
            let scaled = out[eta.row] / eta.pivot;
            out[eta.row] = scaled;
            if scaled != 0.0 {
                for &(i, v) in &eta.column {
                    out[i] -= v * scaled;
                }
            }
        }
    }

    /// BTRAN: solve `B' y = c`, returning `y' = c' B^-1`.
    pub fn btran(&self, c: &[f64], out: &mut Vec<f64>) {
        out.clear();
        out.extend_from_slice(c);
        // Transposed etas, newest first: only the pivot row's entry changes.
        for eta in self.etas.iter().rev() {
            let mut acc = out[eta.row];
            for &(i, v) in &eta.column {
                acc -= v * out[i];
            }
            out[eta.row] = acc / eta.pivot;
        }
        let mut work = std::mem::take(out);
        self.lu.btran(&work, out);
        work.clear();
    }

    /// BTRAN against a unit vector: row `r` of `B^-1`.
    ///
    /// The dual simplex prices a whole pivot row at once and needs this. With an
    /// explicit inverse it was a row copy; here it is an ordinary BTRAN, which is
    /// why it belongs behind the interface rather than being read off the
    /// representation by the caller.
    pub fn btran_unit(&self, r: usize, out: &mut Vec<f64>) {
        let mut unit = vec![0.0; self.m];
        unit[r] = 1.0;
        self.btran(&unit, out);
    }

    /// Record a pivot on row `r`, where `d = B^-1 a_q` for the entering column.
    ///
    /// The caller is responsible for having rejected a pivot too small to be safe.
    pub fn update(&mut self, d: &[f64], r: usize) {
        debug_assert_eq!(d.len(), self.m);
        if self.m == 0 {
            return;
        }
        let pivot = d[r];
        debug_assert!(pivot != 0.0, "pivot on a zero element");
        let column = d
            .iter()
            .enumerate()
            .filter(|&(i, &v)| i != r && v.abs() > ZERO_TOL)
            .map(|(i, &v)| (i, v))
            .collect();
        self.etas.push(Eta {
            row: r,
            pivot,
            column,
        });
    }

    /// Refactorize from the basis columns, given as `(rows, values)` pairs.
    pub fn refactorize(
        &mut self,
        columns: &[(Vec<usize>, Vec<f64>)],
        _pivot_tol: f64,
    ) -> Result<(), BasisError> {
        if self.m == 0 {
            self.etas.clear();
            return Ok(());
        }
        match Lu::factor(self.m, columns, ZERO_TOL) {
            Ok(lu) => {
                self.lu = lu;
                self.etas.clear();
                Ok(())
            }
            Err(Singular { position }) => Err(BasisError::Singular { row: position }),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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

    /// `B^-1 B` must be the identity, checked column by column.
    fn assert_inverts(basis: &Basis, dense: &[Vec<f64>]) {
        let m = basis.dimension();
        let mut out = Vec::new();
        for (j, col) in dense.iter().enumerate() {
            basis.ftran(col, &mut out);
            for (i, &got) in out.iter().enumerate().take(m) {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((got - expected).abs() < 1e-9, "B^-1 B [{i}][{j}] = {got}");
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
        let dense = vec![
            vec![2.0, 1.0, 1.0],
            vec![1.0, 3.0, 2.0],
            vec![1.0, 0.0, 4.0],
        ];
        let mut b = Basis::all_logical(3);
        b.refactorize(&sparse(3, &dense), 1e-9).unwrap();
        assert_inverts(&b, &dense);
        assert_eq!(b.updates(), 0);
    }

    #[test]
    fn refactorize_reports_a_singular_basis() {
        let dense = vec![
            vec![1.0, 0.0, 1.0],
            vec![0.0, 1.0, 1.0],
            vec![1.0, 1.0, 2.0],
        ];
        let mut b = Basis::all_logical(3);
        assert!(matches!(
            b.refactorize(&sparse(3, &dense), 1e-9),
            Err(BasisError::Singular { .. })
        ));
    }

    #[test]
    fn update_matches_a_refactorization_of_the_same_basis() {
        // Replacing a basis column must give the same inverse whether reached by a
        // product-form update or by refactorizing the new basis outright. This is
        // the property the whole eta file rests on.
        let mut dense = vec![
            vec![2.0, 1.0, 1.0],
            vec![1.0, 3.0, 2.0],
            vec![1.0, 0.0, 4.0],
        ];
        let entering = vec![3.0, -1.0, 2.0];

        let mut updated = Basis::all_logical(3);
        updated.refactorize(&sparse(3, &dense), 1e-9).unwrap();
        let mut d = Vec::new();
        updated.ftran(&entering, &mut d);
        updated.update(&d, 1);

        dense[1] = entering;
        assert_inverts(&updated, &dense);
        assert_eq!(updated.updates(), 1);

        let mut fresh = Basis::all_logical(3);
        fresh.refactorize(&sparse(3, &dense), 1e-9).unwrap();
        let probe = [1.0, -2.0, 0.5];
        let (mut a, mut b) = (Vec::new(), Vec::new());
        updated.ftran(&probe, &mut a);
        fresh.ftran(&probe, &mut b);
        for (x, y) in a.iter().zip(&b) {
            assert!((x - y).abs() < 1e-9, "{x} vs {y}");
        }
        updated.btran(&probe, &mut a);
        fresh.btran(&probe, &mut b);
        for (x, y) in a.iter().zip(&b) {
            assert!((x - y).abs() < 1e-9, "{x} vs {y}");
        }
    }

    #[test]
    fn many_updates_stay_consistent() {
        // The eta file must remain correct as it grows, not just for one pivot.
        let mut dense = vec![
            vec![4.0, 1.0, 0.0, 1.0],
            vec![1.0, 3.0, 1.0, 0.0],
            vec![0.0, 1.0, 5.0, 1.0],
            vec![1.0, 0.0, 1.0, 2.0],
        ];
        let mut basis = Basis::all_logical(4);
        basis.refactorize(&sparse(4, &dense), 1e-9).unwrap();

        for (step, replacement) in [
            (0usize, vec![1.0, 2.0, 0.0, 1.0]),
            (2, vec![0.0, 1.0, 3.0, 1.0]),
            (1, vec![2.0, 0.0, 1.0, 1.0]),
            (3, vec![1.0, 1.0, 1.0, 4.0]),
        ] {
            let mut d = Vec::new();
            basis.ftran(&replacement, &mut d);
            basis.update(&d, step);
            dense[step] = replacement;
            assert_inverts(&basis, &dense);
        }
        assert_eq!(basis.updates(), 4);
    }

    #[test]
    fn btran_unit_is_a_row_of_the_inverse() {
        let dense = vec![vec![2.0, 1.0], vec![1.0, 3.0]];
        let mut b = Basis::all_logical(2);
        b.refactorize(&sparse(2, &dense), 1e-9).unwrap();
        for r in 0..2 {
            let mut rho = Vec::new();
            b.btran_unit(r, &mut rho);
            for (j, col) in dense.iter().enumerate() {
                let dot: f64 = rho.iter().zip(col).map(|(a, b)| a * b).sum();
                let expected = if j == r { 1.0 } else { 0.0 };
                assert!((dot - expected).abs() < 1e-9, "row {r}, column {j}: {dot}");
            }
        }
    }

    #[test]
    fn a_zero_row_basis_is_handled() {
        // A model with no rows has an empty basis; nothing here may panic on it.
        let mut b = Basis::all_logical(0);
        assert!(b.refactorize(&[], 1e-9).is_ok());
        let mut out = Vec::new();
        b.ftran(&[], &mut out);
        assert!(out.is_empty());
        b.update(&[], 0);
        assert_eq!(b.updates(), 0);
    }
}
