//! The problem model: a binary integer program.
//!
//! Rows are stored in *range* form (`lb <= a'x <= ub`, with infinite bounds
//! allowed) rather than as a sense plus a right-hand side. That is what HiGHS and
//! Gurobi keep internally, and it pays off twice: `<=`, `>=`, `=` and true range
//! rows all become one case for the simplex, and presolve's bound tightening is
//! just an update to `lb`/`ub` instead of a change of sense.
//!
//! Every column is binary. Bounds are still carried explicitly per column so that
//! presolve fixing and branching can express `x_j = 0` / `x_j = 1` as `lb == ub`
//! without a separate mechanism.

use crate::sparse::SparseMatrix;

/// Optimization direction of the *original* problem as the user stated it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Sense {
    Minimize,
    Maximize,
}

/// A binary integer program.
///
/// The solver works internally in minimization form. [`Problem::objective_value`]
/// converts an internal objective back to the user's original sense and offset, so
/// callers never have to remember whether a negation happened.
#[derive(Clone, Debug)]
pub struct Problem {
    pub name: String,
    /// Direction the user asked for. Coefficients in `obj` are always in
    /// minimization form; for a maximization they were negated on the way in.
    pub sense: Sense,
    /// Objective coefficients, minimization form, one per column.
    pub obj: Vec<f64>,
    /// Constant term of the internal (minimization) objective.
    pub obj_offset: f64,
    /// Constraint matrix, `n_rows x n_cols`, column-major.
    pub matrix: SparseMatrix,
    pub row_lb: Vec<f64>,
    pub row_ub: Vec<f64>,
    pub col_lb: Vec<f64>,
    pub col_ub: Vec<f64>,
    pub col_names: Vec<String>,
    pub row_names: Vec<String>,
}

/// How a row was written before conversion to range form.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RowSense {
    /// `a'x <= rhs`
    Le,
    /// `a'x >= rhs`
    Ge,
    /// `a'x == rhs`
    Eq,
}

impl RowSense {
    /// The `(lb, ub)` range-form bounds for this sense and right-hand side.
    pub fn bounds(self, rhs: f64) -> (f64, f64) {
        match self {
            RowSense::Le => (f64::NEG_INFINITY, rhs),
            RowSense::Ge => (rhs, f64::INFINITY),
            RowSense::Eq => (rhs, rhs),
        }
    }
}

impl Problem {
    pub fn n_cols(&self) -> usize {
        self.obj.len()
    }

    pub fn n_rows(&self) -> usize {
        self.row_lb.len()
    }

    /// Convert an internal (minimization) objective value back to the value in the
    /// user's original sense, including the offset introduced during construction.
    pub fn objective_value(&self, internal: f64) -> f64 {
        let with_offset = internal + self.obj_offset;
        match self.sense {
            Sense::Minimize => with_offset,
            Sense::Maximize => -with_offset,
        }
    }

    /// Check the internal invariants a freshly built problem must satisfy.
    ///
    /// Readers and presolve both go through this, so a malformed model is caught at
    /// the boundary rather than as a confusing failure deep in the simplex.
    pub fn validate(&self) -> Result<(), ModelError> {
        let (n, m) = (self.n_cols(), self.n_rows());
        if self.matrix.n_cols() != n || self.matrix.n_rows() != m {
            return Err(ModelError::ShapeMismatch {
                matrix: (self.matrix.n_rows(), self.matrix.n_cols()),
                vectors: (m, n),
            });
        }
        for (len, what) in [
            (self.col_lb.len(), "col_lb"),
            (self.col_ub.len(), "col_ub"),
            (self.col_names.len(), "col_names"),
        ] {
            if len != n {
                return Err(ModelError::LengthMismatch {
                    what,
                    got: len,
                    expected: n,
                });
            }
        }
        for (len, what) in [
            (self.row_ub.len(), "row_ub"),
            (self.row_names.len(), "row_names"),
        ] {
            if len != m {
                return Err(ModelError::LengthMismatch {
                    what,
                    got: len,
                    expected: m,
                });
            }
        }
        // A NaN bound would silently poison every comparison in the simplex, so it is
        // rejected here rather than allowed to propagate.
        for j in 0..n {
            let (lo, hi) = (self.col_lb[j], self.col_ub[j]);
            if lo.is_nan() || hi.is_nan() || lo > hi {
                return Err(ModelError::InvalidBounds {
                    what: "column",
                    index: j,
                    lb: lo,
                    ub: hi,
                });
            }
            // Binary columns may only be relaxed to [0,1] or fixed within it.
            if lo < 0.0 || hi > 1.0 {
                return Err(ModelError::NotBinary {
                    index: j,
                    lb: lo,
                    ub: hi,
                });
            }
        }
        for i in 0..m {
            let (lo, hi) = (self.row_lb[i], self.row_ub[i]);
            if lo.is_nan() || hi.is_nan() || lo > hi {
                return Err(ModelError::InvalidBounds {
                    what: "row",
                    index: i,
                    lb: lo,
                    ub: hi,
                });
            }
        }
        if self.obj.iter().any(|v| !v.is_finite()) || !self.obj_offset.is_finite() {
            return Err(ModelError::NonFiniteObjective);
        }
        Ok(())
    }
}

/// A malformed [`Problem`].
#[derive(Debug, thiserror::Error)]
pub enum ModelError {
    #[error("matrix is {}x{} but the model describes {}x{}", matrix.0, matrix.1, vectors.0, vectors.1)]
    ShapeMismatch {
        matrix: (usize, usize),
        vectors: (usize, usize),
    },
    #[error("{what} has length {got}, expected {expected}")]
    LengthMismatch {
        what: &'static str,
        got: usize,
        expected: usize,
    },
    #[error("{what} {index} has invalid bounds [{lb}, {ub}]")]
    InvalidBounds {
        what: &'static str,
        index: usize,
        lb: f64,
        ub: f64,
    },
    #[error(
        "column {index} has bounds [{lb}, {ub}], which is not within [0, 1]; ripsolve solves binary programs only"
    )]
    NotBinary { index: usize, lb: f64, ub: f64 },
    #[error("objective contains a non-finite coefficient")]
    NonFiniteObjective,
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tiny(sense: Sense, offset: f64) -> Problem {
        Problem {
            name: "tiny".into(),
            sense,
            obj: vec![1.0, 2.0],
            obj_offset: offset,
            matrix: SparseMatrix::from_triplets(1, 2, [(0, 0, 1.0), (0, 1, 1.0)]),
            row_lb: vec![1.0],
            row_ub: vec![f64::INFINITY],
            col_lb: vec![0.0, 0.0],
            col_ub: vec![1.0, 1.0],
            col_names: vec!["x0".into(), "x1".into()],
            row_names: vec!["c0".into()],
        }
    }

    #[test]
    fn row_sense_to_range_bounds() {
        assert_eq!(RowSense::Le.bounds(3.0), (f64::NEG_INFINITY, 3.0));
        assert_eq!(RowSense::Ge.bounds(3.0), (3.0, f64::INFINITY));
        assert_eq!(RowSense::Eq.bounds(3.0), (3.0, 3.0));
    }

    #[test]
    fn objective_value_undoes_negation_and_offset() {
        // A maximization is stored negated, so the round trip must flip it back.
        assert_eq!(tiny(Sense::Minimize, 5.0).objective_value(2.0), 7.0);
        assert_eq!(tiny(Sense::Maximize, 5.0).objective_value(2.0), -7.0);
    }

    #[test]
    fn validate_accepts_a_well_formed_problem() {
        assert!(tiny(Sense::Minimize, 0.0).validate().is_ok());
    }

    #[test]
    fn validate_rejects_non_binary_columns() {
        let mut p = tiny(Sense::Minimize, 0.0);
        p.col_ub[1] = 4.0;
        assert!(matches!(
            p.validate(),
            Err(ModelError::NotBinary { index: 1, .. })
        ));
    }

    #[test]
    fn validate_rejects_crossed_and_nan_bounds() {
        let mut p = tiny(Sense::Minimize, 0.0);
        p.col_lb[0] = 1.0;
        p.col_ub[0] = 0.0;
        assert!(matches!(
            p.validate(),
            Err(ModelError::InvalidBounds { .. })
        ));

        let mut p = tiny(Sense::Minimize, 0.0);
        p.row_lb[0] = f64::NAN;
        assert!(matches!(
            p.validate(),
            Err(ModelError::InvalidBounds { .. })
        ));
    }

    #[test]
    fn validate_rejects_shape_and_length_errors() {
        let mut p = tiny(Sense::Minimize, 0.0);
        p.obj.push(3.0);
        assert!(matches!(
            p.validate(),
            Err(ModelError::ShapeMismatch { .. })
        ));

        let mut p = tiny(Sense::Minimize, 0.0);
        p.col_names.pop();
        assert!(matches!(
            p.validate(),
            Err(ModelError::LengthMismatch {
                what: "col_names",
                ..
            })
        ));
    }
}
