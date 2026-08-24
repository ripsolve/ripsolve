//! Reading problems from LP and MPS files.
//!
//! Parsing is delegated to `lp_parser_rs`; this module's job is the translation
//! into [`Problem`] — range-form rows and a minimization-form objective.
//!
//! Column and row order follow the order the names appear in the source file.
//! `lp_parser_rs` keys its collections with `IndexMap`, so that order is stable,
//! and stability matters well beyond tidiness: column indices decide branching
//! order and degenerate pivot choices, so an order that varied between runs would
//! make node counts irreproducible on identical input.

use std::collections::HashSet;
use std::path::Path;

use lp_parser_rs::model::{ComparisonOp, Constraint, Sense as LpSense, VariableType};
use lp_parser_rs::parser::parse_file;
use lp_parser_rs::problem::LpProblem;

use crate::model::{Problem, RowSense, Sense, VarType};
use crate::sparse::SparseMatrix;

/// Names declared integral by an LP file's `Binary`, `General` or `Integer`
/// sections, recovered from the source text.
///
/// This exists to work around `lp_parser_rs`: a `Bounds` section overwrites a
/// variable's type, so `x` declared under `General` and then bounded to `[0, 10]`
/// comes back as a plain double-bounded *continuous* variable, its integrality
/// gone. The bounds it reports are right; only the type is lost. Rather than
/// re-parse the file, this recovers the one fact the parser drops and leaves
/// everything else to it.
///
/// Section headers end the list, so it stops at `Bounds`, `End`, and the rest.
fn declared_integer(text: &str) -> HashSet<String> {
    // Headers that begin a list of integral variables.
    const INTEGRAL: [&str; 9] = [
        "binaries", "binary", "bin", "generals", "general", "gen", "integers", "integer", "int",
    ];
    // Any other header ends one. `subject` and `such` cover the two spellings of
    // the constraint header without matching a variable named `st`.
    const OTHER: [&str; 12] = [
        "bounds", "bound", "end", "maximize", "maximise", "minimize", "minimise", "max", "min",
        "subject", "such", "sos",
    ];

    let mut names = HashSet::new();
    let mut collecting = false;
    for line in text.lines() {
        // Strip comments, which run from a backslash to the end of the line.
        let line = line.split('\\').next().unwrap_or("").trim();
        if line.is_empty() {
            continue;
        }
        let first = line
            .split_whitespace()
            .next()
            .unwrap_or("")
            .to_ascii_lowercase();
        let first = first.trim_end_matches(':').to_string();

        if INTEGRAL.contains(&first.as_str()) {
            collecting = true;
            // A header may carry names on the same line.
            names.extend(line.split_whitespace().skip(1).map(str::to_string));
            continue;
        }
        if OTHER.contains(&first.as_str()) {
            collecting = false;
            continue;
        }
        if collecting {
            names.extend(line.split_whitespace().map(str::to_string));
        }
    }
    names
}

/// Failure to turn a model file into a [`Problem`].
#[derive(Debug, thiserror::Error)]
pub enum ReadError {
    #[error("could not read or parse model file: {0}")]
    Parse(String),
    #[error("unrecognised model file extension {0:?}; expected .lp or .mps")]
    UnknownFormat(String),
    #[error("the model has no objective function")]
    NoObjective,
    #[error("the model has {0} objectives; ripsolve supports exactly one")]
    MultipleObjectives(usize),
    #[error("variable {0:?} is {1}, which ripsolve does not support")]
    UnsupportedVariable(String, String),
    #[error("integer variable {0:?} has the fractional bound {1}")]
    FractionalIntegerBound(String, f64),
    #[error("constraint {0:?} is a special ordered set, which ripsolve does not support")]
    SosConstraint(String),
}

/// Which parser to use for a model file.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Format {
    Lp,
    Mps,
}

impl Format {
    /// Infer the format from a path's extension, case-insensitively.
    pub fn from_path(path: &Path) -> Result<Format, ReadError> {
        match path.extension().and_then(|e| e.to_str()) {
            Some(e) if e.eq_ignore_ascii_case("lp") => Ok(Format::Lp),
            Some(e) if e.eq_ignore_ascii_case("mps") => Ok(Format::Mps),
            other => Err(ReadError::UnknownFormat(other.unwrap_or("").to_string())),
        }
    }
}

impl Problem {
    /// Read a binary integer program from a file, choosing the parser by extension.
    pub fn from_file(path: &Path) -> Result<Problem, ReadError> {
        Problem::from_file_as(path, Format::from_path(path)?)
    }

    /// Read a binary integer program from a file in an explicitly given format.
    pub fn from_file_as(path: &Path, format: Format) -> Result<Problem, ReadError> {
        let content = parse_file(path).map_err(|e| ReadError::Parse(e.to_string()))?;
        let lp = match format {
            Format::Lp => LpProblem::parse(&content),
            Format::Mps => LpProblem::parse_mps(&content),
        }
        .map_err(|e| ReadError::Parse(e.to_string()))?;
        // The LP source is consulted only to recover integrality the parser drops;
        // see `declared_integer`.
        let integral = match format {
            Format::Lp => declared_integer(&content),
            Format::Mps => HashSet::new(),
        };
        Problem::from_lp_with_integrality(&lp, &integral)
    }

    /// Translate an already-parsed [`LpProblem`].
    ///
    /// A `Maximize` objective is negated so the solver always works in
    /// minimization form; [`Problem::sense`] records what was asked for and
    /// [`Problem::objective_value`] converts results back.
    pub fn from_lp(lp: &LpProblem) -> Result<Problem, ReadError> {
        Problem::from_lp_with_integrality(lp, &HashSet::new())
    }

    /// As [`Problem::from_lp`], with extra names known to be integral.
    pub fn from_lp_with_integrality(
        lp: &LpProblem,
        integral: &HashSet<String>,
    ) -> Result<Problem, ReadError> {
        if lp.objectives.len() > 1 {
            return Err(ReadError::MultipleObjectives(lp.objectives.len()));
        }
        let objective = lp
            .objectives
            .values()
            .next()
            .ok_or(ReadError::NoObjective)?;
        let name = |id| lp.interner.resolve(id).to_string();

        // Column types and bounds, from whichever section declared them.
        //
        // A variable's *type* and its *bounds* are separate in both formats: an
        // integer declared with no bounds is `[0, inf)`, a binary one is an integer
        // pinned to `[0, 1]`, and a bounds section can narrow either. Reading them
        // as one thing is how a model silently becomes the wrong model.
        let mut col_names = Vec::with_capacity(lp.variables.len());
        let mut col_type = Vec::with_capacity(lp.variables.len());
        let mut col_lb = Vec::with_capacity(lp.variables.len());
        let mut col_ub = Vec::with_capacity(lp.variables.len());

        for var in lp.variables.values() {
            let label = name(var.name);
            let (mut kind, lo, hi) = match &var.var_type {
                VariableType::Binary => (VarType::Integer, 0.0, 1.0),
                VariableType::Integer | VariableType::General => {
                    (VarType::Integer, 0.0, f64::INFINITY)
                }
                VariableType::Free => (VarType::Continuous, f64::NEG_INFINITY, f64::INFINITY),
                VariableType::LowerBound(lo) => (VarType::Continuous, *lo, f64::INFINITY),
                VariableType::UpperBound(hi) => (VarType::Continuous, 0.0, *hi),
                VariableType::DoubleBound(lo, hi) => (VarType::Continuous, *lo, *hi),
                other => {
                    return Err(ReadError::UnsupportedVariable(label, other.to_string()));
                }
            };
            // A `Bounds` entry replaces the declared type, so the section lists are
            // the only surviving evidence that the column is integral.
            if integral.contains(&label) {
                kind = VarType::Integer;
            }
            if kind == VarType::Integer {
                for bound in [lo, hi] {
                    if bound.is_finite() && (bound - bound.round()).abs() > 1e-9 {
                        return Err(ReadError::FractionalIntegerBound(label, bound));
                    }
                }
            }
            col_names.push(label);
            col_type.push(kind);
            col_lb.push(lo);
            col_ub.push(hi);
        }

        // Position within the IndexMap is the column index, so a name lookup is just
        // `get_index_of` -- no side table to build or keep consistent.
        let col_of = |id| {
            lp.variables
                .get_index_of(&id)
                .expect("variable is interned and present")
        };

        let negate = lp.sense == LpSense::Maximize;
        let flip = |v: f64| if negate { -v } else { v };

        let n_cols = col_names.len();
        let mut obj = vec![0.0f64; n_cols];
        for coeff in &objective.coefficients {
            // `+=`: both formats permit a variable to appear more than once in an
            // expression, with the terms accumulating.
            obj[col_of(coeff.name)] += flip(coeff.value);
        }

        let n_rows = lp.constraints.len();
        let mut triplets: Vec<(usize, usize, f64)> = Vec::new();
        let mut row_lb = Vec::with_capacity(n_rows);
        let mut row_ub = Vec::with_capacity(n_rows);
        let mut row_names = Vec::with_capacity(n_rows);

        for (i, constraint) in lp.constraints.values().enumerate() {
            match constraint {
                Constraint::Standard {
                    name: row_name,
                    coefficients,
                    operator,
                    rhs,
                    ..
                } => {
                    let sense = match operator {
                        // Neither format has strict inequalities; some writers emit
                        // them anyway and every solver reads them as non-strict.
                        ComparisonOp::GTE | ComparisonOp::GT => RowSense::Ge,
                        ComparisonOp::LTE | ComparisonOp::LT => RowSense::Le,
                        ComparisonOp::EQ => RowSense::Eq,
                    };
                    let (lo, hi) = sense.bounds(*rhs);
                    row_lb.push(lo);
                    row_ub.push(hi);
                    row_names.push(name(*row_name));
                    for coeff in coefficients {
                        triplets.push((i, col_of(coeff.name), coeff.value));
                    }
                }
                Constraint::SOS { name: row_name, .. } => {
                    return Err(ReadError::SosConstraint(name(*row_name)));
                }
            }
        }

        Ok(Problem {
            name: lp.name.clone().unwrap_or_else(|| "unnamed".to_string()),
            sense: if negate {
                Sense::Maximize
            } else {
                Sense::Minimize
            },
            obj,
            obj_offset: flip(objective.constant),
            matrix: SparseMatrix::from_triplets(n_rows, n_cols, triplets),
            row_lb,
            row_ub,
            col_lb,
            col_ub,
            col_type,
            col_names,
            row_names,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(text: &str) -> Result<Problem, ReadError> {
        let lp = LpProblem::parse(text).map_err(|e| ReadError::Parse(e.to_string()))?;
        Problem::from_lp(&lp)
    }

    const SIMPLE: &str = "\
Minimize
 obj: 3 x1 + 2 x2 + 4 x10
Subject To
 c1: x1 + x2 >= 1
 c2: x2 + x10 <= 1
 c3: x1 + x10 = 1
Binary
 x1 x2 x10
End
";

    #[test]
    fn reads_columns_in_file_order() {
        let p = parse(SIMPLE).unwrap();
        assert_eq!(p.col_names, ["x1", "x2", "x10"]);
        assert_eq!(p.obj, [3.0, 2.0, 4.0]);
        assert!(p.validate().is_ok());
    }

    #[test]
    fn converts_each_sense_to_range_bounds() {
        let p = parse(SIMPLE).unwrap();
        let row = |name: &str| {
            let i = p.row_names.iter().position(|n| n == name).unwrap();
            (p.row_lb[i], p.row_ub[i])
        };
        assert_eq!(row("c1"), (1.0, f64::INFINITY));
        assert_eq!(row("c2"), (f64::NEG_INFINITY, 1.0));
        assert_eq!(row("c3"), (1.0, 1.0));
    }

    #[test]
    fn maximize_is_negated_and_reported_back() {
        let p =
            parse("Maximize\n obj: 3 x1\nSubject To\n c1: x1 <= 1\nBinary\n x1\nEnd\n").unwrap();
        assert_eq!(p.sense, Sense::Maximize);
        // Stored negated for the minimizing solver...
        assert_eq!(p.obj, [-3.0]);
        // ...and converted back on the way out.
        assert_eq!(p.objective_value(-3.0), 3.0);
    }

    #[test]
    fn repeated_terms_accumulate() {
        let p =
            parse("Minimize\n obj: 2 x1 + 3 x1\nSubject To\n c1: x1 + x1 >= 1\nBinary\n x1\nEnd\n")
                .unwrap();
        assert_eq!(p.obj, [5.0]);
        assert_eq!(p.matrix.column(0), (&[0usize][..], &[2.0][..]));
    }

    #[test]
    fn order_is_stable_across_repeated_parses() {
        // Column index drives branching order, so an unstable order would make node
        // counts irreproducible on identical input.
        let first = parse(SIMPLE).unwrap();
        for _ in 0..16 {
            let again = parse(SIMPLE).unwrap();
            assert_eq!(first.col_names, again.col_names);
            assert_eq!(first.row_names, again.row_names);
            assert_eq!(first.matrix, again.matrix);
        }
    }

    #[test]
    fn reads_a_general_integer_variable() {
        // No bounds section, so a General integer is [0, inf) -- not binary.
        let p = parse("Minimize\n obj: x1\nSubject To\n c1: x1 >= 1\nGeneral\n x1\nEnd\n").unwrap();
        assert!(p.is_integer(0));
        assert!(!p.is_binary(0), "an unbounded integer is not binary");
        assert_eq!((p.col_lb[0], p.col_ub[0]), (0.0, f64::INFINITY));
    }

    #[test]
    fn a_binary_declaration_is_an_integer_bounded_to_one() {
        let p = parse(SIMPLE).unwrap();
        assert!((0..p.n_cols()).all(|j| p.is_binary(j)));
        assert_eq!(p.col_ub, vec![1.0; p.n_cols()]);
    }

    #[test]
    fn rejects_a_model_with_no_objective() {
        // `Subject To` alone is not a complete LP file, so build the case directly.
        let empty = LpProblem::default();
        assert!(matches!(
            Problem::from_lp(&empty),
            Err(ReadError::NoObjective)
        ));
    }

    #[test]
    fn format_is_inferred_from_extension_case_insensitively() {
        assert_eq!(Format::from_path(Path::new("m.lp")).unwrap(), Format::Lp);
        assert_eq!(Format::from_path(Path::new("m.MPS")).unwrap(), Format::Mps);
        assert!(matches!(
            Format::from_path(Path::new("m.txt")),
            Err(ReadError::UnknownFormat(_))
        ));
    }
}
