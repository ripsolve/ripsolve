//! A branch-and-cut solver for mixed-integer programs, in pure Rust.
//!
//! Columns may be binary, general integer, or continuous. Every node of the search
//! solves a bounded LP relaxation with the simplex method, and the dual bound that
//! comes back, strengthened by presolve and cutting planes, is what prunes the tree.
//!
//! The solver targets small and medium models, up to roughly a thousand rows. In that
//! range it is competitive with the open-source solvers and with commercial ones.
//! Above it, reach for something else.
//!
//! # Solving a model from a file
//!
//! ```no_run
//! use ripsolve::{Problem, search};
//! use std::path::Path;
//!
//! let problem = Problem::from_file(Path::new("model.lp"))?;
//! problem.validate()?;
//!
//! let solution = search::solve(&problem, search::Options::default());
//! println!("{:?} {:?}", solution.status, solution.objective);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! Both LP and MPS are read, chosen by file extension, or explicitly with
//! [`Problem::from_file_as`].
//!
//! # Building a model directly
//!
//! [`Problem`] is a plain struct with public fields, so a model is assembled rather
//! than built through a builder. Rows are held in range form, `lb <= a'x <= ub`, which
//! is what [`RowSense::bounds`] converts a written sense into. Objective coefficients
//! are always stored in minimization form; a maximization negates them on the way in,
//! and [`Problem::objective_value`] converts a solved value back.
//!
//! Maximize `3b + 2n` subject to `2b + n <= 12`, with `b` binary and `n` integer in
//! `[0, 10]`:
//!
//! ```
//! use ripsolve::model::{Problem, RowSense, Sense, VarType};
//! use ripsolve::{SparseMatrix, search};
//!
//! let (lb, ub) = RowSense::Le.bounds(12.0);
//! let problem = Problem {
//!     name: "example".into(),
//!     sense: Sense::Maximize,
//!     // Minimization form: the maximization objective negated.
//!     obj: vec![-3.0, -2.0],
//!     obj_offset: 0.0,
//!     matrix: SparseMatrix::from_triplets(1, 2, [(0, 0, 2.0), (0, 1, 1.0)]),
//!     row_lb: vec![lb],
//!     row_ub: vec![ub],
//!     col_lb: vec![0.0, 0.0],
//!     col_ub: vec![1.0, 10.0],
//!     col_type: vec![VarType::Integer, VarType::Integer],
//!     col_names: vec!["b".into(), "n".into()],
//!     row_names: vec!["c0".into()],
//! };
//! problem.validate().unwrap();
//!
//! let solution = search::solve(&problem, search::Options::default());
//! assert_eq!(solution.objective, Some(23.0));
//! ```
//!
//! A binary column is an integer column bounded to `[0, 1]`. Nothing in the solver
//! treats that as a distinct case, because branching splits a range and degenerates to
//! fixing at 0 or 1 on its own.
//!
//! # Controlling the search
//!
//! [`search::Options`] carries the limits and the tuning. The library defaults to a
//! single thread so that behaviour is reproducible; the command-line application
//! defaults to the machine's parallelism.
//!
//! ```no_run
//! use ripsolve::search::Options;
//! use std::time::Duration;
//!
//! let options = Options {
//!     threads: 8,
//!     time_limit: Some(Duration::from_secs(60)),
//!     gap_tolerance: 0.01,
//!     ..Options::default()
//! };
//! ```
//!
//! A run that stops early reports [`search::Status::NodeLimit`] or
//! [`search::Status::TimeLimit`] with the best solution found and the remaining gap,
//! rather than failing. [`search::Status::Optimal`] is the only status that claims
//! proof.
//!
//! # Where the pieces live
//!
//! [`model`] holds the problem representation, [`reader`] the LP and MPS readers, and
//! [`search`] the branch-and-cut driver. [`lp`] is the simplex: [`lp::Lp`] is a
//! relaxation and the solves against it. [`presolve`], [`cuts`], [`branch`] and
//! [`heuristic`] are the components the search drives, exposed because they are useful
//! to measure on their own.
//!
//! The reasoning behind the defaults, including the approaches that were measured and
//! rejected, is in [design-notes.md](https://github.com/ripsolve/ripsolve/blob/main/docs/design-notes.md).

pub mod branch;
pub mod cuts;
pub mod generate;
pub mod heuristic;
pub mod lp;
pub mod model;
pub mod presolve;
pub mod reader;
pub mod search;
pub mod sparse;

pub use model::{ModelError, Problem, RowSense, Sense};
pub use reader::{Format, ReadError};
pub use sparse::SparseMatrix;
