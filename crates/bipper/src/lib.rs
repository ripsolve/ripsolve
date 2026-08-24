//! bipper — a branch-and-cut solver for binary integer programs.
//!
//! Every variable is binary. The LP relaxations solved inside the search are
//! ordinary bounded continuous LPs, but the model exposed here is pure BIP: see
//! [`Problem`], which rejects anything else at construction.

pub mod generate;
pub mod model;
pub mod reader;
pub mod sparse;

pub use model::{ModelError, Problem, RowSense, Sense};
pub use reader::{Format, ReadError};
pub use sparse::SparseMatrix;
