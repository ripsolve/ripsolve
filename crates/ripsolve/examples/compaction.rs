//! Report what presolve and compaction together remove from a model.
use ripsolve::Problem;
use ripsolve::presolve::{self, Outcome};

fn main() {
    for path in std::env::args().skip(1) {
        let name = std::path::Path::new(&path)
            .file_stem()
            .map_or_else(|| path.clone(), |s| s.to_string_lossy().into_owned());
        let Ok(original) = Problem::from_file(std::path::Path::new(&path)) else {
            continue;
        };
        let mut reduced = original.clone();
        if presolve::presolve(&mut reduced, 20) == Outcome::Infeasible {
            println!("{name}\tinfeasible");
            continue;
        }
        match presolve::compact(&reduced) {
            Ok((small, _)) => println!(
                "{name}\t{} -> {} rows\t{} -> {} cols\t{} -> {} nnz",
                original.n_rows(), small.n_rows(),
                original.n_cols(), small.n_cols(),
                original.matrix.nnz(), small.matrix.nnz(),
            ),
            Err(_) => println!("{name}\tinfeasible on compaction"),
        }
    }
}
