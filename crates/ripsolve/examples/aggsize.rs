//! What doubleton aggregation removes from a model, before anything else runs.
use ripsolve::Problem;
use ripsolve::aggregate;

fn main() {
    println!("{:<24} {:>8} {:>8} {:>9}   {:>8} {:>8} {:>9}", "model", "rows", "cols", "nnz", "rows", "cols", "nnz");
    for path in std::env::args().skip(1) {
        let name = std::path::Path::new(&path)
            .file_stem()
            .map_or_else(String::new, |s| s.to_string_lossy().into_owned());
        let problem = Problem::from_file(std::path::Path::new(&path)).unwrap();
        match aggregate::aggregate(&problem, 1e-9) {
            Err(_) => println!("{name:<24} proved infeasible by aggregation"),
            Ok(None) => println!(
                "{name:<24} {:>8} {:>8} {:>9}   unchanged",
                problem.n_rows(), problem.n_cols(), problem.matrix.nnz()
            ),
            Ok(Some((smaller, map))) => println!(
                "{name:<24} {:>8} {:>8} {:>9} ->{:>8} {:>8} {:>9}   {} eliminated",
                problem.n_rows(), problem.n_cols(), problem.matrix.nnz(),
                smaller.n_rows(), smaller.n_cols(), smaller.matrix.nnz(), map.eliminated()
            ),
        }
    }
}
