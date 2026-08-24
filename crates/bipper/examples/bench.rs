use bipper::Problem;
use bipper::generate::{Kind, Spec};
use bipper::search::{self, Options};
use lp_parser_rs::problem::LpProblem;
use std::time::Instant;

fn main() {
    let specs = [
        (Kind::Knapsack, 20, 10, 1u64),
        (Kind::Knapsack, 30, 15, 2),
        (Kind::Knapsack, 45, 20, 3),
        (Kind::Covering, 25, 30, 1),
        (Kind::Covering, 40, 50, 2),
        (Kind::Covering, 60, 80, 3),
        (Kind::Signed, 20, 20, 1),
        (Kind::Signed, 32, 32, 2),
        (Kind::Signed, 48, 48, 3),
    ];
    for (kind, c, r, seed) in specs {
        let spec = Spec {
            kind,
            n_cols: c,
            n_rows: r,
            seed,
        };
        let p = Problem::from_lp(&LpProblem::parse(&spec.to_lp()).unwrap()).unwrap();
        let t = Instant::now();
        let s = search::solve(&p, Options::default());
        println!(
            "{:22} obj={:9.2} nodes={:8} simplex={:9} {:>11.3?}",
            spec.name(),
            s.objective.unwrap_or(f64::NAN),
            s.nodes,
            s.simplex_iterations,
            t.elapsed()
        );
    }
    for f in ["v048c048.lp", "v064c064.lp", "v064c200.mps"] {
        let p = Problem::from_file(std::path::Path::new("samples").join(f).as_path()).unwrap();
        let t = Instant::now();
        let s = search::solve(&p, Options::default());
        println!(
            "{:22} obj={:9.2} nodes={:8} simplex={:9} {:>11.3?}",
            f,
            s.objective.unwrap_or(f64::NAN),
            s.nodes,
            s.simplex_iterations,
            t.elapsed()
        );
    }
}
