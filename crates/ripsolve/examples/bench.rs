use lp_parser_rs::problem::LpProblem;
use ripsolve::Problem;
use ripsolve::generate::{Kind, Spec};
use ripsolve::search::{self, Options};
use std::time::Instant;

fn run(p: &Problem, label: &str) {
    let mut line = format!("{label:22}");
    for cuts in [0usize, 20] {
        let options = Options {
            cut_rounds: cuts,
            ..Options::default()
        };
        let t = Instant::now();
        let s = search::solve(p, options);
        line.push_str(&format!(" |{:>8} nodes {:>9.3?}", s.nodes, t.elapsed()));
        if cuts > 0 {
            line.push_str(&format!(
                " {:>3} cuts  root {:.2}->{:.2} opt {}",
                s.cuts_added,
                s.root_bound,
                s.root_bound_after_cuts,
                s.objective.unwrap_or(f64::NAN)
            ));
        }
    }
    println!("{line}");
}

fn main() {
    println!("{:22} |{:^26}|{:^26}", "", "no cuts", "with cuts");
    for (kind, c, r, seed) in [
        (Kind::Knapsack, 30, 15, 2u64),
        (Kind::Knapsack, 45, 20, 3),
        (Kind::Covering, 40, 50, 2),
        (Kind::Covering, 60, 80, 3),
        (Kind::Signed, 32, 32, 2),
    ] {
        let spec = Spec {
            kind,
            n_cols: c,
            n_rows: r,
            seed,
        };
        let p = Problem::from_lp(&LpProblem::parse(&spec.to_lp()).unwrap()).unwrap();
        run(&p, &spec.name());
    }
    for f in ["v048c048.lp", "v064c064.lp", "v064c200.mps"] {
        let p = Problem::from_file(std::path::Path::new("samples").join(f).as_path()).unwrap();
        run(&p, f);
    }
}
