use bipper::Problem;
use bipper::generate::{Kind, Spec};
use bipper::search::{self, Options};
use lp_parser_rs::problem::LpProblem;
use std::time::Instant;

fn run(p: &Problem, label: &str) {
    let mut line = format!("{label:22}");
    for presolve in [false, true] {
        let options = Options {
            presolve,
            ..Options::default()
        };
        let t = Instant::now();
        let s = search::solve(p, options);
        line.push_str(&format!(" | {:>8} nodes {:>9.3?}", s.nodes, t.elapsed()));
        if presolve {
            let st = s.presolve.unwrap_or_default();
            line.push_str(&format!(
                "  fixed={:<4} rows={:<4} coef={:<5} obj={}",
                st.fixed_columns,
                st.redundant_rows,
                st.tightened_coefficients,
                s.objective.unwrap_or(f64::NAN)
            ));
        }
    }
    println!("{line}");
}

fn main() {
    println!(
        "{:22} | {:^28} | {:^28}",
        "", "no presolve", "with presolve"
    );
    for (kind, c, r, seed) in [
        (Kind::Knapsack, 30, 15, 2u64),
        (Kind::Knapsack, 45, 20, 3),
        (Kind::Covering, 40, 50, 2),
        (Kind::Covering, 60, 80, 3),
        (Kind::Signed, 32, 32, 2),
        (Kind::Signed, 48, 48, 3),
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
