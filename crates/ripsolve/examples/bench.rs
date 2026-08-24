use ripsolve::Problem;
use ripsolve::search::{self, Options};
use std::time::Instant;

fn main() {
    let home = std::env::var("HOME").unwrap();
    for (name, path) in [
        ("v064c200", "samples/v064c200.mps".to_string()),
        (
            "v081c162n018",
            format!("{home}/repos/bip-gen/v081c162n018.lp"),
        ),
        (
            "v128c1000n100",
            format!("{home}/repos/bip-gen/v128c1000n100.lp"),
        ),
        (
            "v256c256n100",
            format!("{home}/repos/bip-gen/v256c256n100.lp"),
        ),
    ] {
        let p = Problem::from_file(std::path::Path::new(&path)).unwrap();
        println!("--- {name}");
        for (rounds, per) in [(1usize, 8usize), (2, 12), (3, 32), (5, 32)] {
            let o = Options {
                cut_rounds: rounds,
                cuts_per_round: per,
                time_limit: Some(std::time::Duration::from_secs(60)),
                ..Options::default()
            };
            let t = Instant::now();
            let s = search::solve(&p, o);
            println!(
                "  rounds={rounds} per={per:<3} cuts={:<4} root {:.1}->{:.1} {:>8} nodes {:>8.2?} {}",
                s.cuts_added,
                s.root_bound,
                s.root_bound_after_cuts,
                s.nodes,
                t.elapsed(),
                if s.status == search::Status::Optimal {
                    "optimal"
                } else {
                    "LIMIT"
                }
            );
        }
    }
}
