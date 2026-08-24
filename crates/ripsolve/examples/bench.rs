use ripsolve::Problem;
use ripsolve::search::{self, Options};
use std::time::{Duration, Instant};

fn main() {
    let home = std::env::var("HOME").unwrap();
    let cases = [
        ("v048c048", "samples/v048c048.lp".to_string()),
        ("v064c064", "samples/v064c064.lp".to_string()),
        ("v064c200", "samples/v064c200.mps".to_string()),
        (
            "v081c162n009",
            format!("{home}/repos/bip-gen/v081c162n009.lp"),
        ),
        (
            "v081c162n018",
            format!("{home}/repos/bip-gen/v081c162n018.lp"),
        ),
        (
            "v256c256n100",
            format!("{home}/repos/bip-gen/v256c256n100.lp"),
        ),
        (
            "v064c1000n100",
            format!("{home}/repos/bip-gen/v064c1000n100.lp"),
        ),
        (
            "v128c1000n100",
            format!("{home}/repos/bip-gen/v128c1000n100.lp"),
        ),
    ];
    print!("{:16}", "instance");
    for b in ["off", "100", "1000", "10000"] {
        print!(" | {:>18}", format!("budget {b}"));
    }
    println!();
    for (name, path) in cases {
        let p = Problem::from_file(std::path::Path::new(&path)).unwrap();
        print!("{name:16}");
        for budget in [0usize, 100, 1000, 10000] {
            let o = Options {
                strong_branching_budget: budget,
                threads: 1,
                time_limit: Some(Duration::from_secs(120)),
                ..Options::default()
            };
            let t = Instant::now();
            let s = search::solve(&p, o);
            let tag = if s.status == search::Status::Optimal {
                format!("{:>7}n {:>6.2}s", s.nodes, t.elapsed().as_secs_f64())
            } else {
                format!("{:>7}n g{:>4.0}%", s.nodes, s.gap() * 100.0)
            };
            print!(" | {tag:>18}");
        }
        println!();
    }
}
