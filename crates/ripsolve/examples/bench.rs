use ripsolve::Problem;
use ripsolve::search::{self, Options};
use std::time::{Duration, Instant};

fn main() {
    let home = std::env::var("HOME").unwrap();
    for (name, path) in [
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
    ] {
        let p = Problem::from_file(std::path::Path::new(&path)).unwrap();
        print!("{name:16}");
        for limit in [0usize, 5, 15, 50, usize::MAX] {
            let o = Options {
                plunge_limit: limit,
                time_limit: Some(Duration::from_secs(60)),
                ..Options::default()
            };
            let t = Instant::now();
            let s = search::solve(&p, o);
            let tag = if s.status == search::Status::Optimal {
                format!("{:>6.1}s", t.elapsed().as_secs_f64())
            } else {
                format!("g{:>4.0}%", s.gap() * 100.0)
            };
            print!(
                " | {}:{:>7}n {tag}",
                if limit == usize::MAX {
                    "dfs".into()
                } else {
                    limit.to_string()
                },
                s.nodes
            );
        }
        println!();
    }
}
