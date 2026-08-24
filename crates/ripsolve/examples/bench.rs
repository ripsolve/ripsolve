use ripsolve::Problem;
use ripsolve::search::{self, Options};
use std::time::Instant;

fn main() {
    let home = std::env::var("HOME").unwrap();
    let cases: Vec<(String, String)> = vec![
        (
            "v081c162n009",
            format!("{home}/repos/bip-gen/v081c162n009.lp"),
        ),
        (
            "v081c162n018",
            format!("{home}/repos/bip-gen/v081c162n018.lp"),
        ),
        (
            "v128c256n100",
            format!("{home}/repos/bip-gen/v128c256n100.lp"),
        ),
        (
            "v256c256n100",
            format!("{home}/repos/bip-gen/v256c256n100.lp"),
        ),
        ("v064c064", "samples/v064c064.lp".into()),
        ("v064c200", "samples/v064c200.mps".into()),
        ("v048c048", "samples/v048c048.lp".into()),
    ]
    .into_iter()
    .map(|(a, b)| (a.to_string(), b))
    .collect();

    println!(
        "{:16} {:>9} {:>10}   objective",
        "instance", "nodes", "time"
    );
    for (name, path) in cases {
        let p = Problem::from_file(std::path::Path::new(&path)).unwrap();
        let t = Instant::now();
        let s = search::solve(&p, Options::default());
        println!(
            "{name:16} {:>9} {:>10.2?}   {}",
            s.nodes,
            t.elapsed(),
            s.objective.unwrap_or(f64::NAN)
        );
    }
}
