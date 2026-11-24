//! Solves the reduced model and reports how far the point the map gives back is from
//! feasible in the original, row by row, to find a postsolve record that loses it.
use ripsolve::{Problem, reduce, search};

fn worst(problem: &Problem, x: &[f64]) -> Vec<(f64, String)> {
    let csr = problem.matrix.to_csr();
    let mut out = Vec::new();
    for (j, &v) in x.iter().enumerate() {
        let over = (problem.col_lb[j] - v).max(v - problem.col_ub[j]);
        if over > 1e-6 {
            out.push((over, format!("col {j} = {v} in [{}, {}]", problem.col_lb[j], problem.col_ub[j])));
        }
        if problem.is_integer(j) && (v - v.round()).abs() > 1e-6 {
            out.push(((v - v.round()).abs(), format!("col {j} = {v} fractional")));
        }
    }
    for i in 0..problem.n_rows() {
        let (cols, vals) = csr.column(i);
        let activity: f64 = cols.iter().zip(vals).map(|(&j, &a)| a * x[j]).sum();
        let over = (problem.row_lb[i] - activity).max(activity - problem.row_ub[i]);
        if over > 1e-6 {
            let terms: Vec<String> = cols.iter().zip(vals).map(|(&j, &a)| format!("{a}*x{j}({})", x[j])).collect();
            out.push((over, format!("row {i} activity {activity} in [{}, {}]: {}", problem.row_lb[i], problem.row_ub[i], terms.join(" "))));
        }
    }
    out.sort_by(|a, b| b.0.total_cmp(&a.0));
    out
}

fn main() {
    let mut args = std::env::args().skip(1);
    let input = args.next().unwrap();
    let seconds: f64 = args.next().map_or(20.0, |s| s.parse().unwrap());
    let original = Problem::from_file(std::path::Path::new(&input)).unwrap();
    // A solution file from `ripsolve solve --solution` in place of a time: check it.
    if let Some(path) = args.next() {
        let text = std::fs::read_to_string(path).unwrap();
        let values: std::collections::HashMap<&str, f64> =
            text.lines().filter_map(|l| l.split_once(' ')).map(|(n, v)| (n, v.parse().unwrap())).collect();
        let x: Vec<f64> = original.col_names.iter().map(|n| values[n.as_str()]).collect();
        let score: f64 = original.obj.iter().zip(&x).map(|(c, v)| c * v).sum();
        println!("objective {}", original.objective_value(score));
        for (over, what) in worst(&original, &x).iter().take(15) {
            println!("{over:.3e} {what}");
        }
        return;
    }
    let Ok(Some((smaller, map, stats))) = reduce::reduce_with_presolve(&original, 1e-9) else {
        println!("nothing reduced");
        return;
    };
    println!("{stats:?}");
    let options = search::Options {
        aggregate: false,
        time_limit: Some(std::time::Duration::from_secs_f64(seconds)),
        ..Default::default()
    };
    let solution = search::solve(&smaller, options);
    let Some(objective) = solution.objective else {
        println!("no solution: {:?}", solution.status);
        return;
    };
    println!("reduced objective {objective}, violations in reduced: {:?}", worst(&smaller, &solution.x).iter().take(3).collect::<Vec<_>>());
    let x = map.expand(&solution.x);
    let score: f64 = original.obj.iter().zip(&x).map(|(c, v)| c * v).sum();
    println!("expanded objective {}", original.objective_value(score));
    for (over, what) in worst(&original, &x).iter().take(15) {
        println!("{over:.3e} {what}");
    }
}
