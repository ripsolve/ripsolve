//! How many columns probing could fix, using the conflict graph and row activities.
//!
//! Tentatively fixes each free binary each way and propagates. A value whose
//! propagation proves infeasible cannot be taken, so the other one is forced. This
//! counts what that would find without changing the model.
use ripsolve::Problem;
use ripsolve::cuts::Conflicts;

fn main() {
    for path in std::env::args().skip(1) {
        let problem = match Problem::from_file(std::path::Path::new(&path)) {
            Ok(p) => p,
            Err(e) => {
                println!("{path}: {e}");
                continue;
            }
        };
        let name = std::path::Path::new(&path)
            .file_stem()
            .unwrap()
            .to_string_lossy()
            .to_string();
        let n = problem.n_cols();
        let conflicts = Conflicts::of(&problem);
        let csr = problem.matrix.to_csr();
        let binary = |j: usize| {
            problem.is_integer(j) && problem.col_lb[j] >= 0.0 && problem.col_ub[j] <= 1.0
        };

        // One propagation sweep from a single tentative fixing.
        let probe = |j: usize, value: f64| -> bool {
            let mut lb = problem.col_lb.clone();
            let mut ub = problem.col_ub.clone();
            lb[j] = value;
            ub[j] = value;
            let mut queue = vec![j];
            let mut steps = 0usize;
            while let Some(k) = queue.pop() {
                steps += 1;
                if steps > 50_000 {
                    return true;
                }
                let node = if lb[k] > 0.5 { 2 * k } else { 2 * k + 1 };
                for excluded in conflicts.adjacent(node as u32) {
                    let c = excluded as usize / 2;
                    let forced = if excluded.is_multiple_of(2) { 0.0 } else { 1.0 };
                    if lb[c] >= ub[c] - 1e-9 {
                        if (lb[c] - forced).abs() > 1e-9 {
                            return false;
                        }
                        continue;
                    }
                    lb[c] = forced;
                    ub[c] = forced;
                    queue.push(c);
                }
            }
            // Any row whose activity range no longer meets its bounds refutes the trial.
            for i in 0..problem.n_rows() {
                let (cols, vals) = csr.column(i);
                let (mut lo, mut hi) = (0.0f64, 0.0f64);
                for (&c, &a) in cols.iter().zip(vals) {
                    let (x, y) = (a * lb[c], a * ub[c]);
                    lo += x.min(y);
                    hi += x.max(y);
                }
                if lo > problem.row_ub[i] + 1e-6 || hi < problem.row_lb[i] - 1e-6 {
                    return false;
                }
            }
            true
        };

        let mut fixable = 0usize;
        let mut free = 0usize;
        for j in 0..n {
            if !binary(j) || problem.col_lb[j] >= problem.col_ub[j] {
                continue;
            }
            free += 1;
            let up = probe(j, 1.0);
            let down = probe(j, 0.0);
            if !up || !down {
                fixable += 1;
            }
        }
        println!("{name:<20} {free:>7} free binaries, probing fixes {fixable:>7}");
    }
}
