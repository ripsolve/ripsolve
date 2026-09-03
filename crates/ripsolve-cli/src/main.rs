use std::path::PathBuf;

/// Which column values `solve` should print.
#[derive(Clone, Copy, Debug, PartialEq, Eq, clap::ValueEnum)]
enum Values {
    /// Print none of them.
    None,
    /// Print only the columns away from zero, which is usually what is interesting.
    Nonzero,
    /// Print every column, including the zeros.
    All,
}

use anyhow::{Context, Result};
use ripsolve::Problem;
use ripsolve::generate::{Kind, Spec};
use std::time::Duration;

use clap::{Parser, Subcommand, ValueEnum};
use ripsolve::lp::{Lp, LpStatus};
use ripsolve::search::{self, Options, Status as SearchStatus};

#[derive(Parser)]
#[command(
    name = "ripsolve",
    version,
    about = "A branch-and-cut solver for mixed-integer programs"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Read a model and report its dimensions.
    Info {
        /// Model file, `.lp` or `.mps`.
        path: PathBuf,
    },
    /// Solve a model to proven optimality.
    Solve {
        /// Model file, `.lp` or `.mps`.
        path: PathBuf,
        /// Stop after this many nodes.
        #[arg(long)]
        max_nodes: Option<usize>,
        /// Stop after this many seconds.
        #[arg(long)]
        time_limit: Option<f64>,
        /// Relative optimality gap to stop at. 0 demands a proof of optimality.
        #[arg(long, default_value_t = 1e-4)]
        gap: f64,
        /// Which column values to print: none, nonzero, or all.
        #[arg(long, value_enum, default_value_t = Values::None)]
        values: Values,
        /// Shorthand for `--values nonzero`.
        #[arg(short, long)]
        verbose: bool,
        /// Also write the solution to this file, one `name value` pair per line.
        #[arg(long, value_name = "PATH")]
        solution: Option<PathBuf>,
        /// Skip presolve.
        #[arg(long)]
        no_presolve: bool,
        /// Rounds of root cut separation. Zero disables cuts.
        #[arg(long, default_value_t = 50)]
        cut_rounds: usize,
        /// Most cuts to keep per separation round.
        #[arg(long, default_value_t = 64)]
        cuts_per_round: usize,
        /// Separate cuts at one node in every N, not only at the root. 0 disables.
        #[arg(long, default_value_t = 10)]
        local_cut_frequency: usize,
        /// Worker threads; defaults to the machine's parallelism.
        #[arg(short, long)]
        threads: Option<usize>,
        /// Flips per column the LP-free feasibility search may make. 0 disables it.
        #[arg(long, default_value_t = 200)]
        jump_moves: usize,
    },
    /// Solve a model's LP relaxation and report the bound.
    Relax {
        /// Model file, `.lp` or `.mps`.
        path: PathBuf,
        /// Maximum simplex iterations.
        #[arg(long, default_value_t = 25_000)]
        max_iterations: usize,
        /// Enter through the dual method instead of the primal.
        #[arg(long)]
        dual: bool,
    },
    /// Write a reproducible random instance in LP format.
    Gen {
        #[arg(long, value_enum)]
        kind: GenKind,
        /// Number of binary variables.
        #[arg(long)]
        cols: usize,
        /// Number of constraints.
        #[arg(long)]
        rows: usize,
        #[arg(long, default_value_t = 0)]
        seed: u64,
        /// Output file; writes to stdout when omitted.
        #[arg(short, long)]
        out: Option<PathBuf>,
    },
}

#[derive(Clone, Copy, ValueEnum)]
enum GenKind {
    Knapsack,
    Covering,
    Signed,
}

impl From<GenKind> for Kind {
    fn from(k: GenKind) -> Kind {
        match k {
            GenKind::Knapsack => Kind::Knapsack,
            GenKind::Covering => Kind::Covering,
            GenKind::Signed => Kind::Signed,
        }
    }
}

fn main() -> Result<()> {
    match Cli::parse().command {
        Command::Info { path } => {
            let problem =
                Problem::from_file(&path).with_context(|| format!("reading {}", path.display()))?;
            problem.validate()?;
            println!(
                "{}: {} columns, {} rows, {} nonzeros ({:.2}% dense), sense {:?}",
                problem.name,
                problem.n_cols(),
                problem.n_rows(),
                problem.matrix.nnz(),
                problem.matrix.density() * 100.0,
                problem.sense,
            );
        }
        Command::Solve {
            path,
            max_nodes,
            time_limit,
            gap,
            values,
            verbose,
            solution: solution_path,
            no_presolve,
            cut_rounds,
            cuts_per_round,
            local_cut_frequency,
            threads,
            jump_moves,
        } => {
            let problem =
                Problem::from_file(&path).with_context(|| format!("reading {}", path.display()))?;
            problem.validate()?;

            let options = Options {
                max_nodes: max_nodes.unwrap_or(usize::MAX),
                time_limit: time_limit.map(Duration::from_secs_f64),
                gap_tolerance: gap,
                presolve: !no_presolve,
                threads: threads
                    .unwrap_or_else(|| std::thread::available_parallelism().map_or(1, |n| n.get())),
                cut_rounds,
                cuts_per_round,
                local_cut_frequency,
                jump_moves,
                ..Options::default()
            };
            let started = std::time::Instant::now();
            let solution = search::solve(&problem, options);
            let elapsed = started.elapsed();

            match solution.objective {
                Some(objective) => {
                    println!("objective: {objective}");
                    // Only a search that ran to exhaustion has proven anything; a
                    // limit-terminated run reports its remaining gap instead.
                    if solution.status == SearchStatus::Optimal {
                        println!("status:    optimal");
                    } else if solution.bound.is_finite() {
                        println!(
                            "status:    {:?} (bound {}, gap {:.4}%)",
                            solution.status,
                            solution.bound,
                            solution.gap() * 100.0
                        );
                    } else {
                        // Nothing was proven about this point, and saying so in words
                        // beats printing a bound of NaN and leaving the reader to work
                        // out what that means.
                        println!("status:    {:?} (no bound proven)", solution.status);
                    }
                    // `-v` is the old spelling of `--values nonzero` and still works.
                    let wanted = if verbose && values == Values::None {
                        Values::Nonzero
                    } else {
                        values
                    };
                    if wanted != Values::None {
                        // One per line rather than one long line: a model with
                        // thousands of columns is unreadable otherwise, and a line per
                        // column is what `grep` and `awk` expect.
                        let width = problem
                            .col_names
                            .iter()
                            .map(|n| n.len())
                            .max()
                            .unwrap_or(0)
                            .min(40);
                        println!("values:");
                        for (name, value) in problem.col_names.iter().zip(&solution.x) {
                            // A general integer can be 3, not just 0 or 1, so the value
                            // is printed rather than the name alone.
                            if wanted == Values::All || value.abs() > 1e-9 {
                                println!("  {name:<width$} {value}");
                            }
                        }
                    }
                    if let Some(path) = &solution_path {
                        let mut out = String::new();
                        for (name, value) in problem.col_names.iter().zip(&solution.x) {
                            out.push_str(&format!("{name} {value}\n"));
                        }
                        std::fs::write(path, out)
                            .with_context(|| format!("writing {}", path.display()))?;
                    }
                }
                None => println!("status:    {:?}", solution.status),
            }
            if let Some(stats) = solution.presolve {
                println!(
                    "presolve:  {} columns fixed, {} rows removed, {} coefficients tightened",
                    stats.fixed_columns, stats.redundant_rows, stats.tightened_coefficients
                );
            }
            if solution.cuts_added > 0 {
                println!(
                    "cuts:      {} added, root bound {:.6} -> {:.6}",
                    solution.cuts_added, solution.root_bound, solution.root_bound_after_cuts
                );
            }
            println!("heuristic: {} incumbents", solution.heuristic_solutions);
            println!(
                "{} nodes, {} simplex iterations, {elapsed:.3?}",
                solution.nodes, solution.simplex_iterations
            );
        }
        Command::Relax {
            path,
            max_iterations,
            dual,
        } => {
            let problem =
                Problem::from_file(&path).with_context(|| format!("reading {}", path.display()))?;
            problem.validate()?;

            let started = std::time::Instant::now();
            let mut lp = Lp::relaxation(&problem);
            let solution = if dual {
                lp.solve_cold_dual(max_iterations)
            } else {
                lp.solve_with_limit(max_iterations)
            };
            let elapsed = started.elapsed();

            match solution.status {
                LpStatus::Optimal => {
                    // Reported in the user's original sense, not the internal
                    // minimization form the simplex works in.
                    println!(
                        "relaxation: {:.9}",
                        problem.objective_value(solution.objective)
                    );
                    // How far the relaxation is from integral is the first thing that
                    // matters once branching exists: zero fractional columns means the
                    // node is already integer-feasible and needs no branching at all.
                    let fractional = solution
                        .x
                        .iter()
                        .filter(|v| {
                            let f = *v - v.floor();
                            f > 1e-7 && f < 1.0 - 1e-7
                        })
                        .count();
                    println!("{fractional} of {} columns fractional", problem.n_cols());
                }
                other => println!("relaxation: {other:?}"),
            }
            println!(
                "{} simplex iterations in {elapsed:.3?}",
                solution.iterations
            );
        }
        Command::Gen {
            kind,
            cols,
            rows,
            seed,
            out,
        } => {
            let spec = Spec {
                kind: kind.into(),
                n_cols: cols,
                n_rows: rows,
                seed,
            };
            let text = spec.to_lp();
            match out {
                Some(path) => std::fs::write(&path, text)
                    .with_context(|| format!("writing {}", path.display()))?,
                None => print!("{text}"),
            }
        }
    }
    Ok(())
}
