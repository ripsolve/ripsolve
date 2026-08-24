use std::path::PathBuf;

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
    about = "A branch-and-cut solver for binary integer programs"
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
        /// Stop once the relative optimality gap reaches this.
        #[arg(long, default_value_t = 0.0)]
        gap: f64,
        /// Print the value of every column.
        #[arg(short, long)]
        verbose: bool,
        /// Skip presolve.
        #[arg(long)]
        no_presolve: bool,
    },
    /// Solve a model's LP relaxation and report the bound.
    Relax {
        /// Model file, `.lp` or `.mps`.
        path: PathBuf,
        /// Maximum simplex iterations.
        #[arg(long, default_value_t = 200_000)]
        max_iterations: usize,
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
            verbose,
            no_presolve,
        } => {
            let problem =
                Problem::from_file(&path).with_context(|| format!("reading {}", path.display()))?;
            problem.validate()?;

            let options = Options {
                max_nodes: max_nodes.unwrap_or(usize::MAX),
                time_limit: time_limit.map(Duration::from_secs_f64),
                gap_tolerance: gap,
                presolve: !no_presolve,
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
                    } else {
                        println!(
                            "status:    {:?} (bound {}, gap {:.4}%)",
                            solution.status,
                            solution.bound,
                            solution.gap() * 100.0
                        );
                    }
                    if verbose {
                        let ones: Vec<&str> = problem
                            .col_names
                            .iter()
                            .zip(&solution.x)
                            .filter(|(_, v)| **v == 1)
                            .map(|(n, _)| n.as_str())
                            .collect();
                        println!("1: {}", ones.join(" "));
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
            println!(
                "{} nodes, {} simplex iterations, {elapsed:.3?}",
                solution.nodes, solution.simplex_iterations
            );
        }
        Command::Relax {
            path,
            max_iterations,
        } => {
            let problem =
                Problem::from_file(&path).with_context(|| format!("reading {}", path.display()))?;
            problem.validate()?;

            let started = std::time::Instant::now();
            let solution = Lp::relaxation(&problem).solve_with_limit(max_iterations);
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
