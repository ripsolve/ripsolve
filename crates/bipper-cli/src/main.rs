use std::path::PathBuf;

use anyhow::{Context, Result};
use bipper::Problem;
use bipper::generate::{Kind, Spec};
use bipper::lp::{Lp, LpStatus};
use clap::{Parser, Subcommand, ValueEnum};

#[derive(Parser)]
#[command(
    name = "bipper",
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
