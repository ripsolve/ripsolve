//! The search must reproduce the reference MIP optima.

mod fixtures;

use fixtures::{fixtures, samples_dir, spec_of};
use lp_parser_rs::problem::LpProblem;
use ripsolve::Problem;
use ripsolve::search::{self, Options, Status};

fn check(problem: &Problem, expected: f64, name: &str) -> search::Solution {
    let solution = search::solve(problem, Options::default());
    assert_eq!(solution.status, Status::Optimal, "{name}");
    let got = solution.objective.unwrap_or(f64::NAN);
    let scale = expected.abs().max(1.0);
    assert!(
        (got - expected).abs() <= 1e-6 * scale,
        "{name}: optimum {got}, expected {expected} ({} nodes)",
        solution.nodes
    );
    // A proven-optimal search must close the gap.
    assert!(
        solution.gap() <= 1e-6,
        "{name}: gap {} not closed",
        solution.gap()
    );
    solution
}

/// Check the reported assignment is genuinely feasible and attains the objective.
fn assert_solution_is_valid(problem: &Problem, x: &[f64], objective: f64, name: &str) {
    assert_eq!(x.len(), problem.n_cols(), "{name}: solution length");
    let values: Vec<f64> = x.to_vec();
    // Integer columns must actually be integral.
    for (j, &v) in values.iter().enumerate() {
        if problem.is_integer(j) {
            assert!(
                (v - v.round()).abs() < 1e-6,
                "{name}: integer column {j} is {v}"
            );
        }
    }

    let csr = problem.matrix.to_csr();
    for i in 0..problem.n_rows() {
        let (cols, vals) = csr.column(i);
        let activity: f64 = cols.iter().zip(vals).map(|(&j, &a)| a * values[j]).sum();
        assert!(
            activity >= problem.row_lb[i] - 1e-6 && activity <= problem.row_ub[i] + 1e-6,
            "{name}: row {i} activity {activity} outside [{}, {}]",
            problem.row_lb[i],
            problem.row_ub[i]
        );
    }

    let internal: f64 = problem.obj.iter().zip(&values).map(|(c, v)| c * v).sum();
    let recomputed = problem.objective_value(internal);
    assert!(
        (recomputed - objective).abs() < 1e-6,
        "{name}: reported {objective} but the assignment scores {recomputed}"
    );
}

#[test]
fn generated_optima_match_the_reference_oracle() {
    let data = fixtures();
    for entry in data["instances"].as_array().unwrap() {
        let spec = spec_of(entry);
        let name = entry["name"].as_str().unwrap();
        let problem = Problem::from_lp(&LpProblem::parse(&spec.to_lp()).unwrap()).unwrap();
        let expected = entry["mip_optimum"].as_f64().unwrap();
        let solution = check(&problem, expected, name);
        assert_solution_is_valid(&problem, &solution.x, solution.objective.unwrap(), name);
    }
}

#[test]
fn sample_optima_match_the_reference_oracle() {
    let data = fixtures();
    for entry in data["samples"].as_array().unwrap() {
        let file = entry["file"].as_str().unwrap();
        let problem = Problem::from_file(&samples_dir().join(file)).unwrap();
        let expected = entry["mip_optimum"].as_f64().unwrap();
        let solution = check(&problem, expected, file);
        assert_solution_is_valid(&problem, &solution.x, solution.objective.unwrap(), file);
    }
}

#[test]
fn reports_infeasibility() {
    let text = "Minimize\n obj: x0 + x1\nSubject To\n c0: x0 + x1 >= 1\n \
                c1: x0 + x1 <= 0\nBinary\n x0 x1\nEnd\n";
    let problem = Problem::from_lp(&LpProblem::parse(text).unwrap()).unwrap();
    let solution = search::solve(&problem, Options::default());
    assert_eq!(solution.status, Status::Infeasible);
    assert!(solution.objective.is_none());
}

#[test]
fn a_node_limit_stops_early_without_claiming_optimality() {
    let spec = ripsolve::generate::Spec {
        kind: ripsolve::generate::Kind::Knapsack,
        n_cols: 45,
        n_rows: 20,
        seed: 3,
    };
    let problem = Problem::from_lp(&LpProblem::parse(&spec.to_lp()).unwrap()).unwrap();
    let options = Options {
        max_nodes: 5,
        ..Options::default()
    };
    let solution = search::solve(&problem, options);
    assert_eq!(solution.status, Status::NodeLimit);
    assert!(solution.nodes <= 6, "{} nodes", solution.nodes);
}

#[test]
fn presolve_does_not_change_any_optimum() {
    // Presolve is allowed to reduce the model however it likes, provided the answer
    // is identical. Checked on every sample, both ways.
    let data = fixtures();
    for entry in data["samples"].as_array().unwrap() {
        let file = entry["file"].as_str().unwrap();
        let problem = Problem::from_file(&samples_dir().join(file)).unwrap();

        let with = search::solve(
            &problem,
            Options {
                presolve: true,
                ..Options::default()
            },
        );
        let without = search::solve(
            &problem,
            Options {
                presolve: false,
                ..Options::default()
            },
        );

        assert_eq!(with.status, without.status, "{file}");
        match (with.objective, without.objective) {
            (Some(a), Some(b)) => {
                assert!((a - b).abs() < 1e-6, "{file}: presolve {a}, plain {b}")
            }
            (None, None) => {}
            (a, b) => panic!("{file}: presolve {a:?}, plain {b:?}"),
        }
    }
}

#[test]
fn cuts_do_not_change_any_optimum() {
    // Cuts may only remove fractional points. If one ever removes an integer
    // solution the answer changes silently, so every sample is solved both ways.
    let data = fixtures();
    for entry in data["samples"].as_array().unwrap() {
        let file = entry["file"].as_str().unwrap();
        let problem = Problem::from_file(&samples_dir().join(file)).unwrap();

        let with = search::solve(&problem, Options::default());
        let without = search::solve(
            &problem,
            Options {
                cut_rounds: 0,
                ..Options::default()
            },
        );

        assert_eq!(with.status, without.status, "{file}");
        match (with.objective, without.objective) {
            (Some(a), Some(b)) => assert!((a - b).abs() < 1e-6, "{file}: cuts {a}, plain {b}"),
            (None, None) => {}
            (a, b) => panic!("{file}: cuts {a:?}, plain {b:?}"),
        }
        // Cuts can only tighten the root relaxation, never weaken it. "Tighter"
        // means closer to the optimum, which is a *higher* bound when minimizing and
        // a *lower* one when maximizing -- so the direction-free statement is that
        // the gap to the optimum must not widen.
        if let Some(optimum) = with.objective {
            let gap_before = (optimum - with.root_bound).abs();
            let gap_after = (optimum - with.root_bound_after_cuts).abs();
            assert!(
                gap_after <= gap_before + 1e-6,
                "{file}: cuts widened the root gap from {gap_before} to {gap_after} \
                 (bound {} -> {}, optimum {optimum})",
                with.root_bound,
                with.root_bound_after_cuts
            );
        }
    }
}

#[test]
fn parallel_search_reaches_the_same_answer() {
    // The one property a parallel solver must have. Which node a worker takes
    // depends on timing, so node counts vary run to run; the proven optimum must
    // not, because every bound and cut is globally valid and every worker prunes
    // against the same shared incumbent.
    let data = fixtures();
    for entry in data["samples"].as_array().unwrap() {
        let file = entry["file"].as_str().unwrap();
        let problem = Problem::from_file(&samples_dir().join(file)).unwrap();
        let expected = entry["mip_optimum"].as_f64().unwrap();

        for threads in [2usize, 4, 8] {
            let solution = search::solve(
                &problem,
                Options {
                    threads,
                    ..Options::default()
                },
            );
            assert_eq!(
                solution.status,
                Status::Optimal,
                "{file} on {threads} threads"
            );
            let got = solution.objective.unwrap_or(f64::NAN);
            let scale = expected.abs().max(1.0);
            assert!(
                (got - expected).abs() <= 1e-6 * scale,
                "{file} on {threads} threads: got {got}, expected {expected}"
            );
            assert!(
                solution.gap() <= 1e-6,
                "{file} on {threads} threads: gap {} not closed",
                solution.gap()
            );
        }
    }
}

#[test]
fn parallel_search_respects_a_node_limit() {
    // A limit must stop every worker, not just the one that noticed it.
    let spec = ripsolve::generate::Spec {
        kind: ripsolve::generate::Kind::Knapsack,
        n_cols: 45,
        n_rows: 20,
        seed: 3,
    };
    let problem = Problem::from_lp(&LpProblem::parse(&spec.to_lp()).unwrap()).unwrap();
    let solution = search::solve(
        &problem,
        Options {
            threads: 4,
            max_nodes: 20,
            heuristic_frequency: 0,
            ..Options::default()
        },
    );
    assert_ne!(solution.status, Status::Optimal);
    // Workers in flight may each finish their node, so allow slack over the limit
    // but not an unbounded overshoot.
    assert!(
        solution.nodes <= 40,
        "ran {} nodes past a limit of 20",
        solution.nodes
    );
}
