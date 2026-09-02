//! The simplex must reproduce the reference LP relaxation values.
//!
//! This is the correctness bar for the LP core: for every generated instance and
//! bundled sample, the relaxation solved here must agree with the value an
//! independent solver recorded, in the user's original sense.

mod fixtures;

use std::path::Path;

use fixtures::{fixtures, samples_dir, spec_of};
use lp_parser_rs::problem::LpProblem;
use ripsolve::Problem;
use ripsolve::lp::{Lp, LpStatus};

/// Solve a problem's relaxation and return its value in the original sense.
fn relaxation_value(problem: &Problem, name: &str) -> f64 {
    let solution = Lp::relaxation(problem).solve();
    assert_eq!(
        solution.status,
        LpStatus::Optimal,
        "{name}: {:?}",
        solution.status
    );
    problem.objective_value(solution.objective)
}

fn assert_matches(got: f64, expected: f64, name: &str) {
    // Relative tolerance, since sample optima span 0 to ~1800.
    let scale = expected.abs().max(1.0);
    assert!(
        (got - expected).abs() <= 1e-6 * scale,
        "{name}: relaxation {got}, expected {expected} (diff {:.3e})",
        (got - expected).abs()
    );
}

#[test]
fn generated_relaxations_match_the_reference_oracle() {
    let data = fixtures();
    let instances = data["instances"].as_array().unwrap();
    assert!(!instances.is_empty());

    for entry in instances {
        let spec = spec_of(entry);
        let name = entry["name"].as_str().unwrap();
        let parsed = LpProblem::parse(&spec.to_lp()).unwrap();
        let problem = Problem::from_lp(&parsed).unwrap();
        assert_matches(
            relaxation_value(&problem, name),
            entry["lp_relaxation"].as_f64().unwrap(),
            name,
        );
    }
}

#[test]
fn sample_relaxations_match_the_reference_oracle() {
    let data = fixtures();
    let samples = data["samples"].as_array().unwrap();
    assert!(!samples.is_empty());

    for entry in samples {
        let file = entry["file"].as_str().unwrap();
        let problem = Problem::from_file(&samples_dir().join(file)).unwrap();
        assert_matches(
            relaxation_value(&problem, file),
            entry["lp_relaxation"].as_f64().unwrap(),
            file,
        );
    }
}

#[test]
fn the_dual_cold_entry_reaches_the_same_relaxation_value() {
    // The dual entry is a fallback, reached only where the primal method returned no
    // vertex, so nothing in the search set exercises it. It still has to be right, and
    // it is the kind of thing that goes wrong quietly: on the branch this came from, a
    // basis that had drifted out of dual feasibility was reported as the optimum and
    // returned -1.0 for a relaxation whose value is 2087. Only a comparison against a
    // value the code under test did not produce catches that.
    //
    // Every sample is required to be answered or declined, and a declined one is
    // required to have done no work. Skipping whatever came back non-optimal is the
    // obvious way to write this and is a test that cannot fail: breaking the dual ratio
    // test deliberately stops it reaching the optimum on eighteen of the nineteen
    // samples, and a skipping version passes, slowly and silently.
    //
    // Objectives rather than solutions, because the two methods legitimately end on
    // different optimal vertices. That difference is the whole reason this is a
    // fallback rather than the default.
    let data = fixtures();
    let samples = data["samples"].as_array().unwrap();
    assert!(!samples.is_empty());

    for entry in samples {
        let file = entry["file"].as_str().unwrap();
        let problem = Problem::from_file(&samples_dir().join(file)).unwrap();
        let solution = Lp::relaxation(&problem).solve_cold_dual(200_000);
        if solution.status != LpStatus::Optimal {
            // The one legitimate outcome short of an answer: a model with a column
            // whose cost points at a bound it does not have cannot be parked dual
            // feasibly, and the entry declines it without spending anything.
            assert_eq!(
                solution.iterations, 0,
                "{file}: dual entry ran {} iterations and gave up at {:?}",
                solution.iterations, solution.status
            );
            continue;
        }
        assert_matches(
            problem.objective_value(solution.objective),
            entry["lp_relaxation"].as_f64().unwrap(),
            file,
        );
    }
}

#[test]
fn the_relaxation_of_a_fixed_problem_matches_its_reference_solution() {
    // Fixing every column to the reference integer solution should reproduce the
    // reference optimum exactly -- the check that bounds narrowed by branching are
    // honoured, on real instances rather than toy ones.
    let data = fixtures();
    for entry in data["samples"].as_array().unwrap() {
        let file = entry["file"].as_str().unwrap();
        let path = samples_dir().join(file);
        let mut problem = Problem::from_file(&path).unwrap();

        let solution: Vec<f64> = entry["solution"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_i64().unwrap() as f64)
            .collect();
        if solution.len() != problem.n_cols() {
            // Gurobi's column order need not match ours; skip rather than mis-pair.
            continue;
        }
        for (j, &v) in solution.iter().enumerate() {
            problem.col_lb[j] = v;
            problem.col_ub[j] = v;
        }

        let value = relaxation_value(&problem, file);
        assert_matches(value, entry["mip_optimum"].as_f64().unwrap(), file);
    }
    let _ = Path::new("");
}
