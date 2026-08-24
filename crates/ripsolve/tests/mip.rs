//! Mixed-integer models: general integers and continuous columns together.
//!
//! Everything else in this suite is binary. These instances mix binary,
//! bounded-integer and continuous columns in the same rows, which is where a
//! surviving binary-only assumption shows up — a reduction that is sound for 0/1
//! columns and wrong for the others does not fail loudly, it proves a worse
//! solution optimal.

use std::path::{Path, PathBuf};

use ripsolve::Problem;
use ripsolve::model::VarType;
use ripsolve::search::{self, Options, Status};

fn samples_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../samples")
}

fn fixtures() -> serde_json::Value {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/mip.json");
    serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap()
}

/// The reported point must be feasible, integral where required, and score what
/// was claimed.
fn assert_valid(problem: &Problem, x: &[f64], objective: f64, name: &str) {
    assert_eq!(x.len(), problem.n_cols(), "{name}: solution length");
    for (j, &v) in x.iter().enumerate() {
        assert!(
            v >= problem.col_lb[j] - 1e-6 && v <= problem.col_ub[j] + 1e-6,
            "{name}: column {j} = {v} outside [{}, {}]",
            problem.col_lb[j],
            problem.col_ub[j]
        );
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
        let activity: f64 = cols.iter().zip(vals).map(|(&j, &a)| a * x[j]).sum();
        assert!(
            activity >= problem.row_lb[i] - 1e-6 && activity <= problem.row_ub[i] + 1e-6,
            "{name}: row {i} activity {activity} outside [{}, {}]",
            problem.row_lb[i],
            problem.row_ub[i]
        );
    }

    let internal: f64 = problem.obj.iter().zip(x).map(|(c, v)| c * v).sum();
    assert!(
        (problem.objective_value(internal) - objective).abs() < 1e-6,
        "{name}: reported {objective} but the point scores {}",
        problem.objective_value(internal)
    );
}

#[test]
fn mixed_integer_optima_match_the_reference_oracle() {
    let data = fixtures();
    let instances = data["instances"].as_array().unwrap();
    assert!(!instances.is_empty());

    for entry in instances {
        let file = entry["file"].as_str().unwrap();
        let expected = entry["mip_optimum"].as_f64().unwrap();
        let problem = Problem::from_file(&samples_dir().join(file)).unwrap();

        let solution = search::solve(
            &problem,
            Options {
                threads: 1,
                ..Options::default()
            },
        );
        assert_eq!(solution.status, Status::Optimal, "{file}");
        let got = solution.objective.unwrap();
        assert!(
            (got - expected).abs() <= 1e-6 * expected.abs().max(1.0),
            "{file}: got {got}, expected {expected}"
        );
        assert_valid(&problem, &solution.x, got, file);
    }
}

#[test]
fn column_types_survive_a_bounds_section() {
    // `lp_parser_rs` lets a Bounds entry overwrite a variable's declared type, so a
    // General integer that is also bounded comes back continuous. The reader
    // recovers integrality from the source sections; without that, these models
    // would silently solve as pure LPs.
    let data = fixtures();
    let mut saw_continuous = false;
    for entry in data["instances"].as_array().unwrap() {
        let file = entry["file"].as_str().unwrap();
        let problem = Problem::from_file(&samples_dir().join(file)).unwrap();
        let expected: Vec<&str> = entry["integer_columns"]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| v.as_str().unwrap())
            .collect();

        let got: Vec<&str> = problem
            .col_names
            .iter()
            .enumerate()
            .filter(|&(j, _)| problem.is_integer(j))
            .map(|(_, n)| n.as_str())
            .collect();
        assert_eq!(got, expected, "{file}: integer columns disagree");
        saw_continuous |= problem.col_type.contains(&VarType::Continuous);
    }
    // Not every generated instance happens to have one, but the set must, or this
    // test is not exercising the mixed case at all.
    assert!(
        saw_continuous,
        "no instance in the fixture set has a continuous column"
    );
}

#[test]
fn a_relaxation_integral_in_its_integer_columns_is_accepted() {
    // The distinguishing MIP case: the LP optimum is fractional, but only in
    // continuous columns, so it is already a valid answer and needs no branching.
    let text = "\
Minimize
 obj: 2 a + 3 b
Subject To
 c0: a + b >= 1.5
Bounds
 0 <= a <= 1
 0 <= b <= 5
General
 a
End
";
    let path = std::env::temp_dir().join("ripsolve_mip_relaxation.lp");
    std::fs::write(&path, text).unwrap();
    let problem = Problem::from_file(&path).unwrap();
    assert!(problem.is_integer(0) && !problem.is_integer(1));

    let solution = search::solve(
        &problem,
        Options {
            threads: 1,
            ..Options::default()
        },
    );
    assert_eq!(solution.status, Status::Optimal);
    // a = 1 costs 2 and leaves b = 0.5 costing 1.5; total 3.5.
    assert!(
        (solution.objective.unwrap() - 3.5).abs() < 1e-6,
        "{:?}",
        solution.objective
    );
    assert_valid(
        &problem,
        &solution.x,
        solution.objective.unwrap(),
        "relaxation",
    );
}

#[test]
fn a_general_integer_branches_over_its_whole_range() {
    // Branching splits a range rather than fixing to 0/1, so an optimum at 4 must
    // be reachable.
    let text = "\
Minimize
 obj: - x
Subject To
 c0: 3 x <= 13
Bounds
 0 <= x <= 9
General
 x
End
";
    let path = std::env::temp_dir().join("ripsolve_mip_general.lp");
    std::fs::write(&path, text).unwrap();
    let problem = Problem::from_file(&path).unwrap();

    let solution = search::solve(
        &problem,
        Options {
            threads: 1,
            ..Options::default()
        },
    );
    assert_eq!(solution.status, Status::Optimal);
    assert!((solution.x[0] - 4.0).abs() < 1e-6, "x = {}", solution.x[0]);
    assert!((solution.objective.unwrap() + 4.0).abs() < 1e-6);
}
