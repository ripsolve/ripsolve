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
fn a_search_that_skipped_a_node_only_claims_optimality_when_it_holds() {
    // A node whose LP runs out of iterations leaves its subtree unexamined, which
    // normally forfeits the optimality claim. It does not forfeit it when the incumbent
    // has since overtaken the bound that node inherited, because that bound holds over
    // the whole subtree and the search prunes on exactly that test everywhere else.
    //
    // Squeezing the per-node iteration limit is what makes nodes run out. Whatever the
    // search then reports, `Optimal` has to mean the reference optimum: relaxing the
    // rule to say "optimal" more often is only safe if it is never saying it wrongly.
    let data = fixtures();
    let mut claimed = 0;
    for entry in data["samples"].as_array().unwrap() {
        let file = entry["file"].as_str().unwrap();
        let problem = Problem::from_file(&samples_dir().join(file)).unwrap();
        let expected = entry["mip_optimum"].as_f64().unwrap();
        for limit in [1usize, 4, 16, 64] {
            let solution = search::solve(
                &problem,
                Options {
                    max_iterations_per_node: limit,
                    ..Options::default()
                },
            );
            if solution.status != Status::Optimal {
                continue;
            }
            claimed += 1;
            let got = solution.objective.unwrap_or(f64::NAN);
            let scale = expected.abs().max(1.0);
            assert!(
                (got - expected).abs() <= 1e-6 * scale,
                "{file} at {limit} iterations per node: claimed optimal {got},                  reference {expected}"
            );
        }
    }
    assert!(
        claimed > 0,
        "no run claimed optimality, so the claim was never checked"
    );
}

#[test]
fn reduced_cost_fixing_does_not_change_any_optimum() {
    // Reduced cost fixing narrows columns on the strength of a proof about the whole
    // tree: that no solution better than the incumbent has them anywhere else. A proof
    // that is a little too strong removes the optimum and the search reports the
    // second best answer with every appearance of certainty, which is why this is
    // checked against the same search with the claim turned off rather than against a
    // recorded value that would also have to be trusted.
    let data = fixtures();
    let samples = data["samples"].as_array().unwrap().iter().map(|entry| {
        let file = entry["file"].as_str().unwrap().to_string();
        (
            Problem::from_file(&samples_dir().join(&file)).unwrap(),
            file,
        )
    });
    // The generated families as well as the samples, and it is worth knowing how much
    // this catches. Setting the travel cap to zero, so that every nonbasic column is
    // pinned where it sits, is caught on `v064c200`, which answers 1039 against 225.
    // Merely halving the cap is not caught by either set, so this is a guard against a
    // proof that is wrong in kind rather than a proof that is wrong at the margin.
    let generated = data["instances"].as_array().unwrap().iter().map(|entry| {
        let spec = spec_of(entry);
        let name = entry["name"].as_str().unwrap().to_string();
        (
            Problem::from_lp(&LpProblem::parse(&spec.to_lp()).unwrap()).unwrap(),
            name,
        )
    });
    for (problem, file) in samples.chain(generated) {
        let with = search::solve(
            &problem,
            Options {
                fix_by_reduced_cost: true,
                ..Options::default()
            },
        );
        let without = search::solve(
            &problem,
            Options {
                fix_by_reduced_cost: false,
                ..Options::default()
            },
        );

        assert_eq!(with.status, without.status, "{file}");
        match (with.objective, without.objective) {
            (Some(a), Some(b)) => {
                assert!((a - b).abs() < 1e-6, "{file}: fixing {a}, plain {b}")
            }
            (None, None) => {}
            (a, b) => panic!("{file}: fixing {a:?}, plain {b:?}"),
        }
    }
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

        // The assignment, not just its value. An objective is computed on whatever
        // model the search actually ran, so it stays self-consistent even when the
        // point handed back does not belong to the original: a reduction that
        // renumbered columns, or moved a constant it should have kept, would report a
        // plausible number attached to the wrong vector. Comparing objectives alone
        // cannot see that; asking the original model what it makes of the point can.
        if let Some(reported) = with.objective {
            assert_eq!(
                with.x.len(),
                problem.n_cols(),
                "{file}: presolve returned {} values for {} columns",
                with.x.len(),
                problem.n_cols()
            );
            assert!(
                ripsolve::heuristic::is_feasible(&problem, &with.x, 1e-6),
                "{file}: the point presolve returned is not feasible for the original"
            );
            let restated = problem.objective_value(
                problem.obj.iter().zip(&with.x).map(|(c, v)| c * v).sum::<f64>(),
            );
            assert!(
                (restated - reported).abs() < 1e-6,
                "{file}: reported {reported} but the point is worth {restated}"
            );
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
