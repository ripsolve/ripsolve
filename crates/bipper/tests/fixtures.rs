//! The reference fixtures must describe the instances the generator produces now.
//!
//! `bench/refresh_fixtures.py` records an external solver's LP-relaxation and
//! MIP-optimal values for each generated instance. Those numbers are only
//! meaningful for the exact instance they were computed from, so every entry
//! carries a digest of its LP text. If the generator changes, this test fails and
//! says to refresh — which is the intended outcome, and much better than quietly
//! checking a new instance against stale reference values.

use std::path::Path;

use bipper::generate::{lp_digest, Kind, Spec};

fn fixtures() -> serde_json::Value {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/reference.json");
    let text = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("reading {}: {e}", path.display()));
    serde_json::from_str(&text).expect("fixture file is valid JSON")
}

fn kind_from_str(s: &str) -> Kind {
    match s {
        "knapsack" => Kind::Knapsack,
        "covering" => Kind::Covering,
        "signed" => Kind::Signed,
        other => panic!("unknown kind {other:?} in fixture file"),
    }
}

fn spec_of(entry: &serde_json::Value) -> Spec {
    Spec {
        kind: kind_from_str(entry["kind"].as_str().unwrap()),
        n_cols: entry["n_cols"].as_u64().unwrap() as usize,
        n_rows: entry["n_rows"].as_u64().unwrap() as usize,
        seed: entry["seed"].as_u64().unwrap(),
    }
}

#[test]
fn fixtures_match_the_current_generator() {
    let data = fixtures();
    let instances = data["instances"].as_array().expect("instances array");
    assert!(!instances.is_empty(), "fixture file has no instances");

    for entry in instances {
        let spec = spec_of(entry);
        let name = entry["name"].as_str().unwrap();
        assert_eq!(spec.name(), name, "fixture name disagrees with the spec");

        let expected = entry["digest"].as_str().expect("entry has a digest");
        let actual = format!("{:016x}", lp_digest(&spec.to_lp()));
        assert_eq!(
            actual, expected,
            "{name} has changed since its reference values were recorded; \
             re-run `python3 bench/refresh_fixtures.py`"
        );
    }
}

#[test]
fn fixture_values_are_self_consistent() {
    for entry in fixtures()["instances"].as_array().unwrap() {
        let name = entry["name"].as_str().unwrap();
        let lp = entry["lp_relaxation"].as_f64().unwrap();
        let mip = entry["mip_optimum"].as_f64().unwrap();
        // The relaxation drops the integrality restriction, so for a minimization it
        // can never be worse than the integer optimum.
        assert!(lp <= mip + 1e-6, "{name}: relaxation {lp} exceeds optimum {mip}");

        let solution = entry["solution"].as_array().unwrap();
        assert_eq!(
            solution.len(),
            entry["n_cols"].as_u64().unwrap() as usize,
            "{name}: solution length disagrees with n_cols"
        );
        assert!(
            solution.iter().all(|v| matches!(v.as_i64(), Some(0) | Some(1))),
            "{name}: reference solution is not binary"
        );
    }
}
