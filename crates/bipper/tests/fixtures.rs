//! The reference fixtures must describe the instances being solved *now*.
//!
//! `bench/refresh_fixtures.py` records an external solver's LP-relaxation and
//! MIP-optimal values for each generated instance and bundled sample. Those numbers
//! are only meaningful for the exact instance they came from, so every entry
//! carries a digest of its LP text. If a generated instance or a sample file
//! changes, these tests fail and say to refresh — which is the intended outcome,
//! and much better than quietly checking new instances against stale values.

use std::path::{Path, PathBuf};

use bipper::generate::{Kind, Spec, lp_digest};

pub fn fixture_path() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/reference.json")
}

pub fn samples_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../samples")
}

pub fn fixtures() -> serde_json::Value {
    let path = fixture_path();
    let text = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("reading {}: {e}", path.display()));
    serde_json::from_str(&text).expect("fixture file is valid JSON")
}

pub fn kind_from_str(s: &str) -> Kind {
    match s {
        "knapsack" => Kind::Knapsack,
        "covering" => Kind::Covering,
        "signed" => Kind::Signed,
        other => panic!("unknown kind {other:?} in fixture file"),
    }
}

pub fn spec_of(entry: &serde_json::Value) -> Spec {
    Spec {
        kind: kind_from_str(entry["kind"].as_str().unwrap()),
        n_cols: entry["n_cols"].as_u64().unwrap() as usize,
        n_rows: entry["n_rows"].as_u64().unwrap() as usize,
        seed: entry["seed"].as_u64().unwrap(),
    }
}

/// An entry's display name: generated instances carry `name`, samples carry `file`.
fn label(entry: &serde_json::Value) -> &str {
    entry["name"]
        .as_str()
        .or_else(|| entry["file"].as_str())
        .expect("entry has a name or a file")
}

#[test]
fn generated_fixtures_match_the_current_generator() {
    let data = fixtures();
    let instances = data["instances"].as_array().expect("instances array");
    assert!(!instances.is_empty(), "fixture file has no instances");

    for entry in instances {
        let spec = spec_of(entry);
        let name = label(entry);
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
fn sample_fixtures_match_the_files_on_disk() {
    let data = fixtures();
    let samples = data["samples"].as_array().expect("samples array");
    assert!(!samples.is_empty(), "fixture file has no samples");

    for entry in samples {
        let file = entry["file"].as_str().expect("sample entry has a file");
        let path = samples_dir().join(file);
        let text = std::fs::read_to_string(&path)
            .unwrap_or_else(|e| panic!("reading {}: {e}", path.display()));
        let expected = entry["digest"].as_str().expect("entry has a digest");
        assert_eq!(
            format!("{:016x}", lp_digest(&text)),
            expected,
            "{file} has changed since its reference values were recorded; \
             re-run `python3 bench/refresh_fixtures.py`"
        );
    }
}

#[test]
fn every_sample_file_has_a_fixture_entry() {
    // Adding a sample without refreshing would otherwise leave it silently untested.
    let data = fixtures();
    let recorded: Vec<&str> = data["samples"]
        .as_array()
        .unwrap()
        .iter()
        .map(|e| e["file"].as_str().unwrap())
        .collect();

    for dir_entry in std::fs::read_dir(samples_dir()).expect("samples directory") {
        let path = dir_entry.expect("dir entry").path();
        let is_model = path
            .extension()
            .and_then(|e| e.to_str())
            .is_some_and(|e| e.eq_ignore_ascii_case("lp") || e.eq_ignore_ascii_case("mps"));
        if !is_model {
            continue;
        }
        let name = path.file_name().unwrap().to_str().unwrap();
        assert!(
            recorded.contains(&name),
            "{name} has no fixture entry; re-run `python3 bench/refresh_fixtures.py`"
        );
    }
}

#[test]
fn fixture_values_are_self_consistent() {
    let data = fixtures();
    let all = data["instances"]
        .as_array()
        .unwrap()
        .iter()
        .chain(data["samples"].as_array().unwrap());

    for entry in all {
        let name = label(entry);
        let lp = entry["lp_relaxation"].as_f64().unwrap();
        let mip = entry["mip_optimum"].as_f64().unwrap();
        // Relaxing integrality can only help, so the relaxation bounds the optimum:
        // from below when minimizing, from above when maximizing.
        match entry["sense"].as_str().unwrap_or("minimize") {
            "minimize" => {
                assert!(
                    lp <= mip + 1e-6,
                    "{name}: relaxation {lp} exceeds optimum {mip}"
                )
            }
            "maximize" => {
                assert!(
                    lp >= mip - 1e-6,
                    "{name}: relaxation {lp} below optimum {mip}"
                )
            }
            other => panic!("{name}: unknown sense {other:?}"),
        }
        assert!(
            entry["solution"]
                .as_array()
                .unwrap()
                .iter()
                .all(|v| matches!(v.as_i64(), Some(0) | Some(1))),
            "{name}: reference solution is not binary"
        );
    }
}
