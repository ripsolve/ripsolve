//! Every bundled sample must read cleanly, and the two formats must agree.

use std::path::{Path, PathBuf};

use ripsolve::Problem;

fn samples_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../samples")
}

fn read(name: &str) -> Problem {
    let path = samples_dir().join(name);
    let p = Problem::from_file(&path).unwrap_or_else(|e| panic!("reading {name}: {e}"));
    p.validate()
        .unwrap_or_else(|e| panic!("validating {name}: {e}"));
    p
}

#[test]
fn every_sample_reads_and_validates() {
    let mut seen = 0;
    for entry in std::fs::read_dir(samples_dir()).expect("samples directory") {
        let path = entry.expect("dir entry").path();
        let is_model = path
            .extension()
            .and_then(|e| e.to_str())
            .is_some_and(|e| e.eq_ignore_ascii_case("lp") || e.eq_ignore_ascii_case("mps"));
        if !is_model {
            continue;
        }
        read(path.file_name().unwrap().to_str().unwrap());
        seen += 1;
    }
    assert!(seen >= 18, "expected the bundled samples, found {seen}");
}

#[test]
fn lp_and_mps_of_the_same_model_agree() {
    // v064c064 ships in both formats, which cross-checks the two parser paths
    // against each other without needing a reference solver.
    let lp = read("v064c064.lp");
    let mps = read("v064c064.mps");
    assert_eq!(lp.n_cols(), mps.n_cols());
    assert_eq!(lp.n_rows(), mps.n_rows());
    assert_eq!(lp.sense, mps.sense);
    assert_eq!(lp.obj, mps.obj);
    assert_eq!(lp.col_names, mps.col_names);
    assert_eq!(lp.matrix, mps.matrix);
    assert_eq!(lp.row_lb, mps.row_lb);
    assert_eq!(lp.row_ub, mps.row_ub);
}
