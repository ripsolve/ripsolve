#!/usr/bin/env python3
"""Regenerate the reference-value fixtures used by ripsolve's tests.

Gurobi is a *test oracle* here, not a dependency of the solver: this script runs
by hand, writes the expected LP-relaxation and MIP-optimal values into a JSON
fixture, and `cargo test` afterwards reads only that file. Nothing in ripsolve
links against Gurobi, and the test suite runs on machines without it.

Usage:  python3 bench/refresh_fixtures.py [--out PATH]
"""

import argparse
import json
import pathlib
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
DEFAULT_OUT = ROOT / "crates" / "ripsolve" / "tests" / "fixtures" / "reference.json"

# Kept small on purpose: these are correctness fixtures, not a performance suite,
# so every instance must be solvable by Gurobi in well under a second.
SPECS = [
    ("knapsack", 20, 10, 1), ("knapsack", 30, 15, 2), ("knapsack", 45, 20, 3),
    ("covering", 25, 30, 1), ("covering", 40, 50, 2), ("covering", 60, 80, 3),
    ("signed",   20, 20, 1), ("signed",   32, 32, 2), ("signed",   48, 48, 3),
]


def lp_digest(text):
    """FNV-1a (64-bit), matching `ripsolve::generate::lp_digest`.

    The Rust test regenerates each instance and compares this digest, so a change
    to the generator is caught as a stale fixture rather than silently pairing new
    instances with old reference values.
    """
    h = 0xCBF29CE484222325
    for b in text.encode():
        h = ((h ^ b) * 0x100000001B3) & 0xFFFFFFFFFFFFFFFF
    return h


def ripsolve_gen(kind, cols, rows, seed, path):
    exe = ROOT / "target" / "release" / "ripsolve"
    if not exe.exists():
        sys.exit(f"{exe} not found; run `cargo build --release` first")
    subprocess.run(
        [str(exe), "gen", "--kind", kind, "--cols", str(cols),
         "--rows", str(rows), "--seed", str(seed), "-o", str(path)],
        check=True,
    )


def reference_values(lp_path):
    import gurobipy as gp

    env = gp.Env(params={"OutputFlag": 0})
    model = gp.read(str(lp_path), env=env)
    # For a maximization the relaxation is an upper bound, not a lower one, so the
    # sense has to travel with the values for them to be checkable.
    sense = "maximize" if model.ModelSense == gp.GRB.MAXIMIZE else "minimize"

    # Root LP relaxation, with presolve and cuts off so the value is the honest
    # relaxation of the model as written -- that is what ripsolve's simplex must
    # reproduce, before any of its own presolve or cuts exist.
    relaxed = model.relax()
    relaxed.setParam("Presolve", 0)
    relaxed.setParam("Method", 1)  # dual simplex
    relaxed.optimize()
    if relaxed.Status != gp.GRB.OPTIMAL:
        sys.exit(f"{lp_path.name}: relaxation status {relaxed.Status}")
    lp_value = relaxed.ObjVal

    model.setParam("MIPGap", 0.0)
    model.optimize()
    if model.Status != gp.GRB.OPTIMAL:
        sys.exit(f"{lp_path.name}: MIP status {model.Status}")

    return lp_value, model.ObjVal, [v.X for v in model.getVars()], sense


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=pathlib.Path, default=DEFAULT_OUT)
    args = ap.parse_args()

    scratch = ROOT / "bench" / "out"
    scratch.mkdir(parents=True, exist_ok=True)

    entries = []
    for kind, cols, rows, seed in SPECS:
        name = f"{kind}_c{cols}_r{rows}_s{seed}"
        lp_path = scratch / f"{name}.lp"
        ripsolve_gen(kind, cols, rows, seed, lp_path)
        lp_value, mip_value, solution, sense = reference_values(lp_path)
        entries.append({
            "kind": kind, "n_cols": cols, "n_rows": rows, "seed": seed,
            "name": name,
            "sense": sense,
            "digest": f"{lp_digest(lp_path.read_text()):016x}",
            "lp_relaxation": lp_value,
            "mip_optimum": mip_value,
            # Rounded: Gurobi returns binaries as 0.9999999997 and friends.
            "solution": [int(round(x)) for x in solution],
        })
        print(f"{name:28s} lp={lp_value:12.6f}  mip={mip_value:12.6f}")

    # The bundled samples cost nothing extra to cover and exercise the readers on
    # hand-written files rather than only on generator output.
    samples = []
    for path in sorted((ROOT / "samples").iterdir()):
        if path.suffix.lower() not in (".lp", ".mps"):
            continue
        lp_value, mip_value, solution, sense = reference_values(path)
        samples.append({
            "file": path.name,
            "sense": sense,
            "digest": f"{lp_digest(path.read_text()):016x}",
            "lp_relaxation": lp_value,
            "mip_optimum": mip_value,
            "solution": [int(round(x)) for x in solution],
        })
        print(f"{path.name:28s} lp={lp_value:12.6f}  mip={mip_value:12.6f}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(
        {"generator": "ripsolve gen", "oracle": "gurobi",
         "instances": entries, "samples": samples},
        indent=2,
    ) + "\n")
    print(f"\nwrote {len(entries)} instances and {len(samples)} samples "
          f"to {args.out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
