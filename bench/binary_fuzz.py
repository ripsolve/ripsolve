#!/usr/bin/env python3
"""Differential test of ripsolve against HiGHS on pure binary models.

`mip_fuzz.py` mixes binary, integer and continuous columns, which is right for most of
the solver and wrong for the parts that only act on a pure binary model: the
conflict-driven search before the relaxation reports a point, and proves one optimal,
only when every column is binary. This generates models of that kind with the solver's
own generator, at sizes where that search does the whole job on some and none of it on
others, and requires the objectives to agree.

Usage:  bench/binary_fuzz.py [trials] [first seed]
"""

import pathlib
import subprocess
import sys

import highspy

ROOT = pathlib.Path(__file__).resolve().parent.parent
EXE = ROOT / "target" / "release" / "ripsolve"
OUT = ROOT / "bench" / "out" / "binaryfuzz"
OUT.mkdir(parents=True, exist_ok=True)


def ours(path):
    out = subprocess.run([str(EXE), "solve", "-t", "2", "--gap", "0", str(path)],
                         capture_output=True, text=True, timeout=600).stdout
    if "Infeasible" in out:
        return "infeasible"
    for line in out.splitlines():
        if line.startswith("objective:"):
            assert "status:    optimal" in out, out
            return float(line.split()[1])
    raise RuntimeError(out)


def reference(path):
    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    h.setOptionValue("mip_rel_gap", 0.0)
    h.readModel(str(path))
    h.run()
    if h.getModelStatus() == highspy.HighsModelStatus.kInfeasible:
        return "infeasible"
    return h.getInfo().objective_function_value


def main():
    trials = int(sys.argv[1]) if len(sys.argv) > 1 else 90
    kinds = ["knapsack", "covering", "signed"]
    failures = 0
    first = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    for t in range(first, first + trials):
        kind = kinds[t % 3]
        cols, rows = [(20, 10), (40, 25), (60, 30)][(t // 3) % 3]
        path = OUT / f"{kind}_{cols}_{rows}_{t}.lp"
        subprocess.run([str(EXE), "gen", "--kind", kind, "--cols", str(cols), "--rows", str(rows),
                        "--seed", str(t), "-o", str(path)], check=True)
        got, want = ours(path), reference(path)
        same = got == want if isinstance(want, str) or isinstance(got, str) \
            else abs(got - want) <= 1e-6 * max(1.0, abs(want))
        if not same:
            failures += 1
            print(f"MISMATCH {path.name}: ripsolve {got}, HiGHS {want}")
    print(f"{trials - failures}/{trials} agree")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
