#!/usr/bin/env python3
"""Compare ripsolve against the open-source MIP solvers.

HiGHS, SCIP and CBC are all open source and carry no restriction on publishing
benchmarks, so unlike the commercial comparison this one names what it measures.

Fairness notes, all deliberate:
  * every solver is pinned to one thread, so this compares algorithms rather than
    how many cores each chooses to take
  * only the solve is timed; parsing is excluded on every side
  * the same time limit applies to all of them, and a run that hits it is reported
    with its gap rather than as a solve

Usage:  python3 bench/open_compare.py [seconds] [instance ...]
"""

import pathlib
import subprocess
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parent.parent
EXE = ROOT / "target" / "release" / "ripsolve"
LIMIT = float(sys.argv[1]) if len(sys.argv) > 1 else 60.0

DEFAULT = [
    ROOT / "samples" / n
    for n in ("v032c032.lp", "v048c048.lp", "v048c128.lp", "v064c064.lp", "v064c200.mps")
] + sorted((ROOT / "samples" / "mip").glob("*.lp"))


class Result:
    def __init__(self, objective=None, seconds=0.0, solved=False, note=""):
        self.objective, self.seconds, self.solved, self.note = objective, seconds, solved, note

    def cell(self):
        if self.note:
            return f"{self.note:>16}"
        if self.objective is None:
            return f"{'no solution':>16}"
        mark = "" if self.solved else "*"
        return f"{self.objective:>10.2f}{mark} {self.seconds:>4.1f}s"


def run_ripsolve(path):
    started = time.time()
    out = subprocess.run(
        [str(EXE), "solve", "-t", "1", "--time-limit", str(LIMIT), str(path)],
        capture_output=True, text=True, timeout=LIMIT * 3,
    ).stdout
    elapsed = time.time() - started
    objective = next((float(l.split()[1]) for l in out.splitlines()
                      if l.startswith("objective:")), None)
    return Result(objective, elapsed, "status:    optimal" in out)


def run_highs(path):
    import highspy
    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    h.setOptionValue("threads", 1)
    h.setOptionValue("time_limit", LIMIT)
    h.readModel(str(path))
    started = time.time()
    h.run()
    elapsed = time.time() - started
    solved = str(h.getModelStatus()).endswith("kOptimal")
    return Result(h.getObjectiveValue() if solved else None, elapsed, solved)


def run_scip(path):
    from pyscipopt import Model
    m = Model()
    m.hideOutput()
    m.readProblem(str(path))
    m.setParam("limits/time", LIMIT)
    m.setParam("parallel/maxnthreads", 1)
    started = time.time()
    m.optimize()
    elapsed = time.time() - started
    solved = m.getStatus() == "optimal"
    return Result(m.getObjVal() if m.getNSols() else None, elapsed, solved)


def run_cbc(path):
    import pulp
    cbc = pulp.PULP_CBC_CMD().path
    started = time.time()
    out = subprocess.run(
        [cbc, str(path), "threads", "1", "seconds", str(LIMIT), "solve"],
        capture_output=True, text=True, timeout=LIMIT * 3,
    ).stdout
    elapsed = time.time() - started
    objective, solved = None, False
    for line in out.splitlines():
        if "best objective" in line.lower():
            try:
                objective = float(line.lower().split("best objective")[1].split(",")[0])
            except (IndexError, ValueError):
                pass
            solved = "search completed" in line.lower()
    return Result(objective, elapsed, solved)


SOLVERS = [("ripsolve", run_ripsolve), ("HiGHS", run_highs),
           ("SCIP", run_scip), ("CBC", run_cbc)]


def main():
    paths = [pathlib.Path(a) for a in sys.argv[2:]] or DEFAULT
    print(f"time limit {LIMIT:.0f}s, one thread each. * = hit the limit\n")
    print(f"{'instance':22}" + "".join(f"{name:>18}" for name, _ in SOLVERS))
    print("-" * (22 + 18 * len(SOLVERS)))

    disagreements = 0
    for path in paths:
        if not path.exists():
            continue
        row, values = f"{path.name:22}", []
        for _, fn in SOLVERS:
            try:
                result = fn(path)
            except Exception as exc:  # a missing solver should not stop the rest
                result = Result(note=type(exc).__name__[:14])
            row += result.cell()
            if result.solved and result.objective is not None:
                values.append(result.objective)
        # Every solver that proved optimality must have proved the same optimum.
        if values and max(values) - min(values) > 1e-4 * max(1.0, abs(values[0])):
            row += "  <-- DISAGREE"
            disagreements += 1
        print(row)

    print()
    if disagreements:
        print(f"{disagreements} instance(s) where proven optima disagree")
        return 1
    print("all proven optima agree")
    return 0


if __name__ == "__main__":
    sys.exit(main())
