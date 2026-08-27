#!/usr/bin/env python3
"""Run every solver over the tractable set and write the comparison as a web page.

The set comes from `bench/tractable.py`: MIPLIB instances at least two of HiGHS, SCIP
and CBC close within the budget. On those, a failure is this solver's own and is known
to be reachable, which is what makes the comparison worth drawing.

Usage:  bench/tractable_chart.py [seconds] [threads] [out.html]
"""

import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import tractable as screen  # noqa: E402  - reuses its cached measurement layer

# Set explicitly: the budget is part of every cache key, so measuring at one budget and
# reading at another silently answers a different question.
screen.LIMIT = float(sys.argv[1]) if len(sys.argv) > 1 else 60.0
screen.THREADS = int(sys.argv[2]) if len(sys.argv) > 2 else 16

ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = ROOT / "bench" / "out"
SOLVERS = ["ripsolve", "HiGHS", "SCIP", "CBC", "commercial"]
DEST = pathlib.Path(sys.argv[3]) if len(sys.argv) > 3 else OUT / "tractable.html"


def main():
    names = json.loads((OUT / "tractable.json").read_text())
    cache = screen.load(screen.SCREEN)
    rows = []
    for name in names:
        path = screen.CACHE / f"{name}.mps"
        if not path.exists():
            continue
        result = {k: screen.measure(k, name, path, cache) for k in SOLVERS}
        rows.append((name, result))
        summary = "  ".join(
            f"{k}={result[k]['status'][:3]}/{result[k]['seconds']}s" for k in SOLVERS)
        print(f"{name:30} {summary}", flush=True)

    tally = {k: sum(1 for _, r in rows if r[k]["status"] == "optimal") for k in SOLVERS}
    print("\nsolved:", ", ".join(f"{k} {tally[k]}/{len(rows)}" for k in SOLVERS))
    (OUT / "tractable_results.json").write_text(
        json.dumps({"rows": [[n, r] for n, r in rows], "tally": tally,
                    "limit": screen.LIMIT, "threads": screen.THREADS}, indent=1))
    print(f"wrote {OUT / 'tractable_results.json'}")


main()
