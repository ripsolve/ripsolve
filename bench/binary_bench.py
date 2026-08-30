#!/usr/bin/env python3
"""Benchmark every solver over the pure binary instances held locally.

A pure binary program is one whose every column is an integer confined to [0, 1].
That is a property of the model rather than of what any solver manages on it, so a
set chosen this way cannot flatter the solver being measured, unlike a set screened
by what other solvers already close.

Membership is decided by `cargo run --example classify`, which reads the models with
this solver's own reader, and not by MIPLIB's `binary` tag. The tag is looser than it
sounds: it means an instance has no *general* integers, and admits continuous columns
alongside the binaries. Of the 148 instances it tags binary and easy, eight are mixed
under the definition here, from neos-4382714-ruvuma with a single continuous column to
h80x6320 with 6320 of them against 6320 binaries.

Usage:  bench/binary_bench.py [seconds] [threads] [--refresh]

`--refresh` drops this solver's cached measurements before running, which is needed
whenever its binary has changed; the other solvers' entries are left alone because
theirs have not.
"""

import json
import pathlib
import subprocess
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import tractable as screen  # noqa: E402  - reuses its cached measurement layer

ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = ROOT / "bench" / "out"
SOLVERS = ["ripsolve", "HiGHS", "SCIP", "CBC", "commercial"]

argv = [a for a in sys.argv[1:] if a != "--refresh"]
REFRESH = "--refresh" in sys.argv
screen.LIMIT = float(argv[0]) if len(argv) > 0 else 60.0
screen.THREADS = int(argv[1]) if len(argv) > 1 else 16


def easy_names():
    """MIPLIB's easy classification, which every model here is also required to be in."""
    path = OUT / "easy-v18.test"
    return {l.strip().removesuffix(".mps.gz") for l in path.read_text().splitlines() if l.strip()}


def pure_binary():
    """Names of the locally held easy models whose every column is binary."""
    easy = easy_names()
    models = sorted(p for p in screen.CACHE.glob("*.mps") if p.stem in easy)
    classify = ROOT / "target" / "release" / "examples" / "classify"
    if not classify.exists():
        sys.exit("build the classifier first: cargo build --release --example classify")
    out = subprocess.run([str(classify), *[str(p) for p in models]],
                         capture_output=True, text=True, timeout=1800).stdout
    chosen = []
    for line in out.splitlines():
        parts = line.split("\t")
        if len(parts) < 7 or parts[1] == "unreadable":
            continue
        name, rows, cols, nnz, binary, general, continuous = parts
        if int(general) == 0 and int(continuous) == 0:
            chosen.append((name, int(rows), int(cols), int(nnz)))
    return sorted(chosen, key=lambda t: t[3])


def main():
    models = pure_binary()
    cache = screen.load(screen.SCREEN)
    if REFRESH:
        stale = [k for k in cache if k.split("|")[1] == "ripsolve"]
        for key in stale:
            del cache[key]
        print(f"dropped {len(stale)} cached ripsolve measurements\n")

    print(f"{len(models)} pure binary models, {screen.LIMIT:g}s, "
          f"{screen.THREADS} threads\n")
    rows = []
    for name, r, c, nz in models:
        result = {k: screen.measure(k, name, screen.CACHE / f"{name}.mps", cache)
                  for k in SOLVERS}
        rows.append((name, r, c, nz, result))
        print("%-24s %s" % (name, "  ".join(
            "%s=%s/%ss" % (k, result[k]["status"][:3], result[k]["seconds"])
            for k in SOLVERS)), flush=True)

    tally = {k: sum(1 for row in rows if screen.closed(row[4][k])) for k in SOLVERS}
    print("\nclosed (optimal or proved infeasible): "
          + ", ".join(f"{k} {tally[k]}/{len(rows)}" for k in SOLVERS))
    (OUT / "binary_results.json").write_text(json.dumps(
        {"rows": [[n, r, c, nz, res] for n, r, c, nz, res in rows], "tally": tally,
         "limit": screen.LIMIT, "threads": screen.THREADS}, indent=1))
    print(f"wrote {OUT / 'binary_results.json'}")


main()
