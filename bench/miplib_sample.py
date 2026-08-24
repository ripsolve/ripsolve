#!/usr/bin/env python3
"""Run ripsolve over a seeded random sample of the MIPLIB 2017 benchmark set.

The point of sampling rather than choosing is that a hand-picked set says nothing.
An earlier pass here used ten instances recalled from memory; every one of them was
smaller than 83% of the benchmark set, so the results could not support any claim
about MIPLIB generally.

Every instance in the sample is reported — including ones that fail to download,
fail to parse, or time out. Dropping those is how a benchmark flatters itself.

HiGHS runs alongside as an open-source reference, both to check answers and to say
whether a timeout is ripsolve's problem or the instance's.

Usage:  bench/miplib_sample.py [count] [seconds] [seed]
"""

import csv
import gzip
import pathlib
import random
import re
import shutil
import subprocess
import sys
import time
import urllib.request

ROOT = pathlib.Path(__file__).resolve().parent.parent
EXE = ROOT / "target" / "release" / "ripsolve"
OUT = ROOT / "bench" / "out"
CACHE = OUT / "miplib"
LIST_URL = "https://miplib.zib.de/downloads/benchmark-v2.test"
DATA_URL = "https://miplib.zib.de/WebData/instances/{}.mps.gz"

COUNT = int(sys.argv[1]) if len(sys.argv) > 1 else 25
LIMIT = float(sys.argv[2]) if len(sys.argv) > 2 else 60.0
SEED = int(sys.argv[3]) if len(sys.argv) > 3 else 20260824
# Instances larger than this are recorded as skipped rather than downloaded; the
# benchmark set reaches 40 MB compressed and the point is a survey, not a stress
# test of the network.
MAX_DOWNLOAD = 25 * 1024 * 1024


def instance_names():
    path = OUT / "benchmark-v2.test"
    if not path.exists():
        OUT.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(LIST_URL, path)
    return [l.strip().removesuffix(".mps.gz") for l in path.read_text().splitlines() if l.strip()]


def fetch(name):
    """Return the local .mps path, or a reason it is unavailable."""
    CACHE.mkdir(parents=True, exist_ok=True)
    mps = CACHE / f"{name}.mps"
    if mps.exists():
        return mps, ""
    gz = CACHE / f"{name}.mps.gz"
    try:
        with urllib.request.urlopen(DATA_URL.format(name), timeout=120) as response:
            size = int(response.headers.get("Content-Length", 0))
            if size > MAX_DOWNLOAD:
                return None, f"skipped ({size // 1024 // 1024} MB)"
            gz.write_bytes(response.read())
        with gzip.open(gz, "rb") as src, open(mps, "wb") as dst:
            shutil.copyfileobj(src, dst)
        gz.unlink(missing_ok=True)
        return mps, ""
    except Exception as exc:
        return None, f"download failed: {type(exc).__name__}"


def shape(path):
    out = subprocess.run([str(EXE), "info", str(path)], capture_output=True, text=True)
    m = re.search(r"(\d+) columns, (\d+) rows, (\d+) nonzeros", out.stdout)
    if not m:
        reason = out.stderr.strip().splitlines()[-1] if out.stderr.strip() else "unreadable"
        return None, reason.removeprefix("Caused by:").strip()
    return tuple(int(g) for g in m.groups()), ""


def run_ripsolve(path):
    started = time.time()
    try:
        out = subprocess.run(
            [str(EXE), "solve", "-t", "1", "--time-limit", str(LIMIT), str(path)],
            capture_output=True, text=True, timeout=LIMIT * 4,
        ).stdout
    except subprocess.TimeoutExpired:
        return {"status": "killed", "obj": "", "bound": "", "gap": "", "nodes": "",
                "seconds": round(time.time() - started, 1)}
    seconds = round(time.time() - started, 1)
    grab = lambda pat: (re.search(pat, out).group(1) if re.search(pat, out) else "")
    status = "optimal" if "status:    optimal" in out else grab(r"status:\s+(\w+)") or "none"
    return {"status": status, "obj": grab(r"objective: ([-\d.e+]+)"),
            "bound": grab(r"bound ([-\d.e+]+)"), "gap": grab(r"gap ([\d.]+)%"),
            "nodes": grab(r"(\d+) nodes"), "seconds": seconds}


def run_highs(path):
    import highspy
    h = highspy.Highs()
    h.setOptionValue("output_flag", False)
    h.setOptionValue("threads", 1)
    h.setOptionValue("time_limit", LIMIT)
    started = time.time()
    try:
        h.readModel(str(path))
        h.run()
    except Exception:
        return {"status": "error", "obj": "", "seconds": round(time.time() - started, 1)}
    seconds = round(time.time() - started, 1)
    solved = str(h.getModelStatus()).endswith("kOptimal")
    info = h.getInfo()
    got = getattr(info, "primal_solution_status", 0)
    obj = h.getObjectiveValue() if got else None
    return {"status": "optimal" if solved else "limit",
            "obj": f"{obj:.6g}" if obj is not None else "", "seconds": seconds}


FIELDS = ["instance", "columns", "rows", "nonzeros", "integer_columns",
          "ripsolve_status", "ripsolve_objective", "ripsolve_bound", "ripsolve_gap_pct",
          "ripsolve_nodes", "ripsolve_seconds",
          "highs_status", "highs_objective", "highs_seconds", "agree", "note"]


def main():
    names = instance_names()
    random.seed(SEED)
    sample = sorted(random.sample(names, min(COUNT, len(names))))

    OUT.mkdir(parents=True, exist_ok=True)
    csv_path = OUT / "miplib_sample.csv"
    print(f"{len(sample)} of {len(names)} benchmark instances, seed {SEED}, {LIMIT:.0f}s limit")
    print(f"writing {csv_path}\n")

    with open(csv_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        for i, name in enumerate(sample, 1):
            row = dict.fromkeys(FIELDS, "")
            row["instance"] = name
            path, note = fetch(name)
            if path is None:
                row["note"] = note
                print(f"{i:3}/{len(sample)} {name:28} {note}")
                writer.writerow(row); handle.flush()
                continue

            dims, reason = shape(path)
            if dims is None:
                row["note"] = f"parse failed: {reason}"[:90]
                print(f"{i:3}/{len(sample)} {name:28} parse failed: {reason[:50]}")
                writer.writerow(row); handle.flush()
                continue
            cols, rows_, nnz = dims
            row.update(columns=cols, rows=rows_, nonzeros=nnz)

            r = run_ripsolve(path)
            h = run_highs(path)
            row.update(ripsolve_status=r["status"], ripsolve_objective=r["obj"],
                       ripsolve_bound=r["bound"], ripsolve_gap_pct=r["gap"],
                       ripsolve_nodes=r["nodes"], ripsolve_seconds=r["seconds"],
                       highs_status=h["status"], highs_objective=h["obj"],
                       highs_seconds=h["seconds"])
            # Only meaningful when both proved optimality.
            if r["status"] == "optimal" and h["status"] == "optimal" and r["obj"] and h["obj"]:
                a, b = float(r["obj"]), float(h["obj"])
                row["agree"] = "yes" if abs(a - b) <= 1e-6 * max(1.0, abs(b)) else "NO"
            print(f"{i:3}/{len(sample)} {name:28} {cols:>7}c {rows_:>7}r  "
                  f"ripsolve {r['status']:<10} {r['seconds']:>6}s   highs {h['status']:<8} {h['seconds']:>6}s")
            writer.writerow(row); handle.flush()

    print(f"\nwrote {csv_path}")


if __name__ == "__main__":
    sys.exit(main())
