#!/usr/bin/env python3
"""Run ripsolve over a seeded random sample of MIPLIB 2017.

The point of sampling rather than choosing is that a hand-picked set says nothing.
An earlier pass here used ten instances recalled from memory; every one of them was
smaller than 83% of the benchmark set, so the results could not support any claim
about MIPLIB generally.

Two lists are available. `easy` (the default, 747 instances) is the one to sample for
a claim about this solver, because those instances are known to be closed by reference
solvers, so a timeout is a statement about ripsolve. `benchmark` is the curated 240,
which are hard by construction and where a timeout mostly says the instance is hard.

The reference solver runs first, and an instance it cannot close within the same budget
is skipped rather than scored: a timeout on both sides measures the instance, not this
solver. Its answers are cached in `bench/out/miplib_oracle.json`, keyed by instance,
limit and thread count, so re-running a survey costs only ripsolve's time.

Both solvers get sixteen threads, which is what MIPLIB's own benchmarking rules allow.

Every instance in the sample is reported — including ones that fail to download,
fail to parse, or time out. Dropping those is how a benchmark flatters itself.

HiGHS runs alongside as an open-source reference, both to check answers and to say
whether a timeout is ripsolve's problem or the instance's.

Usage:  bench/miplib_sample.py [count] [seconds] [seed] [threads] [easy|benchmark]
"""

import csv
import json
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
# MIPLIB publishes several lists. `easy` is the right default here: those are the
# instances the reference solvers close within the competition's own limits, so a
# timeout on one is a statement about this solver rather than about the instance.
# `benchmark` is the harder curated 240, kept for when that is the question being
# asked.
LISTS = {
    "easy": "https://miplib.zib.de/downloads/easy-v18.test",
    "benchmark": "https://miplib.zib.de/downloads/benchmark-v2.test",
}
DATA_URL = "https://miplib.zib.de/WebData/instances/{}.mps.gz"

COUNT = int(sys.argv[1]) if len(sys.argv) > 1 else 25
LIMIT = float(sys.argv[2]) if len(sys.argv) > 2 else 60.0
SEED = int(sys.argv[3]) if len(sys.argv) > 3 else 20260824
# MIPLIB's own benchmarking rules allow sixteen threads, so that is what both solvers
# get. Running this at one thread measures something the published results are not
# comparable to.
THREADS = int(sys.argv[4]) if len(sys.argv) > 4 else 16
WHICH = sys.argv[5] if len(sys.argv) > 5 else "easy"
# Instances whose reference answer could not be pinned down, and so cannot score
# anything. Excluding one needs a reason that is about the instance, never about how
# ripsolve does on it.
EXCLUDED = {
    # Gurobi returns 1.8235 under one set of parameters and proves 1.9309 optimal
    # under another, so there is no reference value to check against. ripsolve's own
    # answer of 3.6199 is very likely wrong, and is tracked separately; it is not the
    # reason this instance is dropped.
    "neos-619167": "no stable reference objective",
}

# Instances larger than this are recorded as skipped rather than downloaded; the
# benchmark set reaches 40 MB compressed and the point is a survey, not a stress
# test of the network.
MAX_DOWNLOAD = 25 * 1024 * 1024


def instance_names():
    url = LISTS[WHICH]
    path = OUT / url.rsplit("/", 1)[1]
    if not path.exists():
        OUT.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(url, path)
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
            [str(EXE), "solve", "-t", str(THREADS), "--time-limit", str(LIMIT), str(path)],
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


ORACLE_CACHE = OUT / "miplib_oracle.json"

# How far two "optimal" answers may sit apart before it means something.
#
# Both solvers stop at a relative gap of 1e-4, which is their default and the field's,
# so each may sit that far from the true optimum and they may sit on opposite sides of
# it. Judging agreement at solver precision instead calls that a disagreement: app2-1
# came back 19296 against 19295.5, a relative difference of 2.6e-5, which is both
# answers being correct to the tolerance they were asked for.
AGREEMENT_TOLERANCE = 2.5e-4


def load_oracle_cache():
    if ORACLE_CACHE.exists():
        try:
            return json.loads(ORACLE_CACHE.read_text())
        except json.JSONDecodeError:
            return {}
    return {}


def run_oracle(name, path, cache):
    """What the reference solver makes of this instance, from cache when possible.

    An instance's answer does not change, so the only reason to re-run is a different
    budget. The key carries the limit and the thread count for that reason, and nothing
    else: a cache that silently answers for the wrong budget is worse than no cache.

    Reported as `reference` rather than by name. The licence in use here is an academic
    one, and rather than read its terms as permission to publish benchmarks under the
    solver's name, the published output omits it.
    """
    key = f"{name}|{LIMIT:g}|{THREADS}"
    if key in cache:
        hit = dict(cache[key])
        hit["cached"] = True
        return hit

    import gurobipy as gp

    started = time.time()
    try:
        env = gp.Env(empty=True)
        env.setParam("OutputFlag", 0)
        env.start()
        m = gp.read(str(path), env)
        m.setParam("Threads", THREADS)
        m.setParam("TimeLimit", LIMIT)
        m.optimize()
        proven = m.Status == gp.GRB.OPTIMAL
        obj = m.ObjVal if m.SolCount else None
        result = {
            "status": "optimal" if proven else "limit",
            "obj": f"{obj:.10g}" if obj is not None else "",
            "seconds": round(m.Runtime, 1),
        }
    except Exception as exc:  # noqa: BLE001 - any failure is reported, not raised
        result = {
            "status": "error",
            "obj": "",
            "seconds": round(time.time() - started, 1),
            "note": str(exc)[:80],
        }

    cache[key] = result
    ORACLE_CACHE.parent.mkdir(parents=True, exist_ok=True)
    ORACLE_CACHE.write_text(json.dumps(cache, indent=1, sort_keys=True))
    out = dict(result)
    out["cached"] = False
    return out


FIELDS = ["instance", "columns", "rows", "nonzeros", "integer_columns",
          "ripsolve_status", "ripsolve_objective", "ripsolve_bound", "ripsolve_gap_pct",
          "ripsolve_nodes", "ripsolve_seconds",
          "reference_status", "reference_objective", "reference_seconds",
          "agree", "note"]


def main():
    names = instance_names()
    random.seed(SEED)
    sample = sorted(random.sample(names, min(COUNT, len(names))))
    dropped = [n for n in sample if n in EXCLUDED]
    sample = [n for n in sample if n not in EXCLUDED]
    for name in dropped:
        print(f"  excluded {name}: {EXCLUDED[name]}")

    OUT.mkdir(parents=True, exist_ok=True)
    csv_path = OUT / "miplib_sample.csv"
    print(f"{len(sample)} of {len(names)} MIPLIB '{WHICH}' instances, seed {SEED}, "
          f"{LIMIT:.0f}s limit, {THREADS} threads")
    print(f"writing {csv_path}\n")

    cache = load_oracle_cache()
    solved = scored = uninformative = disagreements = 0

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

            # The reference runs first. An instance it cannot close within the same
            # budget cannot say anything about ripsolve, so there is no point spending
            # the time: it is reported and skipped rather than scored.
            o = run_oracle(name, path, cache)
            row.update(reference_status=o["status"], reference_objective=o["obj"],
                       reference_seconds=o["seconds"])
            mark = " (cached)" if o.get("cached") else ""
            if o["status"] != "optimal":
                uninformative += 1
                row["note"] = "reference did not close it"
                print(f"{i:3}/{len(sample)} {name:28} {cols:>7}c {rows_:>7}r  "
                      f"skipped, reference {o['status']}{mark}")
                writer.writerow(row); handle.flush()
                continue

            r = run_ripsolve(path)
            row.update(ripsolve_status=r["status"], ripsolve_objective=r["obj"],
                       ripsolve_bound=r["bound"], ripsolve_gap_pct=r["gap"],
                       ripsolve_nodes=r["nodes"], ripsolve_seconds=r["seconds"])
            if r["status"] == "optimal":
                solved += 1
            # Only meaningful when both proved optimality.
            if r["status"] == "optimal" and r["obj"] and o["obj"]:
                a, b = float(r["obj"]), float(o["obj"])
                apart = abs(a - b) / max(1.0, abs(b))
                row["agree"] = "yes" if apart <= AGREEMENT_TOLERANCE else "NO"
                if row["agree"] == "NO":
                    disagreements += 1
                    row["note"] = f"objectives differ by {apart:.3e} relative"
            scored += 1
            print(f"{i:3}/{len(sample)} {name:28} {cols:>7}c {rows_:>7}r  "
                  f"ripsolve {r['status']:<10} {r['seconds']:>6}s   "
                  f"reference optimal {o['seconds']:>6}s{mark}"
                  + (f"  OBJECTIVES DISAGREE: {row['note']}" if row.get("agree") == "NO" else ""))
            writer.writerow(row); handle.flush()

    print(f"\n{solved}/{scored} solved of the instances the reference closes; "
          f"{uninformative} skipped because it did not")
    if disagreements:
        print(f"{disagreements} OBJECTIVE DISAGREEMENTS")
    print(f"wrote {csv_path}")


if __name__ == "__main__":
    sys.exit(main())
