#!/usr/bin/env python3
"""This solver alone over MIPLIB 2017's 240 benchmark instances, checked against MIPLIB.

The benchmark set is the one Mittelmann's MILP comparison runs, and unlike the pure
binary set it is mostly general MIP: integers and continuous columns as well as
binaries. Every claimed answer is checked against MIPLIB's published solution file, so
a run reports wrong answers separately from slow ones: an optimum that disagrees with
MIPLIB's beyond the relative gap the solve stopped at, an infeasibility proof for an
instance MIPLIB has solved, or a model that does not read at all.

Usage:  bench/benchmark_set.py [seconds] [threads]

Results are cached per instance and budget in `bench/out/benchmark_cache.json`, so an
interrupted run resumes where it stopped; delete the file after changing the binary.
"""

import gzip
import json
import pathlib
import re
import shutil
import subprocess
import sys
import time
import urllib.request

ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = ROOT / "bench" / "out"
MODELS = OUT / "miplib"
LIST = OUT / "benchmark-v2.test"
SOLU = OUT / "miplib2017-v35.solu"
CACHE = OUT / "benchmark_cache.json"
DATA_URL = "https://miplib.zib.de/WebData/instances/{}.mps.gz"
SOLU_URL = "https://miplib.zib.de/downloads/miplib2017-v35.solu"
LIMIT = float(sys.argv[1]) if len(sys.argv) > 1 else 60.0
THREADS = int(sys.argv[2]) if len(sys.argv) > 2 else 12
# The CLI's default relative gap. An optimum within it of MIPLIB's is the answer the
# solve was asked for, as it would be from any solver run at its default tolerance.
GAP = 1e-4


def fetch(name):
    """The uncompressed model, downloaded on first use. Unlike the screening scripts
    this takes models of any size: the set is fixed, and skipping its large instances
    would make the count look better than it is."""
    mps = MODELS / f"{name}.mps"
    if mps.exists():
        return mps
    MODELS.mkdir(parents=True, exist_ok=True)
    gz = MODELS / f"{name}.mps.gz"
    with urllib.request.urlopen(DATA_URL.format(name), timeout=600) as response:
        with open(gz, "wb") as dst:
            shutil.copyfileobj(response, dst)
    with gzip.open(gz, "rb") as src, open(mps, "wb") as dst:
        shutil.copyfileobj(src, dst)
    gz.unlink()
    return mps


def known():
    """MIPLIB's status and objective for every instance it has one for."""
    if not SOLU.exists():
        urllib.request.urlretrieve(SOLU_URL, SOLU)
    table = {}
    for line in SOLU.read_text().splitlines():
        parts = line.split()
        if len(parts) >= 2:
            table[parts[1]] = (parts[0], float(parts[2]) if len(parts) > 2 else None)
    return table


def run(path):
    started = time.time()
    try:
        done = subprocess.run(
            [str(ROOT / "target" / "release" / "ripsolve"), "solve",
             "-t", str(THREADS), "--time-limit", str(LIMIT), str(path)],
            capture_output=True, text=True, timeout=LIMIT * 4 + 300)
    except subprocess.TimeoutExpired:
        return {"status": "killed", "seconds": round(time.time() - started, 1)}
    elapsed = round(time.time() - started, 1)
    if done.returncode != 0:
        last = (done.stderr.strip().splitlines() or ["no message"])[-1]
        return {"status": "error", "seconds": elapsed, "message": last[:300]}
    status = re.search(r"^status:\s+(\S+)", done.stdout, re.M)
    objective = re.search(r"^objective: (\S+)", done.stdout, re.M)
    bound = re.search(r"\(bound (\S+?),", done.stdout)
    return {"status": status.group(1) if status else "unparsed", "seconds": elapsed,
            "objective": float(objective.group(1)) if objective else None,
            "bound": float(bound.group(1)) if bound else None}


def verdict(result, reference):
    """closed, open, or wrong: an answer that contradicts MIPLIB's."""
    tag, value = reference or (None, None)
    status = result["status"]
    if status in ("error", "killed", "unparsed"):
        return "wrong" if status == "error" else "open"
    if status == "Infeasible":
        return "wrong" if tag in ("=opt=", "=best=") else "closed"
    if tag == "=inf=" and result.get("objective") is not None:
        return "wrong"
    if status == "optimal":
        if tag != "=opt=" or result.get("objective") is None:
            return "closed"
        scale = max(1.0, abs(value))
        return "closed" if abs(result["objective"] - value) <= GAP * scale + 1e-6 else "wrong"
    # An open run's bound must not pass the optimum. The sense is read off the side of
    # the optimum the incumbent is on; with no incumbent, or one at the optimum, there
    # is nothing to read it from.
    objective, bound = result.get("objective"), result.get("bound")
    if tag == "=opt=" and objective is not None and bound is not None:
        slack = GAP * max(1.0, abs(value)) + 1e-6
        if objective > value + slack and bound > value + slack:
            return "wrong"
        if objective < value - slack and bound < value - slack:
            return "wrong"
    return "open"


def main():
    names = [l.strip().removesuffix(".mps.gz") for l in LIST.read_text().splitlines()
             if l.strip()]
    reference = known()
    cache = json.loads(CACHE.read_text()) if CACHE.exists() else {}
    print(f"{len(names)} benchmark instances, {LIMIT:g}s, {THREADS} threads\n", flush=True)
    rows = []
    for name in names:
        key = f"{name}|{LIMIT:g}|{THREADS}"
        if key not in cache:
            try:
                path = fetch(name)
            except Exception as error:
                cache[key] = {"status": "unfetched", "seconds": 0.0, "message": str(error)}
            else:
                cache[key] = run(path)
            CACHE.write_text(json.dumps(cache, indent=1, sort_keys=True))
        result = dict(cache[key])
        result["verdict"] = ("open" if result["status"] == "unfetched"
                             else verdict(result, reference.get(name)))
        result["miplib"] = reference.get(name)
        rows.append((name, result))
        print("%-32s %-8s %-12s %7.1fs %s" % (
            name, result["verdict"], result["status"], result["seconds"],
            result.get("message", "")), flush=True)

    tally = {v: sum(1 for _, r in rows if r["verdict"] == v)
             for v in ("closed", "open", "wrong")}
    print("\n" + ", ".join(f"{v} {n}" for v, n in tally.items()))
    out = OUT / f"benchmark_results_{LIMIT:g}s.json"
    out.write_text(json.dumps({"rows": rows, "tally": tally, "limit": LIMIT,
                               "threads": THREADS}, indent=1))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
