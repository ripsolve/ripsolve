#!/usr/bin/env python3
"""Build a benchmark of MIPLIB instances the open-source solvers can actually close.

An instance nobody solves measures the instance. An instance only a commercial solver
solves measures the gap between commercial and open-source work, which is not this
solver's problem to answer for. What is left, the instances at least two of HiGHS, SCIP
and CBC close within the budget, is a set where every failure is ours and every one is
demonstrably reachable.

Screening is cached per (instance, solver, limit, threads), so a re-run costs only what
has not been measured. Every external solve runs in a fresh interpreter: the bindings
carry process-wide state, and highspy in particular returns kNotset for a second solve
in the same process.

Usage:  bench/tractable.py [count] [seconds] [threads] [seed]
        count is how many instances to screen, not how many qualify
"""

import gzip
import json
import pathlib
import random
import re
import shutil
import subprocess
import sys
import time
import urllib.request

ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = ROOT / "bench" / "out"
CACHE = OUT / "miplib"
SCREEN = OUT / "tractable_screen.json"
LIST_URL = "https://miplib.zib.de/downloads/easy-v18.test"
DATA_URL = "https://miplib.zib.de/WebData/instances/{}.mps.gz"

# Read only when run directly. Importing this module must not consume another
# program's arguments, and must not quietly screen at a different budget: the cache is
# keyed by budget, so it would measure honestly and answer the wrong question.
RUN_DIRECTLY = __name__ == "__main__"
COUNT = int(sys.argv[1]) if RUN_DIRECTLY and len(sys.argv) > 1 else 60
LIMIT = float(sys.argv[2]) if RUN_DIRECTLY and len(sys.argv) > 2 else 60.0
THREADS = int(sys.argv[3]) if RUN_DIRECTLY and len(sys.argv) > 3 else 16
SEED = int(sys.argv[4]) if RUN_DIRECTLY and len(sys.argv) > 4 else 20260826
MAX_DOWNLOAD = 25 * 1024 * 1024

WORKERS = {
    "HiGHS": """
import highspy, time
h = highspy.Highs()
h.setOptionValue("output_flag", False)
h.setOptionValue("threads", {threads})
h.setOptionValue("parallel", "on" if {threads} > 1 else "off")
h.setOptionValue("time_limit", {limit})
h.readModel({path!r})
t = time.time(); h.run(); e = time.time() - t
print("RESULT", "optimal" if str(h.getModelStatus()).endswith("kOptimal") else "limit", e)
""",
    "SCIP": """
from pyscipopt import Model
import time
m = Model(); m.hideOutput(); m.readProblem({path!r})
m.setParam("limits/time", {limit})
m.setParam("parallel/maxnthreads", {threads})
t = time.time(); m.optimize(); e = time.time() - t
print("RESULT", "optimal" if m.getStatus() == "optimal" else m.getStatus(), e)
""",
    "commercial": """
import gurobipy as gp
env = gp.Env(empty=True); env.setParam("OutputFlag", 0); env.start()
m = gp.read({path!r}, env)
m.setParam("Threads", {threads}); m.setParam("TimeLimit", {limit})
m.optimize()
print("RESULT", "optimal" if m.Status == gp.GRB.OPTIMAL else "limit", m.Runtime)
""",
}


def instance_names():
    path = OUT / "easy-v18.test"
    if not path.exists():
        urllib.request.urlretrieve(LIST_URL, path)
    return [l.strip().removesuffix(".mps.gz") for l in path.read_text().splitlines() if l.strip()]


def fetch(name):
    CACHE.mkdir(parents=True, exist_ok=True)
    mps = CACHE / f"{name}.mps"
    if mps.exists():
        return mps
    gz = CACHE / f"{name}.mps.gz"
    try:
        with urllib.request.urlopen(DATA_URL.format(name), timeout=120) as response:
            size = int(response.headers.get("content-length", 0))
            if size > MAX_DOWNLOAD:
                return None
            gz.write_bytes(response.read())
        with gzip.open(gz, "rb") as src, open(mps, "wb") as dst:
            shutil.copyfileobj(src, dst)
        gz.unlink(missing_ok=True)
        return mps
    except Exception:
        return None


def load(path):
    if path.exists():
        try:
            return json.loads(path.read_text())
        except json.JSONDecodeError:
            return {}
    return {}


def measure(kind, name, path, cache):
    key = f"{name}|{kind}|{LIMIT:g}|{THREADS}"
    if key in cache:
        return cache[key]
    if kind == "CBC":
        result = run_cbc(path)
    elif kind == "ripsolve":
        result = run_ripsolve(path)
    else:
        script = WORKERS[kind].format(path=str(path), threads=THREADS, limit=LIMIT)
        try:
            out = subprocess.run([sys.executable, "-c", script], capture_output=True,
                                 text=True, timeout=LIMIT * 3 + 120).stdout
            found = re.search(r"RESULT (\S+) (\S+)", out)
            result = ({"status": found.group(1), "seconds": round(float(found.group(2)), 1)}
                      if found else {"status": "error", "seconds": 0.0})
        except subprocess.TimeoutExpired:
            result = {"status": "killed", "seconds": LIMIT * 3}
    cache[key] = result
    SCREEN.write_text(json.dumps(cache, indent=1, sort_keys=True))
    return result


def run_cbc(path):
    import pulp
    exe = pulp.PULP_CBC_CMD().path
    started = time.time()
    try:
        out = subprocess.run([exe, str(path), "threads", str(THREADS),
                              "seconds", str(LIMIT), "solve"],
                             capture_output=True, text=True, timeout=LIMIT * 3 + 120).stdout
    except subprocess.TimeoutExpired:
        return {"status": "killed", "seconds": LIMIT * 3}
    proven = any("search completed" in l.lower() for l in out.splitlines())
    return {"status": "optimal" if proven else "limit",
            "seconds": round(time.time() - started, 1)}


def run_ripsolve(path):
    started = time.time()
    try:
        out = subprocess.run(
            [str(ROOT / "target" / "release" / "ripsolve"), "solve",
             "-t", str(THREADS), "--time-limit", str(LIMIT), str(path)],
            capture_output=True, text=True, timeout=LIMIT * 4 + 120).stdout
    except subprocess.TimeoutExpired:
        return {"status": "killed", "seconds": LIMIT * 4}
    elapsed = round(time.time() - started, 1)
    if "status:    optimal" in out:
        return {"status": "optimal", "seconds": elapsed}
    if "status:    Infeasible" in out:
        return {"status": "infeasible", "seconds": elapsed}
    return {"status": "limit", "seconds": elapsed}


def main():
    names = instance_names()
    random.seed(SEED)
    sample = sorted(random.sample(names, min(COUNT, len(names))))
    cache = load(SCREEN)

    print(f"screening {len(sample)} of {len(names)} MIPLIB 'easy' instances, "
          f"{LIMIT:g}s, {THREADS} threads, seed {SEED}")
    print("keeping those at least two of HiGHS, SCIP and CBC close\n")

    qualifying = []
    for i, name in enumerate(sample, 1):
        path = fetch(name)
        if path is None:
            print(f"{i:3}/{len(sample)} {name:28} unavailable")
            continue
        open_source = {k: measure(k, name, path, cache) for k in ("HiGHS", "SCIP", "CBC")}
        closed = [k for k, v in open_source.items() if v["status"] == "optimal"]
        mark = "KEEP" if len(closed) >= 2 else "drop"
        print(f"{i:3}/{len(sample)} {name:28} {mark}  closed by {len(closed)}: "
              f"{', '.join(closed) if closed else 'none'}", flush=True)
        if len(closed) >= 2:
            qualifying.append(name)

    print(f"\n{len(qualifying)} instances qualify:")
    for name in qualifying:
        print(f"  {name}")
    (OUT / "tractable.json").write_text(json.dumps(qualifying, indent=1))
    print(f"\nwrote {OUT / 'tractable.json'}")


if __name__ == "__main__":
    main()
