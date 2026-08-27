#!/usr/bin/env python3
"""Every solver on the sampled MIPLIB instances, in one table.

The per-instance surveys compare ripsolve against one reference at a time, which says
whether an instance is lost but not how far off the open-source field it is. This runs
the lot so the shape of the deficit is visible.

ripsolve and the commercial solver are read from the survey's own outputs, since both
have already been measured at this budget; HiGHS, SCIP and CBC are run here. Every
external solve runs in a fresh interpreter, because the bindings carry process-wide
state: solving with highspy at one thread and then at eight in the same process makes
the second return kNotset.

Usage:  bench/miplib_table.py [seconds] [threads]
"""

import csv
import json
import pathlib
import re
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = ROOT / "bench" / "out"
CACHE = OUT / "miplib"
LIMIT = float(sys.argv[1]) if len(sys.argv) > 1 else 60.0
THREADS = int(sys.argv[2]) if len(sys.argv) > 2 else 16

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
}


def run_external(kind, path):
    script = WORKERS[kind].format(path=str(path), threads=THREADS, limit=LIMIT)
    try:
        out = subprocess.run([sys.executable, "-c", script], capture_output=True,
                             text=True, timeout=LIMIT * 3 + 120).stdout
    except subprocess.TimeoutExpired:
        return "killed", LIMIT * 3
    found = re.search(r"RESULT (\S+) (\S+)", out)
    if not found:
        return "error", 0.0
    return found.group(1), float(found.group(2))


def run_ripsolve(path):
    import time
    started = time.time()
    try:
        out = subprocess.run(
            [str(ROOT / "target" / "release" / "ripsolve"), "solve",
             "-t", str(THREADS), "--time-limit", str(LIMIT), str(path)],
            capture_output=True, text=True, timeout=LIMIT * 4 + 120).stdout
    except subprocess.TimeoutExpired:
        return "killed", LIMIT * 4
    elapsed = time.time() - started
    if "status:    optimal" in out:
        return "optimal", elapsed
    if "status:    Infeasible" in out:
        return "Infeasible", elapsed
    return "limit", elapsed


def run_cbc(path):
    import pulp
    exe = pulp.PULP_CBC_CMD().path
    import time
    started = time.time()
    try:
        out = subprocess.run([exe, str(path), "threads", str(THREADS),
                              "seconds", str(LIMIT), "solve"],
                             capture_output=True, text=True, timeout=LIMIT * 3 + 120).stdout
    except subprocess.TimeoutExpired:
        return "killed", LIMIT * 3
    elapsed = time.time() - started
    proven = any("search completed" in l.lower() for l in out.splitlines())
    return ("optimal" if proven else "limit"), elapsed


def main():
    survey = {r["instance"]: r for r in csv.DictReader(open(OUT / "miplib_sample.csv"))}
    oracle = json.loads((OUT / "miplib_oracle.json").read_text())

    def short(status):
        return {"optimal": "opt", "TimeLimit": "lim", "limit": "lim",
                "NodeLimit": "node", "Infeasible": "infeas", "killed": "kill",
                "timelimit": "lim", "userinterrupt": "lim"}.get(status, status[:6])

    print(f"MIPLIB easy sample, {LIMIT:g}s limit, {THREADS} threads. "
          f"opt = proved optimality.\n")
    header = f"{'instance':22}{'rows':>7}{'cols':>7}"
    for name in ("ripsolve", "HiGHS", "SCIP", "CBC", "commercial"):
        header += f"{name:>13}"
    print(header)
    print("-" * len(header))

    tally = {n: 0 for n in ("ripsolve", "HiGHS", "SCIP", "CBC", "commercial")}
    counted = 0
    for name in sorted(survey):
        row = survey[name]
        path = CACHE / f"{name}.mps"
        if not path.exists() or not row["rows"]:
            print(f"{name:22}  {row.get('note', 'unavailable')[:60]}")
            continue

        cells = {}
        if row["ripsolve_status"]:
            cells["ripsolve"] = (short(row["ripsolve_status"]),
                                 float(row["ripsolve_seconds"] or 0))
        else:
            # The survey skips instances the reference cannot close, so those rows
            # carry no ripsolve result. Reporting them as a failure would be a lie in
            # our favour or against us depending on the reader, so they are run here.
            status, seconds = run_ripsolve(path)
            cells["ripsolve"] = (short(status), seconds)
        key = f"{name}|{LIMIT:g}|{THREADS}"
        if key in oracle:
            cells["commercial"] = (short(oracle[key]["status"]), oracle[key]["seconds"])
        else:
            cells["commercial"] = ("?", 0.0)
        for kind in ("HiGHS", "SCIP"):
            status, seconds = run_external(kind, path)
            cells[kind] = (short(status), seconds)
        status, seconds = run_cbc(path)
        cells["CBC"] = (short(status), seconds)

        counted += 1
        line = f"{name:22}{row['rows']:>7}{row['columns']:>7}"
        for solver in ("ripsolve", "HiGHS", "SCIP", "CBC", "commercial"):
            status, seconds = cells[solver]
            if status == "opt":
                tally[solver] += 1
            line += f"{status + ' ' + format(seconds, '.1f'):>13}"
        print(line, flush=True)

    print("-" * len(header))
    summary = f"{'solved of ' + str(counted):22}{'':>14}"
    for solver in ("ripsolve", "HiGHS", "SCIP", "CBC", "commercial"):
        summary += f"{tally[solver]:>13}"
    print(summary)


main()
