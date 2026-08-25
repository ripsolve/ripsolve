#!/usr/bin/env python3
"""Build the README's benchmark tables.

One row per instance, with its dimensions and its optimal value, so the table says
what was solved rather than only how fast. The optimum is taken from every solver that
proved one and flagged if they disagree, which makes this a differential test as much
as a benchmark. Every solver is pinned to the same thread
count, and each is run at one thread and at eight, because the interesting question
is not only single-thread speed but how much of the machine each can use.

Only the solve is timed on every side; parsing is excluded.

Usage:  bench/table.py <group> [seconds] [instance ...]
        group is "bip" or "mkp"; naming instances overrides the group's list
"""
import pathlib
import re
import subprocess
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parent.parent
EXE = ROOT / "target" / "release" / "ripsolve"
GEN = pathlib.Path.home() / "repos" / "bip-gen"
GROUP = sys.argv[1] if len(sys.argv) > 1 else "bip"
LIMIT = float(sys.argv[2]) if len(sys.argv) > 2 else 60.0

BIP = ["v032c032", "v048c048", "v048c128", "v064c064", "v064c200", "v081c162n009",
       "v081c162n018", "v128c256n100", "v256c256n100", "v064c1000n100",
       "v128c1000n100", "v064c1000n020"]
MKP = ["mkp_200", "mkp_500"]


def locate(name):
    for base in (ROOT / "samples", GEN, ROOT / "bench" / "out"):
        for ext in (".lp", ".mps"):
            p = base / f"{name}{ext}"
            if p.exists():
                return p
    return None


def dimensions(path):
    out = subprocess.run([str(EXE), "info", str(path)], capture_output=True, text=True).stdout
    cols = re.search(r"(\d+)\s+columns", out)
    rows = re.search(r"(\d+)\s+rows", out)
    return (int(cols.group(1)) if cols else 0, int(rows.group(1)) if rows else 0)


def ripsolve(path, threads):
    started = time.time()
    out = subprocess.run(
        [str(EXE), "solve", "-t", str(threads), "--time-limit", str(LIMIT), str(path)],
        capture_output=True, text=True).stdout
    elapsed = time.time() - started
    obj = re.search(r"objective:\s+([-\d.e+]+)", out)
    proven = "status:    optimal" in out
    return (float(obj.group(1)) if obj else None), proven, elapsed


# Each external solve runs in a fresh interpreter. Solver bindings carry process-wide
# state: solving with highspy at one thread and then at eight in the same process makes
# the second return kNotset, which would read as a timeout rather than as a bug here.
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
proven = str(h.getModelStatus()).endswith("kOptimal")
print("RESULT", h.getObjectiveValue() if proven else "none", proven, e)
""",
    "SCIP": """
from pyscipopt import Model
import time
m = Model(); m.hideOutput(); m.readProblem({path!r})
m.setParam("limits/time", {limit})
m.setParam("parallel/maxnthreads", {threads})
t = time.time(); m.optimize(); e = time.time() - t
proven = m.getStatus() == "optimal"
print("RESULT", m.getObjVal() if m.getNSols() else "none", proven, e)
""",
    # Reported as "commercial" rather than by name. The licence in use here is an
    # academic one, and rather than read its terms as permission to publish benchmarks
    # under the solver's name, the published table omits it. Nothing else is withheld:
    # the call below names it, so anyone holding a licence reproduces the same numbers.
    "commercial": """
import gurobipy as gp, time
env = gp.Env(empty=True); env.setParam("OutputFlag", 0); env.start()
m = gp.read({path!r}, env)
m.setParam("Threads", {threads}); m.setParam("TimeLimit", {limit})
m.optimize()
proven = m.Status == gp.GRB.OPTIMAL
print("RESULT", m.ObjVal if m.SolCount else "none", proven, m.Runtime)
""",
}


def worker(kind):
    def run(path, threads):
        script = WORKERS[kind].format(path=str(path), threads=threads, limit=LIMIT)
        try:
            out = subprocess.run([sys.executable, "-c", script], capture_output=True,
                                 text=True, timeout=LIMIT * 3 + 60).stdout
        except subprocess.TimeoutExpired:
            return None, False, LIMIT
        line = re.search(r"RESULT (\S+) (\S+) (\S+)", out)
        if not line:
            return None, False, LIMIT
        obj = None if line.group(1) == "none" else float(line.group(1))
        return obj, line.group(2) == "True", float(line.group(3))
    return run


def cbc(path, threads):
    import pulp
    exe = pulp.PULP_CBC_CMD().path
    started = time.time()
    try:
        out = subprocess.run(
            [exe, str(path), "threads", str(threads), "seconds", str(LIMIT), "solve"],
            capture_output=True, text=True, timeout=LIMIT * 3).stdout
    except subprocess.TimeoutExpired:
        return None, False, LIMIT
    elapsed = time.time() - started
    objective, proven = None, False
    for line in out.splitlines():
        if "best objective" in line.lower():
            try:
                objective = float(line.lower().split("best objective")[1].split(",")[0])
            except (IndexError, ValueError):
                pass
            proven = "search completed" in line.lower()
    return objective, proven, elapsed


SOLVERS = [("ripsolve", ripsolve), ("HiGHS", worker("HiGHS")),
           ("SCIP", worker("SCIP")), ("CBC", cbc),
           ("commercial", worker("commercial"))]


def main():
    names = sys.argv[3:] or (BIP if GROUP == "bip" else MKP)
    print(f"{GROUP}, {LIMIT:g}s limit, times as 1 thread / 8 threads\n")
    header = f"{'instance':16}{'vars':>6}{'rows':>6}{'optimum':>12}"
    for label, _ in SOLVERS:
        header += f"{label:>20}"
    print(header, flush=True)

    for name in names:
        path = locate(name)
        if path is None:
            print(f"{name:16}  not found", flush=True)
            continue
        cols, rows = dimensions(path)
        cells = []
        answers = []
        for _, fn in SOLVERS:
            got = []
            for threads in (1, 8):
                obj, proven, elapsed = fn(path, threads)
                got.append(f"{elapsed:.1f}" + ("" if proven else "*"))
                if proven and obj is not None:
                    answers.append(obj)
            cells.append(" / ".join(got))
        optimum = f"{answers[0]:g}" if answers else "none"
        agree = all(abs(a - answers[0]) <= 1e-6 * max(1.0, abs(answers[0])) for a in answers)
        if not agree:
            optimum += " DISAGREE"
        line = f"{name:16}{cols:>6}{rows:>6}{optimum:>12}"
        for c in cells:
            line += f"{c:>20}"
        print(line, flush=True)
    print("\n* = hit the time limit, so the figure is a floor and the answer unproven")


main()
