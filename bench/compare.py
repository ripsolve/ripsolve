#!/usr/bin/env python3
"""Head-to-head comparison of ripsolve against Gurobi.

Fairness notes, all deliberate:
  * Both solvers get the same thread count, passed as the second argument. At one
    thread this compares algorithms; at the machine's full width it compares what a
    user would actually experience.
  * Only the solve is timed on both sides -- gurobipy's Runtime excludes parsing,
    and ripsolve's reported elapsed covers search only.
  * The same time limit applies to both, and a run that hits it is reported with
    its remaining gap rather than as a solve.
"""
import re
import subprocess
import sys
import pathlib

ROOT = pathlib.Path(__file__).resolve().parent.parent
EXE = ROOT / "target" / "release" / "ripsolve"
LIMIT = float(sys.argv[1]) if len(sys.argv) > 1 else 60.0
THREADS = int(sys.argv[2]) if len(sys.argv) > 2 else 1

BIPGEN = pathlib.Path.home() / "repos" / "bip-gen"
SCRATCH = pathlib.Path("/tmp/claude-1000/-home-andy-repos/784cbd5f-e5c1-4a72-b4b8-e191a18aa741/scratchpad/bench")

INSTANCES = [
    ("v032c032",      ROOT / "samples/v032c032.lp"),
    ("v048c048",      ROOT / "samples/v048c048.lp"),
    ("v048c128",      ROOT / "samples/v048c128.lp"),
    ("v064c064",      ROOT / "samples/v064c064.lp"),
    ("v064c200",      ROOT / "samples/v064c200.mps"),
    ("v081c162n009",  BIPGEN / "v081c162n009.lp"),
    ("v081c162n018",  BIPGEN / "v081c162n018.lp"),
    ("v128c256n100",  BIPGEN / "v128c256n100.lp"),
    ("v256c256n100",  BIPGEN / "v256c256n100.lp"),
    ("v064c1000n100", BIPGEN / "v064c1000n100.lp"),
    ("v064c1000n020", BIPGEN / "v064c1000n020.lp"),
    ("v128c1000n100", BIPGEN / "v128c1000n100.lp"),
    ("mkp_200",       SCRATCH / "mkp_200.lp"),
    ("mkp_500",       SCRATCH / "mkp_500.lp"),
]


def gurobi(path):
    import gurobipy as gp
    env = gp.Env(params={"OutputFlag": 0})
    m = gp.read(str(path), env=env)
    m.setParam("Threads", THREADS)
    m.setParam("TimeLimit", LIMIT)
    m.setParam("MIPGap", 0.0)
    m.optimize()
    solved = m.Status == gp.GRB.OPTIMAL
    obj = m.ObjVal if m.SolCount > 0 else None
    gap = m.MIPGap if m.SolCount > 0 else float("inf")
    return dict(solved=solved, obj=obj, nodes=int(m.NodeCount),
                time=m.Runtime, gap=gap)


def ripsolve(path):
    out = subprocess.run(
        [str(EXE), "solve", "--time-limit", str(LIMIT), "-t", str(THREADS), str(path)],
        capture_output=True, text=True, timeout=LIMIT * 3,
    ).stdout
    obj = re.search(r"objective: ([-\d.e+]+)", out)
    nodes = re.search(r"(\d+) nodes", out)
    time = re.search(r"([\d.]+)(m?s)\b(?!.*nodes)", out.splitlines()[-1])
    gap = re.search(r"gap ([\d.]+)%", out)
    solved = "status:    optimal" in out
    secs = None
    if time:
        secs = float(time.group(1)) / (1000.0 if time.group(2) == "ms" else 1.0)
    return dict(solved=solved,
                obj=float(obj.group(1)) if obj else None,
                nodes=int(nodes.group(1)) if nodes else 0,
                time=secs if secs is not None else LIMIT,
                gap=float(gap.group(1)) / 100.0 if gap else (0.0 if solved else float("inf")))


def cell(r):
    if r["obj"] is None:
        return "no solution"
    if r["solved"]:
        return f"{r['obj']:.0f}"
    return f"{r['obj']:.0f} (gap {r['gap']*100:.0f}%)"


print(f"time limit {LIMIT:.0f}s, both solvers on {THREADS} thread(s)\n")
hdr = f"{'instance':16} | {'ripsolve':>26} | {'Gurobi':>26} | ratio"
print(hdr)
print("-" * len(hdr))
for name, path in INSTANCES:
    if not pathlib.Path(path).exists():
        continue
    r = ripsolve(path)
    g = gurobi(path)
    ratio = ""
    if r["solved"] and g["solved"]:
        ratio = f"{r['time'] / max(g['time'], 1e-3):.0f}x"
    elif g["solved"] and not r["solved"]:
        ratio = ">limit"
    elif r["solved"] and not g["solved"]:
        ratio = "ripsolve wins"
    print(f"{name:16} | {cell(r):>14} {r['nodes']:>7}n {r['time']:>4.0f}s"
          f" | {cell(g):>14} {g['nodes']:>7}n {g['time']:>4.0f}s | {ratio}")
