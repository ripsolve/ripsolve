#!/usr/bin/env python3
"""Randomized differential test of ripsolve against Gurobi on mixed-integer models.

Generates models mixing binary, general-integer and continuous columns, solves each
with both, and requires the objectives to agree. Structure is randomized rather than
fixed because the interesting failures are in the combinations -- an integer column
next to a continuous one in the same row is exactly where a binary-only assumption
would survive unnoticed.
"""
import random
import subprocess
import sys
import pathlib

ROOT = pathlib.Path(__file__).resolve().parent.parent
EXE = ROOT / "target" / "release" / "ripsolve"
OUT = pathlib.Path("/tmp/claude-1000/-home-andy-repos/784cbd5f-e5c1-4a72-b4b8-e191a18aa741/scratchpad/mipfuzz")
OUT.mkdir(parents=True, exist_ok=True)


def generate(seed, path):
    rnd = random.Random(seed)
    n = rnd.randint(4, 12)
    m = rnd.randint(2, 6)
    # Each column is binary, a bounded general integer, or continuous.
    kinds = [rnd.choice(["B", "I", "C"]) for _ in range(n)]
    ubs = [1 if k == "B" else rnd.randint(2, 8) for k in kinds]
    obj = [rnd.randint(1, 9) + (0.5 if k == "C" else 0) for k in kinds]

    lines = ["Minimize", " obj: " + " + ".join(f"{obj[j]} x{j}" for j in range(n)), "Subject To"]
    for i in range(m):
        coeffs = [rnd.randint(1, 6) for _ in range(n)]
        # A right-hand side reachable but not free.
        cap = sum(coeffs[j] * ubs[j] for j in range(n))
        lines.append(f" c{i}: " + " + ".join(f"{coeffs[j]} x{j}" for j in range(n))
                     + f" >= {rnd.randint(cap // 5, cap // 2)}")
    lines.append("Bounds")
    for j in range(n):
        lines.append(f" 0 <= x{j} <= {ubs[j]}")
    integral = [f"x{j}" for j in range(n) if kinds[j] in ("B", "I")]
    if integral:
        lines.append("General")
        lines.append(" " + " ".join(integral))
    lines.append("End")
    path.write_text("\n".join(lines) + "\n")
    return kinds


def gurobi(path):
    import gurobipy as gp
    env = gp.Env(params={"OutputFlag": 0})
    m = gp.read(str(path), env=env)
    m.setParam("Threads", 1)
    m.optimize()
    return (m.ObjVal if m.SolCount else None), m.Status


def ripsolve(path):
    out = subprocess.run([str(EXE), "solve", "-t", "1", str(path)],
                         capture_output=True, text=True).stdout
    if "status:    infeasible" in out.lower() or "Infeasible" in out:
        return None, "INFEASIBLE"
    for line in out.splitlines():
        if line.startswith("objective:"):
            return float(line.split()[1]), "OPTIMAL"
    return None, out.strip()[:60]


def main():
    trials = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    mismatches = 0
    for seed in range(trials):
        path = OUT / f"m{seed}.lp"
        generate(seed, path)
        mine, mine_status = ripsolve(path)
        theirs, _ = gurobi(path)
        if mine is None and theirs is None:
            continue
        if mine is None or theirs is None or abs(mine - theirs) > 1e-6 * max(1, abs(theirs)):
            mismatches += 1
            print(f"MISMATCH seed {seed}: ripsolve {mine} ({mine_status}) vs gurobi {theirs}")
            print(path.read_text())
            if mismatches >= 3:
                break
    print(f"{trials - mismatches}/{trials} agree")
    return 1 if mismatches else 0


if __name__ == "__main__":
    sys.exit(main())
