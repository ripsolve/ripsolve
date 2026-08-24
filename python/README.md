# ripsolve for Python

A `gurobipy`-shaped interface for binary integer programs. A gurobipy script that
uses only binary variables should run unchanged after swapping the import.

```python
import ripsolve as gp
from ripsolve import GRB

m = gp.Model("knapsack")
x = m.addVars(n, vtype=GRB.BINARY, name="x")
m.setObjective(gp.quicksum(value[j] * x[j] for j in range(n)), GRB.MAXIMIZE)
m.addConstr(gp.quicksum(weight[j] * x[j] for j in range(n)) <= capacity)
m.optimize()

print(m.ObjVal, [x[j].X for j in range(n)])
```

## Building

```sh
./build.sh                      # writes ripsolve.so beside this file
python3 -m pytest tests -q
```

Needs a Python with development headers. No maturin required.

## What is supported

| | |
|---|---|
| Model | `Model(name)`, `optimize`, `getVars`, `write`, `setParam` |
| Variables | `addVar`, `addVars`, `.X`, `.VarName`, `.Obj`, `.LB`, `.UB` |
| Expressions | `+ - *`, unary `-`, `<= >= ==`, `quicksum` |
| Constraints | `addConstr`, `addConstrs`, `.ConstrName` |
| Objective | `setObjective(expr, GRB.MINIMIZE / GRB.MAXIMIZE)` |
| Attributes | `ObjVal`, `ObjBound`, `Status`, `SolCount`, `NodeCount`, `Runtime`, `MIPGap`, `NumVars`, `NumConstrs`, `ModelName`, `ModelSense` |
| Parameters | `TimeLimit`, `Threads`, `MIPGap`, `OutputFlag` |
| Files | `read(path)` for LP and MPS |

## What is not

Every variable is binary. `addVar` takes `vtype` so that source stays portable, but
rejects anything other than `GRB.BINARY` rather than silently relaxing a model the
solver cannot honour. There are no callbacks, no quadratic terms, no SOS
constraints, no lazy constraints, and no multi-objective support.

Unknown parameter names raise `KeyError` rather than being ignored, so a
misspelling fails loudly.

## Threads

`optimize()` releases the GIL, so a solve does not block other Python threads.
`Threads` defaults to the machine's parallelism, as Gurobi's does.
