"""The Python API, checked against gurobipy where it is available.

The bar is that a gurobipy script using only binary variables runs unchanged
after swapping the import, and returns the same answer. Where gurobipy is not
installed the comparisons skip and the rest still runs.
"""

import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import ripsolve
from ripsolve import GRB

try:
    import gurobipy
    from gurobipy import GRB as GGRB

    HAVE_GUROBI = True
except ImportError:  # pragma: no cover - environment dependent
    HAVE_GUROBI = False

needs_gurobi = pytest.mark.skipif(not HAVE_GUROBI, reason="gurobipy not installed")


def knapsack(gp, grb, values, weights, capacity):
    """A model written the way a gurobipy user would write it."""
    n = len(values)
    m = gp.Model("knapsack")
    m.setParam("OutputFlag", 0)
    x = m.addVars(n, vtype=grb.BINARY, name="x")
    m.setObjective(gp.quicksum(values[j] * x[j] for j in range(n)), grb.MAXIMIZE)
    m.addConstr(gp.quicksum(weights[j] * x[j] for j in range(n)) <= capacity, name="cap")
    m.optimize()
    return m.Status, m.ObjVal, sorted(j for j in range(n) if x[j].X > 0.5)


def cover(gp, grb, sets, costs, universe):
    n = len(sets)
    m = gp.Model("cover")
    m.setParam("OutputFlag", 0)
    x = m.addVars(n, vtype=grb.BINARY, name="x")
    m.setObjective(gp.quicksum(costs[j] * x[j] for j in range(n)), grb.MINIMIZE)
    for element in range(universe):
        m.addConstr(gp.quicksum(x[j] for j in range(n) if element in sets[j]) >= 1)
    m.optimize()
    return m.Status, m.ObjVal, sorted(j for j in range(n) if x[j].X > 0.5)


KNAPSACKS = [
    ([12, 2, 7, 9, 4, 11, 3, 8, 6, 5], [4, 2, 3, 5, 1, 6, 2, 4, 3, 2], 12),
    ([5, 4, 3, 2, 1], [3, 3, 3, 3, 3], 7),
    ([10, 10, 10], [5, 5, 5], 11),
    ([1] * 12, [2] * 12, 13),
]


@pytest.mark.parametrize("values,weights,capacity", KNAPSACKS)
def test_knapsack_is_solved(values, weights, capacity):
    status, objective, chosen = knapsack(ripsolve, GRB, values, weights, capacity)
    assert status == GRB.OPTIMAL
    # The answer must be feasible and score what was reported.
    assert sum(weights[j] for j in chosen) <= capacity
    assert objective == pytest.approx(sum(values[j] for j in chosen))


@needs_gurobi
@pytest.mark.parametrize("values,weights,capacity", KNAPSACKS)
def test_knapsack_matches_gurobi(values, weights, capacity):
    mine = knapsack(ripsolve, GRB, values, weights, capacity)
    theirs = knapsack(gurobipy, GGRB, values, weights, capacity)
    assert mine[0] == theirs[0]
    assert mine[1] == pytest.approx(theirs[1])


@needs_gurobi
def test_set_cover_matches_gurobi():
    sets = [{0, 1}, {1, 2}, {2, 3}, {3, 4}, {0, 4}, {0, 1, 2, 3, 4}]
    costs = [1, 1, 1, 1, 1, 4]
    mine = cover(ripsolve, GRB, sets, costs, 5)
    theirs = cover(gurobipy, GGRB, sets, costs, 5)
    assert mine[0] == theirs[0]
    assert mine[1] == pytest.approx(theirs[1])


def test_an_integer_objective_is_exact():
    # Reported values must compare equal to the integer they are, not to within
    # floating point noise: users write `m.ObjVal == 31`.
    status, objective, _ = knapsack(ripsolve, GRB, *KNAPSACKS[0])
    assert status == GRB.OPTIMAL
    assert objective == 31.0


def test_expression_arithmetic():
    m = ripsolve.Model()
    x = m.addVar(vtype=GRB.BINARY, name="x")
    y = m.addVar(vtype=GRB.BINARY, name="y")
    # Every form a gurobipy user might reasonably write.
    m.setObjective(2 * x + 3 * y - 1, GRB.MINIMIZE)
    m.addConstr(x + y >= 1)
    m.addConstr(-x + 2 * y <= 1)
    m.addConstr(x - y == 0)
    m.optimize()
    assert m.Status == GRB.OPTIMAL
    # x == y and x + y >= 1 forces both to 1, costing 2 + 3 - 1.
    assert m.ObjVal == pytest.approx(4.0)
    assert x.X == pytest.approx(1.0)
    assert y.X == pytest.approx(1.0)


def test_repeated_terms_accumulate():
    m = ripsolve.Model()
    x = m.addVar(vtype=GRB.BINARY)
    m.setObjective(x + x, GRB.MINIMIZE)
    m.addConstr(x + x >= 2)
    m.optimize()
    assert m.ObjVal == pytest.approx(2.0)


def test_infeasible_is_reported():
    m = ripsolve.Model()
    x = m.addVar(vtype=GRB.BINARY)
    m.addConstr(x >= 1)
    m.addConstr(x <= 0)
    m.optimize()
    assert m.Status == GRB.INFEASIBLE
    assert m.SolCount == 0
    with pytest.raises(ValueError):
        _ = m.ObjVal


def test_continuous_is_the_default_vtype():
    # As in gurobipy. Defaulting to binary would be more convenient for this
    # solver's history and would silently change the meaning of a ported script.
    m = ripsolve.Model()
    v = m.addVar()
    assert v.VType == "C"


def test_unknown_parameter_is_rejected():
    m = ripsolve.Model()
    with pytest.raises(KeyError):
        m.setParam("NoSuchParam", 1)


def test_attributes_before_optimize_raise():
    m = ripsolve.Model()
    m.addVar(vtype=GRB.BINARY)
    with pytest.raises(ValueError):
        _ = m.Status


def test_model_shape_attributes():
    m = ripsolve.Model("shape")
    a, b = m.addVar(vtype=GRB.BINARY, name="a"), m.addVar(vtype=GRB.BINARY, name="b")
    m.addConstr(a + b <= 1, name="together")
    assert m.ModelName == "shape"
    assert m.NumVars == 2
    assert m.NumConstrs == 1
    assert a.VarName == "a"
    assert [v.VarName for v in m.getVars()] == ["a", "b"]


def test_read_lp_file():
    here = os.path.dirname(__file__)
    path = os.path.join(here, "..", "..", "samples", "bin_10var_5con.lp")
    m = ripsolve.read(path)
    assert m.NumVars == 10
    m.optimize()
    assert m.Status == GRB.OPTIMAL
    assert m.ObjVal == pytest.approx(25.0)


def test_write_round_trips(tmp_path):
    m = ripsolve.Model("rt")
    x = m.addVars(3, vtype=GRB.BINARY, name="x")
    m.setObjective(ripsolve.quicksum(x[j] for j in range(3)), GRB.MINIMIZE)
    m.addConstr(ripsolve.quicksum(x[j] for j in range(3)) >= 2)
    path = str(tmp_path / "model.lp")
    m.write(path)

    again = ripsolve.read(path)
    again.optimize()
    m.optimize()
    assert again.ObjVal == pytest.approx(m.ObjVal)


def test_time_limit_is_honoured():
    here = os.path.dirname(__file__)
    path = os.path.join(here, "..", "..", "samples", "v064c064.lp")
    m = ripsolve.read(path)
    m.setParam("TimeLimit", 0.05)
    m.setParam("Threads", 1)
    m.optimize()
    assert m.Status in (GRB.OPTIMAL, GRB.TIME_LIMIT)
    assert m.Runtime < 5.0


def test_bounds_can_fix_a_variable():
    m = ripsolve.Model()
    x, y = m.addVar(vtype=GRB.BINARY, name="x"), m.addVar(vtype=GRB.BINARY, name="y")
    m.setObjective(x + 5 * y, GRB.MINIMIZE)
    m.addConstr(x + y >= 1)
    x.UB = 0.0  # forces the expensive column in
    m.optimize()
    assert m.ObjVal == pytest.approx(5.0)


# --- mixed-integer models -------------------------------------------------


def mixed(gp, grb):
    """Binary, general-integer and continuous columns in the same model."""
    m = gp.Model("mixed")
    m.setParam("OutputFlag", 0)
    b = m.addVar(vtype=grb.BINARY, name="b")
    i = m.addVar(vtype=grb.INTEGER, ub=6, name="i")
    c = m.addVar(vtype=grb.CONTINUOUS, ub=4.0, name="c")
    m.setObjective(4 * b + 3 * i + 1.5 * c, grb.MINIMIZE)
    m.addConstr(2 * b + 3 * i + c >= 7.5)
    m.addConstr(i + c >= 3)
    m.optimize()
    return m.Status, m.ObjVal, (b.X, i.X, c.X)


def test_mixed_integer_is_solved():
    status, objective, (b, i, c) = mixed(ripsolve, GRB)
    assert status == GRB.OPTIMAL
    # The integer columns must land on integers; the continuous one need not.
    assert b == pytest.approx(round(b))
    assert i == pytest.approx(round(i))
    assert 2 * b + 3 * i + c >= 7.5 - 1e-6
    assert objective == pytest.approx(4 * b + 3 * i + 1.5 * c)


@needs_gurobi
def test_mixed_integer_matches_gurobi():
    mine = mixed(ripsolve, GRB)
    theirs = mixed(gurobipy, GGRB)
    assert mine[0] == theirs[0]
    assert mine[1] == pytest.approx(theirs[1])


def test_general_integer_takes_a_value_above_one():
    # Branching splits a range rather than fixing to 0/1, so 4 must be reachable.
    m = ripsolve.Model()
    x = m.addVar(vtype=GRB.INTEGER, ub=9, name="x")
    m.setObjective(-x, GRB.MINIMIZE)
    m.addConstr(3 * x <= 13)
    m.optimize()
    assert x.X == pytest.approx(4.0)


def test_continuous_variables_stay_fractional():
    # The distinguishing MIP case: an optimum that is fractional in a continuous
    # column is a valid answer, not something to branch away.
    m = ripsolve.Model()
    a = m.addVar(vtype=GRB.BINARY, name="a")
    b = m.addVar(vtype=GRB.CONTINUOUS, ub=5.0, name="b")
    m.setObjective(2 * a + 3 * b, GRB.MINIMIZE)
    m.addConstr(a + b >= 1.5)
    m.optimize()
    assert m.ObjVal == pytest.approx(3.5)
    assert b.X == pytest.approx(0.5)


def test_vtype_is_reported_back():
    m = ripsolve.Model()
    assert m.addVar(vtype=GRB.BINARY).VType == "B"
    assert m.addVar(vtype=GRB.INTEGER, ub=5).VType == "I"
    assert m.addVar(vtype=GRB.CONTINUOUS).VType == "C"


def test_unsupported_vtype_is_rejected():
    m = ripsolve.Model()
    with pytest.raises(ValueError):
        m.addVar(vtype="S")  # semi-continuous


def test_fractional_integer_bound_is_rejected():
    # Branching splits at an integer, so a bound of 2.5 is unreachable from either
    # side; better to refuse than to solve a subtly different model.
    m = ripsolve.Model()
    with pytest.raises(ValueError):
        m.addVar(vtype=GRB.INTEGER, ub=2.5)
    # The same bound on a continuous column is fine.
    m.addVar(vtype=GRB.CONTINUOUS, ub=2.5)


def test_mixed_model_round_trips_through_a_file(tmp_path):
    m = ripsolve.Model("rt")
    a = m.addVar(vtype=GRB.INTEGER, ub=5, name="a")
    b = m.addVar(vtype=GRB.CONTINUOUS, ub=3.0, name="b")
    m.setObjective(2 * a + b, GRB.MINIMIZE)
    m.addConstr(a + b >= 2.5)
    path = str(tmp_path / "mixed.lp")
    m.write(path)

    again = ripsolve.read(path)
    # Types must survive the trip; a Bounds section must not silently relax `a`.
    assert [v.VType for v in again.getVars()] == ["I", "C"]
    again.optimize()
    m.optimize()
    assert again.ObjVal == pytest.approx(m.ObjVal)
