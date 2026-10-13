"""Callbacks: lazy constraints and cuts added while the search runs.

The model is clique partitioning, written the way a gurobipy user writes a lazy
constraint callback: a binary per pair of points, a degree row per point, and the
triangle rows, y[a,x] + y[x,c] - y[a,c] <= 1, left to the callback. The same function
builds and solves it under both modules, and both must reach the optimum of the model
with every triangle written out.
"""

import itertools
import os
import random
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import ripsolve
from ripsolve import GRB

try:
    import gurobipy

    HAVE_GUROBI = True
except ImportError:  # pragma: no cover - environment dependent
    HAVE_GUROBI = False

needs_gurobi = pytest.mark.skipif(not HAVE_GUROBI, reason="gurobipy not installed")


def points(n, seed):
    rng = random.Random(seed)
    return [(rng.randrange(100), rng.randrange(100)) for _ in range(n)]


def distance(p, q):
    return round(((p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2) ** 0.5)


def triangles(n):
    for x in range(n):
        for a, c in itertools.combinations([i for i in range(n) if i != x], 2):
            yield a, x, c


def y_of(y, i, j):
    return y[min(i, j), max(i, j)]


def partition(gp, grb, at, size, lazy=True):
    """Group the points into groups of `size` with the least total distance within
    groups. Returns the model's status, objective and the groups."""
    n = len(at)
    m = gp.Model("partition")
    m.setParam("OutputFlag", 0)
    pairs = list(itertools.combinations(range(n), 2))
    y = m.addVars(pairs, vtype=grb.BINARY, name="y")
    m.setObjective(gp.quicksum(distance(at[i], at[j]) * y[i, j] for i, j in pairs))
    for i in range(n):
        m.addConstr(gp.quicksum(y_of(y, i, j) for j in range(n) if j != i) == size - 1)
    if not lazy:
        for a, x, c in triangles(n):
            m.addConstr(y_of(y, a, x) + y_of(y, x, c) - y_of(y, a, c) <= 1)
        m.optimize()
    else:
        m._y = y
        m._n = n

        def callback(model, where):
            if where == grb.Callback.MIPSOL:
                value = model.cbGetSolution(model._y)
                for a, x, c in triangles(model._n):
                    if (
                        y_of(value, a, x) + y_of(value, x, c) - y_of(value, a, c)
                        > 1.5
                    ):
                        yy = model._y
                        model.cbLazy(y_of(yy, a, x) + y_of(yy, x, c) - y_of(yy, a, c) <= 1)

        m.setParam("LazyConstraints", 1)
        m.optimize(callback)
    groups = sorted(
        sorted({i} | {j for j in range(n) if j != i and y_of(y, i, j).X > 0.5})
        for i in range(n)
    )
    unique = [g for k, g in enumerate(groups) if k == 0 or g != groups[k - 1]]
    return m.Status, m.ObjVal, unique


def is_partition(groups, n, size):
    return sorted(i for g in groups for i in g) == list(range(n)) and all(
        len(g) == size for g in groups
    )


@pytest.mark.parametrize("seed", [1, 2, 3])
def test_lazy_triangles_reach_the_optimum_of_the_full_model(seed):
    at = points(12, seed)
    status, lazy, groups = partition(ripsolve, GRB, at, 4)
    _, full, _ = partition(ripsolve, GRB, at, 4, lazy=False)
    assert status == GRB.OPTIMAL
    assert lazy == pytest.approx(full)
    assert is_partition(groups, 12, 4)


@needs_gurobi
@pytest.mark.parametrize("seed", [1, 2])
def test_the_same_callback_gives_gurobis_answer(seed):
    at = points(12, seed)
    _, ours, groups = partition(ripsolve, GRB, at, 4)
    _, theirs, _ = partition(gurobipy, gurobipy.GRB, at, 4)
    assert ours == pytest.approx(theirs)
    assert is_partition(groups, 12, 4)


def small_model():
    m = ripsolve.Model()
    m.setParam("OutputFlag", 0)
    x = m.addVars(3, vtype=GRB.BINARY, name="x")
    m.setObjective(x[0] + 2 * x[1] + 3 * x[2], GRB.MAXIMIZE)
    m.addConstr(x[0] + x[1] + x[2] <= 2)
    return m, x


def test_lazy_constraints_need_the_parameter():
    m, x = small_model()

    def callback(model, where):
        if where == GRB.Callback.MIPSOL:
            model.cbLazy(x[2] <= 0)

    with pytest.raises(RuntimeError, match="LazyConstraints"):
        m.optimize(callback)


def test_a_lazy_row_changes_the_answer():
    m, x = small_model()
    m.Params.LazyConstraints = 1
    assert m.Params.LazyConstraints == 1

    def callback(model, where):
        if where == GRB.Callback.MIPSOL and model.cbGetSolution(x[2]) > 0.5:
            model.cbLazy(x[2], "<=", 0)

    m.optimize(callback)
    assert m.Status == GRB.OPTIMAL
    assert m.ObjVal == pytest.approx(3.0)
    assert [x[j].X for j in range(3)] == pytest.approx([1.0, 1.0, 0.0])


def test_cuts_belong_at_nodes():
    m, x = small_model()

    def callback(model, where):
        if where == GRB.Callback.MIPSOL:
            model.cbCut(x[0] <= 1)

    with pytest.raises(RuntimeError, match="MIPNODE"):
        m.optimize(callback)


def test_an_exception_in_the_callback_comes_back_from_optimize():
    m, _ = small_model()

    def callback(model, where):
        raise ZeroDivisionError("from the callback")

    with pytest.raises(ZeroDivisionError, match="from the callback"):
        m.optimize(callback)


def test_values_come_back_shaped_like_what_was_asked():
    m, x = small_model()
    m.addConstr(x[0] + x[1] + x[2] >= 0.5)
    seen = {}

    def callback(model, where):
        if where == GRB.Callback.MIPSOL:
            seen["dict"] = model.cbGetSolution(x)
            seen["list"] = model.cbGetSolution([x[0], x[2]])
            seen["one"] = model.cbGetSolution(x[1])
            seen["obj"] = model.cbGet(GRB.Callback.MIPSOL_OBJ)
        elif where == GRB.Callback.MIPNODE:
            assert model.cbGet(GRB.Callback.MIPNODE_STATUS) == GRB.OPTIMAL
            assert len(model.cbGetNodeRel([x[0], x[1]])) == 2

    m.optimize(callback)
    assert set(seen["dict"]) == {0, 1, 2}
    assert len(seen["list"]) == 2
    assert isinstance(seen["one"], float)
    assert seen["obj"] == pytest.approx(
        sum(c * seen["dict"][j] for j, c in enumerate([1, 2, 3]))
    )


def test_callback_methods_refuse_outside_a_callback():
    m, x = small_model()
    with pytest.raises(RuntimeError):
        m.cbGetSolution(x[0])
