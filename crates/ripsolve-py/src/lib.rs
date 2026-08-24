//! A Python interface shaped like `gurobipy`, for binary integer programs.
//!
//! The goal is that a `gurobipy` script that only uses binary variables runs
//! unchanged after swapping the import. So the names, the attribute spellings, and
//! the operator overloading all follow Gurobi's rather than anything more Rust-like:
//!
//! ```python
//! import ripsolve as gp
//! from ripsolve import GRB
//!
//! m = gp.Model("knapsack")
//! x = m.addVars(n, vtype=GRB.BINARY, name="x")
//! m.setObjective(gp.quicksum(value[j] * x[j] for j in range(n)), GRB.MAXIMIZE)
//! m.addConstr(gp.quicksum(weight[j] * x[j] for j in range(n)) <= capacity)
//! m.optimize()
//! print(m.ObjVal, [x[j].X for j in range(n)])
//! ```
//!
//! Binary, general-integer and continuous variables are all supported, so the
//! models this accepts are ordinary MIPs rather than pure BIPs. Quadratic terms,
//! SOS constraints, callbacks and lazy constraints are not.

use std::collections::BTreeMap;
use std::sync::{Arc, Mutex};

use pyo3::exceptions::{PyIndexError, PyKeyError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyTuple};

use ripsolve::model::{Problem, RowSense, Sense, VarType};
use ripsolve::search::{self, Options, Status as SearchStatus};
use ripsolve::sparse::SparseMatrix;

/// A row as the user built it, before translation into a [`Problem`].
struct Row {
    terms: Vec<(usize, f64)>,
    sense: RowSense,
    rhs: f64,
    name: String,
}

/// Everything a model holds, shared between the `Model` and the `Var`s that refer
/// back into it for their solution values.
struct ModelData {
    name: String,
    var_names: Vec<String>,
    var_types: Vec<VarType>,
    objective: Vec<f64>,
    objective_constant: f64,
    sense: Sense,
    rows: Vec<Row>,
    /// Bounds, so a caller can fix a variable the way Gurobi's `.LB`/`.UB` do.
    lower: Vec<f64>,
    upper: Vec<f64>,
    solution: Option<Solved>,
    params: Params,
}

#[derive(Clone)]
struct Solved {
    status: i32,
    objective: Option<f64>,
    bound: f64,
    values: Vec<f64>,
    nodes: usize,
    runtime: f64,
    gap: f64,
}

#[derive(Clone)]
struct Params {
    time_limit: Option<f64>,
    threads: usize,
    mip_gap: f64,
    output_flag: i32,
}

impl Default for Params {
    fn default() -> Self {
        Self {
            time_limit: None,
            threads: 0,
            mip_gap: 0.0,
            output_flag: 1,
        }
    }
}

impl ModelData {
    fn new(name: String) -> Self {
        Self {
            name,
            var_names: Vec::new(),
            var_types: Vec::new(),
            objective: Vec::new(),
            objective_constant: 0.0,
            // Gurobi minimizes unless told otherwise.
            sense: Sense::Minimize,
            rows: Vec::new(),
            lower: Vec::new(),
            upper: Vec::new(),
            solution: None,
            params: Params::default(),
        }
    }
}

type Shared = Arc<Mutex<ModelData>>;

/// Gurobi's status codes, so `m.Status == GRB.OPTIMAL` compares correctly.
const STATUS_OPTIMAL: i32 = 2;
const STATUS_INFEASIBLE: i32 = 3;
const STATUS_NODE_LIMIT: i32 = 7;
const STATUS_TIME_LIMIT: i32 = 9;

/// The `GRB` namespace of constants.
#[pyclass(name = "GRB", module = "ripsolve", from_py_object)]
#[derive(Clone)]
struct Grb;

#[pymethods]
impl Grb {
    #[classattr]
    const BINARY: char = 'B';
    #[classattr]
    const CONTINUOUS: char = 'C';
    #[classattr]
    const INTEGER: char = 'I';
    #[classattr]
    const MINIMIZE: i32 = 1;
    #[classattr]
    const MAXIMIZE: i32 = -1;
    #[classattr]
    const OPTIMAL: i32 = STATUS_OPTIMAL;
    #[classattr]
    const INFEASIBLE: i32 = STATUS_INFEASIBLE;
    #[classattr]
    const NODE_LIMIT: i32 = STATUS_NODE_LIMIT;
    #[classattr]
    const TIME_LIMIT: i32 = STATUS_TIME_LIMIT;
    #[classattr]
    const INFINITY: f64 = f64::INFINITY;
}

/// A linear expression: a weighted sum of variables plus a constant.
///
/// Terms are kept in a `BTreeMap` keyed by column so that repeated mentions of a
/// variable accumulate — `x + x` is `2 x`, as in Gurobi — and so the resulting row
/// is deterministic regardless of how the expression was assembled.
#[pyclass(name = "LinExpr", module = "ripsolve", from_py_object)]
#[derive(Clone, Default)]
struct LinExpr {
    terms: BTreeMap<usize, f64>,
    constant: f64,
}

impl LinExpr {
    fn of_var(index: usize) -> Self {
        let mut terms = BTreeMap::new();
        terms.insert(index, 1.0);
        Self {
            terms,
            constant: 0.0,
        }
    }

    fn scaled(&self, factor: f64) -> Self {
        Self {
            terms: self.terms.iter().map(|(&j, &v)| (j, v * factor)).collect(),
            constant: self.constant * factor,
        }
    }

    fn add(&self, other: &LinExpr) -> Self {
        let mut out = self.clone();
        for (&j, &v) in &other.terms {
            *out.terms.entry(j).or_insert(0.0) += v;
        }
        out.constant += other.constant;
        out
    }

    /// Interpret a Python object as an expression: a `Var`, a `LinExpr`, or a
    /// number.
    fn extract(value: &Bound<'_, PyAny>) -> PyResult<LinExpr> {
        if let Ok(var) = value.extract::<PyRef<Var>>() {
            return Ok(LinExpr::of_var(var.index));
        }
        if let Ok(expr) = value.extract::<PyRef<LinExpr>>() {
            return Ok(expr.clone());
        }
        if let Ok(number) = value.extract::<f64>() {
            return Ok(LinExpr {
                terms: BTreeMap::new(),
                constant: number,
            });
        }
        Err(PyTypeError::new_err(format!(
            "expected a Var, LinExpr or number, got {}",
            value.get_type().name()?
        )))
    }

    /// Build the constraint `self <sense> other`, moving everything to the left.
    fn compare(&self, other: &Bound<'_, PyAny>, sense: RowSense) -> PyResult<TempConstr> {
        let rhs = LinExpr::extract(other)?;
        let combined = self.add(&rhs.scaled(-1.0));
        Ok(TempConstr {
            terms: combined
                .terms
                .into_iter()
                .filter(|&(_, v)| v != 0.0)
                .collect(),
            sense,
            rhs: -combined.constant,
        })
    }
}

#[pymethods]
impl LinExpr {
    #[new]
    #[pyo3(signature = (value = None))]
    fn new(value: Option<&Bound<'_, PyAny>>) -> PyResult<Self> {
        match value {
            Some(value) => LinExpr::extract(value),
            None => Ok(LinExpr::default()),
        }
    }

    /// The expression's value at the current solution.
    #[getter(getValue)]
    fn get_value_attr(&self) -> f64 {
        self.constant
    }

    fn __add__(&self, other: &Bound<'_, PyAny>) -> PyResult<LinExpr> {
        Ok(self.add(&LinExpr::extract(other)?))
    }

    fn __radd__(&self, other: &Bound<'_, PyAny>) -> PyResult<LinExpr> {
        Ok(self.add(&LinExpr::extract(other)?))
    }

    fn __sub__(&self, other: &Bound<'_, PyAny>) -> PyResult<LinExpr> {
        Ok(self.add(&LinExpr::extract(other)?.scaled(-1.0)))
    }

    fn __rsub__(&self, other: &Bound<'_, PyAny>) -> PyResult<LinExpr> {
        Ok(LinExpr::extract(other)?.add(&self.scaled(-1.0)))
    }

    fn __mul__(&self, factor: f64) -> LinExpr {
        self.scaled(factor)
    }

    fn __rmul__(&self, factor: f64) -> LinExpr {
        self.scaled(factor)
    }

    fn __neg__(&self) -> LinExpr {
        self.scaled(-1.0)
    }

    fn __le__(&self, other: &Bound<'_, PyAny>) -> PyResult<TempConstr> {
        self.compare(other, RowSense::Le)
    }

    fn __ge__(&self, other: &Bound<'_, PyAny>) -> PyResult<TempConstr> {
        self.compare(other, RowSense::Ge)
    }

    fn __eq__(&self, other: &Bound<'_, PyAny>) -> PyResult<TempConstr> {
        self.compare(other, RowSense::Eq)
    }

    fn __repr__(&self) -> String {
        let mut out = String::new();
        for (j, v) in &self.terms {
            out.push_str(&format!(
                "{}{v} x{j} ",
                if out.is_empty() { "" } else { "+ " }
            ));
        }
        if self.constant != 0.0 || out.is_empty() {
            out.push_str(&format!("+ {}", self.constant));
        }
        format!("<LinExpr: {}>", out.trim())
    }
}

/// A comparison that has not yet been added to a model, as produced by `x + y <= 1`.
#[pyclass(name = "TempConstr", module = "ripsolve", from_py_object)]
#[derive(Clone)]
struct TempConstr {
    terms: Vec<(usize, f64)>,
    sense: RowSense,
    rhs: f64,
}

/// A decision variable.
#[pyclass(name = "Var", module = "ripsolve")]
struct Var {
    data: Shared,
    index: usize,
}

impl Var {
    fn expr(&self) -> LinExpr {
        LinExpr::of_var(self.index)
    }
}

#[pymethods]
impl Var {
    /// The variable's value in the current solution, as Gurobi's `.X`.
    #[getter(X)]
    fn x(&self) -> PyResult<f64> {
        let data = self.data.lock().expect("model lock");
        match &data.solution {
            Some(solved) if solved.objective.is_some() => Ok(solved.values[self.index]),
            _ => Err(PyValueError::new_err(
                "no solution available; call optimize() first and check Status",
            )),
        }
    }

    #[getter(VarName)]
    fn var_name(&self) -> String {
        self.data.lock().expect("model lock").var_names[self.index].clone()
    }

    #[getter(VType)]
    fn vtype(&self) -> char {
        let data = self.data.lock().expect("model lock");
        match data.var_types[self.index] {
            VarType::Continuous => 'C',
            VarType::Integer if data.upper[self.index] <= 1.0 && data.lower[self.index] >= 0.0 => {
                'B'
            }
            VarType::Integer => 'I',
        }
    }

    #[getter(Obj)]
    fn obj(&self) -> f64 {
        self.data.lock().expect("model lock").objective[self.index]
    }

    #[getter(LB)]
    fn lb(&self) -> f64 {
        self.data.lock().expect("model lock").lower[self.index]
    }

    #[setter(LB)]
    fn set_lb(&self, value: f64) {
        self.data.lock().expect("model lock").lower[self.index] = value;
    }

    #[getter(UB)]
    fn ub(&self) -> f64 {
        self.data.lock().expect("model lock").upper[self.index]
    }

    #[setter(UB)]
    fn set_ub(&self, value: f64) {
        self.data.lock().expect("model lock").upper[self.index] = value;
    }

    fn __add__(&self, other: &Bound<'_, PyAny>) -> PyResult<LinExpr> {
        Ok(self.expr().add(&LinExpr::extract(other)?))
    }

    fn __radd__(&self, other: &Bound<'_, PyAny>) -> PyResult<LinExpr> {
        Ok(self.expr().add(&LinExpr::extract(other)?))
    }

    fn __sub__(&self, other: &Bound<'_, PyAny>) -> PyResult<LinExpr> {
        Ok(self.expr().add(&LinExpr::extract(other)?.scaled(-1.0)))
    }

    fn __rsub__(&self, other: &Bound<'_, PyAny>) -> PyResult<LinExpr> {
        Ok(LinExpr::extract(other)?.add(&self.expr().scaled(-1.0)))
    }

    fn __mul__(&self, factor: f64) -> LinExpr {
        self.expr().scaled(factor)
    }

    fn __rmul__(&self, factor: f64) -> LinExpr {
        self.expr().scaled(factor)
    }

    fn __neg__(&self) -> LinExpr {
        self.expr().scaled(-1.0)
    }

    fn __le__(&self, other: &Bound<'_, PyAny>) -> PyResult<TempConstr> {
        self.expr().compare(other, RowSense::Le)
    }

    fn __ge__(&self, other: &Bound<'_, PyAny>) -> PyResult<TempConstr> {
        self.expr().compare(other, RowSense::Ge)
    }

    fn __eq__(&self, other: &Bound<'_, PyAny>) -> PyResult<TempConstr> {
        self.expr().compare(other, RowSense::Eq)
    }

    fn __hash__(&self) -> u64 {
        self.index as u64
    }

    fn __repr__(&self) -> String {
        format!("<Var {}>", self.var_name())
    }
}

/// A constraint that has been added to a model.
#[pyclass(name = "Constr", module = "ripsolve")]
struct Constr {
    data: Shared,
    index: usize,
}

#[pymethods]
impl Constr {
    #[getter(ConstrName)]
    fn constr_name(&self) -> String {
        self.data.lock().expect("model lock").rows[self.index]
            .name
            .clone()
    }

    fn __repr__(&self) -> String {
        format!("<Constr {}>", self.constr_name())
    }
}

/// An optimization model.
#[pyclass(name = "Model", module = "ripsolve")]
struct Model {
    data: Shared,
}

impl Model {
    /// Translate the model as built into the solver's own representation.
    fn to_problem(data: &ModelData) -> PyResult<Problem> {
        let n = data.var_names.len();
        if n == 0 {
            return Err(PyValueError::new_err("model has no variables"));
        }
        let m = data.rows.len();
        let mut triplets = Vec::new();
        let (mut row_lb, mut row_ub, mut row_names) = (Vec::new(), Vec::new(), Vec::new());
        for (i, row) in data.rows.iter().enumerate() {
            for &(j, v) in &row.terms {
                triplets.push((i, j, v));
            }
            let (lo, hi) = row.sense.bounds(row.rhs);
            row_lb.push(lo);
            row_ub.push(hi);
            row_names.push(row.name.clone());
        }

        // The solver always minimizes, so a maximization is negated on the way in
        // and `Problem::objective_value` converts results back.
        let negate = data.sense == Sense::Maximize;
        let flip = |v: f64| if negate { -v } else { v };

        let problem = Problem {
            name: data.name.clone(),
            sense: data.sense,
            obj: data.objective.iter().map(|&c| flip(c)).collect(),
            obj_offset: flip(data.objective_constant),
            matrix: SparseMatrix::from_triplets(m, n, triplets),
            row_lb,
            row_ub,
            col_lb: data.lower.clone(),
            col_ub: data.upper.clone(),
            col_type: data.var_types.clone(),
            col_names: data.var_names.clone(),
            row_names,
        };
        problem
            .validate()
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(problem)
    }

    fn solved(&self) -> PyResult<Solved> {
        self.data
            .lock()
            .expect("model lock")
            .solution
            .clone()
            .ok_or_else(|| PyValueError::new_err("call optimize() first"))
    }
}

// gurobipy's method names, deliberately. Renaming them to Rust convention would
// defeat the entire point of the binding, which is that an existing script runs
// unchanged.
#[allow(non_snake_case)]
#[pymethods]
impl Model {
    #[new]
    #[pyo3(signature = (name = ""))]
    fn new(name: &str) -> Self {
        Self {
            data: Arc::new(Mutex::new(ModelData::new(name.to_string()))),
        }
    }

    /// Add one variable.
    ///
    /// `vtype` follows Gurobi: `GRB.BINARY`, `GRB.INTEGER` or `GRB.CONTINUOUS`,
    /// defaulting to continuous as Gurobi's does. Defaulting to binary instead
    /// would be more convenient here and is exactly the wrong trade: a ported
    /// script calling `addVar()` would silently get a different model.
    ///
    /// Binary is an integer pinned to `[0, 1]`, so the default bounds depend on the
    /// type -- as in Gurobi, where an unbounded integer is `[0, inf)`.
    #[pyo3(signature = (obj = 0.0, vtype = 'C', name = "", lb = None, ub = None))]
    fn addVar(
        &self,
        obj: f64,
        vtype: char,
        name: &str,
        lb: Option<f64>,
        ub: Option<f64>,
    ) -> PyResult<Var> {
        let (kind, default_lb, default_ub) = match vtype {
            'B' => (VarType::Integer, 0.0, 1.0),
            'I' => (VarType::Integer, 0.0, f64::INFINITY),
            'C' => (VarType::Continuous, 0.0, f64::INFINITY),
            other => {
                return Err(PyValueError::new_err(format!(
                    "unsupported vtype {other:?}; expected GRB.BINARY, GRB.INTEGER or GRB.CONTINUOUS"
                )));
            }
        };
        let (lb, ub) = (lb.unwrap_or(default_lb), ub.unwrap_or(default_ub));
        if kind == VarType::Integer {
            for bound in [lb, ub] {
                if bound.is_finite() && (bound - bound.round()).abs() > 1e-9 {
                    return Err(PyValueError::new_err(format!(
                        "integer variable has the fractional bound {bound}"
                    )));
                }
            }
        }
        let mut data = self.data.lock().expect("model lock");
        let index = data.var_names.len();
        let name = if name.is_empty() {
            format!("C{index}")
        } else {
            name.to_string()
        };
        data.var_names.push(name);
        data.var_types.push(kind);
        data.objective.push(obj);
        data.lower.push(lb);
        data.upper.push(ub);
        data.solution = None;
        Ok(Var {
            data: Arc::clone(&self.data),
            index,
        })
    }

    /// Add many binary variables, returning a dict keyed as Gurobi's `tupledict` is.
    ///
    /// An integer count gives integer keys `0..n`; an iterable of keys gives those
    /// keys, so `addVars(range(n))` and `addVars(pairs)` both behave as expected.
    // Mirrors gurobipy's signature; splitting it into a struct would make the
    // binding read nothing like the API it exists to imitate.
    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature = (keys, obj = 0.0, vtype = 'C', name = "", lb = None, ub = None))]
    fn addVars<'py>(
        &self,
        py: Python<'py>,
        keys: &Bound<'py, PyAny>,
        obj: f64,
        vtype: char,
        name: &str,
        lb: Option<f64>,
        ub: Option<f64>,
    ) -> PyResult<Bound<'py, PyDict>> {
        let key_list: Vec<Py<PyAny>> = if let Ok(count) = keys.extract::<usize>() {
            (0..count)
                .map(|i| i.into_pyobject(py).unwrap().unbind().into())
                .collect()
        } else {
            keys.try_iter()?
                .map(|k| k.map(|k| k.unbind()))
                .collect::<PyResult<_>>()?
        };

        let out = PyDict::new(py);
        for key in key_list {
            let bound = key.bind(py);
            let label = if name.is_empty() {
                String::new()
            } else {
                format!("{name}[{}]", bound.str()?)
            };
            let var = self.addVar(obj, vtype, &label, lb, ub)?;
            out.set_item(bound, Py::new(py, var)?)?;
        }
        Ok(out)
    }

    /// Add a constraint built by comparing expressions, as in `x + y <= 1`.
    #[pyo3(signature = (constraint, name = ""))]
    fn addConstr(&self, constraint: PyRef<TempConstr>, name: &str) -> PyResult<Constr> {
        let mut data = self.data.lock().expect("model lock");
        let index = data.rows.len();
        let name = if name.is_empty() {
            format!("R{index}")
        } else {
            name.to_string()
        };
        data.rows.push(Row {
            terms: constraint.terms.clone(),
            sense: constraint.sense,
            rhs: constraint.rhs,
            name,
        });
        data.solution = None;
        Ok(Constr {
            data: Arc::clone(&self.data),
            index,
        })
    }

    /// Add several constraints from an iterable of comparisons.
    #[pyo3(signature = (constraints, name = ""))]
    fn addConstrs(&self, constraints: &Bound<'_, PyAny>, name: &str) -> PyResult<usize> {
        let mut added = 0;
        for item in constraints.try_iter()? {
            let item = item?;
            let constraint = item.extract::<PyRef<TempConstr>>()?;
            self.addConstr(constraint, name)?;
            added += 1;
        }
        Ok(added)
    }

    /// Set the objective, and optionally the direction.
    #[pyo3(signature = (expression, sense = None))]
    fn setObjective(&self, expression: &Bound<'_, PyAny>, sense: Option<i32>) -> PyResult<()> {
        let expr = LinExpr::extract(expression)?;
        let mut data = self.data.lock().expect("model lock");
        let n = data.var_names.len();
        data.objective = vec![0.0; n];
        for (&j, &v) in &expr.terms {
            if j >= n {
                return Err(PyIndexError::new_err(
                    "objective refers to an unknown variable",
                ));
            }
            data.objective[j] = v;
        }
        data.objective_constant = expr.constant;
        if let Some(sense) = sense {
            data.sense = if sense < 0 {
                Sense::Maximize
            } else {
                Sense::Minimize
            };
        }
        data.solution = None;
        Ok(())
    }

    /// Set a solver parameter. Unknown names are rejected rather than ignored, so a
    /// misspelling does not silently do nothing.
    fn setParam(&self, name: &str, value: &Bound<'_, PyAny>) -> PyResult<()> {
        let mut data = self.data.lock().expect("model lock");
        match name.to_ascii_lowercase().as_str() {
            "timelimit" => data.params.time_limit = Some(value.extract()?),
            "threads" => data.params.threads = value.extract()?,
            "mipgap" => data.params.mip_gap = value.extract()?,
            "outputflag" | "logtoconsole" => data.params.output_flag = value.extract()?,
            other => {
                return Err(PyKeyError::new_err(format!("unknown parameter {other:?}")));
            }
        }
        Ok(())
    }

    /// Solve the model.
    fn optimize(&self, py: Python<'_>) -> PyResult<()> {
        let (problem, params) = {
            let data = self.data.lock().expect("model lock");
            (Model::to_problem(&data)?, data.params.clone())
        };

        let options = Options {
            time_limit: params.time_limit.map(std::time::Duration::from_secs_f64),
            gap_tolerance: params.mip_gap,
            threads: if params.threads == 0 {
                std::thread::available_parallelism().map_or(1, |n| n.get())
            } else {
                params.threads
            },
            ..Options::default()
        };

        // Release the GIL: the solve is pure Rust and can take minutes, and holding
        // the GIL would freeze every other Python thread for its duration.
        let started = std::time::Instant::now();
        let solution = py.detach(|| search::solve(&problem, options));
        let runtime = started.elapsed().as_secs_f64();

        let status = match solution.status {
            SearchStatus::Optimal => STATUS_OPTIMAL,
            SearchStatus::Infeasible => STATUS_INFEASIBLE,
            SearchStatus::NodeLimit => STATUS_NODE_LIMIT,
            SearchStatus::TimeLimit => STATUS_TIME_LIMIT,
        };

        let mut data = self.data.lock().expect("model lock");
        data.solution = Some(Solved {
            status,
            objective: solution.objective,
            bound: solution.bound,
            values: solution.x.clone(),
            nodes: solution.nodes,
            runtime,
            gap: solution.gap(),
        });
        Ok(())
    }

    /// Every variable, in the order they were added.
    fn getVars(&self, py: Python<'_>) -> PyResult<Vec<Py<Var>>> {
        let count = self.data.lock().expect("model lock").var_names.len();
        (0..count)
            .map(|index| {
                Py::new(
                    py,
                    Var {
                        data: Arc::clone(&self.data),
                        index,
                    },
                )
            })
            .collect()
    }

    /// Write the model to an LP file.
    fn write(&self, path: &str) -> PyResult<()> {
        let data = self.data.lock().expect("model lock");
        let problem = Model::to_problem(&data)?;
        std::fs::write(path, problem_to_lp(&problem))
            .map_err(|e| PyValueError::new_err(format!("writing {path}: {e}")))
    }

    #[getter(ObjVal)]
    fn obj_val(&self) -> PyResult<f64> {
        self.solved()?
            .objective
            .ok_or_else(|| PyValueError::new_err("no solution was found"))
    }

    #[getter(ObjBound)]
    fn obj_bound(&self) -> PyResult<f64> {
        Ok(self.solved()?.bound)
    }

    #[getter(Status)]
    fn status(&self) -> PyResult<i32> {
        Ok(self.solved()?.status)
    }

    #[getter(SolCount)]
    fn sol_count(&self) -> PyResult<usize> {
        Ok(usize::from(self.solved()?.objective.is_some()))
    }

    #[getter(NodeCount)]
    fn node_count(&self) -> PyResult<f64> {
        Ok(self.solved()?.nodes as f64)
    }

    #[getter(Runtime)]
    fn runtime(&self) -> PyResult<f64> {
        Ok(self.solved()?.runtime)
    }

    #[getter(MIPGap)]
    fn mip_gap(&self) -> PyResult<f64> {
        Ok(self.solved()?.gap)
    }

    #[getter(NumVars)]
    fn num_vars(&self) -> usize {
        self.data.lock().expect("model lock").var_names.len()
    }

    #[getter(NumConstrs)]
    fn num_constrs(&self) -> usize {
        self.data.lock().expect("model lock").rows.len()
    }

    #[getter(ModelName)]
    fn model_name(&self) -> String {
        self.data.lock().expect("model lock").name.clone()
    }

    #[getter(ModelSense)]
    fn model_sense(&self) -> i32 {
        match self.data.lock().expect("model lock").sense {
            Sense::Minimize => 1,
            Sense::Maximize => -1,
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "<ripsolve.Model {:?}: {} vars, {} constrs>",
            self.model_name(),
            self.num_vars(),
            self.num_constrs()
        )
    }
}

/// Render a problem back to LP format, for `Model.write`.
fn problem_to_lp(problem: &Problem) -> String {
    use std::fmt::Write as _;
    let mut out = String::new();
    let _ = writeln!(out, "\\ {}", problem.name);
    out.push_str(if problem.sense == Sense::Maximize {
        "Maximize\n obj:"
    } else {
        "Minimize\n obj:"
    });
    // The stored objective is in minimization form; undo that for the file.
    let flip = |v: f64| {
        if problem.sense == Sense::Maximize {
            -v
        } else {
            v
        }
    };
    for (j, &c) in problem.obj.iter().enumerate() {
        if c != 0.0 {
            let _ = write!(out, " + {} {}", flip(c), problem.col_names[j]);
        }
    }
    out.push_str("\nSubject To\n");
    let csr = problem.matrix.to_csr();
    for i in 0..problem.n_rows() {
        let (cols, vals) = csr.column(i);
        if cols.is_empty() {
            continue;
        }
        // A range row needs one line per finite side.
        for (bound, op) in [(problem.row_lb[i], ">="), (problem.row_ub[i], "<=")] {
            if !bound.is_finite() {
                continue;
            }
            let _ = write!(out, " {}: ", problem.row_names[i]);
            for (&j, &v) in cols.iter().zip(vals) {
                let _ = write!(out, "+ {v} {} ", problem.col_names[j]);
            }
            let _ = writeln!(out, "{op} {bound}");
            if problem.row_lb[i] == problem.row_ub[i] {
                break;
            }
        }
    }
    // Bounds first, then integrality: every column needs its range stated, and only
    // the integer ones get declared.
    out.push_str("Bounds\n");
    for (j, name) in problem.col_names.iter().enumerate() {
        let lo = problem.col_lb[j];
        let hi = problem.col_ub[j];
        match (lo.is_finite(), hi.is_finite()) {
            (true, true) => {
                let _ = writeln!(out, " {lo} <= {name} <= {hi}");
            }
            (true, false) => {
                let _ = writeln!(out, " {name} >= {lo}");
            }
            (false, true) => {
                let _ = writeln!(out, " -inf <= {name} <= {hi}");
            }
            (false, false) => {
                let _ = writeln!(out, " {name} free");
            }
        }
    }
    let integers: Vec<&String> = problem
        .col_names
        .iter()
        .enumerate()
        .filter(|&(j, _)| problem.is_integer(j))
        .map(|(_, n)| n)
        .collect();
    if !integers.is_empty() {
        out.push_str("General\n");
        for name in integers {
            let _ = write!(out, " {name}");
        }
        out.push('\n');
    }
    out.push_str("End\n");
    out
}

/// Sum an iterable of variables or expressions.
#[pyfunction]
fn quicksum(terms: &Bound<'_, PyAny>) -> PyResult<LinExpr> {
    let mut total = LinExpr::default();
    for item in terms.try_iter()? {
        total = total.add(&LinExpr::extract(&item?)?);
    }
    Ok(total)
}

/// Read a model from an LP or MPS file.
#[pyfunction]
fn read(py: Python<'_>, path: &str) -> PyResult<Py<Model>> {
    let problem = Problem::from_file(std::path::Path::new(path))
        .map_err(|e| PyValueError::new_err(format!("reading {path}: {e}")))?;

    let model = Model::new(&problem.name);
    {
        let mut data = model.data.lock().expect("model lock");
        data.sense = problem.sense;
        data.var_names = problem.col_names.clone();
        // `Problem` holds the objective in minimization form; the Python side keeps
        // it as the user would have written it.
        let flip = |v: f64| {
            if problem.sense == Sense::Maximize {
                -v
            } else {
                v
            }
        };
        data.objective = problem.obj.iter().map(|&c| flip(c)).collect();
        data.objective_constant = flip(problem.obj_offset);
        data.lower = problem.col_lb.clone();
        data.upper = problem.col_ub.clone();
        data.var_types = problem.col_type.clone();

        let csr = problem.matrix.to_csr();
        for i in 0..problem.n_rows() {
            let (cols, vals) = csr.column(i);
            let terms: Vec<(usize, f64)> = cols.iter().copied().zip(vals.iter().copied()).collect();
            let (lo, hi) = (problem.row_lb[i], problem.row_ub[i]);
            let (sense, rhs) = if lo == hi {
                (RowSense::Eq, lo)
            } else if hi.is_finite() {
                (RowSense::Le, hi)
            } else {
                (RowSense::Ge, lo)
            };
            data.rows.push(Row {
                terms,
                sense,
                rhs,
                name: problem.row_names[i].clone(),
            });
        }
    }
    Py::new(py, model)
}

// Named differently from the crate so the import does not shadow it; `#[pyo3(name)]`
// is what Python sees.
#[pymodule]
#[pyo3(name = "ripsolve")]
fn ripsolve_module(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<Model>()?;
    module.add_class::<Var>()?;
    module.add_class::<Constr>()?;
    module.add_class::<LinExpr>()?;
    module.add_class::<TempConstr>()?;
    module.add_class::<Grb>()?;
    module.add_function(wrap_pyfunction!(quicksum, module)?)?;
    module.add_function(wrap_pyfunction!(read, module)?)?;
    module.add("GRB", module.py().get_type::<Grb>())?;
    let _ = PyTuple::empty(module.py());
    Ok(())
}
