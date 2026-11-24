//! An interior-point method for the LP relaxation: Mehrotra's predictor-corrector.
//!
//! The simplex methods walk vertices, and on some relaxations they walk for a very
//! long time: over MIPLIB's benchmark set, 47 roots of 240 finish with neither of them
//! in a minute, and on 18 of those HiGHS's interior-point code with crossover beats its
//! simplex by more than half or is the only one to finish (`dano3_3` 54s against 3,
//! `cryptanalysiskb128n5obj14` 40 against 3). The cost of an interior-point iteration
//! is one sparse Cholesky factorization, and the count of iterations barely grows with
//! the model, which is what makes it the method for the large, degenerate ones.
//!
//! # Form
//!
//! The same one the simplex uses: `n + m` variables `v = (x, s)`, the rows `A x - s = 0`,
//! and bounds `l <= v <= u`, any of them infinite. A finite lower bound gets a slack
//! `xl = v - l >= 0` with dual `zl >= 0`, a finite upper `xu = u - v >= 0` with `zu`;
//! fixed variables stay where they are and free ones take a primal regularization.
//! Every Newton step reduces to the normal equations `(A Θ_x A' + Θ_s + δ I) dy = r`,
//! where the logicals contribute only a diagonal.
//!
//! The method is infeasible-start: residuals in every equation shrink with the
//! duality measure rather than being kept at zero.
//!
//! What comes out is a point near the middle of the optimal face, not a vertex; the
//! search needs a basis, which [`crate::lp::Lp::solve_interior`] recovers from it.

use crate::lp::cholesky::{Abandoned, Cholesky, Stop};
use crate::sparse::SparseMatrix;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::Instant;

/// Relative primal and dual infeasibility, and relative duality gap, to stop at.
const TOLERANCE: f64 = 1e-8;
const MAX_ITERATIONS: usize = 200;
/// Infeasibility and gap at which a stalled solve is still handed to the crossover.
const APPROXIMATE: f64 = 1e-5;
const APPROXIMATE_GAP: f64 = 1e-3;
/// Bounds at or beyond this are infinite.
const INFINITE_BOUND: f64 = 1e20;
/// Regularizations: primal on the diagonal `D` of each variable, dual on the normal
/// matrix.
const PRIMAL_REG: f64 = 1e-10;
const DUAL_REG: f64 = 1e-10;
/// The fraction of the way to the boundary each step goes.
const STEP_FRACTION: f64 = 0.9995;
/// The factor may hold at most this many entries, 800MB of them.
const MAX_FACTOR_NONZEROS: f64 = 1e8;
/// And cost at most this much arithmetic to compute, `sum` over its columns of the
/// count squared: around ten seconds a factorization, and an interior-point solve is
/// twenty to eighty of them.
const MAX_FACTOR_FLOPS: f64 = 3e11;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum IpmStatus {
    Optimal,
    /// Stalled short of the tolerance but close enough for a crossover to finish from:
    /// the simplex method that ends the crossover repairs what is left.
    Approximate,
    /// Stopped without converging: out of iterations, time, or numerical room.
    NotConverged,
    /// The normal equations would factor too densely to be worth trying.
    TooDense,
}

#[derive(Clone, Debug)]
pub struct IpmResult {
    pub status: IpmStatus,
    /// Values of all `n + m` variables, structurals then logicals.
    pub v: Vec<f64>,
    /// Row duals.
    pub y: Vec<f64>,
    /// Reduced costs of all `n + m` variables, `zl - zu`.
    pub z: Vec<f64>,
    pub iterations: usize,
    pub objective: f64,
}

/// Which form of the normal equations the factor holds.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Mode {
    /// `A Θ_x A' + D`, one row and column per row of `A`.
    Rows,
    /// `Θ_x^-1 + A' D^-1 A`, one per column; its inverse gives the row form's by
    /// `M^-1 = D^-1 - D^-1 A N^-1 A' D^-1`.
    Columns,
}

/// What kind of bounds a variable has.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Kind {
    Fixed,
    Free,
    Lower,
    Upper,
    Boxed,
}

impl Kind {
    fn has_lower(self) -> bool {
        matches!(self, Kind::Lower | Kind::Boxed)
    }
    fn has_upper(self) -> bool {
        matches!(self, Kind::Upper | Kind::Boxed)
    }
}

struct Ipm<'a> {
    n: usize,
    m: usize,
    /// The scaled matrix, by column and by row.
    a: SparseMatrix,
    at: SparseMatrix,
    cost: Vec<f64>,
    lower: Vec<f64>,
    upper: Vec<f64>,
    kind: Vec<Kind>,
    /// Structural columns left out of the factor and handled by conjugate gradients.
    dense: Vec<bool>,
    dense_list: Vec<usize>,
    chol: Cholesky,
    mode: Mode,
    /// Rows left out of the column form's factor: none.
    row_skip: Vec<bool>,
    /// `D^-1` of the last factorization, for the column form.
    d_inv: Vec<f64>,
    /// `Θ_x^-1` of the last factorization, for the column form.
    col_diag: Vec<f64>,
    work: Vec<f64>,
    deadline: Option<Instant>,
    interrupt: Option<&'a Arc<AtomicBool>>,
}

/// Search directions.
#[derive(Clone, Default)]
struct Direction {
    v: Vec<f64>,
    y: Vec<f64>,
    xl: Vec<f64>,
    xu: Vec<f64>,
    zl: Vec<f64>,
    zu: Vec<f64>,
}

impl Direction {
    fn new(total: usize, m: usize) -> Self {
        Direction {
            v: vec![0.0; total],
            y: vec![0.0; m],
            xl: vec![0.0; total],
            xu: vec![0.0; total],
            zl: vec![0.0; total],
            zu: vec![0.0; total],
        }
    }
}

/// Solve `min c'v, [A -I] v = 0, lower <= v <= upper`.
pub fn solve(
    matrix: &SparseMatrix,
    cost: &[f64],
    lower: &[f64],
    upper: &[f64],
    deadline: Option<Instant>,
    interrupt: Option<&Arc<AtomicBool>>,
) -> IpmResult {
    let n = matrix.n_cols();
    let m = matrix.n_rows();
    let total = n + m;
    let failed = |status| IpmResult {
        status,
        v: Vec::new(),
        y: Vec::new(),
        z: Vec::new(),
        iterations: 0,
        objective: f64::NAN,
    };

    // Geometric scaling, rounded to powers of two so it is exact.
    let (row_scale, col_scale) = geometric_scaling(matrix);
    // v = sigma * scaled v: structural sigma = col scale, logical sigma = 1 / row scale.
    let sigma: Vec<f64> = (0..total)
        .map(|j| if j < n { col_scale[j] } else { 1.0 / row_scale[j - n] })
        .collect();
    let mut triplets = Vec::with_capacity(matrix.nnz());
    for j in 0..n {
        let (rs, vs) = matrix.column(j);
        for (&i, &v) in rs.iter().zip(vs) {
            triplets.push((i, j, v * row_scale[i] * col_scale[j]));
        }
    }
    let a = SparseMatrix::from_triplets(m, n, triplets);
    let at = a.to_csr();

    let mut kind = Vec::with_capacity(total);
    let mut lo = vec![0.0; total];
    let mut hi = vec![0.0; total];
    let mut c = vec![0.0; total];
    for j in 0..total {
        let l = lower[j] / sigma[j];
        let u = upper[j] / sigma[j];
        let finite_l = lower[j] > -INFINITE_BOUND;
        let finite_u = upper[j] < INFINITE_BOUND;
        lo[j] = if finite_l { l } else { f64::NEG_INFINITY };
        hi[j] = if finite_u { u } else { f64::INFINITY };
        c[j] = cost[j] * sigma[j];
        kind.push(match (finite_l, finite_u) {
            (true, true) if lower[j] == upper[j] => Kind::Fixed,
            (true, true) => Kind::Boxed,
            (true, false) => Kind::Lower,
            (false, true) => Kind::Upper,
            (false, false) => Kind::Free,
        });
        if lower[j] > upper[j] {
            return failed(IpmStatus::NotConverged);
        }
    }

    // Dense columns: their cliques would fill the factor.
    let mean = a.nnz() as f64 / n.max(1) as f64;
    let cut = (64.0f64).max(10.0 * mean);
    let mut by_size: Vec<usize> = (0..n)
        .filter(|&j| a.column(j).0.len() as f64 > cut && kind[j] != Kind::Fixed)
        .collect();
    by_size.sort_by_key(|&j| std::cmp::Reverse(a.column(j).0.len()));
    by_size.truncate(100);
    let mut dense = vec![false; n];
    for &j in &by_size {
        dense[j] = true;
    }
    let mut skip = dense.clone();
    for j in 0..n {
        if kind[j] == Kind::Fixed {
            skip[j] = true;
        }
    }

    // The normal equations in row space, `A Θ_x A' + D`, m by m, or through the same
    // regularized system's other form, `Θ_x^-1 + A' D^-1 A`, n by n, whichever factors
    // smaller. A model with many more rows than columns can fill the first past any
    // budget and the second not at all: `opm2-z10-s4` is 160633 rows by 6250 columns.
    let analysed = Instant::now();
    let stop = Stop { deadline, interrupt: interrupt.map(|f| &**f) };
    let by_rows = Cholesky::analyse(&a, &at, &skip, MAX_FACTOR_NONZEROS, MAX_FACTOR_FLOPS, stop);
    if matches!(by_rows, Err(Abandoned::Stopped)) {
        return failed(IpmStatus::NotConverged);
    }
    let row_size = by_rows.as_ref().map_or(f64::INFINITY, |c| c.nonzeros() as f64);
    let row_skip = vec![false; m];
    // The column form holds at least its diagonal and, for each row, the clique that
    // row makes among its columns. When that alone outgrows the row factor the column
    // ordering is wasted work: on `bab6` and `s250r10`, with 1.3M and 2.4M entry row
    // factors, it took 115s and 76s.
    let widest_row = (0..m).map(|i| at.column(i).0.len()).max().unwrap_or(0) as f64;
    let column_floor = n as f64 + widest_row * (widest_row - 1.0) / 2.0;
    let by_columns = if column_floor < row_size
        && (row_size > (n as f64) * (n as f64) * 0.5 || row_size > 1e6)
    {
        Cholesky::analyse(&at, &a, &row_skip, MAX_FACTOR_NONZEROS.min(row_size), MAX_FACTOR_FLOPS, stop).ok()
    } else {
        None
    };
    let (chol, mode) = match (by_rows, by_columns) {
        (_, Some(chol)) if (chol.nonzeros() as f64) < row_size => (chol, Mode::Columns),
        (Ok(chol), _) => (chol, Mode::Rows),
        (Err(why), _) => {
            if std::env::var_os("RIPSOLVE_TRACE").is_some() {
                eprintln!(
                    "ipm gave up after {:.2}s: m {m} n {n} nnz(A) {}, {why:?}, dense cols {}",
                    analysed.elapsed().as_secs_f64(),
                    a.nnz(),
                    by_size.len()
                );
            }
            return failed(match why {
                Abandoned::Stopped => IpmStatus::NotConverged,
                Abandoned::TooDense { .. } => IpmStatus::TooDense,
            });
        }
    };
    if mode == Mode::Columns {
        // Every column is in the factor of the column form.
        by_size.clear();
        dense.fill(false);
    }
    if std::env::var_os("RIPSOLVE_TRACE").is_some() {
        eprintln!(
            "ipm analyse {:.2}s: m {m} n {n} nnz(A) {} {mode:?} factor nnz {} dense cols {}",
            analysed.elapsed().as_secs_f64(),
            a.nnz(),
            chol.nonzeros(),
            by_size.len()
        );
    }
    let mut ipm = Ipm {
        n,
        m,
        a,
        at,
        cost: c,
        lower: lo,
        upper: hi,
        kind,
        dense,
        dense_list: by_size,
        chol,
        mode,
        row_skip,
        d_inv: vec![0.0; m],
        col_diag: vec![0.0; n],
        work: Vec::new(),
        deadline,
        interrupt,
    };
    let _ = &ipm.dense;
    let mut result = ipm.run(&skip);
    // Unscale.
    if !result.v.is_empty() {
        for j in 0..total {
            result.v[j] *= sigma[j];
            result.z[j] /= sigma[j];
        }
        for i in 0..m {
            result.y[i] *= row_scale[i];
        }
        result.objective = (0..total).map(|j| cost[j] * result.v[j]).sum();
    }
    result
}

/// Row and column scale factors making `R A C`'s entries close to one, by a few
/// passes of geometric-mean scaling.
fn geometric_scaling(a: &SparseMatrix) -> (Vec<f64>, Vec<f64>) {
    let m = a.n_rows();
    let n = a.n_cols();
    let mut r = vec![1.0; m];
    let mut c = vec![1.0; n];
    let pow2 = |x: f64| if x.is_finite() && x > 0.0 { 2f64.powi(x.log2().round() as i32) } else { 1.0 };
    for _ in 0..6 {
        let mut rmin = vec![f64::INFINITY; m];
        let mut rmax = vec![0.0f64; m];
        for j in 0..n {
            let (rs, vs) = a.column(j);
            for (&i, &v) in rs.iter().zip(vs) {
                let x = (v * c[j]).abs();
                rmin[i] = rmin[i].min(x);
                rmax[i] = rmax[i].max(x);
            }
        }
        for i in 0..m {
            if rmax[i] > 0.0 {
                r[i] = 1.0 / (rmin[i] * rmax[i]).sqrt();
            }
        }
        for j in 0..n {
            let (rs, vs) = a.column(j);
            let mut lo = f64::INFINITY;
            let mut hi = 0.0f64;
            for (&i, &v) in rs.iter().zip(vs) {
                let x = (v * r[i]).abs();
                lo = lo.min(x);
                hi = hi.max(x);
            }
            if hi > 0.0 {
                c[j] = 1.0 / (lo * hi).sqrt();
            }
        }
    }
    (r.into_iter().map(pow2).collect(), c.into_iter().map(pow2).collect())
}

impl Ipm<'_> {
    fn total(&self) -> usize {
        self.n + self.m
    }

    fn out_of_time(&self) -> bool {
        self.deadline.is_some_and(|d| Instant::now() >= d)
            || self.interrupt.is_some_and(|f| f.load(Ordering::Relaxed))
    }

    /// `[A -I] v`.
    fn k_mul(&self, v: &[f64], out: &mut [f64]) {
        for i in 0..self.m {
            out[i] = -v[self.n + i];
        }
        for j in 0..self.n {
            let x = v[j];
            if x == 0.0 {
                continue;
            }
            let (rs, vs) = self.a.column(j);
            for (&i, &a) in rs.iter().zip(vs) {
                out[i] += a * x;
            }
        }
    }

    /// `[A -I]' y`.
    fn kt_mul(&self, y: &[f64], out: &mut [f64]) {
        for j in 0..self.n {
            let (rs, vs) = self.a.column(j);
            out[j] = rs.iter().zip(vs).map(|(&i, &a)| a * y[i]).sum();
        }
        for i in 0..self.m {
            out[self.n + i] = -y[i];
        }
    }

    /// `(A Θ_x A' + diag(extra)) p`, every column included.
    fn normal_mul(&self, theta: &[f64], extra: &[f64], p: &[f64], out: &mut [f64]) {
        for i in 0..self.m {
            out[i] = extra[i] * p[i];
        }
        for j in 0..self.n {
            if theta[j] == 0.0 {
                continue;
            }
            let (rs, vs) = self.a.column(j);
            let dot: f64 = rs.iter().zip(vs).map(|(&i, &a)| a * p[i]).sum();
            if dot == 0.0 {
                continue;
            }
            let t = dot * theta[j];
            for (&i, &a) in rs.iter().zip(vs) {
                out[i] += a * t;
            }
        }
    }

    /// Solve the normal equations, by the factor alone when nothing was left out of
    /// it, and otherwise by conjugate gradients preconditioned with it.
    fn normal_solve(&mut self, theta: &[f64], extra: &[f64], rhs: &[f64], skip: &[bool]) -> Vec<f64> {
        let m = self.m;
        let mut x = rhs.to_vec();
        self.precondition(&mut x);
        let _ = skip;
        // Refine: one or more steps of preconditioned conjugate gradients on the full
        // system. With no dense columns the factor is the exact inverse and this
        // converges at once; each dense column costs about one more step.
        let max_steps = match self.mode {
            Mode::Rows if self.dense_list.is_empty() => 2,
            Mode::Rows => self.dense_list.len() + 20,
            Mode::Columns => 5,
        };
        let norm_b = rhs.iter().fold(0.0f64, |acc, &v| acc.max(v.abs())).max(1e-300);
        let mut r = vec![0.0; m];
        self.normal_mul(theta, extra, &x, &mut r);
        for i in 0..m {
            r[i] = rhs[i] - r[i];
        }
        let converged = |r: &[f64]| r.iter().fold(0.0f64, |acc, &v| acc.max(v.abs())) <= 1e-12 * norm_b;
        if converged(&r) {
            return x;
        }
        let mut zv = r.clone();
        self.precondition(&mut zv);
        let mut p = zv.clone();
        let mut rz: f64 = r.iter().zip(&zv).map(|(a, b)| a * b).sum();
        let mut q = vec![0.0; m];
        for _ in 0..max_steps {
            if converged(&r) || rz == 0.0 {
                break;
            }
            self.normal_mul(theta, extra, &p, &mut q);
            let pq: f64 = p.iter().zip(&q).map(|(a, b)| a * b).sum();
            if pq <= 0.0 || !pq.is_finite() {
                break;
            }
            let alpha = rz / pq;
            for i in 0..m {
                x[i] += alpha * p[i];
                r[i] -= alpha * q[i];
            }
            zv.copy_from_slice(&r);
            self.precondition(&mut zv);
            let rz_new: f64 = r.iter().zip(&zv).map(|(a, b)| a * b).sum();
            let beta = rz_new / rz;
            rz = rz_new;
            for i in 0..m {
                p[i] = zv[i] + beta * p[i];
            }
        }
        x
    }

    /// Apply the row form's factor to `r` in place.
    fn precondition(&mut self, r: &mut [f64]) {
        self.chol.solve(r, &mut self.work);
    }

    /// The Newton step's `dy` and `dv` for primal residual `rp` and reduced right-hand
    /// side `rhat`: `(K Θ K' + δ) dy = rp + K Θ r̂` and `dv = Θ (K' dy - r̂)`, fixed
    /// variables held still.
    ///
    /// In the column form the same system is solved for `dx` first,
    /// `(Θ_x^-1 + A' D^-1 A) dx = A' D^-1 (rp - Θ_s r̂_s) - r̂_x`, with
    /// `dy = D^-1 (rp - A dx - Θ_s r̂_s)`. `dx` is used as solved, not rebuilt from `dy`:
    /// going through `dy` multiplies its error by `Θ_x`, which reaches 1e10 on a column
    /// strictly inside its bounds.
    fn newton(&mut self, theta: &[f64], extra: &[f64], rp: &[f64], rhat: &[f64], skip: &[bool]) -> (Vec<f64>, Vec<f64>) {
        let total = self.total();
        let n = self.n;
        let m = self.m;
        let mut dv = vec![0.0; total];
        match self.mode {
            Mode::Rows => {
                let tr: Vec<f64> = (0..total).map(|j| theta[j] * rhat[j]).collect();
                let mut rhs = vec![0.0; m];
                self.k_mul(&tr, &mut rhs);
                for i in 0..m {
                    rhs[i] += rp[i];
                }
                let dy = self.normal_solve(theta, extra, &rhs, skip);
                let mut kty = vec![0.0; total];
                self.kt_mul(&dy, &mut kty);
                for j in 0..total {
                    if self.kind[j] != Kind::Fixed {
                        dv[j] = theta[j] * (kty[j] - rhat[j]);
                    }
                }
                (dy, dv)
            }
            Mode::Columns => {
                let b: Vec<f64> = (0..m).map(|i| self.d_inv[i] * (rp[i] - theta[n + i] * rhat[n + i])).collect();
                let mut bx: Vec<f64> = (0..n)
                    .map(|j| {
                        if self.kind[j] == Kind::Fixed {
                            return 0.0;
                        }
                        let (rs, vs) = self.a.column(j);
                        rs.iter().zip(vs).map(|(&i, &a)| a * b[i]).sum::<f64>() - rhat[j]
                    })
                    .collect();
                let dx = self.column_solve(&mut bx);
                let mut adx = vec![0.0; m];
                for j in 0..n {
                    let x = dx[j];
                    if x == 0.0 {
                        continue;
                    }
                    let (rs, vs) = self.a.column(j);
                    for (&i, &a) in rs.iter().zip(vs) {
                        adx[i] += a * x;
                    }
                }
                let dy: Vec<f64> =
                    (0..m).map(|i| self.d_inv[i] * (rp[i] - adx[i] - theta[n + i] * rhat[n + i])).collect();
                for j in 0..n {
                    if self.kind[j] != Kind::Fixed {
                        dv[j] = dx[j];
                    }
                }
                for i in 0..m {
                    if self.kind[n + i] != Kind::Fixed {
                        dv[n + i] = theta[n + i] * (-dy[i] - rhat[n + i]);
                    }
                }
                (dy, dv)
            }
        }
    }

    /// `(Θ_x^-1 + A' D^-1 A) x = b` by the column form's factor, refined by conjugate
    /// gradients on the same matrix.
    fn column_solve(&mut self, b: &mut [f64]) -> Vec<f64> {
        let n = self.n;
        let m = self.m;
        let apply = |ipm: &Ipm, x: &[f64], out: &mut [f64]| {
            let mut ax = vec![0.0; m];
            for j in 0..n {
                let v = x[j];
                if v == 0.0 {
                    continue;
                }
                let (rs, vs) = ipm.a.column(j);
                for (&i, &a) in rs.iter().zip(vs) {
                    ax[i] += a * v;
                }
            }
            for i in 0..m {
                ax[i] *= ipm.d_inv[i];
            }
            for j in 0..n {
                let (rs, vs) = ipm.a.column(j);
                out[j] = ipm.col_diag[j] * x[j] + rs.iter().zip(vs).map(|(&i, &a)| a * ax[i]).sum::<f64>();
            }
        };
        let mut x = b.to_vec();
        self.chol.solve(&mut x, &mut self.work);
        let norm_b = b.iter().fold(0.0f64, |acc, &v| acc.max(v.abs())).max(1e-300);
        let mut r = vec![0.0; n];
        apply(self, &x, &mut r);
        for j in 0..n {
            r[j] = b[j] - r[j];
        }
        let converged = |r: &[f64]| r.iter().fold(0.0f64, |acc, &v| acc.max(v.abs())) <= 1e-12 * norm_b;
        if converged(&r) {
            return x;
        }
        let mut zv = r.clone();
        self.chol.solve(&mut zv, &mut self.work);
        let mut p = zv.clone();
        let mut rz: f64 = r.iter().zip(&zv).map(|(a, b)| a * b).sum();
        let mut q = vec![0.0; n];
        for _ in 0..3 {
            if converged(&r) || rz == 0.0 {
                break;
            }
            apply(self, &p, &mut q);
            let pq: f64 = p.iter().zip(&q).map(|(a, b)| a * b).sum();
            if pq <= 0.0 || !pq.is_finite() {
                break;
            }
            let alpha = rz / pq;
            for j in 0..n {
                x[j] += alpha * p[j];
                r[j] -= alpha * q[j];
            }
            zv.copy_from_slice(&r);
            self.chol.solve(&mut zv, &mut self.work);
            let rz_new: f64 = r.iter().zip(&zv).map(|(a, b)| a * b).sum();
            let beta = rz_new / rz;
            rz = rz_new;
            for j in 0..n {
                p[j] = zv[j] + beta * p[j];
            }
        }
        x
    }

    /// Factor the normal matrix for `theta`, returning the row form's diagonal `D`.
    fn factor(&mut self, theta: &[f64], skip: &[bool]) -> Option<Vec<f64>> {
        let stop = Stop { deadline: self.deadline, interrupt: self.interrupt.map(|f| &**f) };
        let n = self.n;
        let extra: Vec<f64> = (0..self.m).map(|i| theta[n + i] + DUAL_REG).collect();
        let began = Instant::now();
        match self.mode {
            Mode::Rows => {
                let (a, at) = (&self.a, &self.at);
                if !self.chol.factor(a, at, skip, theta, &extra, stop) {
                    return None;
                }
            }
            Mode::Columns => {
                for i in 0..self.m {
                    self.d_inv[i] = 1.0 / extra[i];
                }
                // A fixed column does not move: an all but infinite diagonal.
                for j in 0..n {
                    self.col_diag[j] = if theta[j] > 0.0 { 1.0 / theta[j] } else { 1e20 };
                }
                let (a, at) = (&self.a, &self.at);
                if !self.chol.factor(at, a, &self.row_skip, &self.d_inv, &self.col_diag, stop) {
                    return None;
                }
            }
        }
        if std::env::var_os("RIPSOLVE_TRACE").is_some() {
            eprintln!("ipm factor {:.3}s tiny {}", began.elapsed().as_secs_f64(), self.chol.tiny_pivots);
        }
        Some(extra)
    }

    /// The Newton direction for complementarity targets `cl`, `cu`.
    #[allow(clippy::too_many_arguments)]
    fn direction(
        &mut self,
        st: &State,
        res: &Residuals,
        theta: &[f64],
        extra: &[f64],
        cl: &[f64],
        cu: &[f64],
        skip: &[bool],
        d: &mut Direction,
    ) {
        let total = self.total();
        // r̂ = rd - cl/xl - (zl/xl) rl + cu/xu - (zu/xu) ru
        let mut rhat = vec![0.0; total];
        for j in 0..total {
            let k = self.kind[j];
            if k == Kind::Fixed {
                continue;
            }
            let mut r = res.rd[j];
            if k.has_lower() {
                r -= (cl[j] + st.zl[j] * res.rl[j]) / st.xl[j];
            }
            if k.has_upper() {
                r += (cu[j] - st.zu[j] * res.ru[j]) / st.xu[j];
            }
            rhat[j] = r;
        }
        let (dy, dv) = self.newton(theta, extra, &res.rp, &rhat, skip);
        d.y = dy;
        d.v = dv;
        for j in 0..total {
            let k = self.kind[j];
            if k.has_lower() {
                d.xl[j] = d.v[j] - res.rl[j];
                d.zl[j] = (cl[j] - st.zl[j] * d.xl[j]) / st.xl[j];
            } else {
                d.xl[j] = 0.0;
                d.zl[j] = 0.0;
            }
            if k.has_upper() {
                d.xu[j] = res.ru[j] - d.v[j];
                d.zu[j] = (cu[j] - st.zu[j] * d.xu[j]) / st.xu[j];
            } else {
                d.xu[j] = 0.0;
                d.zu[j] = 0.0;
            }
        }
    }

    fn step_lengths(&self, st: &State, d: &Direction) -> (f64, f64) {
        let mut ap = 1.0f64;
        let mut ad = 1.0f64;
        for j in 0..self.total() {
            let k = self.kind[j];
            if k.has_lower() {
                if d.xl[j] < 0.0 {
                    ap = ap.min(-st.xl[j] / d.xl[j]);
                }
                if d.zl[j] < 0.0 {
                    ad = ad.min(-st.zl[j] / d.zl[j]);
                }
            }
            if k.has_upper() {
                if d.xu[j] < 0.0 {
                    ap = ap.min(-st.xu[j] / d.xu[j]);
                }
                if d.zu[j] < 0.0 {
                    ad = ad.min(-st.zu[j] / d.zu[j]);
                }
            }
        }
        (ap, ad)
    }

    fn residuals(&self, st: &State) -> Residuals {
        let total = self.total();
        let mut rp = vec![0.0; self.m];
        self.k_mul(&st.v, &mut rp);
        for r in &mut rp {
            *r = -*r;
        }
        let mut kty = vec![0.0; total];
        self.kt_mul(&st.y, &mut kty);
        let mut rd = vec![0.0; total];
        let mut rl = vec![0.0; total];
        let mut ru = vec![0.0; total];
        for j in 0..total {
            let k = self.kind[j];
            if k == Kind::Fixed {
                continue;
            }
            rd[j] = self.cost[j] - kty[j] - st.zl[j] + st.zu[j];
            if k.has_lower() {
                rl[j] = self.lower[j] - st.v[j] + st.xl[j];
            }
            if k.has_upper() {
                ru[j] = self.upper[j] - st.v[j] - st.xu[j];
            }
        }
        // Objectives: the dual's uses the fixed variables' reduced costs, since their
        // columns sit on the right-hand side.
        let pobj: f64 = (0..total).map(|j| self.cost[j] * st.v[j]).sum();
        let mut dobj = 0.0;
        for j in 0..total {
            let k = self.kind[j];
            if k == Kind::Fixed {
                dobj += (self.cost[j] - kty[j]) * st.v[j];
            }
            if k.has_lower() {
                dobj += self.lower[j] * st.zl[j];
            }
            if k.has_upper() {
                dobj -= self.upper[j] * st.zu[j];
            }
        }
        Residuals { rp, rd, rl, ru, pobj, dobj }
    }

    fn mu(&self, st: &State) -> f64 {
        let mut sum = 0.0;
        let mut count = 0usize;
        for j in 0..self.total() {
            let k = self.kind[j];
            if k.has_lower() {
                sum += st.xl[j] * st.zl[j];
                count += 1;
            }
            if k.has_upper() {
                sum += st.xu[j] * st.zu[j];
                count += 1;
            }
        }
        if count == 0 { 0.0 } else { sum / count as f64 }
    }

    fn theta(&self, st: &State) -> Vec<f64> {
        (0..self.total())
            .map(|j| {
                let k = self.kind[j];
                if k == Kind::Fixed {
                    return 0.0;
                }
                let mut d = PRIMAL_REG;
                if k.has_lower() {
                    d += st.zl[j] / st.xl[j];
                }
                if k.has_upper() {
                    d += st.zu[j] / st.xu[j];
                }
                1.0 / d
            })
            .collect()
    }

    /// Mehrotra's starting point: least-squares primal and dual estimates, shifted
    /// into the interior.
    fn start(&mut self, skip: &[bool]) -> Option<State> {
        let total = self.total();
        let n = self.n;
        // A first primal guess inside the bounds, logicals at their rows' activity.
        let mut v = vec![0.0; total];
        for j in 0..n {
            v[j] = clamp_into(0.0, self.lower[j], self.upper[j]);
        }
        let mut act = vec![0.0; self.m];
        let mut xonly = v.clone();
        for i in 0..self.m {
            xonly[n + i] = 0.0;
        }
        self.k_mul(&xonly, &mut act);
        for i in 0..self.m {
            v[n + i] = if self.kind[n + i] == Kind::Fixed {
                self.lower[n + i]
            } else {
                clamp_into(act[i], self.lower[n + i], self.upper[n + i])
            };
        }
        // Least-norm correction to [A -I] v = 0 with Θ = 1 off the fixed variables.
        let theta: Vec<f64> = (0..total).map(|j| if self.kind[j] == Kind::Fixed { 0.0 } else { 1.0 }).collect();
        let extra = self.factor(&theta, skip)?;
        let mut rp = vec![0.0; self.m];
        self.k_mul(&v, &mut rp);
        for r in &mut rp {
            *r = -*r;
        }
        let zero_total = vec![0.0; total];
        let (_, dv) = self.newton(&theta, &extra, &rp, &zero_total, skip);
        for j in 0..total {
            v[j] += dv[j];
        }
        // Dual: y minimizing |c - K'y| over the free part, z the rest.
        let zero_rows = vec![0.0; self.m];
        let (y, _) = self.newton(&theta, &extra, &zero_rows, &self.cost.clone(), skip);
        let mut kty = vec![0.0; total];
        self.kt_mul(&y, &mut kty);

        let mut xl = vec![0.0; total];
        let mut xu = vec![0.0; total];
        let mut zl = vec![0.0; total];
        let mut zu = vec![0.0; total];
        let mut min_x = f64::INFINITY;
        let mut min_z = f64::INFINITY;
        for j in 0..total {
            let k = self.kind[j];
            if k == Kind::Fixed {
                continue;
            }
            let z = self.cost[j] - kty[j];
            if k.has_lower() {
                xl[j] = v[j] - self.lower[j];
                min_x = min_x.min(xl[j]);
            }
            if k.has_upper() {
                xu[j] = self.upper[j] - v[j];
                min_x = min_x.min(xu[j]);
            }
            match k {
                Kind::Boxed => {
                    zl[j] = z.max(0.0);
                    zu[j] = (-z).max(0.0);
                    min_z = min_z.min(zl[j]).min(zu[j]);
                }
                Kind::Lower => {
                    zl[j] = z;
                    min_z = min_z.min(z);
                }
                Kind::Upper => {
                    zu[j] = -z;
                    min_z = min_z.min(-z);
                }
                _ => {}
            }
        }
        let dx = if min_x.is_finite() { (-1.5 * min_x).max(0.0) } else { 0.0 };
        let dz = if min_z.is_finite() { (-1.5 * min_z).max(0.0) } else { 0.0 };
        let (mut xz, mut sx, mut sz) = (0.0, 0.0, 0.0);
        for j in 0..total {
            let k = self.kind[j];
            if k.has_lower() {
                xl[j] += dx;
                zl[j] += dz;
                xz += xl[j] * zl[j];
                sx += xl[j];
                sz += zl[j];
            }
            if k.has_upper() {
                xu[j] += dx;
                zu[j] += dz;
                xz += xu[j] * zu[j];
                sx += xu[j];
                sz += zu[j];
            }
        }
        let dx2 = if sz > 0.0 { 0.5 * xz / sz } else { 1.0 };
        let dz2 = if sx > 0.0 { 0.5 * xz / sx } else { 1.0 };
        for j in 0..total {
            let k = self.kind[j];
            if k.has_lower() {
                xl[j] = (xl[j] + dx2).max(1e-4);
                zl[j] = (zl[j] + dz2).max(1e-4);
            }
            if k.has_upper() {
                xu[j] = (xu[j] + dx2).max(1e-4);
                zu[j] = (zu[j] + dz2).max(1e-4);
            }
        }
        Some(State { v, y, xl, xu, zl, zu })
    }

    fn run(&mut self, skip: &[bool]) -> IpmResult {
        let total = self.total();
        let Some(mut st) = self.start(skip) else {
            return IpmResult {
                status: IpmStatus::NotConverged,
                v: Vec::new(),
                y: Vec::new(),
                z: Vec::new(),
                iterations: 0,
                objective: f64::NAN,
            };
        };
        let norm_c = self.cost.iter().fold(0.0f64, |acc, &c| acc.max(c.abs()));
        let norm_b = (0..total)
            .flat_map(|j| [self.lower[j], self.upper[j]])
            .filter(|b| b.is_finite())
            .fold(0.0f64, |acc, b| acc.max(b.abs()));
        let mut aff = Direction::new(total, self.m);
        let mut cor = Direction::new(total, self.m);
        let mut best_mu = f64::INFINITY;
        let mut stalled = 0;
        let mut iterations = 0;
        let mut status = IpmStatus::NotConverged;
        while iterations < MAX_ITERATIONS {
            let res = self.residuals(&st);
            let mu = self.mu(&st);
            let pinf = res.rp.iter().chain(&res.rl).chain(&res.ru).fold(0.0f64, |a, &r| a.max(r.abs())) / (1.0 + norm_b);
            let dinf = res.rd.iter().fold(0.0f64, |a, &r| a.max(r.abs())) / (1.0 + norm_c);
            let gap = (res.pobj - res.dobj).abs() / (1.0 + res.pobj.abs());
            if std::env::var_os("RIPSOLVE_TRACE").is_some() {
                eprintln!(
                    "ipm {iterations:3} pobj {:+.8e} dobj {:+.8e} pinf {pinf:.1e} dinf {dinf:.1e} gap {gap:.1e} mu {mu:.1e}",
                    res.pobj, res.dobj
                );
            }
            if pinf < TOLERANCE && dinf < TOLERANCE && gap < TOLERANCE {
                status = IpmStatus::Optimal;
                break;
            }
            if !mu.is_finite() || self.out_of_time() {
                break;
            }
            let close = pinf < APPROXIMATE && dinf < APPROXIMATE && gap < APPROXIMATE_GAP;
            if mu < best_mu * 0.99 {
                best_mu = mu;
                stalled = 0;
            } else {
                stalled += 1;
                if stalled > 10 {
                    if close {
                        status = IpmStatus::Approximate;
                    }
                    break;
                }
            }
            if iterations + 1 >= MAX_ITERATIONS && close {
                status = IpmStatus::Approximate;
            }
            iterations += 1;

            let theta = self.theta(&st);
            let Some(extra) = self.factor(&theta, skip) else {
                break;
            };

            // Predictor.
            let cl: Vec<f64> = (0..total).map(|j| -st.xl[j] * st.zl[j]).collect();
            let cu: Vec<f64> = (0..total).map(|j| -st.xu[j] * st.zu[j]).collect();
            self.direction(&st, &res, &theta, &extra, &cl, &cu, skip, &mut aff);
            let (ap, ad) = self.step_lengths(&st, &aff);
            let mut sum = 0.0;
            let mut count = 0usize;
            for j in 0..total {
                let k = self.kind[j];
                if k.has_lower() {
                    sum += (st.xl[j] + ap * aff.xl[j]) * (st.zl[j] + ad * aff.zl[j]);
                    count += 1;
                }
                if k.has_upper() {
                    sum += (st.xu[j] + ap * aff.xu[j]) * (st.zu[j] + ad * aff.zu[j]);
                    count += 1;
                }
            }
            let mu_aff = if count == 0 { 0.0 } else { sum / count as f64 };
            let sigma = if mu > 0.0 { (mu_aff / mu).clamp(0.0, 1.0).powi(3) } else { 0.0 };

            // Corrector.
            let target = sigma * mu;
            let cl: Vec<f64> = (0..total)
                .map(|j| target - st.xl[j] * st.zl[j] - aff.xl[j] * aff.zl[j])
                .collect();
            let cu: Vec<f64> = (0..total)
                .map(|j| target - st.xu[j] * st.zu[j] - aff.xu[j] * aff.zu[j])
                .collect();
            self.direction(&st, &res, &theta, &extra, &cl, &cu, skip, &mut cor);
            let (ap, ad) = self.step_lengths(&st, &cor);
            let ap = (STEP_FRACTION * ap).min(1.0);
            let ad = (STEP_FRACTION * ad).min(1.0);
            for j in 0..total {
                st.v[j] += ap * cor.v[j];
                st.xl[j] += ap * cor.xl[j];
                st.xu[j] += ap * cor.xu[j];
                st.zl[j] += ad * cor.zl[j];
                st.zu[j] += ad * cor.zu[j];
            }
            for i in 0..self.m {
                st.y[i] += ad * cor.y[i];
            }
            if st.v.iter().any(|x| !x.is_finite()) || st.y.iter().any(|x| !x.is_finite()) {
                break;
            }
        }
        let z: Vec<f64> = (0..total).map(|j| st.zl[j] - st.zu[j]).collect();
        let objective = (0..total).map(|j| self.cost[j] * st.v[j]).sum();
        IpmResult { status, v: st.v, y: st.y, z, iterations, objective }
    }
}

struct State {
    v: Vec<f64>,
    y: Vec<f64>,
    xl: Vec<f64>,
    xu: Vec<f64>,
    zl: Vec<f64>,
    zu: Vec<f64>,
}

struct Residuals {
    rp: Vec<f64>,
    rd: Vec<f64>,
    rl: Vec<f64>,
    ru: Vec<f64>,
    pobj: f64,
    dobj: f64,
}

/// `x` moved into `[lo, hi]`, with some room off a finite bound.
fn clamp_into(x: f64, lo: f64, hi: f64) -> f64 {
    if lo.is_finite() && hi.is_finite() {
        if x <= lo || x >= hi { 0.5 * (lo + hi) } else { x }
    } else if lo.is_finite() {
        x.max(lo + 1.0)
    } else if hi.is_finite() {
        x.min(hi - 1.0)
    } else {
        x
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// min -x - y, x + 2y <= 4, 3x + y <= 6, 0 <= x, y <= 10. Optimum (1.6, 1.2), -2.8.
    #[test]
    fn a_small_lp_reaches_its_optimum() {
        let a = SparseMatrix::from_triplets(2, 2, [(0, 0, 1.0), (0, 1, 2.0), (1, 0, 3.0), (1, 1, 1.0)]);
        let cost = [-1.0, -1.0, 0.0, 0.0];
        let lower = [0.0, 0.0, f64::NEG_INFINITY, f64::NEG_INFINITY];
        let upper = [10.0, 10.0, 4.0, 6.0];
        let r = solve(&a, &cost, &lower, &upper, None, None);
        assert_eq!(r.status, IpmStatus::Optimal);
        assert!((r.objective + 2.8).abs() < 1e-6, "{}", r.objective);
        assert!((r.v[0] - 1.6).abs() < 1e-5 && (r.v[1] - 1.2).abs() < 1e-5);
    }

    /// Free variable, equality row, fixed column.
    #[test]
    fn free_fixed_and_equality_columns() {
        // min x0 + 2 x1 + 0 x2, x0 + x1 + x2 = 5 (row 0), x0 - x1 >= -1 (row 1),
        // x0 free but x0 >= -2 via row 2: x0 >= -2; x1 >= 0, x2 fixed at 1.
        let a = SparseMatrix::from_triplets(
            3,
            3,
            [(0, 0, 1.0), (0, 1, 1.0), (0, 2, 1.0), (1, 0, 1.0), (1, 1, -1.0), (2, 0, 1.0)],
        );
        let inf = f64::INFINITY;
        let cost = [1.0, 2.0, 0.0, 0.0, 0.0, 0.0];
        let lower = [-inf, 0.0, 1.0, 5.0, -1.0, -2.0];
        let upper = [inf, inf, 1.0, 5.0, inf, inf];
        // x2 = 1 so x0 + x1 = 4, x0 - x1 >= -1 -> x0 >= 1.5; min x0 + 2 x1 = 4 + x1
        // -> x1 as small as possible: x1 = 0, x0 = 4. Objective 4.
        let r = solve(&a, &cost, &lower, &upper, None, None);
        assert_eq!(r.status, IpmStatus::Optimal);
        assert!((r.objective - 4.0).abs() < 1e-6, "{}", r.objective);
    }
}
