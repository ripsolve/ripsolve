//! Sparse Cholesky factorization of the interior-point method's normal equations.
//!
//! Every iteration of the interior-point method solves `(A Θ A' + D) y = r` for a new
//! diagonal `Θ`. The pattern of that matrix never changes, so the ordering and the
//! pattern of the factor are worked out once ([`Cholesky::analyse`]) and each
//! iteration only refills the numbers ([`Cholesky::factor`]).
//!
//! # Ordering
//!
//! Approximate minimum degree on the quotient graph, with the columns of `A` as the
//! starting elements: column `j` of `A` puts a clique on its rows into `A A'`, which is
//! exactly what an element is, so `A A'` is never formed to order it. Supervariables
//! (rows that have become indistinguishable) are merged, and an element whose
//! variables all sit inside the new pivot's element is absorbed. Rows of the normal
//! matrix that are dense go last.
//!
//! # Factor
//!
//! Supernodal and left-looking. Columns whose structures nest form a supernode stored
//! as one dense column-major block, so the inner loops are dense and vectorize where a
//! column-at-a-time factor would chase indices.
//!
//! A pivot that comes out tiny or negative is replaced by a huge one, which zeroes that
//! component of the solution instead of failing. The interior-point method produces
//! these near its end, when `Θ` spans thirty orders of magnitude; it is the standard
//! treatment and what keeps the method running to its tolerance.

use crate::sparse::SparseMatrix;

/// Pivots at or below this times their own column's diagonal before elimination are
/// treated as zero: what is left after cancelling all but this share is rounding.
/// Relative to the column's own diagonal and not the largest anywhere, because the
/// column form gives a fixed column a diagonal of 1e20.
const TINY_PIVOT: f64 = 1e-15;
/// What replaces a tiny pivot: large enough that the solve gives that component zero.
const HUGE_PIVOT: f64 = 1e128;
/// Columns factored together inside a supernode's diagonal block.
const PANEL: usize = 16;

/// The symbolic and numeric factor `P (A Θ A' + D) P' = L L'`.
#[derive(Clone, Debug)]
pub struct Cholesky {
    m: usize,
    /// Step `k` eliminates original row `perm[k]`; `pinv` is the inverse.
    perm: Vec<usize>,
    pinv: Vec<usize>,
    /// Supernode `s` holds factor columns `sn_start[s]..sn_start[s + 1]`.
    sn_start: Vec<usize>,
    /// Supernode of each factor column.
    sn_of: Vec<usize>,
    /// Row structure of each supernode, in factor (permuted) indices, sorted: the
    /// supernode's own columns first, then the rows below.
    rows_start: Vec<usize>,
    rows: Vec<usize>,
    /// Offset of each supernode's dense block in `values`.
    block_start: Vec<usize>,
    values: Vec<f64>,
    /// Pivots the last factorization replaced by [`HUGE_PIVOT`].
    pub tiny_pivots: usize,
}

/// Why an analysis returned no factor.
#[derive(Debug)]
pub enum Abandoned {
    /// The factor would hold more entries or cost more arithmetic than allowed.
    TooDense { nonzeros: f64, flops: f64 },
    /// The caller's deadline passed or its interrupt was raised.
    Stopped,
}

/// When to give up partway: a deadline and an interrupt flag, both optional.
#[derive(Clone, Copy, Default)]
pub struct Stop<'a> {
    pub deadline: Option<std::time::Instant>,
    pub interrupt: Option<&'a std::sync::atomic::AtomicBool>,
}

impl Stop<'_> {
    pub fn hit(&self) -> bool {
        self.deadline.is_some_and(|d| std::time::Instant::now() >= d)
            || self.interrupt.is_some_and(|f| f.load(std::sync::atomic::Ordering::Relaxed))
    }
}

/// How many steps of a long loop run between looks at the clock.
const CLOCK_EVERY: usize = 1024;

impl Cholesky {
    /// Order and analyse the normal matrix of `a`, leaving out the columns in `skip`.
    ///
    /// Gives up when the factor would hold more than `max_nonzeros` entries.
    pub fn analyse(
        a: &SparseMatrix,
        at: &SparseMatrix,
        skip: &[bool],
        max_nonzeros: f64,
        max_flops: f64,
        stop: Stop,
    ) -> Result<Cholesky, Abandoned> {
        let m = a.n_rows();
        let order = amd(a, at, skip, stop).ok_or(Abandoned::Stopped)?;
        let mut pinv = vec![0usize; m];
        for (k, &r) in order.iter().enumerate() {
            pinv[r] = k;
        }

        // Elimination tree and column counts on that order, then a postorder, which
        // makes every fundamental supernode a run of consecutive columns.
        let (parent, _) = etree_and_counts(a, at, skip, &order, &pinv, None, stop);
        if stop.hit() {
            return Err(Abandoned::Stopped);
        }
        let post = postorder(&parent);
        let perm: Vec<usize> = post.iter().map(|&k| order[k]).collect();
        for (k, &r) in perm.iter().enumerate() {
            pinv[r] = k;
        }
        let (parent, counts) = etree_and_counts(a, at, skip, &perm, &pinv, Some(max_nonzeros), stop);
        if stop.hit() {
            return Err(Abandoned::Stopped);
        }
        let Some(counts) = counts else {
            return Err(Abandoned::TooDense {
                nonzeros: f64::INFINITY,
                flops: f64::INFINITY,
            });
        };
        let nonzeros: f64 = counts.iter().map(|&c| c as f64).sum();
        let flops: f64 = counts.iter().map(|&c| (c as f64) * (c as f64)).sum();
        if nonzeros > max_nonzeros || flops > max_flops {
            return Err(Abandoned::TooDense { nonzeros, flops });
        }

        // Fundamental supernodes: column k + 1 continues k's supernode when it is k's
        // parent and k's structure is exactly k plus k + 1's.
        let mut sn_start = vec![0usize];
        for k in 1..m {
            let continues = parent[k - 1] == k && counts[k - 1] == counts[k] + 1 && k - sn_start[sn_start.len() - 1] < 128;
            if !continues {
                sn_start.push(k);
            }
        }
        if m > 0 {
            sn_start.push(m);
        } else {
            sn_start = vec![0];
        }
        let ns = sn_start.len() - 1;
        let mut sn_of = vec![0usize; m];
        for s in 0..ns {
            for k in sn_start[s]..sn_start[s + 1] {
                sn_of[k] = s;
            }
        }

        // Row structure of each supernode's first column, by walking each row's
        // subtree: row i of the factor has a nonzero in column t exactly when t is
        // on the path up the tree from some j < i with (P M P')[i][j] nonzero.
        let mut rows_start = vec![0usize; ns + 1];
        for s in 0..ns {
            rows_start[s + 1] = rows_start[s] + counts[sn_start[s]];
        }
        let mut rows = vec![0usize; rows_start[ns]];
        let mut fill: Vec<usize> = rows_start[..ns].to_vec();
        let mut mark = vec![usize::MAX; m];
        for i in 0..m {
            if i % CLOCK_EVERY == 0 && stop.hit() {
                return Err(Abandoned::Stopped);
            }
            mark[i] = i;
            let s = sn_of[i];
            if sn_start[s] == i {
                rows[fill[s]] = i;
                fill[s] += 1;
            }
            let r = perm[i];
            let (cols, _) = at.column(r);
            for &j in cols {
                if skip[j] {
                    continue;
                }
                let (rs, _) = a.column(j);
                for &r2 in rs {
                    let mut t = pinv[r2];
                    if t >= i {
                        continue;
                    }
                    while mark[t] != i {
                        mark[t] = i;
                        let s = sn_of[t];
                        if sn_start[s] == t {
                            rows[fill[s]] = i;
                            fill[s] += 1;
                        }
                        t = parent[t];
                        if t == usize::MAX {
                            break;
                        }
                    }
                }
            }
        }
        debug_assert!((0..ns).all(|s| fill[s] == rows_start[s + 1]));

        let mut block_start = vec![0usize; ns + 1];
        for s in 0..ns {
            let width = sn_start[s + 1] - sn_start[s];
            block_start[s + 1] = block_start[s] + width * (rows_start[s + 1] - rows_start[s]);
        }
        Ok(Cholesky {
            m,
            perm,
            pinv,
            sn_start,
            sn_of,
            rows_start,
            rows,
            values: vec![0.0; block_start[ns]],
            block_start,
            tiny_pivots: 0,
        })
    }

    /// Entries in the factor.
    pub fn nonzeros(&self) -> usize {
        self.values.len()
    }

    /// Factor `A diag(theta) A' + diag(extra)`, the columns `skip` left out.
    pub fn factor(
        &mut self,
        a: &SparseMatrix,
        at: &SparseMatrix,
        skip: &[bool],
        theta: &[f64],
        extra: &[f64],
        stop: Stop,
    ) -> bool {
        let m = self.m;
        let ns = self.sn_start.len() - 1;
        self.tiny_pivots = 0;
        // Each column's diagonal as assembled, before any update reaches it.
        let mut diag0: Vec<f64> = Vec::new();
        // Where each factor row sits in the current supernode's block.
        let mut position = vec![usize::MAX; m];
        // Supernodes waiting to update the supernode their next row falls in, as linked
        // lists: `head[s]` and `next[d]`. `ptr[d]` is the first row of `d`'s structure
        // that has not yet been applied.
        let mut head = vec![usize::MAX; ns];
        let mut next = vec![usize::MAX; ns];
        let mut ptr = vec![0usize; ns];
        let mut update: Vec<f64> = Vec::new();

        for s in 0..ns {
            if s % 64 == 0 && stop.hit() {
                return false;
            }
            let f = self.sn_start[s];
            let l = self.sn_start[s + 1];
            let w = l - f;
            let rs = &self.rows[self.rows_start[s]..self.rows_start[s + 1]];
            let ld = rs.len();
            for (p, &i) in rs.iter().enumerate() {
                position[i] = p;
            }
            let base = self.block_start[s];
            diag0.clear();
            diag0.resize(w, 0.0);
            {
                let block = &mut self.values[base..base + w * ld];
                block.fill(0.0);
                // Assemble the lower part of the normal matrix's columns f..l.
                for c in f..l {
                    let col = &mut block[(c - f) * ld..(c - f + 1) * ld];
                    let r = self.perm[c];
                    col[position[c]] += extra[r];
                    let (cols, vals) = at.column(r);
                    for (&j, &arj) in cols.iter().zip(vals) {
                        if skip[j] || theta[j] == 0.0 {
                            continue;
                        }
                        let scale = arj * theta[j];
                        let (rs2, vs2) = a.column(j);
                        for (&r2, &a2) in rs2.iter().zip(vs2) {
                            let i = self.pinv[r2];
                            if i >= c {
                                col[position[i]] += scale * a2;
                            }
                        }
                    }
                    diag0[c - f] = col[position[c]];
                }
            }

            // Updates from every earlier supernode with rows in f..l.
            let mut d = head[s];
            while d != usize::MAX {
                let after = next[d];
                let drs = &self.rows[self.rows_start[d]..self.rows_start[d + 1]];
                let dw = self.sn_start[d + 1] - self.sn_start[d];
                let dld = drs.len();
                let p = ptr[d];
                let mut q = p;
                while q < dld && drs[q] < l {
                    q += 1;
                }
                let len = dld - p;
                let cols = q - p;
                let (before, from_s) = self.values.split_at_mut(base);
                let dblock = &before[self.block_start[d]..self.block_start[d] + dw * dld];
                let block = &mut from_s[..w * ld];
                // Rows p.. of d landing on consecutive rows of s, the usual case once
                // the factor turns dense, take the update in place; otherwise it goes
                // through a buffer and is scattered.
                let pos0 = position[drs[p]];
                if position[drs[dld - 1]] - pos0 == len - 1 {
                    let col0 = drs[p] - f;
                    lower_product(&mut block[col0 * ld..], ld, pos0, len, cols, dblock, dld, dw, p);
                } else {
                    update.clear();
                    update.resize(len * cols, 0.0);
                    lower_product(&mut update, len, 0, len, cols, dblock, dld, dw, p);
                    for c in 0..cols {
                        let col_in_s = drs[p + c] - f;
                        let col = &mut block[col_in_s * ld..(col_in_s + 1) * ld];
                        for (t, &i) in update[c * len + c..(c + 1) * len].iter().zip(&drs[p + c..]) {
                            col[position[i]] += *t;
                        }
                    }
                }
                ptr[d] = q;
                if q < dld {
                    let target = self.sn_of[drs[q]];
                    next[d] = head[target];
                    head[target] = d;
                }
                d = after;
            }

            // Dense factor of the diagonal block, then the rows below it, in panels: each
            // panel left-looking within itself, then the columns to its right updated
            // from it by the same kernel as above.
            let block = &mut self.values[base..base + w * ld];
            let mut ps = 0;
            while ps < w {
                let pe = (ps + PANEL).min(w);
                for c in ps..pe {
                    for k in ps..c {
                        let lck = block[k * ld + c];
                        if lck == 0.0 {
                            continue;
                        }
                        let (left, right) = block.split_at_mut(c * ld);
                        let src = &left[k * ld + c..k * ld + ld];
                        let dst = &mut right[c..ld];
                        for (t, &v) in dst.iter_mut().zip(src) {
                            *t -= v * lck;
                        }
                    }
                    let col = &mut block[c * ld..(c + 1) * ld];
                    let mut pivot = col[c];
                    if !(pivot > TINY_PIVOT * diag0[c] && pivot > 1e-300) {
                        pivot = HUGE_PIVOT;
                        self.tiny_pivots += 1;
                    }
                    let root = pivot.sqrt();
                    col[c] = root;
                    let inv = 1.0 / root;
                    for v in &mut col[c + 1..] {
                        *v *= inv;
                    }
                }
                if pe < w {
                    let (left, right) = block.split_at_mut(pe * ld);
                    let panel = &left[ps * ld..];
                    lower_product(right, ld, pe, ld - pe, w - pe, panel, ld, pe - ps, pe);
                }
                ps = pe;
            }
            for &i in rs {
                position[i] = usize::MAX;
            }
            ptr[s] = w;
            if w < ld {
                let target = self.sn_of[rs[w]];
                next[s] = head[target];
                head[target] = s;
            }
        }
        true
    }

    /// Solve `(A Θ A' + D) x = b` in place with the last factorization.
    pub fn solve(&self, b: &mut [f64], work: &mut Vec<f64>) {
        let m = self.m;
        work.clear();
        work.extend(self.perm.iter().map(|&r| b[r]));
        let x = work;
        let ns = self.sn_start.len() - 1;
        // L z = P b.
        for s in 0..ns {
            let f = self.sn_start[s];
            let w = self.sn_start[s + 1] - f;
            let rs = &self.rows[self.rows_start[s]..self.rows_start[s + 1]];
            let ld = rs.len();
            let block = &self.values[self.block_start[s]..self.block_start[s] + w * ld];
            for c in 0..w {
                let col = &block[c * ld..(c + 1) * ld];
                let v = x[f + c] / col[c];
                x[f + c] = v;
                if v != 0.0 {
                    for (&i, &lic) in rs[c + 1..].iter().zip(&col[c + 1..]) {
                        x[i] -= lic * v;
                    }
                }
            }
        }
        // L' y = z.
        for s in (0..ns).rev() {
            let f = self.sn_start[s];
            let w = self.sn_start[s + 1] - f;
            let rs = &self.rows[self.rows_start[s]..self.rows_start[s + 1]];
            let ld = rs.len();
            let block = &self.values[self.block_start[s]..self.block_start[s] + w * ld];
            for c in (0..w).rev() {
                let col = &block[c * ld..(c + 1) * ld];
                let mut v = x[f + c];
                for (&i, &lic) in rs[c + 1..].iter().zip(&col[c + 1..]) {
                    v -= lic * x[i];
                }
                x[f + c] = v / col[c];
            }
        }
        for k in 0..m {
            b[self.perm[k]] = x[k];
        }
    }

    /// The factor's column `k`'s supernode, for tests.
    #[cfg(test)]
    fn supernodes(&self) -> usize {
        self.sn_start.len() - 1
    }

    #[cfg(test)]
    fn widest(&self) -> usize {
        self.sn_start.windows(2).map(|w| w[1] - w[0]).max().unwrap_or(0)
    }

    #[cfg(test)]
    fn sn_of_len(&self) -> usize {
        self.sn_of.len()
    }
}

/// `out[c * out_ld + row0 + i] -= sum_k L[p + i][k] * L[p + c][k]` for `i >= c`,
/// over the `dw` columns of the block `l` (leading dimension `ld`), `i < len` and
/// `c < cols`. Rows a little above the diagonal of each column may be written too;
/// they sit in the unused upper triangle and are never read.
///
/// The update is most of a factorization's work. It runs as 8-row by 4-column tiles
/// held in registers across the whole sum, so each entry loaded serves several
/// multiply-adds; one column at a time, it was bound by memory traffic. The AVX2 build
/// of the same loop is picked at run time where the processor has it.
#[allow(clippy::too_many_arguments)]
fn lower_product(out: &mut [f64], out_ld: usize, row0: usize, len: usize, cols: usize, l: &[f64], ld: usize, dw: usize, p: usize) {
    #[cfg(target_arch = "x86_64")]
    {
        if std::is_x86_feature_detected!("avx2") && std::is_x86_feature_detected!("fma") {
            // SAFETY: the features the function is compiled for were just detected.
            unsafe { lower_product_avx2(out, out_ld, row0, len, cols, l, ld, dw, p) };
            return;
        }
    }
    lower_product_body(out, out_ld, row0, len, cols, l, ld, dw, p, |a, b, c| a * b + c);
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2,fma")]
#[allow(clippy::too_many_arguments)]
fn lower_product_avx2(out: &mut [f64], out_ld: usize, row0: usize, len: usize, cols: usize, l: &[f64], ld: usize, dw: usize, p: usize) {
    lower_product_body(out, out_ld, row0, len, cols, l, ld, dw, p, |a, b, c| a.mul_add(b, c));
}

#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn lower_product_body(
    out: &mut [f64],
    out_ld: usize,
    row0: usize,
    len: usize,
    cols: usize,
    l: &[f64],
    ld: usize,
    dw: usize,
    p: usize,
    fma: impl Fn(f64, f64, f64) -> f64 + Copy,
) {
    const R: usize = 8;
    const C: usize = 4;
    let mut c0 = 0;
    while c0 + C <= cols {
        let mut i0 = c0;
        while i0 + R <= len {
            let mut acc = [[0.0f64; R]; C];
            for k in 0..dw {
                let col = &l[k * ld + p..k * ld + ld];
                let b: [f64; C] = col[c0..c0 + C].try_into().unwrap();
                let a: [f64; R] = col[i0..i0 + R].try_into().unwrap();
                for q in 0..C {
                    for r in 0..R {
                        acc[q][r] = fma(a[r], b[q], acc[q][r]);
                    }
                }
            }
            for q in 0..C {
                let at = (c0 + q) * out_ld + row0 + i0;
                let t = &mut out[at..at + R];
                for r in 0..R {
                    t[r] -= acc[q][r];
                }
            }
            i0 += R;
        }
        // Rows left over at the bottom.
        for i in i0..len {
            let mut acc = [0.0f64; C];
            for k in 0..dw {
                let col = &l[k * ld + p..k * ld + ld];
                for q in 0..C {
                    acc[q] = fma(col[i], col[c0 + q], acc[q]);
                }
            }
            for q in 0..C {
                out[(c0 + q) * out_ld + row0 + i] -= acc[q];
            }
        }
        c0 += C;
    }
    for c in c0..cols {
        let at = c * out_ld + row0;
        let target = &mut out[at + c..at + len];
        for k in 0..dw {
            let lk = &l[k * ld + p + c..k * ld + ld];
            let a = -lk[0];
            if a == 0.0 {
                continue;
            }
            for (t, &v) in target.iter_mut().zip(lk) {
                *t = fma(v, a, *t);
            }
        }
    }
}

/// Elimination tree of `P A A' P'` (columns `skip` left out), and when `limit` is
/// given, the count of every factor column including its diagonal, `None` if those
/// counts pass `limit` in total.
fn etree_and_counts(
    a: &SparseMatrix,
    at: &SparseMatrix,
    skip: &[bool],
    perm: &[usize],
    pinv: &[usize],
    limit: Option<f64>,
    stop: Stop,
) -> (Vec<usize>, Option<Vec<usize>>) {
    let m = perm.len();
    let mut parent = vec![usize::MAX; m];
    let mut ancestor = vec![usize::MAX; m];
    for k in 0..m {
        if k % CLOCK_EVERY == 0 && stop.hit() {
            return (parent, None);
        }
        let (cols, _) = at.column(perm[k]);
        for &j in cols {
            if skip[j] {
                continue;
            }
            let (rs, _) = a.column(j);
            for &r2 in rs {
                let mut t = pinv[r2];
                if t >= k {
                    continue;
                }
                // Liu's algorithm with path compression.
                while t != usize::MAX && t != k {
                    let up = ancestor[t];
                    ancestor[t] = k;
                    if up == usize::MAX {
                        parent[t] = k;
                        break;
                    }
                    t = up;
                }
            }
        }
    }
    let Some(limit) = limit else {
        return (parent, None);
    };
    let mut counts = vec![1usize; m];
    let mut mark = vec![usize::MAX; m];
    let mut total = m as f64;
    for i in 0..m {
        if i % CLOCK_EVERY == 0 && stop.hit() {
            return (parent, None);
        }
        mark[i] = i;
        let (cols, _) = at.column(perm[i]);
        for &j in cols {
            if skip[j] {
                continue;
            }
            let (rs, _) = a.column(j);
            for &r2 in rs {
                let mut t = pinv[r2];
                if t >= i {
                    continue;
                }
                while t != usize::MAX && mark[t] != i {
                    mark[t] = i;
                    counts[t] += 1;
                    total += 1.0;
                    t = parent[t];
                }
            }
        }
        if total > limit {
            return (parent, None);
        }
    }
    (parent, Some(counts))
}

/// A postorder of the forest `parent`: `post[k]` is the node visited `k`-th.
fn postorder(parent: &[usize]) -> Vec<usize> {
    let m = parent.len();
    let mut first_child = vec![usize::MAX; m];
    let mut sibling = vec![usize::MAX; m];
    for k in (0..m).rev() {
        let p = parent[k];
        if p != usize::MAX {
            sibling[k] = first_child[p];
            first_child[p] = k;
        }
    }
    let mut post = Vec::with_capacity(m);
    let mut stack = Vec::new();
    for root in 0..m {
        if parent[root] != usize::MAX {
            continue;
        }
        stack.push(root);
        while let Some(&top) = stack.last() {
            let child = first_child[top];
            if child == usize::MAX {
                post.push(top);
                stack.pop();
            } else {
                // Detach the child so the node is finished once all are visited.
                first_child[top] = sibling[child];
                stack.push(child);
            }
        }
    }
    post
}

/// Approximate minimum degree ordering of the rows of `a` for `A A'`, the columns
/// `skip` left out. Returns the rows in elimination order.
fn amd(a: &SparseMatrix, at: &SparseMatrix, skip: &[bool], stop: Stop) -> Option<Vec<usize>> {
    let m = a.n_rows();
    let n = a.n_cols();
    // Node ids: variables (rows) 0..m, which become elements when eliminated, and the
    // starting elements, one per usable column, m + j.
    let total = m + n;
    let mut elems: Vec<Vec<usize>> = vec![Vec::new(); m];
    let mut members: Vec<Vec<usize>> = vec![Vec::new(); total];
    let mut alive_elem = vec![false; total];
    let mut elem_size = vec![0usize; total];
    for j in 0..n {
        if skip[j] {
            continue;
        }
        let (rs, _) = a.column(j);
        if rs.len() < 2 {
            continue;
        }
        members[m + j] = rs.to_vec();
        alive_elem[m + j] = true;
        elem_size[m + j] = rs.len();
        for &r in rs {
            elems[r].push(m + j);
        }
    }
    let _ = at;
    let mut nv = vec![1usize; m];
    let mut eliminated = vec![false; m];
    // Supervariable members, in the order they will be numbered.
    let mut merged: Vec<Vec<usize>> = vec![Vec::new(); m];

    // Exact starting degrees.
    let mut degree = vec![0usize; m];
    let mut stamp = vec![usize::MAX; total.max(m)];
    for i in 0..m {
        let mut d = 0;
        stamp[i] = i;
        for &e in &elems[i] {
            for &r in &members[e] {
                if stamp[r] != i {
                    stamp[r] = i;
                    d += 1;
                }
            }
        }
        degree[i] = d;
    }

    // Rows dense in A A' go last, out of the graph.
    let dense_cut = (16.0f64).max(10.0 * (m as f64).sqrt()) as usize;
    let mut dense = vec![false; m];
    let mut last = Vec::new();
    for i in 0..m {
        if degree[i] > dense_cut && m > 100 {
            dense[i] = true;
            eliminated[i] = true;
            last.push(i);
        }
    }
    if !last.is_empty() {
        for e in m..total {
            if alive_elem[e] {
                members[e].retain(|&r| !dense[r]);
                elem_size[e] = members[e].len();
            }
        }
        for i in 0..m {
            if dense[i] {
                continue;
            }
            let mut d = 0usize;
            stamp[i] = usize::MAX - 1 - i;
            let tag = usize::MAX - 1 - i;
            for &e in &elems[i] {
                for &r in &members[e] {
                    if stamp[r] != tag {
                        stamp[r] = tag;
                        d += 1;
                    }
                }
            }
            degree[i] = d.saturating_sub(1);
        }
    } else {
        for d in &mut degree {
            *d = d.saturating_sub(1);
        }
    }
    stamp.fill(usize::MAX);

    // Degree buckets as doubly linked lists.
    let mut head = vec![usize::MAX; m + 1];
    let mut next = vec![usize::MAX; m];
    let mut prev = vec![usize::MAX; m];
    let insert = |i: usize, d: usize, head: &mut [usize], next: &mut [usize], prev: &mut [usize]| {
        next[i] = head[d];
        prev[i] = usize::MAX;
        if head[d] != usize::MAX {
            prev[head[d]] = i;
        }
        head[d] = i;
    };
    let remove = |i: usize, d: usize, head: &mut [usize], next: &mut [usize], prev: &mut [usize]| {
        if prev[i] != usize::MAX {
            next[prev[i]] = next[i];
        } else {
            head[d] = next[i];
        }
        if next[i] != usize::MAX {
            prev[next[i]] = prev[i];
        }
    };
    for i in 0..m {
        if !eliminated[i] {
            degree[i] = degree[i].min(m);
            insert(i, degree[i], &mut head, &mut next, &mut prev);
        }
    }

    let mut order = Vec::with_capacity(m);
    let mut left: usize = (0..m).filter(|&i| !eliminated[i]).count();
    let mut min_degree = 0;
    let mut we = vec![0isize; total];
    let mut we_stamp = vec![usize::MAX; total];
    let mut in_pivot = vec![usize::MAX; m];
    let mut step = 0usize;
    let mut lp: Vec<usize> = Vec::new();
    let mut hashes: Vec<(u64, usize)> = Vec::new();

    while left > 0 {
        if step.is_multiple_of(CLOCK_EVERY) && stop.hit() {
            return None;
        }
        while min_degree <= m && head[min_degree] == usize::MAX {
            min_degree += 1;
        }
        let p = head[min_degree];
        remove(p, min_degree, &mut head, &mut next, &mut prev);
        step += 1;

        // The new element: everything p's elements touch.
        lp.clear();
        in_pivot[p] = step;
        let p_elems = std::mem::take(&mut elems[p]);
        for &e in &p_elems {
            if !alive_elem[e] {
                continue;
            }
            for &r in &members[e] {
                if nv[r] > 0 && !eliminated[r] && in_pivot[r] != step {
                    in_pivot[r] = step;
                    lp.push(r);
                }
            }
            alive_elem[e] = false;
            members[e] = Vec::new();
        }
        eliminated[p] = true;
        order.push(p);
        order.extend(merged[p].iter().copied());
        left -= nv[p];

        // |L_e \ Lp| for every element next to Lp.
        for &i in &lp {
            remove(i, degree[i], &mut head, &mut next, &mut prev);
            for &e in &elems[i] {
                if !alive_elem[e] {
                    continue;
                }
                if we_stamp[e] != step {
                    we_stamp[e] = step;
                    we[e] = elem_size[e] as isize;
                }
                we[e] -= nv[i] as isize;
            }
        }
        let degme: usize = lp.iter().map(|&i| nv[i]).sum();

        // New element lists and approximate degrees.
        hashes.clear();
        for &i in &lp {
            let mut external = 0usize;
            let list = &mut elems[i];
            list.retain(|&e| {
                if !alive_elem[e] {
                    return false;
                }
                let w = we[e];
                if w <= 0 {
                    // Inside the new element: absorbed.
                    alive_elem[e] = false;
                    return false;
                }
                external += w as usize;
                true
            });
            list.push(p);
            let bound = degree[i] + degme - nv[i];
            let d = (external + degme - nv[i]).min(bound).min(left - nv[i]);
            degree[i] = d;
            list.sort_unstable();
            let mut h: u64 = 0;
            for &e in list.iter() {
                h = h.wrapping_mul(1_000_003).wrapping_add(e as u64);
            }
            hashes.push((h, i));
        }
        // Supervariables: rows with identical element lists.
        hashes.sort_unstable();
        let mut a = 0;
        while a < hashes.len() {
            let mut b = a + 1;
            while b < hashes.len() && hashes[b].0 == hashes[a].0 {
                b += 1;
            }
            for x in a..b {
                let i = hashes[x].1;
                if nv[i] == 0 {
                    continue;
                }
                for y in x + 1..b {
                    let j = hashes[y].1;
                    if nv[j] == 0 || elems[i] != elems[j] {
                        continue;
                    }
                    nv[i] += nv[j];
                    degree[i] = degree[i].saturating_sub(nv[j]);
                    nv[j] = 0;
                    let moved = std::mem::take(&mut merged[j]);
                    merged[i].push(j);
                    merged[i].extend(moved);
                    elems[j] = Vec::new();
                }
            }
            a = b;
        }

        lp.retain(|&i| nv[i] > 0);
        for &i in &lp {
            let d = degree[i].min(m);
            degree[i] = d;
            insert(i, d, &mut head, &mut next, &mut prev);
            min_degree = min_degree.min(d);
        }
        members[p] = lp.clone();
        alive_elem[p] = true;
        elem_size[p] = degme;
        // Elements absorbed above were only flagged; drop their lists.
        for &e in &p_elems {
            members[e] = Vec::new();
        }
    }
    order.extend(last);
    debug_assert_eq!(order.len(), m);
    Some(order)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn random_matrix(m: usize, n: usize, per_col: usize, seed: u64) -> SparseMatrix {
        let mut state = seed;
        let mut rand = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        let mut triplets = Vec::new();
        for j in 0..n {
            for _ in 0..per_col {
                let i = (rand() % m as u64) as usize;
                let v = ((rand() % 2000) as f64 - 1000.0) / 250.0;
                triplets.push((i, j, v));
            }
        }
        SparseMatrix::from_triplets(m, n, triplets)
    }

    #[test]
    fn the_factor_solves_the_normal_equations() {
        for seed in 1..46u64 {
            // The last few are large and fill heavily, which makes wide supernodes:
            // the paths through the panel factor and the in-place updates.
            let (m, n, per_col) = if seed >= 40 {
                (150 + 40 * (seed as usize - 40), 400, 4 + (seed as usize % 3))
            } else {
                let m = 5 + (seed as usize * 7) % 60;
                (m, m + (seed as usize * 3) % 50, 1 + (seed as usize) % 4)
            };
            let a = random_matrix(m, n, per_col, seed);
            let at = a.to_csr();
            let skip = vec![false; n];
            let mut chol = Cholesky::analyse(&a, &at, &skip, 1e9, 1e15, Stop::default()).unwrap();
            assert_eq!(chol.sn_of_len(), m);
            assert!(chol.supernodes() >= 1);
            if seed >= 40 {
                assert!(chol.widest() > 2 * PANEL, "seed {seed}: widest supernode {}", chol.widest());
            }
            let theta: Vec<f64> = (0..n).map(|j| 0.5 + (j % 5) as f64).collect();
            let extra: Vec<f64> = (0..m).map(|i| 1e-3 + (i % 3) as f64 * 0.1).collect();
            assert!(chol.factor(&a, &at, &skip, &theta, &extra, Stop::default()));
            assert_eq!(chol.tiny_pivots, 0);
            let x: Vec<f64> = (0..m).map(|i| (i as f64 * 0.37).sin()).collect();
            // b = (A Θ A' + D) x
            let mut b: Vec<f64> = (0..m).map(|i| extra[i] * x[i]).collect();
            for j in 0..n {
                let (rs, vs) = a.column(j);
                let dot: f64 = rs.iter().zip(vs).map(|(&i, &v)| v * x[i]).sum();
                for (&i, &v) in rs.iter().zip(vs) {
                    b[i] += theta[j] * dot * v;
                }
            }
            let mut work = Vec::new();
            chol.solve(&mut b, &mut work);
            for i in 0..m {
                assert!((b[i] - x[i]).abs() < 1e-6, "seed {seed} row {i}: {} vs {}", b[i], x[i]);
            }
        }
    }

    #[test]
    fn the_ordering_is_a_permutation() {
        for seed in 1..30u64 {
            let a = random_matrix(200, 300, 3, seed);
            let at = a.to_csr();
            let skip: Vec<bool> = (0..300).map(|j| j % 17 == 0).collect();
            let order = amd(&a, &at, &skip, Stop::default()).unwrap();
            let mut seen = vec![false; 200];
            for &r in &order {
                assert!(!seen[r]);
                seen[r] = true;
            }
            assert_eq!(order.len(), 200);
        }
    }
}
