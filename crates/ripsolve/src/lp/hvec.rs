//! A dense vector that also lists where its nonzeros are.
//!
//! The vectors a simplex solves for are mostly zero: 46 nonzeros in a right-hand side
//! of 131991 on `ns1116954`. Held as a plain dense array, every pass over one costs
//! the full length whatever it holds, and a pivot makes a few dozen such passes, so on
//! a model of fifty thousand rows the pivot costs milliseconds where the work it does
//! is microseconds. Carrying the positions alongside the values lets each pass read
//! only those.
//!
//! The list is a superset: a position on it may hold a zero, after cancellation, but
//! every nonzero is on it, and none is on it twice. A vector too full for the list to
//! be worth keeping drops it and is read densely, which is what `dense` says.

/// See the module notes.
#[derive(Clone, Debug, Default)]
pub struct HVec {
    /// The values, by position, the full length always.
    pub val: Vec<f64>,
    /// Positions that may be nonzero, when `dense` is false.
    pub idx: Vec<usize>,
    /// Whether `idx` has been given up and `val` must be read whole.
    pub dense: bool,
}

impl HVec {
    pub fn new(m: usize) -> Self {
        HVec {
            val: vec![0.0; m],
            idx: Vec::new(),
            dense: false,
        }
    }

    pub fn len(&self) -> usize {
        self.val.len()
    }

    pub fn is_empty(&self) -> bool {
        self.val.is_empty()
    }

    /// Zero the vector at length `m`, reading only what the list says was written.
    pub fn reset(&mut self, m: usize) {
        if self.val.len() != m || self.dense {
            self.val.clear();
            self.val.resize(m, 0.0);
        } else {
            for &i in &self.idx {
                self.val[i] = 0.0;
            }
        }
        self.idx.clear();
        self.dense = false;
    }

    /// The unit vector `e_r` at length `m`.
    pub fn set_unit(&mut self, m: usize, r: usize) {
        self.reset(m);
        self.val[r] = 1.0;
        self.idx.push(r);
    }

    /// Copy a dense vector in, listing its nonzeros.
    pub fn set_dense(&mut self, values: &[f64]) {
        self.val.clear();
        self.val.extend_from_slice(values);
        self.reindex();
    }

    /// Make this a copy of `other`, reading only what its list holds.
    pub fn copy_from(&mut self, other: &HVec) {
        if other.dense {
            self.val.clear();
            self.val.extend_from_slice(&other.val);
            self.idx.clear();
            self.dense = true;
            return;
        }
        self.reset(other.len());
        for &i in &other.idx {
            self.val[i] = other.val[i];
        }
        self.idx.extend_from_slice(&other.idx);
    }

    /// Rebuild the list from the values, by a pass over all of them.
    pub fn reindex(&mut self) {
        self.idx.clear();
        self.idx.extend((0..self.val.len()).filter(|&i| self.val[i] != 0.0));
        self.dense = false;
    }

    /// Add `v` at position `i`, listing it if it is not listed already.
    ///
    /// A sum that cancels to exactly zero and is then added to again must not list its
    /// position twice, and nothing marks it as listed but the value. So that is where
    /// the list is given up instead. Listed twice, it was counted twice by every reader:
    /// the long step's flips summed on `n2seq36q`, whose entries are all ones, moved
    /// the basic values by thousands, and the dual method proved infeasible node LPs
    /// whose optimum was the root's.
    pub fn add(&mut self, i: usize, v: f64) {
        let was = self.val[i];
        if !self.dense && was == 0.0 {
            self.idx.push(i);
        }
        self.val[i] += v;
        if was != 0.0 && self.val[i] == 0.0 {
            self.dense = true;
        }
    }

    /// How many positions the list holds, or the length when there is no list.
    pub fn count(&self) -> usize {
        if self.dense { self.val.len() } else { self.idx.len() }
    }

    /// Call `f` on every position that may be nonzero, with its value.
    #[inline]
    pub fn for_each(&self, mut f: impl FnMut(usize, f64)) {
        if self.dense {
            for (i, &v) in self.val.iter().enumerate() {
                if v != 0.0 {
                    f(i, v);
                }
            }
        } else {
            for &i in &self.idx {
                let v = self.val[i];
                if v != 0.0 {
                    f(i, v);
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_sum_that_cancels_is_never_listed_twice() {
        let mut v = HVec::new(4);
        v.add(2, 1.0);
        v.add(2, -1.0);
        v.add(2, 1.0);
        v.add(1, 3.0);
        let mut seen = Vec::new();
        v.for_each(|i, x| seen.push((i, x)));
        seen.sort_by_key(|&(i, _)| i);
        assert_eq!(seen, vec![(1, 3.0), (2, 1.0)]);
    }
}
