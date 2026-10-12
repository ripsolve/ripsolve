//! Rows supplied during the search, against the same rows written into the model.
//!
//! Clique partitioning is the model a separator exists for: one binary per pair, "these
//! two share a group", a degree row per point, and a transitivity row per ordered triple,
//! `y[a,x] + y[x,c] - y[a,c] <= 1`. The triples are cubic in the point count and nearly
//! all slack, so they are the textbook lazy constraint.

use ripsolve::cuts::Cut;
use ripsolve::model::{Builder, RowSense, Sense};
use ripsolve::search::{self, Options, Point, Status};
use ripsolve::{Problem, VarType};

struct Partition {
    n: usize,
    cost: Vec<Vec<f64>>,
    size: usize,
}

impl Partition {
    fn new(n: usize, size: usize, seed: u64) -> Self {
        let mut state = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        // Scattered points in a square, which leaves the relaxation a crowd of
        // half-groups to choose between and the search a tree to walk.
        let at: Vec<(f64, f64)> = (0..n)
            .map(|_| ((next() % 100) as f64, (next() % 100) as f64))
            .collect();
        let cost = (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| ((at[i].0 - at[j].0).hypot(at[i].1 - at[j].1)).round())
                    .collect()
            })
            .collect();
        Self { n, cost, size }
    }

    fn pair(&self, i: usize, j: usize) -> usize {
        let (i, j) = if i < j { (i, j) } else { (j, i) };
        i * self.n - i * (i + 1) / 2 + (j - i - 1)
    }

    fn triangle(&self, a: usize, x: usize, c: usize) -> Cut {
        let mut coefficients = vec![
            (self.pair(a, x), 1.0),
            (self.pair(x, c), 1.0),
            (self.pair(a, c), -1.0),
        ];
        coefficients.sort_by_key(|&(j, _)| j);
        Cut {
            coefficients,
            lb: f64::NEG_INFINITY,
            ub: 1.0,
        }
    }

    fn triangles(&self) -> Vec<Cut> {
        let n = self.n;
        let mut out = Vec::new();
        for x in 0..n {
            for a in (0..n).filter(|&a| a != x) {
                for c in (a + 1..n).filter(|&c| c != x) {
                    out.push(self.triangle(a, x, c));
                }
            }
        }
        out
    }

    fn model(&self, full: bool) -> Problem {
        let n = self.n;
        let mut b = Builder::new(Sense::Minimize);
        let mut objective = Vec::new();
        for i in 0..n {
            for j in i + 1..n {
                let c = b.column(&format!("y{i}_{j}"), 0.0, 1.0, VarType::Integer);
                objective.push((c, self.cost[i][j]));
            }
        }
        b.objective(&objective);
        for i in 0..n {
            let row: Vec<_> = (0..n)
                .filter(|&j| j != i)
                .map(|j| (self.pair(i, j), 1.0))
                .collect();
            b.row(&row, RowSense::Eq, (self.size - 1) as f64);
        }
        if full {
            for t in self.triangles() {
                b.row(&t.coefficients, RowSense::Le, 1.0);
            }
        }
        b.build()
    }

    fn separate(&self, x: &[f64], _: Point) -> Vec<Cut> {
        self.triangles()
            .into_iter()
            .filter(|t| t.violation(x) > 1e-6)
            .collect()
    }

    /// Whether `x` describes a partition into groups of the right size.
    fn is_partition(&self, x: &[f64]) -> bool {
        let n = self.n;
        let together = |i: usize, j: usize| x[self.pair(i, j)] > 0.5;
        (0..n).all(|i| {
            let group: Vec<usize> = (0..n).filter(|&j| j != i && together(i, j)).collect();
            group.len() == self.size - 1
                && group
                    .iter()
                    .all(|&j| group.iter().all(|&k| k == j || together(j, k)))
        })
    }
}

fn options(threads: usize) -> Options {
    Options {
        threads,
        gap_tolerance: 0.0,
        ..Options::default()
    }
}

#[test]
fn lazy_rows_reach_the_optimum_of_the_full_model() {
    for seed in 1..=6 {
        for threads in [1, 4] {
            let partition = Partition::new(16, 4, seed);
            let full = search::solve(&partition.model(true), options(threads));
            let separator = |x: &[f64], point: Point| partition.separate(x, point);
            let lazy = search::solve_with(
                &partition.model(false),
                options(threads),
                None,
                &separator,
            );
            assert_eq!(full.status, Status::Optimal, "seed {seed}");
            assert_eq!(lazy.status, Status::Optimal, "seed {seed}, {threads} threads");
            assert_eq!(lazy.objective, full.objective, "seed {seed}, {threads} threads");
            assert!(partition.is_partition(&lazy.x), "seed {seed}, {threads} threads");
        }
    }
}

#[test]
fn without_its_rows_the_model_has_a_cheaper_answer_that_is_not_a_partition() {
    // The premise of the test above: the degree rows alone admit something the triangles
    // forbid, so the separator really is doing the work.
    let partition = Partition::new(16, 4, 3);
    let bare = search::solve(&partition.model(false), options(1));
    let full = search::solve(&partition.model(true), options(1));
    assert!(bare.objective < full.objective);
    assert!(!partition.is_partition(&bare.x));
}

#[test]
fn a_start_that_breaks_a_supplied_row_is_not_believed() {
    let partition = Partition::new(16, 4, 3);
    let model = partition.model(false);
    let bare = search::solve(&model, options(1));
    let separator = |x: &[f64], point: Point| partition.separate(x, point);
    let lazy = search::solve_with(&model, options(1), Some(&bare.x), &separator);
    assert!(lazy.objective > bare.objective);
    assert!(partition.is_partition(&lazy.x));
}

#[test]
fn cuts_offered_beside_a_complete_model_change_nothing_but_the_work() {
    for seed in 1..=4 {
        let partition = Partition::new(16, 4, seed);
        let model = partition.model(true);
        let plain = search::solve(&model, options(1));
        let separator = |x: &[f64], point: Point| partition.separate(x, point);
        let offered = search::solve_with(&model, options(1), None, &separator);
        assert_eq!(offered.objective, plain.objective, "seed {seed}");
    }
}
