# Where this is, and what to do next

Written at the end of a long session so the next one does not start by re-deriving it.
The reasoning behind every claim here is in `design-notes.md`; this is the index.

## Standing

Measured over the 140 pure binary MIPLIB instances, 60s, 16 threads:

```text
ripsolve  26 + 2 flips     HiGHS 46     SCIP 45     CBC 30     commercial 94
```

**Do not quote a single run's count.** Two instances sit on the sixty second line, so the
headline reads 26, 27 or 28 depending on where they land, and one benchmark run cannot
tell a real gain of one from noise of one. Quote the set that closes every time, with the
flips named beside it. Measured over five consecutive runs each:

```text
closes every run        26
nw04                     4 runs in 5, 32 to 38s
irp                      2 runs in 5, 46 to 56s
```

Against the same figures at the start of the session: **25 every run, `nw04` never, `irp`
a coin flip.** So the session gained `mitre` outright and turned `nw04` from an instance
that had never closed in any run of any version into one that closes four times in five.
The full benchmark run that produced the current standing shows neither, because it caught
both flips on the wrong side; it is kept as `docs/baselines/binary-2026-09-02c.json`.

**The set that matters is 31, not 113.** Of the instances this solver misses, most are
missed by every open source solver too. The 31 that some open source solver closes and
this one does not are the whole of the gap.

## Where the 31 are, measured rather than inherited

Re-derived this session by running all 31 and reading nodes, cuts and whether any feasible
point was found; the full table is in `design-notes.md` under "The addressable set,
re-derived". The breakdown the last three sessions worked from was stale in two ways, so
do not inherit this one either without checking it.

| group | how many | what they need |
|---|---|---|
| stop at 1-2 nodes, no cuts separated | 11 | the root relaxation, which never finishes |
| search and never find a point | 7 | feasibility |
| hold a point and lose on the bound | 10 | the bound, and **four are within 4%** |
| finish their root, separate cuts, still stop at 2 nodes | 3 | not yet understood |

The last two rows are where the cheap conversions are and are not where the last three
sessions have been looking:

```text
n2seq36f          0.38%   158581 nodes    HiGHS closes it in 3.6s
nw04              0.40%     1574 nodes    CBC closes it in 15.4s
neos-1516309      1.43%   468406 nodes    HiGHS closes it in 0.3s
neos-1599274      3.85%    33419 nodes    HiGHS closes it in 0.6s
```

## Presolve is finished, and it was not the lever it looked like

Probing had been parked twice for cost. The cost was not in the idea: a sweep rebuilt its
state before doing any reasoning, three times over, and a probe does two sweeps per
column. With the rebuilding gone it is 30 to 100 times cheaper, finds *more* than it used
to, and is bounded by a budget measured in matrix entries read rather than in probes.

What that bought is one instance. Ten times the budget finds several thousand more
columns on `ex9`, `ex10`, `air04`, `rail01` and `neos-4754521-awarau` and closes none of
them, because all of them stop at one node with the root relaxation unfinished. `air04`
fixes 1343 columns of 8904 — a reference solver fixes 1400 — and still spends sixty
seconds on 37994 simplex iterations without leaving the root.

So the twelve instances that lose a fifth of their columns are not twelve conversions
waiting on presolve. They are waiting on the same thing the twenty root bound instances
are waiting on, and that is where the next work goes.

## The root relaxation was attacked this session and is answered as far as it goes

Two defects and one capability, none of which converted anything, all of which are in.

The rescue that already existed for a stalled root was given the caller's entire
remaining clock, so on all eight models that need it, it spent the whole run and returned
nothing — and nothing after it could ever run. It is now bounded by what the attempt it
repeats spent. A second rescue, the dual method entered cold, was taken from
`dual-cold-start` and gated so warm starts keep the row selection they had. `air04`'s
relaxation goes from not finishing in 130 seconds to 1.5, `tanglegram6`'s to 0.5, and in
the search `air04` goes from one node to 380 and `tanglegram6` from no incumbent to one.
Five of the seven still do not finish either way.

Nothing closes, and the reason is worth carrying forward: at a 60 second limit the first
attempt is given 40% of the run before anything else may be tried, so `air04`'s root is
answered at 25 seconds. But do not go straight at `ROOT_LP_FIRST_SHARE`. `air04` needs
well over 60 seconds even from a free root — the branch reached a 0.62% gap only at 120
seconds — and `tanglegram6`'s bound is 0 against an incumbent of 8856. This group is not
where the next conversion is.

## Reduced cost fixing is built, at the root only, and that is where to continue

Built this session and it is what moved `nw04`. At an optimal basis a nonbasic column's
reduced cost `d` says the objective cannot fall by leaving the bound it sits on, so moving
`t` off it raises the objective by at least `|d| t`; anything better than the incumbent
`u` needs `root + |d| t < u`, capping `t` at `(u - root) / |d|`, which for an integer
column rounds inwards and on a binary usually decides it. `nw04` fixes 25471 columns of
87482 this way, `neos-1516309` 350 of 4500.

**The obvious continuation is to re-run it as the incumbent improves.** It currently runs
once, at the root, against whatever incumbent the root heuristics found, and that
incumbent is far weaker than the one the search ends with: `n2seq36f` is at a 39.7% gap at
the root and 0.38% at the end. The room `(u - root)` is widest exactly when the fixing runs
and narrowest when it would pay most, so most of this reduction is still on the table.

Two obstacles, found by scoping it and worth having before starting rather than after:

- **Globally**, the model is shared immutably across the search's threads, so a
  mid-search tightening has nowhere to live that a running worker can read. A restart is
  the shape that fits, and restarts are separately on the list for nine of the 31.
- **Per node**, where it would otherwise be the textbook answer, `Node.fixings` is not
  free to grow: `search.rs` uses `node.fixings.len()` as the node's *depth* in the
  best-bound pool's tie-break. Appending implied bounds to it would silently corrupt node
  ordering, and on `nw04` it would append 25471 entries to a vector that is cloned per
  child. Implied bounds need their own field, kept out of the depth measure, and a cap.

The cost side is already dealt with: `Lp::reduced_costs` reuses the factorization the
basis was left with, so asked right after a node's own solve it is one BTRAN and a pass
over the columns, not a refactorization.

Two things not to get wrong, both learned the hard way here and both written up in
`design-notes.md`:

- **Any margin belongs on the generous side.** A smaller room caps the travel harder and
  fixes more, which is the direction that removes a solution nobody has seen. Subtracting
  the feasibility tolerance "to be safe" fixed 1566 columns of `irp` where the sound
  version fixes 3, and it passed the tests.
- **The invariance test only catches a proof wrong in kind.** Solving every sample with
  the fixing on and off catches the travel cap set to zero, and only once the generated
  families were added; halving the cap is caught by neither fixture set.

## And after that: the seven that search and never find a point

Untouched this session and the second largest group. `acc-tight4`, `acc-tight5`,
`graph20-20-1rand`, `air04`, `air05`, `neos-3045796-mogo` and `neos-820879` reach between
twelve and 8720 nodes and never find a feasible assignment, with rounding, diving, the
pump, fixing with propagation and the LP-free jump all running and all returning nothing.

## Parked, with the evidence for each in design-notes.md

- `stash@{0}` — row equilibration. 6-13% more LP iterations per second where
  coefficients vary, no solves, nothing to scale on the five all-ones matrices that
  most need help.
- `longstep-phase1`, `dual-cold-start` — older LP work. The second is now the next thing
  to try rather than a dead end; see above.
- `probing-parked` — **superseded**, and can be deleted. Its reduction is on the main
  line and larger.

## Three lessons this session kept re-learning

1. **A share of the time limit is not a budget.** It prices work against how long the
   caller happens to be willing to wait. Made three times before this session, and the
   fix each time was to measure against the model. Probing now adds the other half of
   that lesson: a multiple of the matrix is not a budget either, because it hands a
   model with eight million nonzeros a hundred times the work it gives one with eighty
   thousand, for the same reduction. What a budget encodes is how much reasoning is
   worth doing, which is a property of neither the clock nor the matrix.
2. **Cheapest first is wrong when a heuristic is cheap because its answers are bad.**
   Made three times; such a heuristic belongs last, where it runs only if everything
   else failed.
3. **The parallel search is too noisy to read a regression from, and now too noisy to
   read the headline count from.** Incumbent counts at 16 threads swing enough to invent
   both losses and gains. Anything that matters is confirmed at one thread, or measured
   with the search taken out of the picture entirely — `cargo run --example probecost`
   does that for presolve, and finding the real cost of probing took two sessions longer
   than it should have because every earlier measurement of it had a whole solve sitting
   on top. This session the noise ran the other way for the first time: a single
   benchmark run reported that reduced cost fixing changed nothing, when five runs of the
   instance it was aimed at show four closes where there had never been one.

## Measuring here without being misled

- `bench/binary_bench.py [seconds] [threads] [--refresh]` is the standing. `--refresh`
  drops this solver's cached rows only. `bench/out/` is gitignored and every run writes
  over it, so the run behind the figures above is kept as
  `docs/baselines/binary-2026-09-02c.json`, with the run before reduced cost fixing kept
  as `binary-2026-09-02b.json`, the one before the root work as `binary-2026-09-02.json`
  and the one before probing as `binary-2026-09-01.json`. Diff against those rather than
  against a remembered number, and copy the current one aside before a run that matters.
- `cargo run --release --example probecost -- <models>` times presolve with and without
  probing on the same model and reports both reductions. Use it in preference to a solve
  for anything about presolve: it is seconds rather than an hour, and it has no search
  noise on top.
- Never run two measurements at once. A benchmark and an A/B on the same machine corrupt
  each other, and the tell is an instance that suddenly cannot solve something it solves
  in isolation. A full 140 instance run takes around two hours, so plan what goes into it
  before starting one.
- A suspended laptop writes impossible wall times into the cache. Entries above the time
  limit by a wide margin are that, and should be dropped and re-measured.
