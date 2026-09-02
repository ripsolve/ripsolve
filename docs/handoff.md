# Where this is, and what to do next

Written at the end of a long session so the next one does not start by re-deriving it.
The reasoning behind every claim here is in `design-notes.md`; this is the index.

## Standing

Measured over the 140 pure binary MIPLIB instances, 60s, 16 threads, one consistent run:

```text
ripsolve  26 or 27     HiGHS 46     SCIP 45     CBC 30     commercial 94
```

Started the session at 26. The one gained is `mitre`, at about 18 seconds, from presolve
probing.

**The figure is 26 or 27 because `irp` is a coin flip.** It closes in 57 to 60 seconds
against a 60 second limit, so it lands on either side of the line: three consecutive runs
of one build gave optimal, TimeLimit, TimeLimit. Two full benchmark runs of this session's
work came back 27 and then 26, differing in `irp` and nothing else, and the code that
changed between them provably never executes on it. Say which side of the flip a figure
came from, or the next session will chase a regression that is not there. This is the same
trap `nw04` set last session.

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

## Next: reduced cost fixing, for the ten that lose on the bound

The only group with instances within a percent of closing, and the obvious untried thing
for them. With an incumbent and the root's duals in hand, a nonbasic binary whose reduced
cost exceeds the remaining gap cannot take its other value in any better solution and can
be fixed for good. At `nw04`'s 0.4% that is most of the model.

Nothing here reads a reduced cost outside the simplex: `LpSolution` carries the status,
the objective, the primal values and the basis, and no duals. Exposing them is the first
step and is small — the solver already computes `y` by BTRAN of the basic costs every
iteration and has `Solver::reduced_cost`.

Two things to get right, both of which decide whether it is sound:

- It is only valid against a *proven* bound and a *feasible* incumbent, so it belongs
  after the root LP has actually reached Optimal, not after whatever the rescues left.
- The obvious place to apply it is the root, and at the root there is usually no
  incumbent yet: on this set the incumbent arrives during the search. So the version
  worth building re-runs it when the incumbent improves, not once at the start.

Measure it on the ten first — `ripsolve solve` on each at 60s and 16 threads — before
spending two hours on a full benchmark run. A full run is the confirmation, not the
experiment.

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
3. **The parallel search is too noisy to read a regression from.** Incumbent counts at
   16 threads swing enough to invent both losses and gains. Anything that matters is
   confirmed at one thread, or measured with the search taken out of the picture
   entirely — `cargo run --example probecost` does that for presolve, and finding the
   real cost of probing took two sessions longer than it should have because every
   earlier measurement of it had a whole solve sitting on top.

## Measuring here without being misled

- `bench/binary_bench.py [seconds] [threads] [--refresh]` is the standing. `--refresh`
  drops this solver's cached rows only. `bench/out/` is gitignored and every run writes
  over it, so the run behind the figures above is kept as
  `docs/baselines/binary-2026-09-02b.json`, with the run before the root work kept as
  `binary-2026-09-02.json` and the one before probing as `binary-2026-09-01.json`. Diff
  against those rather than against a remembered number, and copy the current one aside
  before a run that might matter.
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
