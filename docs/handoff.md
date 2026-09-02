# Where this is, and what to do next

Written at the end of a long session so the next one does not start by re-deriving it.
The reasoning behind every claim here is in `design-notes.md`; this is the index.

## Standing

Measured over the 140 pure binary MIPLIB instances, 60s, 16 threads, one consistent run:

```text
ripsolve     27        HiGHS 46     SCIP 45     CBC 30     commercial 94
```

Started the session at 26. The one gained is `mitre`, at 18.1 seconds, from presolve
probing; nothing was lost. Four instances came back one to five seconds slower and two
came back one to four seconds faster, all of them multi-second parallel solves, which is
the sixteen-thread spread this file already warns about and not a change in the solver:
the probing every one of them pays for was measured separately and is under a fifth of a
second on each.

**The set that matters is 31, not 113.** Of the instances this solver misses, most are
missed by every open source solver too. The 31 that some open source solver closes and
this one does not are the whole of the gap.

## What the comparison says to build

Taken from a reference solver's own logs across all 31, not from first principles.

| what | reach | state |
|---|---|---|
| presolve | 12 of 31 lose >20% of columns | **done, and it converted one of them.** See below |
| root relaxation | 7 of 31 never finish it; 20 reach only 1-3 nodes | **the blocker now** |
| sub-MIP | fires on 21 of 31 | blocked: every construction needs a relaxation, and the instances that need it have none |
| restarts | 9 of 31 | not built |
| symmetry | 6 of 31 | not built |

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

## Next: the root relaxation, through the dual method as a fallback

The largest identified group, and the one the last three sessions have kept arriving at
from different directions. `dual-cold-start` is kept unmerged because the dual method
reaches a *different* optimal vertex and the search downstream reads that vertex —
Gomory cuts come off its tableau, branching reads its fractional values — so on this set
most gaps widened even though four relaxations improved sharply.

That objection applies only where the primal method reaches a vertex at all. On the seven
where it returns nothing, there is no vertex to be worse than. The shape to try is a
second rescue beside `perturbed_root`, which already sits at exactly that point in
`search.rs` and already rescued `neos-1324574` and `tanglegram6` this way: when the
primal has not finished the root, re-enter cold through the dual method rather than
returning nothing. By construction it cannot reach the instances the branch regressed.

Cheap first step before building anything: check out `dual-cold-start`, build, and run
`ripsolve relax` on `air04`, `bley_xl1`, `cod105`, `neos-3226448-wkra`, `supportcase4`,
`ex9` and `ex10`. If those relaxations finish there and do not on main, the fallback is
worth the work; if they do not, this whole line is answered and the notes should say so.

Two things to carry into that work if it happens. The dual entry needs gating so that
warm starts — every node of the search — keep the row selection they have now; the
branch changes them too, and that is a second variable in a measurement that already has
enough. And the branch's steepest edge row selection is most of why its cold starts work
(`drayage-100-23`, 487230 iterations to 2374), so it has to come along.

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
  `docs/baselines/binary-2026-09-02.json`. Diff against that rather than against a
  remembered number, and copy the current one aside before a run that might matter.
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
