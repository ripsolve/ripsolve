# Where this is, and what to do next

Written at the end of a long session so the next one does not start by re-deriving it.
The reasoning behind every claim here is in `design-notes.md`; this is the index.

## Standing

Measured over the 140 pure binary MIPLIB instances, 60s, 16 threads, one consistent run:

```text
ripsolve     26        HiGHS 46     SCIP 45     CBC 30     commercial 94
```

Started the session at 24. The three gained are `acc-tight2`, `disctom` and
`neos-913984`, all from the LP-free feasibility search, each closing in one node.
`nw04` left the list and is a coin flip in both versions, timing out in three runs of
each: it was a lucky parallel run in the old count, not a capability.

**The set that matters is 31, not 114.** Of the instances this solver misses, most are
missed by every open source solver too. The 31 that some open source solver closes and
this one does not are the whole of the gap, and passing the best of them means
converting about 23 of those 31.

## What the comparison says to build

Taken from a reference solver's own logs across all 31, not from first principles. Every
component invented from first principles in this session produced no solves; everything
that worked came from reading those logs or from a measurement already written down and
not followed up.

| what | reach | state |
|---|---|---|
| presolve | 12 of 31 lose >20% of columns; `ex9` and `ex10` are solved outright by it | **best lever**, partly built on `probing-parked` |
| sub-MIP | fires on 21 of 31 | blocked: every construction needs a relaxation, and the instances that need it have none |
| restarts | 9 of 31 | not built |
| symmetry | 6 of 31 | not built |

## Parked, with the evidence for each in design-notes.md

- `probing-parked` — presolve probing. Reduces `mitre` by 3677 columns of 10724 against
  a reference solver's 4693, `air04` by 787, and is the only thing tried that closes
  `mitre`. Costs `air03` 3.7x its solve. Six budgets tried, all trading the reduction
  against the cost. Needs propagation that does not rebuild per probe, not a seventh
  budget.
- `stash@{0}` — row equilibration. 6-13% more LP iterations per second where
  coefficients vary, no solves, nothing to scale on the five all-ones matrices that
  most need help.
- `longstep-phase1`, `dual-cold-start` — older LP work.

## Three lessons this session kept re-learning

1. **A share of the time limit is not a budget.** It prices work against how long the
   caller happens to be willing to wait. Made three times: the root heuristic budget,
   the jump budget, and probing. The fix each time was to measure against the model, the
   relaxation's own cost or the reach of a column.
2. **Cheapest first is wrong when a heuristic is cheap because its answers are bad.**
   Made three times: fix-and-propagate ahead of diving cost `cap6000` 500 nodes to 1500,
   and the same was nearly done with the jump and with the corner points. Such a
   heuristic belongs last, where it runs only if everything else failed.
3. **The parallel search is too noisy to read a regression from.** Incumbent counts at 16
   threads swing enough to invent both losses and gains. Anything that matters is
   confirmed at one thread.

## Measuring here without being misled

- `bench/binary_bench.py [seconds] [threads] [--refresh]` is the standing. `--refresh`
  drops this solver's cached rows only. `bench/out/` is gitignored and every run writes
  over it, so the run behind the figures above is kept as
  `docs/baselines/binary-2026-09-01.json`. Diff against that rather than against a
  remembered number, and copy the current one aside before a run that might matter.
- Never run two measurements at once. A benchmark and an A/B on the same machine
  corrupt each other, and the tell is an instance that suddenly cannot solve something
  it solves in isolation.
- A suspended laptop writes impossible wall times into the cache. Entries above the time
  limit by a wide margin are that, and should be dropped and re-measured.
