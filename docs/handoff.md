# Where this is, and what to do next

Written at the end of a long session so the next one does not start by re-deriving it.
The reasoning behind every claim here is in `design-notes.md`; this is the index.

## Standing

Measured over the 140 pure binary MIPLIB instances, 60s, 16 threads:

```text
ripsolve     33        HiGHS 46     SCIP 45     CBC 30     commercial 94
```

Confirmed by a full run, `docs/baselines/binary-2026-09-07.json`, which gained
`neos-3045796-mogo` against the previous one and lost nothing. Ahead of CBC, and
thirteen behind HiGHS.

The 33 all close on every run, with no instance left sitting on the sixty second line.
That matters, because for a long stretch two of them were: quote the set that closes
every time rather than one run's count, and check a new one against repeated runs before
believing it.

Against the start of these sessions, **25 every run, with `nw04` never closing in any run
of any version and `irp` a coin flip**. The six gained:

```text
mitre           presolve probing, once the propagator stopped rebuilding itself
nw04            reduced cost fixing at the root
irp             those fixings collected as the incumbent earns them
neos-1516309    the Gomory density filter
neos-1599274    the same, plus not forfeiting optimality to a prunable skipped node
n2seq36f        mod-2 cuts
neos-3226448-wkra  four budgets between the feasibility search and a point it could reach
neos-3045796-mogo  the neighbourhood the root's reduced costs would fix
```

Five of the 33 are instances HiGHS does not close: `disctom`, `eil33-2`, `decomp2`,
`nw04` and `neos-4754521-awarau`.

## Speed, on the instances both solvers close

Worth knowing before optimising anything, and it is not one number. Measured over the 25
instances this solver and HiGHS both close, with this solver's time net of MPS parsing,
since the benchmark times HiGHS around `run()` with reading outside the timer and this
solver as a whole subprocess:

```text
                      ripsolve   HiGHS
shifted geomean (10)      4.37    1.92
shifted geomean (1)       1.91    1.17
median ratio              0.94x
faster on                14 of 25
```

**The median instance is a dead heat and the distribution is bimodal.** This solver is
two to ten times *faster* on the easy end -- `app2-1` 0.29s against 3.10, `p0201` 0.20
against 1.20, `neos-3437289-erdre` 0.09 against 1.10 -- and far behind on a handful of
hard ones, which is what drags the geometric mean to 1.28x and the shifted mean to 2.3x:

```text
acc-tight2        37.4s vs 0.1     250x
neos-1599274      22.8  vs 0.6      35x
mitre             29.0  vs 1.4      20x
neos-913984       35.9  vs 4.3       8x
```

So per-node overhead is not the problem and is not where to look. The gap is entirely in
the heavy machinery, which is the same conclusion the missing-instance analysis reaches
from the other direction. `acc-tight2` is the one worth a look on its own account:
HiGHS closes it in a tenth of a second where this solver needs 37 and two nodes, which
smells like a reduction rather than a search.

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

## Where the cut work got to

`n2seq36f` closes. The section below is kept because the *reasoning* in it was wrong in
an instructive way and the measurements in it stand.

**The answer was not aggregation.** `HighsPathSeparator` aggregates along continuous
columns, and `n2seq36f` has none, so it produces nothing there; naming it as the target
was reasoning from a name rather than from the file. `HighsModkSeparator` is what moves
that bound, and mod-2 separation is now built: take tight rows with integer coefficients,
add a subset whose coefficient sums are all even and whose right-hand side is odd, halve
and round down. A Python prototype settled it in an afternoon -- 52000 to 52200 in two
rounds and eighteen cuts -- before any Rust was written.

Aggregation over continuous columns is still untried and is still the obvious next cut
family, but it is aimed at *mixed* models and there is no measurement here saying it
would pay on this set, which is pure binary by construction.

## The reasoning that led there, with one wrong turn

The first specific target this file has had for the bound group, and it comes from
watching HiGHS solve `n2seq36f` rather than from first principles.

Both solvers start from the same relaxation bound of 52000 against an optimum of 52200.
HiGHS moves the bound to 52200 at two seconds with 1120 cuts, 101 of them held in the LP,
and that closes the whole gap from the bound side. This solver generates 1019 cuts in the
same place and moves the bound by **nothing at all**, to the unit.

Four explanations were tested and all four are wrong:

- Not the *number* of cuts. 4036 of them move it zero.
- Not how long they are kept. The loop ages a cut out after two rounds sitting slack;
  holding every one of them instead still moves it zero.
- Not reduced cost fixing, and not a restart. `n2seq36f` narrows to 78% of its columns
  decided and the bound is still 52000. `cargo run --example restartsim` simulates the
  restart in a few minutes rather than building it, and says not to build it.
- Not the density filter, though that one was a real defect and fixing it closed two
  other instances. With it out of the way `n2seq36f` produces twenty GMI cuts, each
  violated by a full unit, and the bound stays at 52000 to the digit.

What is left is what the cuts *are*. `n2seq36f` has 285 rows and 8100 columns and an
enormous optimal face; every cut here removes one vertex of it and the LP steps to another
vertex worth the same 52000. All four families in `cuts.rs` -- cover, MIR, clique, GMI --
separate from a single row or a single tableau row. Cutting a face needs aggregation: MIR
over *combinations* of rows, which in HiGHS is `HighsPathSeparator` over `HighsLpAggregator`.

That is the build. It is the largest identified group (ten instances lose on the bound,
four of them within 4%), and it is the only lever on that group that has not been tried
and measured this session.

## Reduced cost fixing is built, root and beyond, and is not the bound group's answer

Built this session and it is what moved `nw04`. At an optimal basis a nonbasic column's
reduced cost `d` says the objective cannot fall by leaving the bound it sits on, so moving
`t` off it raises the objective by at least `|d| t`; anything better than the incumbent
`u` needs `root + |d| t < u`, capping `t` at `(u - root) / |d|`, which for an integer
column rounds inwards and on a binary usually decides it. `nw04` fixes 25471 columns of
87482 this way, `neos-1516309` 350 of 4500.

It now runs beyond the root as well, as a table rather than a mutation, which is how the
"shared immutably across threads" obstacle dissolves: nothing in the derivation changes
except the room, so the bound each column will take and the incumbent at which it takes it
are both worked out at the root. Entries are ordered by the incumbent they need, so the
ones in force are a prefix and a worker recomputes its length only when the incumbent
moves. On `n2seq36f`, 5910 of 6642 come into force and the node count over a minute
doubles.

**It is not what the bound group needs, and that is now measured rather than assumed.**
`n2seq36f`'s bound never leaves 52000 however many columns are fixed underneath it. See
the section above for what does move it.

Still genuinely open, if someone wants the rest of the reduction: node-local fixing, using
each node's own LP duals against its own bound, which is stronger than the root's. One
obstacle to know about first -- `search.rs` uses `node.fixings.len()` as the node's *depth*
in the best-bound pool's tie-break, so implied bounds cannot be appended to that vector
without silently corrupting node ordering, and on `nw04` it would append 25471 entries to
something cloned once per child. They need their own field and a cap. The cost side is
already dealt with: `Lp::reduced_costs` reuses the factorization the basis was left with,
so asked right after a node's own solve it is one BTRAN and a pass over the columns.

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
  `docs/baselines/binary-2026-09-07.json`, with the run before the reduced cost
  neighbourhood kept as `binary-2026-09-06.json`, the one before the feasibility work as
  `binary-2026-09-05.json`, the one before mod-2 cuts as
  `binary-2026-09-04.json`, the one before the cut density fix as
  `binary-2026-09-03.json`, the one before this round's reduced cost work as
  `binary-2026-09-02c.json`, the one before the root work as `binary-2026-09-02.json` and
  the one before probing as `binary-2026-09-01.json`. Diff against those rather than
  against a remembered number, and copy the current one aside before a run that matters.
- `cargo run --release --example probecost -- <models>` times presolve with and without
  probing on the same model and reports both reductions. Use it in preference to a solve
  for anything about presolve: it is seconds rather than an hour, and it has no search
  noise on top.
- Never run two measurements at once. A benchmark and an A/B on the same machine corrupt
  each other, and the tell is an instance that suddenly cannot solve something it solves
  in isolation. A full 140 instance run takes around two hours, so plan what goes into it
  before starting one.
- `pb-fit2d` runs for around 300 seconds against a sixty second limit, and has in every
  run recorded here, so it is a standing defect rather than a suspended laptop. Something
  below the search does not check the clock on that model. `n3seq24` overruns to 87
  seconds; nothing else in the set is out by more than a few.
- A suspended laptop writes impossible wall times into the cache. Entries above the time
  limit by a wide margin are that, and should be dropped and re-measured.
