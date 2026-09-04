# Where this is, and what to do next

Written at the end of a long session so the next one does not start by re-deriving it.
The reasoning behind every claim here is in `design-notes.md`; this is the index.

## Where the dual work had got to when it stopped

Nothing uncommitted; this is a note on the last measurement so it is not re-run. The
primal side is settled and answered — see "The primal side is not the gap" in
`design-notes.md` — and what is left on the three closest instances is bound. Running
each cut family alone against the same relaxation, with `cargo run --release --example
cutfamilies`:

```text
neos-953928     root -99.9200   cover 19, mir 65, clique 62, gomory 12, mod2 0
                                158 cuts between them and the bound does not move at all
neos-820879     root 24874.27   gomory 82 -> 24959.6, mod2 6 -> 24887, all -> 24969.8
                                optimum 25468
air05           root 25877.61   gomory generates *nothing*, clique 17 -> 25906.4
                                optimum 26374
```

Two of the three are now answered, and neither is short of a filter setting; see
"`air05` and `neos-953928`: the filter is right and the face is large" in
`design-notes.md`.

`air05` was indeed a filter discarding cuts silently — all 217 of them, every one exactly
6557 terms of 7195 columns — and this time the filter is right: one such cut is an eighth
of the matrix, two hundred of them will not finish an LP solve in 200 seconds, and taken a
few at a time they are worth less bound than the clique cuts already getting in. It needs
five hundred of bound and the families here argue about thirty.

`neos-953928` moves by nothing at any density, across every family. Large optimal face,
and mod-2 — which answered that on `n2seq36f` — finds no cut on it at all.

`neos-820879` is answered too, and the answer is no — see "`neos-820879`: the cut families
are exhausted" in `design-notes.md`. Its cut loop converges in three rounds and then no
family finds a violated cut at all, at 25114 against a needed 25468, and raising the round
limit from 10 to 200 changes nothing. Two real defects were found and reverted on the way
(the restart trigger cannot see the root's own fixing; every root budget is a share
measured from the solve's start, so a restarted pass gets none of it), because with both
fixed the restarted pass offers zero cuts at the same bound: the families are exhausted
whatever the model's size.

That lead is now closed too. `compact` is built and tested — see "Compaction, and the
question it finally answers" — and simulating the restart sequence round by round shrinks
`neos-820879` to 2638 columns, within a couple of hundred of what a reference solver
restarts on, with the bound converging to 25108: the same place the uncompacted loop
stops, and no family finding a single cut from the third round on. The exhaustion is a
property of the families, not of the representation.

The module is kept and nothing calls it. It is the prerequisite for any restart that does
pay, and the general version was measured before as "correct, and 14% slower for nothing",
so wiring it in without a reason would repeat that.

**Three explanations for `neos-820879` — restarts, sub-MIPs, compaction — have each been
built and each turned out to be downstream of the same thing: cut families that stop
finding anything at 25114 against a needed 25468.** The next attempt on it should be a cut
family, and there is no evidence here about which.

## Standing

Measured over the 140 pure binary MIPLIB instances, 60s, 16 threads:

```text
ripsolve     34        HiGHS 46     SCIP 45     CBC 30     commercial 94
```

Confirmed by a full run, `docs/baselines/binary-2026-09-10.json`, which gained
`neos-787933` and lost the coin flip on `neos-3045796-mogo`. Ahead of CBC, and twelve
behind HiGHS. Six of the 34 are instances HiGHS does not close.

`ex9` is the widest margin of the 34: 24.9 seconds against HiGHS's 3.5 and SCIP's 23.1,
with CBC not closing it at all. The instances the probing budgets were most likely to
cost -- `mitre` at 27.5s, `eil33-2` at 32.3s, `irp` at 40.7s, `nw04` at 56.2s, `air03` at
1.5s -- all still close, and `irp` has more room than the note above credits it with.

**`neos-3045796-mogo` is a coin flip and this file used to claim otherwise.** It has
closed at 38 to 41 seconds of a 60 second budget in every baseline since it was gained,
and three runs of the current build close one of three, four runs without the newest
heuristic two of four. It is marginal on its own account. The rule this file states --
quote the set that closes every time rather than one run's count -- is right, and was not
applied to `mogo` when it was written. `neos-787933` was checked that way: three runs,
three optimal.

Against the start of these sessions, **25 every run, with `nw04` never closing in any run
of any version and `irp` a coin flip**. The nine gained:

```text
mitre           presolve probing, once the propagator stopped rebuilding itself
nw04            reduced cost fixing at the root
irp             those fixings collected as the incumbent earns them
neos-1516309    the Gomory density filter
neos-1599274    the same, plus not forfeiting optimality to a prunable skipped node
n2seq36f        mod-2 cuts
neos-3226448-wkra  four budgets between the feasibility search and a point it could reach
neos-3045796-mogo  the neighbourhood the root's reduced costs would fix
ex9             probing, once it was allowed to finish
neos-787933     every implication a big-M row states, aggregated and then rewritten
```

Five of the 34 are instances HiGHS does not close: `disctom`, `eil33-2`, `decomp2`,
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

## Where the presolve work got to, and where it goes next

**`ex9` closes.** It is the first instance gained from presolve rather than from cuts,
heuristics or fixing, and the reasoning is in "What was actually stopping `ex9`" in
`design-notes.md`. The short version is that nothing needed building: probing already
reduces `ex9` from 10404 columns to 328 and the model then closes in 25 seconds against
a 60 second limit. What stopped it was `PROBE_PATIENCE`, which cut probing off inside a
57.8 million entry gap between two proofs at 63 probes of 8560. Twice in this history
`ex9` has been recorded as a model presolve could not reduce; it was a model presolve was
not allowed to finish reducing.

Two changes went with it. A column that *both* of a probe's suppositions force is forced
outright -- worth 8244 columns against 10076 on `ex9` -- and `ABANDON_RUN` stops probing
on a model whose sweeps never finish, which is what `irp`, `nw04`, `air03` and `eil33-2`
needed and what raising the patience alone would otherwise have cost them 9 to 23 seconds
each to discover.

Three things to know before pulling on this further.

**The presolve lead is otherwise closed.** HiGHS's presolve collapses `bley_xl1` to 9.7%
of its rows, `neos-780889` to 33%, `neos-633273` to 25%, `neos18` to 23% of its columns.
Handing this solver those presolved models directly converts **only `ex10`**; the others
still reach one or two nodes on a fifth of the model. Do not build presolve reductions
for that group on the strength of the size comparison -- the experiment has been run.

**`ex10` is not the second half of `ex9`.** Given the patience it wants it does close --
objective 100, two nodes -- but in **161 seconds**, of which 78 are probing and 83 are one
root relaxation. Compacting it to 16760 rows by 408 columns changes that 83 by two
seconds. Parallelising probing would buy the 78 and leave the 119. Do not start there.

## What SCIP and CBC close that HiGHS does not

A better-posed question than "what does HiGHS do", and it was worth asking. Ten of the 140
qualify; four this solver already closes. Putting the other six to SCIP's own parameters:

```text
cod105, graph20-20-1rand   symmetry      (SCIP times out without it)
neos-787933                presolve
tanglegram6                nothing in particular
```

**The symmetric pair needs a real automorphism group** -- neither has a single duplicate
column, so no cheap detection reaches it. Note also that turning SCIP's cuts *off* makes
both faster, `cod105` 57.7s to 11.5. Cuts are not what that group wants.

**`neos-787933` closes, and HiGHS does not close it.** Its relaxation was weak by
construction: 1764 big-M linking rows, LP 3.0 against an optimum of 30. Three things
between them: the conflict graph now reads all 133 implications a linking row states
rather than one, a cut family aggregates `>=` rows through those implied bounds (root
bound 9 -> 30, exact), and the model rewritten over the gates is a restriction small
enough to answer outright -- 1764 columns, a fifth of a second. Compaction is chained
into that rewriting and is the first thing in the solver to use it.

## The one place the remaining conversions are

After everything above, the sixty-second losses that are neither bound-limited nor
already answered come down to a single cause, and all four are many rows against few
columns:

```text
ex10           root relaxation 119s on 16760 x 408 after everything
bley_xl1       root relaxation 337s on 17039 x 751, 83581 iterations at 247/s
supportcase4   root relaxation does not finish in 600s
neos-633273    two nodes in sixty seconds
```

Presolve has been tested against this group and converts `ex9` alone; compaction converts
none of it; cut families do not apply, since a bound needs a relaxation to bound. **The
next conversion at this limit is a faster relaxation on a row-heavy model.**

The specific hypothesis worth testing first, and it is a hypothesis: the hyper-sparse
triangular solve was implemented and reverted on the measured grounds that `B^-1 a_q`
runs 12 to 33% dense, and that was measured on column-rich models. A basis of 17039 rows
holding at most 751 structural columns is nearly all logicals, which is the shape that
measurement did not cover. One attempt at measuring it here was wasted by counting every
`ftran` including the dense `B^-1 b` -- it read 95% everywhere, which means nothing. The
measurement has to isolate the entering column's solve.

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
