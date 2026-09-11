# Where this is, and what to do next

Written at the end of a long session so the next one does not start by re-deriving it.
The reasoning behind every claim here is in `design-notes.md`; this is the index.

## The automorphism group is built, and it is answered

**Do not start here. This was the previous file's "next build" and it is done, measured
and closed.** `crates/ripsolve/src/symmetry.rs` finds generators of the automorphism
group of the coefficient matrix, in pure Rust -- no dependency was taken and the README's
"nothing links against another solver" still holds. `Options::symmetry` uses them when
branching. It is **off by default because it converts nothing**, and the full reasoning
is under "Symmetry" in `design-notes.md`.

The detection works and is not the problem:

```text
cod105              11 generators, all 1024 columns in one orbit, 0.07s
graph20-20-1rand    36 generators, 58 orbits, widest 37, 0.18s
acc-tight2/4         7 generators, 90 orbits, widest 56, 0.02s
drayage-25-23       78 generators, widest orbit 1600, 0.58s
```

33 of the 140 pure binary instances have an orbit at least four columns wide, so this is
not a two-instance phenomenon either.

**What it ran into is that the models with the widest orbits never reach a second node.**
`cod105`, whose entire column set is one orbit, reaches one node in sixty seconds and
proves no bound. So do `neos-953928`, `neos-960392`, `neos-957143`, `neos-780889`,
`ex10`, `cvs16r106-72`, `cvs16r128-89`, `cvs16r70-62` and `drayage-25-23`. Eight of the
nine instances that have symmetry *and* that HiGHS closes and this does not are in that
state. There is no symmetric subtree to collapse because there is no subtree.

Over 33 symmetric instances at sixty seconds and one thread: nothing converted, nothing
lost, and a detection cost of 0.1 to 2.6 seconds on every model. Where a tree does exist
the mechanism does what it should -- `cvs16r89-60` reaches the same bound in 13 nodes
against 193 -- which is what the option is kept for.

**No full benchmark was run for this, and none was needed.** With the option off the
search takes the same path it took before -- detection is not called, and the node depth
that the open set's tie-break reads is the bound-change count it always was, because only
an orbital branch appends more than one. The standing below is therefore still
`binary-2026-09-11.json` and still means what it said.

Three things it rules out, so they are not attempted again: a better automorphism search
(the group is found in a tenth of a second and finding more of it changes nothing),
orbitopal fixing or a Schreier-Sims tower on top of it (stronger versions of a reduction
with nothing to reduce), and the idea that `cod105` and `graph20-20-1rand` are symmetry
instances *for this solver*. They are root relaxation instances. SCIP needs symmetry for
them because SCIP gets past the root.

Two faults were found and fixed on the way and are worth knowing about, both in
`design-notes.md`: the textbook orbital branching rule overrules pseudocost branching and
loses `neos-1599274` and `irp`, and taking the orbit *of the column pseudocosts chose* is
strictly better; and a group whose widest orbit is below the usable floor still cost
`irp` its solve by being re-examined at every node, so it is now dropped once at the root.

## The root relaxation was profiled, and it is not what this file predicted

**`Lu::factor` is the largest single item on every model that spends sixty seconds inside
its first relaxation.** 92% of `cod105`, 47% of `bley_xl1`, 32 to 40% of `ex10`,
`neos-953928`, `neos-633273` and `cvs16r106-72`. The triangular solves this file has
pointed at for three sessions are 7 to 11%.

Inside it was a quadratic: the fill-in test searched the column it was updating instead of
asking a marker, so it cost `|pivot column| * |updated column|` per column touched. On
`cod105` that was **182 billion operations against 822 million** in the elimination it
belongs to. Fixed, and the same 8000 iterations go from 57.8 seconds to 9.0. See "The
fill-in test was quadratic" in `design-notes.md`.

**It converts nothing at sixty seconds and it is in anyway**, on the same grounds the root
cut loop's warm start was: it is a real fault in the hottest function in the profile.
`cod105` still reaches one node.

### And what that leaves is the iteration count

Measured against HiGHS on the same relaxations rather than inferred, through `highspy`
with every column relaxed to continuous and one thread:

```text
                  HiGHS iters   HiGHS s      here
cod105                   4090      7.47   >100,000, does not converge
cvs16r70-62              7636      0.88   >100,000, does not converge
cvs16r106-72            12650      1.66   >100,000, does not converge
cvs16r128-89            14516      2.34   >100,000, does not converge
ex10                    25974     71.17   >100,000, does not converge
neos-953928             12710      5.48       10,751 in 16.5s
```

**These are two different gaps.** On `cod105` and the `cvs16r*` family it is the iteration
*count*, eight to twenty-four times too many -- and per iteration this solver is now the
faster of the two there, 823 a second against HiGHS's 547. On `neos-953928` it is the
reverse: the counts agree within a fifth and the time is three times worse.

The count looked like a Dantzig pricing symptom -- `Solver::price` scanned every column
and took the largest reduced cost, where the dual method already carried steepest edge
weights. It was half that and half something else, and the something else was larger.

### It was Bland's rule, and both fixes are in

Devex was built and on its own it changed nothing. So the iterations were counted by what
produced them, with `cargo run --release --example pivotmix`:

```text
                 under Bland   zero-step
    cod105             89.8%       93.9%
    cvs16r70-62        89.1%       99.0%
    cvs16r106-72       87.3%       98.7%
```

**Nine iterations in ten ran under Bland's rule**, entered after 100 degenerate steps in a
row. 100 in a row is not a cycle -- `cvs16r70-62` is 94% zero-step *while converging
normally*. The primal's threshold is now 2000 (the dual's is untouched at 100). With Bland
out of the way devex looked worth having, which it was not before:

```text
                    before             now      HiGHS
    cod105          no convergence     3801      4090
    cvs16r70-62     no convergence     7180      7636
    cvs16r106-72    no convergence     9914     12650
    neos-953928     10751 / 16.5s      7827     12710
```

Three of those are now *below* HiGHS's iteration count.

**It was benchmarked, and half of it was wrong. Devex is now off by default.** The full
run the previous session asked for scored **32 against the standing 34**: no gain, and two
instances lost. Both reproduce three runs of three, and the parent commit closes both, so
this was not variance.

```text
    neos-787933     optimal 38.4s  ->  limit
    neos-913984     optimal 51.3s  ->  limit
```

The commit changed two things at once, so they were separated, and the whole cost is on
devex. The threshold is exonerated:

```text
    pricing    stall limit        neos-787933      neos-913984
    Dantzig    100 (parent)       optimal 40s      optimal 47-50s
    Dantzig    2000               optimal 40-42s   optimal 48.7s
    devex      100                limit            optimal, 1 run of 2
    devex      2000 (as shipped)  limit            limit
```

**The case for devex was iteration counts on relaxations that do not converge, and none of
those instances close at sixty seconds either way.** `cod105`, `cvs16r70-62`,
`cvs16r106-72` and `cvs16r128-89` are all still `limit` in both runs. That is a real gain
against HiGHS's iteration count and it converts nothing, while the cost is two instances
that were closing. `Tolerances::primal_devex` is kept and defaults to `false`; the reasoning
is under "Devex is off, and what the measurement was worth" in `design-notes.md`.

**The stall threshold is the half worth keeping, and it is the whole of the improvement.**
Re-run with devex off: 34 again, none lost, none gained, and three instances materially
faster.

```text
    nw04             47.8s -> 37.5s
    neos-1599274     25.3s -> 15.1s
    irp              55.4s -> 47.8s
```

`irp` is the one to notice. This file has called it borderline in every run recorded, and
it now closes with twelve seconds to spare rather than five. `nw04` has never closed in any
run of any version before reduced-cost fixing was added. `mitre` and `decomp1`, which devex
had cost 4.7s and 8.5s, are back at baseline.

`docs/baselines/binary-2026-09-13.json` is that run.

**A benchmark artifact worth knowing about, because it will happen again.** Two instances
reported impossible times -- `neos-831188` at 240.3s and `pb-fit2d` at 15122.1s against a
sixty second limit -- and the machine suspending mid-run is why, confirmed in the journal at
21:21:28. The harness times a solve with `time.time()`, which counts suspended time, while
both `subprocess.run`'s timeout and the solver's own `Instant::now()` deadline use
CLOCK_MONOTONIC, which does not. **Statuses are unaffected** -- the solver still gets its
real sixty seconds -- but the reported seconds are not, and a contaminated entry is cached
like any other. `neos-831188` re-ran at 60.09s standalone. If a time looks impossible, check
`journalctl | grep -i suspend` before believing it, and drop that one key from
`bench/out/tractable_screen.json` rather than re-running the set.

**`pb-fit2d` is a real time-limit overrun underneath that artifact, and it is not new.** Its
suspend-free time is **458.8s against a sixty second limit**, no bound proven, which is why
it crosses the harness's `LIMIT * 4 + 120` cap and lands in the 09-13 baseline as `killed`
where 09-12 has `limit` at 294.2s. Neither counts as solved and the standing is untouched
either way. It is not this change: the parent build was still inside the same instance after
twenty minutes when this was written. Something in that model defeats the deadline check --
the checks are every 256 iterations in each simplex loop, so the time is going somewhere
that is not those loops. Nobody has looked; `pb-fit2d` appears nowhere else in these notes.

## Node propagation is built, and it converts `neos-3045796-mogo`

`crates/ripsolve/src/propagate.rs`. Two rules run to a fixpoint from the columns a branch
just moved: the conflict graph (a fixed binary excludes every literal it is joined to) and
row activity (a row's reachable range bounds every column in it, rounding inwards on
integers). A node it proves empty needs no LP at all. It is bounded by work rather than by
reaching a fixpoint -- `NODE_PROPAGATION_BUDGET`, in matrix entries -- because it runs on
every node, and stopping early is sound.

**What it converts.** `neos-3045796-mogo`, the coin flip this file said not to spend more
on without a plan for the primal side:

```text
    propagation OFF   optimal 1 of 5   (47.5s)
    propagation ON    optimal 5 of 5   (23.2 - 26.4s)
```

That is the rule this file states -- quote the set that closes every run -- and it passes
it. Over a 20 instance A/B, both arms back to back, **nothing was lost** and the instances
that already closed mostly got faster: `irp` 50.8s to 41.4 on fewer nodes, `nw04` 53.7 to
50.3, `acc-tight2` 44.5 to 41.2, `neos-787933` 44.3 to 40.8. On instances that do not close
it buys nodes: `chromaticindex32-8` 6652 to 13739, `neos18` 350 to 486, `air05` 948 to 1282.
`docs/baselines/node-propagation-ab-2026-09-09.txt`.

**The full 140 run scored 33 and it is not usable.** It reported losing `ex9` and
`neos-913984` -- and both fail **three runs of three with propagation off as well as on**.
Thirty of its thirty-two closing instances were slower than the baseline. The machine could
not reproduce the 09-14 conditions: under sustained sixteen-thread load this laptop settles
at 1.87 to 2.26 GHz against a 2.5 GHz base, where a short run gets turbo. `ex9` closed at
24.5s earlier the same day and at 23.6/24.4s in the A/B an hour before that.

**So the standing is not re-measured and 34 still stands from `binary-2026-09-14.json`,
with `mogo` a verified gain on top of it.** Before quoting a new number, get a full run
where the closing instances are not uniformly slower than the baseline -- that is the
signal that the machine is comparable, and it is worth checking on the first twenty rows
rather than at the end.

## The conflict-driven line is finished. Do not extend it; the next lever is elsewhere.

Five mechanisms built over this and the previous session, all sound, all tested, all
measured. **One converts.**

```text
    node propagation     neos-3045796-mogo, 1 of 5 runs -> 5 of 5          ON
    conflict analysis    no conflicts to work from at 16 threads           off
    Farkas certificates  196 clauses on chromaticindex32-8, no change      off
    root probing         20000 clauses on acc-tight5, no change            off
    rapid learning dive  19915 clauses, zero gained and zero lost          off
```

A/Bs on disk, not to be re-derived: `node-propagation-ab-2026-09-09.txt`,
`conflict-analysis-ab-2026-09-10.txt`, `root-probing-ab-2026-09-10.txt`,
`rapid-learning-ab-2026-09-11.txt`, and the unbounded dive as a cautionary one.

**The closing argument, measured rather than argued.** The reference solver's
`rapidlearning` returns 200 conflict constraints *and* 620 global domain reductions on
`chromaticindex32-8`. We now produce the clauses and they do not move the bound --
`examples/divebound` puts the best 500 into the relaxation as rows and the root stays at 3
against its 4. So it is the fixings that matter, and `examples/cutoffprobe` shows where
they cannot come from: probing every column against the objective cutoff decides **zero**,
with or without the cutoff row. **Propagation cannot derive them** -- "objective at most 3
and this column is one is infeasible" is an LP statement.

So anything further here needs a sub-search that *solves LPs* -- the solver run recursively
under a node limit on the cutoff model, importing the bounds it proves. That is a large
build and it is the only remaining form of this idea. Weigh it against the rest of the
board before starting, because five mechanisms on this axis have now returned one instance.

**Where the board actually stands.** The standing is 34 of 140, held across
`binary-2026-09-12/13/14` and unchanged by any of the above except `mogo`. The instances
two open solvers close and this one does not are nine, and their failure modes are in
"Where the nine instances two open solvers close actually fail". Of those, the ones nothing
here has touched are the two whose **root relaxation never finishes** (`supportcase4`,
`neos-780889`) and `neos-633273`, whose gap is presolve: **10648 of its 21781 rows are the
identical doubleton `x_p - x_q = 0`**, `aggregate.rs` already removes them, and it is
neutral because the model that remains is still too hard. Those are different problems from
this one and none of them wants another clause.

## The whole conflict-driven line is built and none of it converts

Four mechanisms, in the order they were built, each aimed at the instances two open solvers
close and this one does not. All are correct, tested and measured. **One converts.**

```text
    node propagation     converts neos-3045796-mogo, 1 of 5 runs -> 5 of 5.   ON
    conflict analysis    fires only on propagation conflicts, of which there   off
                         are none: zero over eleven instances.
    Farkas certificates  gives it the conflicts -- 196 clauses on              off
                         chromaticindex32-8 -- and changes no outcome.
    root probing         20000 clauses on acc-tight5 in 0.09s, 13279 on        off
                         neos18, and changes no outcome.
```

The A/Bs are on disk and should not be re-derived:
`docs/baselines/node-propagation-ab-2026-09-09.txt`,
`conflict-analysis-ab-2026-09-10.txt`, `root-probing-ab-2026-09-10.txt`.

**The finding that matters most is the last one, because it closes off the explanation the
other three rested on.** Root probing learns **zero clauses on `chromaticindex32-8`**. That
model is pure set packing, so every implication is already an edge of the conflict graph --
which holds **12642 edges against the reference solver's 8438**. We are not short of
implications there, and we never were. Whatever gives that solver 56x on it, it is not
information about which pairs of columns exclude one another.

So do not build: more propagation strength, a stronger clause learner, a better conflict
source, or a richer implication graph. Each has now been tried and measured.

**What is left unexplained, precisely.** That solver closes `chromaticindex32-8` in **one
node** with a bound of 4, having reduced nothing and separated no cuts. We hold the same
clique information, the optimum, and a root bound of 3, and 12000 nodes will not lift it.
The remaining difference has to be in what it does *inside* that one node -- some form of
repeated propagate-and-learn over a dive or a probing tree that raises the root bound
itself, rather than anything the tree search does. That is the thing to go and read, and
`conflict/enable` being worth 56x says the learning is part of it.

**Ablate the reference solver before building anything.** Three of these four builds were
aimed one level shallower than the answer, and each time ten minutes of
`m.setParam('conflict/enable', False)` and its neighbours would have said so. The table is
under "What the other solvers do that this one does not".

## What propagation did NOT do, which was the reason it was built

The argument for building it was `chromaticindex32-8`: a reference solver reduces nothing
on it, separates no cuts, and closes it in 0.71s where this spends sixty seconds with the
optimum already in hand. **Propagation doubles its node throughput and the bound stays at
3.** The prediction was wrong.

What is still unexplained is how that solver *proves* the bound there with no presolve and
no cuts. The candidate not yet accounted for is **conflict analysis** -- learning a clause
from each infeasible node -- which is the standard way a solver proves a combinatorial
statement like "no 3-colouring exists", and which this solver does not have. That is the
next thing to test on that instance, and it is a different mechanism from more propagation
strength.

## Doubleton aggregation is built, measured, and neutral

`crates/ripsolve/src/aggregate.rs` substitutes out equality rows on two columns, removing
a row *and* a column, and carries the first real postsolve in this solver. It is wired
around everything in `solve_from`, behind `Options::aggregate` and the CLI's `--aggregate`,
**off by default**.

**Measured properly and it converts nothing and loses nothing.** Aggregation touches 26 of
the 140 pure binary instances; the other 114 have no doubleton equality and `aggregate`
returns `None`, so they cannot be affected. Both arms run back to back on mains in
`docs/baselines/aggregate-ab-2026-09-09.txt`:

```text
    26 instances, every status identical on and off
    the only consistent time difference is neos-4754521-awarau, 4.2/4.9s -> 6.6/7.1s
```

The reduction itself is large where it applies -- `rail02` loses 39086 columns,
`neos-780889` 26600, `rail01` 20791, `neos-633273` half its rows with its root relaxation
going 4.32s to 1.58 -- and it still changes no outcome. Six of those (`rail01`, `rail02`,
`neos-780889`, `supportcase10`, `bley_xl1`, `bab6`) reach **one node either way**: their
root relaxations do not finish, which no amount of column removal fixes.

**Ignore the earlier run that scored 33.** It was taken on battery, 1.06 to 1.19x slower
than the mains baseline it was compared against, and every instance it "lost" is one
aggregation never touches. `ex9`, the reported loss, fails three runs of three on the
committed build `6fa6233` itself and closes at 24.5s once plugged in. **Check
`cat /sys/class/power_supply/AC/online` before trusting a timing.**

## The density cap is NOT what stops aggregation paying. Three levers tested, three refuted

`neos-633273` is the instance aggregation was built for, and with cutting off it goes from
727 nodes to **1964** with a better bound. Turn cutting on and it is 2 nodes. The obvious
explanation was the Gomory density cap: aggregation deletes the model's shortest rows, the
mean row length goes 2.93 to 3.74, and under integer division the cap widens from 48 terms
to 72. **That explanation is wrong, and so are the two after it.**

**The cap level is not the lever.** Forcing it, at fifty rounds:

```text
    cap    aggregation ON        aggregation OFF
     16    1894 nodes (0 cuts)   482 nodes (0 cuts)
     30       2 nodes (645)      117 nodes (665)
     48       2 nodes (774)      111 nodes (757)
     72       2 nodes (803)        2 nodes (968)
    120       2 nodes (977)        2 nodes (1020)
```

With aggregation on, **every cap that admits any cut at all collapses to two nodes**. Only
a cap so tight it admits none (16) leaves a search. The same ~650 cuts at the same cap
leave 117 nodes on the unaggregated model. So it is not the density of the cuts.

**The ageing trigger is not the lever either.** `active.len() > base.n_rows()` never fires
on a large model, which is a real defect, but forcing it to fire changes nothing here:
dividing the threshold by 20, 100 and 400 gives 2 nodes at every setting, down to 190
cuts, because each purge it triggers costs the cold rebuild the loop otherwise avoids.

**Column density explains the wrong instance.** Aggregation merges columns, so a survivor
carries the union of both columns' entries -- but measured, `neos-633273`'s columns get
*sparser*, mean 5.7 to 3.9 with the median and max unmoved. It does explain
`neos-4754521-awarau`, whose columns double (mean 29.7 to 58.3, median 12 to 26, max 9801
to 19600) and which is the one instance consistently ~2s slower with aggregation on.

**So what is left is the question nobody has answered: why do a few hundred cuts of
bounded density cost twenty seconds a node on the aggregated `neos-633273` and a fraction
of a second on the unaggregated one?** Phase timing says root 1.58s, five rounds 1.3s,
purge 12s -- about 15 seconds of 60 -- and the remaining 45 go into *two* node LPs. That
is the thing to profile, and `examples/aggsize.rs` now reports column statistics to help.
Do not spend another session tuning `DENSITY_FACTOR`; it has been swept.

## The root cut loop's guard was taking the run, and it is fixed## The root cut loop's guard was taking the run, and it is fixed

**A third instance of the fault this file has now recorded three times: a guard given the
caller's whole remaining clock.** The root cut loop is budgeted by `ROOT_CUT_SHARE`; the
purge that runs after it was not. Phase timing on `neos-633273`:

```text
    separate 0.10s -> resolve 0.08s warm     (x5 rounds, ~1.5s in total)
    final-purge cold 54.43s, 11901 iterations
```

The loop costs a second and a half. Its guard took **54.4 of a 60 second budget**, leaving
two nodes and no incumbent where the same model reaches 732 nodes with no cutting at all.
`ROOT_PURGE_SHARE` now bounds it at a fifth, chosen by measurement (the table is on the
constant). Bounding it is safe by construction: a purge that does not finish leaves the
untrimmed model, which already solves and already holds the bound.

Measured over the 140 at 60s and 16 threads: **34, none gained, none lost**, with `irp`
47.8s to 41.1, `cap6000` 18.3 to 15.4, `neos-787933` 41.3 to 39.2, `acc-tight2` 38.7 to
37.3, against `nw04` 3 to 4 seconds slower over three runs each way.

**Read this before touching a deadline anywhere.** Writing the fix introduced a bug the
whole test suite passed with: the shortened deadline was set on `candidate`, which is
conditionally promoted to the LP the *entire search* runs on, so every node inherited it
and gave up at a fifth of the run -- `p2m2p1m1p0n100` reported the time limit after 12.1
seconds of 60. 213 tests did not see it, and neither did a targeted A/B over ten models.
The full benchmark's partial diff caught it at 31 instances. If a bounded deadline is set
on an `Lp` that outlives the solve it was for, restore it immediately afterwards.

## Where the nine instances two open solvers close actually fail, and the two defects left

The set this solver misses that at least two of HiGHS, SCIP and CBC close is nine. Two
were already answered (`air05`, `neos-3045796-mogo`), two never finish their root
relaxation (`supportcase4`, `neos-780889`), and the full table is in `design-notes.md`
under "Where the nine instances two open solvers close actually fail".

Two things were established there that the next session should not re-derive:

**The Gomory density filter is not a defect, and loosening it is a trap.**
`chromaticindex32-8` separates one cut where its neighbours separate hundreds, and the
filter explains it -- rows averaging 4 columns of 2304 give an allowance of 72 terms and
the one surviving cut is exactly 72. Disabling it admits 1459 more cuts and the bound does
not move at all, while `neos18`'s bound gets *worse*, 8.5 to 7.6, because dense cuts crowd
out the sparse ones that pay. This was tested, not assumed.

**`neos-633273` is a presolve instance, and it is the argument for building postsolve.**
SCIP closes it in 0.47s, reducing 11154 x 21781 to 5079 x 5498; our presolve reports zero
of everything. **10648 of its 21781 rows are the identical doubleton equality
`x_p - x_q = 0`** on binaries. Merging equal columns removes a row and a column each time
and cannot be expressed in place, because `x_p` has to be recovered from `x_q` --
`presolve.rs` is deliberately postsolve-free and says so. Doubleton aggregation is the
first reduction that needs postsolve, and this is the instance that justifies it.

A second, smaller defect is diagnosed and unfixed: **the in-loop cut ageing never fires on
a large model**, because its trigger is `active.len() > base.n_rows()` -- on a 21781-row
model the pool would need 21781 cuts before a slack one ages out. That is why fifty rounds
accumulate 969 cuts and the node LPs crawl even with the purge bounded: `neos-633273`
reaches 290 nodes at five rounds and still 2 at the default fifty. It needs its own A/B.

## The order after that, and this is the second measurement saying so

The order this file used to state was LP throughput, then symmetry, then `mogo`'s primal
side. LP throughput was closed as a *conversion* lever and symmetry is now closed as one
too, and **both closed onto the same thing**: a group of instances that spend sixty
seconds inside their first relaxation.

**That hypothesis has now been tested and is closed.** The entering column's solve was
isolated, as the note demanded, and being row-heavy does not make it sparse: `ex10` 79.0%
dense, `supportcase4` 70.2%, `bley_xl1` 36.1%, all with under 4% of solves under the
threshold. The only hyper-sparse models are `neos-633273` and `neos-555001`, both of which
already solve their relaxation. See "The hyper-sparse triangular solve is refuted" in
`design-notes.md`. Do not build it.

What is new since that was written is the size of the group it would serve. It is not the
four models named there. Adding the symmetric instances that stop at one node, the list
of models whose sixty seconds go into a relaxation that does not finish is at least:

```text
cod105  ex10  neos-953928  neos-960392  neos-957143  neos-780889  cvs16r106-72
cvs16r128-89  cvs16r70-62  drayage-25-23  drayage-100-23  chromaticindex128-5
bley_xl1  supportcase4  neos-633273  air04  tanglegram6
```

Six of those are instances HiGHS closes and this does not.

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

`docs/baselines/binary-2026-09-14.json`, the run that validated the cut-purge bound below:
the same 34 instances as `binary-2026-09-13.json` and `binary-2026-09-12.json`, none gained
and none lost across all three. The paragraph below was written
against the 09-12 run and still holds instance for instance. **34 is the
number this file has been saying the standing really is**: the run before it reported 35
and the 35th was `neos-3045796-mogo`, a coin flip at about two in five, which landed the
other way this time. Nothing else changed -- no instance gained, none lost -- and the
factorization fix shows up as speed on the instances that were already closing: `disctom`
16.2s to 13.1, `neos-1599274` 28.8 to 25.3, `acc-tight2` 40.1 to 36.8, `neos-3226448-wkra`
8.9 to 7.8, eight of the ten slowest between 0.81x and 0.96x. `irp` went the other way,
47.2 to 55.4, which is inside what sixteen threads swing by on an instance that has been
borderline in every run recorded here -- and the 09-13 run puts it back at 47.8, which
settles that: it was swing, not the factorization.

Ahead of CBC, and twelve behind HiGHS. Six of the 34 are instances HiGHS does not close.

`ex9` is the widest margin of the 34: 24.9 seconds against HiGHS's 3.5 and SCIP's 23.1,
with CBC not closing it at all. The instances the probing budgets were most likely to
cost -- `mitre` at 27.5s, `eil33-2` at 32.3s, `irp` at 40.7s, `nw04` at 56.2s, `air03` at
1.5s -- all still close, and `irp` has more room than the note above credits it with.

**`neos-3045796-mogo` is a coin flip at about two in five, and no cut budget changes
that.** Its root bound is the optimum before any cut, so it is purely primal, and what
decides it is whether the parallel search reaches -175. Three settings were measured over
six runs each -- default 2 of 6, cut pool capped at half the model's rows 3 of 6,
cut-rounds 5 **1 of 6** -- after a three-run sweep had made cut-rounds 5 look like three
of three. Do not spend more on it without a plan for the primal side.

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
column, so no cheap detection reaches it. It was built and it is not the answer *here*;
see the top of this file. `cod105`'s whole column set comes back as one orbit in 0.07
seconds and the model still reaches one node in sixty. Note also that turning SCIP's cuts
*off* makes both faster, `cod105` 57.7s to 11.5. Cuts are not what that group wants
either, which is a second sign that what it wants is a relaxation it can finish.

**`neos-787933` closes, and HiGHS does not close it.** Its relaxation was weak by
construction: 1764 big-M linking rows, LP 3.0 against an optimum of 30. Three things
between them: the conflict graph now reads all 133 implications a linking row states
rather than one, a cut family aggregates `>=` rows through those implied bounds (root
bound 9 -> 30, exact), and the model rewritten over the gates is a restriction small
enough to answer outright -- 1764 columns, a fifth of a second. Compaction is chained
into that rewriting and is the first thing in the solver to use it.

## The root cut loop was solving the model cold every round

Fixed, and it is the largest throughput change here: `decomp1` 32.2s to 1.4, `decomp2`
42.3 to 2.5, `mitre` 27.4 to 3.6, `f2gap801600` 7.0 to 1.0. Adding rows to a solved model
is a dual warm start, `BasisState::extend_for_rows` already existed for it, and
`Lp::solve_with_rows` was built on it for node-local cutting; the root loop rebuilt from
scratch and solved cold instead. Ageing a cut out forces the rebuild, so it now happens
when the pool passes the model's own row count rather than every round.

**If you are looking for more of the same, look for cold solves.** This one was in plain
sight behind a comment explaining why the rebuild was necessary, and the machinery to
avoid it was two files away with a comment saying what it was for.

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

**Since measured, and mostly closed as a lead.** Warm-starting the cut loop was found by
chasing this group and is worth far more than the group is: `neos-633273` still reaches
two nodes. Given 300 seconds instead of 60, *none* of `neos-633273`, `acc-tight4`,
`neos-957143` or `neos-960392` closes, so LP speed is not what they are short of. And the
hyper-sparse hypothesis is refuted -- the row-heavy models are the densest, not the
sparsest. Node warm starts on `neos-633273` take 1 to 9 pivots and 0.04 seconds against a
17 second cold solve, so the node LPs were never the problem either.

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
  `docs/baselines/binary-2026-09-12.json`, with the run before the symmetry work and the
  factorization fix kept as `binary-2026-09-11.json`, the run before the reduced cost
  neighbourhood as `binary-2026-09-06.json`, the one before the feasibility work as
  `binary-2026-09-05.json`, the one before mod-2 cuts as
  `binary-2026-09-04.json`, the one before the cut density fix as
  `binary-2026-09-03.json`, the one before this round's reduced cost work as
  `binary-2026-09-02c.json`, the one before the root work as `binary-2026-09-02.json` and
  the one before probing as `binary-2026-09-01.json`. Those filenames are a sequence
  rather than calendar dates, and the next one continues it. Diff against them rather than
  against a remembered number, and copy the current one aside before a run that matters.
  A full 140 instance run takes a little over two hours.
- `cargo run --release --example pivotmix -- <models>` reports, per relaxation, the share
  of iterations that ran under Bland's rule and the share whose step moved nothing.
  `--no-devex` prices by the largest reduced cost instead. Reach for it before a profiler
  on any model whose LP is slow: an iteration count says a solve is slow and says nothing
  about what the iterations are, and not asking cost three sessions here.
- `cargo run --release --example symmetry -- <models>` reports each model's automorphism
  group and what finding it cost, with `--seconds=`, `--stall-work=` and
  `--max-vertices=` setting the budget. Use it in preference to a solve for anything
  about symmetry: the whole local set takes minutes and has no search noise on it.
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
