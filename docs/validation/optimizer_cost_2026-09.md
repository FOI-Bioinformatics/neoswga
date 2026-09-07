# Optimizer cost and the dimer criterion: what changed on three real designs

This records what the September 2026 optimizer work (branch
`optimizer-cost-and-dimer-criterion`, commits `dd083c7` through `b3a00ba`) did to
the three GC-tier designs. Six changes went in: an exact vectorised form of the
pairwise dimer relation, a bounded-memory heterodimer screen with a cheap
pre-screen, a screen computed once per pool instead of once per alternative set,
the configured `max_dimer_bp` threaded into that screen in place of a hardcoded
free-energy cutoff, a dimer rejection guard in the greedy selection, and an
incremental coverage counter for background pruning. One of them (the guard) is expected to change
which primers are selected. The others are not. Both claims are tested below,
and one headline claim the plan was written around is not supported.

## What was measured, and against what

The `before` column is the same code path run at `0f1a560`, the commit this
branch starts from, in a `git worktree` on the identical `step3_df.csv`
candidate pool. The `after` column is `b3a00ba`.

The pre-existing designs stored in `runs/gc_tiers/` are **not** used as the
before image, for a reason worth recording. Re-running `low_saureus` at
`0f1a560` from its own `params.json` under its own recorded method delivers 96
primers; the stored design holds 200, against a `target_set_size` of 96. The
same mismatch appears on the other two (stored 160 against a target of 24, and
36 against 16). Whatever produced the stored designs on 2026-09-05 is older than
this branch's base or used different arguments, so comparing against it would
attribute earlier work to this plan. Those files were nonetheless preserved
before being overwritten, in `runs/gc_tiers_pre_plan2_backup/`, and the
comparison against them is reported separately at the end.

Inputs, all with `max_dimer_bp: 3`:

| design | genome | candidates in step3_df.csv | target set size | shipped method |
|---|---|---|---|---|
| `low_saureus` | *S. aureus* | 1215 | 96 | `dominating-set` |
| `mid_ecoli` | *E. coli* | 449 | 24 | `dominating-set` |
| `high_mtb` | *M. tuberculosis* | 319 | 16 | `network` |

Each run produces `max_sets: 5` alternative sets; the tables below describe set
0, which is the set the metrics and the summary describe. Wall times are the
`real` line from `/usr/bin/time -p` and cover the whole invocation, all five
sets included. One machine, runs taken one at a time.

The *S. aureus* pool is 1215 candidates, not the 1222 the plan text carried; the
optimizer reduces it to 972 with a background pre-filter before selection.

## Wall time

| design | method | before (s) | after (s) | ratio |
|---|---|---|---|---|
| *S. aureus* | `dominating-set` | 11.58 | 13.05 | 1.13 |
| *E. coli* | `dominating-set` | 10.96 | 12.30 | 1.12 |
| *M. tuberculosis* | `network` | 27.83 | 29.02 | 1.04 |
| *M. tuberculosis* | `dominating-set` | 11.26 | 10.69 | 0.95 |
| *S. aureus* | `hybrid` | 2499.91 | 2392.31 | 0.96 |

These are the figures as measured on 2026-09-06. The 2026-09-07 fix to the
relaxation scope (see "What the log said, and what it says now" below) adds
about 21% to the two `dominating-set` S. aureus columns, because the dimer
matrix is consulted for the whole run rather than being switched off partway.
The other rows are unaffected: the relaxation never fires on the *E. coli* or
*M. tuberculosis* pools.

**The plan's headline speed claim is not met.** It was written expecting
`hybrid` to finish on this pool in time comparable to `dominating-set`. It does
not: 2392 s against 13.05 s on the same input, a factor of 183. Nor is the
change from the work in this branch large - `hybrid` is 4.3% faster than at the
branch base, which is within the range of run-to-run variation on this machine
and should not be read as a speedup at all.

The per-set optimizer times say why. Summing the `Total runtime` lines the
hybrid optimizer prints, one per alternative set:

| set | before (s) | after (s) |
|---|---|---|
| 0 | 649.54 | 637.37 |
| 1 | 504.65 | 436.48 |
| 2 | 485.97 | 409.68 |
| 3 | 428.15 | 493.85 |
| 4 | 416.51 | 400.99 |
| total | 2484.82 | 2378.37 |

The thermodynamic screen this branch made cheaper and stopped repeating is a
small part of hybrid's cost on this pool. The dominant cost is Stage-2 network
refinement, which this branch did not touch. The screen work is still doing what
it was built to do, and it is visible in the logs:

- Before, `PRE-STAGE: Thermodynamic Filtering` appears five times in one
  `max_sets: 5` run, once per alternative set, screening 972, 876, 780, 684 and
  588 primers in turn.
- After, it appears once, followed by four lines reading `PRE-STAGE: reusing the
  thermodynamic screen computed over 972 candidates`.
- The screen also reports its work, which it did not before. The old line reads
  `Checking heterodimers between 972 primers...`; the new one reads `Checking
  heterodimers: 188170 of 471906 pairs among 972 primers`, so the cheap
  pre-screen is discarding 60% of the pairs before the thermodynamic test.

The three runs under their shipped methods are 4 to 13% slower after the change
(*S. aureus* `dominating-set` +13%, *E. coli* `dominating-set` +12%,
*M. tuberculosis* `network` +4%); the extra *M. tuberculosis* `dominating-set`
comparison is 5% faster. That spread is the cost of building and consulting the
dimer matrix during selection set against ordinary run-to-run variation, and on
runs of 11 to 29 seconds it is not a practical concern either way.

## Delivered sets

Set 0 in each case. `Jaccard` is between the before and after sets of primers.

| design | method | n | Jaccard | fg_coverage | selectivity_density | worst heterodimer (bp) |
|---|---|---|---|---|---|---|
| *S. aureus* | `dominating-set` | 96 -> 96 | 0.306 | 0.9104 -> 0.8803 | 10.2 -> 10.7 | 11 -> 11 |
| *E. coli* | `dominating-set` | 24 -> 24 | 0.043 | 0.5529 -> 0.4895 | 42.3 -> 28.1 | 10 -> 3 |
| *M. tuberculosis* | `network` | 16 -> 16 | **1.000** | 0.5282 -> 0.5282 | 732.0 -> 732.0 | 9 -> 9 |
| *M. tuberculosis* | `dominating-set` | 16 -> 16 | 0.143 | 0.6918 -> 0.5570 | 412.7 -> 318.7 | 10 -> 3 |
| *S. aureus* | `hybrid` | 96 -> 96 | 0.352 | 0.9104 -> 0.9011 | 10.2 -> 10.7 | 11 -> 11 |

The *M. tuberculosis* `network` row is the control, and it is the cleanest
result here. `network` has its own greedy and never reaches the guard, so it
exercises the other five changes alone. It returns the identical 16 primers at
identical coverage and identical selectivity density. Tasks 1, 2, 3, 4 and 6 do
not change selection.

Every method that does reach the guard returns a substantially different set,
which is what the guard is for.

One further observation, not part of the acceptance condition. At the branch
base, `hybrid` and `dominating-set` returned the *identical* 96-primer set on
*S. aureus* (Jaccard 1.000, same coverage, same density) - the behaviour Known
Issue #8 in CLAUDE.md records, where hybrid costs far more for the same answer.
After the change they differ (Jaccard 0.455), because hybrid's Stage-2 network
refinement now starts from a dimer-screened Stage-1 selection. Hybrid still
costs 183 times what `dominating-set` costs, and now buys 2.1 percentage points
of foreground coverage for it (0.9011 against 0.8803).

### Coverage fell by more than the plan allowed for

The plan set two percentage points as the level above which a coverage fall
should be reported rather than accepted. Three of the four guarded runs are
above it:

- *S. aureus*, `dominating-set`: -3.0 pp
- *E. coli*, `dominating-set`: -6.3 pp
- *M. tuberculosis*, `dominating-set`: -13.5 pp
- *S. aureus*, `hybrid`: -0.9 pp

Selectivity density moves in both directions: it improves slightly on
*S. aureus* (10.2 to 10.7) and falls on the other two (42.3 to 28.1, and 412.7
to 318.7). So the constraint is not buying specificity in exchange for the
coverage it costs on those two designs; it is buying a set with no dimerising
pair in it, and paying coverage and selectivity density for that.

This is reported, not resolved. Whether a 13.5-point coverage loss on
*M. tuberculosis* is worth a dimer-free set is a design decision that depends on
the application, and `max_dimer_bp: 3` at 12-mer primer length is a strict
setting - the constraint forbids any 4 bp complementary run between any pair. A
follow-up should measure the coverage/dimer trade-off across `max_dimer_bp`
values rather than leaving 3 as an unexamined default.

## Is `max_dimer_bp` respected in each delivered pool?

The acceptance condition binds only pools delivered by `dominating-set`,
`hybrid` or `background-aware`, the three methods that route through the guard
in `DominatingSetOptimizer.optimize_greedy`. It is met in every such pool:

- ***E. coli*, `dominating-set`: met outright.** worst heterodimer 3 bp against
  a configured `max_dimer_bp` of 3, down from 10. No relaxation warning in
  set 0.
- ***M. tuberculosis*, `dominating-set`: met outright.** 3 bp, down from 10. No
  relaxation warning in set 0. This is the run that answers whether the guard
  works on this genome, and it does.
- ***S. aureus*, `dominating-set`: not met, and the escape is doing most of the
  work.** The worst heterodimer stays at 11 bp. Selection holds the constraint
  for the first 29 primers; from primer 30 to primer 96 every pick stalls, so
  every one of those 67 primers is admitted without being screened against the
  set. The delivered 96-primer panel carries **1911 pairs above `max_dimer_bp`
  of 3, out of 4560** (41.9%). At `-n 200`, the size of the shipped design, it
  is **9172 of 19900** (46.1%), with 171 unscreened admissions.
- ***S. aureus*, `hybrid`: the same, through the same Stage-1 greedy.**

So on the two designs where a set of the requested size can be built inside the
constraint, it is. On *S. aureus* it is not: a 96-primer panel is large relative
to a 972-candidate pool and the constraint cannot be held that far, so from
primer 30 onward the run is effectively unguarded.

### What the log said, and what it says now (2026-09-07)

This section previously read: "the constraint was lifted to admit
`ATTTTCGCAAAA` and the warning fired ... the same relaxation fires once in each
of the five alternative sets", and concluded that "the relaxation reports that
in the log rather than delivering a silently non-conforming pool". Both
statements described the log, and the log was wrong.

`_run_greedy_selection` cleared its local dimer matrix on the first stall and
never restored it, while the flag that gates the warning was cleared after the
first admission. So the constraint was off for the whole remainder of the run
and exactly one primer was named. On the 200-primer panel: five warnings
against 9172 violating pairs.

With the constraint restored after each admission, the counts line up:

| run | unscreened admissions | warnings | pairs above 3 bp | pairs with neither primer named |
|---|---|---|---|---|
| `-n 96`, set 0 | 67 | 67 | 1911 of 4560 | 0 |
| `-n 200`, set 0 | 171 | 171 | 9172 of 19900 | 0 |

**The delivered panel does not change.** `-n 96` returns the same 96 primers at
the same `fg_coverage` (0.8803), the same `selectivity_density` (10.7) and the
same 11 bp worst heterodimer as the pre-fix run; `-n 200` likewise (Jaccard
1.000, `fg_coverage` 0.9899). On this pool, relaxing per pick and relaxing
permanently select identically, because the stall recurs on essentially every
pick once it starts. What changes is that the log now names every primer
admitted unscreened, and every violating pair in the delivered panel involves at
least one named primer - so the list is complete, which is the property the
"swap it by hand at ordering time" rationale depends on.

It costs time, because the matrix is now consulted for the whole run instead of
being switched off at primer 29. Paired runs alternating the two versions on the
same machine, `-n 96`: 15.50 / 15.60 / 16.69 s restored against 12.59 / 12.74 /
14.18 s unrestored, about +21%. Load average was 6 to 8 throughout, so these are
inflated against the 13.05 s in the wall-time table above, which was measured on
a quiet machine; the ratio is the part to read, not the absolute values.

**This does not make the S. aureus panel fit for ordering at
`max_dimer_bp: 3`.** 171 named primers is not a hand-swap list. The honest
reading is that the escape hatch's justification holds for a handful of
admissions and not for this design, where the constraint is inoperative over
two thirds of the panel. Either the pool is too small for a 96- or 200-primer
set at this threshold, or `-m clique` (the only method that constrains dimers
structurally) is the right method here. Both are decisions for a follow-up, not
for this branch.

### Two gaps, both out of scope here

**`network` is unguarded.** It has its own `optimize_greedy`
(`network_optimizer.py:824`), which applies `calculate_dimer_score` as a soft
penalty and never builds the dimer matrix. Its *M. tuberculosis* set keeps a
9 bp worst heterodimer against a configured 3, and no warning is emitted,
because none is due - the guard is not on this path. The shipped
*M. tuberculosis* design came from `network`, so this is the method a user of
that design is actually running.

**`background-aware` can reintroduce a dimerising pair.** Its Stage 1 calls the
guarded greedy (`background_aware_optimizer.py:202`), but its Stage 3 is the
`network` refinement, which can put back a pair the screened Stage 1 excluded.
Its delivered pool therefore carries the acceptance condition without a
mechanism that guarantees it.

Neither was in this plan's scope. Both are findings for the next one.

## Component measurements

Two numbers below come from the Task 6 review, which loaded the pre-change
module alongside the new one and ran both on the same input. They are
**component measurements, not whole-run ones** - neither is visible as a wall
time in the tables above, and the background pruning figure in particular is a
measurement of one loop, not of an optimizer run.

- Background pruning, reducing 50 primers to 20: 4.14 s before, 0.01 s after.
  The `after` figure carries one significant figure, so the ratio is about 400x
  rather than the 414x the two numbers divide to. Same set, same coverage.
- The incremental counter is exact, not an approximation. On a 14-primer pool
  over 267 bins across 8 successive removals, the counter's coverage fraction
  and `HybridOptimizer._calculate_coverage` differed by 0.000e+00 at every step,
  as did the loop's predicted coverage against a real rebuild of the reduced
  set.

## The stored designs, for the record

The comparison the plan originally specified, against the designs stored in
`runs/gc_tiers/` before this task overwrote them (preserved in
`runs/gc_tiers_pre_plan2_backup/`):

| design | n | Jaccard | fg_coverage | selectivity_density | worst heterodimer (bp) |
|---|---|---|---|---|---|
| `low_saureus` | 200 -> 96 | 0.254 | 0.9932 -> 0.8803 | 10.4 -> 10.7 | 11 -> 11 |
| `mid_ecoli` | 160 -> 24 | 0.095 | 0.9431 -> 0.4895 | 33.1 -> 28.1 | 11 -> 3 |
| `high_mtb` | 36 -> 16 | 0.444 | 0.7560 -> 0.5282 | 628.7 -> 732.0 | 11 -> 9 |

These differences are dominated by set size, not by this branch. The stored
designs hold 2.1, 6.7 and 2.2 times the primers their `target_set_size` asks
for, and a larger set covers more; the coverage figures are not comparable. The
branch-base re-runs in the tables above deliver exactly the requested sizes, so
the size difference is already present at `0f1a560` and predates this work.
Reporting the numbers in this last table as the effect of this plan would be
wrong, which is why the base re-runs were made.

## Reproducing

```bash
# after
neoswga optimize -j runs/gc_tiers/low_saureus/params.json -m dominating-set
# the shipped 200-primer S. aureus design, from its run manifest:
# neoswga optimize -j runs/gc_tiers/low_saureus/params.json -m dominating-set -n 200
neoswga optimize -j runs/gc_tiers/mid_ecoli/params.json   -m dominating-set
neoswga optimize -j runs/gc_tiers/high_mtb/params.json    -m network

# before: same pools, same methods, at the branch base
git worktree add /tmp/optcost/base 0f1a560
# then run `python3 -m neoswga.cli_unified optimize` from that worktree with a
# params.json whose data_dir points outside runs/
```

`worst_heterodimer` is a top-level object in `step4_improved_df_summary.json`,
not a field of `metrics`, shaped `{"length_bp": int, "pair": [primer, primer]}`.
It is written to the file only and never logged.

Test suite at `b3a00ba`: 4100 passed, 24 skipped, 0 failed
(`pytest tests/ -n 8 -q -p no:randomly`, 213 s).
