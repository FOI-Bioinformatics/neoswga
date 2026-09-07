# Optimizer cost and the dimer criterion: what changed on three real designs

This records what the September 2026 optimizer work (branch
`optimizer-cost-and-dimer-criterion`, commits `dd083c7` through `8194721`) did to
the three GC-tier designs. Six changes went in: an exact vectorised form of the
pairwise dimer relation, a bounded-memory heterodimer screen with a cheap
pre-screen, a screen computed once per pool instead of once per alternative set,
the configured `max_dimer_bp` threaded into that screen in place of a hardcoded
free-energy cutoff, a dimer rejection guard in the greedy selection, and an
incremental coverage counter for background pruning. One of them (the guard) is expected to change
which primers are selected. The others are not.

The short version: the control holds, the guard changes selection as intended,
and the acceptance condition is met in every guarded pool - but met by the
relaxation escape in every one of them, never outright. At the panel sizes these
designs actually use, the configured `max_dimer_bp` of 3 cannot be held, and the
log says so on every run. Two claims the plan was written around are not
supported, and both are recorded below.

## The invocation is part of the measurement

Every number in this document names the full command that produced it, primer
count included. That is not decoration. An earlier revision of this measurement
omitted `-n`, the runs fell back to `num_primers` in params.json, and the three
designs came out at 96, 24 and 16 primers against shipped panels of 200, 160 and
36. Every metric moves with panel size, so nothing measured that way could be
attributed to this branch. The *E. coli* worst heterodimer read 11 bp to 3 bp
under that mistake and looked like the headline result; at the correct panel
size of 160 it does not move at all. A 24-primer panel has 276 pairs where a
160-primer one has 12,720.

The counts below are the invocations that produced the shipped designs,
recovered from the last `optimize` entry of `cli_invocation` in each run's
`run_manifest.json`.

## What was measured, and against what

The `before` column is the same code path run at `0f1a560`, the commit this
branch starts from, in a `git worktree`, on the identical `step3_df.csv`
candidate pool, under the identical method and the identical `-n`. Those
before-runs reproduce the shipped designs exactly: Jaccard 1.000 against the
set-0 primers of all three backed-up `step4_improved_df.csv` files. The
`after` column is `b3a00ba`.

Inputs, all with `max_dimer_bp: 3`:

| design | genome | candidates in step3_df.csv | delivered panel | method |
|---|---|---|---|---|
| `low_saureus` | *S. aureus* | 1215 | 200 | `dominating-set` |
| `mid_ecoli` | *E. coli* | 449 | 160 | `dominating-set` |
| `high_mtb` | *M. tuberculosis* | 319 | 36 | `network` |

Note that the delivered panel is set by `-n`, not by `target_set_size` in
params.json, which reads 96, 24 and 16 respectively and does not describe any of
the three shipped designs.

Each run produces `max_sets: 5` alternative sets; the tables describe set 0,
which is the set the metrics and the summary describe. Wall times are the `real`
line from `/usr/bin/time -p` and cover the whole invocation, all five sets
included. One machine, runs taken one at a time.

The *S. aureus* pool is 1215 candidates, not the 1222 the plan text carried; the
optimizer reduces it to 972 with a background pre-filter before selection.

## Wall time

```
neoswga optimize -j runs/gc_tiers/low_saureus/params.json -m dominating-set -n 200
neoswga optimize -j runs/gc_tiers/mid_ecoli/params.json   -m dominating-set -n 160
neoswga optimize -j runs/gc_tiers/high_mtb/params.json    -m network        -n 36
neoswga optimize -j runs/gc_tiers/high_mtb/params.json    -m dominating-set -n 36
neoswga optimize -j runs/gc_tiers/low_saureus/params.json -m hybrid         -n 200
```

| design | method | -n | before (s) | after (s) | ratio |
|---|---|---|---|---|---|
| *S. aureus* | `dominating-set` | 200 | 17.80 | 15.00 | 0.84 |
| *E. coli* | `dominating-set` | 160 | 12.02 | 11.56 | 0.96 |
| *M. tuberculosis* | `network` | 36 | 124.61 | 132.82 | 1.07 |
| *M. tuberculosis* | `dominating-set` | 36 | 15.05 | 16.73 | 1.11 |
| *S. aureus* | `hybrid` | 200 | not run | 2649.43 | - |
| *S. aureus* | `hybrid` | 96 | 2499.91 | 2392.31 | 0.96 |

The four fast runs move between 16% faster and 11% slower, in both directions,
on runs of 12 to 133 seconds. That spread is run-to-run variation on this
machine, not a measured effect. Building and consulting the dimer matrix during
selection does not cost anything visible at this scale.

**The plan's headline speed claim is not met.** It was written expecting
`hybrid` to finish on this pool in time comparable to `dominating-set`. It does
not: 2649.43 s against 15.00 s at the delivered panel size of 200, a factor of
177. Nor is the change from this branch large - the size-matched before/after
pair is the 96-primer row, 2392.31 s after against 2499.91 s before, 4.3% apart,
which is inside run-to-run variation and should not be read as a speedup.

`hybrid` at `-n 200` also returns the **identical** 200-primer set that
`dominating-set` returns (Jaccard 1.000, same 0.9899 coverage, same 10.6
selectivity density, same 11 bp worst heterodimer) for that 177x. This is the
behaviour Known Issue #8 in CLAUDE.md records, measured here on the current code
at the panel size the design actually ships.

The per-set optimizer times say why. Summing the `Total runtime` lines the
hybrid optimizer prints, one per alternative set, for the 96-primer pair:

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
  pre-screen discards 60.1% of the pairs before the thermodynamic test.

## Delivered sets

Set 0 in each case, before and after at the same `-n`. `Jaccard` is between the
before and after sets of primers.

| design | method | -n | n | Jaccard | fg_coverage | selectivity_density | worst heterodimer (bp) |
|---|---|---|---|---|---|---|---|
| *S. aureus* | `dominating-set` | 200 | 200 -> 200 | 0.476 | 0.9932 -> 0.9899 | 10.4 -> 10.6 | 11 -> 11 |
| *E. coli* | `dominating-set` | 160 | 160 -> 160 | 0.808 | 0.9431 -> 0.9388 | 32.6 -> 32.2 | 11 -> 11 |
| *M. tuberculosis* | `network` | 36 | 36 -> 36 | **1.000** | 0.7560 -> 0.7560 | 628.7 -> 628.7 | 11 -> 11 |
| *M. tuberculosis* | `dominating-set` | 36 | 36 -> 36 | 0.241 | 0.9066 -> 0.8052 | 418.2 -> 336.3 | 11 -> 10 |

The *M. tuberculosis* `network` row is the control, and it is the cleanest
result here. `network` has its own greedy and never reaches the guard, so it
exercises the other five changes alone. It returns the identical 36 primers at
identical coverage and identical selectivity density. **Tasks 1, 2, 3, 4 and 6
do not change selection.**

Both `dominating-set` rows return a different set, which is what the guard is
for. The *E. coli* panel changes least (Jaccard 0.808) and the
*M. tuberculosis* one most (0.241).

### Coverage

The plan set two percentage points as the level above which a coverage fall
should be reported rather than accepted. Three of the four runs are inside it
and one is well outside:

- *S. aureus*, `dominating-set` at 200: -0.3 pp
- *E. coli*, `dominating-set` at 160: -0.4 pp
- *M. tuberculosis*, `network` at 36: 0.0 pp
- *M. tuberculosis*, `dominating-set` at 36: **-10.1 pp**, with selectivity
  density also falling from 418.2 to 336.3

On the two large panels the guard is close to free in coverage terms. On the
36-primer *M. tuberculosis* panel it is not, and it is not buying specificity in
exchange either - both coverage and selectivity density fall. That design is
also the one where the guard delivers least: 10 bp against a configured 3.

## Is `max_dimer_bp` respected in each delivered pool?

The acceptance condition binds only pools delivered by `dominating-set`,
`hybrid` or `background-aware`, the three methods that route through the guard
in `DominatingSetOptimizer.optimize_greedy`. **It is met in every such pool, and
in every case by the relaxation escape rather than outright.** No guarded pool
here reaches 3 bp.

- ***S. aureus*, `dominating-set -n 200`: 11 bp, relaxation fired repeatedly.**
  After 29 primers, no remaining candidate was both dimer-free against them and
  able to add coverage, so the constraint was lifted to admit `ATTTTCGCAAAA`.
  It was then lifted again for nearly every subsequent pick: set 0 records 171
  unscreened admissions in a 200-primer panel, and 640 across all five
  alternative sets.

  An earlier revision of this section said the relaxation "fires once in each of
  the five alternative sets" and quoted a warning reading `the delivered pool
  contains at least one pair above max_dimer_bp=3`. Both were wrong, and the
  second string no longer exists in `neoswga/`. The relaxation used to lift the
  constraint permanently on the first stall while suppressing only the log, so
  one primer was named and the rest were admitted unscreened and unrecorded.
  That was fixed in `9ed10f3`: the constraint is now restored after each
  admission, so one warning is emitted per unscreened primer.

  The accounting is now complete, and that is checkable rather than asserted. Of
  the 9172 pairs in the delivered panel that bind above `max_dimer_bp=3`, the
  number with neither primer named in the log is **zero**.

  The fix changes the record, not the panel. The delivered set is identical
  before and after it (Jaccard 1.000, identical coverage), and that is a
  property of the algorithm rather than of this pool: after a stall, every
  candidate with positive marginal coverage dimerises with the selected set,
  and since both coverage and the selected set only grow, the dimer-free
  candidates that still add coverage can never come back. So every later pick
  stalls too, and admitting one at a time reaches the same set as lifting the
  constraint permanently. Checked over 400 random pools against the pre-fix
  loop: zero divergences, with 200 runs relaxing and 132 of those making two or
  more admissions.

  It costs wall time, because the matrix is now consulted for the whole run
  rather than abandoned at the first stall: 19.68 s against about 15.5 s on
  this pool, roughly 21 percent. An independent re-run gave 21.84 s under
  higher load. Read the ratio rather than the seconds.
- ***E. coli*, `dominating-set -n 160`: 11 bp, relaxation fired.** Log line 45,
  inside set 0 (closes at line 72), after 31 primers, admitting `AGGCCGGATAAG`.
- ***M. tuberculosis*, `dominating-set -n 36`: 10 bp, relaxation fired.** Log
  line 48, inside set 0 (closes at line 51), after 26 primers, admitting
  `CGACGCCGACGA`. This is the run that answers whether the guard changes
  anything on this genome. It moves the worst heterodimer by one base pair and
  costs 10.1 points of coverage.
- ***S. aureus*, `hybrid -n 200`: 11 bp, relaxation fired.** Log line 88, inside
  set 0 (closes at line 161), after 29 primers, admitting `ATTTTCGCAAAA` - the
  same primer at the same point as the `dominating-set` run, which is expected,
  since hybrid's Stage 1 is that greedy.

So the guard holds for the first 24 to 31 primers of each panel and then cannot
hold any further. That is the honest result. It is the behaviour the relaxation
was designed for - it reports rather than silently delivering a non-conforming
pool - but it means **no delivered design in this repository currently respects
its configured `max_dimer_bp`**, and the earlier reading that *E. coli* had been
brought to 3 bp was an artefact of measuring a 24-primer panel.

Whether that is a problem depends on the panel size a user wants. A dimer-free
panel of roughly 26 to 31 primers is available from these pools today; a
dimer-free panel of 160 or 200 is not. `max_dimer_bp: 3` at 12-mer primer length forbids
any 4 bp complementary run between any pair, which is strict, and the right
follow-up is to measure the achievable panel size across `max_dimer_bp` values
rather than leaving 3 as an unexamined default.

### Two gaps, both out of scope here

**`network` is unguarded.** It has its own `optimize_greedy`
(`network_optimizer.py:824`), which applies `calculate_dimer_score` as a soft
penalty and never builds the dimer matrix. Its *M. tuberculosis* set keeps an
11 bp worst heterodimer against a configured 3, and no warning is emitted,
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

- Background pruning, reducing 50 primers to 20: 4.14 s before, 0.01 s after,
  about 414x (4.14 divided by 0.01), returning the same set at the
  same coverage.
- The incremental counter is exact, not an approximation. On a 14-primer pool
  over 267 bins across 8 successive removals, the counter's coverage fraction
  and `HybridOptimizer._calculate_coverage` differed by 0.000e+00 at every step,
  as did the loop's predicted coverage against a real rebuild of the reduced
  set.

## Reproducing

```bash
# after
neoswga optimize -j runs/gc_tiers/low_saureus/params.json -m dominating-set -n 200
neoswga optimize -j runs/gc_tiers/mid_ecoli/params.json   -m dominating-set -n 160
neoswga optimize -j runs/gc_tiers/high_mtb/params.json    -m network        -n 36

# before: same pools, same methods, same -n, at the branch base
git worktree add /tmp/optcost/base 0f1a560
# then run `python3 -m neoswga.cli_unified optimize ... -n <count>` from that
# worktree with a params.json whose data_dir points outside runs/
```

The designs as shipped on 2026-09-05 are preserved in
`runs/gc_tiers_pre_plan2_backup/`, which `.gitignore:179` excludes from version
control along with the rest of `runs/`. The before-runs above reproduce them
exactly.

`worst_heterodimer` is a top-level object in `step4_improved_df_summary.json`,
not a field of `metrics`, shaped `{"length_bp": int, "pair": [primer, primer]}`.
It is written to the file only and never logged.

Test suite at `b3a00ba`: 4100 passed, 24 skipped, 0 failed
(`pytest tests/ -n 8 -q -p no:randomly`, 213 s).
