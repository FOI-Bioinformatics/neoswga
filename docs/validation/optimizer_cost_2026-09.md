# Optimizer cost and the dimer criterion: what changed on three real designs

This records what the September 2026 optimizer work (branch
`optimizer-cost-and-dimer-criterion`, `dd083c7` through `4cdd3a1`) did to the
three GC-tier designs. The work added an exact vectorised form of the pairwise
dimer relation, a bounded-memory heterodimer screen with a cheap pre-screen, a
screen computed once per pool instead of once per alternative set, the
configured `max_dimer_bp` threaded into that screen and into every stage of the
hybrid optimizer, a dimer rejection guard in the greedy selection, and an
incremental coverage counter for background pruning. The guard is expected to
change which primers are selected. Nothing else is.

The short version, and it is mostly a negative result:

- The control holds. The one method that never reaches the guard returns a
  byte-identical panel, so the other changes do not perturb selection.
- The guard works on the small panel and not on the large ones. On the
  36-primer *M. tuberculosis* design it cuts dimerising pairs from 72.4% to
  27.6% of the panel. On the 160- and 200-primer designs it changes essentially
  nothing: 41.2% to 41.7% and 47.3% to 46.1%.
- No delivered design reaches its configured `max_dimer_bp` of 3. The
  constraint holds for the first 26 to 31 primers of each panel and is then
  inoperative for the rest, which on the *S. aureus* panel is 171 of 200
  primers.
- The branch does not make `hybrid` faster, and `hybrid` is not comparable in
  cost to `dominating-set`. The screen it made cheaper was not the bottleneck.

## The invocation is part of the measurement

Every number here names the command that produced it, primer count included.
That is not decoration. An earlier revision of this measurement omitted `-n`,
the runs fell back to `num_primers` in params.json, and the three designs came
out at 96, 24 and 16 primers against shipped panels of 200, 160 and 36. Every
metric moves with panel size, so nothing measured that way could be attributed
to this branch. The *E. coli* worst heterodimer read 11 bp to 3 bp under that
mistake and looked like the plan's headline result; at the correct panel size of
160 it does not move at all. A 24-primer panel has 276 pairs where a 160-primer
one has 12,720.

The counts used below are the invocations that produced the shipped designs,
recovered from the last `optimize` entry of `cli_invocation` in each run's
`run_manifest.json`.

## What was measured, and against what

The `before` column is the same code path run at `0f1a560`, the commit this
branch starts from, in a `git worktree`, on the identical `step3_df.csv`
candidate pool, under the identical method and the identical `-n`. Those
before-runs reproduce the shipped designs exactly: Jaccard 1.000 against the
set-0 primers of all three backed-up `step4_improved_df.csv` files. The `after`
column is the branch tip.

The delivered sets were first measured at `b3a00ba` and re-measured at the tip
six commits later, after the fixes to the relaxation scope, the fixed-primer
screen, the `uint8` accumulation in the dimer matrix, and the Stage-1 dimer
limit landed. **All four re-runs return Jaccard 1.000 against their `b3a00ba`
counterparts**, with identical worst heterodimer, foreground coverage and
selectivity density in every case. Those six commits corrected real defects in
what the run reports and in how the constraint is applied, and they changed no
delivered panel on these three designs. The `after` numbers below are the tip
ones; the earlier measurement would have given the same table.

Inputs, all with `max_dimer_bp: 3`:

| design | genome | candidates in step3_df.csv | delivered panel | method |
|---|---|---|---|---|
| `low_saureus` | *S. aureus* | 1215 | 200 | `dominating-set` |
| `mid_ecoli` | *E. coli* | 449 | 160 | `dominating-set` |
| `high_mtb` | *M. tuberculosis* | 319 | 36 | `network` |

The delivered panel is set by `-n`, not by `target_set_size` in params.json,
which reads 96, 24 and 16 respectively and describes none of the three designs.

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
| *S. aureus* | `dominating-set` | 200 | 17.80 | 20.14 | 1.13 |
| *E. coli* | `dominating-set` | 160 | 12.02 | 10.83 | 0.90 |
| *M. tuberculosis* | `network` | 36 | 124.61 | 126.58 | 1.02 |
| *M. tuberculosis* | `dominating-set` | 36 | 15.05 | 14.49 | 0.96 |
| *S. aureus* | `hybrid` | 200 | not run | 2770 (CPU) | - |

The four fast runs move between 10% faster and 13% slower, in both directions,
on runs of 11 to 127 seconds. Building and consulting the dimer matrix for the
whole run costs something real, but how much is not settled, and the two
measurements of it disagree:

| measurement | *S. aureus* `dominating-set` cost |
|---|---|
| paired runs alternating the two versions, `-n 96` | about +21% |
| the table above, 17.80 s before against 20.14 s after, `-n 200` | +34% |
| an independent re-run of the same case under higher load | 21.84 s |

These are not reconcilable into one figure and are not averaged here. The
paired-run design controls for load by alternating the versions, which the
single before/after pair in the table does not, so +21% is the better-controlled
estimate; but it was taken at `-n 96` and the table row at `-n 200`, so they do
not measure quite the same thing either. The 21.84 s re-run sits between the
two and was taken under a loaded machine. The most defensible statement is that
the cost is somewhere around a fifth to a third of the run, that the spread is
dominated by machine load rather than by anything in the code, and that on runs
of 11 to 127 seconds it is not a practical concern at any of these values.

**The plan's headline speed claim is not met.** It was written expecting
`hybrid` to finish on this pool in time comparable to `dominating-set`. It does
not: about 2770 s against 20.14 s at the delivered panel size of 200, a factor
of roughly 138. Nor does this branch make `hybrid` faster. The size-matched
before/after pair available for `hybrid` is a 96-primer run measured earlier on
the branch, 2392.31 s after against 2499.91 s before, 4.3% apart, which is
inside run-to-run variation.

The `hybrid` figure is CPU time, not elapsed time, and the distinction matters
here. `/usr/bin/time -p` reported `real 9832.82` against `user 2770.46` and
`sys 53.42`, because that run shared the machine with a test suite and several
other optimize runs. Reporting 9832 s would overstate the cost by more than
three times. The other rows in this table were taken one at a time and their
real and CPU figures agree, so they are reported as elapsed.

**What `hybrid` buys for that time, at this panel size, is nothing.** At the tip
it returns the *identical* 200-primer set to `dominating-set`: Jaccard 1.000,
the same `fg_coverage` of 0.9899, the same 9172 dimerising pairs. This is the
behaviour Known Issue 8 in CLAUDE.md already records, and an earlier revision of
this document reported that the dimer guard had broken the tie (Jaccard 0.455).
That was measured on the 96-primer runs. At the panel size these designs
actually use, the tie is intact, and `hybrid` costs about 138 times what
`dominating-set` costs for a byte-identical answer.

The per-set optimizer times say why. Summing the `Total runtime` lines the
hybrid optimizer prints, one per alternative set, for that 96-primer pair:

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
exercises every other change alone. It returns the identical 36 primers at
identical coverage and identical selectivity density. **The non-selection
changes do not change selection.**

Coverage falls by 0.3 and 0.4 points on the two large panels, which is inside
the two-point level the plan set as a reporting threshold, and by 10.1 points on
the 36-primer *M. tuberculosis* panel, which is not. Selectivity density falls
with it there, from 418.2 to 336.3, so on that design the constraint is not
buying specificity in exchange for the coverage it costs. It is buying a real
reduction in dimerising pairs, which the next section quantifies.

## Does the guard deliver a less dimerising panel?

`worst_heterodimer` is a single worst pair and is a poor summary: one admitted
primer pins it. Counting every pair in the delivered panel above the configured
`max_dimer_bp` of 3 is the informative measure.

| design | method | -n | pairs above 3 bp, before | after |
|---|---|---|---|---|
| *S. aureus* | `dominating-set` | 200 | 9420 of 19900 (47.3%) | 9172 of 19900 (46.1%) |
| *E. coli* | `dominating-set` | 160 | 5242 of 12720 (41.2%) | 5299 of 12720 (41.7%) |
| *M. tuberculosis* | `dominating-set` | 36 | 456 of 630 (72.4%) | 174 of 630 (27.6%) |
| *M. tuberculosis* | `network` | 36 | 450 of 630 (71.4%) | 450 of 630 (71.4%) |

**The guard works on the 36-primer panel and not on the two large ones.** On
*M. tuberculosis* it removes about three fifths of the dimerising pairs. On
*S. aureus* it removes 1.2 percentage points and on *E. coli* it adds 0.5 - a
different panel of the same quality, not a better one.

The reason is in the next section: on the large panels the constraint is
inoperative over most of the selection.

## Is `max_dimer_bp` respected in each delivered pool?

No. It is not met outright in any guarded pool, and the escape hatch carries all
three.

The acceptance condition binds only pools delivered by `dominating-set`,
`hybrid` or `background-aware`, the three methods routing through the guard in
`DominatingSetOptimizer.optimize_greedy`. It is satisfied in the narrow sense
that every such pool's log carries the relaxation warning explaining why the
worst heterodimer sits above 3 bp. It is not satisfied in the sense the plan
intended.

Counted on set 0 of each run at the branch tip:

| design | method | -n | constraint holds through | admitted unscreened | pairs above 3 bp | violating pairs naming no admitted primer |
|---|---|---|---|---|---|---|
| *S. aureus* | `dominating-set` | 200 | 29 primers | 171 | 9172 of 19900 | 0 |
| *E. coli* | `dominating-set` | 160 | 31 primers | 129 | 5299 of 12720 | 0 |
| *M. tuberculosis* | `dominating-set` | 36 | 26 primers | 10 | 174 of 630 | 0 |

So the guard holds for the first 26 to 31 picks of every panel and then stalls:
no remaining candidate is both dimer-free against what is already selected and
able to add coverage. From that point each further pick is admitted unscreened,
with the constraint restored before the next one, so every unscreened admission
is named in the log. The final column confirms the list is complete - every
violating pair in every delivered panel involves at least one named primer.

That completeness is the property the escape hatch's rationale depends on: a
short list of named primers can be swapped by hand at ordering time. **171 named
primers is not a hand-swap list.** For *M. tuberculosis* at 10 admissions out of
36 it arguably is, which is the same design where the pair count actually
improves. The honest reading is that the escape hatch's justification holds for
a handful of admissions and not for the two large panels, where the constraint
is inoperative across most of the selection.

A fifth run answers the question the four above raise, and it is the clearest
result in this document. The same *M. tuberculosis* pool at `-n 16` under
`dominating-set` delivers a **completely dimer-free panel: 0 violating pairs of
120**. The guard is not weak. The constraint is simply unsatisfiable once the
requested panel outgrows what the candidate pool can support, and past that
point the guard degrades from enforcing the threshold to reporting that it
cannot hold it. Panel size, not the guard, is the variable that decides whether
a design conforms:

| design | method | -n | pairs above 3 bp | unscreened admissions |
|---|---|---|---|---|
| *M. tuberculosis* | `dominating-set` | 16 | 0 of 120 | 7 |
| *M. tuberculosis* | `dominating-set` | 36 | 174 of 630 | 10 |
| *E. coli* | `dominating-set` | 160 | 5299 of 12720 | 129 |
| *S. aureus* | `dominating-set` | 200 | 9172 of 19900 | 171 |

The first row is worth reading carefully, because it corrects an obvious
reading of the other three. The relaxation fires seven times at `-n 16` and the
delivered panel is still completely dimer-free. An unscreened admission is not a
violation. The warning fires when the guard is lifted for a pick, and says the
admitted primer *may* pair above the threshold with something already selected;
whether it does is a separate question, answered by counting pairs. Treat the
admission count as an upper bound on how much of the panel went unchecked, and
the pair count as what was actually delivered.

Either these pools are too small to yield a 160- or 200-primer panel at this
threshold, or `-m clique`, the only method that constrains dimers structurally,
is the right method for these designs. `max_dimer_bp: 3` at 12-mer primer length
forbids any 4 bp complementary run between any pair, which is strict, and the
achievable panel size across `max_dimer_bp` values has not been measured. All
three are questions for a follow-up, not for this branch.

### A note on the warning text

The relaxation warning was reworded on this branch when the relaxation was
scoped to a single pick (it previously stayed lifted for the rest of the run and
named only the first primer). Anything grepping for the old string `without the
dimer constraint` will silently find nothing on current logs and report zero
relaxations against a log full of them. The current strings are `unscreened
against the already-selected set`, `Dimer constraint lifted`, and `Dimer-aware
selection disabled`.

### Two gaps, both out of scope here

**`network` is unguarded.** It has its own `optimize_greedy`
(`network_optimizer.py:824`), which applies `calculate_dimer_score` as a soft
penalty and never builds the dimer matrix. Its *M. tuberculosis* panel keeps
450 of 630 pairs above the configured 3 bp, unchanged by this branch, and emits
no warning, because none is due - the guard is not on this path. The shipped
*M. tuberculosis* design came from `network`, so this is the method a user of
that design is actually running, and it is the one panel here that gets no
benefit at all.

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
time in the tables above, and the background pruning figure in particular
measures one loop, not an optimizer run.

- Background pruning, reducing 50 primers to 20: 4.14 s before, 0.01 s after,
  about 404x, returning the same set at the same coverage.
- The incremental counter is exact, not an approximation. On a 14-primer pool
  over 267 bins across 8 successive removals, the counter's coverage fraction
  and `HybridOptimizer._calculate_coverage` differed by 0.000e+00 at every step,
  as did the loop's predicted coverage against a real rebuild of the reduced
  set.

## Where the evidence is kept

Nothing this document cites is in version control. `.gitignore` excludes
`runs/`, so the run outputs and logs behind every number here live only on disk:

- `runs/optcost_evidence_2026-09/` holds the run directories and the 24 logs for
  the before, after and branch-tip measurements, including the relaxation
  warnings the accounting column counts.
- `runs/gc_tiers_pre_plan2_backup/` holds the three shipped designs as they stood
  before any of this work overwrote them, verified by SHA-256 at the time.
- `runs/gc_tiers/` holds the re-derived designs.

These were originally written under `/tmp`, which does not survive a reboot, and
were copied to `runs/` for that reason. If a number here is ever disputed, those
directories are the record; if they are gone, the runs have to be repeated.

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
