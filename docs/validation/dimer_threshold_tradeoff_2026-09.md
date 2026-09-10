# What panel size each candidate pool supports at each max_dimer_bp

`max_dimer_bp` forbids any complementary run longer than its value between any
two primers in a delivered panel. All three shipped designs set it to 3 and none
of them meets it. `docs/validation/optimizer_cost_2026-09.md` established why:
the guard in the greedy selection works, but the constraint becomes
unsatisfiable once the requested panel outgrows what the candidate pool can
support, and past that point the guard degrades from enforcing the threshold to
reporting that it cannot hold it.

The value 3 had never been examined. This document measures, for each of the
three GC-tier pools and each `max_dimer_bp` from 3 to 8, the largest panel size
that delivers zero pairs binding above the configured threshold, and what
coverage that costs.

## The table

`largest` is the largest `-n` whose delivered set 0 contains no pair binding
above the configured threshold. `cov` is `metrics.fg_coverage` from
`step4_improved_df_summary.json`, the measured base-by-base figure. The last
three columns describe the run at the shipped panel size.

| pool | max_dimer_bp | largest conforming n | cov at that n | cov at shipped n | violating pairs at shipped n | unscreened admissions at shipped n |
|---|---|---|---|---|---|---|
| *S. aureus*, 1215 candidates, shipped 200 | 3 | 29 | 0.3751 | 0.9899 | 9172 of 19900 | 171 |
| | 4 | 83 | 0.7413 | 0.9834 | 3097 of 19900 | 117 |
| | 5 | 126 | 0.8740 | 0.9818 | 791 of 19900 | 74 |
| | 6 | 190 | 0.9701 | 0.9803 | 39 of 19900 | 10 |
| | 7 | **218** | 0.9843 | 0.9828 | **0 of 19900** | 0 |
| | 8 | (guard off) | - | 0.9932 | 55 of 19900 | 0 |
| *E. coli*, 449 candidates, shipped 160 | 3 | 31 | 0.4979 | 0.9388 | 5299 of 12720 | 129 |
| | 4 | 72 | 0.6978 | 0.9351 | 1884 of 12720 | 88 |
| | 5 | 104 | 0.7955 | 0.9287 | 848 of 12720 | 56 |
| | 6 | 132 | 0.8519 | 0.9231 | 188 of 12720 | 28 |
| | 7 | **158** | 0.9012 | 0.9056 | 10 of 12720 | 2 |
| | 8 | (guard off) | - | 0.9431 | 63 of 12720 | 0 |
| *M. tuberculosis*, 319 candidates, shipped 36 | 3 | 26 | 0.6063 | 0.8052 | 174 of 630 | 10 |
| | 4 | **55** | 0.8485 | 0.8203 | **0 of 630** | 0 |
| | 5 | 78 | 0.9248 | 0.8559 | 0 of 630 | 0 |
| | 6 | 99 | 0.9680 | 0.8719 | 0 of 630 | 0 |
| | 7 | 118 | 0.9827 | 0.8877 | 0 of 630 | 0 |
| | 8 | (guard off) | - | 0.9066 | 6 of 630 | 0 |

The `max_dimer_bp: 8` rows are not a threshold measurement. At 8 the guard is
disabled outright; see the section below.

## max_dimer_bp = 8 turns the dimer guard off

`dimer_matrix.build` represents the pairwise relation as two indicator matrices
over the `4**(max_dimer_bp + 1)` possible t-mers, and refuses above
`MAX_CODES = 4**8` (`neoswga/core/dimer_matrix.py:51`). At `max_dimer_bp: 8`
that is `4**9 = 262144` codes, so `build` raises and
`DominatingSetOptimizer._build_dimer_matrix_for_greedy` catches it:

```
WARNING: Dimer-aware selection disabled: max_dimer_bp=8 needs 262144 t-mer
codes, above the 65536 this representation allocates. Use dimer.is_dimer_fast
pairwise for a threshold this loose.
```

That warning appears in every `max_dimer_bp: 8` run in this sweep, in all three
pools, and in no run at 3 through 7. It explains the otherwise contradictory
shape of those rows: zero unscreened admissions (the guard never stalls, because
it is not running) alongside violating pairs (nothing screened them). The
"largest conforming n" the search returns at 8 - 3 for *S. aureus*, 17 for
*E. coli*, 9 for *M. tuberculosis* - is just where an unguarded panel happens to
first contain a 9 bp complementary run. It says nothing about what the pool
supports at that threshold.

**7 is the highest value of `max_dimer_bp` the optimizer can actually enforce.**
Configuring 8 or more silently drops the constraint. It is logged at WARNING and
not otherwise reported.

## What this means per genome

**The threshold decides whether the delivered panel honours the constraint. It
barely changes what the panel covers.**

Read the two coverage columns against each other. On *S. aureus* the shipped
200-primer panel covers 0.9899 of the foreground at `max_dimer_bp: 3` and 0.9828
at 7 - 0.7 percentage points apart. The panel at 3 contains 9172 dimerising
pairs of 19900 and the panel at 7 contains none. The high coverage at 3 is not
bought by the tight threshold; it is bought by abandoning it, 171 unscreened
admissions into a 200-primer panel.

The same reading holds on the other two pools, with a different sign on
*M. tuberculosis*:

- *E. coli* at 160 primers: 0.9388 at threshold 3 with 5299 violating pairs,
  0.9056 at threshold 7 with 10. Cost of conforming: 3.3 coverage points.
- *M. tuberculosis* at 36 primers: 0.8052 at threshold 3 with 174 violating
  pairs, 0.8203 at threshold 4 with none. Loosening to 4 is better on **both**
  axes here. There is no trade-off to make on this design.

What is expensive is holding threshold 3 and conforming, because that caps the
panel at 26 to 31 primers on all three pools:

| pool | conforming panel at threshold 3 | its coverage | shipped panel coverage |
|---|---|---|---|
| *S. aureus* | 29 primers | 0.3751 | 0.9899 at 200 primers |
| *E. coli* | 31 primers | 0.4979 | 0.9388 at 160 primers |
| *M. tuberculosis* | 26 primers | 0.6063 | 0.8052 at 36 primers |

A conforming *S. aureus* design at threshold 3 covers 37.5% of the genome. The
same pool at threshold 7 delivers 218 primers covering 98.4%, with no dimerising
pair at all.

At the loose end the pool, not the constraint, becomes the limit. On
*S. aureus* at threshold 7 both `-n 300` and `-n 400` deliver 255 primers, so
218 of a reachable 255 conform. Whether that ceiling moves at other thresholds
was not measured.

## Is 3 the right default?

No. On these three pools it is unreachable at any panel size a user would ship,
and the pipeline's response to that is to abandon the constraint rather than
refuse. What a user gets from `max_dimer_bp: 3` today is a panel with 27% to 46%
of its pairs dimerising and a log full of warnings, not a stricter design.

The measurement supports **7** as the default, on three grounds:

1. It is the only value at which all three pools deliver a conforming panel at
   or near the shipped size: 218 against a shipped 200 on *S. aureus*, 158
   against 160 on *E. coli*, 118 against 36 on *M. tuberculosis*.
2. It costs between 0.7 and 3.3 coverage points against the shipped panels, and
   nothing at all on *M. tuberculosis*, where it gains 8.3 points.
3. It is the highest value the guard can enforce; 8 disables it.

6 is the more conservative choice and is close behind: 190, 132 and 99 primers
conforming, and the shipped panels carry 39, 188 and 6 violating pairs rather
than 0, 10 and 0.

What this measurement does not settle is whether a 7 bp complementary run
between two 12-mers is acceptable chemistry. That is over half the primer, and
the choice between 3 and 7 is a reaction-design judgement, not a
selection-algorithm one. What the measurement does settle is that 3 is not
currently a choice: no shipped design meets it, and the designs that would meet
it cover 37% to 61% of their genomes. Whichever value is chosen, it should be
one the pools can actually support, so that the guard enforces rather than
reports.

## Method

All runs use `-m dominating-set`, including *M. tuberculosis*, whose shipped
design came from `network`. The question here is what the candidate pool
supports, and `network` cannot answer it: it has its own `optimize_greedy`
(`neoswga/core/network_optimizer.py:824`) which applies `calculate_dimer_score`
as a soft penalty, never builds the dimer matrix and never reaches the guard.
The *M. tuberculosis* rows therefore describe the pool, not the shipped design.
The shipped design's own dimer load is in `optimizer_cost_2026-09.md`: 450 of
630 pairs above 3 bp under `network`, unchanged by the guard work.

Nothing in `runs/gc_tiers/` or `runs/gc_tiers_pre_plan2_backup/` was written to.
Each measurement ran in a scratch directory holding a copy of the pool's
params.json with `data_dir` repointed and `max_dimer_bp` set, plus a copy of
`step3_df.csv`. The position HDF5 files live under `fg_prefixes`
(`runs/gc_tiers/kmers/`) rather than under `data_dir`; their mtimes were
compared before and after a trial run and none changed.

### Invocations

Each row of the table comes from runs of this form, with `<dir>` the scratch
copy for that (pool, threshold) pair and `<n>` the panel size:

```
python3 -m neoswga.cli_unified optimize -j <dir>/params.json -m dominating-set -n <n>
```

The panel size is always given explicitly. Omitting `-n` falls back to
`num_primers` in params.json, which reads 96, 24 and 16 on these three pools and
describes none of the three designs; an earlier measurement in this series made
that mistake and produced a clean, plausible, wrong result.

Panel sizes measured per pool and threshold, in the order the binary search
visited them:

```
low_saureus  d3: 200 100 50 25 37 31 28 29 30
low_saureus  d4: 200 100 50 75 87 81 84 82 83
low_saureus  d5: 200 100 150 125 137 131 128 126 127
low_saureus  d6: 200 100 150 175 187 193 190 191
low_saureus  d7: 200 400 300 250 225 212 218 221 219
low_saureus  d8: 200 100 50 25 13 7 4 2 3
mid_ecoli    d3: 160 80 40 20 30 35 32 31
mid_ecoli    d4: 160 80 40 60 70 75 72 73
mid_ecoli    d5: 160 80 120 100 110 105 102 103 104
mid_ecoli    d6: 160 80 120 140 130 135 132 133
mid_ecoli    d7: 160 80 120 140 150 155 157 158 159
mid_ecoli    d8: 160 80 40 20 10 15 17 18
high_mtb     d3: 36 18 27 22 24 25 26
high_mtb     d4: 36 72 54 63 58 56 55
high_mtb     d5: 36 72 144 108 90 81 76 78 79
high_mtb     d6: 36 72 144 108 90 99 103 101 100
high_mtb     d7: 36 72 144 108 126 117 121 119 118
high_mtb     d8: 36 18 9 13 11 10
```

Each search starts at the shipped size. Where that size already conforms - all
of `high_mtb` at 4 and above, and `low_saureus` at 7 - the search continues
upward by doubling until it fails, then bisects; those are the runs at 72, 144,
400 and so on. Where it does not conform the search bisects downward from the
shipped size. `low_saureus` at 7 additionally shows the pool ceiling: `-n 400`
and `-n 300` both deliver 255 primers, so 218 of a possible 255 conform.

The monotonicity scans below add 12 further sizes on `high_mtb_d3` (16, 20, 21,
23, 28 to 35) and 16 on `mid_ecoli_d3` (25 to 29, 33, 34, 36 to 39, 41 to 45).

178 distinct optimize runs, 2744.6 s of run time in total (the sum of the `real`
line from `/usr/bin/time -p` over all of them; several ran concurrently, so
elapsed wall time was less). Mean 15.4 s a run.

### Counting

Violating pairs are counted directly over set 0 of `step4_improved_df.csv` with
`neoswga.core.dimer.is_dimer_fast(a, b, max_dimer_bp)`, which returns True when
binding *exceeds* its third argument, so the configured value is passed and not
one more. `worst_heterodimer` in the summary is a single worst pair and is not
used here.

Unscreened admissions are counted by the log line containing `unscreened against
the already-selected`, which is the current wording. The older string `without
the dimer constraint` no longer exists; grepping for it returns zero, which
reads exactly like a clean run.

The count must also be scoped to set 0. Each of the five alternative sets
restarts the optimizer, and a whole-log count gives 640, 229 and 96 admissions
at the three shipped sizes where set 0 alone has 171, 129 and 10 - up to 3.7
times too high. Set 0 is the log region from the first `Dominating set
optimization:` line to the third, since that line is emitted twice per set.

### Harness validation

Before measuring anything new, the harness was checked against every number in
`optimizer_cost_2026-09.md` it could reproduce. All eleven agree exactly:

| case | prior document | here |
|---|---|---|
| mtb, d3, `-n 16`, violating pairs | 0 of 120 | 0 of 120 |
| mtb, d3, `-n 16`, unscreened admissions | 7 | 7 |
| mtb, d3, `-n 36`, violating pairs | 174 of 630 | 174 of 630 |
| mtb, d3, `-n 36`, fg_coverage | 0.8052 | 0.8052 |
| ecoli, d3, `-n 160`, violating pairs | 5299 of 12720 | 5299 of 12720 |
| ecoli, d3, `-n 160`, fg_coverage | 0.9388 | 0.9388 |
| saureus, d3, `-n 200`, violating pairs | 9172 of 19900 | 9172 of 19900 |
| saureus, d3, `-n 200`, fg_coverage | 0.9899 | 0.9899 |
| set-0 admissions at shipped size, three pools | 10 / 129 / 171 | 10 / 129 / 171 |

The prior document's "constraint holds through 26 / 31 / 29 primers" figures are
also exactly the largest conforming sizes found here at threshold 3, arrived at
independently by binary search.

### Monotonicity in n

A panel at one size is not a subset of the panel at a larger size, so
conformance is not guaranteed to be monotone in `-n` and a binary search could
in principle miss a conforming size above the boundary it finds. This was
checked by exhaustive scan on two pools at threshold 3:

```
high_mtb  d3, n = 20..36: conforming through 26, violating from 27 on
                          (0,0,0,0,0,0,0, then 22,33,49,62,73,98,116,126,160,174)
mid_ecoli d3, n = 25..45: conforming through 31, violating from 32 on
                          (0 x7, then 12,40,63,86,89,94,118,124,134,159,187,211,237,251)
```

Both are cleanly monotone: a single boundary, no conforming size above it, and
the violating count rises steadily with `n`. The binary search is sound on these
pools. This is an empirical result over two pools at one threshold, not a
guarantee.

## Where the evidence is kept

The 178 run directories, logs and per-size step-4 outputs are under the session
scratch directory and are not in version control. `.gitignore` excludes `runs/`,
so nothing cited here could be committed in any case. If a number is disputed,
the invocations above reproduce it; each takes about 15 s.
