# Panel limits under `worst-case`: one design, one threshold, two verdicts

Measured 8 November 2026, on `main` at `5343bc8` (plus the script this record
names). Against the first "Not measured" item of
[`2026-11-08-background-aggregate-measured.md`](2026-11-08-background-aggregate-measured.md):

> Panel limits under `worst-case` on real genomes. Covered by synthetic tests
> only (`tests/test_per_host_figures.py`).

Reproduce with:

```
PYTHONPATH=. python docs/validation/aggregate_limits_2026-11-08/measure.py \
    --sum-work <dir> --worst-work <dir> --out results.json \
    --sizes 6,12,24 --limit min_selectivity_density=100
```

Results, both pools' configurations and the arms in
[`aggregate_limits_2026-11-08/`](aggregate_limits_2026-11-08/).
**No code changed.**

## Why this run exists

Step 4's entire difference between `sum` and `worst-case` is how
`min_selectivity_density`, `max_background_sites` and `max_host_coverage` are
read: against the pooled figure, or against the per-host reduction. The
published comparison on D3 ran with **no panel limits configured** and
concluded that `worst-case` "changes little". It could not have concluded
otherwise: with no limit set, the only thing that differs is the step-2 gate.

This runs the arm that can tell them apart, on real genomes.

## Setup

*M. tuberculosis* H37Rv (4,411,532 bp, GC 65.6%, 1 record, no N) as the
target, against two real hosts **88x apart in length**: *S. aureus* NCTC 8325
(2,821,361 bp) and human **chr1** (hg38, 248,956,422 bp). A bacterial target in
a human background with a co-colonising commensal, which is also D3's
small-plus-large shape -- the case Known Issue 6 says the pooled reading gets
wrong.

Design parameters are D3's (`DESIGN_PARAMS` in
`scripts/benchmarking/diversity_baseline.py`): phi29 at 30 C, 12-mers,
`max_bg_freq` 5e-6, `max_primer` 2000, `hybrid`, sizes 6/12/24. GC bounds are
left to the pipeline rather than carried over, since D3's were for a 35%-GC
target and this one is 65.6%. Both genomes came from the AWS iGenomes mirror
and are verified by length, record count, GC and SHA-256;
`scripts/fetch_reference_genomes.py` carries the entries.

**One data directory per reading**, each built from `count-kmers`: the step-2
gate differs, so the shortlists differ and one directory cannot serve both.

## The gate difference is real, and does not reach the panel

| after | `sum` | `worst-case` |
|---|---:|---:|
| hard QC | 441,523 | 410,990 |
| Gini | 68,367 | 62,606 |
| shortlist (`max_primer`) | 2,000 | 2,000 |

The per-host gate removes **30,533 more** candidates, 6.9% of those the pooled
gate passed (D3: 11.6%). With no limits set the two readings then deliver
**byte-identical panels at all three sizes** -- same primers, same order, same
`fg_coverage`. D3 differed by one primer at 24; here, nothing. So the step-2
attrition does not reach the delivered panel on this pair either.

## The figures a threshold has to straddle

Unlimited, both readings (the pooled figures are identical by contract):

| n | pooled density | worst-host density | pooled bg sites | worst-host sites | pooled bg_cov | worst-host bg_cov |
|---:|---:|---:|---:|---:|---:|---:|
| 6 | 351.57 | **35.17** | 1 | 1 | 2.3831e-05 | 2.4101e-05 |
| 12 | 275.37 | **28.84** | 3 | 3 | 7.1492e-05 | 7.2302e-05 |
| 24 | 140.96 | **27.11** | 8 | 8 | 1.9064e-04 | 1.9280e-04 |

Three facts about this table decide what can be measured:

- **The density gap is tenfold**, because the pooled denominator is 98.9%
  chr1. A straddling threshold is easy to choose.
- **`max_background_sites` admits no straddle here.** chr1 holds every site, so
  the per-host maximum EQUALS the pooled total and the limit reads identically
  under both readings. The reduction only departs from the pooled figure when
  sites are spread across hosts. That arm was not run, for that reason.
- **The host-coverage gap is 1.1%**, for the same denominator reason. A
  straddle exists but is narrow and size-specific.

**The worst host is the host the panel never matches.** *S. aureus* holds
**zero** exact sites for every delivered panel, yet it is the worst host at
every size: its occupancy-weighted load is 9.9 per Mb against chr1's 0.89. In
the regime this pool reports -- mean matched/mismatched occupancy 1.00 with
100% of candidates saturated, Known Issue 17 -- that load is a mismatch
neighbourhood figure, not an exact-site count. **The limit that binds under
`worst-case` here is set by the occupancy model on a host with no exact
matches.**

The thresholds below are chosen from this table to straddle the two readings.
They are **not recommendations**: no reference for these limits exists
(CLAUDE.md, "Panel limits"), and a threshold picked to separate two readings is
a measurement instrument, not advice.

## 1. `min_selectivity_density: 100` -- delivered under `sum`, refused under `worst-case`

| n | `sum` | `worst-case` |
|---:|---|---|
| 6 | 6 primers, no violation, **qualified** | **`selectivity below minimum`**, repair did not resolve, panel unchanged |
| 12 | 12 primers, no violation, **qualified** | **blocking**, repair moved the panel: worst-host density 28.84 to 30.53 |
| 24 | 24 primers, no violation, **qualified** | **blocking**, repair shrank it to **21 primers**, worst-host density 27.11 to 30.59 |

The same design and the same threshold are **delivered and qualified under
`sum` at every size, and refused under `worst-case` at every size.** The
refusal reaches the export gate: the validator records `panel_limit_not_met`,
which `limit_violation_issue` describes as the only channel by which a
configured limit reaches `export`, so the panel needs `--allow-unqualified`.
Under `sum` the issue list is empty.

The artifact says which length the floor was read on. With the limit
configured, `panel_regime.criteria` gains a `selectivity_density` row
annotated *"configured floor, on the WORST HOST's own length"*, and
`panel_regime.limiting` flips from `coverage` to `selectivity_density` at
n=24.

**The bounded repair is not a no-op, and it is not free.** It optimizes the
axis the limit names and pays on others:

| n | worst-host density | fg_coverage | pooled fg/bg ratio | panel |
|---:|---|---|---|---|
| 12 | 28.84 → **30.53** | 0.2067 → 0.2293 | 4.82 → **1.88** | 12 primers |
| 24 | 27.11 → **30.59** | 0.2628 → 0.2457 | 2.47 → 4.58 | **24 → 21 primers** |

The direction on the pooled axes is not consistent: the ratio collapses at 12
and improves at 24. At 24 the shrink adds its own validator issue
(`set_size_mismatch`) and an advisory violation. Runtime goes from 69-104 s to
**295-345 s** when a `worst-case` limit is unmet, about 3-5x, which is the
repair searching.

It never reaches 100. It cannot: the binding host is one the panel does not
exactly match, so no choice of primers from this pool moves its load far.

## 2. `max_host_coverage: 1.91e-04` -- the reading changes the DELIVERED panel

The narrow straddle gives the other outcome, and it is the more interesting one.

| n | `sum` | `worst-case` |
|---:|---|---|
| 6, 12 | unchanged, no violation | unchanged, no violation (both figures below the limit) |
| 24 | pooled 1.9064e-04 ≤ limit, **unchanged panel** | worst-host 1.9280e-04 > limit, repair **resolved it** |

At 24 primers under `worst-case` the repair delivers a **different panel** --
Jaccard **0.600** against the one `sum` delivers, still 24 primers, no
violation, qualified:

| | `sum` (unchanged) | `worst-case` (repaired) |
|---|---:|---:|
| fg_coverage | 0.2628 | **0.2805** |
| worst-host bg_coverage | 1.9280e-04 | **9.6402e-05** |
| pooled bg_coverage | 1.9064e-04 | 9.5322e-05 |
| pooled density | 140.96 | 103.23 |
| worst-host density | 27.11 | **21.00** |

**This is the first measured instance of `background_aggregate` changing what a
design delivers.** The unlimited comparisons -- D3's and this record's -- could
not see it, because with no limit configured the reading has nothing to act on.

Two things follow:

- **On this arm the `worst-case` panel is better on the axes the user asked
  about**: half the host coverage and 0.018 more target coverage, at the same
  panel size. That is one arm at one size, and it is the repair's doing rather
  than the search's -- the optimizer's own output at 24 was the panel with
  twice the host coverage.
- **The three background-measured limits do not move together.** Repairing
  host coverage pushed the worst-host DENSITY down, 27.11 to 21.00. Both are
  read against per-host reductions under `worst-case`, and
  `per_host_metrics` already says each reduces independently, "because the
  host with the lowest density need not be the one with the most sites". Here
  improving one reduction measurably worsened another, in a single bounded
  repair pass. A design setting both would need the repair to satisfy both at
  once, and §1 shows the density floor is out of reach on this pool, so that
  combination was not run.

## What this changes

- **The first "Not measured" item of the 2026-11-08 record is answered on real
  genomes**: the limits do read differently, the verdict differs at every size
  for a straddling density floor, and the delivered panel differs once a
  repair can act.
- **`worst-case` is still not shown to be the better default.** Its one
  measured panel change was an improvement on two axes and a worsening on a
  third, at one size, on one pair. The default stays `sum`, as the schema's own
  reasoning requires.
- **A user reading these figures needs to know which host binds.** On this
  design it is the small host with no exact matches, through the occupancy
  model. A `min_selectivity_density` that looks generous against the pooled
  figure (100 against 141-352) is unreachable against the worst host (27-35).

## Scope of the claim

One target, one pair of hosts, one optimizer, three panel sizes, one run per
arm, 12-mers, phi29 at 30 C, two of the three background-measured limits. The
specific caveats:

- **The pool cannot discriminate** (occupancy 1.00, 100% saturated -- Known
  Issue 17), inherited from mirroring D3's phi29-at-30 C parameters. Every
  occupancy-weighted figure here is at the saturated end, and the per-host
  density is one of them. A warmer reaction would move the worst-host density,
  and might move which host is worst.
- **`max_background_sites` was not exercised**, for the reason given above: one
  host holds every site on this design, so the reduction equals the pooled
  total.
- **The thresholds are instruments, not references.** Both were read off the
  unlimited arms to straddle the two readings. Nothing here says a design
  should use either.
- **chr1 is 7.4% N** and the denominators include it, identically in every arm.
- No wet-lab outcome is involved, and there is still no BAM in this repository.
