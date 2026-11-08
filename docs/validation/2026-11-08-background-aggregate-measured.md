# `background_aggregate` on a two-host design, measured

Date: 2026-11-08. Code: `d7dfdb0` (main after #149), plus the measurement
option added to the script in the same change as this record.

## Question

Phase 7 of the genomic-diversity plan made `background_aggregate:
"worst-case"` reachable in step 2 (the background gate judges each host on its
own length) and in step 4 (the three background-measured panel limits bind per
host). The plan's verification asks for Phase 3's D3 design to be rerun under
`worst-case` and compared with `sum` on per-host density, target coverage and
panel size. This record is that comparison. Whether `worst-case` should become
the default is the decision it was meant to inform (plan, open decision 6).

## Method

`scripts/benchmarking/diversity_baseline.py --aggregate-control
'D3__wolbachia__vs__lactobacillus+drosophila'` builds the design twice from
`count-kmers` onward, once per reading, each in its own data directory under
`tests/validation/genomes/diversity_baseline/pools_aggregate_<mode>/`. Both
are built on the same code; the Phase 3 record is not used as the `sum` arm,
because it was produced by an earlier version of the optimizer, which has
changed since (#146 among others), and a difference against it would mix two
changes.

- Target: *Wolbachia* wMel. Hosts in the design: *L. plantarum* WCFS1
  (3.35 Mb) and *D. melanogaster* (144 Mb).
- Design parameters as Phase 3 (`DESIGN_PARAMS` in the script): phi29 at
  30 C, 12-mers, `max_bg_freq` 5e-6, `max_primer` 2000, `hybrid`, no
  panel limits configured.
- Sizes requested: 6, 12, 24. Set 0 of each delivered design is evaluated
  against all five *Wolbachia* strains and all three hosts with
  `reference_panel_evaluation.evaluate_reference_panel`, occupancy-weighted
  load included. hg38 is read by k-mer counts only, so its sites and
  densities are exact and its coverage is unavailable.

Command (about 17 minutes wall, largest step `filter` at 3.3 GiB peak RSS):

```
python scripts/benchmarking/diversity_baseline.py \
    --hosts lactobacillus,drosophila --counts-only-hosts human_full \
    --aggregate-control 'D3__wolbachia__vs__lactobacillus+drosophila' \
    --sizes 6,12,24 --start-free-percent 45 --min-free-disk-gb 8
python scripts/benchmarking/aggregate_control_tables.py
```

The record is `tests/validation/genomes/diversity_baseline/aggregate_control.json`
(gitignored).

## Results

### Step 2

| | `sum` | `worst-case` |
|---|---:|---:|
| after the background gate | 707,464 | 625,157 |
| after thermodynamic filter | 492,245 | 435,421 |
| after Gini | 20,649 | 17,737 |
| shortlist (`max_primer` cut) | 2,000 | 2,000 |

The per-host gate removes 82,307 more candidates (11.6% of those that passed
the pooled gate). After the `max_primer` cut the two shortlists share 1,812 of
2,000 candidates.

### Step 4, set 0

| n | shared primers | reading | wMel cov. | worst strain cov. (wBm) | *L. plantarum* density (weighted) | *Drosophila* density (weighted) | hg38 density (weighted) |
|---:|---:|---|---:|---:|---:|---:|---:|
| 6 | 6 of 6 | both | 0.6521 | 0.1630 | 203.91 (24.31) | 437.60 (32.16) | 273.27 (12.45) |
| 12 | 12 of 12 | both | 0.7334 | 0.2571 | 273.11 (19.96) | 393.36 (25.62) | 220.33 (11.40) |
| 24 | 23 of 24 | `sum` | 0.7778 | 0.3042 | 234.75 (16.78) | 289.94 (19.45) | 173.29 (11.70) |
| 24 | | `worst-case` | 0.7792 | 0.2931 | 234.75 (16.95) | 305.32 (19.72) | 175.81 (11.76) |

Density is the exact selectivity density of wMel against each host;
the weighted figure uses the occupancy-weighted load at 30 C.

At 6 and 12 primers the two readings deliver the same panel. At 24 they differ
by one primer: under `worst-case` the *Drosophila* density is 5% higher and
the hg38 density 1.5% higher, the *L. plantarum* figure is unchanged, wMel
coverage is 0.0014 higher and the worst-strain coverage (wBm) 0.011 lower.

Step 4's own `worst_host_selectivity_density` under `worst-case` (24.31,
19.96, 16.95) equals the evaluator's weighted density against *L. plantarum*
at every size, and names *L. plantarum*, the smallest host, as the worst. The
step-4 per-host figure and the independent evaluator agree.

## Interpretation

On this design, with no panel limits configured, `worst-case` changes little.
The per-host gate removes about a tenth of the candidates that pass the pooled
gate, but most of them fall below the `max_primer` cut in either case, and the
optimizer selects the same panel at two of three sizes. The one difference, at
24 primers, is a single primer and moves every figure by a few percent, in
both directions.

This does not show that the two readings are equivalent in general. Their
step-4 difference is in how the configured panel limits are read, and no limit
was configured here, so that part was not exercised. On a design with
`min_selectivity_density`, `max_background_sites` or `max_host_coverage` set,
`worst-case` would hold the panel to the worst host's figure, which on this
design is *L. plantarum* at a weighted density of 17 to 24. No reference
threshold for those limits exists (CLAUDE.md, "Panel limits"), so that arm was
not run.

The size of this claim: one target, one pair of hosts, one optimizer, three
sizes, one run each. No benefit of `worst-case` is demonstrated, and none is
ruled out. The default stays `sum`.

## Not measured

- Panel limits under `worst-case` on real genomes. Covered by synthetic tests
  only (`tests/test_per_host_figures.py`). *Measured 2026-11-08 on a second
  pair in [2026-11-08-limits-under-worst-case.md](2026-11-08-limits-under-worst-case.md):
  they do read differently. A straddling `min_selectivity_density` is
  delivered and qualified under `sum` and refused under `worst-case` at all
  three sizes, and a straddling `max_host_coverage` makes the bounded repair
  deliver a DIFFERENT panel (Jaccard 0.600) -- the first measured instance of
  this key changing what a design delivers. The conclusion above, that
  `worst-case` changes little, holds only for a run with no limits set.*
- The H3 design (hg38 as a design background) under `worst-case`. That route
  needs the host-scale pass and about 10 GiB per step; it was not run.
- C2 random panels for the `worst-case` pool.
- `background-aware` under `worst-case`, which is refused
  (`per_host_metrics.refuse_uncalibrated_host_term`).
