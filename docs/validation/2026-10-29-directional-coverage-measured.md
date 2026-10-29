# The symmetric coverage window inflates `fg_coverage` by 1.3x to 1.9x

Measured 29 October 2026, on `main` at `37f8acd`. Phase 5.1 of
[the plan](../reports/PLAN_2026-10-28_pool_design_trustworthiness.md), which
required a measurement before any fix and did not pre-approve one. Reproduce
with:

```
PYTHONPATH=. python docs/validation/directional_coverage_2026-10-29/measure.py \
    --prefix <work>/target --length 1270000 --panel <work>/step4_improved_df.csv
```

## The defect, restated

`coverage._mark_window` marks `occupied[pos-extension : pos+extension]`
(`coverage.py:234`) -- a window `2r` wide centred on every site -- and both
coverage entry points default to `strand="both"` (`coverage.py:28`, `:117`), so
the position list holds the oligo's own occurrences and its reverse
complement's together. Every occurrence therefore earns coverage in the
direction it cannot extend in.
[`2026-09-22-swga-comparative-review.md`](2026-09-22-swga-comparative-review.md)
put it at its sharpest: `polymerase_extension_reach` returns an EXTENSION
distance and the code spends it as a RADIUS. That entry argued the defect from
source and explicitly declined to quantify it. This quantifies it.

The directional model is one direction per occurrence, which is what the
chemistry allows and what COATswga implements:

- the oligo occurring literally at `i` anneals to the minus strand, the nascent
  strand is plus-sense, extension runs toward INCREASING coordinates, covering
  `[i, i + r]`;
- its reverse complement occurring at `j` means the oligo anneals to the plus
  strand and extension runs toward DECREASING coordinates, covering
  `[j + len - r, j + len]`.

`position_cache` already separates the two: `strand="forward"` is the primer's
own entry and `strand="reverse"` is its reverse complement's
(`position_cache.py:369-386`). Both models are read off the SAME position
lists, so what follows is the geometry alone and not a different set of sites.

## The measurement

Three panels, each selected by the tool itself through the full four-step
pipeline, against a 1.27 Mb target at 35% GC. At the shipped 3 kb reach:

| panel | symmetric (shipped) | directional | inflation | ratio |
|---|---|---|---|---|
| 6 primers | 0.4585 | 0.2484 | **21.0 pp** | **1.85x** |
| 12 primers | 0.7218 | 0.4637 | **25.8 pp** | **1.56x** |
| 24 primers | 0.9083 | 0.7137 | **19.5 pp** | **1.27x** |

So the 12-primer panel the tool reports as covering 72% of the target covers
46% of it under the model the chemistry supports. That is not a rounding error
on the headline figure; it is the difference between a panel that looks adequate
and one that does not.

### How it varies with site density

Sweeping the reach moves the same panels across the density regimes:

| reach | n=6 ratio | n=12 ratio | n=24 ratio |
|---|---|---|---|
| 250 | 1.99 | 1.95 | 1.91 |
| 500 | 1.98 | 1.93 | 1.87 |
| 1000 | 1.94 | 1.85 | 1.72 |
| 2000 | 1.92 | 1.69 | 1.45 |
| 3000 | 1.85 | 1.56 | 1.27 |
| 6000 | 1.67 | 1.31 | 1.07 |

The ratio tends to **2.0** as windows stop overlapping -- an isolated site marks
`2r` symmetrically and `r` directionally -- and to **1.0** as the genome
saturates, where both models mark everything.

**This refines the comparative review's expectation rather than confirming it.**
That entry predicted "a small error on dense panels and a larger one on sparse
panels". True of the RATIO. But the absolute inflation in percentage points
peaks in the MIDDLE and falls off at both ends: at n=24 it rises from 7.9 pp at
250 bp to 24.5 pp at 2000 bp and back down to 6.7 pp at 6000 bp. Very sparse
panels have a large ratio on a small number, and saturated ones have no room to
be wrong. The worst absolute error is at intermediate density -- which is where
a real SWGA design sits, since that is the regime where adding primers still
buys coverage.

### The implementation was validated analytically, not trusted

The whole finding rests on the directional marking being right, so it was
checked against a case with a known answer: a single isolated site at r=1000
marks exactly **1000 bp** directionally and **2000 bp** symmetrically, for a
forward and a reverse occurrence alike. Both strands also carry substantial
site counts in every panel (forward 64-248, reverse 52-209), so the comparison
is not dominated by one strand.

## What this does NOT establish

- **The genome is synthetic**, because no real reference is reachable from this
  container: `tests/validation/genomes/` is gitignored and absent, and
  `ftp.ncbi.nlm.nih.gov` is refused by the egress proxy (403 on CONNECT), so
  `scripts/fetch_reference_genomes.py` cannot run. Random sequence gives
  Poisson-distributed site positions; real genomes cluster their k-mer sites,
  and the review's argument is that the error depends on how forward and
  reverse-complement sites interleave. **The direction and the mechanism
  generalise; the exact percentages are for this site geometry.** The panels and
  the position indexes are real -- selected and built by the pipeline, not
  hand-made.
- **A single record**, so no record-join geometry is exercised. That is the
  other half of the coverage question and
  [`record_geometry_on_drosophila_2026-09-21.md`](record_geometry_on_drosophila_2026-09-21.md)
  measures it separately.
- **The directional figure is not what a fixed tool would deliver.** These
  panels were SELECTED under the symmetric model. A set cover that knew each
  site covers `r` rather than `2r` would choose different primers and do better
  than the directional column here. So the table bounds how much the reported
  figure overstates THIS panel, not how much coverage the tool would lose.
- No wet-lab outcome is involved. Neither column is validated against measured
  sequencing breadth.

## What follows, and what it costs

Phase 5.2 is now a decision with a number attached rather than an open
question. Acting on it is not a small change:

- `fg_coverage` is the authoritative figure in every
  `step4_improved_df_summary.json` and quoted throughout `docs/validation/`.
  Every recorded figure moves.
- `reach_calibration.md` fitted the 3 kb reach **against the symmetric model**,
  so the reach itself would need refitting: a directional model with a refitted
  reach is not the same as this table's directional column.
- Selection would change, not just reporting, if Stage 1's set cover is taught
  the directional geometry -- which is the point of fixing it, and why the
  delivered panels would differ.

The cheapest honest intermediate step is to report both figures and label which
is which, as the hybrid run already does for its two coverage numbers, rather
than silently swapping the definition of the headline metric.
