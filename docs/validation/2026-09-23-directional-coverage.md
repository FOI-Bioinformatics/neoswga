# Symmetric against directional coverage windows

Measured 23 September 2026. This reconciles the figure already recorded in
[reach_calibration.md](reach_calibration.md) under **What is not established**,
and reports what that figure does not: which bases the two models claim, and
that the directional model is already implemented and already used by every
method's Stage-1 selection, while everything that SCORES a panel is symmetric.

## The defect

`coverage._mark_window` marks `occupied[pos - r : pos + r]`
(`coverage.py:242`), and the sites come from `get_positions(..., strand="both")`
(`coverage.py:28`, `:117`). A site can only extend one way:

- An oligo occurring literally at `i` anneals to the MINUS strand there, the
  nascent strand is plus-sense, and extension runs toward increasing
  coordinates.
- An occurrence of its reverse complement at `j` is the mirror: the oligo
  anneals to the PLUS strand and extension runs toward decreasing coordinates.

So each occurrence is credited in the direction it cannot extend.

`PositionCache` already separates the two -- `"forward"` is `db[primer]` and
`"reverse"` is `db[rc(primer)]`, both plus-strand offsets
(`position_cache.py:374-386`) -- so nothing needs re-indexing.

## Two directional variants, and the distinction is the whole reconciliation

| name | window per site | isolates |
|---|---|---|
| `one-sided r` | `r`, on the reachable side | geometry AND magnitude together |
| `one-sided 2r` | `2r`, on the reachable side | geometry alone |

`one-sided 2r` credits the same total per site as the symmetric model and
differs only in where those bases fall.

## Delivered panels

```bash
python scripts/benchmarking/directional_coverage.py panels
```

Prevotella, 3,168,282 bp, reach 3,000. The symmetric column reproduces
production exactly: the out12 panel's recorded `coverage` is
`0.06234830106663485` and this returns `0.0623`.

| panel | k | sites | symmetric | 1-sided r | 1-sided 2r | map agreement |
|---|---|---|---|---|---|---|
| out10 | 10 | 386 | 0.4688 | -37.1% | +8.0% | 72.4% |
| out11 | 11 | 51 | 0.0930 | -48.5% | -2.0% | 53.5% |
| out12 | 12 | 33 | 0.0623 | -49.9% | -1.5% | 51.5% |
| out12big | 12 | 61 | 0.1087 | -46.9% | +3.6% | 55.2% |

**The number is nearly right and the map is about half wrong.** At equal width
the headline moves by -2.0% to +8.0%, while only 51.5% to 72.4% of the bases
the symmetric model claims are also claimed by the directional one. At k=12,
95,784 bp are claimed by symmetric alone and 92,836 bp by directional alone.
The totals nearly cancel; the maps do not.

Gap statistics move modestly: out12 `max_gap` 320,415 -> 326,403, `mean_gap`
90,023 -> 92,928; out10 `max_gap` unchanged at 84,334 with `n_gaps` 189 -> 170.

## Why the totals nearly cancel

Per-primer strand balance is close to automatic. A k-mer and its reverse
complement each occur about `N / 4**k` times, so the counts differ by a
sampling fluctuation. Every primer in the out10 panel has skew
`|fwd - rev| / total` at or below 0.17, with per-primer deltas between -2.2%
and +4.0%. On the sparse k=12 panels the windows do not overlap at all, so the
two models cover an identical NUMBER of bases and the per-primer delta is
exactly 0.0% for five of six primers.

The symmetric model is approximately right because it is approximately
symmetric in aggregate, not because the mechanism is symmetric. Nothing in
selection enforces balance.

## Reconciling the recorded 4-14%

`reach_calibration.md:179-181` records, under **What is not established**:

> **The symmetric window is still an approximation.** Extension from a binding
> site is one-sided and strand-determined; measured on Prevotella, symmetric
> overstates coverage by 4-14% (growing with set size, shrinking with reach).
> Correcting it would move the fitted reach up slightly.

No script accompanies it and no panel is named, so it could not be reproduced
directly. Rebuilding its stated setup -- greedy maximum-coverage panels over
Prevotella 12-mers at the three calibrated reaches -- does reproduce it:

```bash
python scripts/benchmarking/directional_coverage.py trend
```

| reach | n=8 | n=16 | n=32 | n=64 | | n=8 | n=16 | n=32 | n=64 |
|---|---|---|---|---|---|---|---|---|---|
| | **vs one-sided 2r** | | | | | **vs one-sided r** | | | |
| 3.0 kb | +0.1% | -0.0% | +5.6% | +6.6% | | +43.0% | +41.5% | +40.3% | +36.1% |
| 4.5 kb | +2.8% | +4.4% | +6.1% | +5.1% | | +40.0% | +37.8% | +34.5% | +28.0% |
| 6.2 kb | +4.6% | +5.8% | +9.8% | +5.3% | | +38.5% | +36.7% | +33.0% | +20.6% |

**The recorded figure describes the width-preserving comparison.** The
`one-sided 2r` column is 0-10%, against the recorded 4-14%, and grows with set
size as recorded. The `one-sided r` column is 20-43% and SHRINKS with set size,
which matches neither the magnitude nor the trend.

That settles the question the two figures appeared to disagree about: they are
the same quantity, and the recorded one is the width-preserving correction.
It also makes `reach_calibration.md`'s own consequence -- "would move the
fitted reach up slightly" -- consistent, which a 37-50% correction would not
be.

**The "shrinking with reach" half does not reproduce.** At n=32 the
width-preserving gap is 5.6%, 6.1% and 9.8% as reach grows from 3.0 to 6.2 kb.
At n=64 it is roughly flat. Not investigated further.

## The directional model already exists, and selection already uses it

`dominating_set_optimizer.BipartiteGraph.add_primer_coverage:205-224`
implements it, with a docstring stating the mechanism:

> forward sites extend downstream only, reverse sites extend upstream only.
> This models the biological reality that phi29 polymerase extends from the 3'
> end of the primer in one direction along the template.

It takes that branch only when the caller supplies the separated strands, and
`optimize_greedy` (`:885`) does, at `:939` and `:957`. **Every method's Stage-1
set cover therefore selects directionally** -- `hybrid` through
`hybrid_optimizer.py:708`, `background-aware` through the same, and
`dominating-set` through `dominating_set_adapter.py:169`.

**Everything that SCORES is symmetric.** So a run selects on one geometry and
is judged on another:

| stage | function | geometry |
|---|---|---|
| Stage-1 set cover | `optimize_greedy` -> `add_primer_coverage:205-224` | **directional**, binned |
| Stage-2 swap | `occupancy_coverage` -> `merged_window_intervals` | symmetric, base-resolution |
| `fg_coverage`, `bg_coverage` | `base_optimizer._union_coverage` -> `_mark_window` | symmetric, base-resolution |
| hybrid's "(estimated, binned)" figure and the `_prune_background` floor | `hybrid_optimizer._calculate_coverage:1235` | symmetric, binned |
| hybrid's swap bins | `hybrid_optimizer._coverage_bins_by_primer:1140` | symmetric, binned |

`hybrid_optimizer.py:1138` has the separated arrays in hand at the last two and
concatenates them:

```python
positions = np.concatenate([fw, rv])
```

**This is the finding, and it is not the one this document originally
recorded.** A first pass reported that the default method does not use the
directional model at all, having found the two symmetric call sites in
`hybrid_optimizer` and not the `optimize_greedy` call above them. That was
wrong. The defect is not that one method uses the worse model; it is that
every method uses BOTH, at different stages, and nothing reconciles them.

CLAUDE.md's **Coverage reach (important)** records that the realistic reach was
threaded into Stage-1 set cover so that "selection and the reported
`fg_coverage` agree (and ensemble comparisons are fair)". The REACH was
reconciled. The GEOMETRY was not, and the stated goal is not met while it
differs.

## What this does not establish

- **No wet-lab outcome anchors either geometry.** `reach_phi29` is `assumed` in
  `registry/model_evidence.json:278-291`, `calibrate-reach --bam` exists, and
  there is no BAM in this repository.
- **Whether the reach value should change.** The value and the geometry were
  fitted together in `reach_calibration.md`; correcting one without refitting
  the other is unsound in either direction. Not done here.
- **Whether making the scoring directional moves the delivered panel.**
  Everything above scores already-delivered panels. Since Stage-1 already
  selects directionally, aligning the scoring would change what the swap
  refinement optimises and what the result is judged on, and neither was
  measured here. That is the question deciding whether this is a reporting
  change or a design change.
- **Which geometry is right for the BIN-level Stage-1 model.** Its directional
  branch is one-sided at `r`, not at `2r`, so Stage 1 and the width-preserving
  correction measured above are not the same model either.
- **The record-confinement disagreement is untouched and is present in both
  columns.** `compute_per_prefix_coverage` confines windows to record
  boundaries and `merged_window_intervals` does not. Prevotella has two
  records, so a single join sits inside these figures either way.
- **One target, one background, one reach convention.** Prevotella against
  nothing, at three reaches, on four delivered panels and one greedy sweep.

## Eight other coverage implementations exist

Found while tracing, recorded so the next reader does not assume `coverage.py`
is the only one. Most are dead or test-only and none is on the path that
produces `fg_coverage`:

| site | geometry |
|---|---|
| `swga_simulator._calculate_bin_coverage:469` | bin-only; a site claims its whole 10 kb bin |
| `minimal_primer_selector:115` | no window at all -- coverage is `len(sites) / genome_length` |
| `advanced_features._estimate_coverage:457` | symmetric, literal 3000 rather than the polymerase's |
| `amplicon_network.calculate_coverage:203` | connected-component span; an isolated site covers 0 bp |
| `simulation_analysis._analyze_coverage:150` | bin-only |
| `primer_expansion.identify_gaps:268` | no window; raw inter-site distances |

`primer_expansion._recovered_deficit_fraction:573` is the only one in the set
that passes `record_starts`.
