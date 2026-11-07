# `off_gap_gini`, computed for the first time, and what it sees

Measured 29 October 2026, on `main` at `162a699` plus the change this document
records. Phase 7.2 of
[the plan](../reports/PLAN_2026-10-28_pool_design_trustworthiness.md).
Reproduce with:

```
PYTHONPATH=. python docs/validation/off_gap_gini_2026-10-29/measure.py \
    --out results.json
```

The script generates its own genomes from a fixed seed, so it needs nothing on
disk. Results in [`off_gap_gini_2026-10-29/results.json`](off_gap_gini_2026-10-29/results.json);
regressions in `tests/test_the_background_evenness_is_computed.py` (26 tests).

## What it is, and why it was not here

swga 2.0's set-level model has five terms, fitted by ridge regression with
10-fold cross-validation against the proportion of the target genome reaching
1x sequencing coverage, over 46 published sets from *M. tuberculosis* and
*H. sapiens*. `off_gap_gini` carries the second-largest recorded weight, +0.281
behind `freq_ratio` at +0.321: a Gini index of the distances between background
binding sites, computed per strand, averaged across strands and then across
genomes. A POSITIVE weight means an UNEVEN host is better, because clustered
host sites leave most of the host out of a polymerase's reach.

[`tool_comparison.md`](tool_comparison.md) recorded three times that nothing
here computed it. It is computed now, inside
`PositionCache.compute_strand_alternation_stats` as a sixth figure per genome
(`gap_gini_per_strand`) and aggregated over the background genomes by
`core/background_evenness.off_gap_gini`.

**It is reported and never scored.** The sample behind the weight is 46 sets
with one response variable, read from a table nobody here could open, and
COATswga computes no gap statistic at all on the stated ground that a
per-primer Gini cannot speak for a set. `tool_comparison.md` already took that
position and Known Issue 18 keeps `max_gap` out of every score for the same
reason. `test_it_does_not_enter_normalized_score` pins it.

## It sees a geometry the site counts cannot

The direction check, through the shipped function, on two synthetic hosts with
the **same site count** -- 400 per strand over 4 Mb:

| host geometry | sites per strand | per-strand gap Gini |
|---|---|---|
| dispersed (uniform at random) | 400 | 0.497 |
| clustered (ten clusters inside 2% of the host) | 400 | **0.966** |

`total_bg_sites` is 800 in both. `selectivity_density` and
`selectivity_ratio` are identical in both, since both are additive in
per-primer counts and no position enters them. So this is the one background
quantity besides `bg_coverage` that distinguishes the two, and it is the shape
the fitted weight claims to reward.

## On a real pipeline run

One synthetic pool, equiphi29 at 43 C, k 10-12, 1.27 Mb target at 35% GC,
4 Mb host at 42% GC, 500 candidates. One `optimize` per panel size.

| primers | off_gap_gini | target per-strand | bg_coverage | bg sites | sel_dens | fg_cov |
|---|---|---|---|---|---|---|
| 6 | 0.388 | 0.514 | 0.0240 | 12 | 4.152 | 0.2615 |
| 12 | 0.518 | 0.427 | 0.0393 | 21 | 4.294 | 0.4749 |
| 18 | 0.527 | 0.384 | 0.0413 | 22 | 4.311 | 0.5964 |
| 24 | 0.491 | 0.362 | 0.0433 | 23 | 4.302 | 0.6912 |

**It moves, and not with the size.** The range is 0.14 over the four panels and
the trajectory is non-monotone, peaking at 18 primers. The target's own
per-strand figure falls monotonically as the panel grows, which is what more
sites should do to a spacing Gini; the host's does not, because adding a primer
adds host sites wherever that primer happens to bind.

**No correlation is claimed.** The rank correlations the script computes are
over FOUR panels, which is not a sample: `selectivity_density` reads 0.8 and
the other three read 0.4, and on four points those numbers are what random
orderings produce. The direction check above, not this table, is the evidence
that the figure carries information the counts do not. Whether it predicts
anything would need the response variable swga 2.0 fitted against, which is
sequencing breadth, and this repository has no BAM or CRAM anywhere.

## Cost

Computed from the per-strand position lists
`compute_strand_alternation_stats` already builds, so it costs arithmetic and
not a second pass through `get_positions`, which is the hottest method in the
package. On 20 primers against a warm cache: 0.07 ms with the figure against
0.05 ms without, so about 0.02 ms per call, against roughly 40 ms for one
`compute_metrics` call on the prepared Wolbachia example.

## What the run turned up on the way: the panels are one-stranded

On every panel this pool delivered, `strand_coverage_ratio` is 0.0 and
`strand_alternation_score` is 0.0 on BOTH genomes. Counted directly on the
24-primer panel:

| genome | forward sites | reverse-complement sites | index entries absent |
|---|---|---|---|
| target | 140 | 1 | 0 |
| host | 23 | 0 | 0 |

**These are measured zeros, not missing measurements.** Every reverse-complement
entry is present in the position index and empty, so the figures mean what they
say, and the 7.1 strand section flags the target correctly rather than
reporting an artifact. The distinction was checked explicitly because a 0.0
there is exactly what Known Issues 5, 6 and 13 are about.

**The mechanism is selection, in a low-multiplicity regime.** The genome is
generated i.i.d. with p(A)=p(T) and p(G)=p(C), under which a k-mer and its
reverse complement have EQUAL expected counts. At k=11-12 on 1.27 Mb that
expectation is a fraction of one occurrence, so a candidate that passes
`min_fg_freq` is a large upper-tail draw on its own count; its reverse
complement is an independent draw with the same small mean, and lands at zero.
Filtering on a k-mer's own frequency therefore selects for one strand whenever
the expected multiplicity is near zero.

How much of that survives on a real genome is **not determined here**: real
sequence obeys Chargaff's second parity rule and carries far higher k-mer
multiplicity at these lengths, so the asymmetry should be milder, and none of
the five reference genomes `CLAUDE.md` lists is reachable from this container.
The observation is recorded because it is what a one-strand panel looks like
and because it says which pools the strand figures will read as degenerate --
not as a claim about any real design.

## Scope

One synthetic pool, four panel sizes, one reaction, one refinement method; one
direction check at one site count and two geometries. Nothing is selected,
scored or constrained differently: `normalized_score`, every panel limit and
every delivered panel are unchanged. No wet-lab outcome is implied, and the
coefficient the term carries in swga 2.0 is recorded rather than verified --
[`pool_selection_audit_2026-09-18.md`](pool_selection_audit_2026-09-18.md)
explains why.
