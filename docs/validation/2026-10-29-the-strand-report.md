# The strand report: ten figures measured, two reported, one of them wrong

Changed 29 October 2026, on `main` at `c4f15f0` plus the change this document
records. Phase 7.1 of
[the plan](../reports/PLAN_2026-10-28_pool_design_trustworthiness.md).
Reproduce with:

```
PYTHONPATH=. python docs/validation/strand_report_2026-10-29/probe.py
```

The regression is `tests/test_the_strand_report_reads_every_genome.py`
(26 tests).

## What was measured and what was read

`base_optimizer._compute_metrics` calls `strand_metrics.collect_strand_stats`,
which returns five figures for every foreground genome and five for the host,
and `PrimerSetMetrics.to_dict` serialises all of them under
`metrics.strand_stats` in `step4_improved_df_summary.json`. Ten figures on a
one-host run.

The quality report rendered two: the alternation score and the site-count
balance of the FIRST target genome. The host's strand structure, both
convergent-gap figures, and the longest same-strand run were computed,
written to disk, and read by nobody — Known Issue 8's shape inside the report
rather than inside the search.

## The one that was wrong, not merely absent

Alternation is the fraction of ADJACENT site pairs on opposite strands, so a
panel with fewer than two sites on a genome is 0/0.
`PositionCache.compute_strand_alternation_stats` writes `null` for that case on
purpose, and `tests/test_strand_scores_say_when_they_are_unmeasurable.py` pins
it. The report read it through `_safe_float`, whose default is the literal
`0.0`, and rendered:

```
Strand alternation score   0.00   fraction of adjacent sites that alternate
                                  strands (higher = better exponential
                                  amplification)
```

"Could not be measured" shown as "measured: no alternation at all", in the
direction that looks like a finding. Known Issues 5, 6 and 13. The same
module already reads the gap figures with `default=None` for exactly this
reason; the strand pair was missed.

## What the section renders now

From the probe, on a panel with 25 forward sites to 1 reverse on the target:

```
Strand alternation score   not measurable   ... 0/0 and so not measurable
                                            below two sites
Strand coverage ratio      0.04             min(forward, reverse) /
                                            max(forward, reverse) over binding
                                            site COUNTS, not per-strand
                                            coverage (1.0 = balanced)

Genome  Role    Site balance  Alternation      Widest gap  Mean gap  Longest run
host    host    0.95          0.71             48,000      9,100     3
target  target  0.04          not measurable   412,000     61,000    25

- host (host): widest stretch with no convergent pair is 48,000 bp, 16.0x the
  3,000 bp reach. Larger is better here. No reference.
- target (target): widest stretch with no convergent pair is 412,000 bp,
  137.3x the 3,000 bp reach, so no two sites face each other across it.
  No reference.
```

With the target's reverse strand empty instead:

```
- target (target) has binding sites on one strand only, so it holds no
  convergent pair anywhere and amplifies linearly at best. This is the
  imbalance a mean of the two strand coverages hides.
```

## Two readings the section has to keep apart

**A starved strand on the target is bad and on the host is good.** Exponential
amplification needs two sites in convergent orientation within the
polymerase's reach; one site primes linearly at best. The role comes from
`fg_prefixes`/`bg_prefixes` in params.json, and a prefix params.json does not
place is labelled `unknown` rather than defaulted to either.

**An absent genome is not a zero one.** A prefix missing from `strand_stats`
gets no row. `collect_strand_stats` omits a prefix the cache could not answer
for, and a row of zeros would read as a measurement.

## What is flagged, and what is deliberately not

Known Issue 18 keeps the strand scalars unconstrainable, and the reason given
there is that an unmeasured zero cannot be told from a measured one. That
reason does not extend to reporting, which is what this closes. But it does
mean **no threshold is introduced**, because none has a reference. Two things
are flagged:

- **A fact.** `strand_coverage_ratio` is min/max over the site counts, so 0.0
  with any site at all means the minority strand carries none.
- **A mechanism.** The widest stretch with no convergent pair, in units of the
  polymerase's reach. The reach is the distance the enzyme can span, so a gap
  wider than it is a statement about geometry, not a calibrated limit. Every
  such note ends "No reference".

At most one note per genome, and the starved-strand note wins: the gap figures
hold the GENOME LENGTH when no opposite-strand pair exists at all (a
conservative encoding, which `worst_convergent_gap` also reads that way), and
printing that as "423x the reach" would dress the encoding up as a
measurement.

## The gap this closes, and the one it does not

[`2026-09-22-swga-comparative-review.md`](2026-09-22-swga-comparative-review.md)
records that COATswga ranks on the MEAN of its two strand coverages, so
0.99/0.30 scores 0.645 and beats a balanced 0.60/0.60, and nothing in it flags
the imbalance — then says NeoSWGA "has the data to make the equivalent check
and does not make it either". It does now, as a report.

**It is not the same quantity.** COATswga's two numbers are per-strand
COVERAGE fractions; `strand_coverage_ratio` is a balance of site COUNTS. The
section says so in the cell it renders, because a reader who takes one for the
other would think a 0.95 ratio meant 95% of the coverage on the weaker strand.
A per-strand coverage figure would need the directional coverage model
[`2026-10-29-directional-coverage-measured.md`](2026-10-29-directional-coverage-measured.md)
measures, which is a separate and larger decision.

## Scope

Nothing is selected, scored or constrained differently. `normalized_score`,
every panel limit and every delivered panel are unchanged; the whole change is
in what the report reads off a summary that already carried it, plus one
conversion default. No wet-lab outcome is implied by any figure here, and no
threshold on any of them is claimed.
