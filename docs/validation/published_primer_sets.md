# Validation against published primer sets with wet-lab outcomes

**Date:** 2026-07-31
**Data:** Dwivedi-Yu et al. (2023) *PLOS Comput Biol* 19(4):e1010137, S2 Table
(six *Prevotella melaninogenica* primer sets with their in-silico statistics) and
S3 Table (percent of reads mapping to target after amplification).
**Test suite:** `tests/validation/test_published_primer_sets.py`
**Fixture:** `tests/validation/data/dwivedi_yu_2023_prevotella.json`

This is a rare dataset: designed primer sets, their computed metrics, and a
measured outcome for each. It is the only wet-lab-validated benchmark we have.

## The data

| Set | n | mean spacing | max gap | Gini | fg/bg | on-target reads | vs control |
|---|---|---|---|---|---|---|---|
| Prev06 | 3 | 4980 bp | **31.3 kb** | **0.552** | **0.0980** | 72.22 % | **120x** |
| Prev03 | 4 | 4100 bp | **31.5 kb** | **0.533** | **0.0916** | 57.34 % | **96x** |
| Prev04 | 9 | **1980 bp** | 33.3 kb | 0.623 | 0.0664 | 7.70 % | 13x |
| Prev02 | 3 | 5030 bp | 42.2 kb | 0.566 | 0.0494 | 1.24 % | 2.1x |
| Prev01 | 2 | 7390 bp | 65.3 kb | 0.600 | 0.0538 | 1.23 % | 2.0x |
| Prev05 | 6 | 2760 bp | 37.7 kb | 0.638 | 0.0804 | 1.20 % | 2.0x |

Unamplified control: 0.60 % on-target.

## Finding 1: mean binding distance does not predict success

Spearman rank correlation against fold enrichment (n = 6):

| Metric | rho | direction |
|---|---|---|
| **max gap** (largest coverage hole) | **-0.83** | lower is better |
| **Gini** (evenness) | **-0.77** | lower is better |
| fg/bg ratio | +0.60 | higher is better |
| **mean binding distance** | **-0.09** | effectively none |
| set size | +0.09 | none |

Ranking the six sets by each metric and taking the top two:

| Ranked by | Top 2 | Matches the wet-lab winners? |
|---|---|---|
| max gap | Prev06, Prev03 | yes |
| Gini | Prev03, Prev06 | yes |
| fg/bg ratio | Prev06, Prev03 | yes |
| mean binding distance | Prev04, Prev05 | **no** |

Ranking by max gap alone reproduces the measured top three in exact order
(Prev06, Prev03, Prev04).

This contradicts the widely repeated summary that "the set with the lowest mean
binding distance wins", which had been carried into this project's guidance.
Prev04 has the best mean spacing of all six and finished third; Prev06 has the
second-worst and won outright; Prev05 is denser than both winners and failed.

The mechanism is straightforward: a region no primer can reach stays unamplified
however tightly spaced the sites are elsewhere. Average density hides holes;
`max_gap` and `gap_gini` measure them.

## Finding 2 (CORRECTED): the shipped scoring does rank the winners first

An earlier version of this report claimed `normalized_score` ranked the two
winners 3rd and 4th and put a 2.1x failure first. **That was wrong.** The cause
was a bug in the reconstruction used to build `PrimerSetMetrics` from the
published statistics, not in the scoring: background sites were counted over the
3.2 Mb *target* genome instead of the 3.1 Gb *human background*, understating
them roughly 1000-fold and inverting the selectivity ordering.

With the background sized correctly:

| Rank by our score | Set | Score | Actual |
|---|---|---|---|
| 1 | Prev03 | 0.542 | **96x** |
| 2 | Prev06 | 0.540 | **120x** |
| 3 | Prev02 | 0.533 | 2.1x |
| 4 | Prev04 | 0.532 | 13x |
| 5 | Prev05 | 0.528 | 2.0x |
| 6 | Prev01 | 0.479 | 2.0x |

The top two are the wet-lab winners, with no change to the scoring weights.

Two caveats remain, and they are real:

- **The margin is thin.** The winners lead the next set by about 0.008 of score,
  against a 60-fold spread in measured enrichment. The ordering is right; it is
  not a confident separation.
- **`max_gap` is still unused.** `normalized_score` computes it and never reads
  it, even though it is the strongest single predictor in this dataset. A
  coverage-hole term was implemented and tested against this benchmark: it
  changed nothing in the ranking, so it was **not shipped**. A change that
  re-ranks all optimizer output has to demonstrate a benefit first, and on the
  only validation data available this one does not.

## Caveat

**n = 6.** These correlations are suggestive, not established, and any
re-weighting tuned on six points risks fitting noise — with this many sets it is
trivial to find weights that rank perfectly and mean nothing.

What the data supports firmly is Finding 1: the metric previously treated as
dominant (mean binding distance) has no correlation with outcome here, while the
worst coverage hole and evenness both track it. Those are computed directly from
the published tables and do not depend on any reconstruction.

Finding 2 is weaker in both directions: the scoring gets the order right, but by
a margin too small to call a real separation on six sets.

## Recommended next step

Extend the benchmark before changing any weights. The Clarke et al. (2017)
*M. tuberculosis* sets are the obvious second source, and any in-house sets with
measured outcomes would be worth more than either.

A `max_gap` coverage-hole term remains the most plausible improvement --
mechanistically motivated, already computed, and the strongest single predictor
here. It was implemented and measured against this benchmark and made no
difference to the ranking, so it is not in the codebase. Revisit it when there
is enough outcome data to tell whether it helps.

## Follow-up: not yet done

- Clarke et al. (2017) *M. tuberculosis* sets are not yet in the fixture.
- Our own metrics are validated against the *published* statistics, not
  recomputed from the genome. Recomputing would need GCF_000144405.1 plus a
  human background, which is too large to fetch during a test run; a separate
  opt-in benchmark would close that loop.
