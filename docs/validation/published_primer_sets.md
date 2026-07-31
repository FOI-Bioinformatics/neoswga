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

## Finding 2: `normalized_score` does not reproduce the outcome

`PrimerSetMetrics.normalized_score` weights coverage, selectivity, dimer risk,
evenness and Tm. It **computes `max_gap` but never uses it**.

Scored on the published statistics:

| Rank by our score | Set | Score | Actual |
|---|---|---|---|
| 1 | Prev02 | 0.684 | 2.1x |
| 2 | Prev04 | 0.663 | 13x |
| 3 | Prev03 | 0.659 | **96x** |
| 4 | Prev06 | 0.656 | **120x** |
| 5 | Prev05 | 0.654 | 2.0x |
| 6 | Prev01 | 0.610 | 2.0x |

The two wet-lab winners rank 3rd and 4th, and a 2.1x failure ranks first. The
same ordering holds for the `balanced`, `clinical` and `enrichment` application
profiles.

The score is also barely discriminating: the full spread is 0.610-0.684, about
7 %, across sets whose measured performance spans 60-fold.

Tracked as a strict `xfail` in
`tests/validation/test_published_primer_sets.py::test_normalized_score_ranks_the_wet_lab_winners_first`
so the gap cannot be silently tolerated. It is **not** fixed here: changing the
scoring re-ranks all optimizer output and needs its own change with integration
baselines regenerated deliberately.

## Caveat

**n = 6.** These correlations are suggestive, not established, and any
re-weighting tuned on six points risks fitting noise — with this many sets it is
trivial to find weights that rank perfectly and mean nothing.

What the data does support firmly is the negative result: the current scoring
ranks the winners 3rd and 4th on the only real outcome data available, and the
metric previously treated as dominant has no correlation here. That is enough to
justify re-examining the weights; it is not enough to fix them by fitting.

## Recommended next step

Include `max_gap` in `normalized_score` as a coverage-hole penalty. It is
already computed, it is the strongest single predictor here, and it is
mechanistically motivated independently of this dataset.

Before committing to weights, extend the benchmark: the Clarke et al. (2017)
*M. tuberculosis* sets are the obvious second source, and any in-house sets with
measured outcomes would help more than tuning on these six.

## Follow-up: not yet done

- Clarke et al. (2017) *M. tuberculosis* sets are not yet in the fixture.
- Our own metrics are validated against the *published* statistics, not
  recomputed from the genome. Recomputing would need GCF_000144405.1 plus a
  human background, which is too large to fetch during a test run; a separate
  opt-in benchmark would close that loop.
