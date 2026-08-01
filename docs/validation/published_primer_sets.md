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

---

# Second benchmark: per-primer plasmid amplification

**Data:** same paper, S1 Table — 310 qPCR measurements of fold amplification for
284 primers on the pcDNA and pLTR plasmids. Those are the two plasmids shipped in
`examples/plasmid_example`, so binding sites are recomputed with this codebase
rather than taken on trust.
**Test suite:** `tests/validation/test_plasmid_amplification.py`
**Fixture:** `tests/validation/data/dwivedi_yu_2023_plasmid_amplification.json`

## What this benchmark can and cannot test

Every primer has **exactly one on-target site and zero off-target sites** — they
were selected that way. Site geometry is constant, so this dataset cannot speak
to whether foreground/background site ratio predicts selectivity: the predictor
has no variance.

That is worth recording because the obvious analyses are actively misleading:

- Spearman of site ratio against measured selectivity returns about **-0.01**.
  That reads as "site ratio does not predict selectivity" and would be a false
  and consequential conclusion. It actually means "no variance".
- Splitting into deciles by site ratio appears to show a **four-fold** effect
  (median selectivity 0.53 vs 2.12). That is entirely an artefact of how ties get
  ordered.

Both were computed during this work before the constant geometry was noticed.

## What it does test: sequence properties predict amplification

With site geometry held fixed, measured on-target amplification still spans
0.0x to 161x. Something other than binding-site count drives it, and
sequence-level properties track it:

| Property | Spearman vs measured on-target amplification |
|---|---|
| **primer length** | **+0.52** |
| **melting temperature** | **+0.43** |
| GC fraction | +0.12 |

Length and Tm are correlated with each other, so these are not independent
signals. Still, this is a direct check on the premise behind the pipeline's
thermodynamic filtering, and it holds at n = 310 — a far better-powered result
than the six-set benchmark above.

None of these predict measured *selectivity* (all rho about 0.1), which is
expected: with identical site geometry there is little selectivity to predict.

## A real bug this benchmark caught

Recomputing the site counts surfaced a defect in the scan fallback added for
`evaluate-set`: `string_search.get_all_positions_per_k` scans the **forward
strand only**, so the reverse complement has to be queried as a separate k-mer —
which is what `string_search.get_positions` does when building the HDF5 index.
The scan fallback omitted it.

The effect was that every reverse-strand binding site was silently lost: **147 of
568** primer/genome counts on this dataset came back zero, each of which would
read as "this primer does not bind" rather than "we only looked at one strand".
146 of 284 primers appeared to bind neither plasmid.

Fixed, and `test_scanner_agrees_with_a_naive_search` now compares our scanner
against an independent search for all 284 primers on both plasmids. The single
remaining difference is a genuine circular wrap-around that the naive search does
not model.

---

# Third benchmark: Clarke et al. (2017) M. tuberculosis sets

**Data:** Clarke EL, Sundararaman SA, Seifert SN, Bushman FD, Hahn BH, Brisson D
(2017) *Bioinformatics* 33(14):2071-2077. Table 2 (set characteristics, with the
four most effective sets footnoted) and Supplementary Table 2 (primer sequences),
obtained via EuropePMC (PMC5870857) since the publisher PDF is not directly
retrievable.
**Test suite:** `tests/validation/test_clarke_2017_mtb.py`
**Fixture:** `tests/validation/data/clarke_2017_mtb.json`

Twelve sets against *M. tuberculosis* in a human background, including two
chosen deliberately as negative controls (`MtbUneven`, highest Gini;
`MtbSparse`, sparsest binding).

## This benchmark contradicts the Prevotella one

| Metric | Clarke Mtb (n=12) | Dwivedi-Yu Prevotella (n=6) |
|---|---|---|
| mean binding distance | **separates winners** | does not |
| worst gap (`max_gap`) | does not | **separates winners** |
| Gini (evenness) | does not (runs slightly backwards) | **separates winners** |

The four most effective Mtb sets are the four densest, cleanly separated with no
overlap. Gini runs the *wrong* way: the winners' median Gini is 0.522 against
0.482 for the sets they beat.

## Why, and what to take from it

The papers' methods account for the difference. Clarke pre-filtered the primer
pool for even binding before assembling sets, "constrain[ing] the range of set
binding evenness ... by removing primers that cluster on repeat regions". With
evenness already controlled it had little residual variance to explain the
outcome, so density dominated. Dwivedi-Yu did not control evenness, so it varied
across their sets and the coverage holes dominated.

**Neither "density predicts" nor "evenness predicts" is universally true. What
predicts is whichever property is currently limiting.** That is a more useful
rule than either dataset alone yields, and it is only visible because there are
now two benchmarks that disagree.

This revises the conclusion drawn from the Prevotella data alone, which had been
carried into the swga-chemistry skill as a general claim. Both the skill and
`references/troubleshooting.md` are corrected.

## A trap worth recording

`fg_bg_ratio` is defined identically in both papers (fg_mean_distance /
bg_mean_distance), and it appears to separate winners in *both* datasets — but in
**opposite directions**: lower is better for Clarke, higher for Dwivedi-Yu.

An analysis that supplies the expected direction per dataset will report it as a
consistent predictor. It is not one, and that circularity was present in a first
pass of this comparison before being caught.

## Status of the `max_gap` question

The original motivation for a third benchmark was to decide whether to add a
`max_gap` coverage-hole term to `normalized_score`. The answer is now clearer,
and it is *not* a straightforward yes: `max_gap` separates winners on Prevotella
and does not on Clarke. A term weighted for one dataset would be wrong for the
other.

The remaining honest position is unchanged from before: do not tune scoring
weights on this much data. What both benchmarks support is reporting `mean_gap`,
`max_gap` and `gap_gini` alongside any composite score so the user can see which
regime they are in, rather than compressing them into one number.

---

# Fourth benchmark: Clarke et al. (2017) Wolbachia sets

**Data:** same paper, Table 1 / Supplementary Table 1, with outcomes quoted from
Results 3.1. Five sets against *W. pipientis* in a *D. melanogaster* background,
including the original hand-designed set from Leichty and Brisson (2014).
**Test suite:** `tests/validation/test_clarke_2017_wolbachia.py`

This one is unusually clean because the authors built it as a controlled
contrast: within each of two primer melting-temperature ranges they took the
**most selective** set and the **most even** set, pitting the two design criteria
directly against each other.

| Set | selected for | Tm range | Gini | mean gap | max gap | outcome |
|---|---|---|---|---|---|---|
| `TmL/Even` | evenness | low | **0.537** | 6853 | 34.7 kb | **effective** — 10x over 60–75 % of genome |
| `TmL/Selective` | selectivity | low | 0.654 | 5327 | 31.1 kb | partial — uneven coverage |
| `Leichty 2014` | hand-designed | — | 0.712 | **5305** | 112.8 kb | partial — uneven coverage |
| `TmH/Even` | evenness | high | 0.537 | 13070 | 112.8 kb | ineffective |
| `TmH/Selective` | selectivity | high | 0.660 | 12074 | 128.2 kb | ineffective |

Unamplified control reached 10x over 2.8 % of the genome.

**Evenness wins, and density runs backwards.** `TmL/Even` was effective with
*worse* mean spacing than both sets it beat. The hand-designed Leichty set has
the best density of all five and the worst evenness, and did not improve
sequencing efficiency.

A second, independent failure mode also shows up: both high-Tm sets failed
regardless of which criterion they were selected for. Raising the Tm window
shrank the usable primer pool, roughly doubling mean spacing and pushing max gaps
past 110 kb. Tm range set the outcome before either criterion could apply.

---

# Synthesis across the three set-level benchmarks

| Dataset | n | density predicts? | evenness predicts? |
|---|---|---|---|
| Prevotella (Dwivedi-Yu 2023) | 6 | no | **yes** (with worst gap) |
| Wolbachia (Clarke 2017) | 5 | no (backwards) | **yes** |
| *M. tuberculosis* (Clarke 2017) | 12 | **yes** | no (backwards) |

Two of three favour evenness. The odd one out is *M. tuberculosis*, and the paper
states why: there the primer **pool** was pre-filtered for even binding, which
removed most of the variance in evenness and left density as the limiting
property.

All three are consistent with a single rule: **whichever property is currently
limiting is the one that predicts.** That is what the reporting change is built
around — `mean_gap`, `max_gap` and `gap_gini` are now shown together with their
benchmark context in both HTML reports and in `evaluate-set`, rather than folded
into a composite that would have to commit to one regime.

## Why the scoring weights were still not changed

With 6, 5 and 12 sets, any weighting fitted to one dataset is contradicted by
another. A `max_gap` term helps on Prevotella, is neutral-to-harmful on
*M. tuberculosis*, and does not cleanly separate the Wolbachia sets either
(`TmL/Selective` has a *better* max gap than the winner). The composite score is
left alone; the three statistics are surfaced so the reader can see the regime.

## Datasets sought but not obtained

- **Oyola et al. (2016)** *Plasmodium falciparum* sWGA (Malaria Journal
  12936-2016-1641). Three probe pools (10, 20, 28 primers) with measured yield
  and >18-fold enrichment, but the primer sequences sit in Additional file 1
  Table S2, which is not retrievable through EuropePMC or the Springer static
  endpoints. Its outcome axis is also different — yield by pool size rather than
  per-set coverage — so it would test set *size*, not gap structure.
