# Validation against published primer sets with wet-lab outcomes

**Date:** 2026-07-31
**Data:** Dwivedi-Yu et al. (2023) *PLOS Comput Biol* 19(4):e1010137, S2 Table
(six *Prevotella melaninogenica* primer sets with their in-silico statistics) and
S3 Table (percent of reads mapping to target after amplification).
**Test suite:** `tests/validation/test_published_primer_sets.py`
**Fixture:** `tests/validation/data/dwivedi_yu_2023_prevotella.json`

See [tool_comparison.md](tool_comparison.md) for how NeoSWGA compares with
swga 1.0, swga 2.0 and COATswga, including why the coverage reach and the
set-level scoring differ deliberately.

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

## The default Tm window costs most of the measured amplification

Running the 141 pcDNA-designed primers through the shipped phi29 filter
settings (`min_tm` 15, `max_tm` 45) splits them three ways:

| | n | median measured amplification |
|---|---|---|
| below the window | 14 | 0.80x |
| inside the window | 51 | 6.45x |
| **above the window** | **76** | **19.95x** |

The lower bound earns its place: the primers it rejects genuinely amplified
badly. The upper bound does not behave like a filter on quality. The 76 primers
above it amplified about three times better than those inside, and carry **74%
of all amplification measured in the experiment**. Only 2 of the 141 survive
the complete filter.

This is recorded rather than acted on, because the benchmark cannot settle the
question it raises. Every primer here has exactly one on-target site and zero
off-target sites, so there is no selectivity to lose — the same constant
geometry noted above. At genome scale a longer, higher-Tm primer binds more
sites including background ones, and an upper Tm bound is defensible as a
selectivity constraint rather than an amplification one.

What was previously unmeasured is the cost side of that trade. Raising `max_tm`
is not consequence-free in either direction, and the decision needs data this
benchmark cannot provide: per-primer amplification measured against a genome
with real background. Pinned by
`tests/validation/test_tm_window_vs_measured_amplification.py`.

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

# Fifth benchmark: Leichty and Brisson (2014), the original SWGA paper

Data: *Genetics* 198:473, Table S1 (sequences) and Results (outcomes). Three sets:
`B31-BL21` (20 primers, *Borrelia burgdorferi* against *E. coli*) and `SR1` / `SR2`
(10 and 2 primers, *Wolbachia pipientis* against *D. melanogaster*).

Fixture `tests/validation/data/leichty_brisson_2014.json`, tests
`tests/validation/test_leichty_brisson_2014.py`.

## A nested-set experiment

`SR2` is a strict subset of `SR1` — the same two primers, nothing else. It was the
*more* selective of the two: with restriction digest it reached 1264-fold target
amplification against 9-fold background (~140x), and up to 70 % of reads on target
against 2.4–8.8 % unamplified. Fewer, more selective primers beat more primers on
enrichment. This is the only place in the suite where set size varies with
composition held fixed.

## It documents the other half of a trade-off

Clarke et al. (2017) later benchmarked this same `SR2` set — it is their
"Leichty and Brisson 2014" row — and found it gave *uneven* coverage and did not
improve sequencing efficiency, despite competitive binding density. It has the
worst Gini of the five Wolbachia sets.

Both results are correct. They measure different things: `SR2` maximises
enrichment and sacrifices evenness. Evenness was never a design criterion here —
the 2014 selection scored target/background k-mer frequency ratio with a Tm cutoff,
predating the Gini index in the `swga` toolkit entirely.

This is the clearest evidence in the suite that **a single composite score cannot
serve both goals**, and it is only visible by combining two papers. Which end of
the trade-off you want depends on whether you need target DNA or even coverage.

---

# Sixth benchmark: Dwivedi-Yu et al. (2023) Leishmania braziliensis sets

Data: *PLOS Pathogens* 19(3):e1011230, S1 Table (both sheets) for the sets,
Results for the outcomes. Four sets of ten 8-mers drawn from a shared pool of 23,
designed with `swga2.0` against **three** backgrounds (human, *S. aureus*,
*S. pyogenes*). PS1 and PS4 reached ≥60 % parasite-mapping reads from 1 % and
0.1 % starting parasite DNA; PS2 and PS3 "performed more poorly on synthetic
samples" and were relegated to nested second-round reactions.

Fixture `tests/validation/data/dwivedi_yu_2023_leishmania.json`, tests
`tests/validation/test_leishmania_2023.py`.

## The closest thing here to a controlled A/B

PS1 and PS2 are **entirely disjoint** — twenty primers, no overlap — and PS3/PS4
are hybrids. Membership in the PS1 pool tracks the outcome exactly:

| Set | from PS1 | from PS2 | outcome |
|---|---|---|---|
| PS1 | 10 | 0 | effective |
| PS4 | 6 | 2 | effective |
| PS3 | 2 | 6 | poor |
| PS2 | 0 | 10 | poor |

## Composition does not explain it, and points the wrong way

Neither mean Tm, Tm spread, GC content nor 3' clamp separates the winners:

| Set | effective | mean Tm | Tm spread | mean GC |
|---|---|---|---|---|
| PS1 | yes | 28.98 | 0.83 | 0.613 |
| PS2 | no | 28.88 | 0.95 | 0.600 |
| PS3 | no | 28.49 | 1.65 | 0.588 |
| PS4 | **yes** | 27.85 | **1.71** | 0.575 |

PS4 won with the *widest* Tm spread and the *lowest* GC of the four. The composite
score carries a Tm-uniformity term; on this dataset that term points backwards.
One counter-example is not grounds for removing it, but it is grounds for not
treating uniform within-set Tm as established. Recorded as a test rather than
acted on.

What this dataset cannot test is binding geometry — that would need the 32 Mb
*L. braziliensis* assembly, too large to carry as a fixture. The narrow and honest
reading is: whatever separated these sets, it was not primer composition.

It also fills two gaps in the suite. It is the only **multi-background** design,
and the only one against a **GC-rich** target (~60 %), where the others are AT-rich
bacterial genomes — `B31-BL21`'s primers are near-zero GC. Both facts are pinned by
tests so the suite does not quietly narrow back to the easy case.

---

# Seventh benchmark: Oyola et al. (2016) Plasmodium falciparum

Data: *Malaria Journal* 15:597, Table S1 (Additional file 1, a PowerPoint whose
slide 5 carries the table) for sequences, Results for outcomes. One published set,
`Probe_10`, applied to 156 patients: >18-fold enrichment, 73.2 % (sd 4.4, N=5) of
reads on target against <1 % for standard WGA; parasitaemia threshold 0.03 %.

Fixture `tests/validation/data/oyola_2016_pfalciparum.json`, tests
`tests/validation/test_oyola_2016_pfalciparum.py`.

## It replicates the nested-set result

Three nested pools were compared on mock samples — `Probe_10` (first 10 primers),
`Probe_20` (first 20), `Probe_28` (all 28) — and the **smallest won**, on yield and
cost. Leichty and Brisson (2014) found the same with `SR2` beating the `SR1` it was
drawn from. Two groups, different organisms, twelve years apart, same direction:
adding primers to a working pool made it worse.

Neither is conclusive alone, and only `Probe_10`'s sequences were published so the
losing pools cannot be scored. But it is now two independent observations pointing
against "more primers means more coverage", which is worth not designing against.

## It was a set this tool could not design

*P. falciparum* is ~19 % GC and the working primers are **pure A/T** — one is
twelve consecutive adenines. Every GC rule rejected them, in two different ways.

| Rule | Before | After |
|---|---|---|
| GC clamp (`adaptive_filters`) | 0 / 10 | **10 / 10** |
| Adaptive GC window at genome GC 0.19 | 0 / 10 | **10 / 10** |
| GC content window (`parameter`, production path) | 0 / 10 | **10 / 10** |
| `filter_extra` overall | 0 / 10 | 3 / 10 |

Neither failure was a matter of threshold tuning:

- The **GC clamp** cannot be satisfied by a zero-GC primer at any setting. It is
  standard PCR primer-design practice, imported into a regime it does not fit —
  SWGA against an AT-rich genome at 30 °C, where AT-rich primers are the point.
- The **adaptive GC window** floored its lower bound at a positive value (0.20 in
  `adaptive_filters`, `max(0.15, genome_gc - tol)` in `parameter`). For a 0.19-GC
  genome the "adapted" window was 0.20–0.34 — it excluded the target's own average
  composition. A filter that cannot admit sequences at the target's GC is not
  adapting.

An important correction to an earlier version of this section: the clamp was **not**
the binding constraint in the production path. `filter.filter_extra` already made
the clamp conditional on `parameter.genome_gc`, which `get_params` does set from the
foreground genome. What blocked the set there was the GC *content* floor. Both the
clamp and the window are now conditional, on shared thresholds.

## What changed

`parameter.EXTREME_AT_GENOME_GC` (0.30) and `EXTREME_GC_GENOME_GC` (0.70) are now
the single definition of a compositionally extreme target, read by
`parameter.adaptive_gc_window`, `filter.filter_extra` and both filters in
`adaptive_filters` — the three had each carried their own literal, which is how
they came to disagree. Beyond a threshold the bound on the *matching* side is
released rather than centred on the genome mean:

| Target | genome GC | primer GC window |
|---|---|---|
| *P. falciparum* | 0.19 | **0.00**–0.34 |
| *B. burgdorferi* | 0.28 | **0.00**–0.43 |
| *Francisella* | 0.33 | 0.18–0.48 |
| *E. coli* | 0.50 | 0.35–0.65 |
| *Burkholderia* | 0.67 | 0.52–0.82 |
| GC-rich | 0.72 | 0.57–**1.00** |

The rationale is that SWGA selects on differential target/background k-mer
frequency, which on an AT-rich target drives primers to pure A/T — and short
8–12mers quantise GC coarsely enough that "near the genome mean" is not reachable
anyway.

### Scope

Deliberately narrow. Across all 31 published sets, exactly two were rejected
outright and both are very AT-rich targets — this one and `B31-BL21`
(*B. burgdorferi*, 0.00 mean GC, 0/20 before). Every other set already passed at
n−1 of n or better, and GC-normal behaviour is unchanged (a 0.50-GC target still
gets 0.35–0.65, clamp minimum still 1). Tests pin both directions.

### Necessary but not sufficient

`filter_extra` now keeps **3 of 10**, not 10. The rest are held back by rules this
change did not touch:

- the **homopolymer** rule rejects the twelve-adenine primer;
- the **self-dimer** rule rejects the six pure (AT)n primers. These are genuinely
  self-complementary — but at 30 °C a 10 bp AT duplex melts well below the reaction
  temperature, so counting paired bases without a thermodynamic threshold
  over-rejects in this regime.

Whether to relax either for AT-rich targets is a separate question with its own
evidence, not decided here. The boundary is pinned by test so that answering it
later is a visible change.

---

# Synthesis across the set-level benchmarks

| Dataset | n | density predicts? | evenness predicts? |
|---|---|---|---|
| Prevotella (Dwivedi-Yu 2023) | 6 | no | **yes** (with worst gap) |
| Wolbachia (Clarke 2017) | 5 | no (backwards) | **yes** |
| *M. tuberculosis* (Clarke 2017) | 12 | **yes** | no (backwards) |
| Leishmania (Dwivedi-Yu 2023) | 4 | not measurable | not measurable |
| Borrelia / Wolbachia (Leichty 2014) | 3 | — | — (not a design criterion) |
| *P. falciparum* (Oyola 2016) | 1 (of 3 nested) | — | — (only the winner published) |

Among the three where both are measurable, two favour evenness. The odd one out is
*M. tuberculosis*, and the paper states why: there the primer **pool** was
pre-filtered for even binding, which removed most of the variance in evenness and
left density as the limiting property.

All three are consistent with a single rule: **whichever property is currently
limiting is the one that predicts.** That is what the reporting change is built
around — `mean_gap`, `max_gap` and `gap_gini` are now shown together with their
benchmark context in both HTML reports and in `evaluate-set`, rather than folded
into a composite that would have to commit to one regime.

The two datasets added later sharpen rather than change this:

- **Leishmania** rules out primer composition as the explanation in at least one
  case, which pushes the signal further toward binding geometry.
- **Leichty 2014** shows the objective itself is not single-valued. `SR2` is
  simultaneously the most selective set on record and the least even. Any score
  that ranks sets on one axis is silently choosing for the user.
- **Oyola 2016** moves the question upstream of scoring entirely. Before any set
  can be ranked, its primers have to survive filtering — and for very AT-rich
  targets ours do not. A candidate that is never generated cannot be scored well
  or badly.

Two of the seven datasets (Leichty's `SR2` ⊂ `SR1`, Oyola's `Probe_10` ⊂
`Probe_20` ⊂ `Probe_28`) also independently found the **smaller** of two nested
pools to perform better. `target_set_size` defaults to 6 and `--auto-size` scales
with application profile; neither dataset contradicts those defaults, but both
argue against treating a larger set as safer.

## Why the scoring weights were still not changed

With 6, 5, 12 and 4 sets, any weighting fitted to one dataset is contradicted by
another. A `max_gap` term helps on Prevotella, is neutral-to-harmful on
*M. tuberculosis*, and does not cleanly separate the Wolbachia sets either
(`TmL/Selective` has a *better* max gap than the winner). The Tm-uniformity term
already in the score runs backwards on Leishmania. The composite score is left
alone; the statistics are surfaced so the reader can see the regime they are in.

## Datasets sought but not obtained

- **Cowell et al. (2017)** *Plasmodium vivax* `pvset1` / `pvset1920`. Sequences are
  recoverable, but the published outcome is a per-sample enrichment series rather
  than a comparison between sets, so it cannot separate a winner from a loser.

- **Thesis-derived sets** (*P. falciparum* 6A/8A; *Toxoplasma gondii*). The
  *T. gondii* sets are `swga` tool output with no wet-lab result, so they measure
  the tool's own agreement with itself rather than an outcome. Excluded on that
  basis.

## Survey method

The starting point was a local collection of 19 SWGA papers, screened
programmatically by scoring each PDF for the density of DNA-like tokens alongside
outcome keywords (fold, enrichment, coverage). That flagged five candidates, of
which one — Leichty and Brisson (2014) — carried both complete sequences and
quantified outcomes. The remaining benchmark, Leishmania, came from a subsequent
literature search rather than the local set.

The screen's value was mostly negative: most SWGA papers either report sequences
without a set-level outcome, or an outcome without sequences. Requiring both is
what keeps the suite at six datasets.
