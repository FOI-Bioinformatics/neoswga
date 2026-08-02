# NeoSWGA compared with published SWGA design tools

**Date:** 2026-08-02
**Scope:** capability comparison against the three published SWGA primer-set
designers, and an audit of whether NeoSWGA's own design options change the design.

Related: [published_primer_sets.md](published_primer_sets.md), which validates the
scoring against wet-lab outcomes rather than against other tools.

## The tools

- **swga 1.0** — Clarke et al. 2017, *Bioinformatics*. The original. Selects sets by
  finding a maximum clique in a primer compatibility graph, using the `cliquer` C
  library, with vertices weighted by mean background binding distance.
- **swga 2.0 / soapswga** — Dwivedi-Yu et al. 2023, *PLOS Comput Biol* 19(4):e1010137.
  Adds a random forest primer filter trained on 396 primers from rolling-circle
  amplification experiments, and a set-level score fitted by ridge regression against
  observed 1x coverage. Roughly 12x faster than swga 1.0.
- **COATswga** — bioRxiv 2025.11.26.688640. Interval-based tiling for coverage
  uniformity, with gap filling and refinement of a pre-existing set.

## Comparison

| | swga 1.0 | swga 2.0 | COATswga | NeoSWGA |
|---|---|---|---|---|
| K-mer counting | DSK (disk-based, low memory) | jellyfish + h5py cache | k-mer counting | jellyfish (hard dependency) + Bloom filter |
| Set formation | max-clique on compatibility graph | breadth-first greedy with drop-out | interval tiling | greedy set cover + network refinement; **clique also available** |
| Set scoring | binding metrics | ridge regression fitted to observed coverage | coverage uniformity | `normalized_score`, hand-chosen weights |
| Primer scoring | filters only | RF, 1500 trees / depth 50 | ML | RF, 100 trees / depth 15 |
| Coverage reach | — | 70 kb (phi29 processivity) | — | 3 kb (realistic per-primer) |
| Heterodimer-free set | guaranteed by clique | pairwise check while growing | — | guaranteed by `--optimization-method clique`; penalised in the other four |
| Refine an existing set | no | no | yes | yes (`expand-primers`, `analyze-coverage --bam`) |
| Multi-genome / pan-target | no | no | no | yes |
| Wet-lab validation suite | — | — | — | yes, seven published datasets |

## Two deliberate divergences

These are positions, not gaps. They are recorded here so the disagreement is explicit
rather than silent.

**Coverage reach: 3 kb, not 70 kb.** swga 2.0's `coverage_ratio` uses phi29's
single-molecule processivity of ~70 kb as the per-primer reach. NeoSWGA separates the
two quantities: `processivity_bp` (~70 kb) governs amplification-*network*
connectivity, while `coverage_reach` (~3 kb, `coverage.polymerase_extension_reach`)
governs set-cover selection and reported `fg_coverage`. One molecule can be extended
70 kb; that is not the reach at which a primer reliably covers a genome in a
multi-primer reaction. Conflating them overstates coverage by more than an order of
magnitude.

**No fitted set-level score.** swga 2.0 fits five set-level terms by ridge regression
against measured 1x coverage (`freq_ratio` +0.321, `off_gap_gini` +0.281,
`mean_gap_ratio` -0.0368, `coverage_ratio` -0.0318, `on_gap_gini` -0.0131). NeoSWGA
uses hand-chosen weights instead. This is not an oversight: the validation work in
[published_primer_sets.md](published_primer_sets.md) covers seven datasets, which is
far too small a sample to fit weights on without overfitting to them. Hand-chosen
weights that can be argued from mechanism are more honest at this sample size than
fitted weights that cannot be validated out of sample.

## What this audit changed

**Coverage was binned more coarsely than the polymerase reaches.**
`BipartiteGraph.add_primer_coverage` marks a bin covered when any part of it is within
reach of a binding site, which is only sound while a bin is no larger than the reach.
At the shipped defaults it was not — 10 kb bins against the realistic ~3 kb reach — so
a primer covering 30% of a bin claimed all of it. Two consequences, both silent:
Stage-1 greedy set cover saw the genome covered after two primers and stopped, making
`--num-primers` inoperative on the default path; and the coverage reported for that set
was 1.000 where the honest figure was 0.433. `coverage_bin_size()` now caps the bin at
a quarter of the reach, bounding the overstatement at 25%. Coverage is now a property
of the genome rather than of the bin size, and at 100 kb and 1 Mb the requested count
is honoured exactly.

**A clique optimizer is available again.** `--optimization-method clique` builds the
swga 1.0 compatibility graph — an edge joins two primers only when they do not
dimerise — and selects a clique, so a dimer-free set follows from the representation
rather than from a threshold. The four other methods penalise dimers but cannot
guarantee their absence, because a greedy search can always judge the coverage worth
the penalty. The module had been deleted in commit 08044d9 for being unused rather
than for being wrong; two defects were fixed on the way back in (an enumeration cap
that could not stop a single large clique from expanding to millions of subsets, and
a selection objective that disagreed with the metric the result was then judged on).
It is limited to pools of roughly 200 candidates and is not in the default ensemble.

**Five options did nothing.** `--use-cooperative-binding`, `--primer-strategy` and
`--scoring-weights` were accepted, logged and stored, and read by no module under
`neoswga/core/`; all three were removed (`--scoring-weights` duplicated
`--application`, which is genuinely wired). `--max-extension` was ignored in favour
of a hardcoded 70000 and is now wired. `--minimize-primers` and `--target-coverage`
were swallowed by `**kwargs` — read off an `OptimizerConfig` that has no such fields —
and now drive a real post-processing pass.

The first attempt to guard this was a grep for the option name under `neoswga/core/`,
and it was too weak: a name can appear as a config field that is set, forwarded and
never read, or only as a substring of an unrelated identifier
(`target_coverage` vs `per_target_coverage`). Stored is not consumed.
`tests/test_design_options_have_effect.py` now runs the dispatch path twice with the
option varied and requires the selected set or its score to differ, which is the only
check that separates a working option from decoration.

One near-miss worth recording: `--uniformity-weight` looks dead by the same test and
is not. It is consumed by `NetworkOptimizer._evaluate_primer_addition`, but on a
genome without strong spacing differences the uniformity term clamps to the same value
for every candidate, so the base score breaks ties in the same order whatever the
weight. An end-to-end assertion would have failed for the wrong reason and condemned a
working option.

**Reaction conditions now reach the chemistry.** `ReactionConditions` takes 20
parameters and nothing constructed it with all of them — each of nineteen call sites
hand-listed its own subset, from two of twenty in `cli/evaluate.py` to eleven in
`core/filter.py`. No site passed `k_conc`, `nh4_conc`, `dntp_conc` or `dtt_mm`, so the
phi29 buffer species reached no calculation at all. `build_reaction_conditions()` now
reads the field list off the constructor, so a new field is forwarded everywhere
without anyone updating call sites. Effects that previously reached nothing: K⁺ +1.07 °C,
NH₄⁺ +0.46 °C, dNTP −1.16 °C (Mg²⁺ chelation), glycerol changing predicted amplification,
SSB doubling the binding rate. The last two act through `MechanisticModel` rather than
Tm, and the `--auto-size` path that builds that model was among the sites dropping them.

**A third coverage semantics was found and not adopted.**
`minimal_primer_selector.MinimalPrimerSelector` — never called from anywhere — counts
`covered_positions` as the set of binding-site *coordinates*, so a primer covers as
many bases as it has sites and extension is ignored. On a 30 kb genome, 30 sites reads
as 0.1% coverage. The minimisation pass uses `optimizer.compute_metrics` instead, so
selection, trimming and reporting all share the figure written to
`step4_improved_df_summary.json`.

## Gaps not closed

- **No interval tiling.** COATswga's uniformity approach has no counterpart here.
  `tiling_optimizer.py` was deleted in the same commit as the clique optimizer and
  could be restored the same way.
- **No background gap evenness.** `off_gap_gini` carries the second-largest weight in
  swga 2.0's fitted model and is not computed anywhere in NeoSWGA. `PrimerSetMetrics`
  has `gap_gini` for the foreground only.
- **jellyfish is a hard dependency.** There is no in-process counting fallback;
  `kmer_counter.count_kmers_in_sequence` exists but has no callers.
- **RF training provenance is undocumented.** swga 2.0 states its training set (396
  primers from RCA experiments); the model shipped in
  `neoswga/core/models/random_forest_filter.skops` does not have an equivalent record.
- **`peg_percent`, `bsa_ug_ml` and `dtt_mm` have no model behind them.** They are
  accepted, range-validated and now correctly forwarded, but nothing consumes them:
  no Tm coefficient and no mechanistic pathway. This was previously masked by the
  forwarding bug — they could not be seen to do nothing while they were not arriving.
  Either give them a pathway or drop them; PEG is a crowding agent and BSA a
  stabiliser, so neither belongs in a Tm correction, and a mechanistic term would be
  the honest place. Note `glycerol_percent` and `ssb` are in the same shape but *do*
  have pathways (enzyme stability/speed, and kon respectively).
