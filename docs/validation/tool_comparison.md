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
| Coverage reach | — | 70 kb (phi29 processivity) | — | **3 kb, fitted to a measured outcome (3.0-6.2 kb band)**; configurable, also fittable from BAM depth; reported at 3/10/70 kb |
| Heterodimer-free set | guaranteed by clique | pairwise check while growing | — | guaranteed by `--optimization-method clique`; penalised in the other four |
| Refine an existing set | no | no | yes | yes (`expand-primers`, `analyze-coverage --bam`) |
| Multi-genome / pan-target | no | no | no | yes |
| Wet-lab validation suite | — | — | — | yes, seven published datasets |

## Two deliberate divergences

These are positions, not gaps. They are recorded here so the disagreement is explicit
rather than silent.

**Coverage reach: 3 kb by default, not 70 kb — and now configurable and
measurable.** swga 2.0's `coverage_ratio` uses phi29's single-molecule processivity of
~70 kb as the per-primer reach. NeoSWGA separates the two quantities:
`processivity_bp` (~70 kb) governs amplification-*network* connectivity, while
`coverage_reach` (~3 kb, `coverage.resolve_coverage_reach`) governs set-cover selection
and reported `fg_coverage`. One molecule can be extended 70 kb; that is not the reach at
which a primer reliably covers a genome in a multi-primer reaction. Conflating them
overstates coverage by more than an order of magnitude.

That position stands, and one half of it is now measured -- see
[reach_calibration.md](reach_calibration.md). Fitting the reach to Clarke et
al.'s (2017) *Wolbachia* outcome, the one published set with a quantitative
breadth figure, puts it at **3.0-6.2 kb**. The 3 kb selection default sits at
the bottom of that range and is corroborated; the 10 kb reporting reach is
outside it. What follows was written before that measurement and is left
because the reasoning still holds for the reporting convention:

**Neither convention was measured.** 70 kb is what one molecule *can* do, not what a product reaches before a
neighbour's strand displacement truncates it; and the ~3 kb default is derived from the
inter-primer *spacings* Clarke et al. (2017) and Dwivedi-Yu et al. (2023) designed for,
which is a designer's choice rather than an observation. The difference is not academic
— the same set measures:

| reach | coverage |
|---|---|
| 3 kb | 0.418 |
| 10 kb | 0.836 |
| 70 kb | ~1.0 |

Three changes follow. The reach is **configurable** (`coverage_reach` in params.json,
`--coverage-reach`), because the right value is a property of the reaction rather than
of the enzyme alone. Coverage is **reported at 3 / 10 / 70 kb** plus whichever reach was
used for selection (`coverage_by_reach`), so a NeoSWGA figure can be placed beside a
published one without either convention having to move. And `neoswga calibrate-reach
--bam` **fits the reach to real sequencing depth**: if products reach R bases, depth
falls away from binding sites on that length scale. It reports the whole grid and marks
the result uninformative below a correlation floor, because a library with no
primer-position signal would otherwise still yield a confident-looking number. No other
SWGA tool estimates this from data.

**No fitted set-level score.** swga 2.0 fits five set-level terms by ridge regression
against measured 1x coverage (`freq_ratio` +0.321, `off_gap_gini` +0.281,
`mean_gap_ratio` -0.0368, `coverage_ratio` -0.0318, `on_gap_gini` -0.0131). NeoSWGA
uses hand-chosen weights instead. This is not an oversight: the validation work in
[published_primer_sets.md](published_primer_sets.md) covers seven datasets, which is
far too small a sample to fit weights on without overfitting to them. Hand-chosen
weights that can be argued from mechanism are more honest at this sample size than
fitted weights that cannot be validated out of sample.

## What the 2026-08-02 audit changed

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

**PEG, BSA and DTT act on the enzyme, and now do so.** All three were accepted,
range-validated and consumed by nothing. They were given mechanistic pathways rather
than Tm coefficients, because none is a duplex-stability agent at these
concentrations: PEG crowding raises kon up to a plateau and then costs speed through
viscosity (so the model reports an interior optimum rather than "more is better"),
BSA is a saturating stabiliser, DTT a deficiency penalty. This also required folding
`stability_factor` into `predicted_amplification_factor` — it was computed, reported,
and left out of the product, so every stability-mediated term was invisible in the
headline number.

**A third coverage semantics was found and not adopted.**
`minimal_primer_selector.MinimalPrimerSelector` — never called from anywhere — counts
`covered_positions` as the set of binding-site *coordinates*, so a primer covers as
many bases as it has sites and extension is ignored. On a 30 kb genome, 30 sites reads
as 0.1% coverage. The minimisation pass uses `optimizer.compute_metrics` instead, so
selection, trimming and reporting all share the figure written to
`step4_improved_df_summary.json`.

## What the 2026-08-19 audit changed

**The reach is calibrated.** Previously an assumption on both sides; now fitted
to a measured outcome at 3.0-6.2 kb, with the site-finding validated against
the paper's own statistics for three of five sets. Full write-up and the two
sets that do not reconcile: [reach_calibration.md](reach_calibration.md). No
other SWGA tool publishes a fitted reach, and swga 2.0's `coverage_ratio` uses
70 kb -- at which almost any set scores ~1.0, so a comparison at that reach
carries little information.

**Coverage was six quantities.** Three of them were compared against a 0.95
threshold as if they meant the same thing. `DominatingSetOptimizer` divided
covered bins by the bins the *candidate pool* could reach, so a design covering
a tenth of the target scored 1.000 and reported SUCCESS; that number was also
its `score` and its stop condition. `expand-primers` measured bin occupancy
rather than coverage, its value set by `bin_size` rather than the polymerase.
Neither has a counterpart in swga 1.0 or 2.0, which compute one coverage
definition each.

**The set-size ceiling was below the measured requirement.** swga 1.0 and 2.0
impose no cap; NeoSWGA's schema capped `num_primers` at 50 while a >95% design
on a 3.2 Mb target needs 31-96 primers inside the calibrated reach band. Raised
to 200. The `--auto-size` ceiling of 20 is deliberately not raised: the estimator
behind it is genome-independent (`genome_length` cancels out), so at every SWGA
primer length the clamp is the whole answer.

**The quality grades were uncalibrated and self-inconsistent.** `report` and
`interpret` carried separate threshold tables that had drifted -- enrichment
"excellent" was 500x in one and 200x in the other -- and both were set above
anything ever published. Every set in the seven-dataset benchmark suite,
including the 120-fold *Prevotella* winner, graded "poor" or "critical" on
binding evenness. Both now read one table calibrated against that suite, and the
evenness grade is documented as descriptive rather than predictive, consistent
with the position below on not folding these statistics into a single score.

Neither swga 1.0 nor swga 2.0 grades a design at all, so this is a NeoSWGA-only
surface and correspondingly the one with least external calibration.

## Gaps not closed

- **`off_gap_gini` is still not computed**, and the reach calibration does not
  change that: it carries the second-largest weight in swga 2.0's fitted model
  and is the one predictor NeoSWGA cannot produce. Reporting it is worthwhile;
  folding it into `normalized_score` is not, on the sample size available.
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
- **Additive coverage is broader than the published tools', and correspondingly less
  validated.** swga 1.0 and 2.0 model no additives at all. NeoSWGA models thirteen
  plus four buffer species, but the Tm coefficients are literature-derived while the
  enzyme- and kinetics-pathway magnitudes are EMPIRICAL, calibrated to reproduce
  qualitative protocol behaviour rather than measured. `docs/SCIENCE_CITATIONS.md`
  marks which is which. Predictions of *relative* effect are the defensible use;
  absolute yield needs calibration against local wet-lab data.
