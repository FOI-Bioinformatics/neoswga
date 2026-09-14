# Pipeline candidate retention and additive-aware optimization

Follow-up audit, 2026-09-14, commit `8c2428c`. This extends the
[additive/coverage audit](README.md). Production behavior was not changed.

## Answers

1. QC-passing candidates are **not all retained until optimization**. There are
   single-primer distribution filters, ranking shortlists, a hard candidate cap,
   and additional optimizer-dependent screens.
2. Additives influence candidate admission and ranking, and the final metrics.
   They do **not consistently influence the actual coverage search objective**.
3. Longer oligos plus additives are a plausible way to explore specificity,
   but the current implementation is not sufficiently consistent or calibrated
   to establish that such a design will improve experimental stringency.

## Stage-by-stage trace

| Stage | Candidate handling | Additive influence | Audit assessment |
|---|---|---|---|
| Initialization | Selects k-mer lengths, polymerase defaults, GC/Tm windows and optional long-primer settings. Explicit settings can override defaults. | Condition presets and adaptive settings affect which candidates are even considered. | A comparison must record resolved settings, not just CLI differences. |
| `count-kmers` | Enumerates/counts motifs at configured lengths for target/background references. | None on the sequence counts themselves. | Longer lengths not counted here cannot appear later. Existing record-boundary issues remain relevant to downstream position scans. |
| Candidate loading | Reads count files; rejects ambiguous bases, applies a broad GC prefilter and an additive-aware Tm window with margin. | Uses `conditions.calculate_effective_tm`. | This is already a condition-specific candidate universe. Its broad GC prefilter is distinct from subsequent GC QC. |
| Frequency filtering | Applies foreground/background frequency thresholds. For k > 10, thresholds scale by `4 ** (10-k)`. | Exact-count gates do not use additives. | Same JSON thresholds are not equivalent site-count thresholds across lengths. Approximate background-index modes need separate resolution checks. |
| Sequence QC | Applies effective Tm, GC range, homopolymers, repeats, GC clamp, end-composition and self-dimer checks. Optional exclusion/blacklist filters also run. | Tm includes additives. Sequence-complementarity thresholds are not a chemical dimer-equilibrium calculation. | Changing chemistry requires re-evaluating condition-dependent QC. The global primer-concentration propagation bug remains. |
| Distribution filter | Calculates foreground site-spacing Gini; rejects missing values and values above `max_gini`. | No direct additive term. | A primer useful only in a difficult region can be lost before its contribution to a pool is considered. |
| Ranking and cap | Sorts by exact background/foreground ratio and foreground abundance to form an occupancy shortlist; condition-dependent ranking runs within it; then only `max_primer` remain. | Occupancy ranking incorporates corrected Tm and mismatch-class loads. | A quality/ranking cut is treated as permanent exclusion. No guarantee exists that it preserves the best pool. |
| Background position indexes | Generated for the capped pool, not all Gini survivors. | None on positions. | The capped-away pool cannot simply be supplied later without indexing its background sites too. |
| `score` | Default carries stage-2 candidates forward; `--amp-model` restores the legacy efficacy score/filter. | Default adds no new chemistry-based selection. | Stage 3 is not the main loss point. The opt-in efficacy gate can reject additional candidates. Sampling settings in that legacy path are not evidence of sampling the default candidate inventory. |
| Standard `optimize` entry | Loads stage 3; checks indexes; normally applies a background prefilter unless disabled. | Prefilter uses exact counts, `fg/(bg+1)`, not corrected occupancy or per-base density. | Default ratio threshold 1 and removal cap 20% add another possible loss after stage 3. |
| Hybrid/background-aware prescreen | Polymerase-dependent secondary Tm/GC/structure screen; may remove primers with many pool-wide dimer conflicts. | Inconsistent Tm handling, described below. | A candidate's compatibility with the eventual small pool differs from its number of conflicts with every candidate. If no candidate passes, this screen returns the unfiltered pool, so it is not a uniformly enforced hard QC gate. |
| Coverage-greedy selection | Chooses compatible primers that cover the most new geometric regions. | No per-primer occupancy term; no full conditions object passed to the dominating-set engine. | Chemistry can change the input pool/order, but is not priced directly in marginal coverage gains. |
| Network refinement | Scores network structure; optional Tm and mechanistic terms can affect choices. Background-aware pruning uses background costs. | Some condition dependence, with a separate inconsistent Tm path. | Optional terms do not repair the mismatch between coverage selection and the reported effective coverage metric. |
| Swap refinement | Lexicographic optimization: raw covered bases first, exact background-site reduction second, under dimer compatibility. | No corrected occupancy in these two objectives. | A coverage improvement may increase background; background only resolves equal-coverage alternatives. |
| Final evaluation/report | Recomputes actual panel metrics. `plan-pool` checks effective coverage, specificity and dimer constraints and recommends the smallest qualifying found. | Effective coverage/selectivity respond to additives. | These are acceptance checks after search. They do not make the search itself condition-aware or constrained by the requested density threshold. |

`plan-pool` uses a distinct entry path: it selects one oligo length, validates
sequences, removes self-dimers, then constructs an optimizer directly. It skips
the standard `optimize_step4` background prefilter. It still inherits the
earlier `max_primer` cut and any internal optimizer screens. It fixes chemistry
for each run; it does not jointly optimize oligo length and additive doses.

## Measured Wolbachia funnel

Saved `work/filter_stats.json`, confirmed against the two candidate CSVs:

| Point | Candidates remaining |
|---|---:|
| Loaded after early candidate prefilters | 874,596 |
| Foreground frequency | 874,596 |
| Background frequency | 706,756 |
| Sequence/thermodynamic QC | 491,836 |
| Gini filter | 20,670 |
| Ranking cap | 2,000 |
| Default stage 3 | 2,000 |

Stage-2 and stage-3 primer sets are identical in this example. The cap removes
18,670 of the 20,670 Gini survivors, about 90.3%. A later optimizer cannot
select them. This does not prove that those candidates improve coverage; it
does mean the saved result cannot establish what is possible with all survivors.
The funnel's `total_kmers` label is also misleading: it starts after early
loading prefilters, not at every motif initially counted.

## New findings about chemistry during optimization

### P1: secondary Tm screening can undo additive-aware QC

`thermodynamic_filter.py:189` calls `calculate_tm_with_salt` using Na and Mg,
then applies its Tm limits before building conditions for the structure checks.
It does not apply additives to that Tm value. `hybrid_thermo_screen.py` invokes
this screen for configured polymerases and uses polymerase GC/temperature
criteria, not necessarily the full condition-specific settings already used by
stage 2. In particular, a longer primer made acceptable by an additive-induced
Tm shift can be rejected again using its uncorrected Tm.

### P1: network scoring uses yet another Tm

`network_optimizer.py:746` forwards Na concentration to `calculate_tm_with_salt`
but leaves Mg at the low-level default of zero, then adds the additive shift.
It also omits K, NH4, dNTP correction and explicit primer concentration. The
standard condition object's calculation includes additional buffer species.

Numerical probe for `ATCGATCGATCG`, Na 50 mM, Mg 10 mM, temperature 30 C,
DMSO 5%, betaine 1 M (a code diagnostic, not a proposed recipe):

| Calculation | Tm |
|---|---:|
| Main condition-aware calculation | 43.85 C |
| Secondary screen | 47.72 C |
| Network scoring | 34.92 C |

These differences can change filtering or rankings without any experimental
condition changing. They should be fixed before interpreting differences between
long-oligo/additive designs.

### P1: selection and acceptance use different objectives

The default `plan-pool` background-aware/swap path optimizes raw reach and exact
background sites, then evaluates occupancy-weighted coverage and density.
Requested minimum density is not passed into the inner search as a constraint.
Therefore `no qualifying panel found` may reflect a search-objective mismatch,
not an infeasible design. Conversely, a high raw-coverage panel can pass one
search stage while having poor binding at the chosen reaction temperature.

## Longer oligos, additives and stringency

Sequence specificity and reaction stringency are related but distinct.
Increasing length often reduces exact off-target occurrences, while also
changing target site density and duplex stability. Longer oligos can have more
stable near-matched off-target duplexes too. Additives that lower Tm can shift
the operating point so weaker duplexes bind less while target binding remains
adequate. If both target and off-target duplexes are strongly bound, a small
shift may change little; if target binding is already marginal, the same shift
may mainly reduce useful priming.

The mechanism is plausible, with direct PCR evidence for using DMSO and Tm to
adjust specificity ([Chester & Marshak 1993](https://pubmed.ncbi.nlm.nih.gov/8470801/)).
PCR evidence does not validate a specific SWGA dose, oligo length or predicted
recovery improvement. Betaine also affects GC dependence and secondary
structure, so it is not simply a universal specificity enhancer. The earlier
audit documents the limitations of the coefficients used here.

Within the current model, an exact match for the same oligo receives the same
occupancy in target and background. An additive cannot distinguish those two
sites simply from their organism labels. Relative selectivity changes through
mismatch classes or changes in which primers dominate the summed load. The
fixed mismatch-count penalty also ignores mismatch identity, position and
extension competence. A better density ratio can coincide with reduced useful
target binding and does not imply improved genome recovery.

Longer, more unique oligos may have too few target sites to cover the genome
with a small pool. Adding more unique oligos and using fewer repeated oligos
are different tradeoffs. The appropriate comparison is a pool-size/coverage/
specificity frontier across lengths and conditions, with per-oligo or total
pool concentration explicitly fixed and actual panel metrics recomputed.

## Recommended design changes

1. **Preserve the candidate inventory.** Store every enumerated candidate with
   QC outcomes/reasons and metrics. Retain every candidate passing the declared
   hard constraints in a reusable inventory. Label any search shortlist as a
   computational limit, not QC rejection. Make it possible to expand that
   shortlist as uncovered regions or dimer conflicts are encountered.
2. **Move pool properties into pool optimization.** Single-primer Gini, abundance
   and pool-wide conflict counts are useful ranking features, but need not be
   permanent exclusion gates. A localized primer may complement the rest of
   a pool. Keep explicit biological exclusion/background limits enforceable.
3. **Unify chemistry.** Use the same resolved conditions, concentration and Tm
   implementation in loading, QC, secondary screening, selection and reporting.
   Cache entries need condition identity when conditions can change.
4. **Match selection to acceptance.** Optimize marginal effective coverage with
   background limits enforced during search. Keep dimer constraints on the
   selected panel. Report the achieved tradeoff rather than only a final
   pass/fail outcome from a differently optimized panel.
5. **Compare conditions from the reusable inventory.** Recompute condition-based
   QC and scores for every chemistry/length comparison. If any admissible tested
   condition permits a candidate, preserve it in the shared inventory, while
   allowing its use only in panels whose own conditions pass QC. Do not reuse
   one chemistry's already-capped CSV as the universe for all other chemistries.
6. **Validate promising frontiers experimentally.** Resolve the earlier coverage
   and coefficient issues before claiming improved stringency or recovery.

The immediate priority is retaining useful candidates and making the search
use consistent chemistry. Increasing additive levels alone would not fix these
software limitations.

## Verification and entry points

- Candidate cuts: `core/pipeline.py`, `_rank_and_cut_candidates`,
  `_rank_by_occupancy`, `step2`, `step3`.
- Frequency scaling: `core/filter.py`, `_scale_freq_threshold`.
- Default prefilter: `core/unified_optimizer.py`, `_prefilter_by_background`.
- Secondary screen: `core/hybrid_thermo_screen.py`, `core/thermodynamic_filter.py`.
- Greedy coverage: `core/dominating_set_optimizer.py`, `_select_next_primer`.
- Swap objective: `core/swap_refinement.py`, `refine_hybrid_stage2`.
- Final acceptance: `core/pool_planner.py`, `plan_pool`.

Numerical Tm comparisons and the saved candidate counts were checked directly
against the current implementation. No experimental comparison of longer
oligos/additives was performed, and no production code was changed.

The 60 targeted tests for candidate ranking, loading, stage-3 retention, swap
refinement and background-prefilter configuration passed. Four warnings came
from synthetic non-DNA candidate labels in existing tests.

Reproduce the numerical checks with
`PYTHONPATH=. python docs/validation/additives_coverage_audit_2026-09-14/pipeline_probes.py`.
The script writes [pipeline_probe_results.json](pipeline_probe_results.json).
