# NeoSWGA compared with SWGA, SOAPswga and COATswga

Review date: 22 September 2026. Verdict: NeoSWGA has a useful, comparatively broad
design and diagnosis framework, but the inspected evidence does **not** establish
that it produces the best experimentally performing pools. Its immediate priority
is connecting and validating existing safeguards, not adding more optimizers.

## Scope and evidence

Inspected snapshots:

| Application | Revision |
|---|---|
| NeoSWGA | `1c4f20cbd2d7a396fb2b61667c7c60b08bd2d6d9` |
| Original SWGA | [`b36cefa`](https://github.com/eclarke/swga/tree/b36cefaeb50d595e34068a86f466c3c9fde5aab7) |
| SOAPswga / swga2.0 | [`c25667e`](https://github.com/songlab-cal/swga2/tree/c25667ef87a856a923bde5bea353fe2608c59492) |
| COATswga | [`e2e0ff4`](https://github.com/bailey-lab/coatswga/tree/e2e0ff4dae1957f24908bdebafa2d32763f3b996) |

The working tree was clean at the start. The no-fallback implementation is merged;
the previous audit's unfinished-code and size failures are not assumed current.
Competitor repositories were downloaded for static inspection. Their programs were
not installed or benchmarked against NeoSWGA, and no new amplification experiments
were performed. Therefore this review makes no measured cross-tool speed or
biological-superiority claim. No production implementation was changed.

Fresh NeoSWGA verification: **325 focused tests passed**, including request,
failure, model-evidence, panel-assessment, budget, feedback, shared-contract and
size checks. **13 strict-pipeline integration tests passed**. These include real
command execution and refusal/export checks. **101 published-data validation tests
passed**. Total: **439 passed**, with no skips in these runs. This is targeted
verification, not a full-suite rerun or a prospective biological validation.
The commands and results are saved in the
[verification record](2026-09-22-comparison/verification.json).

The read-only [probe script](2026-09-22-comparison/probe_contracts.py) and its
[recorded output](2026-09-22-comparison/probe_results.json) reproduce the contract
findings. Named-call inspection was combined with manual tracing; absence of an
AST call alone is not proof against every possible dynamic invocation.

## What the external evidence supports

**Original SWGA has directly relevant experimental evidence.** Clarke et al.
tested designs for Wolbachia against Drosophila and for M. tuberculosis against
human material. That makes the original tool an essential baseline for this
project's Wolbachia objective, despite its older implementation. The paper does
not establish that any arbitrary compatible primer set will work.
[Original study](https://pmc.ncbi.nlm.nih.gov/articles/PMC5870857/).

**SOAPswga and swga2.0 are names for the same linked tool.** Its study trained an
individual-primer model on experimental amplification measurements and describes
greater usefulness for rejecting poor primers than accurately separating the
best ones. It reports experimental Prevotella results and a substantial runtime
improvement over original SWGA in its own benchmark. Those results are evidence
for that tested setting, not a present-day runtime comparison with NeoSWGA.
[Paper](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010137),
[repository](https://github.com/songlab-cal/swga2).

**COATswga provides a useful coverage-oriented comparator.** Its repository
documents KMC counting, interval tiling, directional coverage, and existing-primer
support. Its associated work is identifiable as a November 2025 bioRxiv preprint.
The publisher's full-text endpoint was inaccessible in this review, so the
reported experimental superiority in the supplied summary was not independently
assessed against methods, replication and sequencing-depth controls. No confirmed
later peer-reviewed version was established here.
[Repository](https://github.com/bailey-lab/coatswga),
[publisher record](https://connect.biorxiv.org/qr/2025.11.26.688640),
[preprint DOI](https://doi.org/10.1101/2025.11.26.688640).

The summary supplied with this review needs one qualification: geometric tiling
does not guarantee uniform observed amplification. Likewise, a modelled coverage
fraction is not sequencing breadth at a stated depth. These distinctions apply
to all four tools.

## Comparative strengths and weaknesses

| Dimension | Original SWGA | SOAPswga / swga2.0 | COATswga | NeoSWGA at this revision |
|---|---|---|---|---|
| Main selection strategy | Compatible sets through clique search and binding statistics | Learned efficacy filtering and guided, bounded set search | Parallel starts with interval-based additions | Multiple proposal methods plus shared repair, swaps and contraction |
| Candidate management | Frequency/Tm/dimer filters and an active candidate set | Filtering and top-ranked candidate restriction | Frequency/ratio filters and ordered starts | Persisted QC inventory and refillable frontiers; strongest potential to avoid permanent shortlist loss |
| Spatial design | Binding spacing/evenness | Strand/gap features in set scoring | Explicit forward/reverse interval coverage | Geometry, weighted-coverage proxies, gap measures and per-target constraints |
| Chemistry | Standard Tm/buffer and sequence compatibility criteria | Thermodynamic features and learned efficacy | Tm filtering and nearest-neighbor dimer alignment | Broader additive/polymerase controls; evidence/domain and concentration integration remain incomplete |
| Number of oligos | Evaluates compatible sets under configured search | Configured iterative set construction | Coverage goal and minimum set size | Size curves, constrained deletion and fixed-primer expansion; no global minimum guarantee |
| Sequencing iteration | Not a central workflow in inspected code | Not a central workflow in inspected code | Existing-primer/gap-filling support | BAM-derived gaps, reach fitting and deficit-oriented expansion; full feedback provenance/calibration not complete |
| Experimental evidence | Published design-and-test study, including Wolbachia | Published experimental model and pool evaluation | Associated preprint; full methods not assessed here | Retrospective published-data fixtures and diagnostics; no prospective head-to-head evidence identified in reviewed artifacts |
| Engineering | Compact older architecture; Python 2 installation barrier | Reusable HDF5/caching; aging pinned dependencies/model serialization | Small focused Python application with external KMC/bedtools | Modern Python, broad tests, CI, strict references and failure records; considerably more integration surface |

Competitor algorithm/package observations come from the pinned
[SWGA implementation](https://github.com/eclarke/swga/tree/b36cefaeb50d595e34068a86f466c3c9fde5aab7/swga),
[SOAP implementation](https://github.com/songlab-cal/swga2/tree/c25667ef87a856a923bde5bea353fe2608c59492/src),
and [COAT implementation](https://github.com/bailey-lab/coatswga/tree/e2e0ff4dae1957f24908bdebafa2d32763f3b996/src).
An absent feature means not found in the inspected main path, not a claim that
users cannot implement it externally.

### What NeoSWGA does well

- `candidate_inventory.py` and `candidate_source.py` separate eligibility from a
  memory/search frontier. This is a useful improvement over treating the top-N
  candidates as the complete design universe.
- `pool_objective.py`, `optimization_service.py`, `panel_refinement.py` and
  `panel_contraction.py` provide reusable constrained improvement. Fixed/excluded
  primers and explicit background limits are more useful than a single opaque score.
- Reference identity/geometry checks, structured failures and blocked exports
  address a consequential scientific failure mode: a plausible result computed
  against missing or wrong reference data. The integration tests exercised these.
- `occupancy_coverage.py` uses interval events rather than whole-genome arrays
  for every primer. This is a sensible computational design, although its speed
  advantage over competitors has not been measured here.
- BAM-driven gaps and the existing blocked reach-calibration diagnostics offer a
  practical route from a first design to a better second design.
- Modern packaging and Linux/macOS/Python-version CI are advantages over original
  SWGA's explicit Python 3 rejection and SOAP's older dependency pins.
  [SWGA setup](https://github.com/eclarke/swga/blob/b36cefaeb50d595e34068a86f466c3c9fde5aab7/setup.py#L16),
  [SOAP requirements](https://github.com/songlab-cal/swga2/blob/c25667ef87a856a923bde5bea353fe2608c59492/requirements.txt).

### What should be learned from competitors

**From SWGA:** transparent binding/evenness criteria and experimentally tested
Wolbachia designs. Its age is not a reason to discard its biological baseline.
Clique compatibility is a modelled property, not a laboratory guarantee; a search
over prefiltered candidates is not exhaustive over every possible oligo pool.

**From SOAP:** experimental training data and strand/gap features. NeoSWGA's retired
synthetic-data amplification predictor is not equivalent to SOAP's experimentally
trained filter. Retiring the former does not demonstrate that a properly validated
efficacy model is useless. The inspected SOAP evaluator loads a serialized learned
set model and uses gap/strand/background features; benchmark it as shipped before
replacing its objective with NeoSWGA's own.
[SOAP evaluator](https://github.com/songlab-cal/swga2/blob/c25667ef87a856a923bde5bea353fe2608c59492/src/optimize.py#L926).

**From COAT:** explicit directional interval coverage and a small, focused design
path. However, `sets.py` implements ordered additions from multiple starts, not
an exhaustive search over all subsets. Its outer loop uses CPU count to determine
start batches, with a low stopping threshold when forced coverage is disabled;
resource settings can therefore alter the search performed. Existing primers are
included in the returned panel/dimer context while initial interval coverage is
built from the seed. That limits treating this path as a fully assessed incumbent
pool refinement. These are static observations requiring runtime fixtures before
quantifying their impact. Its minimum-size padding also does not prove a minimum
pool, despite the printed wording.
[COAT set construction](https://github.com/bailey-lab/coatswga/blob/e2e0ff4dae1957f24908bdebafa2d32763f3b996/src/sets.py).

## Priority findings in NeoSWGA

### P1: fixed-total concentration is a declared policy, not an evaluator policy

`design_request.py:185` calculates per-member concentrations, but there are no
production calls to `concentrations_molar`. `base_optimizer.py:1441` evaluates
weighted coverage with its existing reaction conditions. The probe's fixed total
of 12 micromolar across two oligos gives a declared allocation of 6 micromolar
each, while `request.conditions.primer_conc` remains 0.5 micromolar.

The concentration keys also do not appear in `params.schema.json`; accepting them
in the request resolver alone does not establish usable command-wide support.
Consequently, this feature must not yet be advertised as concentration-aware
pool-size optimization. Wire actual panel allocation into the evaluator, model
cache keys, filtering policy and report, or refuse the unsupported mode.

### P1: the new authoritative assessment is not authoritative in execution

`panel_evaluation.evaluate_panel` has no production caller. Its implementation
does not include effective coverage in the returned named metrics and delegates
violations only to an optional objective. It does not itself enforce dimer,
composition or requested size. A controlled one-primer assessment without an
objective qualifies under a request for twelve primers.

This is an interface/protection gap, not proof that the existing CLI currently
exports that panel: existing stages have separate acceptance checks. It does mean
the new record cannot yet be trusted as the sole export gate. Integrate a complete
assessment and remove duplicate acceptance logic, rather than simply calling the
current incomplete function.

### P1: recorded requests do not uniquely describe executed searches

`design_request_for_run` (`design_request.py:555`) reads the parameter file but
ignores command-line overrides. A probe supplying temperature 42 and seed 17 on
the args object records temperature 30 and no seed from the file. `plan-pool`
similarly resolves its request before applying command-specific choices.

Changing `bg_circular`, `stage1_objective_width`, or `swap_max_evaluations` produces
the same request hash in the probe. Nested `ReactionConditions` remains mutable:
changing temperature succeeds and changes the supposed immutable request's hash.
Explicit `target_set_size=0` is replaced by another size instead of rejected.

Build one effective request after all overrides, freeze nested data, validate
values without truthiness defaults and hash the settings/artifacts actually used.

### P1: evidence annotation does not yet enforce a quantitative support domain

`model_evidence.py` explicitly says the cited primary literature was not reread
when compiling the registry. This is honest and useful, but its temperature,
buffer and sequence domains are text. `require_model_support` checks polymerase,
enzyme-associated length range and whether an additive record is absent; it does
not establish mixture interactions, modification support or a numerical domain
test for each requested calculation. Enzyme-recommended primer lengths are also
not a substitute for the domain of a thermodynamic model.

Use machine-checkable domains and versioned parameter provenance. Keep supported
proxy calculations available with an evidence label, but fail unsupported required
calculations. Do not infer that additives make longer oligos more selective merely
because they change computed Tm.

### P2: budget and frontier semantics still constrain quality inconsistently

`search_control.py:145` stops at first qualification. This is suitable for a fast
feasible mode, not a claim that the best small pool has been sought across the
inventory. `panel_acceptance.py:219` calls repair outside the service budget;
`unified_optimizer.py:232` calls the optimizer directly for alternatives and catches
all exceptions as a stopped alternative search. Stage-one objective use remains
conditional in `swap_refinement.py:83`.

Provide explicit fast-feasible and quality-until-budget modes. Enforce one ledger
at evaluation boundaries, including alternatives and late repairs, and distinguish
algorithm approximations from full objective-based acceptance. A missed threshold
under a bounded search is not evidence that no qualifying pool exists.

### P2: feedback building blocks are ahead of feedback validation

`sequencing_feedback.py` currently contains a disjoint-experiment-ID guard; it has
no production caller. This does not negate the working BAM-gap and reach-fitting
paths, but it is not a complete provenance-aware, held-out learning workflow.

Require alignment/reference identity, read-filter/depth policy, pool composition,
reaction metadata and a frozen calibration artifact before a learned model affects
a new design. Fit and assess on distinct experiments; within-run blocked validation
alone does not establish transfer to a different reaction or sample.

### P2: the strongest scientific claims still exceed the validation available

The occupancy coverage model combines primer contributions under an explicit
independence approximation and a chosen reach. Correct arithmetic does not establish
calibrated sequencing-recovery probabilities. Published-data tests are valuable
retrospective checks, but ranking a few already-known pools is not prospective
proof that new NeoSWGA designs outperform competitors. Several fixtures use
published statistics rather than recalculating every metric from raw input.

Evaluate independently held-out outcomes, separating target read fraction, breadth
at fixed sequencing effort, uniformity and pool size. Prevent fixture-informed
weight tuning from being counted as independent validation.

## Engineering assessment

| Dimension | Assessment | Consequence |
|---|---|---|
| Correctness | Strong local regression coverage; incomplete cross-component contracts | Require production-path tests for every protection, not merely helper tests |
| Performance | Promising caches, intervals and staged search; no fair cross-tool measurements here | Benchmark complete pipelines and search-only runs separately |
| Maintainability | Many useful modules, but duplicated settings and acceptance paths remain | Consolidate ownership and retire redundant routes before adding algorithms |
| Data integrity/security | Stronger reference checks and failure/export controls; security not exhaustively audited | Make artifact identities mandatory and keep dependency/security checks actionable |
| Scientific validity | Useful documented proxies and retrospective checks; empirical generalization unestablished | Restrict claims to evidence tier and perform blinded/held-out comparisons |

CI currently treats lint, static typing and dependency audit as informational.
That is a concrete enforcement limitation, not evidence of a discovered security
vulnerability. Prioritize blocking checks for newly changed critical modules.

## What “best” should mean

For a declared target/background and experiment domain, seek the smallest pool
that achieves the required observed breadth and specificity at an acceptable
design and sequencing cost. Measure success on a Pareto frontier: a tool need
not win every metric to be the better choice for a stated use case.

NeoSWGA should aim to lead in reproducibility, trustworthy failure behavior,
small-pool tradeoffs and sequencing-informed iteration. Superiority in recovery
or speed remains a hypothesis until a pinned, fair comparison supports it.
The [implementation and evaluation roadmap](../design/2026-09-22-swga-leadership-plan.md)
sets the work and acceptance gates needed to test that hypothesis.

---

# Addendum, 23 September 2026: COATswga re-examined, and one defect it exposes here

The review above inspected COATswga statically and says so
(**Scope and evidence**: competitor programs "were not installed or benchmarked").
This addendum re-examines the same revision, `e2e0ff4`, and reports the first
execution of COATswga code recorded in this repository. It is still narrow: two
calls into its dimer model, no pipeline run, no amplification experiment.

Three of the observations below were already recorded elsewhere here and are
re-verified rather than new. They are marked as such, because a finding that
restates the repository's own notes is worth less than one that does not, and
conflating the two is how a review inflates its own contribution.

## What is new

### COATswga does not screen self-dimers, and says it does

Its README ("without forming primer-primer dimers or self-dimers"), the
docstring of its dimer routine, and the preprint ("checking for potential self-
or primer-primer dimers") all claim the screen. The shipped code screens pairs
only.

Measured by calling its own code: a perfectly self-complementary 12-mer scores
**-24.29** against itself, against its -2.79 threshold, and its `is_dimer`
returns `False`. Searching the 40 COATswga references across this repository
found no prior mention of this, so it appears to be new here.

NeoSWGA screens self-dimers in three places, against a separate threshold
`max_self_dimer_bp` (default 4, distinct from `max_dimer_bp` default 3):
`filter.py:560-563` inside `filter_extra`, whose own docstring lists five rules
and omits this one; `optimization_service.py:404-406`, where a screen that
empties a non-empty pool raises `NoCandidatesError` naming the threshold; and
`clique_optimizer.py:117-118`.

### Its variable-length support is row pooling, not length-aware design

COATswga's stated distinction is being "the only currently available pipeline
that supports the generation of variable-length primer sets". The capability is
real and the selection is not informed by it.

After `make_df` the frame carries `primer, fg_count, bg_count, ratio`
(`filter.py:206-208`), later `cov_len` and `sort_val`. **Length is not a
column.** The search reads `primer`, `fg_count` and `bg_count`
(`sets.py:109-110`, `:173-174`) and calls `len()` at exactly two sites
(`sets.py:103`, `:189`), both to compute an interval endpoint. Nothing
recomputes Tm, weights by duplex stability, or treats a 10-mer differently from
a 14-mer. A 6-mer and a 22-mer compete on `cov_len / ratio**2` alone.

This compounds with a Tm floor that does nothing. The ceiling is applied per k
before pooling (`filter.py:96`); `min_tm` is in the defaults, exposed as a flag
and documented, and read nowhere. (That last point is NOT new: it is recorded at
`pool_selection_audit_2026-09-18.md:362`.) So a design declaring EquiPhi29 at
42 C admits short oligos at any stability and ranks them against long ones on a
criterion blind to the difference. This repository's own measurement says what
that admits: at equiphi29 42 C, median occupancy is 0.001 at k=7, 0.205 at k=10
and 0.917 at k=12 (`variable_oligo_length_2026-09-21.md:78-81`). COATswga's
default k range starts at 8 and its shipped `params.json` starts at 6.

NeoSWGA's contrast is `length_occupancy.occupancy_by_length`, reported from
`filter` and from `optimize` and enforced nowhere. It reports and does not gate,
which is the resolution Known Issue 17 reached for the same quantity -- but it
is the measurement COATswga cannot make.

### Ranking on the mean of the two strand coverages hides a starved strand

COATswga runs two greedy passes, forward and reverse (`sets.py:112`, `:171`).
Reaching `target_coverage` on both is a **stopping heuristic, not a
requirement**: each loop breaks when the novelty ladder bottoms out
(`sets.py:164-165`, `:218-219`) whether or not the target was met; the two are
not symmetric, the reverse loop carrying an extra `coverage_change >= 0.01`
guard that stops it a rung early; and the padding loop to `min_set_size`
(`sets.py:222-228`) appends primers with a dimer check and no coverage check at
all.

Ranking is then on the **mean** of the two (`sets.py:262`, `:273`), so a set at
0.99 forward and 0.30 reverse scores 0.645 and beats a balanced 0.60/0.60 at
0.60. That is the failure `min_per_target_coverage` exists to name in NeoSWGA,
applied across strands instead of across targets. Two numbers are printed
(`sets.py:233-234`) and nothing flags the imbalance.

NeoSWGA has the data to make the equivalent check and does not make it either:
`strand_metrics` collects all five strand figures per genome onto
`PrimerSetMetrics.strand_stats`, and Known Issue 18 records that the two
headline strand scalars are deliberately not constrainable, because a zero there
cannot be told from an unmeasured one. That reason does not extend to a
balance REPORT, which is a gap on this side worth closing.

## What COATswga does better, and the defect it exposes here

**Its coverage model is directional, and that is mechanistically correct.**
A forward site covers `[pos, min(chr_len, pos + frag)]` (`sets.py:93`); a
reverse-complement site covers `[max(0, pos + len - frag), pos + len]`
(`sets.py:103`). The derivation is short. If oligo P occurs literally at `i`,
P anneals to the minus strand there, the nascent strand is plus-sense, and
extension runs toward increasing coordinates. If `rc(P)` occurs at `j`, P
anneals to the plus strand and extension runs toward decreasing coordinates.
One direction per occurrence, never both.

**NeoSWGA credits both.** Verified in source: both coverage entry points default
to `strand="both"` (`coverage.py:28`, `:117`), so the position list holds
occurrences of the oligo and of its reverse complement together; and
`_mark_window` marks `occupied[pos-extension : pos+extension]`
(`coverage.py:242`), a symmetric window `2r` wide centred on the site. So every
occurrence earns coverage in the direction it cannot extend in.

Put at its sharpest: `polymerase_extension_reach` returns an EXTENSION distance
and the code spends it as a RADIUS.

**The magnitude is not measured, and this entry does not claim one.** The
argument is structural: it predicts inflation, not how much. Where forward and
reverse-complement sites interleave, the wrong-direction half is often covered
anyway from the neighbouring opposite-orientation site -- which is the same
effect behind the record-geometry measurement reading exactly zero on
two-chromosome Prevotella while reading +4.788% on *Drosophila*
(`record_geometry_on_drosophila_2026-09-21.md`). The expectation is therefore a
small error on dense panels and a larger one on sparse panels, and sparse is
where SWGA designs live. Nothing here establishes that.

Anyone acting on this should note what it touches: `fg_coverage` is the
authoritative figure in every `step4_improved_df_summary.json` and in the
validation documents throughout this directory. Measure before changing.

**COATswga's own model collapses in one place.** In the reverse loop,
`sets.py:189` builds the candidate's reverse intervals from the FORWARD
occurrence list extended leftward, where `sets.py:103` correctly uses the
reverse-complement list for the same job; `sets.py:206` then builds forward
intervals from that identical list extended rightward. One occurrence earns both
sides there. Every other interval construction in the tool is one-sided
(`sets.py:93`, `:103`, `:134`, `:153`, `:206`, `filter.py:296`).

**Its record handling avoids this repository's record-join class by
construction, for a representation reason rather than a modelling insight.**
Intervals are `(chr, start, end)` triples keyed by chromosome
(`filter.py:256-260`), merged within chromosome, and totalled per record
(`sets.py:60-63`); both endpoints clamp to the record. A window spanning a join
is not expressible. NeoSWGA concatenates records into one coordinate space and
must then re-impose the boundaries, which `compute_per_prefix_coverage` does and
`merged_window_intervals` does not -- the disagreement recorded under **The two
production coverage paths** in CLAUDE.md. The same representation also removes
the scanner-fabrication class: `filter.py:142-150` scans each record separately,
so a k-mer spanning a join is never formed.

## What this addendum does NOT establish

- No COATswga pipeline was run. Its reported experimental results are neither
  reproduced nor disputed here, and the scope disclaimer above still governs.
- The directional over-credit in NeoSWGA is argued from source, not measured.
- Two calls into COATswga's dimer model are the whole of the execution evidence.
- The chemistry comparison rests on reading both codebases, not on any reaction.
  NeoSWGA has a reaction model where COATswga has none, and that model is
  largely extrapolated: of 23 constants in `registry/model_evidence.json`, 7 are
  `measured`, both coverage reaches are `assumed`, `mismatch_penalty` is
  `assumed` and uniform, and the betaine, trehalose and urea coefficients each
  carry an in-code retraction in `mechanistic_params.py`.
  `additive_specificity.md:24-30` states the same caveat for its own contents:
  "none of this is validated against a measured wet-lab outcome." Having a model
  and having evidence for it are different claims, and only the first is
  established.
