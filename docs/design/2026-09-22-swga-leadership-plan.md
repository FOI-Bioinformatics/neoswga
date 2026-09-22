# NeoSWGA comparative improvement and validation plan

Date: 22 September 2026. Input: the
[pinned comparative review](../validation/2026-09-22-swga-comparative-review.md).
This is a prioritized follow-up to the valid-design plan, not a claim that all
earlier planned capabilities are now connected. Backward compatibility is not a
requirement. Execute inline unless delegation is requested.

## Objective and decision rule

Make NeoSWGA a well-supported choice for designing small, specific SWGA pools by
demonstrating correctness, useful search tradeoffs and experimental performance.
Do not define success as having the most methods, the highest internal score or
the largest number of passing tests.

Compare observed callable breadth at fixed sequencing effort, target read fraction,
uniformity, oligo count/cost, runtime, peak memory and reproducibility. Predeclare
the use-case priorities and acceptable tradeoffs; publish failures and cases where
another tool is better. Restrict superiority claims to tested domains.

## 1. Connect the contracts before adding features — release blocking

**Files:** `core/design_request.py`, `design_context.py`, `panel_evaluation.py`,
`panel_acceptance.py`, `optimization_service.py`, `unified_optimizer.py`,
`core/schema/params.schema.json`, `cli/pipeline.py`, `cli/plan_pool.py`,
`cli/iterate.py`, `core/pool_plan_report.py` (all under `neoswga/`).

- [ ] Merge parameter-file values, defaults and explicit CLI overrides once.
  Pass the resulting deeply immutable request to every design path. Include
  geometry, objectives, seeds, all active search controls and artifact/model
  identities in its recorded form. Remove duplicated global-parameter reads.
- [ ] Reject zero/negative sizes and unsupported concentration policies rather
  than substituting defaults or accepting settings with no execution path.
- [ ] Complete `PanelAssessment`: actual membership/count, fixed/excluded oligos,
  QC eligibility, hard dimer restrictions, raw/effective coverage, per-target
  limits, background policy, chemistry identity and unavailability reasons.
  Make this the required search acceptance and recommendation export gate.
- [ ] Connect concentration policy to every actual panel evaluation and cache
  key. At fixed total concentration, reevaluate all remaining oligos after an
  addition, swap or deletion. Keep concentration-dependent QC separate from
  safe candidate-level rejection. Expose unsupported modes only after this works.
- [ ] Add real-command tests proving CLI overrides reach both evaluator and
  manifest; fixed-total size changes affect evaluation; required assessment
  failure blocks export; and no alternative command bypasses the contract.
- [ ] Turn the review probes into regression tests for corrected behavior.
  Existing helper tests remain useful but cannot satisfy this gate alone.

**Gate:** all entry points produce equivalent effective requests and final
assessments for equivalent inputs; injected model/reference failures export no
recommendation; changing any active setting changes the recorded request when
appropriate. Preserve current module/function ceilings through extraction.

## 2. Establish a defensible scientific model — release blocking for model claims

**Files:** `core/model_evidence.py`, `core/registry/model_evidence.json`,
`core/reaction_conditions.py`, `core/thermodynamics.py`, `core/occupancy.py`,
`core/occupancy_coverage.py`, `core/coverage.py`, `core/lazy_dimer.py`;
create `docs/validation/model_domain_review.md`.

- [ ] Verify parameter sources against primary tables/equations, with units and
  source-specific conditions. Record what was measured, extrapolated or assumed.
  Review salt, magnesium, concentration and additive interactions separately.
- [ ] Replace free-text-only model domains with checked ranges/capabilities for
  length, temperature, buffer, mixture and modification support. Distinguish
  recommended enzyme operating ranges from thermodynamic validity ranges.
- [ ] Support explicit geometry and chemistry-proxy analyses with appropriately
  limited claims. Unsupported required calculations fail; no catch-and-substitute
  model selection. An explicit approximation is a declared model, not a fallback.
- [ ] Audit directional versus symmetric windows with an independent small-genome
  oracle. Do not treat COAT fragment length and NeoSWGA window radius as equivalent
  settings. Report sensitivity to reach and occupancy assumptions.
- [ ] Keep exact-match specificity distinct from mismatch-aware analysis and
  observed enrichment. Require identical model treatment for identical target
  and background duplexes; assess mismatches only with a supported model/index.
- [ ] Test identical inputs across scalar/batch/cached implementations and all
  supported chemistry paths. Check numerical validity and interval boundaries.

**Gate:** every computed scientific quantity has a model/version/domain and
evidence label; independent oracle tests agree within declared numeric tolerance;
no report equates a geometric/occupancy proxy with observed sequencing recovery.
Unverified chemistry blocks the relevant claim, not unrelated geometry analysis.

## 3. Improve search quality with a small set of accountable methods

**Files:** `core/search_control.py`, `optimization_service.py`, `pool_planner.py`,
`panel_refinement.py`, `panel_contraction.py`, `panel_beam.py`,
`swap_refinement.py`, `optimizer_factory.py`, `candidate_source.py`.

- [ ] Keep one fast baseline and one quality search as the primary user choices.
  Retain other algorithms only as measured proposal generators or explicit
  research options; method count is not a quality metric.
- [ ] Provide `first_feasible` and `improve_until_budget` policies. The quality
  policy must continue examining candidates after initial feasibility when it
  can improve pool size, coverage, specificity or robustness under the request.
- [ ] Preserve all hard-QC survivors in the inventory. Show the counts admitted,
  ranked, examined and not yet examined. Never make a top-N shortlist permanent
  exclusion or relax QC/constraints to fill a requested pool.
- [ ] Integrate swaps, additions, deletion and existing bounded beam repair under
  one objective. Keep feasible incumbents separately from exploratory states.
  Retain accepted within-stage improvements when a budget expires.
- [ ] Account for all search evaluations, including alternatives, proposal scoring
  and late repairs. Track cached and uncached work separately. Share command-wide
  allowances across size/chemistry sweeps with an explicit allocation policy.
- [ ] Maintain size/coverage/specificity tradeoffs across multiple seeds. Do not
  assume larger panels improve performance under fixed-total concentration.
  Label results as the smallest qualifying pool found unless a valid exact
  certificate covers the full declared universe and objective.
- [ ] Add exhaustive small-instance oracles and difficult synthetic cases where
  local greedy deletion fails. Benchmark ablations: refill, swaps, beam, deletion,
  chemistry weighting and directional coverage, at equal declared resources.

**Gate:** no hard-constraint regressions; caps are respected under their stated
scope; known toy optima are recovered by the exact reference path; quality mode
retains a no-worse feasible incumbent. Choose default budgets from benchmark
quality/runtime curves rather than arbitrary expansion of existing limits.

## 4. Build a neutral four-tool benchmark

**Create:** `benchmarks/swga_comparison/manifest.json`, `adapters/`,
`evaluate_outputs.py`, `run_comparison.py`, and a documented result schema;
extend `scripts/benchmarking/sequential_panel_search.py` only if it can remain
focused. Pin the three external revisions recorded in the review.

- [ ] Supply isolated, reproducible environments. Original SWGA needs a suitable
  legacy runtime; dependency incompatibility must be reported rather than silently
  substituting NeoSWGA's reimplementation for the competitor.
- [ ] Start with Wolbachia/Drosophila, synthetic multi-record/circular references,
  and nonpathogenic references spanning GC content, related backgrounds, genome
  size and strain variation. Use published outcome datasets for retrospective
  external validation; do not tune and validate on the same datasets.
- [ ] Run two separate comparisons: native complete pipelines, and selection-only
  comparisons where a common eligible inventory/constraint set is representable.
  Report unsupported settings. Record the candidate universe for each run.
- [ ] For all runs fix input checksums, CPU/memory allowance, wall time, tool
  version and seed policy. Report cold and warm starts separately. COAT CPU count
  can alter starts searched, so show it as both a resource and search setting.
- [ ] Evaluate exported oligos using an independent evaluator with explicit
  common geometry and reference definitions. Include each tool's native score
  for interpretation but never rank tools solely by NeoSWGA's own objective.
- [ ] Compare multiple fixed seeds and include failed/timeout runs. Report interval
  uncertainty for aggregate outcomes and quality-versus-runtime curves, not only
  the best run. Separate counting/index cost from selection cost without hiding
  either in an end-to-end comparison.
- [ ] Publish oligo counts, constraint violations, estimated breadth, largest
  holes, per-target coverage, background load, wall time, peak memory and all
  deviations needed to run old software. Correct comparator defects in separate
  labelled runs; never silently patch the baseline.

**Gate:** another researcher can reproduce the benchmark from pinned manifests.
An engineering advantage claim requires the declared quality target at lower
resource cost, or better independently assessed quality within the same allowance,
without hiding a specificity/size tradeoff. No claimed win from incomparable
coverage definitions or unequal candidate pools.

## 5. Validate efficacy modelling and sequencing-informed redesign

**Files:** `core/sequencing_feedback.py`, `bam_coverage.py`, `reach_calibration.py`,
`deficit_objective.py`, `primer_expansion.py`; create a versioned experiment and
calibration schema and `tests/integration/test_feedback_redesign.py`.

- [ ] Preserve working depth-policy and reference-layout checks. Add experiment,
  pool, concentration, chemistry and mapping provenance to feedback artifacts.
  Uninformative or out-of-domain fits must not silently alter a new design.
- [ ] Connect the disjoint-experiment guard to fitting and application. Use guarded
  spatial folds for within-run diagnostics and separate experiments for transfer
  validation. Record callable masks, mapping filters and sequencing effort.
- [ ] Generate constrained gap-targeted additions/swaps and compare them with the
  original pool and equal-size alternatives. Do not infer each oligo's causal
  efficacy from one mixed-pool run.
- [ ] Reassess an optional efficacy model using experimental training data and
  appropriate held-out sequences/experiments. Include SOAP's poor-primer filtering
  concept as a baseline; do not substitute synthetic rule labels for empirical
  amplification measurements. Promote it only if it adds out-of-sample value over
  simpler models and does not silently remove useful gap-filling candidates.
- [ ] Freeze the model before prospective evaluation. Predeclare breadth at fixed
  total sequencing bases as a primary endpoint; target read fraction, uniformity,
  failure frequency and pool cost as secondary endpoints. Equal target-read
  subsampling is useful diagnostically but cannot establish enrichment efficiency.
- [ ] Evaluate Wolbachia designs with independent amplification replicates and
  held-out material, using published/competitor pools as baselines. Choose sample
  size and numerical noninferiority/superiority margins from application needs and
  pilot-based precision estimates before unblinding results. Avoid attributing a
  combined chemistry/length/pool change to one causal factor.

**Gate:** redesigned pools improve predeclared held-out outcomes or demonstrate
acceptable breadth/specificity with fewer oligos. Claims remain specific to the
validated sample, chemistry and organism domain; negative results remain visible.
Prospective measurements require laboratory data and cannot be completed by code
refactoring alone.

## 6. Make the application easier to trust and use

**Files:** CLI request builders, reports, `README.md`, guides and
`.github/workflows/ci.yml`; add installation and export smoke tests.

- [ ] Provide a guided path: references/backgrounds, supported chemistry, coverage
  and specificity goals, pool-cost limits, then fast or quality search. Offer a
  small set of interpretable outputs rather than requiring users to choose among
  many algorithms.
- [ ] Show a coverage-versus-oligo-count frontier, largest unresolved gaps,
  background tradeoffs and the exact reason a request could not be satisfied.
  Keep measured sequencing results visibly distinct from estimates.
- [ ] Make each experimental capability's state clear: implemented, connected,
  tested against an oracle, retrospectively evaluated, or prospectively validated.
  Generate model support tables from checked metadata.
- [ ] Test installed wheels and example workflows, not only editable checkouts.
  Make lint/type checks blocking for new critical modules and use a tracked,
  shrinking allowlist for inherited issues. Treat dependency-audit findings
  explicitly rather than always suppressing their exit status.
- [ ] Publish a limitations page and a regression dashboard from the independent
  benchmark. Retire options with no measured benefit or reliable support.

**Gate:** a fresh installation can run a documented small design and interpret
its output without implementation knowledge. Every advertised capability has a
production-path test and an accurate evidence tier.

## Order and completion policy

Implement 1 first. Establish model support in 2 and benchmark infrastructure in 4
before selecting default search improvements in 3. Complete feedback integration
and empirical evaluation in 5; deliver usability and quality gates in 6 alongside
each functional increment. Keep all necessary scientific failures explicit.

Update the earlier plan's partial tasks with test evidence, not merely the
presence of new classes. A software release can legitimately offer well-tested
proxy design before prospective validation, but cannot claim experimental
superiority until the appropriate held-out comparison has passed.
