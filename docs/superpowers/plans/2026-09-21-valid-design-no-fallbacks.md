# Valid Oligo-Pool Design Without Silent Fallbacks Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task, inline. Steps use checkbox (`- [ ]`) syntax for tracking. Do not delegate without a separate user request.

**Goal:** Produce reproducible, constraint-satisfying oligo designs with explicit model limitations and failures, and evaluate how well they predict sequencing recovery.

**Architecture:** Replace compatibility-driven orchestration with one immutable design request and one authoritative evaluator/acceptance contract. Treat optimization algorithms as proposal generators, retain the complete QC-passing inventory, and use staged search with shared accounting. Introduce evidence-qualified chemistry and coverage models before enabling calibrated sequencing predictions.

**Tech Stack:** Existing Python 3.11+, dataclasses, NumPy, pandas, HDF5, Jellyfish, pytest and Hypothesis; existing BAM/CRAM readers. Add no new runtime dependency until a measured need or model requirement justifies it.

**Spec:** [Valid design contract](../../design/2026-09-21-valid-design-contract.md). Read this specification and the current code before implementing.

## Global Constraints

- Backward compatibility is not a requirement.
- Missing data, stale provenance, unsupported requested chemistry, non-finite calculations and software failures must not be replaced with another model, zero binding, an empty index or a successful result.
- Every QC-passing candidate remains available to search.
- No stage changes chemistry or thresholds.
- Software correctness does not establish experimental validity.
- Documentation uses modest scientific language. Do not use Unicode in Nextflow files.
- Keep existing module/function size ceilings; extract code instead of raising them.
- Preserve the current uncommitted work while revising it. Do not reset the workspace or commit unrelated edits.

## Starting point and known unfinished work

The shared `OptimizationRequest`/`run_panel_search`, objective, repair, swap and
contraction code already exists. Extend or replace these units deliberately;
do not add a second parallel service. Beam repair already exists in `panel_beam.py`.
The historical `score` stage now prepares ordered candidates; do not reintroduce
the retired amplification predictor under its old name.

The latest refill/budget changes are **unfinished**, not a validated baseline.
On 21 September the focused service, contract, refill, end-to-end, ensemble and
size checks produced **191 passed and 2 failed**. Failures are
`cli/iterate.py::run_expand_primers` at 247 lines versus 245, and
`core/parameter.py` at 1,605 lines versus 1,601. The earlier 5,599-test pass predates
these changes. Do not describe it as validation of the current worktree.

Current concrete replacement targets:

| Location | Behavior to remove or restructure |
|---|---|
| `candidate_source.open_source_or_list` | Inventory-open `ValueError` changes the search to the CSV/list; a condition mismatch can hide candidates |
| `coverage._record_starts_for` | Missing getter or getter exception removes record boundaries |
| `coverage.polymerase_extension_reach`, `product_reach` | Broad exceptions return a default reach |
| `thermodynamic_filter._check_heterodimer_pair` | Calculation error returns zero free-energy contribution |
| `thermodynamics.calculate_tm_batch` | Unexpected calculation failure becomes NaN without a required failed-run decision |
| `occupancy.discrimination_profile` | Failed primer calculations can be skipped in an aggregate diagnostic |
| `unified_optimizer` | Broad catches around proposal methods, validation and alternatives need explicit result/error handling |
| `search_control`, `optimization_service` | Budget currently counts only uncached shared-objective calls; reporting/proposal work is not a hard total-compute bound |
| `panel_acceptance.apply_configured_limits`, `collect_alternative_sets` | Late repair/alternative search can bypass the newly introduced budget |
| `optimization_service.run_panel_search` | Frontier summary needs a complete acceptance record; an empty filtered frontier must not prevent examination of later eligible candidates |

These are code-inspection findings, not claims that every path is reachable in
every command. Task 1 records reachability and tests each production path before
removing it. Allocation choices that preserve the exact computation, such as
dense versus lazy dimer checks, are implementation strategies rather than model
fallbacks; retain them only with equivalence tests.

## Delivery order and file boundaries

Tasks 1-4 establish trustworthy inputs and model support. Tasks 5-8 establish
consistent evaluation and search. Tasks 9-10 establish sequencing feedback and
release evidence. Each task has its own test/review gate; the experimental gate
does not block shipping honestly labelled proxy analysis.

| Unit | Responsibility |
|---|---|
| Existing `design_context.py`, new `design_request.py` | Resolve and validate immutable configuration, references, models and concentration policy |
| Existing `exceptions.py`, new `design_result.py` | Lightweight errors, run/qualification states and output eligibility |
| New `model_evidence.py`, `registry/model_evidence.json` | Versioned parameter provenance and supported domains |
| Existing `position_cache.py`, `candidate_inventory.py`, `candidate_source.py` | Verified reference answers and complete eligible candidate access |
| Existing `pool_objective.py`, new `panel_evaluation.py` | Authoritative measurements and hard-constraint assessment |
| Existing `optimization_service.py`, `search_control.py` | Search stages, incumbents, candidate traversal and budget ledger |
| Existing `pool_planner.py`, `panel_contraction.py`, `primer_expansion.py` | Request-specific size and feedback policies using the shared service |
| New `sequencing_feedback.py`, existing `bam_coverage.py`, `reach_calibration.py` | Experiment validation, measured coverage and versioned calibration |
| Existing reports and benchmark script | Render saved assessments; evaluate correctness, search quality and evidence limits |

Avoid expanding `unified_optimizer.py` or `parameter.py`. Remove their duplicate
resolution paths as consumers move to `design_request.py`.

## Task 1: Explicit failures and honest completion states

**Files:** modify `neoswga/core/exceptions.py`, `optimization_service.py`,
`unified_optimizer.py`, `thermodynamic_filter.py`, `thermodynamics.py`,
`neoswga/cli_unified.py`; create `neoswga/core/design_result.py` and
`tests/test_design_failure_contract.py`.

**Interface:** add dependency-light `DesignError` subclasses
`InvalidDesignRequest`, `ReferenceDataError`, `UnsupportedModelError`,
`ModelEvaluationError`. Preserve exception chaining. `SearchBudgetExhausted`
remains control flow distinct from these errors. Define independent run status,
termination reason, and panel qualification; do not overload `SUCCESS`.

- [ ] Add this output-eligibility test and tests that inject a thermodynamic error into a production design call and assert no recommended FASTA is written:

```python
import pytest
from neoswga.core.design_result import recommendation_allowed

@pytest.mark.parametrize("state,qualified,expected", [
    ("finished", True, True), ("finished", False, False),
    ("failed", True, False), ("interrupted", True, False),
])
def test_recommendation_requires_finished_qualified_run(state, qualified, expected):
    assert recommendation_allowed(state, qualified) is expected
```

- [ ] Run `python -m pytest tests/test_design_failure_contract.py -q`; confirm the new interface/behavior fails before implementation.
- [ ] Implement `recommendation_allowed(run_state: str, qualified: bool) -> bool` as `run_state == "finished" and qualified`. At the CLI boundary catch `DesignError`, write a structured failure with stage/input/model identifiers, and return a nonzero exit code. Unexpected exceptions also fail and retain a traceback. A numerical error is never classified as candidate QC rejection.
- [ ] Classify each broad catch on the design path: required calculation fails the run; an explicitly selected method failing fails the requested search; an optional renderer can fail independently after validated JSON is committed. Do not automatically try another algorithm or substitute zero/NaN. For bad user sequences, return named QC rejections rather than fabricated metrics.
- [ ] Extract the new source-setup code from `run_expand_primers` and search-setting resolution from `parameter.py`; restore both size checks without raising ceilings.
- [ ] Run the new tests, `tests/test_reaction_conditions_init_errors.py`, `tests/test_function_length_ratchet.py`, and `tests/test_module_size_ratchet.py`. Review the error-path diff and record this task as complete only when all pass.

## Task 2: One immutable request with complete provenance

**Files:** create `neoswga/core/design_request.py`,
`tests/test_resolved_design_request.py`; modify `design_context.py`,
`parameter.py`, `base_optimizer.py`, `schema/params.schema.json`,
`neoswga/cli/pipeline.py`, `neoswga/cli/plan_pool.py`, `neoswga/cli/iterate.py`.

**Interface:** `resolve_design_request(params: Mapping[str, object]) -> DesignRequest`.
The frozen request contains resolved conditions; reference/record manifests;
candidate-source identity; fixed and excluded oligos; objective and constraint
definitions; size policy; search budgets; seed; model and calibration identifiers;
and concentration policy. Nested content must also be immutable. Store a
canonical serialized request hash with results. Keep runtime caches out of it.

- [ ] Add tests for unknown keys, non-finite values, `coverage_reach=0`, negative budgets, contradictory fixed/excluded oligos, incompatible model requests and missing required specificity backgrounds. Explicit zero must not become a default through `or`.
- [ ] Include this production resolver regression:

```python
import pytest
from neoswga.core.design_request import resolve_design_request
from neoswga.core.exceptions import InvalidDesignRequest

def test_explicit_zero_reach_is_not_replaced_by_a_default():
    with pytest.raises(InvalidDesignRequest, match="coverage_reach"):
        resolve_design_request({"coverage_reach": 0})
```

- [ ] Run `python -m pytest tests/test_resolved_design_request.py -q` and observe the failures.
- [ ] Resolve defaults once, record their source, reject unknown/retired settings, and pass the result to every command. Remove read-time dependence on mutable `parameter` globals from evaluator code. Missing required fields receive a field-specific error, not an implicit model choice.
- [ ] Replace `OptimizationRequest.optimizer: Any` as the owner of scientific settings: the service receives the resolved request and constructs runtime proposal generators/evaluator separately. Require identical request hashes for equivalent CLI and Python requests.
- [ ] Run resolver tests plus `tests/test_no_schema_key_is_inert.py`, `tests/test_shared_optimization_contract.py` and command tests. Check that all supplied parameters either affect the resolved request or are rejected.

## Task 3: Verified references and complete candidate inventories

**Files:** modify `kmer_counter.py`, `position_cache.py`, `candidate_inventory.py`,
`candidate_source.py`, `filter.py`, `pipeline.py`, `coverage.py`; create
`tests/test_strict_reference_answers.py` and
`tests/test_candidate_inventory_completeness.py` under `tests/`.

**Interface:** `PositionCache.require_record_metadata` and `require_entries`
become mandatory for design. Replace `open_source_or_list` with explicit
`open_inventory_source(...)` and `open_explicit_source(...)`; neither changes
source type after an error. Explicit candidate files undergo the same applicable
QC and verified position acquisition as inventory candidates.

- [ ] Add tests distinguishing a recorded empty site array from a missing primer/prefix answer, a corrupt index, an old schema, and a mismatched reference digest. Assert that all but the recorded empty array fail. Inject a metadata getter error into coverage and assert that it propagates as `ReferenceDataError`.
- [ ] Add a candidate-survival test that constructs more eligible candidates than the initial frontier, prepares candidates, then exhausts the source and checks:

```python
def assert_inventory_preserved(qc_eligible, visited, rejection_reasons):
    assert set(visited) == set(qc_eligible)
    assert not set(qc_eligible).intersection(rejection_reasons)
```

  Call this assertion from fixtures covering multiple lengths and changed chemistry. The fixture's expected sequences must be independently enumerated from the small input FASTA, not read back from the inventory being tested.
- [ ] Run the two new test modules and confirm failure for the existing permissive paths.
- [ ] Require full reference digest, record IDs/lengths/circularity, strand convention, k range, counting tool/version and schema metadata. Rebuild incompatible artifacts only through an explicit preparation step; never accept them as empty/current. Verify count/position agreement on exhaustive small references, including reverse complements, palindromes, ambiguous bases and circular origins.
- [ ] Persist candidate sequence identity independently of chemistry-specific eligibility. Preserve every hard-QC survivor; retain reason codes for rejections and model/condition fingerprints for derived eligibility. Ranking, top-N views and memory frontiers cannot change the inventory universe. Reassess stale eligibility or fail with a named preparation command; never use stale CSVs instead.
- [ ] Rename the candidate-preparation command from `score` to `prepare-candidates`, remove the alias and retired score parameters, and update pipeline/docs/tests together. Preserve existing measured features without inventing an amplification probability.
- [ ] Run the new tests and existing inventory, geometry, stale-provenance and score-retirement tests, updated for the intentional API break. Review every removal from the eligible universe against a recorded hard-QC reason.

## Task 4: Evidence-qualified chemistry and concentration models

**Files:** create `neoswga/core/model_evidence.py`,
`neoswga/core/registry/model_evidence.json`, `tests/test_model_evidence_contract.py`;
modify `reaction_conditions.py`, `additives.py`, `thermodynamics.py`,
`occupancy.py`, `thermodynamic_filter.py`, `registry/polymerases.py`,
`registry/views.py`, `design_request.py`, `pyproject.toml`; create
`docs/validation/chemistry_model_evidence.md`.

**Interface:** a versioned evidence record identifies the quantity, units,
model equation/version, source location, measured versus assumed status,
sequence-length domain, buffer/concentration/temperature domain, modification
support, mixture support and uncertainty. `require_model_support(request) -> None`
raises `UnsupportedModelError` for a requested computation outside its supported
domain. An assumption is not promoted to a measurement because it has a citation.

- [ ] Add tests for unknown polymerase, unsupported modified base, unsupported additive mixture, missing required nearest-neighbor parameters, concentration-unit mistakes and unmodelled requested dimer chemistry. Each must fail before search. Test that changing conditions invalidates all derived caches.
- [ ] Add and run a concentration conservation test for the new request method:

```python
import pytest

def assert_concentration_conserved(request, oligos):
    values = request.concentrations_molar(tuple(oligos))
    assert len(values) == len(oligos)
    assert sum(values) == pytest.approx(request.total_primer_molar)
    assert all(value > 0 for value in values)
```

  Use this for a fixed-total request with two and four oligos; separately assert fixed-per-oligo mode keeps each concentration constant while the total changes.
- [ ] Audit original papers and authoritative parameter tables, recording exact supporting tables/equations and experimental domain. Prioritize salt/free-magnesium treatment, additive coefficients, mismatch terms, primer concentration and the interpretation of polymerase reach. Verify whether evidence actually concerns short oligos and the requested reaction rather than assuming transfer from PCR or bulk product sizes.
- [ ] Include the registry JSON in package data and add an installed-package smoke test that loads it without relying on the repository working directory. Fail model initialization if the evidence artifact is missing or its schema/digest is invalid.
- [ ] Remove unsupported default mismatch penalties from strict predictive calculations. Keep exact-match analysis an explicit model with exact-match-only specificity claims. Enable context/position-dependent mismatch or modified-oligo models only where parameters and domain tests exist. Terminal protection, base substitutions, mismatched binding, extension competence and primer-primer dimers require distinct capability records.
- [ ] Implement `DesignRequest.concentrations_molar(oligos: tuple[str, ...]) -> tuple[float, ...]` for declared allocation modes and propagate it into every evaluation/cache key. Recheck concentration-sensitive QC during panel evaluation; use only safe, documented bounds for candidate-level pruning.
- [ ] Record geometry reach as an explicit assumption or calibration, not polymerase identity converted silently to predicted recovery. Consolidate conflicting reach documentation and remove runtime exception-to-default paths.
- [ ] Run chemistry/evidence/concentration tests and existing thermodynamic reference tests. Release only supported computations. Longer oligos plus additives remain a hypothesis to compare under supported models and validate experimentally; lowering Tm alone does not establish improved selectivity or recovery.

## Task 5: One authoritative panel evaluator and independent coverage oracle

**Files:** create `neoswga/core/panel_evaluation.py`,
`tests/test_authoritative_panel_evaluation.py`,
`tests/test_coverage_independent_oracle.py`; modify `pool_objective.py`,
`base_optimizer.py`, `occupancy_coverage.py`, `coverage.py`, `panel_acceptance.py`,
`lazy_dimer.py`, `pool_plan_report.py`.

**Interface:** `evaluate_panel(request, oligos) -> PanelAssessment`, an immutable
assessment holding panel identity, request hash, named metrics and units,
per-target results, evidence/model identity, all hard-constraint violations and
qualification. Internally separate measurements from policy, but return one
complete acceptance record. Numerical failures raise; a missing optional quantity
is explicitly unavailable, and a constraint requiring it fails preflight.

- [ ] Write a separate base-by-base oracle for tiny multi-record references. Do not call production interval-merging code from the oracle. Test linear ends, per-record circular wrap, no sites, overlapping sites, strand convention and every record boundary.
- [ ] Pin the currently implemented occupancy grouping without calling the production evaluator to compute the expected answer:

```python
import pytest
from neoswga.core.occupancy_coverage import occupancy_weighted_coverage

def test_no_sites_is_zero_even_with_a_long_reach():
    assert occupancy_weighted_coverage(
        {}, 100, extension_reach=1000, circular=True, conditions=object()
    ) == 0.0

def assert_two_independent_primers_cover_half_a_window(actual):
    # Two distinct primers, each theta=0.5, both reaching 40 of 100 bases.
    assert actual == pytest.approx(0.4 * (1 - (1 - 0.5) ** 2))
```

  Supply a controlled condition/enthalpy fixture to the second assertion. Test repeated overlapping sites of one primer separately. Label independence as a model assumption rather than establishing its empirical correctness through this test.
- [ ] Run the new oracle/assessment tests and observe failures for missing geometry and disagreement between acceptance/reporting quantities.
- [ ] Route coverage, specificity and configured dimer checks into `PanelAssessment`; use the same record in all stages and exports. Reject NaN/non-finite required quantities. Represent a verified zero-background denominator explicitly; do not serialize infinity as ordinary JSON or interpret it as guaranteed enrichment.
- [ ] Keep exact and mismatch background measurements separate. Identical duplexes under identical local assumptions must get identical occupancy regardless of foreground/background label. A hypothetical mismatch discrimination diagnostic must not become measured specificity when the index contains only exact matches.
- [ ] Report geometric and occupancy-weighted coverage with their declared reach and denominator. Add uniformity, largest holes and per-target floors where requested. Any second reach is a named sensitivity scenario, not a replacement headline metric. Do not call the proxy a recovery probability.
- [ ] Run independent-oracle, occupancy, coverage-boundary, dimer and report tests; require saved JSON and rendered report to agree for the exact exported panel.

## Task 6: Shared staged search and complete budget accounting

**Files:** modify `optimization_service.py`, `search_control.py`,
`panel_refinement.py`, `panel_acceptance.py`, `unified_optimizer.py`,
`optimizer_factory.py`, `pool_planner.py`, `primer_expansion.py`; create
`tests/test_search_budget_contract.py`, `tests/test_frontier_service_contract.py`.

**Interface:** `run_panel_search` becomes the sole orchestrator over the Task 2
request and Task 5 evaluator. Proposal methods advertise capabilities and return
proposals, not authoritative success/acceptance. A shared ledger owns uncached
objective-evaluation limits and cooperative elapsed-time limits for the declared
run scope; separate counters record proposal effort and reporting evaluations.

- [ ] Add and run this exact-cap test, then production-path tests spanning proposal selection, repair, swaps, deletion, refills, ensemble methods and alternatives:

```python
import pytest
from neoswga.core.search_control import SearchBudget, SearchBudgetExhausted

def test_one_allowance_cannot_be_reset_by_the_next_stage():
    budget = SearchBudget(max_evaluations=2)
    budget.consume()
    budget.consume()
    with pytest.raises(SearchBudgetExhausted):
        budget.consume()
    assert budget.evaluations == 2
```

- [ ] Bind the ledger at the evaluator boundary. Direct `compute_metrics` calls made for search must not bypass it. Cache hits remain separately counted; non-search final assessment is explicitly recorded and cannot be used as an uncharged way to explore proposals. Use deterministic fake-clock tests for cooperative time limits.
- [ ] Implement stages: prepare verified frontier, generate proposals, assess, repair, refine, attempt deletions, revisit swaps, and widen the frontier under the selected search policy. Return valid incumbents at budget boundaries, including improvements already accepted within a stage. Never turn a numerical error into a budget stop.
- [ ] Reuse one ledger across requested sizes, ensemble methods, refills, late repair and alternatives. A chemistry grid has a declared command-wide allowance and recorded allocation by condition; no silent resets. Do not claim a hard wall-clock limit for opaque calls. If such a limit is required, introduce cancellable worker execution and verified checkpoints as a separately tested execution policy.
- [ ] Empty filtered frontiers advance when unexamined candidates remain. Prepare their positions and apply eligibility/exclusion checks before gap scoring. Refill stages preserve fixed oligos and assess the best prior incumbent. A refill's final stage record contains full qualification data so ensemble ranking cannot accidentally favor an invalid panel.
- [ ] Remove duplicated late repair and direct alternative optimization. Route standard, planning, grid, expansion, contraction and alternatives through the same contract. Unrequested solver substitutions are errors. Capability-limited algorithms may propose approximate panels, but cannot make unsupported chemistry or acceptance claims.
- [ ] Run the two new test modules, service/ensemble/refill tests and size ratchets. Demonstrate that the counted allowance is never exceeded and every termination has a reason and candidate-examination statistics.

## Task 7: Pool-size search that continues beyond first feasibility

**Files:** modify `pool_planner.py`, `optimization_service.py`,
`panel_refinement.py`, `panel_contraction.py`, `panel_beam.py`;
create `tests/test_smallest_pool_search.py`.

**Interface:** keep the best qualifying incumbent for each actual size and a
Pareto archive of size, coverage, background load and relevant uncertainty.
Support `first_feasible` and `improve_until_budget` policies; default design
quality searches use the latter with an explicitly resolved finite budget.
Define ranking centrally: hard feasibility first; smallest qualifying size;
then the request's declared coverage/specificity tie-breaks.

- [ ] Add a small independently enumerated fixture where a larger frontier supplies a better or smaller qualifying pool even though the initial frontier already meets the target. Add a fixture where greedy deletion gets stuck but swap-then-delete succeeds. Enumerate all subsets for at most ten candidates in the test, not through the production solver.
- [ ] Add the following result assertion to each exact fixture:

```python
def assert_matches_enumerated_optimum(result, qualifying_subsets):
    best_size = min(map(len, qualifying_subsets))
    assert len(result.primers) == best_size
    assert frozenset(result.primers) in {frozenset(p) for p in qualifying_subsets}
```

- [ ] Run `python -m pytest tests/test_smallest_pool_search.py -q` and confirm that early-stop or deletion-only behavior fails the appropriate fixture.
- [ ] Continue candidate examination after feasibility in quality mode using the same specificity-aware assessment. Maintain incumbents instead of letting a wider frontier replace a better panel. Combine additions, swaps, multi-step beam repair and deletion under the ledger; use restarts only with recorded seeds and budget allocation.
- [ ] Search configured sizes with reusable incumbents, but do not assume feasibility is monotone with size, especially under fixed-total concentration or dimer constraints. Do not binary-search size without a proven applicable monotonicity property. Evaluate final panel concentration again after every composition change.
- [ ] Calculate lower bounds only for a documented relaxation that bounds the actual objective. Exhaustive toy cases may produce minimum certificates. Large heuristic runs say `smallest qualifying pool found`, expose the searched size range, unseen candidates and stopping reason, and never equate `no pool found` with infeasibility.
- [ ] Run oracle/minimization/fixed-primer/dimer tests. Report all qualifying count/coverage alternatives so the user can compare synthesis cost and predicted performance rather than receiving an unexplained fixed pool size.

## Task 8: Auditable outputs and reproducible comparisons

**Files:** modify `pool_plan_report.py`, `neoswga/core/report/metrics.py`,
`scripts/benchmarking/sequential_panel_search.py`; create
`tests/test_design_report_provenance.py`,
`docs/validation/strict_design_benchmark.md`; update
`docs/guides/optimization_guide.md` and `CLAUDE.md` to the new contract.

**Interface:** a versioned result schema stores the immutable request, reference
and model digests, inventory/QC/examination counts, actual panel, assessment,
stage changes, budget ledger, run status, termination and claim tier. Rendering
consumes saved results without reinterpreting defaults or recalculating metrics.

- [ ] Write round-trip tests proving JSON, HTML, CSV and recommended FASTA identify the same panel and qualification. Test failed final validation and stale output directories so neither can leave an apparently current recommendation.
- [ ] Run `python -m pytest tests/test_design_report_provenance.py -q` before implementation.
- [ ] Publish validated machine-readable results atomically and generate recommendation exports only through Task 1's gate. Display measured versus predicted coverage, concentration policy, unsupported/unavailable quantities, all constraint results, budget stopping and evidence limitations prominently.
- [ ] Replace the exploratory per-stage benchmark comparison with equal declared run allowances, the same references/inventory/constraints and multiple fixed seeds. Include proposal generation in timings; state precisely what the evaluation count covers. Give each arm a fresh evaluator cache or report shared warm-cache effects explicitly. Keep final verification separate and identical.
- [ ] Use wMel plus the supplied full host reference, synthetic adversarial geometry, multiple lengths supported by the selected model, several coverage targets and fixed-total/per-oligo policies. Compare smallest qualifying size, coverage, specificity, failures, runtime, memory and candidate examination. Show uncertainty across seeds; do not claim superiority from one favorable panel.
- [ ] Run report tests and benchmark smoke fixtures. Document previous 9-oligo/70.07% results as historical proxy results, not evidence that this plan improves experimental recovery.

## Task 9: Sequencing feedback with verified coordinates and held-out evaluation

**Files:** create `neoswga/core/sequencing_feedback.py`,
`tests/test_sequencing_feedback_contract.py`,
`docs/design/sequencing_feedback_contract.md`; modify `bam_coverage.py`,
`reach_calibration.py`, `deficit_objective.py`, `primer_expansion.py`,
`neoswga/cli/iterate.py`.

**Interface:** a feedback artifact records experiment ID, design/request hash,
oligo identities/concentrations, reaction conditions, reference digests, alignment
and filtering provenance, per-record depth summaries, target/background read
fractions, excluded regions and calibration eligibility. A calibration artifact
adds training/validation experiment IDs, model version/domain, fitted parameters,
uncertainty and out-of-sample results.

- [ ] Write fixtures for reference-name collision, same-length but different sequence, missing contigs, coordinate aliases and mismatched design/chemistry. Reject ambiguous matching; do not auto-assign a BAM record by length. Explicit aliases require verified coordinate identity, not just similar names.
- [ ] Add and run this leakage guard for the new artifact validator:

```python
import pytest
from neoswga.core.sequencing_feedback import require_disjoint_experiments

def test_an_experiment_cannot_validate_its_own_fit():
    with pytest.raises(ValueError, match="overlap"):
        require_disjoint_experiments(("run_a", "run_b"), ("run_b",))
```

- [ ] Implement `require_disjoint_experiments(training: tuple[str, ...], validation: tuple[str, ...]) -> None`. Specify depth thresholds, mapping/base quality, duplicate and multimapper policy, callable-reference denominator, and multi-contig handling in the feedback contract. Report observed breadth at several declared depths and sequencing effort; do not compare runs at unequal depth without an explicit adjustment.
- [ ] Reuse the existing blocked reach cross-validation and informative-fit flags. Add guarded spatial holdouts/buffers when reach-sized windows could overlap folds; distinguish within-run tuning diagnostics from independent-run validation. Refuse to apply an uninformative or out-of-domain calibration; retain the original design plus the reason as a diagnostic result.
- [ ] Use low-depth regions to propose targeted additions and swaps under the original hard constraints. Separate likely mapping/repeat problems from design deficits using recorded masks. Avoid assigning per-oligo causal efficacy from one mixed-pool sequencing run; fit such effects only with enough varied, replicated pools to identify them.
- [ ] For each redesigned pool save baseline, changes, predicted deficit improvement, size/cost and constraint assessment. Evaluate the fixed design on subsequent held-out sequencing runs. Do not silently retrain a model during optimization or use validation results to select a panel while still calling them held out.
- [ ] Run feedback, BAM, reach-calibration, deficit-selection and expansion tests. Verify that an intentionally incompatible sequencing run cannot influence a new design.

## Task 10: Separate software release and experimental validation gates

**Files:** create `docs/validation/design_release_gates.md`,
`tests/integration/test_strict_design_pipeline.py`; update benchmark and chemistry
evidence documents and the wMel example instructions.

**Interface:** every released model/claim tier has a validation record. A software
release may support proxy design without claiming validated sequencing recovery.
Only models with applicable empirical evidence may expose calibrated-recovery
claims, and the report cites their frozen validation record.

- [ ] Add an end-to-end small-reference test covering count, filter, preparation, design, contraction and sequencing-informed expansion, plus error injections at each required step. The failure cases must return nonzero status and no qualifying recommendation artifact.
- [ ] Run the new integration test and all focused tests from Tasks 1-9; fix failures before broader testing.
- [ ] Run `python -m pytest tests/ -q`, repository formatting/lint checks and both size ratchets. Record skips and reasons. Require deterministic fixed-seed behavior where promised and preserve provenance for nondeterministic solvers.
- [ ] Predeclare empirical evaluation endpoints: observed callable breadth at specified depths, target read fraction, breadth uniformity and oligo count/cost. Use independent amplification replicates and held-out experiments. Include the previous pool and an appropriate baseline; compare chemistry/length changes without attributing a combined change to one factor. Keep this document a study design, not an unsupported reaction recipe.
- [ ] Define numerical empirical acceptance margins before observing the validation outcomes, based on the intended application's coverage requirements and a pilot-based precision/power assessment. Version and freeze those margins with the study; software cannot invent a universal acceptable recovery threshold.
- [ ] Release supported proxy calculations after the software gate. Promote a calibrated model only after its predeclared empirical gate passes. Report negative findings and revise or retire models that do not transfer; never substitute another model silently to make the validation pass.

## Completion and review checklist

- [ ] No required design calculation catches an error and returns a plausible replacement value.
- [ ] Every accepted setting reaches one resolved request; unsupported settings fail.
- [ ] All required reference answers are verified, and every QC survivor remains searchable.
- [ ] Every delivered panel passes the same complete assessment under its actual concentration and chemistry.
- [ ] Refilling and size reduction share limits, retain incumbents and expose unexamined candidates.
- [ ] Final reports distinguish proxy coverage, calibrated prediction and observed sequencing breadth.
- [ ] Sequencing calibration is versioned, domain-limited and evaluated without data leakage.
- [ ] All release-gate results refer to the actual final code and model versions.

## Recommended next implementation increment

Complete Tasks 1-3 first: explicit failures and result states, a single resolved
request, and strict reference/candidate provenance. These address conditions
that can make a plausible result invalid. Then establish the supported chemistry
and concentration contract before increasing search complexity. Search quality
benchmarks and sequencing calibration follow once the evaluator has a stable,
testable scientific meaning.

This plan replaces the unfinished compatibility-oriented refill/budget work as
the roadmap. Reuse its useful implementation pieces, but do not finish old
fallback behavior merely because it was already started. Implementation remains
inline unless the user requests delegation. Review each task's scoped diff and
record test evidence before marking it complete; commit only the intended files
when committing is part of the active execution request.
