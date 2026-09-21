# Optimization entry points and chemistry audit

## Implementation follow-up

The findings below describe the audited baseline. Subsequent uncommitted changes
address chemistry-error fallbacks, parameter propagation, and a shared panel-stage
service (`OptimizationRequest` / `run_panel_search`). Standard selection, ensemble
proposal comparison, planning/grid rows and configured expansion now use the
shared objective for repair/refinement; reduction preserves limits and can revisit
swaps after a deletion. Expansion receives the design context, and its dominating
set adapter forwards fixed primers with a new-primer budget. Stage history is
saved with results and displayed in pool reports. `contract-set` now uses the
same constrained deletion stage, with read-only FASTA scanning of supplied
oligos and explicit raw/effective coverage and target attainment.

Bounded beam repair already exists in `panel_beam.py` and `repair_panel`; it is
not an outstanding implementation item. Incremental interval performance work
remains conditional on profiling.

Remaining architecture work: consolidate reference/candidate-source setup and
frontier refill across front ends, add a total-run budget spanning proposal and
all improvement stages. Final repair now precedes derived coverage and validation
reporting and records any changed panel in the stage history. A strict equal-total-budget
benchmark and prospective sequencing validation remain separate evidence tasks.
The accompanying sequential comparison is exploratory and uses per-stage budgets.

Validation completed 20 September 2026: the full test suite passed with 5,599
passed and 26 skipped (5 warnings). The focused shared-service, contraction,
chemistry-contract and function-size checks also passed (22 tests). Repository
formatting/import checks and lint checks for the new modules passed.

The [Wolbachia sequential comparison](sequential_panel_service_2026-09-19.md)
found a 9-oligo subset from a 24-oligo proposal at 70.07% predicted effective
coverage, under the explicitly configured density and dimer limits. This is a
bounded-search result, not evidence of minimum size or measured genome recovery.


Original audit: 19 September 2026 at `c3cf80c`, before the implementation
follow-up above.

## Answers

- **One shared optimizer interface exists, but not one complete orchestration
  entry point.** `BaseOptimizer`/`OptimizerFactory` unify construction and result
  metrics. Standard optimization, pool planning and expansion still assemble
  different workflows and configurations.
- **Sequential optimization is appropriate and already partly implemented.**
  Hybrid performs screening, initial coverage selection and refinement. Pool
  planning adds size enumeration, candidate-frontier expansion and repair.
  These should become stages of one service, preserving the best qualifying
  incumbent under one explicit objective.
- **Not every method uses chemistry at every selection stage.** Correct
  condition-aware evaluation does not mean that candidate ranking, proposals,
  refinement and size reduction all optimize that evaluation. Sequence-based
  compatibility rules are intentionally independent of additives. The optional
  free-energy dimer check uses temperature, not full buffer/additive chemistry.

## Validation and scope

75 focused tests passed in 1.43 seconds:

```bash
python -m pytest tests/test_optimizer_chemistry_consistency.py tests/test_stage_one_selects_on_the_judged_quantity.py tests/test_clique_optimizer.py tests/test_clique_fixed_primers_and_defaults.py tests/test_the_deficit_objective_drives_selection.py -q
```

These cover canonical Tm agreement under changed buffer/additive inputs,
stage-one objective routing, clique behavior and deficit-driven selection.
This review also traced the current command handlers and registered methods.
It did not rerun large-genome benchmarks or establish experimental validity of
the chemistry or recovery models. Earlier audit findings are not presumed to
remain open when their implementations have changed.

Additional direct probes confirmed that a params mapping with
`stage1_objective_width=64` resolves to `None` through the planner's context/config
construction, and that the default `stage1_objective` returns `None` even with
an attached objective. A canonical/network Tm comparison under two DMSO inputs
agreed in both cases; the remaining issue is where that calculation influences
selection, not disagreement between those successful Tm calls.

## Current entry points

| Entry | Actual route | Distinction |
|---|---|---|
| `optimize` | CLI pipeline -> `optimize_step4` -> `run_optimization` -> factory -> method | Resolves much configuration through global parameters; default refinement is network; optional limits are applied again after selection |
| Python configuration API | `run_optimization_from_config` -> `run_optimization` | Shares the standard orchestration, but its configuration object is not the complete `DesignContext` |
| `plan-pool` | `run_plan_pool` -> `DesignContext` -> factory -> `plan_pool` | Separate orchestration; attaches `PoolObjective`, defaults to swap refinement, evaluates sizes and can refill candidate frontier |
| `plan-pool --design-grid` | Grid handler -> condition reassessment/sweep -> factory/planner per condition | Now connected, unlike the 16 September audit; settings must still be forwarded consistently |
| `expand-primers` | Iteration handler/core expander -> context -> hybrid or adapter | Shares resolved chemistry and adds a deficit objective; remains its own orchestration |
| `ensemble` | Standard dispatcher -> several factory methods -> winner/combination | Alternative proposal generation, not a sequential chain optimizing one objective |

Key files: `neoswga/core/unified_optimizer.py`, `design_context.py`,
`optimizer_factory.py`, `pool_planner.py`, `primer_expansion.py` and
`neoswga/cli/plan_pool.py`.

## Findings

### 1. P1: configured chemistry can still be replaced after an error

`run_optimization` at `neoswga/core/unified_optimizer.py:919` catches condition
construction errors and continues with `conditions=None`. It logs the loss
of chemistry and later adds a validator message, but still attempts design.
The `DesignContext` path used by pool planning raises instead.

`NetworkOptimizer._get_primer_tm` at `network_optimizer.py:736` correctly uses
canonical Tm and condition-aware cache keys on success. On an exception it
warns and substitutes a fixed-buffer estimate without additives.

**Consequence:** identical invalid chemistry can abort one command and produce
a fallback design through another. A warning does not make that design obey
the requested reaction.

**Recommendation:** fail new design requests when explicitly requested
chemistry cannot be evaluated. If a legacy approximation is needed, require
an explicit mode and record the actual model used. Do not change conditions
halfway through a run.

### 2. P1: stage-one chemistry scoring is optional, and pool planning drops its setting

`stage1_objective` in `neoswga/core/swap_refinement.py:55` returns `None` when
`config.stage1_objective_width` is unset, **before** checking for an attached
pool objective. The default is `None`. Thus the default dominating-set/hybrid
stage-one search still uses unweighted coverage bins. Setting a width enables
full objective evaluation only for a shortlist ranked by geometric gain.

This default is deliberate and benchmarked, not an unimplemented helper:
[stage-one benchmark](stage_one_objective_2026-09-19.md) documents coverage,
specificity and runtime tradeoffs. Preserve that decision in compatibility mode
rather than silently enabling a different objective for everyone.

There is nevertheless a propagation defect: `run_plan_pool` constructs its
config through `context.optimizer_config` with refinement and swap-budget
overrides, but does not forward `params['stage1_objective_width']`.
`DesignContext` does not carry that field. A setting honored by standard
`optimize` therefore remains unset in ordinary `plan-pool` configuration.

**Recommendation:** forward all supported search controls through a shared
request/configuration builder. Describe raw versus effective stage-one
selection explicitly. Offer a consistent objective mode with configured
specificity constraints, and expose the shortlist width as an approximation
budget. A production-command test should assert the value reaching the inner
selector, not only test a factory instance with an explicit config.

### 3. P1: methods do not share one selection/refinement objective

| Method/stage | Chemistry use | Remaining difference |
|---|---|---|
| Dominating-set | Shared final effective metrics; effective stage-one scoring when width is configured | Default selection is geometric; cheap shortlist is geometric even in effective mode |
| Hybrid | Conditions reach secondary screening and network Tm; can use shared stage-one/swap objectives | Standard default uses network refinement, which does not read `pool_objective` |
| Background-aware | Conditions reach the inner hybrid; background terms and optional shared objectives | Inherits hybrid's stage/mode differences; its name does not imply full chemical modeling |
| Network | Canonical Tm and optional mechanistic terms | Tm/mechanistic weights default to zero at the adapter; standard application profiles can override some weights. Merely passing conditions does not enable these terms |
| Clique | Full shortlisted-panel metrics use conditions and normalized score | Candidate cap, compatibility graph and cheap set ranking are sequence/count based; it does not choose by the attached planner objective |
| Planner swaps/repair | Shared coverage/constraint objective is attached to wrapper and inner hybrid | This is substantially improved; it does not make every alternative proposal generator optimize that objective |
| Deficit expansion | Resolved chemistry plus attached deficit objective; swap refinement can use it | A distinct request/acceptance workflow still needs unification with ordinary design |

Clique currently defaults to 200 candidates, 10,000 enumerated sets and 100
fully scored sets (`clique_optimizer.py:58`). Its winner maximizes
`metrics.normalized_score()`, after cheap ranking by exact binding ratio.
That can differ from maximizing effective coverage subject to a density limit
in `plan-pool`. The caps are computational approximations, not completeness
guarantees. Network proposals similarly pursue connectivity and their configured
terms rather than the planner's complete objective.

**Recommendation:** treat clique/network as proposal generators unless they
explicitly implement the common objective. Accept or reject their proposals
through one evaluator and constraint contract. Replace the coarse
`ADDITIVE_AWARE` boolean with capabilities such as condition-aware QC,
condition-aware selection, mismatch model and dimer model. Reporting an
effective metric is different from selecting on it.

### 4. P2: pool minimization changes the coverage definition

`_minimize_primer_count` at `unified_optimizer.py:720` tests and ranks deletions
using `metrics.fg_coverage`, the geometric metric. It does not directly preserve
the effective coverage or the whole constraint set while deleting. A later
configured-limit check can detect/repair some failures, but it does not turn
the deletion search into constraint-preserving minimization.

**Recommendation:** use the same objective/constraints for add, swap and drop.
Retain the previous feasible panel if shrinking fails. Search pool size as a
frontier, and describe the result as the smallest qualifying pool found.

### 5. P2: dimer checks are not a full modified-chemistry model

The strict base-pair complementarity limit is a sequence rule. It is reasonable
for that rule to stay unchanged when additives change; do not imply that
additives automatically permit its relaxation.

The optional stability check goes through `LazyDimerCompatibility` and
`dimer.is_dimer_thermodynamic`. The latter explicitly uses only temperature
from conditions, and evaluates the longest complementary region. It does not
apply the configured salts/additives to a full primer-primer equilibrium model.
Clique's compatibility graph uses sequence-based checks; its search does not
use this optional stability floor in graph construction. Final validation and
proposal compatibility therefore need an explicit capability distinction.

**Recommendation:** label the implemented models accurately and propagate the
same configured constraints to all stages. A chemically detailed dimer model
would require a separately specified and validated implementation, not reuse
of a duplex Tm correction as if it established dimer equilibrium.

## Recommended sequential architecture

Keep the existing CLI commands as convenient front ends, but have them build
one `OptimizationRequest` and call one design service. The request should carry
resolved conditions, reference/index provenance, candidate source, fixed and
excluded primers, objective/constraints, pool-size bounds, stage budgets and
seed. Conditions and model identity must remain fixed within a run.

The service can execute these stages:

1. **Validate and resolve:** chemistry, current QC eligibility, indexes and
   capabilities. Refuse unsupported mandatory constraints before selection.
2. **Build initial proposals:** bounded geometric or effective greedy, with
   optional network/clique proposals. Record which approximation generated
   each proposal and which candidates were examined.
3. **Evaluate under one objective:** recompute every proposal with the shared
   conditions and hard constraints. Maintain feasible incumbents separately
   from promising but infeasible search states.
4. **Improve:** constrained swaps and additions; refill the candidate frontier
   on unmet targets or a declared improvement-search policy. Optional network
   proposals must not overwrite a better feasible incumbent automatically.
5. **Minimize:** remove redundant primers while preserving the same selected
   coverage metric and all constraints. Revisit swaps after removals.
6. **Compare sizes and conditions:** retain nondominated alternatives and the
   smallest qualifying found pool. Each chemistry is a separate validated run;
   changing chemistry invalidates condition-dependent assessments and caches.
7. **Validate and report:** recompute final metrics, check the complete delivered
   pool, and record stage deltas, budgets, stop reasons, model identity and
   raw/effective coverage separately.

For sequencing feedback, substitute a deficit-priority objective at stages
3-5 and preserve the same configuration and acceptance infrastructure. The
observed coverage profile remains evidence from the prior experiment; changes
to a proposed pool remain predictions until measured.

The benefit of a sequence of stages is complementary search behavior, not a
guarantee that every later method improves the previous answer. Use one
acceptance contract and preserve the best feasible incumbent at every stage.

## Next implementation order

1. Fix condition-error fallbacks and the dropped stage-one setting; add CLI
   parity tests across standard optimization, planning and expansion.
2. Introduce the shared request/service while preserving current modes and
   defaults. Make stage objective and method capabilities visible in reports.
3. Route proposal acceptance, refinement and minimization through the same
   objective. Add small oracle tests for chemistry-dependent ranking and
   constraint-preserving add/swap/drop behavior.
4. Benchmark the sequential mode against existing modes at comparable budgets
   and declared specificity limits. Keep coverage, specificity, pool size and
   runtime as separate outcomes; do not choose a winner from coverage alone.
