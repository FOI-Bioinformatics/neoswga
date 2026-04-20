# Phase 11+ — Primer-Set Optimization Quality

**Date:** 2026-04-20 (second plan this day)
**Target:** Production-grade primer-set optimization, not just plumbing
**Prior work:** Phases 1-10 merged; ReactionConditions reaches every optimizer, scientific constants under regression test, `swap-primer` / `contract-set` / `rescore-set` CLIs wired, lab-result feedback loop documented.

## Why a third plan

The user's emphasis — *"focus also on the primer set optimisations"* — flagged that plumbing and ergonomics are done, but the actual **selection logic** across the 16 registered optimizers has qualitative gaps. A third audit pass surfaced the following.

### Audit findings (condensed)

1. **Conditions-awareness in scoring is uneven.** After Phase 7 every optimizer *accepts* `self.conditions`; only `network_optimizer`, `genetic_algorithm`, and `equiphi29_optimizer` actually *use* it in scoring. `greedy`, `dominating-set`, `weighted-set-cover`, `hybrid` (Stage 1), `milp`, `clique`, `tiling`, `normalized`, `moea`, `multi-agent`, `serial-cascade` treat conditions as metadata only. A user running `optimize -m dominating-set` with `betaine_m=1.5` vs `0.0` gets **identical** primer sets.

2. **Cross-optimizer scores are apples-to-oranges.** `hybrid.score` returns network-connectivity values (~5-50), `greedy.score` returns coverage fractions (0-1), `milp.score` returns region counts (100-5000), `moea.score` is a single best from the Pareto front. Only `PrimerSetMetrics.normalized_score()` is comparable — and it isn't populated uniformly.

3. **Pan-genome coverage is not enforced per target.** In multi-target mode (`fg_prefixes=[A,B,C]`), every optimizer aggregates coverage into one number. A set with 99% coverage on target A and 40% on targets B and C passes with "99% fg_coverage". Clinical pan-primer designs fail silently.

4. **`swap-primer` improvement unproven.** Current test asserts JSON shape, not that `after_metrics.mean_severity < before_metrics.mean_severity`. Could silently degrade a set.

5. **`contract-set` is quality-blind.** My Phase 8B implementation removes primers purely by coverage marginals. When two primers contribute equally to coverage, it picks arbitrary — should prefer to remove the one with worse Tm fit / higher dimer risk / lower scoring.

6. **`rescore-set` is half-measure.** Phase 8C computes per-primer quality scores under new conditions but **does not recompute coverage, background hits, or set-level dimer energies**. User asking "how does my validated phi29 set perform at equiphi29 + 1 M betaine?" gets per-primer Tm shifts and that is all.

7. **Reproducibility has MOEA gap.** `--seed` flag exists and is applied globally in `unified_optimizer`, but `pymoo` inside `moea_optimizer.py` does not receive the seed, so same-seed runs of MOEA can diverge.

8. **No post-optimization sanity check.** `OptimizationResult` is a frozen dataclass but nothing validates no-duplicates, `len(primers) == target_size`, `coverage > 0` before returning to the user. A broken optimizer silently returns empty or duplicate-heavy sets.

9. **`--use-mechanistic-model` is a placeholder.** Flag exists at `cli_unified.py:~850`; at `:2114` it has a "not yet integrated" warning. Zero effect on selection today.

10. **Temperature-dependent dimer stability untested.** `StructurePrediction(conditions)` supposedly uses `conditions.temp`, but no regression test asserts `severity(30 C) > severity(45 C)` for the same primer pair.

11. **Grade vs score alignment untested.** `results_interpreter.py` grades A-F; `OptimizationResult.normalized_score()` returns 0-1. No test links them. User might see "GOOD" at `normalized_score = 0.55`.

12. **MOEA Pareto front is collapsed.** NSGA-III computes the non-dominated front internally but `moea_optimizer.optimize()` returns a single best solution. `--show-frontier` exists at CLI level but needs an actual frontier to visualize.

13. **`neoswga doctor` never built.** Listed in the previous plan's Phase 10; deferred.

Nothing above is catastrophically broken, but a production-quality tool for multi-target, multi-background, additive-aware, pan-genome primer design needs each of these to land. Plan phases 11-14 close them.

---

## Phase 11 — Optimizer Scoring Correctness & Coverage Parity

### Task 11A: Every optimizer's scoring honours `self.conditions`

For each of the following optimizers, locate the inner scoring helper and replace any buffer-agnostic Tm / dimer / secondary-structure calculation with its ReactionConditions-aware counterpart (`thermodynamics.calculate_tm_with_salt` + `conditions.calculate_tm_correction`, or `StructurePrediction(conditions)`):

- `neoswga/core/greedy_optimizer.py` — Tm window check + dimer rejection.
- `neoswga/core/dominating_set_optimizer.py` — add Tm window check when conditions provided.
- `neoswga/core/hybrid_optimizer.py` Stage 1 (dominating-set phase) — same.
- `neoswga/core/milp_optimizer.py` — add a Tm-window constraint parameter from conditions.
- `neoswga/core/clique_optimizer.py` — dimer edges should use `StructurePrediction(conditions)`.
- `neoswga/core/tiling_optimizer.py` — Tm filter.
- `neoswga/core/normalized_optimizer.py` / `multi_agent_optimizer.py` — forward to inner optimizer (already done if sub-optimizer is patched).

Acceptance: a parametric ball-bearing test over the 10 CLI-exposed methods shows same-pool + different additives produces measurably different top-K for **every** method (tolerate no-change only when the pool is empty after filter).

### Task 11B: `normalized_score` is populated by every optimizer

`OptimizationResult.normalized_score()` (`base_optimizer.py`) already defines the formula. Ensure every optimizer populates the inputs (fg_coverage, selectivity, dimer_risk, gap_gini, tm_spread). Today several optimizers leave some fields at 0.0, biasing the normalized score. Audit with a test that runs each optimizer on the plasmid example and asserts every field is populated, not a stand-in zero.

### Task 11C: Output column parity — enforced via a single writer

`unified_optimizer.save_results()` already emits the canonical columns. Add a regression test asserting the union of columns across **every** optimizer is identical on the plasmid example. Today the network optimizer emits `strand_alternation_score` non-zero while others fallback to 0.0 — that is acceptable as long as the column exists. Catch regressions where an optimizer skips a column entirely.

### Task 11D: Per-target coverage enforcement for multi-genome

Add a `min_per_target_coverage: float = 0.50` parameter to `OptimizerConfig`. In `OptimizationResult`, add a `per_target_coverage: Dict[str, float]` field populated by every multi-genome-aware optimizer. Post-optimization, if any target is below `min_per_target_coverage`, mark the result `PARTIAL` and log a warning with the per-target breakdown. Wire `--min-per-target-coverage` to `optimize` and `swap-primer`. Update `interpret` so the grading reflects per-target minimums, not only aggregate coverage.

Files: `base_optimizer.py`, `optimizer_factory.py`, `unified_optimizer.py`, every optimizer that iterates `fg_prefixes`, `results_interpreter.py`.

---

## Phase 12 — Primer-Set Quality Assurance

### Task 12A: Post-optimization validator

Add `OptimizationResult.validate() -> ValidationReport`. Checks:

- No duplicate primers (`len(set(primers)) == len(primers)`).
- `len(primers) == target_set_size` or status is `PARTIAL`.
- `fg_coverage > 0` or status is `PARTIAL` / `FAILURE`.
- For multi-target runs, `all(v >= min_per_target_coverage for v in per_target_coverage.values())` or `PARTIAL`.
- If `conditions` provided, every primer's effective Tm is within `[min_tm, max_tm]` as adjusted by additives; otherwise warn.
- If `bl_genomes` present in params, no primer hits a blacklist k-mer (re-check via cache; guards against accidental re-injection via `expand-primers`).

Called automatically at the end of `run_optimization`. Output stored on `OptimizationResult.validation` and written as `step4_improved_df_validation.json`.

### Task 12B: `swap-primer` strict-improvement regression

Augment `tests/test_current_set_improvement_cli.py` with a test that builds a deliberately dimer-heavy set (using known complementary primer pairs), runs swap-primer, and asserts:

- `after_metrics.mean_severity <= before_metrics.mean_severity - 0.01`
- `after_metrics.max_severity <= before_metrics.max_severity`
- No primer in after_set is a duplicate of another.
- After set size matches before set size (net-zero swap).

Also add an inner regression in `dimer_network_analyzer.optimize_set_greedy`: never commit a swap that increases `mean_severity` even transiently.

### Task 12C: Quality-weighted `contract-set`

Rewrite `run_contract_set` (`cli_unified.py`) so that when choosing which primer to drop, it ranks candidates by a **quality deficit** score (higher = first to remove):

```
deficit(primer) = w_cov * (1 - coverage_contribution) +
                  w_dimer * dimer_risk +
                  w_tm * |tm_effective - target_tm| / tm_window +
                  w_bg * bg_freq
```

Weights default to `w_cov=0.40, w_dimer=0.25, w_tm=0.20, w_bg=0.15`. Only remove if post-removal coverage stays >= `--min-coverage`. Keep iterating until no candidate passes the constraint. Add a regression test that shows a set with one obvious low-quality primer has exactly that primer removed first.

### Task 12D: `rescore-set` recomputes coverage + bg + set-level dimer energies

Extend `run_rescore_set` to:

- Compute `fg_coverage` and `per_target_coverage` by querying `PositionCache` for each primer and applying the condition-aware extension reach (polymerase processivity from `ReactionConditions.polymerase`).
- Compute `bg_coverage` from `bg_prefixes` if present.
- Compute set-level dimer energies with `StructurePrediction(conditions)` so a dimer that was stable at 30 C but melts at 42 C gets a low severity under equiphi29 conditions.
- Emit a side-by-side "current-conditions" vs "requested-conditions" table when the user passes both the `-j params.json` and scenario overrides.

### Task 12E: Temperature-dependent dimer regression

New test in `tests/test_dimer_temperature.py`:

- Define a primer pair known to form a stable dimer at 30 C (`AAACCCGGG` + `CCCGGGTTT` or similar with strong complementarity).
- Instantiate `StructurePrediction(ReactionConditions(temp=30, polymerase='phi29'))` and compute severity.
- Repeat at `temp=45, polymerase='equiphi29'`.
- Assert severity drops by a meaningful amount (say ≥ 30%).

If the current implementation does not already show this, fix the temperature usage inside `secondary_structure.py` / `dimer_network_analyzer.py`.

---

## Phase 13 — Reproducibility, Diagnostic, and the `--use-mechanistic-model` promise

### Task 13A: MOEA actually honours `--seed`

`moea_optimizer.py` calls `pymoo.optimize.minimize(...)`. Pass `seed=self.config.seed` into `minimize()`. Add a reproducibility test:

```
run1 = optimize(... method='moea', seed=42)
run2 = optimize(... method='moea', seed=42)
assert run1.primers == run2.primers
```

Do the same for any other optimizer that uses randomness (hybrid Stage 2 network tie-breaking, genetic if not already seeded).

### Task 13B: Decide the fate of `--use-mechanistic-model`

Two options. Pick one, land it, remove ambiguity:

1. **Integrate**: when the flag is set, `unified_optimizer.run_optimization` builds a `MechanisticModel(conditions)`, computes per-primer `predicted_amplification_factor`, and passes a `mech_weight * predicted_amplification` term into every conditions-aware optimizer's scoring sum (weight from `--mechanistic-weight`, default 0.3).
2. **Deprecate**: remove the flag; document that the mechanistic model is available via `predict-efficiency` only; bump the existing `DeprecationWarning` when the flag is set.

Recommendation: **option 1**. The scientific value is real, and the machinery (`mechanistic_params.py`, `mechanistic_model.py`) is already calibrated. Integration effort is modest (~3 h).

### Task 13C: `neoswga doctor`

New CLI `neoswga doctor` that reports:

- Python + neoswga version, jellyfish location and version.
- Which optimizers are registered, with a column "additive-aware: yes/no" computed at runtime by introspecting whether their scoring path references `self.conditions`.
- Current `params.json` summary (resolved polymerase, reaction_temp, additives, Mg²⁺, GC-adaptive state, bl/bg/fg prefix existence on disk, k-mer count file existence).
- If provided `--primers`, run the post-optimization validator against that set plus the current params.
- Non-zero exit if any "blocking" issue found (missing jellyfish, missing k-mer files, invalid params).

### Task 13D: Grade-vs-score alignment test

`tests/test_interpret_vs_score.py`: for a handful of known primer sets (plasmid example top-3, a deliberately poor set, a synthetic ideal set), assert that `results_interpreter.rate_metric(...)` A-grade sets have `normalized_score() > 0.80`, B-grade > 0.65, C-grade > 0.50. Fail if there is crossover (e.g., an A-grade set with normalized_score 0.55).

---

## Phase 14 — Pareto Frontier & Application Consistency

### Task 14A: MOEA returns the full Pareto front

Extend `OptimizationResult` with `pareto_front: Optional[List[List[str]]] = None`. `moea_optimizer.optimize()` populates it when NSGA-III returns non-dominated solutions. Other optimizers leave it `None`.

### Task 14B: `--show-frontier` renders a tradeoff table

Today `--show-frontier` exists but only plots if Plotly is installed. Make it emit a plain-text table even without Plotly:

```
#  coverage  specificity  dimer_risk  tm_spread  normalized_score
1   0.98       87.5          0.12        3.4         0.81
2   0.94       112.3         0.08        2.9         0.83
3   0.91       142.7         0.06        2.6         0.82
...
```

Wire to `moea.pareto_front` when available; otherwise enumerate results from a multi-start `hybrid` run with varying `uniformity_weight`.

### Task 14C: `--application` consistency across optimizers

Today `--application` only affects `--auto-size` (set-size optimizer). Extend it to tune the `normalized_score` weights per application:

| Application | Coverage weight | Selectivity | Dimer | Gini | Tm |
|---|---|---|---|---|---|
| discovery | 0.50 | 0.15 | 0.10 | 0.15 | 0.10 |
| clinical | 0.20 | 0.45 | 0.20 | 0.10 | 0.05 |
| enrichment | 0.35 | 0.30 | 0.15 | 0.10 | 0.10 |
| metagenomics | 0.60 | 0.10 | 0.10 | 0.15 | 0.05 |

Default: `balanced` = current weights (0.35 / 0.30 / 0.15 / 0.10 / 0.10). Store per-application weights in `optimizer_config.normalized_weights`. Pass to `OptimizationResult.normalized_score` to compute application-appropriate rankings.

---

## Critical files

| Path | Phase | Action |
|---|---|---|
| `neoswga/core/base_optimizer.py` | 11B, 12A, 14A | Add `per_target_coverage`, `pareto_front`, `validation`, weighted `normalized_score` |
| `neoswga/core/unified_optimizer.py` | 11A, 11D, 12A, 13A, 14B | Post-opt validator hook, `--application` → weights, seed forwarding |
| `neoswga/core/{greedy,dominating_set,clique,tiling,normalized,multi_agent,moea,milp}_optimizer.py` | 11A | Add condition-aware Tm / dimer in scoring |
| `neoswga/core/moea_optimizer.py` | 13A, 14A | Forward seed; export full Pareto front |
| `neoswga/core/dimer_network_analyzer.py` | 12B, 12E | Strict-improvement guard; temperature-dependent severity |
| `neoswga/cli_unified.py` | 12C, 12D, 13B, 13C, 14B | Rewrite contract-set; fatten rescore-set; wire `doctor`; frontier table |
| `neoswga/core/secondary_structure.py` | 12E | Confirm temperature is applied to dimer severity |
| `tests/test_optimizer_additive_propagation.py` | 11A | Parametric "additives change set" per optimizer |
| `tests/test_per_target_coverage.py` | 11D | Multi-genome per-target assertion |
| `tests/test_swap_primer_improvement.py` | 12B | Strict after < before |
| `tests/test_contract_set_quality_weighted.py` | 12C | Low-quality primer is removed first |
| `tests/test_rescore_set_coverage.py` | 12D | Coverage and bg present in output |
| `tests/test_dimer_temperature.py` | 12E | Severity drops at higher temp |
| `tests/test_reproducibility_moea.py` | 13A | Same seed → same set |
| `tests/test_interpret_vs_score.py` | 13D | Grade vs normalized_score |
| `docs/optimizer-quality.md` | 11-14 | Explain the score semantics, application weights, and assurance matrix |

## Verification

After phases 11-14 are merged:

1. `pytest tests/ -m "not scale" --timeout=300` stays green (target: 2100+).
2. `pytest tests/integration/test_additive_ranking_e2e.py` now parametrizes all CLI-exposed optimizers and each passes.
3. `pytest tests/test_scientific_constants.py tests/test_dimer_temperature.py` green — scientific correctness stays tight.
4. `neoswga doctor` on a fresh checkout with the plasmid example reports "ready".
5. `neoswga optimize -m moea --seed 42 --show-frontier` emits a deterministic Pareto table; rerun produces identical output.
6. `neoswga swap-primer` on a known-bad set produces strictly lower `mean_severity`; contract-set removes the lowest-quality primer first; rescore-set output carries `fg_coverage` and `per_target_coverage`.
7. Cross-optimizer ranking works: running 5 methods on the same pool and sorting by `normalized_score` gives a consistent winner across re-runs with `--seed 42`.

## Risk register

| Risk | Mitigation |
|---|---|
| Post-opt validator's blacklist re-check is slow on large pools | Skip if `bl_prefixes` is empty; reuse in-memory PositionCache |
| Integrating mechanistic model biases scoring away from good practical sets | Keep default `--mechanistic-weight=0.3` (minority share); add off-switch |
| MOEA seed wiring may not expose pymoo's internal RNGs | Verify via same-seed re-run test; if unstable, add a `@pytest.mark.flaky(reruns=2)` escape hatch only for the seed test |
| `--application` weight changes affect grading and cause surprise | Document migration note for 3.7 → 3.8; default = balanced preserves current behaviour |

## Effort

| Phase | Effort |
|---|---|
| 11 (scoring correctness) | 2-3 days |
| 12 (quality assurance) | 2 days |
| 13 (reproducibility + doctor + mech-model) | 1.5 days |
| 14 (Pareto + application) | 1 day |
| **Total** | **6-7 working days** |

Phase 11A is the highest-leverage single task: once every optimizer's scoring actually uses `self.conditions`, the "additives change the selected set" claim holds across the full CLI surface instead of just hybrid/genetic/equiphi29.
