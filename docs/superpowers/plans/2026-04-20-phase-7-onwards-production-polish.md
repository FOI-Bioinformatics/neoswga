# Phase 7+ Production Polish — Optimizer Plumbing, Iterative Improvement, Scientific Validity

**Date:** 2026-04-20
**Relationship to prior plans:** Builds on Phase 1-6 (release blockers through example library). Phases 1-6 are merged (commits 6ca0528...e3f45fc).

## Context

A second round of deep exploration surfaces three dimensions the first round did not cover, each a production blocker:

1. **The critical plumbing bug**: `unified_optimizer.run_optimization()` and `optimize_step4()` never construct `ReactionConditions` from `parameter.*` globals and never pass them to `OptimizerFactory.create()`. The session-1 work wired additive-aware Tm into `integrated_quality_scorer._primer_tm` and `network_optimizer._get_primer_tm`, but with the conditions never reaching those objects at runtime, the user-facing behaviour is unchanged. A user running `neoswga optimize -j params.json` with `betaine_m: 0` vs `betaine_m: 1.0` gets **identical** primer sets today.

2. **No iterative-improvement CLI for existing sets**. `expand-primers`, `predict-efficiency`, and `active-learn` exist and are wired, but there is no CLI for:
   - swapping one primer for a better candidate (`dimer_network_analyzer.optimize_set_greedy` exists but is Python-API only)
   - contracting a set while keeping coverage
   - rescoring an existing set under arbitrary reaction conditions
   - a consolidated `improve-set` meta-command

3. **Scientific validity of empirical constants**: the Owczarzy salt correction, SantaLucia NN parameters, and literature-cited additive coefficients at 37 C are solid (grade A+). But several *empirical* constants are **uncited in code and not regression-tested**:
   - Betaine Arrhenius activation energy 1800 J/mol (`mechanistic_params.py:66`)
   - Trehalose coefficient -3.0 C/M (no multi-temp data)
   - Betaine peak enhancement 1.0 M (literature says 5.2 M for full equalization)
   - TMAC uniform coefficient -0.5 C/M (literature says primarily GC-dependent)
   - DMSO threshold for phi29 5% (empirically tuned)

### The scenario matrix the user cares about

| Axis | Status after phases 1-6 | Status after phase 7+ (target) |
|---|---|---|
| Phi29 + equiphi29 scoring | Additive-aware in scorer, not reaching optimizer | Reaches every optimizer end-to-end |
| Different additive cocktails | Filter Tm honours additives; optimizer ranking does not | Optimizer ranking reflects every additive |
| Bacterial size | Pipeline works; no current-set resize | `contract-set` / `expand-primers` round-trip |
| GC content | Auto-adaptive filter; scorer honours extreme GC | Rescoring works under arbitrary GC targets |
| Multi-target / bg / blacklist | Multi-genome plumbing end-to-end | Rescoring works across multiple targets |
| New-set design | Works | Same, but optimizer now honours additives |
| Current-set improvement | Only `expand-primers` | `improve-set` meta with swap/expand/contract/rescore |
| Lab-result feedback | `active-learn` exists; schema undocumented | Documented schema; end-to-end tutorial; calibration CSV format |
| Scientific validity | A- | Regression-tested against primary literature; citations in code |

---

## Phase 7 — Optimizer Plumbing & Additive-Aware Ranking

### Task 7A (critical): Plumb ReactionConditions through `unified_optimizer`

`neoswga/core/unified_optimizer.py:208-373` (run_optimization) and `:586-677` (optimize_step4) currently pass `**kwargs` to `OptimizerFactory.create()` but never construct `ReactionConditions` from `parameter.*`. Fix:

```python
# Add near the top of run_optimization (~line 240)
from neoswga.core.reaction_conditions import ReactionConditions
from neoswga.core import parameter
conditions = ReactionConditions(
    temp=getattr(parameter, 'reaction_temp', 30.0),
    polymerase=getattr(parameter, 'polymerase', 'phi29'),
    na_conc=getattr(parameter, 'na_conc', 50.0),
    mg_conc=getattr(parameter, 'mg_conc', 10.0),
    dmso_percent=getattr(parameter, 'dmso_percent', 0.0),
    betaine_m=getattr(parameter, 'betaine_m', 0.0),
    trehalose_m=getattr(parameter, 'trehalose_m', 0.0),
    formamide_percent=getattr(parameter, 'formamide_percent', 0.0),
    ethanol_percent=getattr(parameter, 'ethanol_percent', 0.0),
    urea_m=getattr(parameter, 'urea_m', 0.0),
    tmac_m=getattr(parameter, 'tmac_m', 0.0),
)
# Then pass into the factory call:
optimizer = OptimizerFactory.create(
    name=method,
    ...,
    conditions=conditions,
    **kwargs,
)
```

Also reuse `filter.py:_get_reaction_conditions()` if available (DRY).

**Files:** `neoswga/core/unified_optimizer.py`, `neoswga/core/optimizer_factory.py` (accept + forward `conditions`).

### Task 7B: Make every optimizer accept `conditions` (even as no-op)

Goal: consistent signature. Pass-through is acceptable for optimizers that don't yet use additives; the explicit acceptance keeps the API honest and makes the next upgrade trivial.

| Optimizer | Today | Target |
|---|---|---|
| network_optimizer | Uses conditions | no change |
| integrated_quality_scorer | Uses conditions | no change |
| genetic_algorithm | Accepts conditions kwarg | Actually receive them from factory |
| hybrid_optimizer | Ignored | Accept + forward to inner NetworkOptimizer |
| dominating_set | Ignored | Accept (no-op for now) |
| background_aware | Ignored | Accept + forward to stage-3 HybridOptimizer |
| greedy_optimizer | Ignored | Accept + use additive-aware Tm in scoring |
| milp / moea / clique / tiling / normalized | Ignored | Accept as no-op |
| equiphi29_optimizer | Hardcoded 42 C | Read temperature from conditions |

**Files:** all `neoswga/core/*_optimizer.py` + `neoswga/core/optimize.py` (greedy).

### Task 7C: Greedy optimizer additive-aware Tm

`neoswga/core/optimize.py` is the original breadth-first greedy. When a new primer is considered for addition to a set, its Tm window check should use the effective Tm (with additives) rather than raw Wallace or basic Tm. Low-risk change: wrap Tm computation in a helper that falls back to legacy behaviour when `conditions is None`.

### Task 7D: Ball-bearing test — additives change the chosen set

New integration test in `tests/integration/test_additive_ranking_e2e.py`:

1. Use `examples/plasmid_example/` with pre-built k-mer counts.
2. Run `optimize` with `betaine_m: 0.0`, record top-5 primer set.
3. Run again with `betaine_m: 1.5`, record top-5.
4. Assert at least 1 primer differs OR the Tm distribution of the selected set shifts by >1 C.
5. Parametrise across optimizers that *should* honour conditions (network, hybrid, genetic, background-aware, integrated-scorer path).

This is the canonical acceptance test for Phase 7.

### Task 7E: Unit tests on filter.py `_get_reaction_conditions` reuse

Confirm that the `run_optimization` code reads from the same cached `_get_reaction_conditions()` helper as `filter.py`, so filter-time and optimize-time additives cannot drift. Add a regression test that changes `parameter.betaine_m` between runs and confirms the cache resets (which Phase 2 already wired — new test just locks in coverage).

---

## Phase 8 — Iterative Improvement of Existing Primer Sets

### Task 8A: `neoswga swap-primer` CLI

Expose the dormant `dimer_network_analyzer.optimize_set_greedy()` via CLI:

```bash
neoswga swap-primer \
    --primers SEQ1 SEQ2 SEQ3 ... \
    --candidates candidates.csv \
    --max-swaps 3 \
    --objective dimer|coverage|specificity|composite \
    -j params.json \
    -o improved_set.csv
```

Wraps `dimer_network_analyzer.py:335-436`:

- `identify_primers_to_replace(set)` ranks primers by replacement priority.
- `suggest_replacements(primer, candidates)` scores alternatives.
- `optimize_set_greedy(set, candidates, max_swaps)` iteratively swaps.

Output: JSON + CSV with before/after scores, delta-score per swap, and final improved set.

### Task 8B: `neoswga contract-set` CLI (new)

Given a current set `{p1..pN}` and a `min_coverage` threshold, find the smallest subset that still meets coverage. Approach: greedy leave-one-out. Keep a primer only if removing it drops coverage below threshold.

```bash
neoswga contract-set \
    --primers SEQ1 SEQ2 ... \
    -j params.json \
    --min-coverage 0.70 \
    --output reduced_set.csv
```

Reuses `PositionCache` and coverage calculation from `background_aware_optimizer._calculate_coverage()`.

### Task 8C: `neoswga rescore-set` CLI under arbitrary conditions

Given a primer set, compute its scores under any polymerase + additive cocktail:

```bash
neoswga rescore-set \
    --primers SEQ1 SEQ2 ... \
    --polymerase equiphi29 \
    --reaction-temp 43 \
    --betaine-m 1.0 \
    --dmso-percent 5.0 \
    -j params.json \
    -o rescored.json
```

Uses `IntegratedQualityScorer(conditions=...)` + `PositionCache` to produce per-primer and set-level scores with current conditions vs requested conditions side-by-side.

### Task 8D: `neoswga improve-set` meta-command

Thin CLI wrapper that runs `expand-primers` → `swap-primer` → `contract-set` → `rescore-set`, reporting delta at each stage. Honours `--only swap,contract` if user wants a subset.

### Task 8E: Lab-result schema documentation

Create `docs/active-learning-guide.md`:

- CSV format expected by `active-learn --experimental-results`:
  - `primer_set_id` (str), `primers` (semicolon-joined), `enrichment_fold` (float), `uniformity` (float 0-1), `coverage` (float 0-1), `median_depth` (float, optional).
- JSON format (current, keep for backward compatibility).
- End-to-end tutorial: design → synthesize → sequence → log results → re-design.
- Calibration report walkthrough using `ExperimentalTracker.get_calibration_report()`.

Add `scripts/csv_to_experimental_results_json.py` to convert lab-readout CSVs into the JSON format `active-learn` expects.

---

## Phase 9 — Calibration & Lab-Feedback Reality Check

### Task 9A: `neoswga calibrate` CLI

Ingest a lab-result CSV + the predictions that were made for those primer sets, and emit calibration diagnostics:

- Pearson / Spearman correlation between predicted and observed enrichment
- Systematic over/under-prediction per scoring component
- Recommended scaling factor if the model is consistently off

Bridges `ExperimentalTracker.record_outcome()` + `get_calibration_report()` into a user-facing command.

### Task 9B: `predict-efficiency` learns from calibration

When calibration data is available (via `--tracker-path`), apply the recommended scaling factor during prediction. Track via `ExperimentalTracker` so the next calibration is tighter.

### Task 9C: In-silico sanity regression test

Using the plasmid example, run `rescore-set` for the same primer set under:

- phi29 at 30 C, no additives
- equiphi29 at 43 C, betaine 1 M + DMSO 5 %

Assert (1) effective Tm shifts by ≥1 C, (2) accessibility / processivity factors shift per the mechanistic model, (3) the quality ranking within the set rearranges.

This acts as a "non-trivial plumbing" canary without requiring wet-lab data.

---

## Phase 10 — Scientific Validity Hardening

### Task 10A: In-code citations for empirical constants

For each of the five flagged values, add a `# SOURCE:` comment with the primary literature citation (or explicit `# EMPIRICAL: no primary literature, tuned from <run id / dataset>`):

1. `mechanistic_params.py:66` — betaine Ea 1800 J/mol.
2. `mechanistic_params.py:86`, `additives.py:498` — trehalose -3.0 C/M.
3. `mechanistic_params.py:218` — betaine peak enhancement 1.0 M.
4. `additives.py:552`, `mechanistic_params.py:142` — TMAC uniform -0.5 C/M.
5. `mechanistic_params.py:182` — DMSO threshold 5% for phi29.

Where literature doesn't yet back the value, document it is empirical.

### Task 10B: Consolidated citation doc

`docs/SCIENCE_CITATIONS.md` lists every scientific constant in the codebase with its file:line, value, primary citation (DOI), and any deviation from the literature. Acts as a single audit surface for reviewers.

### Task 10C: Numeric regression tests against primary literature

`tests/test_scientific_constants.py`:

- SantaLucia 1998 NN values (`thermodynamics.py:34-80`): assert exact match to Table 1.
- Owczarzy salt correction coefficient 0.368 (`thermodynamics.py:301-341`): assert.
- DMSO -0.55 C/%: test Tm(primer, dmso=5%) - Tm(primer, dmso=0%) ≈ -2.75 C ± 0.2.
- Betaine -1.2 C/M (at 50 % GC): test shift at 1 M.
- Formamide -0.65 C/%: test shift at 10 %.
- Plus sanity checks for the empirically-tuned constants from 10A (ensure they stay within their documented range).

These tests act as alarm bells if someone tweaks a coefficient by accident.

### Task 10D: `neoswga doctor` CLI (diagnostic)

One-shot diagnostic that reports:

- Which optimizer paths are additive-aware (pulls from a runtime capability matrix).
- Which reaction conditions your current `params.json` maps to (effective Tm window, Mg2+, etc.).
- Whether your `bl_genomes` counts exist on disk.
- Whether your genome_gc is auto-calculated or user-provided.
- Warnings if any scientific constant has been overridden from literature default.

Useful for support tickets and reproducibility.

---

## Verification Plan

After executing phases 7-10, run:

1. `pytest tests/ -m "not scale"` green (expect ~2100+ tests).
2. `pytest tests/ -m scale` green with Jellyfish installed.
3. `pytest tests/integration/test_additive_ranking_e2e.py -v` — confirms the Phase 7 ball-bearing.
4. End-to-end sanity: run `neoswga rescore-set` on the plasmid set under 5 scenarios (phi29, phi29+betaine, equiphi29, equiphi29+DMSO+betaine, extreme-GC profile). Log-inspect the output; Tm and ranking should differ per scenario.
5. `neoswga doctor` output documents a clean additive-propagation state.
6. `neoswga improve-set -j params.json --primers X Y Z` round-trips through expand → swap → contract → rescore.
7. `docs/SCIENCE_CITATIONS.md` cross-referenced against `tests/test_scientific_constants.py` — every constant in one appears in the other.

## Risk register

| Risk | Mitigation |
|---|---|
| Plumbing ReactionConditions breaks optimizers that coincidentally rely on kwargs names | Add `**kwargs` safely; drop unknown keys with a debug log. |
| Swap / contract CLIs produce inferior sets in edge cases | Never replace the input set in place; always emit to a new file + provide before/after scores. |
| Scientific regression tests flag noise-level drift | Use ±0.2 C tolerance on Tm deltas, ±0.05 on fraction-valued shifts. |
| Lab-result schema churn | Freeze CSV columns in a `schema_version` field; bump deliberately. |

## Estimated effort

| Phase | Effort |
|---|---|
| 7 (optimizer plumbing) | 1-1.5 days — small code, broad impact |
| 8 (current-set CLIs) | 2-3 days |
| 9 (calibration) | 1-2 days |
| 10 (science hardening) | 1 day |
| **Total** | **5-7 working days** |

Phase 7 Task A alone is half a day and unblocks every previous optimizer-level claim about additive awareness. That should land first.
