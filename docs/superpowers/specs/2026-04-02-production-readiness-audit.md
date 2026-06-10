# NeoSWGA Production Readiness Audit & Gap Analysis

**Date**: 2026-04-02
**Scope**: Full codebase audit for PyPI publication readiness
**Codebase**: 56,341 lines across 91 core modules, 1,914 tests across 70 test files
**Methodology**: 6 specialized audit agents covering packaging, error handling, bugs, tests, security, and API stability

---

## Executive Summary

NeoSWGA is a scientifically capable tool with solid foundations: safe subprocess handling, no shell injection, good optional dependency management, and comprehensive thermodynamics validation. However, it is **not yet ready for public PyPI release**. The audit identified **7 critical**, **17 high**, **30 medium**, and **23 low** severity findings across 6 audit dimensions.

### Release-Blocking Issues (Must Fix)

1. **License mismatch** -- pyproject.toml/code says MIT, LICENSE file is AGPL-3.0
2. **Pickle deserialization fallbacks** -- 4 sites silently fall back to unrestricted `pickle.load()`, defeating the safe_pickle module
3. **No public API definition** -- 91 flat modules with no `__all__`, no public/private separation
4. **Incorrect Tm calculation** in `adaptive_filters.py` -- scientific correctness issue affecting primer filtering
5. **GA optimizer crashes on macOS** -- `ProcessPoolExecutor` with bound methods fails to pickle
6. **Zero real integration tests** -- the integration test file contains `assert True` placeholders
7. **No params.json schema versioning** -- silent behavior changes between versions

### Production Readiness Score: C+

| Dimension | Grade | Summary |
|-----------|-------|---------|
| Packaging & Distribution | C | License mismatch, broken MANIFEST.in, placeholder emails |
| Error Handling | B+ | Good CLI error dispatch, but ParamValidator not wired into pipeline |
| Code Quality | C+ | Scientific core is solid; GA, filters have real bugs |
| Test Coverage | D+ | 1914 tests but critical paths untested; integration tests are stubs |
| Security | C | Safe subprocesses, but pickle fallbacks create exploitable paths |
| API Stability | D | No public API, global mutable state, version confusion |

---

## Critical Findings (7)

### C1. License Mismatch (Packaging)
- `pyproject.toml` declares MIT; `neoswga/__init__.py` declares MIT
- Actual `LICENSE` file contains AGPL-3.0
- **Impact**: Legal liability; PyPI policy violation; user confusion
- **Fix**: Reconcile -- choose one license and update all references

### C2. Pickle Fallback to Unrestricted Deserialization (Security)
- `background_filter.py:294`, `background_filter.py:411`, `efficiency_predictor.py:394` -- fall back to raw `pickle.load()` when `safe_load()` fails
- `rf_preprocessing.py:132-137` -- falls back via `unsafe_load_with_warning()`
- **Impact**: Arbitrary code execution via crafted pickle files that intentionally trigger `RestrictedUnpickler` errors
- **Fix**: Remove all fallbacks. If safe_load fails, raise the error.

### C3. No Public API Definition (API)
- No `__all__` in `neoswga/__init__.py` or `neoswga/core/__init__.py`
- 91 modules in `core/` with no public/private distinction
- **Impact**: Users cannot know which imports are stable; any refactoring risks breaking downstream code
- **Fix**: Define `__all__` with key symbols; prefix internal modules with `_`

### C4. Incorrect Tm Calculation in adaptive_filters.py (Code Quality)
- `adaptive_filters.py:149-185` -- `_nn_tm()` reimplements nearest-neighbor Tm with incorrect entropy values and missing primer concentration term `R * ln(C/4)`
- Used in `ImprovedPipeline` phase 2 filtering
- **Impact**: Systematically incorrect Tm values affect which primers pass filtering
- **Fix**: Replace with `thermodynamics.calculate_tm_basic()` or `calculate_tm_with_salt()`

### C5. GA Optimizer ProcessPoolExecutor Pickling Crash (Code Quality)
- `genetic_algorithm.py:~341` -- submits `self._evaluate_individual` (bound method of complex object) to `ProcessPoolExecutor`
- Bound methods with `position_cache`, `dimer_matrix`, numpy arrays fail to pickle on macOS/Linux `spawn` start method (Python 3.8+ default on macOS)
- **Impact**: `PicklingError` crash for users selecting genetic algorithm optimizer
- **Fix**: Use module-level function or `ThreadPoolExecutor`

### C6. Zero Real Integration Tests (Testing)
- `tests/integration/test_integration.py:94` -- `assert True  # Replace with actual pipeline call when ready`
- The main user workflow (count-kmers -> filter -> score -> optimize) has never been tested end-to-end
- **Impact**: No confidence that the pipeline works as a whole
- **Fix**: Implement real integration test using plasmid_example data

### C7. No params.json Schema Versioning (API)
- No version field in params.json
- When defaults change between versions, old configs silently produce different results
- **Impact**: Reproducibility failures for published experiments
- **Fix**: Add `schema_version` field; emit warning for unversioned configs

---

## High Findings (17)

### Packaging (2)
- **H1. MANIFEST.in references missing files** -- `LICENSE.txt`, `ENHANCEMENTS_QUICKSTART.md`, `setup.py`, etc. Causes build warnings.
- **H2. Model file inclusion fragile** -- `neoswga/core/models/` has no `__init__.py`; `package-data` glob may not traverse it. Must verify with test build.

### Error Handling (2)
- **H3. ParamValidator not wired into pipeline commands** -- Type validation only runs via `validate-params` command, not during `count-kmers`/`filter`/`score`/`optimize`. Users who skip validation get raw Python errors from invalid params.json types.
- **H4. Empty filter results not caught at step2** -- If all primers are filtered out, an empty `step2_df.csv` is written. User must wait until step3 starts to get an error.

### Code Quality (5)
- **H5. `filter.py:170` division by zero** -- `get_bg_rates_via_bloom()` divides by `total` which is 0 when `primer_list` is empty.
- **H6. `genome_io.py:132` ZipFile handle leak** -- `_open_zip()` creates ZipFile but only the TextIOWrapper is closed.
- **H7. `BackgroundAwareOptimizer` only uses first foreground genome** -- `_calculate_coverage()` at line 360 uses `fg_prefixes[0]` only, underestimating coverage for multi-genome targets.
- **H8. Module-level `_reaction_conditions` cache never invalidated** -- `filter.py:24-49` caches conditions at first call; adaptive GC filtering changes params but doesn't reset the cache. Stale Tm calculations.
- **H9. `parameter._json_data` accessed without guard** -- `pipeline.py:441` lacks the `hasattr` check used elsewhere.

### Security (2)
- **H10. Hash verification bypassable** -- `rf_preprocessing.py:108` only verifies hashes for filenames in `_TRUSTED_MODEL_HASHES`. Different filename = no verification.
- **H11. Shipped pickle in PyPI package** -- `random_forest_filter.p` is inherently dangerous in distributed packages. Consider ONNX/skops migration.

### API Stability (4)
- **H12. Version duplicated in 5+ locations** -- `"3.0.0"` hardcoded in `__init__.py`, `pyproject.toml`, 3 report modules, and test fixtures. Version bump will leave stale copies.
- **H13. Version/changelog mismatch** -- `__version__` = 3.0.0 but CHANGELOG.md documents through 3.6.0.
- **H14. Global mutable state for parameters** -- `parameter.py` uses 40+ module-level globals. Breaks concurrent usage, testing, and library consumption.
- **H15. `get_params()` has surprising side effects** -- Creates directories, reads genome files, computes GC. Surprising for library users.

### Testing (2)
- **H16. 91 core modules, ~30 have zero test coverage** -- Including pipeline-critical `primer_attributes.py`, `thermo_background_filter.py`, `unified_optimizer.py`, `greedy_optimizer.py`.
- **H17. No deprecation policy** -- No documented timeline for deprecated features. Users cannot plan migrations.

---

## Medium Findings (30)

### Packaging
- M1. readme field lacks content-type specification
- M2. No upper bounds on optional deps (torch, plotly, tensorflow)
- M3. matplotlib/seaborn/openpyxl as core deps (add ~50MB install weight for optional features)
- M4. Tests/docs not excluded from sdist
- M5. CLAUDE.md included in sdist
- M6. README shows editable install only (no `pip install neoswga`)
- M7. No PyPI badges
- M8. Placeholder email addresses (`andreas.sjodin@example.com`)

### Error Handling
- M9. `get_value_or_default()` returns None for missing required params -- downstream crashes with unhelpful TypeError
- M10. `read_args_from_json()` no error handling for JSONDecodeError/FileNotFoundError
- M11. ZipFile handle leak in `genome_io.py`
- M12. `load_genome()` loads entire genome into memory -- OOM for human genome
- M13. Zero-length k-mer files pass step2 validation (existence check only, not size)
- M14. NaN/Inf values not checked in parameter validation
- M15. Missing timeout in `MultiGenomeKmerCounter.count_kmers_jellyfish()`

### Code Quality
- M16. `pipeline.py:24-101` O(P*K*F) scan -- reads entire k-mer files per primer
- M17. `genetic_algorithm.py:321` self-dimer check (range should start at i+1)
- M18. `pipeline.py:626` DataFrame view vs copy -- potential SettingWithCopyWarning
- M19. `adaptive_filters.py:379` loads entire genome for GC calculation
- M20. `improved_pipeline.py:311` doesn't handle empty FASTA
- M21. `utility.py:142` sigmoid overflow for large values
- M22. Multiple modules serve similar purposes (thermodynamics.py vs thermo_estimation.py vs melting_temp.py)

### Security
- M23. torch.load without `weights_only=True` in `deep_learning.py:607`
- M24. Shipped pickle model -- supply chain risk
- M25. No limits on k-mer range (`max_k=100` would exhaust resources)
- M26. No genome file size limits

### API Stability
- M27. Heavy deps imported at module level in ~20+ core modules
- M28. `thermo_estimation.py` emits DeprecationWarning on import
- M29. 16 optimization methods -- too many for public users
- M30. `--use-mechanistic-model` silently ignored with no deprecation warning

### Testing
- M31. Mechanistic model only tests Pathway 1 of 4
- M32. Network optimizer tests use only mock data
- M33. Excessive mocking in filter tests hides real behavior

---

## Gap Analysis: What's Missing for PyPI Publication

### 1. Legal & Compliance
| Gap | Current State | Required State | Effort |
|-----|---------------|----------------|--------|
| License consistency | MIT in code, AGPL-3.0 in LICENSE | One license everywhere | 1 hour |
| Contact email | Placeholder | Real email | 5 min |
| GitHub URL consistency | Mixed FOI-Bioinformatics/andreassjodin | One canonical URL | 15 min |

### 2. Security
| Gap | Current State | Required State | Effort |
|-----|---------------|----------------|--------|
| Pickle safety | 4 unsafe fallbacks | No fallbacks; fail if safe_load fails | 2 hours |
| Model format | Pickle (.p) | ONNX or skops (ideal) or verified pickle (minimum) | 1 day |
| Resource limits | No k-mer range cap | Validate k range, warn on large genomes | 2 hours |

### 3. API & Architecture
| Gap | Current State | Required State | Effort |
|-----|---------------|----------------|--------|
| Public API | 91 flat modules, no __all__ | Defined __all__, private prefixes | 1 day |
| Parameter handling | Global mutable state | Immutable config objects (or at minimum, documented) | 2-3 days |
| Version management | 5+ hardcoded copies | Single source of truth | 2 hours |
| Schema versioning | None | Version field in params.json | 4 hours |
| Deprecation policy | Ad-hoc | Documented policy with timeline | 2 hours |

### 4. Scientific Correctness
| Gap | Current State | Required State | Effort |
|-----|---------------|----------------|--------|
| Tm calculation | Duplicate, incorrect in adaptive_filters | Single implementation (thermodynamics.py) | 2 hours |
| Multi-genome coverage | Only first prefix used | All foreground genomes | 1 hour |
| Stale reaction conditions cache | Never invalidated | Reset on parameter change | 1 hour |

### 5. Testing
| Gap | Current State | Required State | Effort |
|-----|---------------|----------------|--------|
| Integration tests | `assert True` stubs | Real end-to-end pipeline test | 1 day |
| Critical module coverage | ~30 modules untested | Tests for pipeline-critical modules | 3-5 days |
| Mechanistic model coverage | 1 of 4 pathways tested | All 4 pathways + interactions | 1 day |
| Property-based tests | None | Hypothesis tests for numerical code | 1-2 days |

### 6. User Experience
| Gap | Current State | Required State | Effort |
|-----|---------------|----------------|--------|
| ParamValidator integration | Optional validate-params command | Automatic validation before pipeline | 2 hours |
| Empty results handling | Silent empty CSV | Clear error with remediation hints | 1 hour |
| Install instructions | Editable install only | `pip install neoswga` primary | 15 min |
| 232 print() statements | Pollutes stdout | Replace with logger calls | 4 hours |

---

## Prioritized Remediation Plan

### Phase 1: Release Blockers (Est. 3-4 days)
1. Fix license mismatch (C1)
2. Remove pickle fallbacks (C2)
3. Fix Tm calculation in adaptive_filters.py (C4)
4. Fix GA ProcessPoolExecutor crash (C5)
5. Replace placeholder emails (M8)
6. Fix MANIFEST.in (H1)
7. Verify model file inclusion with test build (H2)
8. Wire ParamValidator into pipeline commands (H3)
9. Add empty-result detection at end of step2 (H4)
10. Fix division by zero in filter.py (H5)

### Phase 2: API & Architecture (Est. 3-5 days)
11. Define `__all__` and public API surface (C3)
12. Add params.json schema versioning (C7)
13. Single-source version number (H12)
14. Fix version/changelog mismatch (H13)
15. Replace all `print()` with `logger` (19 files)
16. Fix stale reaction_conditions cache (H8)
17. Fix BackgroundAwareOptimizer multi-genome (H7)
18. Document deprecation policy (H17)

### Phase 3: Testing (Est. 5-7 days)
19. Real integration test (C6)
20. Tests for untested pipeline-critical modules (H16)
21. Mechanistic model pathways 2-4 (M31)
22. Network optimizer with real data (M32)
23. Filter tests without excessive mocking (M33)
24. Property-based tests for thermodynamics

### Phase 4: Polish (Est. 2-3 days)
25. Move matplotlib/seaborn/openpyxl to optional deps (M3)
26. Update README for PyPI (M6, M7)
27. Add resource limits for k-mer range (M25)
28. Fix ZipFile handle leak (H6)
29. Move experimental modules out of core/ or prefix with _
30. Migrate model from pickle to safer format (H11)

### Total Estimated Effort: 13-19 days

---

## Appendix: Findings by Audit Agent

| Agent | Critical | High | Medium | Low | Info |
|-------|----------|------|--------|-----|------|
| pkg-auditor (Packaging) | 2 | 2 | 7 | 3 | 5 |
| error-auditor (Error Handling) | 0 | 2 | 7 | 6 | 4 |
| bug-hunter (Code Quality) | 3 | 7 | 8 | 4 | 0 |
| test-auditor (Test Coverage) | 1* | 2* | 3* | 0 | 0 |
| security-auditor (Security) | 0 | 2 | 5 | 5 | 5 |
| api-auditor (API Stability) | 2 | 6 | 7 | 6 | 7 |

*Test findings are cross-referenced with other categories for severity.
