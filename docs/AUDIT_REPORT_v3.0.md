# NeoSWGA v3.0.0 — Comprehensive Audit Report

**Date:** 2026-03-30
**Audit Team:** 6 specialists + lead coordinator
**Scope:** Scientific correctness, algorithm quality, UX, testing, production readiness, gap analysis

---

## Executive Summary

NeoSWGA is an ambitious and architecturally sound primer design tool with strong thermodynamic foundations, 13+ optimization algorithms, adaptive filtering for extreme-GC genomes, and a genuinely excellent setup wizard. However, the tool has significant gaps that prevent production deployment: ~30% of core modules lack tests, the genetic algorithm crossover operator is broken, simulation engines are unvalidated against experimental data, and dependency management creates installation fragility.

### Overall Production Readiness: 6.0 / 10 (Functional Beta)

| Dimension | Score | Auditor |
|-----------|-------|---------|
| Thermodynamic Models | 7.2 | Scientific Methods Reviewer |
| Filtering & Scoring | 7.0 | Filter/Scoring Specialist |
| Optimization Algorithms | 6.5 | Algorithm Specialist |
| User Experience | 6.2 | UX Specialist |
| Test Coverage & QA | 5.0 | QA Lead |
| Production Engineering | 5.5 | DevOps Specialist |
| **Weighted Average** | **6.0** | |

---

## Team Structure

| Member | Role | Modules Audited | Score Given |
|--------|------|----------------|-------------|
| **Lead** | Audit Coordinator | Cross-cutting synthesis | — |
| **Member 1** | Scientific Methods Reviewer | 10 thermodynamics/chemistry modules | 7.2/10 |
| **Member 2** | Filtering & Scoring Specialist | 12 filter/scorer/report modules | 7.0/10 |
| **Member 3** | Algorithm & Optimization Specialist | 13 optimizers + 3 simulators + 3 performance modules | 6.5/10 |
| **Member 4** | UX & Onboarding Specialist | CLI, wizard, docs, 5 user scenarios | 6.2/10 |
| **Member 5** | QA & Testing Lead | 63 test files, CI/CD pipeline | 5.0/10 |
| **Member 6** | Production Engineering | Dependencies, security, logging, packaging | 5.5/10 |

---

## I. Scientific Correctness (7.2/10)

### Strengths
- **SantaLucia (1998) NN parameters correctly implemented** with all 10 canonical stacks verified
- **Comprehensive literature citations** throughout (Owczarzy 2008, Rees 1993, Blanco 1989, etc.)
- **Adaptive GC filtering** solves a real problem for extreme-GC organisms (Francisella, Burkholderia)
- **Four-pathway mechanistic model** is novel and scientifically ambitious
- **Additive chemistry** well-modeled with Arrhenius kinetics and sigmoid GC-normalization

### Critical Issues

| # | Finding | Severity | Module |
|---|---------|----------|--------|
| S1 | Jacobson-Stockmayer loop penalty missing RT factor — loop penalties ~60% too high for large loops | HIGH | `secondary_structure.py` |
| S2 | Salt correction is additive to Tm rather than entropy-based (Owczarzy 2004) — 1-3C error for short primers | MEDIUM | `thermodynamics.py` |
| S3 | No mismatch thermodynamics in primary Tm pathway — primers bind with mismatches in practice | MEDIUM | `thermodynamics.py` |
| S4 | Deprecated `thermo_estimation.py` still callable with incorrect coefficients | MEDIUM | `thermo_estimation.py` |
| S5 | Model validation tests can never fail (temperature optimum, binding kinetics) | MEDIUM | `model_validation.py` |
| S6 | Activation energies for all additives are "estimated" without uncertainty ranges | LOW | `mechanistic_params.py` |
| S7 | AT-rich template accessibility assumed perfect (GC <= 50% = 1.0) | LOW | `mechanistic_model.py` |

---

## II. Filtering & Scoring Pipeline (7.0/10)

### Strengths
- Five-rule sequence filtering with genome-GC-aware adaptation
- Thermodynamic filter with polymerase-aware criteria and parallel heterodimer checking
- Multi-genome differential filtering (target/background/blacklist with 5x penalty)
- Three-tier stringency presets for background filtering (lenient/moderate/strict)

### Critical Issues

| # | Finding | Severity | Module |
|---|---------|----------|--------|
| F1 | **Monkey-patching** in adaptive_filters.py replaces filter_extra at runtime — fragile, duplicates logic | HIGH | `adaptive_filters.py` |
| F2 | **RF model is a black box** — no training data docs, no validation metrics, no feature importances published | HIGH | `rf_preprocessing.py` |
| F3 | **Coverage estimation fallback** uses `n_primers * 30000` heuristic — wildly inaccurate | HIGH | `report/metrics.py` |
| F4 | **Scoring weight contradiction** — dimer weight is 35% in integrated_quality_scorer but 5% in report/quality.py | MEDIUM | Multiple |
| F5 | **K-mer sampling enabled by default** (10%) — silently reduces scoring accuracy | MEDIUM | `rf_preprocessing.py` |
| F6 | **Short primer Tm bias** — optimal Tm range set to reaction_temp + 15-30C, but 8bp primers have Tm ~20C | MEDIUM | `integrated_quality_scorer.py` |
| F7 | Primers with missing jellyfish counts pass frequency filter | LOW | `filter.py` |

---

## III. Optimization Algorithms (6.5/10)

### 13 Optimizers Assessed

| Optimizer | Correctness | Reproducible | Scalable | Production Ready |
|-----------|:-----------:|:------------:|:--------:|:----------------:|
| Greedy BFS | 8 | 9 | 7 | YES |
| Network | 8 | 8 | 5 | YES |
| Dominating Set | 9 | 9 | 8 | YES |
| Hybrid (default) | 8 | 8 | 7 | YES |
| Background-Aware | 8 | 8 | 6 | YES |
| Genetic Algorithm | 6 | 5 | 4 | NO |
| MILP | 9 | 10 | 3 | CONDITIONAL |
| MOEA (NSGA-III) | 7 | 7 | 3 | NO |
| EquiPhi29 | 7 | 7 | 7 | YES |
| Tiling | 7 | 9 | 8 | YES |
| Clique | 7 | 8 | 4 | NO |
| Normalized | 8 | 9 | 7 | YES |
| Multi-Agent | 6 | 5 | 6 | NO |

### Critical Issues

| # | Finding | Severity | Module |
|---|---------|----------|--------|
| O1 | **GA crossover is broken** — `_crossover()` samples from entire pool, not parent primers. GA is effectively random search with elitism | CRITICAL | `genetic_algorithm.py` |
| O2 | **GPU acceleration is a no-op** — Python for-loop over primers negates CuPy benefit. "10-100x speedup" claim is misleading | HIGH | `gpu_acceleration.py` |
| O3 | **No experimental validation** of any simulation predictions against real SWGA data | HIGH | All simulators |
| O4 | **Network optimizer amplification prediction uncalibrated** — 2^(n/10) heuristic has no experimental basis | MEDIUM | `network_optimizer.py` |
| O5 | Replication simulator binding model is not a proper Poisson process | MEDIUM | `replication_simulator.py` |
| O6 | MILP optimizer not registered in factory (test failure confirms) | MEDIUM | `milp_optimizer.py` |
| O7 | Dual optimizer interface (native vs BaseOptimizer adapters) creates confusion | LOW | Multiple |

### Simulation Engines

| Engine | Physics Score | Validation Status |
|--------|:------------:|:-----------------:|
| Agent-Based Phi29 | 7/10 | No experimental validation |
| Gillespie Stochastic | 8/10 | No experimental validation |
| SWGA Comprehensive | 6/10 | No experimental validation |

---

## IV. User Experience & Ease of Getting Started (6.2/10)

### Scenario Ratings

| Scenario | Rating | Key Blocker |
|----------|:------:|-------------|
| 1. New user, standard SWGA (phi29, bacterial) | 7/10 | Wizard doesn't mention `neoswga design` single-command pipeline |
| 2. GC-rich organism (equiphi29) | 8/10 | `suggest` and `init` commands disconnected |
| 3. Clinical diagnostic lab | 5/10 | **No reproducibility controls, no audit trail** |
| 4. Multi-genome pan-primer design | 6/10 | No wizard integration, CLI-only, no example data |
| 5. Simulation & validation | 5/10 | **No pipeline-to-simulation bridge** — manual primer extraction |

### Critical Issues

| # | Finding | Severity | Impact |
|---|---------|----------|--------|
| U1 | **No --seed flag anywhere in CLI** — clinical users cannot reproduce results | CRITICAL | Scenario 3 |
| U2 | **No pipeline-to-simulation bridge** — must manually extract primers from CSV | HIGH | Scenario 5 |
| U3 | **No audit trail** — no run metadata (version, timestamp, params hash) saved with results | HIGH | Scenario 3 |
| U4 | **No progress indicators** for long-running steps | MEDIUM | All scenarios |
| U5 | Multi-genome command bypasses params.json entirely — different workflow | MEDIUM | Scenario 4 |
| U6 | Error message quality: 7.5/10 — excellent for most paths, generic catch-all at CLI top level | LOW | All |

### Documentation Assessment
- **Exists and good**: README, QUICK_START, user-guide, optimization_guide, multi-genome-guide, SWGA_SCIENCE, FROM_RESULTS_TO_LAB, CHANGELOG
- **Missing**: Troubleshooting guide, performance tuning, clinical/regulatory workflow, interactive tutorial

---

## V. Test Coverage & QA (5.0/10)

### Test Suite: 1,728 passed, 2 failed, 28 skipped

### Coverage Summary

| Category | Count | Notes |
|----------|:-----:|-------|
| Modules with good tests (quality >= 7) | 12 | Thermodynamics (9/10), filter (8/10), position_cache (8/10) |
| Modules with adequate tests (quality 5-6) | 20 | Mix of real and mock-heavy tests |
| Modules with smoke-only tests (quality 4) | 5 | wizard, workflow_selector, condition_suggester, param_validator |
| **Modules with ZERO tests** | **~30** | Including: hybrid_optimizer, replication_simulator, multi_genome_pipeline, network_optimizer (dedicated), rf_preprocessing (dedicated) |

### Critical Untested Paths

1. **Hybrid optimizer** — the DEFAULT optimization method has no dedicated tests
2. **Simulation engine** — all 3 simulators have zero tests
3. **Multi-genome pipeline** — key differentiating feature, zero tests
4. **RF preprocessing** — scoring pipeline feature engineering, zero dedicated tests

### CI/CD Gaps
- No coverage reporting despite pytest-cov being installed
- No Jellyfish in CI — cannot test core pipeline end-to-end
- No type checking (mypy) despite extensive type hints
- No security scanning (bandit, safety)
- `-x` flag stops at first failure, hiding breadth of issues

---

## VI. Production Readiness (5.5/10)

### Checklist (20 items)

| Status | Count | Key Items |
|--------|:-----:|-----------|
| **Pass** | 6 | Exception hierarchy, logging framework, entry point, no shell=True, Jellyfish check, param validation |
| **Partial** | 8 | print/logger ratio (1,185 vs 1,629), bare except blocks (21), memory limits, timeouts, GPU integration, config schema, path traversal, CI pipeline |
| **Fail** | 6 | Dependency version bounds, pickle security, Dockerfile, PyPI publish, conda recipe, file size limits |

### Critical Production Blockers

| # | Finding | Severity |
|---|---------|----------|
| P1 | **No dependency upper bounds** — `scikit-learn>=0.21.3` permits breaking v2.0 | CRITICAL |
| P2 | **Pickle deserialization** — 4 `pickle.load()` calls are arbitrary code execution vectors | HIGH |
| P3 | **1,185 print() calls** mixed with 1,629 logger calls — cannot control output | MEDIUM |
| P4 | **21 bare `except Exception:` blocks** swallow errors silently | MEDIUM |
| P5 | No Dockerfile for reproducible deployment | MEDIUM |
| P6 | No PyPI publish workflow or conda recipe | MEDIUM |

---

## VII. Consolidated Gap Analysis

### Tier 1 — Must Fix (blocks production use)

| # | Gap | Source | Effort |
|---|-----|--------|--------|
| G1 | GA crossover produces random individuals, not recombinations | Optimizer audit | Low |
| G2 | No --seed CLI flag for reproducible optimization | UX audit | Low |
| G3 | Dependency version bounds (sklearn, numpy, pandas) | DevOps audit | Low |
| G4 | Pickle deserialization security (4 calls) | DevOps audit | Medium |
| G5 | Hybrid optimizer (default method) has zero tests | QA audit | Medium |
| G6 | Jacobson-Stockmayer loop penalty missing RT factor | Science audit | Low |
| G7 | MILP optimizer not registered in factory | QA audit | Low |

### Tier 2 — Should Fix (significant quality gaps)

| # | Gap | Source | Effort |
|---|-----|--------|--------|
| G8 | Simulation engines have zero tests and no experimental validation | QA + Optimizer | High |
| G9 | RF model is undocumented black box (no training data, no metrics) | Filter audit | Medium |
| G10 | Coverage estimation fallback is misleading (n_primers * 30000) | Filter audit | Low |
| G11 | Scoring weight contradiction (dimer 35% vs 5%) | Filter audit | Low |
| G12 | Monkey-patching in adaptive_filters.py | Filter audit | Medium |
| G13 | No pipeline-to-simulation bridge (manual primer extraction) | UX audit | Medium |
| G14 | No audit trail / run metadata saved with results | UX audit | Medium |
| G15 | GPU acceleration claims are misleading (Python for-loop) | Optimizer audit | Low |
| G16 | No coverage reporting in CI | QA audit | Low |
| G17 | Multi-genome pipeline has zero tests | QA audit | Medium |

### Tier 3 — Nice to Have (polish and completeness)

| # | Gap | Source | Effort |
|---|-----|--------|--------|
| G18 | Salt correction should use entropy-based method (Owczarzy 2004) | Science audit | Medium |
| G19 | No mismatch thermodynamics in primary Tm pathway | Science audit | Medium |
| G20 | No property-based testing (hypothesis) | QA audit | Medium |
| G21 | Dockerfile for reproducible deployment | DevOps audit | Low |
| G22 | Troubleshooting documentation | UX audit | Low |
| G23 | Clinical/regulatory workflow guide | UX audit | Medium |
| G24 | No progress indicators for long-running steps | UX audit | Medium |
| G25 | Deprecated thermo_estimation.py still callable | Science audit | Low |
| G26 | Short primer Tm bias in quality scorer | Filter audit | Low |

---

## VIII. What NeoSWGA Does Well

Despite the gaps, the audit team identified substantial strengths worth preserving:

1. **Thermodynamic foundation is solid** — SantaLucia NN implementation with proper caching (1M entries) is production-quality
2. **Adaptive GC filtering is innovative** — Solves a real problem that other tools don't address
3. **Setup wizard is genuinely excellent** — Genome-adaptive, 5 GC classes, smart polymerase recommendations
4. **Optimizer breadth is impressive** — 13+ algorithms covering different tradeoffs with a clean factory/registry pattern
5. **BaseOptimizer abstract interface** — Well-designed with immutable results and normalized scoring
6. **Exception hierarchy** — 20+ specific exception types with structured fields
7. **Multi-genome differential filtering** — Target/background/blacklist with weighted penalties
8. **Mechanistic four-pathway model** — Novel approach for primer design tools
9. **Position cache** — O(1) lookup with 1000x speedup over HDF5 reads
10. **Test suite breadth** — 1,728 tests with excellent thermodynamics coverage (9/10)

---

## IX. Recommended Priority Actions

### Sprint 1 (1 week) — Critical fixes
1. Fix GA crossover to recombine parent primers, not sample from pool
2. Add `--seed` flag to CLI for all stochastic optimizers
3. Add dependency upper bounds to pyproject.toml
4. Fix Jacobson-Stockmayer loop penalty (add RT factor)
5. Fix MILP optimizer factory registration
6. Remove `-x` from CI, add `--cov` with threshold

### Sprint 2 (2 weeks) — Testing gaps
7. Write dedicated tests for hybrid_optimizer (default method)
8. Write tests for replication_simulator and stochastic_simulator
9. Write tests for multi_genome_pipeline
10. Add coverage reporting to CI
11. Fix scoring weight contradiction (standardize dimer weight)

### Sprint 3 (3 weeks) — Production hardening
12. Replace pickle.load() with safer alternatives (skops.io or RestrictedUnpickler)
13. Document RF model (training data, metrics, feature importances)
14. Unify adaptive filtering (eliminate monkey-patching)
15. Add --from-results flag to simulate command
16. Add audit trail (version, timestamp, params hash) to pipeline output
17. Create Dockerfile

### Sprint 4 (4 weeks) — Polish
18. Add property-based tests for thermodynamic invariants
19. Fix GPU acceleration or remove misleading claims
20. Add troubleshooting documentation
21. Add progress indicators for long-running steps
22. Implement entropy-based salt correction (Owczarzy 2004)
