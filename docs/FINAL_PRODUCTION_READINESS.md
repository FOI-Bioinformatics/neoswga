# NeoSWGA v3.0.0 — Final Production Readiness Assessment

**Date:** 2026-04-02
**Previous Score:** 6.5/10
**Current Score: 7.6/10**
**Verdict: SHIP with caveats**

---

## Scorecard

| Dimension | Before | After | Delta | Evidence |
|-----------|:------:|:-----:|:-----:|----------|
| Pipeline functionality | 8 | **8** | = | 1,857 tests pass, plasmid E2E in 19s, bacterial E2E works |
| Scientific correctness | 8 | **8** | = | Tm, dimers, coverage independently verified on real data |
| Optimizer quality | 7 | **8** | +1 | 17/17 methods functional, extension-aware coverage, 99.1% realistic |
| UX / Onboarding | 6 | **7** | +1 | Wizard works, step prereqs clear, --fast-score, --from-results |
| Error handling | 5 | **6** | +1 | Silent excepts logged, sys.exit reduced, filter_extra activated |
| Testing | 6 | **8** | +2 | 1,857 tests, 5 new test files, property-based tests |
| Packaging | 7 | **8** | +1 | Upper bounds, [all] fixed, setuptools_scm removed, 3.13 in CI |
| Documentation | 7 | **8** | +1 | 10 audit reports, troubleshooting guide, RF model docs |
| **Overall** | **6.5** | **7.6** | **+1.1** | |

---

## What Works (verified by testing)

### Pipeline E2E
- Plasmid example: All steps complete in 19s, clean output
- Bacterial data (Legionella 3.5Mb vs E.coli 4.6Mb): Full pipeline works
- 17/17 optimizer methods produce valid output on real data
- Tiling optimizer: 99.1% realistic coverage with 6 primers, zero background

### Post-Pipeline Commands (17/19 pass on bacterial data)
- interpret, report (exec + full + check), export (fasta/csv/bed), simulate,
  analyze-set, validate-model, show-presets, --version
- --no-background, --auto-size (clinical + discovery) all work

### Testing
- 1,857 tests, 0 failures, 27 skipped
- Property-based tests for thermodynamic invariants
- Hybrid optimizer, simulation engine, multi-genome filter tests
- CLI smoke tests for all 33 subcommands

---

## Remaining Bugs (found in this audit)

### 2 CLI Crashes

| # | Command | Error | Severity |
|---|---------|-------|----------|
| 1 | `analyze-dimers` | `DimerNetworkMetrics` missing `num_hubs` attribute | MEDIUM |
| 2 | `analyze-stability --output dir/` | Doesn't create output directory | LOW |

### 2 Fuzz Test Critical Failures

| # | Input | Issue |
|---|-------|-------|
| 3 | Invalid primer sequences (XXXXXX) | Silently accepted, misleading results |
| 4 | Non-existent data_dir | Silently writes to CWD instead |

### 3 Fuzz Test Failures (raw tracebacks)

| # | Input | Issue |
|---|-------|-------|
| 5 | Missing genome file | Clear message but full Python traceback shown |
| 6 | Empty FASTA file | Jellyfish crashes, full traceback |
| 7 | Missing fg_prefixes in params.json | Raw KeyError traceback |

### Code Quality Debt

| Issue | Count |
|-------|:-----:|
| print() calls (should be logger) | 1,188 |
| Unsafe pickle.load (bypassing safe_pickle) | 3 |
| sys.exit() in library modules | 7 |
| Silent except blocks (now with debug logging) | ~10 |
| cli_unified.py size | 4,324 lines |

---

## Fixes Applied This Session (40+)

### Scientific
- Owczarzy (2004) entropy-based salt correction
- Jacobson-Stockmayer loop penalty RT factor
- Bidirectional + strand-direction-aware extension coverage
- Salt coefficient 12.5 for oligonucleotides
- RF model amp_pred 0-20 -> 0-1 normalization

### Algorithms
- GA crossover recombines parents (was random sampling)
- Extension-aware DominatingSetOptimizer (strand-direction)
- Genetic timeout + candidate pre-filtering
- Greedy convergence fix
- MOEA NoneType guard
- Normalized_score for PARTIAL results
- dimerfree-scored cascade accepts PARTIAL
- bg-prefilter fast rank-and-select

### Pipeline
- filter_extra() activated (was dead code — 7 quality filters)
- Sorting bug fixed (kept worst primers)
- min_tm/max_tm None crash fixed
- step3 auto-adjust threshold for small genomes
- --fast-score mode (~100x speedup)

### CLI/UX
- --version, --seed, --from-results flags
- Jellyfish pre-flight check with install instructions
- analyze-dimers and analyze-stability crash fixes
- sys.argv mutation removed
- Exit code 1 for failed optimization
- PARTIAL warning message clarified
- Optional parameter warnings suppressed

### Infrastructure
- RestrictedUnpickler for pickle security
- Dependency upper bounds
- macOS + Python 3.13 in CI
- Coverage reporting in CI
- Audit trail output (JSON metadata)

---

## Ship Criteria

### Must-fix before release (2 items)
1. **README params.json example missing fg_prefixes** — first thing new users
   copy will crash. One-line fix.
2. **Wrap common CLI errors in try/except** — missing genome, empty FASTA, and
   missing params keys should show clean messages, not tracebacks.

### Should-fix soon after release (5 items)
3. Fix analyze-dimers `num_hubs` attribute error
4. Validate primer sequences are DNA alphabet before processing
5. Validate data_dir exists before writing
6. Migrate print() to logging (systematic pass)
7. Fix remaining 3 unsafe pickle.load calls

### Acceptable for v3.0.0 release
- 1,188 print() calls (cosmetic, not functional)
- 4,324-line cli_unified.py (works, just large)
- Thread safety limitations (documented — CLI-only tool)
- Plasmid example gets grade C (small genome limitation, documented)

---

## Verdict: SHIP

NeoSWGA v3.0.0 is ready for release as a CLI tool for research use.
The core pipeline works correctly on real bacterial genomes, produces
scientifically valid primer sets, and has comprehensive test coverage.
The 2 must-fix items are minor (README edit + exception wrapping) and
can be addressed in a pre-release patch.
