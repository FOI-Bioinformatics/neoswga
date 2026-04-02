# NeoSWGA v3.0.0 — Production Readiness Report

**Date:** 2026-04-02
**Team:** 6 specialists (stress test, code quality, packaging, concurrency, science validation, UX gap analysis)

---

## Overall Production Readiness: 6.5 / 10

| Dimension | Score | Assessment |
|-----------|:-----:|-----------|
| Core pipeline functionality | 8/10 | Steps 1-4 work E2E on real bacterial data |
| Scientific quality | 8/10 | Primers pass all thermodynamic checks, lab-usable |
| Packaging & dependencies | 7/10 | Build works, deps bounded, version consistent |
| User experience | 6/10 | Wizard good, but minimal params.json crashes |
| Error handling & logging | 5/10 | 1,188 print(), 13 silent except, 3 analyze-* crashes |
| Thread safety | 4/10 | Global state makes library use unsafe |

---

## Bugs Found (ranked by severity)

### 3 CLI Crashes (HIGH)

| # | Command | Error | Location |
|---|---------|-------|----------|
| 1 | `analyze-set` | `'dict' object has no attribute 'energy'` — hairpin returns dict, code expects object | `cli_unified.py:3071` |
| 2 | `analyze-dimers` | `module has no attribute 'analyze_dimer_network'` — function doesn't exist | `cli_unified.py:3122` |
| 3 | `analyze-stability` | `module has no attribute 'analyze_primers'` — function doesn't exist | `cli_unified.py:3148` |

### 3 Critical UX Blockers

| # | Issue | Impact |
|---|-------|--------|
| 4 | Missing `fg_prefixes` in params.json crashes with raw `KeyError` traceback | Blocks every new user who doesn't use the wizard |
| 5 | README example params.json omits `fg_prefixes`/`bg_prefixes` | Users who copy the example crash immediately |
| 6 | `validate-params` doesn't flag missing `fg_prefixes` as ERROR | Silent until crash at step1 |

### Concurrency & State Issues

| # | Issue | Severity |
|---|-------|----------|
| 7 | `parameter.py` ~40 mutable globals — two pipeline runs in one process corrupt each other | HIGH for library use |
| 8 | `adaptive_filters.py` mutates `sys.argv` during execution | HIGH (process-wide state) |
| 9 | `thermodynamics.py` LRU cache keys ignore reaction conditions — returns wrong Tm if conditions change | MEDIUM |
| 10 | No file locking on CSV/HDF5 writes — concurrent runs to same dir corrupt files | MEDIUM |

### Code Quality Issues

| # | Issue | Count |
|---|-------|:-----:|
| 11 | `print()` calls that should be `logger.*` | ~968 |
| 12 | Silent `except Exception: pass` blocks | 13 |
| 13 | `sys.exit()` in library modules (should raise exceptions) | 6 |
| 14 | Spurious "PARTIAL" warning when primer count matches target | 1 |

### Packaging Issues

| # | Issue | Severity |
|---|-------|----------|
| 15 | `[all]` extra includes dev deps (pytest, black) and omits plotly | MEDIUM |
| 16 | Python 3.13 in classifiers but not in CI matrix | LOW |
| 17 | `setuptools_scm` in build-requires but unused | LOW |

### Minor Issues

| # | Issue |
|---|-------|
| 18 | interpret/export show phi29/30C instead of equiphi29/42C |
| 19 | auto-size uses hardcoded 1M genome length |
| 20 | analyze-set/dimers parse multiple primers as single string |
| 21 | 11 missing-parameter warnings for optional params with known defaults |

---

## What Works Well

1. **Core pipeline is functional** — steps 1-4 complete in ~80s on 3.5Mb genome
2. **Scientific output is valid** — primers pass all thermodynamic checks
3. **17/17 optimizers produce output** — tiling achieves 99.1% realistic coverage
4. **Zero background binding** on Legionella vs E.coli
5. **Wizard produces valid params.json** — best onboarding path
6. **Reports generate correctly** — executive summary and full technical HTML
7. **All export formats work** — FASTA, CSV, BED, BedGraph
8. **Simulation runs** from optimizer output via --from-results
9. **Packaging is sound** — build succeeds, twine passes, version consistent
10. **1,857 tests pass** with good thermodynamics coverage

---

## Risk Assessment

### Safe for: CLI single-user research use
The tool works correctly when run as a CLI by a single user on one genome
at a time. This is the intended use case and it works well.

### Unsafe for: Library import, web service, parallel pipelines
The global mutable state in parameter.py makes NeoSWGA unsuitable for:
- Importing as a Python library in multi-run workflows
- Web service backends (concurrent requests)
- Running multiple pipelines in the same process
- Any environment where reaction conditions might change between calls

### Risk for: New users who don't use the wizard
Users who create params.json manually (following README examples) will crash
on the first command due to missing `fg_prefixes`. The wizard path works but
is not the default suggestion in most documentation.

---

## Priority Fixes for Production Release

### P0 — Must fix (blocks users)
1. **Auto-derive fg_prefixes from fg_genomes** — eliminates #4, #5, #6
2. **Fix 3 analyze-* crashes** — missing functions and wrong return types
3. **Remove sys.argv mutation** in adaptive_filters.py

### P1 — Should fix (quality)
4. **Add reaction conditions to thermodynamics LRU cache keys**
5. **Fix [all] extra** — remove dev deps, add plotly
6. **Fix spurious PARTIAL warning** when count matches target
7. **Fix auto-size hardcoded genome length**
8. **Suppress optional parameter warnings**

### P2 — Should fix (production hardening)
9. **Replace ~968 print() with logger** (systematic pass)
10. **Fix 13 silent except blocks** — at minimum add logging
11. **Move 6 sys.exit() from library to CLI layer**
12. **Add Python 3.13 to CI matrix**
13. **Add file locking for concurrent safety**

### P3 — Nice to have
14. Fix interpret/export polymerase display
15. Fix analyze-set/dimers argument parsing
16. Remove unused setuptools_scm build-require
17. Add progress bars for long-running steps
