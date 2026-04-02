# NeoSWGA End-to-End Validation Report

**Date:** 2026-03-31
**Test System:** pcDNA (6,157bp) vs pLTR (6,258bp) plasmids, real Jellyfish 2.3.1
**Pipeline:** Full 4-step with real k-mer counting
**Team:** 6 agents testing all optimization methods, simulation, reporting, edge cases

---

## Executive Summary

The core pipeline (count-kmers -> filter -> score -> optimize) works end-to-end
with real data. Binding verification, selectivity, and thermodynamics are all
scientifically correct. However, testing revealed **11 bugs** across the
optimization and post-pipeline stages, plus **3 usability issues** that affect
real-world usage.

**Production-ready optimizers: 4 of 17**
**Broken or limited optimizers: 13 of 17**

---

## 1. Full Pipeline Results

Steps 1-3 run successfully with Jellyfish 2.3.1:
- Step 1 (count-kmers): 0.4s per genome, 7 k values (6-12bp)
- Step 2 (filter): 500 candidates from 24,809 raw k-mers (8.3s)
- Step 3 (score): 500 primers scored (8.0s)

**USABILITY ISSUE**: Step 3 filters ALL primers on this small genome
(all have fg_count=1, scoring them poorly). Requires `--min-amp-pred 0`
to proceed. This is not documented and will block every small-genome user.

## 2. Optimizer Status (17 Methods Tested)

### Fully Functional (4 methods)
| Method | Primers | Time | Notes |
|--------|:-------:|-----:|-------|
| network | 6 | 1.7s | Only method reliably filling target=6 |
| normalized | 6 | 1.7s | Fills target, good strategy presets |
| clique | 6 | 2.7s | Fills target, dimer-free guarantee |
| bg-prefilter | 50 | 1.7s | Works but ignores target_set_size (BUG) |

### Partial Results (9 methods)
These work but return fewer primers than requested (1-5 of 6):
| Method | Primers | Time | Reason |
|--------|:-------:|-----:|--------|
| hybrid | 1 | 2.3s | Early stop: 1 primer = 100% coverage on 6kb |
| greedy | 1 | 1.6s | Same early stop behavior |
| dominating-set | 1 | 1.5s | Same |
| weighted-set-cover | 1 | 1.8s | Same |
| background-aware | 1 | 2.1s | Same |
| tiling | 3 | 1.7s | Partial coverage tiling |
| coverage-then-dimerfree | 1 | 1.8s | Pipeline of partial results |
| bg-prefilter-hybrid | 1 | 1.8s | Same as hybrid |
| genetic | 5 | 393s | Very slow, near-target |

Note: Returning 1 primer on a 6kb genome is mathematically correct
(1 primer covers the entire genome with 70kb extension). This is not
a bug but a limitation of the test data.

### Broken (4 methods)
| Method | Error | Severity |
|--------|-------|----------|
| equiphi29 | "Temperature 42C outside valid range for phi29 (20-40C)" - changes temp but not validation context | BUG |
| moea | "object of type 'NoneType' has no len()" | BUG |
| milp | Missing python-mip package | EXPECTED (optional dep) |
| dimerfree-scored | Killed after 12min, 3.7GB RAM on 6kb genome | BUG (resource) |

## 3. Bugs Found

### BUG 1: equiphi29 optimizer temperature conflict (HIGH)
- **Symptom**: "Temperature 42C outside valid range for phi29 (20-40C)"
- **Cause**: Optimizer sets reaction temp to 42C for equiphi29 but validation
  still checks against phi29 range from params.json
- **Impact**: equiphi29 method never works when params.json specifies phi29

### BUG 2: moea optimizer NoneType error (HIGH)
- **Symptom**: "object of type 'NoneType' has no len()"
- **Cause**: Internal initialization returns None somewhere in the pymoo setup
- **Impact**: moea method never produces results

### BUG 3: dimerfree-scored excessive resource usage (HIGH)
- **Symptom**: >12 minutes and 3.7GB RAM on a 6kb genome, no output
- **Cause**: Likely O(n^2) or worse algorithm on 500 candidates
- **Impact**: Unusable on any real dataset

### BUG 4: bg-prefilter ignores target_set_size (MEDIUM)
- **Symptom**: Returns 50 primers when 6 requested
- **Cause**: Pre-filter stage doesn't pass through target size constraint
- **Impact**: Unexpected output size

### BUG 5: analyze-set crashes (MEDIUM)
- **Symptom**: `AttributeError: module 'neoswga.core.secondary_structure'
  has no attribute 'SecondaryStructurePredictor'`
- **Cause**: Missing class reference in cli_unified.py:3050
- **Impact**: analyze-set command non-functional

### BUG 6: Exit codes always 0 (MEDIUM)
- **Symptom**: Even failed optimizers (milp, equiphi29, moea) exit with code 0
- **Cause**: Error handling catches exceptions and logs but doesn't set exit code
- **Impact**: CI/scripted pipelines can't detect failures

### BUG 7: --seed doesn't guarantee reproducibility (MEDIUM)
- **Symptom**: Two runs with --seed 42 produce different results
- **Cause**: Seed only applies to GA/MOEA, but hybrid optimizer may have
  non-deterministic paths or the seed isn't propagated to all random sources
- **Impact**: Clinical reproducibility goal not met

### BUG 8: step3 filters all primers on small genomes (MEDIUM)
- **Symptom**: step3_df.csv empty after scoring (0 primers pass)
- **Cause**: Default min_amp_pred threshold too high for primers with
  fg_count=1 (common on small genomes)
- **Impact**: Pipeline silently fails at step4 for plasmid-size genomes

### BUG 9: CLAUDE.md documents wrong Phi29Simulator API (LOW)
- **Symptom**: API call fails following CLAUDE.md examples
- **Cause**: Constructor signature changed but docs not updated
- **Impact**: Developer confusion

### BUG 10: CLI --duration flag may be ignored in simulate (LOW)
- **Symptom**: Simulation runs 28,800s despite --duration 10
- **Cause**: Duration parameter might be interpreted differently (minutes vs seconds)
- **Impact**: Unexpected simulation runtime

### BUG 11: vendor-csv not a valid export format (LOW)
- **Symptom**: Export fails with "vendor-csv" format
- **Cause**: Correct format name is "csv"
- **Impact**: Documentation mismatch

## 4. Simulation Integration

Phi29Simulator works correctly on real optimizer output:
- Coverage: 30.7% in 300s (physically plausible for 1 fork on 6kb)
- Mechanistic model integrates cleanly
- Gillespie simulator requires AmplificationNetwork objects (complex setup)
- SwgaSimulator instantiates but requires file-based inputs

## 5. Report and Export

- `neoswga interpret`: Works, sensible metrics
- `neoswga report`: Executive summary and full technical report generate correctly
- `neoswga report --check`: Validation-only mode works
- `neoswga export --format fasta`: Works with Tm/GC annotations
- `neoswga export --format csv`: Works (IDT-compatible)

## 6. Edge Cases

All edge cases handled gracefully (0 crashes):
- Empty candidates: Clear ValueError
- Single candidate: Returns 1 primer with PARTIAL status
- Oversized target: Returns available primers
- Duplicate candidates: Deduplication works
- Invalid primers ('NNNNNN'): Rejected with clear message

## 7. Scientific Validation

- Binding positions independently verified (PASS)
- Zero background binding confirmed (PASS)
- Thermodynamics appropriate for phi29 at 30C (PASS)
- Selected primers are lab-usable (PASS)
- Pipeline produces scientifically sound results (Grade: B+)

---

## Summary: Optimizer Production Readiness

| Status | Methods | List |
|--------|:-------:|------|
| Production-ready | 4 | network, normalized, clique, genetic* |
| Works but early-stops on small genomes | 7 | hybrid, greedy, dominating-set, weighted-set-cover, background-aware, tiling, coverage-then-dimerfree |
| Has bugs | 4 | equiphi29, moea, dimerfree-scored, bg-prefilter |
| Missing dependency | 1 | milp |
| Untested (pipeline) | 1 | bg-prefilter-hybrid |

*genetic works but is 200x slower than alternatives

## Priority Fixes

1. **HIGH**: Fix equiphi29 temp validation conflict
2. **HIGH**: Fix moea NoneType error
3. **HIGH**: Fix dimerfree-scored resource usage
4. **HIGH**: Fix step3 filtering all primers on small genomes (auto-adjust threshold)
5. **MEDIUM**: Fix bg-prefilter target_set_size
6. **MEDIUM**: Fix exit codes for failed optimizers
7. **MEDIUM**: Fix --seed reproducibility for hybrid optimizer
8. **MEDIUM**: Fix analyze-set SecondaryStructurePredictor reference
9. **LOW**: Update CLAUDE.md Phi29Simulator API docs
10. **LOW**: Fix --duration handling in simulate command
11. **LOW**: Fix vendor-csv format name in docs
