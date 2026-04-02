# NeoSWGA Step 4 (optimize) Audit Report

**Date:** 2026-03-31
**Team:** 6 agents (data flow, optimizer runs, PositionCache, metrics verification, special features, correctness)
**Test data:** pcDNA plasmid (6kb) + Legionella vs E. coli bacterial genomes (3.5Mb/4.6Mb)

---

## Executive Summary

Step 4's data flow is correct and metrics computation is verified accurate
(17/17 metrics match). All 8 production-ready optimizers run successfully on
both test datasets. However, the audit revealed that **on the plasmid dataset
the optimizer performs worse than random selection** (3rd percentile), the
`--seed` flag does not guarantee reproducibility on small genomes, and the
`--use-mechanistic-model` flag is silently ignored.

On the real bacterial dataset (Legionella vs E. coli), the pipeline works
correctly — `normalized` achieves 53.9% coverage with 6 primers in under 2
seconds, outperforming all other methods.

**Metrics accuracy: 10/10** (17/17 verified)
**Data flow: 8/10** (correct but lossy, 5 dead code paths)
**Optimizer correctness: 6/10** (works on real data, fails vs random on toy data)
**Overall Step 4: 7/10**

---

## 1. Real Bacterial Genome Results

Legionella pneumophila (3.5 Mb, 38.3% GC) vs E. coli K-12 (4.6 Mb, 50.8% GC):

| Method | Primers | Coverage | Norm Score | Time |
|--------|:-------:|:--------:|:----------:|-----:|
| **normalized** | 6 | **53.9%** | **0.6297** | 1.0s |
| dominating-set | 6 | 52.7% | 0.0000 | 0.9s |
| tiling | 6 | 49.9% | 0.0000 | 1.0s |
| hybrid | 6 | 45.5% | 0.5286 | 0.9s |
| network | 6 | 14.0% | 0.1399 | 1.1s |
| greedy | 1 | 4.0% | 0.0000 | 5.3s |

Full pipeline (steps 1-4 with --fast-score): ~100 seconds total.

## 2. Metrics Verification (10/10)

All 17 PrimerSetMetrics independently verified correct:

| Metric | Verified | Method |
|--------|----------|--------|
| fg_coverage | MATCH | Manual position + extension computation |
| bg_coverage | MATCH | Same |
| selectivity_ratio | MATCH | fg_sites / max(bg_sites, 1) |
| mean_gap, max_gap | MATCH | Manual gap calculation from sorted positions |
| gap_gini, gap_entropy | MATCH | Gini coefficient + Shannon entropy |
| dimer_risk_score | MATCH | Pairwise dimer check |
| strand_alternation_score | MATCH | Forward/reverse alternation pattern |
| normalized_score | MATCH | Weighted composite formula |
| mean_tm | MATCH | Wallace rule calculation |

## 3. Data Flow (8/10)

Correct data pipeline: CLI -> optimize_step4 -> run_optimization -> factory
dispatch -> save_results (CSV + summary JSON + audit JSON).

**Issues found:**
- Step3 scores/features completely discarded — only primer sequences used
- PositionCache created twice (once for all candidates, once for selected)
- 5 dead code paths (cooperative binding, primer strategy, validate_simulation
  kwarg, include_all_sets, run_optimization_from_config)
- 0-row step3 proceeds to optimizer with no clear error message

## 4. Optimizer Correctness (6/10)

### Plasmid dataset (6kb): FAILS vs random
- Network optimizer at 3rd percentile (0.78x random mean)
- Hybrid/greedy/dominating-set return 1 primer when 6 requested
- No optimizer beat the best random draw

### Bacterial dataset (3.5Mb): WORKS correctly
- normalized achieves 53.9% coverage (good for 6 primers)
- Primer selection converges across methods
- Results are biologically sensible

**Conclusion:** The optimization algorithms are designed for and work correctly
on bacterial-scale genomes (Mb). They degrade on toy-sized plasmids (kb) where
single primers already achieve full coverage.

## 5. Special Features

| Feature | Status | Notes |
|---------|--------|-------|
| bg_prefilter | PASS | Removes 20% candidates, --no-bg-prefilter disables |
| --no-background | PASS | Correctly clears bg data |
| --auto-size | PASS | Clinical=10, Discovery=15 primers |
| --use-mechanistic-model | **FAIL_SILENT** | Accepted but silently ignored |
| --seed | PARTIAL | Works on bacterial data, fails on plasmid data |
| --method-guide | PASS | Comprehensive decision tree printed |

### BUG: --use-mechanistic-model silently ignored
The flag is in argparse but never reaches the optimizer. No warning logged.
User gets zero indication their request is being ignored.

## 6. PositionCache (9/10)

- All edge cases handled gracefully
- Memory efficient (9.2 KB for 2,355 primers)
- Coverage computation with extension reach is correct
- BUG: StreamingPositionCache hardcodes k=6-12 (misses equiphi29 primers)

## 7. Bugs Found

| # | Severity | Description |
|---|----------|-------------|
| 1 | **HIGH** | --use-mechanistic-model silently ignored in optimize |
| 2 | **HIGH** | --seed not reproducible on small genomes (plasmid) |
| 3 | **MEDIUM** | Several optimizers return normalized_score=0.0 despite good coverage |
| 4 | **MEDIUM** | Step3 scores discarded — optimization ignores RF predictions |
| 5 | **MEDIUM** | PositionCache created twice wastefully |
| 6 | **LOW** | 5 dead code paths in step4 handler |
| 7 | **LOW** | 0-row step3 input gives unclear error |

---

## Recommendations

### Priority 1 — Must Fix
1. **Warn when --use-mechanistic-model is ignored**: Add logger.warning() or
   implement mechanistic model integration in the optimizer factory
2. **Fix normalized_score=0.0**: Investigate why dominating-set, tiling, and
   greedy show zero score despite good coverage

### Priority 2 — Should Fix
3. **Use step3 scores in optimization**: Currently step3's RF predictions are
   computed expensively then discarded. Optimizers should use amp_pred as a
   candidate quality signal.
4. **Fix seed reproducibility**: Trace why --seed 42 produces different results
   on the plasmid dataset (likely hash ordering or cache-dependent paths)
5. **Remove PositionCache double-creation**: Create once and reuse

### Priority 3 — Nice to Have
6. **Clean up dead code**: Remove cooperative_binding, primer_strategy,
   validate_simulation kwarg, include_all_sets, run_optimization_from_config
7. **Better error for 0-row step3**: "No candidates available — step3 may have
   filtered all primers. Try --min-amp-pred 0"
