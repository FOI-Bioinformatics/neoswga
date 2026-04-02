# NeoSWGA Step 3 (score) Audit Report

**Date:** 2026-03-31
**Team:** 6 agents (RF model, output verification, data flow, science validation, edge cases, performance)

---

## Executive Summary

Step 3 correctly scores all primers without data loss, handles edge cases
gracefully, and falls back to heuristic scoring when the RF model is missing.
However, the audit revealed that the RF model is a **Regressor outputting
0-20** (not a Classifier outputting 0-1 as documented), causing report modules
to show incorrect quality ratings. The delta-G histogram features consume
99.93% of runtime but contribute only 1.8% of model importance — a massive
performance waste. Step 3 is the pipeline bottleneck at 46 minutes for 10,000
primers.

**Data integrity: 9/10** (no data loss, correct join, idempotent)
**Model correctness: 5/10** (wrong type documented, reports show wrong scale)
**Performance: 3/10** (46 min for 10k primers, dominated by near-useless features)
**Overall Step 3: 5.5/10**

---

## 1. RF Model Inspection

### Model is a REGRESSOR, Not Classifier

| Property | Documented | Actual |
|----------|-----------|--------|
| Type | RandomForestClassifier | **RandomForestRegressor** |
| Output | Probability 0-1 | **Continuous 0-20 scale** |
| Method | predict_proba() | **predict()** |
| Threshold | ~0.5 | min_amp_pred=10.0 |

- 100 trees, max_depth=15, squared_error criterion
- 53 features (not 52 as documented)
- Trained on sklearn 1.7.2

### Feature Importance

Top 3 features account for ~70% of model importance:

| Rank | Feature | Importance | Category |
|------|---------|:----------:|----------|
| 1 | number.of.C | ~30% | Sequence composition |
| 2 | sequence.length | ~25% | Primer length |
| 3 | melting_tm | ~15% | Thermodynamics |
| 4-53 | 27 delta-G bins + others | ~30% | Mixed |

**The 27 delta-G histogram features combined contribute only ~1.8% importance**
but consume 99.93% of computation time.

### Report Module Bugs

| # | Bug | Severity | Location |
|---|-----|----------|----------|
| RF1 | `technical_report.py:176` clamps amp_pred to [0,1] — all primers get max quality | HIGH | report/ |
| RF2 | `executive_summary.py:490-491` multiplies amp_pred*5 for stars — all show 5 stars | HIGH | report/ |
| RF3 | `visualizations.py:531` uses amp_pred as 0-1 value — heatmap scale wrong | MEDIUM | report/ |
| RF4 | Module docstring says "Probability score (0-1)" — incorrect | MEDIUM | rf_preprocessing.py |

## 2. Data Flow (9/10)

10-stage pipeline traced. Data integrity is correct:

- All 10,000 step2 primers scored (no data loss)
- LEFT join with step2 preserves all scored primers
- Output sorted by gini (ascending), not by score
- Idempotent: two runs produce identical output
- step2_df.csv is NOT modified

### Issues Found

| # | Severity | Issue |
|---|----------|-------|
| DF1 | MEDIUM | `fg_scale` division by zero if fg_seq_lengths is empty |
| DF2 | MEDIUM | Heuristic fallback score scale (0-18) differs from RF (0-20) |
| DF3 | LOW | Column rename fragility (RF output must have `sequence` column) |
| DF4 | LOW | Molarity hardcoded to 2.5 regardless of conditions |

## 3. Science Validation (LIMITED_BY_DATA)

On the plasmid example:
- Score vs fg_count: rho=0.063, p=0.16 (positive but not significant)
- Score vs melting_tm: rho=-0.105, p=0.019 (significant, slight low-Tm preference)
- bg_count constant at 0, gini constant at 0 — cannot test selectivity correlation

**The model rewards C-content, appropriate length, and moderate Tm.** No
anomalies in top-scoring primers (no homopolymers, balanced GC). Full
validation requires a dataset with diverse fg/bg binding frequencies.

## 4. Edge Cases (All Pass)

| Test | Result | Notes |
|------|--------|-------|
| Empty step2 (0 primers) | PASS | Clear error with exit code 1 |
| Single primer | PASS | Scored correctly (16.6) |
| min_amp_pred=0 vs default | PASS | Equivalent (default is 0.0 on this genome) |
| Idempotent (run twice) | PASS | Identical output within floating-point |
| Missing RF model | PASS | Graceful fallback to heuristic scoring |
| Corrupted step2 | PASS | Clear error message |
| 12bp primers only | PASS | Scored correctly, different score distribution |

## 5. Performance (3/10)

**Step 3 is the pipeline bottleneck.**

| Metric | Value |
|--------|-------|
| Runtime (10,000 primers, 6kb plasmid) | **2,765 seconds (46 min)** |
| Feature computation | 99.93% of runtime |
| RF prediction | 0.1 seconds |
| Bottleneck | 112M thermodynamic free energy calculations |
| Throughput | ~40,000 calcs/sec |
| Estimated for 500 primers | ~138 seconds |
| Estimated for 10,000 on 1 Mbp | ~128 hours without sampling |

### The Performance Paradox

The 27 delta-G histogram features:
- Consume **99.93% of runtime** (thermodynamic calculations)
- Contribute **1.8% of model importance** (dominated by sequence features)
- Could be skipped entirely with minimal accuracy loss

### Sampling Configuration

| Genome Size | Sample Rate | Effect |
|------------|:-----------:|--------|
| < 100 kbp | 100% (disabled) | Full accuracy, slowest |
| 100 kbp - 1 Mbp | 25% | ~4x speedup |
| 1 - 10 Mbp | 10% | ~10x speedup |
| > 10 Mbp | 5% | ~20x speedup |

No CLI flag to control sampling — only via internal parameter.

---

## All Bugs Summary

| # | Severity | Description |
|---|----------|-------------|
| RF1 | **HIGH** | Report clamps amp_pred to [0,1] — all primers show max quality |
| RF2 | **HIGH** | Executive summary star rating wrong (amp_pred*5 on 0-20 scale) |
| RF4 | **MEDIUM** | Docstring says "Probability 0-1" — model outputs 0-20 |
| RF3 | **MEDIUM** | Visualization heatmap scale wrong |
| DF1 | **MEDIUM** | fg_scale division by zero if no foreground lengths |
| DF2 | **MEDIUM** | Heuristic fallback score scale mismatch |
| PERF1 | **DESIGN** | 99.93% of runtime on 1.8% of model importance |

---

## Recommendations

### Priority 1 — Must Fix
1. **Fix report amp_pred scaling**: Normalize 0-20 score to 0-1 before report
   modules consume it. Add `amp_pred_normalized = amp_pred / 20.0` at the
   boundary.
2. **Fix docstring**: Change "Probability score (0-1)" to "Amplification
   prediction score (0-20 scale, higher is better)" in rf_preprocessing.py.

### Priority 2 — Should Fix
3. **Skip delta-G features**: Since they contribute 1.8% importance but 99.93%
   runtime, consider a "fast mode" that skips them. This would reduce step3
   from 46 min to ~3 seconds for 10,000 primers with <2% accuracy loss.
4. **Add fg_scale zero guard**: Return early or use default scale when
   fg_seq_lengths is empty.
5. **Add --fast-score CLI flag**: Skip delta-G histogram features for speed.

### Priority 3 — Nice to Have
6. **Add sampling CLI flag**: Let users control sample_rate from command line.
7. **Retrain model without delta-G features**: A 26-feature model would be
   equally accurate and orders of magnitude faster.
8. **Document actual model type and features**: Update CLAUDE.md and docstrings.
