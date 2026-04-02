# NeoSWGA Step 2 (filter) Audit Report

**Date:** 2026-03-31
**Team:** 6 agents (filter rules, output verification, data flow, parameter testing, adaptive filtering, thermodynamic filtering)

---

## Executive Summary

Step 2 has a **critical architectural bug**: the main quality filtering function
`filter_extra()` — containing 7 filters for Tm, homopolymer runs, GC clamp,
dinucleotide repeats, self-dimers, GC content, and 3'-end GC — is **dead code
that is never called** by any pipeline path. This means primers with homopolymer
runs (AAAAAA), self-dimers, and bad GC clamps pass through to optimization
unchecked. Only loose pre-filters (wide GC range, wide Tm margin) and the
frequency/Gini filters are active.

Additionally, a sorting bug keeps the WORST primers (highest bg/fg ratio) when
capping at max_primer, and min_tm/max_tm defaults crash step2 when not in
params.json.

**Core data accuracy: 10/10** (all column values independently verified)
**Filtering completeness: 3/10** (7 of 10 intended filters are dead code)
**Overall Step 2: 5/10**

---

## 1. Data Accuracy (10/10)

All step2_df.csv column values independently verified for 20 random primers:

| Column | Matches | Method |
|--------|:-------:|--------|
| fg_count | 20/20 | Cross-checked against Jellyfish dump files |
| bg_count | 20/20 | Cross-checked against Jellyfish dump files |
| ratio | 20/20 | Confirmed ratio = bg_count / fg_count |
| gini | 20/20 | Recomputed from HDF5 position files |

No NaN, Inf, or invalid sequences found.

## 2. CRITICAL: filter_extra() is Dead Code

The function `filter_extra()` at `filter.py:216` defines 7 quality filters
but is **never called** by any pipeline code path:

| Filter | What It Checks | Would Reject | Status |
|--------|---------------|:------------:|--------|
| Effective Tm | Tm with additive corrections | 1,862 | DEAD CODE |
| Homopolymer | Runs of 5+ identical bases | 155 | DEAD CODE |
| GC content | Strict 37.5-62.5% range | 2,578 | DEAD CODE |
| GC clamp | 1-3 G/C in last 5 bases | 888 | DEAD CODE |
| 3' end GC | Max 2/3 G/C at 3' end | 299 | DEAD CODE |
| Dinucleotide | Repeat patterns (ATATAT...) | 0 | DEAD CODE |
| Self-dimer | >4bp self-complementarity | 415 | DEAD CODE |

**If applied, ~6,197 additional primers would be removed** from the 11,454
post-frequency pool, yielding higher-quality candidates for optimization.

### What IS Active

Only these filters actually execute:

| Stage | Filter | Removes | Active In |
|-------|--------|:-------:|-----------|
| Pre-filter | GC range (28-72%) | 7,712 | kmer_counter.py |
| Pre-filter | Tm range (0-60C, wide margin) | 3,214 | kmer_counter.py |
| Pipeline | Foreground frequency (min_fg_freq) | varies | pipeline.py |
| Pipeline | Background frequency (max_bg_freq) | 12,599 | pipeline.py |
| Pipeline | Gini index (max_gini) | varies | pipeline.py |
| Pipeline | max_primer cap | remainder | pipeline.py |

## 3. Sorting Bug: Keeps Worst Primers

`pipeline.py:691-693` — ratio = bg_count/fg_count sorted DESCENDING then
capped at max_primer. This keeps primers with the HIGHEST background-to-foreground
ratio (the worst candidates). Should sort ASCENDING to keep the best.

Impact: When more primers survive than max_primer (e.g., 10,000 survive but
max_primer=500), the 500 with the WORST selectivity are kept.

## 4. min_tm/max_tm None Crash

`kmer_counter.py:529` — `min_tm - wide_tm_margin` crashes with TypeError when
min_tm/max_tm are not in params.json (default None). Any params.json without
explicit Tm values crashes step2. 3 of 5 parameter tests failed from this.

## 5. Tm Method Disagreement

Two different Tm calculation methods are used:
- **Active** (pre-filter): `melting_temp.temp()` — vendored melt package with
  intentional GC bug (Na=10mM, Mg=20mM defaults)
- **Dead code** (filter_extra): `calculate_effective_tm()` — Owczarzy entropy
  method (Na=50mM, Mg=0mM defaults)

**26C average Tm disagreement** between the two methods for 6bp primers.
User salt settings don't affect the active pre-filter.

## 6. Adaptive Filtering: Working Correctly

The monkey-patching refactor is verified clean:
- parameter.gc_min/gc_max correctly read by filter.py
- Original values restored in finally block
- No remaining monkey-patching code
- Minor: AdaptiveGCFilter caps at 0.20-0.80, parameter.py caps at 0.15-0.85

## 7. Parameter Configuration Tests

| Test | Result | Notes |
|------|--------|-------|
| Very strict | PASS | 0 primers (expected) |
| Very relaxed | PASS | 10,000 primers. GC pre-filter still active (removed 5,719) |
| No background | CRASH | TypeError: None - float |
| Single k=8 | CRASH | TypeError: None - float |
| equiphi29 | CRASH | TypeError: None - float |

Tests 3-5 all crash from the min_tm/max_tm None bug.

---

## All Bugs (Ranked)

| # | Severity | Description | Location |
|---|----------|-------------|----------|
| 1 | **CRITICAL** | filter_extra() is dead code — 7 quality filters never run | filter.py:216 (never called) |
| 2 | **HIGH** | Sorting bug keeps worst primers (highest bg/fg ratio) | pipeline.py:691 |
| 3 | **HIGH** | min_tm/max_tm=None crashes step2 | kmer_counter.py:529 |
| 4 | **HIGH** | 26C Tm disagreement between active and dead-code methods | melting_temp vs thermodynamics |
| 5 | **MEDIUM** | Inconsistent salt defaults (Na=10/Mg=20 vs Na=50/Mg=0) | melting_temp vs thermodynamics |
| 6 | **MEDIUM** | SettingWithCopyWarning in 4 locations | pipeline.py, filter.py |
| 7 | **MEDIUM** | GC pre-filter always active even with relaxed params | kmer_counter.py |
| 8 | **MEDIUM** | Adaptive GC bounds inconsistency (0.20-0.80 vs 0.15-0.85) | adaptive_filters vs parameter |
| 9 | **LOW** | No Tm column in step2_df.csv output | pipeline.py |
| 10 | **LOW** | Duplicated GC calculation in adaptive_filters and filter.py | Multiple |

---

## Recommendations

### Priority 1 — Must Fix
1. **Activate filter_extra()**: Call it from pipeline.py step2() after frequency
   filtering and before position file creation. This is a one-line fix with
   major quality impact.
2. **Fix sorting bug**: Change `ascending=False` to `ascending=True` at
   pipeline.py:691 so the best primers (lowest bg/fg ratio) are kept.
3. **Fix min_tm/max_tm None crash**: Add defaults in kmer_counter.py:
   `min_tm = min_tm if min_tm is not None else 0`

### Priority 2 — Should Fix
4. **Unify Tm methods**: Use calculate_tm_with_salt() (Owczarzy) everywhere
   instead of the vendored melting_temp.temp() with its intentional GC bug.
5. **Fix SettingWithCopyWarning**: Add .copy() after DataFrame slicing.
6. **Add Tm column to step2_df.csv**: Users should see primer Tm values.

### Priority 3 — Nice to Have
7. Unify adaptive GC bounds (0.15-0.85 everywhere)
8. Make GC pre-filter margin configurable
9. Deprecate run_step2_with_adaptive_gc() (parameter.py handles it natively)
