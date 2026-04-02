# NeoSWGA Step 1 (count-kmers) Audit Report

**Date:** 2026-03-31
**Team:** 6 agents (Jellyfish integration, count verification, position files, edge cases, GC-adaptive strategy, code path tracing)

---

## Executive Summary

Step 1 core functionality is scientifically correct: k-mer counts match
100% between Jellyfish and manual Python counting across k=6, 8, 12.
Position files are accurate (verified for 20 primers). Edge cases are
handled gracefully. However, the audit uncovered **4 HIGH-severity bugs**
in the integration layer (CLI blacklist handling, QA flag crash, params
KeyError, GC-adaptive override conflict).

**Core k-mer counting: 9/10**
**Integration and CLI layer: 6/10**
**Overall Step 1: 7.5/10**

---

## 1. K-mer Counting Correctness (9/10)

**VERIFIED: 100% accurate** across all tested k-values.

| k | Unique k-mers | Matches | Mismatches |
|---|:-------------:|:-------:|:----------:|
| 6 | 1,891 | 1,891 | 0 |
| 8 | 5,324 | 5,324 | 0 |
| 12 | 6,023 | 6,023 | 0 |

- Canonical counting (Jellyfish -C): correctly combines k-mer and RC
- NeoSWGA canonicalizes via `min(kmer, rc)` before lookup: correct
- Circular genome: Jellyfish counts linearly (0.11% undercount for
  k-1 wrap-around positions). Position search IS circular-aware.
- Ambiguous bases: Jellyfish correctly skips N-containing k-mers

## 2. Jellyfish Integration (7.8/10)

**Strengths:**
- All subprocess calls use list-form arguments (no shell=True): 9/10
- Version checking rejects v1.x with clear error + install URL
- Output parsing handles blank/malformed lines gracefully

**Bugs found:**

| # | Severity | Issue |
|---|----------|-------|
| J1 | MEDIUM | `MultiGenomeKmerCounter.count_kmers_jellyfish()` has NO subprocess timeout — could hang indefinitely |
| J2 | MEDIUM | Missing `text=True` in class subprocess calls (binary/text mismatch) |
| J3 | LOW | Raw `CalledProcessError` propagates without genome/k-value context |
| J4 | LOW | Progress display out-of-order (`as_completed` vs sequential) |

Code duplication: `MultiGenomeKmerCounter` duplicates logic from
`_run_jellyfish_for_k()` with weaker error handling. Should delegate.

## 3. Position Files (9/10)

**VERIFIED: 100% accurate** for all 20 tested primers (both strands).

- Primary search: Aho-Corasick (5-7x speedup) with naive fallback
- Circular genome handling: correct in both code paths
- HDF5 schema: `{prefix}_{k}mer_positions.h5`, keyed by primer sequence

**Issues:**
- Empty datasets stored as float64, populated as int64 (no functional impact)
- No HDF5 metadata (genome length, source, k-mer size)
- StreamingPositionCache hardcodes k=6-12 (misses equiphi29 12-18bp primers)

## 4. Edge Cases (9/10)

All 6 edge cases passed:

| Genome | Size | Result |
|--------|------|--------|
| Tiny (100bp) | 100bp | Pass |
| Multi-contig (3x500bp) | 1500bp | Pass |
| N bases | 220bp | Pass (N k-mers skipped) |
| Homopolymer (all-A) | 200bp | Pass (1 unique k-mer) |
| AT-rich (90%) | 1000bp | Pass |
| GC-rich (90%) | 1000bp | Pass |

**Bug found:**
- **E1 (HIGH)**: `pipeline.py:372` crashes with `KeyError` on `params["bg_prefixes"]`
  when params.json lacks background genome fields. Should use `.get()`.

## 5. GC-Adaptive Strategy (7/10)

**3 classes** (not 5 as documented): AT_RICH (<0.35), BALANCED (0.35-0.65),
GC_RICH (>0.65). Sub-thresholds at 0.25 and 0.70 act as implicit tiers.

**Bugs found:**

| # | Severity | Issue |
|---|----------|-------|
| GC1 | HIGH | Cannot distinguish user-set phi29 from default — always overrides to equiphi29 for balanced/GC-rich genomes |
| GC2 | HIGH | Step4 reads polymerase from raw JSON, not adapted values — pipeline steps may use different polymerases |
| GC3 | MEDIUM | Sharp discontinuity at GC=0.35 boundary (all params jump simultaneously) |
| GC4 | MEDIUM | Adaptive strategy activates silently without user awareness |
| GC5 | LOW | Confidence scores are cosmetic (hardcoded, not data-derived) |

K-mer range priority is correctly implemented: CLI > params.json > GC-adaptive > defaults.

## 6. Code Path Analysis (6/10)

**14 issues found across 6 code paths:**

### Critical Bugs

| # | Severity | Issue |
|---|----------|-------|
| C1 | **HIGH** | **CLI --blacklist genomes never get k-mer counted.** pipeline.step1() completes BEFORE CLI sets parameter.bl_genomes. Blacklist filtering silently skipped in step2. |
| F1 | **HIGH** | **--enable-qa crashes step1.** Calls `pipeline_qa_integration.run_step1_with_qa()` which does not exist (QA module only covers step2-4). |

### Medium Issues
- B1: Exclusion genomes counted twice if in both JSON and CLI
- C2: CLI blacklist path never computes bl_seq_lengths
- D1: Genome library lookup uses inconsistent keys
- X1: min_k/max_k applied 3 times with GC-adaptive override between
- X6: get_value_or_default returns None for missing required params

### Notable Behaviors
- Genome library can skip Jellyfish entirely via symlinks
- GPU flag accepted but has zero effect on step1
- Foreground genomes NOT resolved from genome library (only bg/bl)
- GC auto-calculated by loading entire genome into memory

---

## All Bugs Summary (Ranked)

| # | Severity | Category | Description |
|---|----------|----------|-------------|
| C1 | **HIGH** | CLI | --blacklist genomes never get k-mer counted |
| F1 | **HIGH** | CLI | --enable-qa crashes step1 (function doesn't exist) |
| GC1 | **HIGH** | Strategy | Can't distinguish user-set phi29 from default |
| GC2 | **HIGH** | Strategy | Step4 reads raw JSON polymerase, not adapted value |
| E1 | **HIGH** | Pipeline | KeyError on missing bg_prefixes in params.json |
| J1 | MEDIUM | Jellyfish | MultiGenomeKmerCounter has no subprocess timeout |
| J2 | MEDIUM | Jellyfish | Missing text=True in class subprocess calls |
| GC3 | MEDIUM | Strategy | Sharp boundary discontinuity at GC=0.35 |
| GC4 | MEDIUM | Strategy | Silent adaptive activation |
| B1 | MEDIUM | CLI | Exclusion genomes may be counted twice |
| C2 | MEDIUM | CLI | bl_seq_lengths never computed for CLI blacklist |
| D1 | MEDIUM | Library | Inconsistent genome library lookup keys |
| X1 | MEDIUM | Pipeline | min_k/max_k applied 3 times with override |
| X6 | MEDIUM | Pipeline | None returned for missing required params |
| J3 | LOW | Jellyfish | Raw CalledProcessError without context |
| J4 | LOW | Jellyfish | Progress display out-of-order |
| P1 | LOW | Position | StreamingPositionCache hardcodes k=6-12 |
| P2 | LOW | Position | No HDF5 metadata attributes |
| GC5 | LOW | Strategy | Cosmetic confidence scores |
| UX1 | LOW | CLI | 20 "Missing parameter" warnings for optional params |

---

## Recommendations

### Priority 1 (Must Fix)
1. Fix CLI blacklist k-mer counting order (move Jellyfish call before step1, or integrate into step1)
2. Remove or implement --enable-qa for step1
3. Fix bg_prefixes KeyError (use .get() with empty list default)
4. Fix GC-adaptive phi29 override (check if user explicitly set polymerase)
5. Fix step4 polymerase reading (use parameter module, not raw JSON)

### Priority 2 (Should Fix)
6. Add subprocess timeout to MultiGenomeKmerCounter
7. Add text=True to class subprocess calls
8. Extend StreamingPositionCache k range for equiphi29/bst primers
9. Suppress warnings for optional params with known defaults

### Priority 3 (Nice to Have)
10. Add HDF5 metadata (genome length, source, k-mer size)
11. Unify Jellyfish calling code (eliminate MultiGenomeKmerCounter duplication)
12. Add transition zone at GC classification boundaries
