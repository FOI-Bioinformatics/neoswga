# Optimization Method Selection Guide

This guide helps you choose the right optimization method for your SWGA primer design.

## Quick Decision Tree

```
                    Start Here
                        |
                        v
              Is this for clinical/
              diagnostic use?
                   /        \
                 Yes         No
                  |           |
                  v           v
           background-    How many
              aware       candidates?
                             |
                       /          \
                    >500         <=500
                     |             |
                     v             v
              dominating-       hybrid
                  set          (default)
```

**Summary:**
1. **Clinical/diagnostic?** Use `background-aware`
2. **Large candidate pool (>500)?** Use `dominating-set`
3. **Otherwise:** Use `hybrid` (default)

## Method Comparison

| Method | Speed | Coverage | Specificity | Best For |
|--------|-------|----------|-------------|----------|
| `hybrid` | Slow, and superlinear in set size | Excellent | Good | General use (default) |
| `dominating-set` | Fast, and flat in set size | Excellent | Fair | Large pools, quick results |
| `background-aware` | Slowest | Good | Excellent | Clinical, low background |
| `network` | Slow | Fair | Good | Tm-balanced sets |
| `clique` | Medium | Fair | Best measured | Dimer-free primer sets |
| `ensemble` | Slowest | Best of the above | Varies | When unsure which to use |

> **Read this table with the measurements below.** On the one target measured
> so far, `hybrid`, `dominating-set` and `background-aware` returned the
> *identical* primer set at every size tested, so the speed column is the only
> column separating them. `genetic`, `moea`, `milp` and `greedy` were removed
> from this guide on 2026-08-31: those optimizers no longer exist and the CLI
> rejects the names.

## Detailed Method Descriptions

### hybrid (Default)

**Two-stage approach combining coverage and connectivity.**

```bash
neoswga optimize -j params.json --optimization-method=hybrid
```

**How it works:**
1. Stage 1: Uses dominating-set to select primers with good coverage
2. Stage 2: Refines selection using network connectivity analysis

**Strengths:**
- Good balance of speed and quality
- Both coverage and amplification optimized
- Suitable for most applications

**Weaknesses:**
- Not the fastest option for very large pools
- Background minimization not explicit

**When to use:**
- General primer design
- When you want reliable results without tuning


### dominating-set

**Graph-based coverage optimization (8x faster).**

```bash
neoswga optimize -j params.json --optimization-method=dominating-set
```

**How it works:**
- Models primer selection as a set cover problem
- Greedy selection of primers covering the most uncovered regions
- Provable ln(n) approximation to optimal

**Strengths:**
- Very fast (8x faster than hybrid)
- Guaranteed coverage bounds
- Works well with large candidate pools

**Weaknesses:**
- Ignores amplification network structure
- Does not consider background binding

**When to use:**
- Large candidate pools (>500 primers)
- Quick screening runs
- When coverage is primary concern


### background-aware

**Three-stage optimizer with explicit background minimization.**

```bash
neoswga optimize -j params.json --optimization-method=background-aware
```

**How it works:**
1. Stage 1: Filter by background binding frequency
2. Stage 2: Score by selectivity ratio (fg/bg)
3. Stage 3: Optimize coverage subject to background constraint

**Strengths:**
- 10-20x reduction in background amplification
- Designed for clinical samples with high background
- Explicit selectivity optimization

**Weaknesses:**
- Slower than other methods
- May sacrifice some coverage for specificity

**When to use:**
- Clinical/diagnostic applications
- Samples with high host DNA
- When background contamination is a concern


### network

**Tm-weighted amplification network optimization.**

```bash
neoswga optimize -j params.json --optimization-method=network
```

**How it works:**
- Builds amplification network from primer binding sites
- Weights edges by Tm compatibility
- Penalizes primer pairs with dimer potential
- Selects primers maximizing network connectivity

**Strengths:**
- Considers thermodynamic compatibility
- Dimer-aware selection
- Good for uniform amplification

**Weaknesses:**
- Slower for large networks
- May not maximize coverage

**When to use:**
- When Tm uniformity is important
- Primer sets with dimer concerns
- Uniform amplification desired


### clique

**Dimer-free primer sets via clique finding.**

```bash
neoswga optimize -j params.json --optimization-method=clique
```

**How it works:**
- Builds a compatibility graph where edges connect dimer-free primer pairs
- Finds maximum clique (largest set with no primer dimers)

**Strengths:**
- Guarantees dimer-free primer sets
- Suitable when dimer avoidance is critical

**Weaknesses:**
- May sacrifice coverage for dimer-freedom
- Computationally intensive for large candidate pools

**When to use:**
- When primer-dimer formation is a primary concern
- Multiplex applications requiring dimer-free sets


### ensemble

**Runs several methods on one shared position cache and keeps the best.**

```bash
neoswga optimize -j params.json --optimization-method=ensemble
neoswga optimize -j params.json --optimization-method=ensemble --ensemble-combine=union
```

Selection is by `normalized_score`, a [0,1] value comparable across optimizers,
weighted by `--application`; the raw `score` is not comparable between methods.
The runner-up table is written to `step4_improved_df_summary.json` as
`ensemble_comparison`. `--ensemble-combine union` additionally re-optimizes over
the pooled primers from all methods, guarded so it never returns a worse set.

Note that on a target where the methods converge — as they did on the one
measured here — the comparison table will show several indistinguishable
runners-up. That is information about the target, not a failure of the ensemble.

### Host-Free Optimization

When no background genome is available, use the `--no-background` flag:

```bash
neoswga optimize -j params.json --no-background
```

This skips all background-related scoring and focuses on target genome coverage and primer compatibility.


## Command Line Examples

### Basic Usage

```bash
# Default (hybrid)
neoswga optimize -j params.json

# Explicit method selection
neoswga optimize -j params.json --optimization-method=dominating-set

# View detailed method comparison
neoswga optimize --method-guide
```

### Clinical Workflow

```bash
# For clinical samples with high host DNA
neoswga optimize -j params.json \
  --optimization-method=background-aware \
  --num-primers 10
```

### Fast Screening

```bash
# Quick screening with large candidate pool
neoswga optimize -j params.json \
  --optimization-method=dominating-set \
  --num-primers 15
```

### High-Quality Design

```bash
# Thorough optimization for important designs
neoswga optimize -j params.json \
  --optimization-method=hybrid \
  --iterations 20
```

## Iterative Design Workflow

For iterative wet-lab optimization:

```bash
# 1. Design initial set
neoswga optimize -j params.json --optimization-method=hybrid

# 2. Predict efficiency before synthesis
neoswga predict-efficiency -j params.json \
  --primers SEQ1 SEQ2 SEQ3 SEQ4 SEQ5 SEQ6

# 3. After testing, expand with additional primers
neoswga expand-primers -j params.json \
  --fixed-primers SEQ1 SEQ2 SEQ3 \
  --failed-primers SEQ4 \
  --num-new 4 \
  --output expanded/
```

## Performance Benchmarks

Measured 2026-08-31 on a 3.17 Mb bacterial target against human chr21, 2000
candidates, seed 42, phi29 at 30 C, one run each (repeat noise about 9%). Driven
through the CLI, so the figures include dispatch, position-cache build and
metric computation. Reproduce with `scripts/benchmarking/sweep_optimize.py`.

| Method | 12 primers | 32 primers | 64 primers | 128 primers | fg_coverage at 32 |
|--------|---:|---:|---:|---:|---:|
| dominating-set | 8.3 s | 8.2 s | 9.0 s | 8.6 s | 0.5376 |
| hybrid (default) | 13.7 s | 63.3 s | 347 s | 2239 s | 0.5376 |
| background-aware | 20.1 s | 133 s | — | — | 0.5376 |
| clique | 13.1 s | 34.7 s | — | — | 0.4499 |
| network | 64.9 s | 377 s | — | — | 0.4283 |

Three things this table is for. **`dominating-set` is flat in set size** while
coverage climbs; every other method's cost curve is superlinear. **The first
three rows returned the identical primer set** at every size tested (Jaccard
1.000), so the extra runtime bought nothing on this target. **`network` is worse
on coverage and costs 46x more** at 32 primers.

`clique` is the one method that is not dominated: it gives up coverage but was
the only one to return `dimer_risk` 0.0000, and it had the best selectivity
density of any method here (42.4 against 31.7).

Caveats that belong with these numbers: one target, one k, one run each. Enough
to show that a method returning an identical set at 38x the cost has a defect;
not enough to change a shipped default on. See
[`AUDIT_2026-08_alternatives_and_scaling.md`](AUDIT_2026-08_alternatives_and_scaling.md)
(F5, F5b) for the full measurement and its limits.

## Troubleshooting

### Poor coverage

Try:
1. Use `dominating-set` for coverage-focused optimization
2. Increase `--num-primers`
3. Check if candidates have sufficient binding sites

### High background

Try:
1. Use `background-aware` method
2. Pre-filter with stricter `max_bg_freq`
3. Check background genome is correctly specified

### Slow optimization

Try:
1. Use `dominating-set` for speed
2. Reduce candidate pool size
3. Decrease `--iterations`

### Poor Tm uniformity

Try:
1. Use `network` method
2. Adjust `min_tm` and `max_tm` parameters
3. Consider polymerase-specific presets
