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
| `hybrid` | Medium | Excellent | Good | General use (default) |
| `dominating-set` | Fast (8x) | Excellent | Fair | Large pools, quick results |
| `background-aware` | Slow | Good | Excellent | Clinical, low background |
| `network` | Medium | Good | Good | Tm-balanced sets |
| `genetic` | Slow | Good | Good | Multi-objective exploration |
| `moea` | Slow | Good | Good | Pareto optimization |
| `milp` | Variable | Optimal | Good | Exact solutions (small sets) |
| `clique` | Medium | Good | Good | Dimer-free primer sets |
| `greedy` | Fast | Fair | Fair | Simple baseline |

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


### genetic

**Multi-objective genetic algorithm.**

```bash
neoswga optimize -j params.json --optimization-method=genetic
```

**How it works:**
- Evolutionary optimization with crossover and mutation
- Balances multiple objectives (coverage, selectivity, Tm)
- Explores diverse solution space

**Strengths:**
- Good for complex multi-objective problems
- Can escape local optima
- Explores diverse solutions

**Weaknesses:**
- Slow (many iterations)
- Non-deterministic results
- May require tuning

**When to use:**
- Complex optimization with many constraints
- When other methods give poor results
- Exploration of solution space


### moea

**Multi-objective evolutionary algorithm (Pareto optimization).**

```bash
neoswga optimize -j params.json --optimization-method=moea
```

**How it works:**
- NSGA-II style Pareto optimization
- Returns set of non-dominated solutions
- Balances coverage vs selectivity trade-offs

**Strengths:**
- Provides multiple Pareto-optimal solutions
- Explicit trade-off analysis
- Good for understanding design space

**Weaknesses:**
- Slow
- Requires choosing from multiple solutions
- Complex interpretation

**When to use:**
- Trade-off analysis
- When multiple objectives matter equally
- Research applications


### milp

**Mixed Integer Linear Programming (exact solution).**

```bash
neoswga optimize -j params.json --optimization-method=milp
```

**How it works:**
- Formulates as integer program
- Solves using branch-and-bound
- Finds provably optimal solution

**Strengths:**
- Optimal solution guaranteed (if feasible)
- Provable bounds
- Good for small problems

**Weaknesses:**
- Very slow for large problems
- Requires `mip` package
- May timeout

**When to use:**
- Small candidate pools (<100)
- When optimality proof needed
- Benchmarking other methods


### greedy

**Simple greedy breadth-first search.**

```bash
neoswga optimize -j params.json --optimization-method=greedy
```

**How it works:**
- Iteratively adds primer with best marginal improvement
- Simple coverage-based scoring

**Strengths:**
- Very fast
- Simple and predictable
- Good baseline

**Weaknesses:**
- Can get stuck in local optima
- No global optimization

**When to use:**
- Quick baseline
- Testing and development
- When speed is critical


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


### Pipeline Methods

These methods chain multiple optimizers together:

- `coverage-then-dimerfree`: Runs dominating-set for coverage, then clique filtering for dimer-free refinement
- `dimerfree-scored`: Clique-based dimer-free selection followed by network scoring
- `bg-prefilter`: Background pre-filtering before optimization
- `bg-prefilter-hybrid`: Background pre-filtering combined with hybrid optimization

```bash
neoswga optimize -j params.json --optimization-method=coverage-then-dimerfree
```

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

Approximate runtimes on typical datasets (500 candidates, 10 target primers):

| Method | Runtime | Memory |
|--------|---------|--------|
| greedy | 1-2s | Low |
| dominating-set | 2-5s | Low |
| hybrid | 10-30s | Medium |
| network | 20-60s | Medium |
| background-aware | 30-120s | Medium |
| genetic | 60-300s | Medium |
| moea | 60-300s | Medium |
| milp | Variable | High |

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
