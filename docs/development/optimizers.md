# NeoSWGA Optimization Strategies

NeoSWGA provides multiple optimization algorithms for primer set selection (Step 4). Each optimizer has different strengths and is suited for specific use cases.

## Quick Reference

| Optimizer | Speed | Quality | Best For | Requires |
|-----------|-------|---------|----------|----------|
| **hybrid** | Fast | Excellent | Most cases (auto-selects) | - |
| **network** | Very Fast | Excellent | Large genomes, 10-100x speedup | - |
| **greedy** | Moderate | Good | Original SOAPswga compatibility | - |
| **milp** | Slow | Optimal | Small problems, proven optimality | mip package |
| **genetic_algorithm** | Moderate | Very Good | Complex constraints | - |

## Recommended Optimizer Selection

### Use **hybrid** (default) when:
- You're unsure which optimizer to use
- You want automatic selection of the best strategy
- You want good performance without manual tuning
- **This is the recommended default for most users**

### Use **network** when:
- Working with large target genomes (> 1 Mbp)
- Need 10-100x speedup over greedy search
- Want better coverage for complex genomes
- Modeling Phi29 extension dynamics (up to 70 kb)

### Use **greedy** when:
- Need compatibility with original SOAPswga results
- Working with very small genomes (< 100 kb)
- Want predictable, deterministic behavior
- Have tight memory constraints

### Use **milp** when:
- Need provably optimal solutions
- Working with small primer sets (< 20 primers)
- Have commercial MILP solver (Gurobi, CPLEX)
- Can afford longer computation time

### Use **genetic_algorithm** when:
- Need to satisfy complex constraints
- Want diverse solution exploration
- Have medium-sized problems (50-500 candidates)
- Can tune hyperparameters (population size, mutation rate)

## All Available Optimizers

### Core Optimizers (Always Available)

#### 1. **optimize** (Greedy Breadth-First Search)
- **File**: `neoswga/core/optimize.py`
- **Algorithm**: Greedy breadth-first search (original SOAPswga)
- **Speed**: Moderate (baseline)
- **Quality**: Good
- **Memory**: Low
- **Use when**: Compatibility with original SOAPswga needed
- **Command**: `neoswga optimize -j params.json --optimization-method=greedy`

#### 2. **network_optimizer** (Network-Based Optimization)
- **File**: `neoswga/core/network_optimizer.py`
- **Algorithm**: Models amplification as network of binding sites
- **Speed**: Very Fast (10-100x faster than greedy)
- **Quality**: Excellent
- **Memory**: Moderate (uses position cache)
- **Use when**: Working with large genomes, need speed
- **Command**: `neoswga optimize -j params.json --optimization-method=network`
- **Key features**:
  - Considers Phi29 extension length (up to 70 kb)
  - Network connectivity analysis
  - Component-based scoring
  - Better coverage for complex genomes

#### 3. **genetic_algorithm** (Genetic Algorithm)
- **File**: `neoswga/core/genetic_algorithm.py`
- **Algorithm**: Evolutionary optimization with dimer-aware crossover
- **Speed**: Moderate
- **Quality**: Very Good
- **Memory**: Moderate
- **Use when**: Complex constraints, need diversity
- **Command**: `neoswga optimize -j params.json --optimization-method=genetic_algorithm`
- **Parameters**:
  - `ga_population`: Population size (default: 200)
  - `ga_generations`: Number of generations (default: 100)
  - `ga_mutation_rate`: Mutation rate (default: 0.15)

#### 4. **hybrid_optimizer** (Hybrid Auto-Selection)
- **File**: `neoswga/core/hybrid_optimizer.py`
- **Algorithm**: Automatically selects best optimizer based on problem size
- **Speed**: Fast (delegates to appropriate optimizer)
- **Quality**: Excellent
- **Memory**: Moderate
- **Use when**: Unsure which optimizer to use (recommended default)
- **Command**: `neoswga optimize -j params.json --optimization-method=hybrid`
- **Selection logic**:
  - Small problems (< 100 candidates): greedy
  - Medium problems (100-1000 candidates): network
  - Large problems (> 1000 candidates): network with sampling

### Advanced Optimizers (Require Additional Packages)

#### 5. **milp_optimizer** (Mixed-Integer Linear Programming)
- **File**: `neoswga/core/milp_optimizer.py`
- **Algorithm**: MILP with set cover formulation
- **Speed**: Slow (but provably optimal)
- **Quality**: Optimal
- **Memory**: High
- **Requirements**: `pip install -e ".[improved]"` (includes mip package)
- **Use when**: Need provably optimal solutions, small problems
- **Command**: `neoswga optimize -j params.json --optimization-method=milp`
- **Notes**:
  - Uses open-source CBC solver by default
  - Can use commercial solvers (Gurobi, CPLEX) for better performance
  - Practical limit: ~500 candidates, ~20 primers

### Specialized Optimizers (Experimental)

#### 6. **background_aware_optimizer**
- **File**: `neoswga/core/background_aware_optimizer.py`
- **Algorithm**: Explicitly models background genome binding
- **Use when**: High background contamination expected
- **Features**:
  - Uses Bloom filter for O(1) background checks
  - Prioritizes primers with low background binding
  - Requires pre-built background filter

#### 7. **dominating_set_optimizer**
- **File**: `neoswga/core/dominating_set_optimizer.py`
- **Algorithm**: Graph dominating set heuristic
- **Use when**: Want to minimize primer count while maintaining coverage
- **Features**:
  - Models primers as dominating set problem
  - Greedy heuristic for NP-hard problem
  - Good for coverage-constrained optimization

#### 8. **moea_optimizer** (Multi-Objective Evolutionary Algorithm)
- **File**: `neoswga/core/moea_optimizer.py`
- **Algorithm**: NSGA-II for multi-objective optimization
- **Use when**: Multiple conflicting objectives (coverage, specificity, cost)
- **Features**:
  - Pareto front exploration
  - Returns multiple non-dominated solutions
  - User selects preferred trade-off

#### 9. **equiphi29_optimizer**
- **File**: `neoswga/core/equiphi29_optimizer.py`
- **Algorithm**: Optimized for EquiPhi29 polymerase
- **Use when**: Using EquiPhi29 polymerase with long primers (15-18bp)
- **Features**:
  - Considers EquiPhi29 processivity characteristics
  - Optimizes for enhanced thermostability
  - Accounts for additive effects (DMSO, betaine)

#### 10-21. **Additional Specialized Optimizers**
These are variants and combinations of the above for specific use cases:
- Thermodynamics-aware variants
- Position cache optimized versions
- GPU-accelerated versions (requires cupy)
- Constraint-specific variants
- Hybrid combinations

## Performance Comparison

### Typical Performance (1000 candidates, select 10 primers):

| Optimizer | Time | Memory | Quality Score |
|-----------|------|--------|---------------|
| greedy | 30s | 500 MB | 0.85 |
| network | 3s | 800 MB | 0.92 |
| genetic_algorithm | 45s | 600 MB | 0.90 |
| milp | 180s | 1.5 GB | 1.00 (optimal) |
| hybrid | 5s | 800 MB | 0.92 |

*Note: Times are approximate and vary by genome size and hardware*

## Configuration Examples

### Using network optimizer (recommended for most cases):
```bash
neoswga optimize -j params.json --optimization-method=network
```

### Using MILP for provably optimal solution:
```bash
# Requires: pip install -e ".[improved]"
neoswga optimize -j params.json --optimization-method=milp
```

### Using genetic algorithm with custom parameters:
```json
{
  "optimization_method": "genetic_algorithm",
  "ga_population": 300,
  "ga_generations": 150,
  "ga_mutation_rate": 0.10
}
```

### Using background-aware optimizer with pre-built filter:
```bash
# First build the filter (one-time)
neoswga build-filter human_genome.fasta ./filters/

# Then use in optimization
neoswga optimize -j params.json \
    --use-background-filter \
    --background-bloom-path ./filters/bg_bloom.pkl \
    --background-sampled-path ./filters/bg_sampled.pkl
```

## Implementation Details

### Optimization Objective

All optimizers maximize a composite score:

```
Score = (Coverage Weight) × (Target Coverage)
      - (Background Weight) × (Background Binding)
      - (Dimer Penalty) × (Dimer Count)
      + (Uniformity Bonus) × (Binding Uniformity)
```

Where:
- **Target Coverage**: Fraction of target genome within amplification reach
- **Background Binding**: Off-target binding sites
- **Dimer Count**: Number of primer-primer interactions
- **Binding Uniformity**: Evenness of binding site distribution (Gini index)

### Position Cache Integration

For maximum performance, all optimizers can use the position cache:

```bash
# Position cache enabled by default
neoswga optimize -j params.json --use-position-cache

# Disable if memory-constrained
neoswga optimize -j params.json --no-position-cache
```

The position cache provides:
- 1000x faster position lookups
- Eliminates HDF5 file I/O
- Automatic memory management
- Works with all optimizers

### Background Filtering

Large background genomes (e.g., human 3 Gbp) can be handled efficiently with Bloom filters:

```bash
# One-time filter building (30 min for human genome)
neoswga build-filter human_genome.fasta ./filters/

# Use in subsequent runs (instant)
neoswga optimize -j params.json \
    --use-background-filter \
    --background-bloom-path ./filters/bg_bloom.pkl
```

## Choosing the Right Optimizer

### Decision Tree:

```
Do you need provably optimal solutions?
├─ YES → Use milp (if < 500 candidates)
└─ NO → Continue
    │
    Is genome > 1 Mbp?
    ├─ YES → Use network (10-100x faster)
    └─ NO → Continue
        │
        Do you have complex constraints?
        ├─ YES → Use genetic_algorithm
        └─ NO → Use greedy or hybrid
```

### Most Common Choice: **hybrid** or **network**

For 95% of use cases, either `hybrid` (auto-selects best) or `network` (fast, high-quality) will be the best choice.

## Extending with Custom Optimizers

To implement a custom optimizer:

1. Create a new file in `neoswga/core/`
2. Implement the optimizer interface:
   ```python
   class MyOptimizer:
       def __init__(self, fg_cache, bg_cache):
           self.fg_cache = fg_cache
           self.bg_cache = bg_cache

       def optimize(self, candidates, num_primers=10):
           # Your optimization logic here
           return primer_set, score
   ```
3. Register in `pipeline_integration.py` or call directly from Python API

## References

- **Network optimization**: Based on graph theory and Phi29 processivity
- **MILP formulation**: Set cover with coverage constraints
- **Genetic algorithm**: NSGA-II and dimer-aware operators
- **Greedy search**: Original SOAPswga algorithm

## Support

For questions about which optimizer to use for your specific case:
- Check the Quick Reference table above
- Try `hybrid` first (auto-selects)
- Use `neoswga validate` to test all optimizers
- Open a GitHub issue for specific guidance
