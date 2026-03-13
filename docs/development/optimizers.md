# NeoSWGA Optimization Strategies

NeoSWGA provides multiple optimization algorithms for primer set selection (Step 4). All optimizers inherit from `BaseOptimizer` and are registered with `OptimizerFactory` for uniform instantiation.

## Quick Reference

| Method Name | Speed | Best For | Requires |
|-------------|-------|----------|----------|
| `hybrid` | Medium | General use (default) | - |
| `greedy` | Fast | Simple optimization, small genomes | - |
| `network` | Medium | Tm-weighted, large genomes | - |
| `dominating-set` | Fast | Set-cover, large pools | - |
| `weighted-set-cover` | Fast | Score-weighted coverage | - |
| `background-aware` | Slow | Clinical, background reduction | - |
| `genetic` | Moderate | Evolutionary exploration | - |
| `equiphi29` | Medium | EquiPhi29 at 42-45C | - |
| `tiling` | Fast | Interval-based coverage | - |
| `normalized` | Medium | Strategy-preset scoring | - |
| `clique` | Moderate | Dimer-free via max-clique | networkx |
| `multi-agent` | Slow | Parallel ensemble | - |
| `bg-prefilter` | Medium | Background pruning wrapper | - |
| `coverage-then-dimerfree` | Medium | DS then clique cascade | networkx |
| `dimerfree-scored` | Medium | Clique then network cascade | networkx |
| `bg-prefilter-hybrid` | Medium | BG pre-filter then hybrid | - |
| `moea` | Slow | Pareto multi-objective | pymoo |
| `milp` | Variable | Exact optimal solutions | mip |

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

#### 10. **tiling_optimizer** (Interval-Based Tiling)
- **File**: `neoswga/core/tiling_optimizer.py`
- **Algorithm**: Model primer binding as genomic intervals, tile genome with minimal gaps
- **Command**: `neoswga optimize -j params.json --optimization-method=tiling`

#### 11. **normalized_optimizer** (Strategy Presets)
- **File**: `neoswga/core/normalized_optimizer.py`
- **Algorithm**: Normalized [0,1] scoring with application-specific weight presets
- **Command**: `neoswga optimize -j params.json --optimization-method=normalized`

#### 12. **clique_optimizer** (Dimer-Free Clique)
- **File**: `neoswga/core/clique_optimizer.py`
- **Algorithm**: Build compatibility graph, find maximum clique for dimer-free sets
- **Requirements**: `networkx`
- **Command**: `neoswga optimize -j params.json --optimization-method=clique`

#### 13. **multi_agent_optimizer** (Parallel Ensemble)
- **File**: `neoswga/core/multi_agent_optimizer.py`
- **Algorithm**: Run multiple optimizers in parallel, aggregate via best/union/voting/pareto
- **Command**: `neoswga optimize -j params.json --optimization-method=multi-agent`

#### 14. **background_prefilter** (BG Pruning Wrapper)
- **File**: `neoswga/core/background_prefilter.py`
- **Algorithm**: Prune candidates by fg/bg ratio, then delegate to inner optimizer
- **Command**: `neoswga optimize -j params.json --optimization-method=bg-prefilter`

#### 15-17. **Serial Cascade Optimizers**
- **File**: `neoswga/core/serial_cascade_optimizer.py`
- `coverage-then-dimerfree`: Dominating-set then clique cascade
- `dimerfree-scored`: Clique then network scoring cascade
- `bg-prefilter-hybrid`: Background pre-filter then hybrid
- **Command**: `neoswga optimize -j params.json --optimization-method=coverage-then-dimerfree`

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

## BaseOptimizer Interface

All optimizers inherit from `BaseOptimizer` (in `base_optimizer.py`), which enforces a consistent interface:

```python
from neoswga.core.base_optimizer import BaseOptimizer, OptimizationResult, OptimizerConfig

class MyOptimizer(BaseOptimizer):
    """Custom optimizer example."""

    def optimize(
        self,
        candidates: List[str],
        target_size: int,
        **kwargs
    ) -> OptimizationResult:
        # Implementation here
        ...
```

Key data classes:
- `OptimizerConfig`: Common configuration (target_size, max_iterations, timeout_seconds, min_coverage)
- `OptimizationResult`: Frozen result container (primer_set, metrics, status, message)
- `PrimerSetMetrics`: Comprehensive metrics (coverage, selectivity, dimer risk, gap stats, strand alternation) with a `normalized_score()` method for cross-optimizer comparison
- `OptimizationStatus`: Enum (SUCCESS, PARTIAL, NO_CONVERGENCE, ERROR)
- `CompositeOptimizer`: Chains multiple optimizers sequentially

## OptimizerFactory and Registry

Optimizers are registered via the decorator pattern and instantiated through `OptimizerFactory`:

```python
from neoswga.core.optimizer_factory import OptimizerFactory

# Register a new optimizer
@OptimizerFactory.register('my-optimizer', aliases=['my', 'custom'],
                           description='My custom optimization strategy')
class MyOptimizer(BaseOptimizer):
    ...

# Instantiate by name
optimizer = OptimizerFactory.create(
    'my-optimizer',
    cache=position_cache,
    fg_prefixes=fg_prefixes,
    fg_seq_lengths=fg_seq_lengths
)

# List all registered optimizers
for name, desc in OptimizerFactory.list_optimizers().items():
    print(f"{name}: {desc}")
```

Registration happens automatically when optimizer modules are imported. The `unified_optimizer._ensure_optimizers_registered()` function triggers all imports at first use.

## Normalized Scoring

`PrimerSetMetrics.normalized_score()` computes a weighted [0,1] composite:

| Component | Default Weight | Metric |
|-----------|---------------|--------|
| Coverage | 0.35 | fg_coverage |
| Selectivity | 0.30 | selectivity_ratio (clamped) |
| Dimer safety | 0.15 | 1 - dimer_risk_score |
| Evenness | 0.10 | 1 - gap_gini |
| Tm uniformity | 0.10 | Tm range penalty |

The `NormalizedOptimizer` provides strategy presets (discovery, clinical, enrichment, metagenomics) that adjust these weights.

## Serial Cascade Optimizers

The `SerialCascadeOptimizer` base class enables chaining optimizers in sequence. Three pre-built cascades are registered:

- `coverage-then-dimerfree`: Dominating-set for coverage, then clique for dimer-free refinement
- `dimerfree-scored`: Clique dimer-free selection, then network scoring
- `bg-prefilter-hybrid`: Background pre-filter pruning, then hybrid optimization

## Extending with Custom Optimizers

1. Create a new file in `neoswga/core/`
2. Inherit from `BaseOptimizer` and implement `optimize()`
3. Decorate with `@OptimizerFactory.register('name', aliases=[...])`
4. Add the import to `unified_optimizer._ensure_optimizers_registered()`

## References

- **Network optimization**: Graph connectivity and Phi29 processivity modeling
- **MILP formulation**: Set cover with coverage constraints
- **Genetic algorithm**: Tournament selection with dimer-aware crossover
- **Greedy search**: Original SOAPswga breadth-first algorithm
- **Clique**: Maximum clique enumeration on primer compatibility graph
- **Dominating set**: Greedy set cover with ln(n) approximation guarantee
