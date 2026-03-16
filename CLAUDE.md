# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

NeoSWGA is a command-line tool for selecting primer sets for selective whole-genome amplification (SWGA). See [README.md](README.md) for user-facing documentation and quick start.

**External dependency**: Jellyfish k-mer counter must be in PATH.

## Architecture

### Entry Point

- `neoswga/cli_unified.py`: Main CLI entry point (all commands)
- Entry point defined in `pyproject.toml`: `neoswga = neoswga.cli_unified:main`

### Core Modules (neoswga/core/, 75 modules)

**Pipeline implementation**:
- `pipeline.py`: Main pipeline functions (count-kmers, filter, score, optimize)
- `improved_pipeline.py`: Performance-optimized pipeline with background filtering
- `pipeline_qa_integration.py`: Quality assurance hooks for filtering
- `auto_swga_pipeline.py`: Automatic parameter optimization pipeline
- `multi_genome_pipeline.py`: Pan-genome primer design for multiple targets

**Filtering (filter command)**:
- `filter.py`: Frequency, Gini index, GC content, complexity filters
- `adaptive_filters.py`: Adaptive GC filtering for extreme GC genomes
- `kmer.py`: K-mer counting via jellyfish subprocess
- `kmer_counter.py`: Multi-genome k-mer counter using Jellyfish
- `multi_genome_filter.py`: Filter primers across multiple target genomes
- `thermodynamic_filter.py`: Filter by Tm, secondary structure
- `thermo_background_filter.py`: Thermodynamic-aware background filtering

**Scoring (score command)**:
- `rf_preprocessing.py`: Random forest feature engineering
- `primer_attributes.py`: Tm calculation, self-complementarity
- `integrated_quality_scorer.py`: Multi-criteria quality scoring
- Pre-trained model: `neoswga/core/models/random_forest_filter.p`

**Optimization (optimize command)**:
- `network_optimizer.py`: Network-based optimization (default, includes Tm weighting and dimer penalty)
- `genetic_algorithm.py`: GA-based optimization
- `optimize.py`: Original greedy breadth-first search
- `milp_optimizer.py`: Mixed-integer linear programming (requires `mip` package)
- `hybrid_optimizer.py`: Hybrid approach combining multiple strategies
- `dominating_set_optimizer.py`: Graph-based greedy set cover (8x faster, ln(n) approximation)
- `background_aware_optimizer.py`: Three-stage optimizer with explicit background minimization (10-20x reduction)
- `moea_optimizer.py`: Multi-objective evolutionary algorithm
- `equiphi29_optimizer.py`: EquiPhi29-specific optimization at 42-45C
- `optimal_oligo_generator.py`: Alternative comprehensive primer design
- `minimal_primer_selector.py`: Post-process to minimize primer count

**Thermodynamics**:
- `thermodynamics.py`: SantaLucia nearest-neighbor calculations with LRU caching (1M entries)
- `thermo_estimation.py`: Free energy with salt corrections
- `reaction_conditions.py`: Polymerase presets (phi29, equiphi29, bst, klenow), additive effects
- `secondary_structure.py`: Hairpin and dimer prediction
- `dimer.py`: Primer-dimer calculations
- `three_prime_stability.py`: 3' end stability analysis

**Performance**:
- `position_cache.py`: In-memory position cache (1000x speedup)
- `background_filter.py`: Bloom filter for large background genomes
- `gpu_acceleration.py`: CuPy-based GPU acceleration for thermodynamics (10-100x speedup)

**Simulation**:
- `replication_simulator.py`: Agent-based phi29 DNA replication simulation
- `swga_simulator.py`: SWGA reaction simulation
- `stochastic_simulator.py`: Gillespie algorithm stochastic simulation
- `simulation_fitness.py`: Fitness functions for simulation
- `simulation_analysis.py`: Analysis of simulation results
- `simulation_plots.py`: Visualization of simulation results
- `simulation_report.py`: HTML report generation

**Analysis**:
- `genome_analysis.py`: Genome suitability analysis for SWGA
- `genome_io.py`: Genome file I/O utilities
- `strand_bias_analyzer.py`: Strand bias detection
- `dimer_network_analyzer.py`: Primer dimer interaction network
- `amplicon_network.py`: Amplicon network analysis

**Utilities**:
- `utility.py`: Multiprocessing helpers, file utilities
- `parameter.py`: Global parameter handling
- `string_search.py`: String matching algorithms
- `validation.py`: Installation validation and testing framework

**User Experience (New)**:
- `wizard.py`: Setup wizard for guided params.json creation
- `param_validator.py`: Parameter validation with error/warning/info levels
- `condition_suggester.py`: Reaction condition recommendations
- `results_interpreter.py`: Quality assessment and go/no-go recommendations
- `workflow_selector.py`: Interactive menu for feature discovery

**Mechanistic Modeling**:
- `mechanistic_params.py`: Literature-based parameters for four-pathway model
- `mechanistic_model.py`: Four-pathway mechanistic model (Tm, accessibility, enzyme, kinetics)
- `set_size_optimizer.py`: Automatic primer set size optimization by application profile
- `model_validation.py`: Validation tests against expected literature behavior
- `additives.py`: Additive effects on Tm with sigmoid GC normalization

**Reporting**:
- `report/`: Quality report generation module
  - `metrics.py`: Collect pipeline metrics from results
  - `quality.py`: Quality grading (A-F) with component scoring
  - `executive_summary.py`: One-page HTML summary report
  - `technical_report.py`: Comprehensive technical analysis report
  - `visualizations.py`: Interactive Plotly charts (optional dependency)
  - `validation.py`: Input validation before report generation
  - `utils.py`: Shared utilities and chart color schemes

**Experimental/Advanced** (in `neoswga/core/experimental/`, not fully integrated):
- `cooperative_binding_selector.py`: Network-based cooperative binding model
- `hybrid_primer_strategy.py`: Core + supplemental mixed-length primer strategy
- `simulate_command.py`: Alternative simulation CLI implementation

**Experimental** (in main core/, not fully validated):
- `active_learning.py`: Active learning for iterative primer optimization
- `adaptive_search.py`: Adaptive search strategies
- `advanced_features.py`: Advanced feature engineering
- `deep_learning.py`: Deep learning-based prediction (stub, requires training data)
- `gc_adaptive_strategy.py`: GC-adaptive primer design

### Data Flow

```
count-kmers            filter                 score                  optimize
     |                    |                     |                       |
     v                    v                     v                       v
 *_Xmer_all.txt  -->  step2_df.csv +    -->  step3_df.csv      -->  step4_improved_df.csv
 (k-mer counts)       positions.h5           (with amp_pred)        (final primer sets)
```

**File outputs** (in `data_dir`):
- `step2_df.csv`: Filtered primers with fg_freq, bg_freq, gini, Tm
- `step3_df.csv`: Primers with amplification prediction scores
- `step4_improved_df.csv`: Final optimized primer sets with enrichment scores
- `*_positions.h5`: HDF5 files with primer binding positions

## CLI Commands

### Standard Pipeline
```bash
neoswga count-kmers -j params.json  # Step 1: Generate k-mer counts
neoswga filter -j params.json       # Step 2: Filter candidate primers
neoswga score -j params.json        # Step 3: Score amplification efficacy
neoswga optimize -j params.json     # Step 4: Find optimal primer sets
```

### Optimization Methods
```bash
# Standard (default)
neoswga optimize -j params.json --optimization-method=hybrid

# Fast graph-based (8x faster)
neoswga optimize -j params.json --optimization-method=dominating-set

# Clinical (10-20x background reduction)
neoswga optimize -j params.json --optimization-method=background-aware

# Other options: greedy, milp, network, genetic, moea
```

**Optimization Method Comparison**:

| Method | Speed | Best For | Notes |
|--------|-------|----------|-------|
| `hybrid` | Medium | General use (default) | Combines network + greedy approaches |
| `dominating-set` | Fast (8x) | Large primer pools | Graph-based set cover, ln(n) approximation |
| `background-aware` | Slow | Clinical applications | 10-20x background reduction, three-stage |
| `greedy` | Fast | Simple optimization | Original breadth-first search |
| `network` | Medium | Tm-weighted selection | Dimer penalty aware |
| `genetic` | Slow | Complex multi-objective | GA-based, good for exploration |
| `moea` | Slow | Pareto optimization | Multi-objective evolutionary |
| `milp` | Variable | Exact solutions | Requires `mip` package |

### Utility Commands
```bash
neoswga validate --quick            # Validate installation
neoswga build-filter genome.fna ./  # Build Bloom filter for large background
neoswga show-presets                # Show reaction condition presets
```

### Setup Commands (New)
```bash
# Interactive workflow selector - discover all features
neoswga start

# Setup wizard - create params.json with guided configuration
neoswga init --genome target.fna [--background host.fna] [-o params.json]

# Validate params.json before running pipeline
neoswga validate-params -j params.json

# Suggest optimal reaction conditions
neoswga suggest --genome-gc 0.65 --primer-length 15
neoswga suggest --genome target.fna  # Auto-calculates GC

# Interpret results after pipeline completes
neoswga interpret -d results/

# Generate quality report after pipeline completes
neoswga report -d results/                        # Executive summary (default)
neoswga report -d results/ --level full           # Full technical report
neoswga report -d results/ --interactive          # With interactive Plotly charts
neoswga report -d results/ --level full --interactive  # Full report with charts
neoswga report -d results/ --check                # Validate only, don't generate

# Validate mechanistic model against expected behavior
neoswga validate-model               # Run all validation tests
neoswga validate-model --output-json  # Output results as JSON
```

### Optimization with Mechanistic Model
```bash
# Auto-size primer set based on application profile
neoswga optimize -j params.json --auto-size --application clinical

# Applications: discovery (high coverage), clinical (high specificity),
#              enrichment (balanced), metagenomics (capture diversity)

# Use mechanistic model for primer weighting
neoswga optimize -j params.json --use-mechanistic-model --mechanistic-weight 0.3
```

### Advanced Commands
```bash
# Multi-genome pan-primer design
neoswga multi-genome --genomes target1.fna target2.fna --output results/

# Replication simulation
neoswga simulate --primers SEQ1 SEQ2 --genome target.fna --output sim/

# Analyze existing primer set
neoswga analyze-set --primers SEQ1 SEQ2 --fg target.fna --fg-kmers data/target --output analysis/

# Genome analysis
neoswga analyze-genome --genome target.fna --output analysis/

# Dimer network analysis
neoswga analyze-dimers --primers SEQ1 SEQ2 --output dimers/ --visualize

# Active learning for iterative optimization (experimental)
neoswga active-learn -j params.json --output active_learn/ --num-candidates 10
```

## Key Parameters (params.json)

**Primer filtering**:
- `min_k`, `max_k`: Primer length range (default: 6-12, use 12-18 for longer primers)
- `min_fg_freq`: Minimum foreground frequency (default: 1e-5)
- `max_bg_freq`: Maximum background frequency (default: 5e-6)
- `max_gini`: Maximum Gini index for binding evenness (default: 0.6)
- `max_primer`: Primers to keep after filtering (default: 500)

**Thermodynamics**:
- `polymerase`: "phi29" (30C), "equiphi29" (42-45C), "bst" (60-65C), "klenow" (25-40C)
- `reaction_temp`: Reaction temperature in Celsius
- `na_conc`, `mg_conc`: Salt concentrations (mM)
- `dmso_percent`, `betaine_m`, `trehalose_m`: Common additive concentrations
- `ethanol_percent`, `urea_m`, `tmac_m`, `formamide_percent`: Advanced additives
- `min_tm`, `max_tm`: Melting temperature range

**Polymerase Presets**:

| Polymerase | Temp | Primer Length | Use Case |
|------------|------|---------------|----------|
| `phi29` | 30C | 6-12 bp | Standard SWGA, high processivity |
| `equiphi29` | 42-45C | 12-18 bp | Higher specificity, GC-rich targets |
| `bst` | 60-65C | 15-25 bp | LAMP-like applications, thermostable |
| `klenow` | 25-40C | 8-15 bp | Room temperature, lower processivity |

**Optimization**:
- `optimization_method`: 'hybrid' (default), 'dominating-set' (fast), 'background-aware' (clinical), 'greedy', 'milp', 'network', 'genetic', 'moea'
- `num_primers`, `target_set_size`: Desired primer set size (default: 6)
- `iterations`: Search iterations (default: 8)
- `max_sets`: Parallel primer sets to build (default: 5)

## Key APIs

### PositionCache

In-memory cache for primer binding positions (1000x speedup over HDF5 reads):

```python
from neoswga.core.position_cache import PositionCache

# Initialize with file prefixes and primers to load
cache = PositionCache(fname_prefixes=['data/target'], primers=['ATCGATCG', 'GCTAGCTA'])

# Get positions - strand is 'forward', 'reverse', or 'both' (NOT '+'/'-')
positions = cache.get_positions('data/target', 'ATCGATCG', strand='both')
```

### GPU Acceleration

Optional CuPy-based GPU acceleration:

```python
from neoswga.core.gpu_acceleration import is_gpu_available, get_gpu_info, GPUThermodynamics

if is_gpu_available():
    info = get_gpu_info()
    print(f"GPU: {info['device_name']}")

    gpu_calc = GPUThermodynamics(conditions)
    tms = gpu_calc.batch_calculate_tm(primers)  # 10-100x faster
```

### Reaction Conditions

```python
from neoswga.core.reaction_conditions import get_enhanced_conditions, get_standard_conditions

# Standard Phi29 at 30C
conditions = get_standard_conditions()

# EquiPhi29 with DMSO/betaine for longer primers
conditions = get_enhanced_conditions()  # 42C, 5% DMSO, 1M betaine
```

### Dominating Set Optimizer

Graph-based coverage optimization (8x faster than greedy):

```python
from neoswga.core.dominating_set_optimizer import DominatingSetOptimizer
from neoswga.core.position_cache import PositionCache

cache = PositionCache(fg_prefixes, candidates)
optimizer = DominatingSetOptimizer(cache, fg_prefixes, fg_seq_lengths, bin_size=10000)

result = optimizer.optimize_greedy(candidates, max_primers=10)
print(f"Coverage: {result['coverage']:.1%}")
```

### Mechanistic Model

Four-pathway model for additive effects on SWGA:

```python
from neoswga.core.mechanistic_model import MechanisticModel, MechanisticEffects
from neoswga.core.reaction_conditions import ReactionConditions

# Create conditions with additives
conditions = ReactionConditions(
    temp=42.0,
    polymerase='equiphi29',
    dmso_percent=5.0,
    betaine_m=1.0,
    mg_conc=2.5
)

# Calculate mechanistic effects for a primer
model = MechanisticModel(conditions)
effects = model.calculate_effects('ATCGATCGATCG', template_gc=0.5)

# Access individual pathway factors
print(f"Processivity factor: {effects.processivity_factor:.2f}")
print(f"Accessibility factor: {effects.accessibility_factor:.2f}")
print(f"Binding rate: {effects.effective_binding_rate:.2f}")
print(f"Predicted amplification: {effects.predicted_amplification_factor:.2f}")
```

**Four pathways**:
1. **Tm modification**: DMSO, betaine, formamide effects on primer-template stability
2. **Secondary structure accessibility**: Template melting, GC-dependent structure
3. **Enzyme activity**: Polymerase processivity, speed, stability
4. **Binding kinetics**: Association/dissociation rates (kon/koff)

### Set Size Optimizer

Automatic primer set size based on application profile:

```python
from neoswga.core.set_size_optimizer import recommend_set_size, quick_size_estimate
from neoswga.core.mechanistic_model import MechanisticModel
from neoswga.core.reaction_conditions import ReactionConditions

# Quick estimate without full mechanistic model
size = quick_size_estimate('clinical', genome_length=1_000_000)

# Full recommendation with mechanistic effects
conditions = ReactionConditions(temp=42.0, polymerase='equiphi29')
model = MechanisticModel(conditions)
effects = model.calculate_effects('ATCGATCGATCG', template_gc=0.5)

recommendation = recommend_set_size(
    application='clinical',
    genome_length=1_000_000,
    primer_length=12,
    mech_effects=effects
)
print(f"Recommended: {recommendation['recommended_size']} primers")
print(f"Range: {recommendation['size_range']}")
print(f"Target coverage: {recommendation['target_coverage']:.0%}")
```

**Application profiles**:

| Application | Coverage Target | Specificity | Typical Size | Use Case |
|-------------|-----------------|-------------|--------------|----------|
| `discovery` | 90% | 60% | 10-15 | Pathogen discovery, maximize sensitivity |
| `clinical` | 70% | 90% | 6-10 | Diagnostics, minimize false positives |
| `enrichment` | 80% | 75% | 8-12 | Sequencing enrichment, balanced |
| `metagenomics` | 95% | 50% | 15-20 | Capture diversity |

### Report Visualizations

Interactive Plotly charts for reports (optional dependency):

```python
from neoswga.core.report.visualizations import (
    is_plotly_available,
    render_filtering_funnel,
    render_component_radar,
    render_tm_gc_distribution,
    render_coverage_specificity_scatter,
    render_primer_heatmap,
    render_dimer_network_heatmap,
    render_dimer_network_graph,
)

# Check if Plotly is installed
if is_plotly_available():
    # Render filtering funnel chart
    funnel_html = render_filtering_funnel([
        ("Total k-mers", 100000),
        ("After frequency filter", 50000),
        ("After background filter", 10000),
        ("Final candidates", 100),
    ])

    # Render dimer interaction heatmap
    dimer_html = render_dimer_network_heatmap(primers, max_primers=15)
```

**Available chart functions**:
- `render_filtering_funnel()`: Pipeline filtering stages as funnel chart
- `render_component_radar()`: Quality component scores as radar chart
- `render_tm_gc_distribution()`: Tm histogram and Tm vs GC scatter
- `render_coverage_specificity_scatter()`: Coverage vs specificity with Pareto frontier
- `render_primer_heatmap()`: Primer metrics comparison heatmap
- `render_dimer_network_heatmap()`: Primer-dimer interaction matrix
- `render_dimer_network_graph()`: Network graph of strong dimer interactions

All functions return empty string when Plotly is not installed (graceful degradation).

## Testing

```bash
pytest tests/                         # All unit tests
pytest tests/test_genetic_algorithm_integration.py  # Specific test
neoswga validate --quick              # Quick validation
```

**Integration tests** (`tests/integration/`):
- `phi29_baseline/`, `phi29_ga/`: Phi29 polymerase tests
- `equiphi29_baseline/`, `equiphi29_long/`: EquiPhi29 tests

**Validation tests** (`validation_tests/`):
- `test_gc_filter_fix.py`: Adaptive GC filter validation
- `test_complexity_filter.py`: Sequence complexity filtering

## Development Tasks

**Generate k-mer files for non-standard lengths**:
```bash
neoswga count-kmers -j params.json --min-k 15 --max-k 18
```

**Retrain random forest model** (for sklearn updates):
```bash
python scripts/retrain_rf_model.py --output neoswga/core/models/random_forest_filter.p
```

**Run pipeline on test data**:
```bash
cd tests/integration/equiphi29_baseline
neoswga count-kmers -j params.json
neoswga filter -j params.json
neoswga score -j params.json
neoswga optimize -j params.json
```

## Code Patterns

**Parameter handling**:
```python
from neoswga.core.parameter import get_params
params = get_params('params.json')
```

**Multiprocessing**:
```python
from neoswga.core.utility import create_pool
with create_pool(cpus) as pool:
    results = pool.map(process_func, items)
```

**Position data** (HDF5 format):
```python
import h5py
with h5py.File('positions.h5', 'r') as f:
    positions = f[primer_sequence][:]
```

## Known Issues

1. **Large background genomes**: Use `neoswga build-filter` to pre-build a Bloom filter for human genome.

2. **sklearn compatibility**: Random forest model may need retraining after sklearn updates.

3. **Memory usage**: The filter command loads all background k-mers into memory. Use Bloom filter for large backgrounds.

4. **PositionCache strand parameter**: Uses 'forward', 'reverse', 'both' (not '+' or '-').

## Package Structure

```
neoswga/
  __init__.py              # Package init, version
  cli_unified.py           # Main CLI (all commands)
  core/                    # Core functionality (75 modules)
    # Pipeline
    pipeline.py, improved_pipeline.py, multi_genome_pipeline.py

    # Filtering
    filter.py, adaptive_filters.py, kmer.py, kmer_counter.py

    # Scoring
    rf_preprocessing.py, primer_attributes.py

    # Optimization (9 optimizers)
    network_optimizer.py, genetic_algorithm.py, optimize.py,
    milp_optimizer.py, hybrid_optimizer.py, dominating_set_optimizer.py,
    background_aware_optimizer.py, moea_optimizer.py, equiphi29_optimizer.py

    # Thermodynamics
    thermodynamics.py, thermo_estimation.py, reaction_conditions.py,
    secondary_structure.py, dimer.py, three_prime_stability.py

    # Mechanistic Modeling
    mechanistic_params.py, mechanistic_model.py, set_size_optimizer.py,
    model_validation.py, additives.py

    # Performance
    position_cache.py, background_filter.py, gpu_acceleration.py

    # Simulation
    replication_simulator.py, swga_simulator.py, stochastic_simulator.py
    simulation_fitness.py, simulation_analysis.py, simulation_plots.py

    # Analysis
    genome_analysis.py, strand_bias_analyzer.py, dimer_network_analyzer.py

    # Reporting
    report/                # Quality report generation
      metrics.py, quality.py, executive_summary.py,
      technical_report.py, visualizations.py, validation.py, utils.py

    # Utilities
    utility.py, parameter.py, validation.py, genome_io.py

    # Experimental (not fully integrated)
    experimental/
      cooperative_binding_selector.py, hybrid_primer_strategy.py,
      simulate_command.py

    models/                # ML models
      random_forest_filter.p
tests/
  integration/             # Integration test scenarios
  report/                  # Report module tests
  test_*.py                # Unit tests
examples/
  plasmid_example/         # Self-contained example (pcDNA vs pLTR plasmids)
scripts/                   # Development utilities
validation_tests/          # Validation test suite
docs/                      # Documentation
  archive/                 # Historical documents
  validation/              # Validation reports
```

## Running the Example

The `examples/plasmid_example/` provides a quick test with two small plasmids:

```bash
cd examples/plasmid_example
neoswga count-kmers -j params.json
neoswga filter -j params.json
neoswga score -j params.json
neoswga optimize -j params.json
```
