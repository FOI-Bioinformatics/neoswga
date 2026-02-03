# NeoSWGA User Guide

This guide provides comprehensive documentation for using NeoSWGA, including installation, configuration, command-line usage, and Python API examples.

**Validation Note**: Performance improvements described in this guide are based on theoretical analysis using established biophysical models. Experimental validation of comparative claims is ongoing.

## Table of Contents

1. [Installation](#installation)
2. [Core Concepts](#core-concepts)
3. [Command-Line Usage](#command-line-usage)
4. [Python API Reference](#python-api-reference)
5. [Configuration](#configuration)
6. [Advanced Features](#advanced-features)
7. [Troubleshooting](#troubleshooting)

---

## Installation

### System Requirements

- Python >= 3.11
- Jellyfish k-mer counter (v2.0+)
- 8+ GB RAM (16+ GB recommended for large genomes)
- Optional: CUDA-capable GPU for acceleration

### Basic Installation

```bash
# Clone repository
cd /path/to/swga2  # Note: directory name differs from package name

# Create virtual environment
python3 -m venv neoswga_env
source neoswga_env/bin/activate  # Windows: neoswga_env\Scripts\activate

# Install package
pip install -e .

# Verify installation
python test_install.py
```

### Optional Dependencies

**GPU Acceleration** (requires CUDA):
```bash
# For CUDA 11.x
pip install -e ".[gpu]"

# For CUDA 12.x
pip install -e ".[gpu-cuda12]"
```

**Deep Learning** (optional):
```bash
# PyTorch backend
pip install -e ".[deep-learning-torch]"

# TensorFlow backend
pip install -e ".[deep-learning-tf]"
```

**All Features**:
```bash
pip install -e ".[all]"
```

### Installing Jellyfish

Jellyfish is required for k-mer counting. Install from:
https://www.cbcb.umd.edu/software/jellyfish/

Verify installation:
```bash
jellyfish --version
```

---

## Core Concepts

### Thermodynamic Modeling

NeoSWGA implements the SantaLucia nearest-neighbor model for DNA hybridization thermodynamics:

**Key parameters**:
- Nearest-neighbor dinucleotide stacking energies (16 parameters)
- Temperature-dependent free energy: ΔG = ΔH - TΔS
- Salt correction using Owczarzy unified formula (Na+, Mg2+)
- Helix initiation and terminal penalties

**Example**:
```python
from neoswga.core import thermodynamics as thermo

seq = 'ATCGATCGATCG'
tm = thermo.calculate_tm_with_salt(seq, na_conc=50, mg_conc=0)
dg = thermo.calculate_free_energy(seq, temp=37)
gc = thermo.gc_content(seq)

print(f"Tm: {tm:.1f}°C")
print(f"ΔG at 37°C: {dg:.2f} kcal/mol")
print(f"GC content: {gc:.1%}")
```

**Estimated accuracy**: Based on SantaLucia parameters; validation against experimental melting curves pending.

### Reaction Conditions

Reaction conditions affect primer binding and extension:

**Polymerases**:
- **Phi29**: Optimal 30-40°C, standard enzyme
- **EquiPhi29**: Optimal 42-45°C, thermostable variant

**Additives**:
- **DMSO**: Reduces Tm (~0.6°C per %)
- **Betaine**: Equalizes AT/GC stability (~2.3°C per M)
- **Trehalose**: Enzyme stabilization (~5°C per M)

**Example**:
```python
from neoswga.core import reaction_conditions as rc

# Standard conditions
std = rc.get_standard_conditions()
print(f"Max primer length: {std.max_primer_length()}")  # 12

# Enhanced conditions (enables longer primers)
enhanced = rc.get_enhanced_conditions()
print(f"Max primer length: {enhanced.max_primer_length()}")  # 15

# Custom conditions
custom = rc.ReactionConditions(
    temp=42.0,
    dmso_percent=5.0,
    betaine_m=1.0,
    polymerase='equiphi29'
)
tm_eff = custom.calculate_effective_tm('ATCGATCGATCGATC')
print(f"Effective Tm: {tm_eff:.1f}°C")
```

**Note**: Additive effects are based on literature values; specific effects may vary by formulation.

### Secondary Structure Prediction

Detection of primer dimers and hairpins using dynamic programming:

**Methods**:
- Heterodimer detection (primer-primer interactions)
- Homodimer detection (self-complementarity)
- Hairpin prediction (stem-loop structures)
- Severity scoring (0-1 scale)

**Example**:
```python
from neoswga.core import secondary_structure as ss
from neoswga.core import reaction_conditions as rc

conditions = rc.get_standard_conditions()

# Check for primer dimers
dimer = ss.check_heterodimer('ATCGATCGAT', 'ATCGATCGAT', conditions)
print(f"Dimer ΔG: {dimer['energy']:.1f} kcal/mol")
print(f"Severity: {dimer['severity']:.2f}")
print(f"Forms dimer: {dimer['forms_dimer']}")

# Check for hairpins
hairpins = ss.check_hairpins('GCGATCGCAAAAGCGATCGC')
if hairpins:
    hp = hairpins[0]
    print(f"Hairpin: {hp['stem_length']}bp stem, {hp['loop_size']}nt loop")
    print(f"Stable: {hp['stable']}")
```

**Estimated sensitivity**: Dynamic programming approach may detect additional dimer configurations compared to simple string matching; comprehensive validation pending.

---

## Command-Line Usage

### Quick Design

For rapid primer design with default parameters:

```bash
neoswga quick-design --fg target_genome.fasta --output results/
```

### Pipeline Workflow

#### count-kmers: K-mer Preprocessing

Identify all k-mers (6-12 or 6-18 mers) in target and off-target genomes:

```bash
neoswga count-kmers -j params.json
```

**Output**: `*_Xmer_all.txt` files for k=6 to 12 (or custom range via `--min-k`, `--max-k`)

#### filter: Candidate Filtering

Filter primers based on thermodynamic and binding properties:

```bash
neoswga filter -j params.json
```

**Filtering criteria**:
- Melting temperature range
- Binding frequency (normalized)
- Evenness (Gini index)
- Self-complementarity
- Adaptive GC content

**Output**: `step2_df.csv` and HDF5 files with primer positions

#### score: Amplification Scoring

Predict amplification potential using machine learning:

```bash
neoswga score -j params.json
```

**Output**: `step3_df.csv` with ranked primers by predicted amplification score

#### optimize: Primer Set Optimization

Select optimal primer combinations:

```bash
neoswga optimize -j params.json
```

**Methods**:
- Network-based optimization (default, recommended)
- Greedy breadth-first search
- Genetic algorithm
- MILP (mixed-integer linear programming)

**Output**: `step4_improved_df.csv` with top primer sets and evaluation scores

### Common Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `-j, --json_file` | Configuration file path | None |
| `--min_fg_freq` | Min frequency in target | 1e-5 |
| `--max_bg_freq` | Max frequency in off-target | 5e-5 |
| `--max_gini` | Max Gini index | 0.6 |
| `--min_tm` | Min melting temperature | 15°C |
| `--max_tm` | Max melting temperature | 45°C |
| `--max_dimer_bp` | Max complementary base pairs | 4 |
| `--cpus` | Number of CPU cores | All |

---

## Python API Reference

### Thermodynamics Module

```python
from neoswga.core import thermodynamics as thermo

# Basic calculations
tm = thermo.calculate_tm_with_salt(seq, na_conc=50, mg_conc=0)
dg = thermo.calculate_free_energy(seq, temp=37)
gc = thermo.gc_content(seq)
prob = thermo.calculate_binding_probability(seq, temp=30)

# Returns
# tm: float (°C)
# dg: float (kcal/mol)
# gc: float (0.0-1.0)
# prob: float (0.0-1.0)
```

### Reaction Conditions Module

```python
from neoswga.core import reaction_conditions as rc

# Preset conditions
std = rc.get_standard_conditions()          # Phi29, 37°C
equiphi = rc.get_equiphi_conditions()       # EquiPhi29, 42°C
enhanced = rc.get_enhanced_conditions()     # With additives
high_gc = rc.get_high_gc_conditions()       # For GC-rich genomes
low_temp = rc.get_low_temp_conditions()     # With SSB proteins

# Custom conditions
conditions = rc.ReactionConditions(
    temp=42.0,
    na_conc=50.0,           # mM
    mg_conc=0.0,            # mM
    dmso_percent=5.0,       # %
    betaine_m=1.0,          # M
    trehalose_m=0.0,        # M
    polymerase='equiphi29',
    has_ssb=False
)

# Query properties
max_len = conditions.max_primer_length()
tm_range = conditions.get_polymerase_range()
tm_eff = conditions.calculate_effective_tm(seq)
```

### Secondary Structure Module

```python
from neoswga.core import secondary_structure as ss

# Dimer checking
heterodimer = ss.check_heterodimer(seq1, seq2, conditions)
homodimer = ss.check_homodimer(seq, conditions)

# Returns dict with:
# - energy: float (kcal/mol)
# - tm: float (°C)
# - severity: float (0-1)
# - forms_dimer: bool
# - binding_length: int

# Hairpin detection
hairpins = ss.check_hairpins(seq)
# Returns list of hairpin dicts

# Batch filtering
filtered = ss.filter_primers_by_structure(
    primers,
    max_hairpin_tm=35.0,
    max_self_dimer_severity=0.3,
    conditions=conditions
)
```

### Adaptive Search Module

```python
from neoswga.core import adaptive_search
from neoswga.core import reaction_conditions as rc

conditions = rc.get_enhanced_conditions()

result = adaptive_search.adaptive_search_pipeline(
    fg_genomes=['target.fasta'],
    bg_genomes=['offtarget.fasta'],
    fg_prefixes=['target_kmers'],
    bg_prefixes=['offtarget_kmers'],
    conditions=conditions,
    min_primers=100,
    max_primers=1000
)

# Access results
optimal_k = result.optimal_k           # Selected k-mer length
primers = result.primers               # List of primer sequences
metrics = result.metrics               # Performance metrics
```

### Genetic Algorithm Module

```python
from neoswga.core import genetic_algorithm as ga

config = ga.GAConfig(
    population_size=200,
    generations=100,
    mutation_rate=0.15,
    crossover_rate=0.8,
    elitism_fraction=0.1
)

best_set = ga.optimize_primer_set_ga(
    primer_pool=candidate_primers,
    fg_prefixes=['target_kmers'],
    bg_prefixes=['offtarget_kmers'],
    fg_lengths=[4400000],
    bg_lengths=[3000000000],
    conditions=conditions,
    config=config
)

print(f"Fitness: {best_set.fitness:.4f}")
print(f"Primers: {best_set.primers}")
```

### GPU Acceleration

```python
from neoswga.core import gpu_acceleration as gpu

# Check availability
if gpu.GPU_AVAILABLE:
    print("GPU acceleration active")
    print(f"CuPy version: {gpu.cp.__version__}")
else:
    print("Using NumPy fallback")

# GPU operations are automatic when available
# No API changes required
```

### Deep Learning Module

```python
from neoswga.core import deep_learning as dl

# Check backend
if dl.DL_AVAILABLE:
    print(f"Using {dl.DL_BACKEND} backend")
else:
    print("Using fallback embeddings")

# Enhanced scoring (when available)
# Integrated into amplification scoring step
```

---

## Configuration

### JSON Configuration File

Create `params.json` for workflow configuration:

```json
{
  "fg_genomes": ["target.fasta"],
  "bg_genomes": ["offtarget.fasta"],
  "fg_prefixes": ["target_kmers"],
  "bg_prefixes": ["offtarget_kmers"],
  "data_dir": "./project/",

  "min_fg_freq": 1e-05,
  "max_bg_freq": 5e-05,
  "max_gini": 0.6,
  "max_primer": 500,
  "min_tm": 15,
  "max_tm": 45,
  "max_self_dimer_bp": 4,

  "min_amp_pred": 5,
  "max_dimer_bp": 4,
  "iterations": 8,
  "max_sets": 5,
  "top_sets_count": 10,

  "polymerase": "equiphi29",
  "dmso_percent": 5.0,
  "betaine_m": 1.0,
  "cpus": 8
}
```

---

## Advanced Features

### Network Analysis

Analyze primer set coverage using graph-based methods:

```python
from neoswga.core import amplicon_network

network = amplicon_network.AmpliconNetwork(
    primers=['ATCGAT', 'GCTAGC', 'TTAACC'],
    genome_length=4400000,
    max_amplicon_length=5000
)

network.load_positions_from_hdf5('target_kmers')
network.build_network()

metrics = network.analyze()
print(f"Coverage: {metrics.coverage_fraction:.1%}")
print(f"Number of amplicons: {metrics.num_amplicons}")
print(f"Hub primers: {metrics.num_hubs}")

# Identify critical primers
critical = network.find_critical_primers()
```

**Note**: Network metrics provide theoretical coverage estimates; experimental validation required.

### Replication Simulation

Agent-based simulation of phi29 replication:

```python
from neoswga.core import replication_simulator

positions = {
    'ATCGAT': {'forward': [100, 500], 'reverse': [2000]},
    'GCTAGC': {'forward': [300], 'reverse': [1500]}
}

result = replication_simulator.simulate_primer_set(
    primers=['ATCGAT', 'GCTAGC'],
    primer_positions=positions,
    genome_length=4400000,
    genome_sequence=genome_seq,
    conditions=conditions,
    duration=3600,  # seconds
    n_replicates=5
)

print(f"Mean coverage: {result['mean_coverage']:.1%}")
print(f"Std coverage: {result['std_coverage']:.1%}")
```

**Parameters**:
- Extension rate: 167 bp/s (Phi29), 200 bp/s (EquiPhi29)
- GC penalty: Slower extension in high-GC regions
- Stochastic primer binding based on thermodynamics

**Note**: Simulation results are theoretical predictions based on phi29 kinetics; comparison with experimental data pending.

---

## Troubleshooting

### Import Errors

**Problem**: `ModuleNotFoundError: No module named 'neoswga'`

**Solution**: Ensure package is installed:
```bash
pip install -e .
python -c "import neoswga; print(neoswga.__version__)"
```

### Melting Temperature Issues

**Problem**: Calculated Tm values seem incorrect

**Possible causes**:
- Salt concentration not specified (default: 50 mM Na+)
- Sequence too short (< 6 bases)
- Non-standard bases in sequence

**Check**:
```python
from neoswga.core import thermodynamics as thermo

seq = 'ATCGATCG'
tm_default = thermo.calculate_tm_with_salt(seq)  # 50 mM Na+
tm_100 = thermo.calculate_tm_with_salt(seq, na_conc=100)  # 100 mM Na+
print(f"Default: {tm_default:.1f}°C, 100mM: {tm_100:.1f}°C")
```

### GPU Not Detected

**Problem**: GPU acceleration not working despite CUDA installation

**Check**:
```python
from neoswga.core import gpu_acceleration as gpu
print(f"GPU available: {gpu.GPU_AVAILABLE}")

if not gpu.GPU_AVAILABLE:
    print("Install CuPy: pip install cupy-cuda11x")
```

**Note**: CuPy version must match CUDA version.

### Memory Issues

**Problem**: Out of memory errors with large genomes

**Solutions**:
- Reduce `max_primer` parameter (e.g., 500 → 200)
- Use disk-based processing for k-mer files
- Increase system RAM
- Process chromosomes separately

### Jellyfish Not Found

**Problem**: `jellyfish: command not found`

**Solution**: Ensure Jellyfish is installed and in PATH:
```bash
# Check installation
which jellyfish

# If not found, install from source or package manager
# https://www.cbcb.umd.edu/software/jellyfish/
```

---

## Performance Considerations

### Dataset Size Guidelines

| Genome Size | RAM Required | Time (Approx) |
|-------------|--------------|---------------|
| < 5 Mb | 4 GB | Minutes |
| 5-50 Mb | 8 GB | 10-30 min |
| 50-500 Mb | 16 GB | 1-3 hours |
| > 500 Mb | 32+ GB | Hours |

**Note**: Times are estimates and vary by CPU speed and primer pool size.

### Optimization Tips

1. **Use adaptive search**: Automatically selects optimal k-mer length
2. **Limit primer pool**: Use `max_primer` parameter to control memory usage
3. **Parallelize**: Use `--cpus` parameter for multi-core processing
4. **GPU acceleration**: Install CuPy for large-scale calculations
5. **Filter early**: Strict filtering in step 2 reduces downstream computation

### Known Performance Limitations

- Dimer matrix computation: O(n²) for n primers
- Network analysis: Memory-intensive for > 1000 primers
- Simulation: ~1 minute per replicate per simulated hour
- GPU speedup varies: 2-10x typical, depends on dataset size

---

## Validation Status

| Feature | Implementation | Validation Status |
|---------|---------------|-------------------|
| Thermodynamics (SantaLucia) | Complete | Literature-based, experimental validation pending |
| Salt correction (Owczarzy) | Complete | Literature-based |
| Secondary structure | Complete | Algorithm-based, sensitivity estimation pending |
| Adaptive k-mer search | Complete | Theoretical, comparative testing required |
| Genetic algorithm | Complete | Functional, performance comparison pending |
| Network analysis | Complete | Theoretical coverage prediction |
| Replication simulation | Complete | Physics-based, experimental comparison required |
| GPU acceleration | Complete | Functional, benchmark data pending |
| Deep learning | Complete | Functional, validation ongoing |

**Summary**: Core features are implemented based on established biophysical models. Comparative performance claims require experimental validation against original SOAPswga and laboratory results.

---

## Getting Help

- Review this guide and [README.md](README.md)
- Check [MIGRATION.md](MIGRATION.md) for API changes
- See [DEVELOPMENT.md](DEVELOPMENT.md) for implementation details
- Report issues via GitHub issue tracker

---

## References

**Biophysical Models**:
- SantaLucia, J. (1998). A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics. PNAS, 95(4), 1460-1465.
- Owczarzy, R., et al. (2008). Predicting stability of DNA duplexes in solutions containing magnesium and monovalent cations. Biochemistry, 47(19), 5336-5353.

**Original Method**:
- Dwivedi-Yu, J. A., et al. (2023). A fast machine-learning-guided primer design pipeline for selective whole genome amplification. PLOS Computational Biology, 19(4), e1010137.

**Algorithms**:
- Zuker, M. (2003). Mfold web server for nucleic acid folding and hybridization prediction. Nucleic Acids Research, 31(13), 3406-3415.
