# NeoSWGA Development Documentation

This document provides technical details for developers working on NeoSWGA, including architecture overview, implementation details, packaging information, and contribution guidelines.

## Table of Contents

1. [Project Overview](#project-overview)
2. [Architecture](#architecture)
3. [Module Descriptions](#module-descriptions)
4. [Implementation Status](#implementation-status)
5. [Testing](#testing)
6. [Packaging and Distribution](#packaging-and-distribution)
7. [Contributing](#contributing)
8. [Future Work](#future-work)

---

## Project Overview

### History and Motivation

NeoSWGA extends the original SOAPswga (Dwivedi-Yu et al., 2023) with enhanced thermodynamic modeling, additional optimization algorithms, and optional GPU/deep learning support. The project aims to implement biophysically accurate primer design while maintaining computational efficiency.

**Design Goals**:
- Implement rigorous thermodynamic calculations based on established models
- Provide alternative optimization strategies (genetic algorithms, network analysis)
- Enable longer primer design (13-15mers) through reaction condition modeling
- Maintain compatibility with original SOAPswga workflow

**Development Status**: Beta - Core features implemented, experimental validation ongoing.

### Key Enhancements from SOAPswga

1. **Thermodynamics**: Full SantaLucia nearest-neighbor model with Owczarzy salt corrections
2. **Reaction Modeling**: Support for thermodynamic additives (DMSO, betaine) and multiple polymerases
3. **Secondary Structure**: Dynamic programming algorithms for dimer/hairpin detection
4. **Optimization**: Genetic algorithms as alternative to greedy search
5. **Analysis Tools**: Network-based coverage prediction and replication simulation
6. **Performance**: Optional GPU acceleration and deep learning integration

**Validation Note**: Comparative performance claims (e.g., accuracy improvements, sensitivity gains) are theoretical estimates based on algorithmic analysis. Experimental validation against original SOAPswga and laboratory data is ongoing.

---

## Architecture

### Package Structure

```
swga2/  (repository directory)
├── neoswga/  (Python package)
│   ├── __init__.py
│   ├── cli.py  (unified command-line interface)
│   ├── logging_config.py
│   └── core/  (core modules)
│       ├── __init__.py
│       ├── thermodynamics.py  (459 lines)
│       ├── reaction_conditions.py  (415 lines)
│       ├── secondary_structure.py  (548 lines)
│       ├── adaptive_search.py  (540 lines)
│       ├── genetic_algorithm.py  (480 lines)
│       ├── amplicon_network.py  (430 lines)
│       ├── replication_simulator.py  (470 lines)
│       ├── advanced_features.py
│       ├── gpu_acceleration.py
│       ├── deep_learning.py
│       ├── unified_pipeline.py
│       └── [legacy modules]
├── test_install.py
├── pyproject.toml
├── setup.py
└── README.md
```

### Design Patterns

**Module Organization**:
- Core scientific modules in `neoswga.core/`
- Legacy SOAPswga modules preserved for compatibility
- Optional dependencies (GPU, deep learning) import conditionally

**Configuration**:
- Preset reaction conditions for common scenarios
- JSON-based parameter files for workflow configuration
- Backward compatibility with SOAPswga parameter structure

**Error Handling**:
- Graceful fallback for optional dependencies (CuPy, PyTorch, TensorFlow)
- Warning messages for missing optional features
- Validation of input parameters

---

## Module Descriptions

### thermodynamics.py (459 lines)

**Purpose**: DNA hybridization thermodynamics calculations

**Key Functions**:
- `calculate_tm_with_salt(seq, na_conc, mg_conc)`: Melting temperature with salt correction
- `calculate_free_energy(seq, temp)`: ΔG at specified temperature
- `calculate_binding_probability(seq, temp)`: Boltzmann-weighted binding probability
- `gc_content(seq)`: GC fraction

**Implementation**:
- SantaLucia (1998) nearest-neighbor parameters
- 16 dinucleotide stacking energies (ΔH, ΔS)
- Owczarzy et al. (2008) unified salt correction formula
- Helix initiation penalties
- Symmetry corrections for self-complementary sequences

**Dependencies**: NumPy

**Validation Status**: Parameters from literature; accuracy comparison with experimental melting curves pending.

### reaction_conditions.py (415 lines)

**Purpose**: Model reaction conditions affecting primer behavior

**Key Classes**:
- `ReactionConditions`: Container for temperature, salt, additives, polymerase type
- Preset factory functions: `get_standard_conditions()`, `get_enhanced_conditions()`, etc.

**Features**:
- Polymerase profiles (Phi29, EquiPhi29)
- Thermodynamic additive effects (DMSO, betaine, trehalose, SSB)
- Effective Tm calculation with additive corrections
- Maximum primer length based on conditions

**Additive Effects**:
- DMSO: -0.6°C per % (Jones et al., 2001)
- Betaine: -2.3°C per M (Rees et al., 1993)
- Trehalose: -5°C per M (estimated)
- SSB proteins: -5°C effective (estimated)

**Dependencies**: thermodynamics module

**Validation Status**: Additive effects from literature; specific formulation effects may vary.

### secondary_structure.py (548 lines)

**Purpose**: Predict secondary structures (dimers, hairpins)

**Key Functions**:
- `check_heterodimer(seq1, seq2, conditions)`: Primer-primer interaction
- `check_homodimer(seq, conditions)`: Self-complementarity
- `check_hairpins(seq)`: Stem-loop structures
- `filter_primers_by_structure(primers, conditions)`: Batch filtering

**Algorithm**:
- Dynamic programming (similar to Zuker algorithm)
- Scoring: Watson-Crick base pairing, mismatches, bulges, loops
- Severity scoring: 0-1 scale based on ΔG and Tm
- 3' end extension risk assessment

**Dependencies**: thermodynamics, reaction_conditions, NumPy

**Validation Status**: Dynamic programming approach; comprehensive sensitivity testing against experimental dimer data pending.

### adaptive_search.py (540 lines)

**Purpose**: Automatic selection of optimal k-mer length

**Algorithm**:
1. Test k-mer lengths from 6 to max (12 or 15 based on conditions)
2. Filter by thermodynamic properties (Tm range)
3. Filter by secondary structure (hairpins, self-dimers)
4. Filter by specificity (fg/bg ratio)
5. Rank by composite score
6. Select k with highest score and sufficient candidates

**Composite Score** (weighted):
- Specificity (40%): log10(fg/bg ratio)
- Coverage (25%): Sigmoid of binding site count
- Tm optimality (20%): Gaussian around polymerase midpoint
- GC balance (10%): Penalty for extreme GC
- Free energy (5%): Binding favorability

**Dependencies**: Core modules, HDF5, multiprocessing

**Validation Status**: Theoretical scoring function; performance comparison with fixed k-mer lengths requires experimental testing.

### genetic_algorithm.py (480 lines)

**Purpose**: Evolutionary optimization for primer set selection

**Algorithm**: Standard genetic algorithm with domain-specific operators

**Components**:
- Population: 200 individuals (configurable)
- Selection: Tournament (k=5)
- Crossover: Uniform with dimer checking
- Mutation: Add/remove/replace primers
- Elitism: Preserve top 10%

**Fitness Function**:
```python
fitness = (
    0.35 * coverage_score +
    0.30 * specificity_score +
    0.15 * evenness_score +
    0.10 * thermodynamic_score +
    0.10 * dimer_avoidance_score
)
```

**Dependencies**: Core modules, multiprocessing

**Validation Status**: Functional implementation; comparative performance vs greedy search requires benchmarking.

### amplicon_network.py (430 lines)

**Purpose**: Graph-based analysis of amplification coverage

**Model**:
- Nodes: Primer binding sites (position, strand, primer_id)
- Edges: Potential amplicons (forward site → reverse site, length ≤ max)

**Metrics**:
- Coverage fraction: Connected genome regions
- Hubs: High-degree nodes (amplification hotspots)
- Betweenness centrality: Critical primers
- Network density: Connectivity measure

**Dependencies**: NetworkX, NumPy, HDF5

**Validation Status**: Theoretical coverage prediction; correlation with experimental coverage requires validation.

### replication_simulator.py (470 lines)

**Purpose**: Agent-based simulation of phi29 replication

**Model Components**:
- Agents: Individual replication forks
- States: ACTIVE, PAUSED, TERMINATED, COLLIDED
- Kinetics: Temperature-dependent primer binding (Boltzmann)
- Extension: Polymerase-specific rates with GC penalty
- Collisions: Fork-fork interactions

**Parameters**:
- Extension rate: 167 bp/s (Phi29), 200 bp/s (EquiPhi29)
- Time step: 1 second
- Duration: 3600 seconds (1 hour) typical
- GC penalty: 50% slower in 70% GC regions

**Dependencies**: Core modules, NumPy

**Validation Status**: Physics-based simulation; comparison with time-course experimental data pending.

### gpu_acceleration.py

**Purpose**: Optional GPU acceleration via CuPy

**Implementation**:
- Conditional import of CuPy
- Fallback to NumPy if GPU unavailable
- Transparent GPU usage for compatible operations

**Dependencies**: CuPy (optional), NumPy

**Performance**: Speedup varies by dataset size and hardware (2-10x typical for large datasets).

### deep_learning.py

**Purpose**: Optional deep learning for primer scoring

**Implementation**:
- Conditional import of PyTorch or TensorFlow
- Fallback embeddings if DL unavailable
- Transformer-based primer embeddings (when available)

**Dependencies**: PyTorch or TensorFlow (optional)

**Validation Status**: Functional implementation; comparative scoring accuracy requires validation.

---

## Implementation Status

### Completed Features

| Module | Status | Lines | Dependencies |
|--------|--------|-------|--------------|
| thermodynamics | Complete | 459 | NumPy |
| reaction_conditions | Complete | 415 | thermodynamics |
| secondary_structure | Complete | 548 | thermodynamics, NumPy |
| adaptive_search | Complete | 540 | Core modules, HDF5 |
| genetic_algorithm | Complete | 480 | Core modules |
| amplicon_network | Complete | 430 | NetworkX |
| replication_simulator | Complete | 470 | Core modules |
| gpu_acceleration | Complete | - | CuPy (optional) |
| deep_learning | Complete | - | PyTorch/TF (optional) |
| unified_pipeline | Complete | - | Core modules |

**Total**: ~3,800+ lines of enhanced code

### Known Issues and Limitations

1. **Testing**: Test suite is minimal; only `test_install.py` exists
2. **Validation**: Performance claims require experimental verification
3. **Documentation**: API docstrings need expansion
4. **Benchmarking**: No systematic performance comparisons
5. **Edge Cases**: Limited testing with extreme parameter values
6. **GPU**: Performance varies significantly by hardware and dataset size
7. **Deep Learning**: Model training and validation data not included

---

## Testing

### Current Test Coverage

**Existing Tests**:
- `test_install.py`: Basic import and functionality verification
- Manual testing during development

**Test Status**: Minimal

### Needed Tests

**Unit Tests** (priority):
- `tests/test_thermodynamics.py`: Validate against literature values
- `tests/test_secondary_structure.py`: Test dimer/hairpin detection
- `tests/test_reaction_conditions.py`: Verify additive effects
- `tests/test_adaptive_search.py`: Test k-mer selection logic
- `tests/test_genetic_algorithm.py`: Validate fitness calculations

**Integration Tests**:
- Full pipeline workflow (count-kmers → filter → score → optimize)
- API examples from documentation
- Command-line interface

**Validation Tests**:
- Compare Tm predictions with experimental melting curves
- Compare dimer predictions with gel electrophoresis
- Compare primer sets with experimental amplification data

### Running Tests

```bash
# Current minimal test
python test_install.py

# Future comprehensive testing (when implemented)
pytest tests/ -v --cov=neoswga --cov-report=html
```

---

## Packaging and Distribution

### Package Configuration

**Files**:
- `pyproject.toml`: Modern packaging (PEP 517/518)
- `setup.py`: Backward compatibility
- `MANIFEST.in`: Package data specification
- `VERSION`: Version file (3.0.0)

**Package Name**: `neoswga` (note: directory is `swga2`)

### Building Distribution

```bash
# Clean previous builds
make clean

# Build distribution packages
python -m build

# Output: dist/neoswga-3.0.0.tar.gz, dist/neoswga-3.0.0-py3-none-any.whl
```

### Version Management

**Current Version**: 3.0.0

**Versioning Scheme**: Semantic Versioning (MAJOR.MINOR.PATCH)
- MAJOR: Breaking API changes (current: 3)
- MINOR: New features, backward compatible
- PATCH: Bug fixes

**Version Locations**:
- `VERSION` file
- `neoswga/__init__.py`: `__version__ = "3.0.0"`
- `pyproject.toml`: `version = "3.0.0"`
- `setup.py`: `version = "3.0.0"`

### Dependencies

**Core Requirements**:
- numpy >= 1.18.0
- pandas >= 1.0.0
- scipy >= 1.4.1
- scikit-learn >= 0.21.3
- h5py >= 2.10.0
- networkx >= 2.5
- biopython >= 1.78
- joblib >= 0.14.1
- matplotlib >= 3.3.0
- seaborn >= 0.11.0
- melt >= 1.0.3 (legacy compatibility)
- openpyxl >= 3.0.0

**Optional Dependencies**:
- GPU: `cupy-cuda11x` or `cupy-cuda12x`
- Deep Learning: `torch` or `tensorflow`
- Development: `pytest`, `black`, `mypy`, etc.

### Installation Methods

```bash
# Editable install (development)
pip install -e .

# With optional features
pip install -e ".[gpu]"
pip install -e ".[deep-learning-torch]"
pip install -e ".[all]"

# Development environment
make init-dev  # Installs dev dependencies + pre-commit hooks
```

---

## Contributing

### Development Setup

```bash
# Clone repository
cd /path/to/swga2

# Create development environment
python3 -m venv dev_env
source dev_env/bin/activate

# Install with dev dependencies
pip install -e ".[dev]"

# Install pre-commit hooks (when configured)
pre-commit install
```

### Code Style

**Python Style**:
- PEP 8 compliance
- Black formatter (line length: 100)
- isort for imports
- Type hints where appropriate

**Documentation**:
- Docstrings for all public functions
- NumPy/SciPy docstring format
- Examples in docstrings

**Commands**:
```bash
# Format code
make format

# Lint
make lint

# Type check
make type-check
```

### Contribution Guidelines

1. **Issues**: Report bugs or feature requests via GitHub issues
2. **Pull Requests**: Fork, create feature branch, submit PR
3. **Testing**: Add tests for new features
4. **Documentation**: Update relevant documentation files
5. **Validation**: Include validation data where possible

### Priority Contributions

**High Priority**:
- Comprehensive test suite
- Experimental validation data
- Benchmark comparisons with SOAPswga
- API documentation expansion

**Medium Priority**:
- Additional polymerase profiles
- Improved error messages
- Performance optimizations
- Tutorial examples with real data

**Low Priority**:
- GUI interface
- Additional export formats
- Visualization tools

---

## Future Work

### Validation Studies

**Needed Validation**:
1. **Thermodynamics**: Compare predicted vs measured Tm (qPCR melting curves)
2. **Secondary Structure**: Compare predictions vs gel electrophoresis
3. **Primer Sets**: Compare coverage vs Illumina sequencing
4. **Performance**: Benchmark against original SOAPswga

**Experimental Collaboration**: Validation requires wet-lab experiments or collaboration with groups performing SWGA.

### Algorithm Enhancements

**Potential Improvements**:
- Multi-objective optimization (Pareto front)
- Machine learning for primer set scoring
- Improved dimer detection (consider mismatches)
- Temperature-dependent kinetics in simulation
- Support for additional polymerases (Bst, Klenow)

### Performance Optimization

**Optimization Opportunities**:
- Numba JIT compilation for hot loops
- Cython extensions for critical functions
- Memory-mapped HDF5 for large genomes
- Distributed computing for large-scale screening

### Documentation

**Needed Documentation**:
- Comprehensive API reference (auto-generated)
- Step-by-step tutorials with example data
- Performance tuning guide
- Algorithm descriptions and derivations
- Benchmark results (when available)

### Community

**Building Community**:
- Publish validation studies
- Present at conferences (e.g., RECOMB, ISMB)
- Create example datasets repository
- Establish user forum or discussion board

---

## References

### Scientific Basis

**Thermodynamics**:
- SantaLucia, J. (1998). A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics. PNAS, 95(4), 1460-1465.
- Owczarzy, R., et al. (2008). Predicting stability of DNA duplexes in solutions containing magnesium and monovalent cations. Biochemistry, 47(19), 5336-5353.

**Additives**:
- Rees, W. A., et al. (1993). Betaine can eliminate the base pair composition dependence of DNA melting. Biochemistry, 32(1), 137-144.
- Henke, W., et al. (1997). Betaine improves the PCR amplification of GC-rich DNA sequences. Nucleic Acids Research, 25(19), 3957-3958.

**Algorithms**:
- Zuker, M. (2003). Mfold web server for nucleic acid folding and hybridization prediction. Nucleic Acids Research, 31(13), 3406-3415.
- Goldberg, D. E. (1989). Genetic Algorithms in Search, Optimization and Machine Learning.

**Original Method**:
- Dwivedi-Yu, J. A., et al. (2023). A fast machine-learning-guided primer design pipeline for selective whole genome amplification. PLOS Computational Biology, 19(4), e1010137.

---

## Contact and Support

**Maintainer**: Andreas Sjodin

**Issues**: Use GitHub issue tracker for bug reports and feature requests

**Contributions**: See [Contributing](#contributing) section above

**Original SOAPswga**: Jane Dwivedi-Yu and collaborators

---

## License

MIT License - See LICENSE file for details

---

## Changelog

### Version 3.0.0 (Current)

**Major Changes**:
- Complete rebrand from swga2 to NeoSWGA
- Restructured package: `src/` → `neoswga.core/`
- Single CLI command: `neoswga`
- Python 3.11+ requirement

**Features**:
- Enhanced thermodynamic modeling
- Reaction condition support
- Genetic algorithm optimization
- Network analysis tools
- Optional GPU/deep learning

**Status**: Beta - Core features implemented, validation ongoing

### Version 2.0.0 (Previous)

- Original enhanced implementation
- Multiple CLI commands
- Python 3.8+ support

### Version 1.0.0 (SOAPswga)

- Original SOAPswga by Dwivedi-Yu et al. (2023)
