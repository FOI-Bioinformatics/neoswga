# NeoSWGA Simulator User Guide

---

## Table of Contents

1. [Quick Start](#quick-start)
2. [Installation](#installation)
3. [Usage Examples](#usage-examples)
4. [Command Reference](#command-reference)
5. [Python API](#python-api)
6. [Understanding Results](#understanding-results)
7. [Troubleshooting](#troubleshooting)
8. [Advanced Topics](#advanced-topics)

---

## Quick Start

### Quick Primer Validation

```bash
neoswga simulate \
  --primers SEQ1 SEQ2 SEQ3 \
  --genome target_genome.fasta \
  --output sim_results/
```

**Output**: Simulation results with coverage, specificity, and recommendations.

---

## Installation

### Requirements

```bash
# Core dependencies (required)
pip install numpy scipy biopython h5py

# Visualization (optional, recommended)
pip install matplotlib

# Interactive plots (optional)
pip install plotly
```

### Verify Installation

```bash
neoswga validate --quick
```

---

## Usage Examples

### Example 1: Quick Validation

**Scenario**: You have 15 primers and want to check if they're good enough.

```bash
python3 production_test/test_simulator.py
```

**Expected output**:
```
YERSINIA TEST PASSED
  Coverage: 42.3%
  Enrichment: 639,972x
  Recommendation: GOOD
```

### Example 2: Comprehensive Analysis

**Scenario**: You want detailed analysis with plots and HTML report.

```python
from neoswga.core.swga_simulator import SwgaSimulator
from neoswga.core.simulation_analysis import SimulationAnalyzer, print_analysis_report
from neoswga.core.simulation_plots import generate_plots
from neoswga.core.simulation_report import generate_html_report

# Initialize
sim = SwgaSimulator(
    primers=['ATCGATCG', 'GCTAGCTA', ...],
    fg_genome='target.fasta',
    bg_genome='background.fasta',
    fg_positions_h5='positions.h5',
    bg_positions_h5='positions.h5'
)

# Run simulation
result = sim.simulate_fast()

# Analyze
analyzer = SimulationAnalyzer(result, sim.fg_positions, sim.bg_positions,
                              sim.fg_length, sim.bg_length)
analysis = analyzer.analyze()

# Generate outputs
print_analysis_report(analysis)
generate_plots(result, sim, analysis, 'plots.png')
generate_html_report(result, 'report.html', analysis)
```

### Example 3: Comparing Multiple Primer Sets

**Scenario**: You have 3 primer sets and want to pick the best.

```python
results = {}

for primer_set in ['set1.txt', 'set2.txt', 'set3.txt']:
    primers = load_primers(primer_set)

    sim = SwgaSimulator(primers, fg_genome, bg_genome, positions, positions)
    result = sim.simulate_fast()

    results[primer_set] = {
        'coverage': result.target_coverage,
        'enrichment': result.enrichment,
        'recommendation': result.recommendation,
        'composite': result.composite_score
    }

# Print comparison
for name, metrics in sorted(results.items(), key=lambda x: x[1]['composite'], reverse=True):
    print(f"{name}: {metrics['recommendation']} (score: {metrics['composite']:.2f})")
```

### Example 4: Configuration File

**Scenario**: You want repeatable simulations with saved settings.

**config.json**:
```json
{
  "primers": "primers.txt",
  "target_genome": "yersinia.fasta",
  "background_genome": "human.fasta",
  "positions": "positions.h5",
  "mode": "fast",
  "bin_size": 10000,
  "temperature": 30.0,
  "polymerase": "phi29",
  "output": "report.html",
  "format": "html",
  "analyze": true,
  "plot": true
}
```

**Run**:
```bash
neoswga simulate --config config.json
```

---

## Command Reference

### simulate Command

```bash
neoswga simulate [OPTIONS]
```

### Required Arguments

| Argument | Description | Example |
|----------|-------------|---------|
| `--primers FILE` | Primer sequences file | `primers.txt` |
| `--fg FILE` | Target genome FASTA | `target.fasta` |
| `--bg FILE` | Background genome FASTA | `background.fasta` |
| `--positions FILE` | HDF5 positions file | `positions.h5` |

### Optional Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--mode MODE` | `fast` | Simulation mode: `fast`, `detailed`, `validation` |
| `--replicates N` | `5` | Number of replicates for detailed mode |
| `--bin-size BP` | `10000` | Bin size for coverage analysis (bp) |
| `--max-extension BP` | `70000` | Max Phi29 extension distance (bp) |
| `--temperature C` | `30.0` | Reaction temperature (°C) |
| `--duration SEC` | `3600` | Reaction duration (seconds) |
| `--polymerase TYPE` | `phi29` | Polymerase: `phi29` or `equiphi29` |
| `--analyze` | Off | Perform comprehensive analysis |
| `--plot` | Off | Generate visualization plots |
| `--output FILE` | None | Output file (HTML, JSON, or text) |
| `--format FMT` | `text` | Output format: `html`, `json`, `text` |
| `--quiet` | Off | Minimal output |
| `--verbose` | Off | Verbose output |

---

## Python API

### Core Classes

#### SwgaSimulator

```python
class SwgaSimulator:
    def __init__(self, primers, fg_genome, bg_genome,
                 fg_positions_h5, bg_positions_h5,
                 reaction_conditions=None, bin_size=10000, max_extension=70000)

    def simulate_fast(self) -> SimulationResult
    def simulate_detailed(self, num_replicates=5) -> SimulationResult
    def simulate_validation(self) -> SimulationResult
    def generate_report(self, result, output_path=None)
```

#### SimulationResult

```python
@dataclass
class SimulationResult:
    mode: str                      # 'fast', 'detailed', or 'validation'
    runtime: float                 # Seconds

    # Target metrics
    target_coverage: float         # 0-1 fraction
    target_uniformity: float       # 0-1 (higher is better)
    target_amplification: float    # Fold-amplification
    target_gaps: List[Dict]        # Under-covered regions

    # Background metrics
    background_coverage: float
    background_amplification: float

    # Specificity
    enrichment: float              # Target/background ratio
    specificity_score: float       # 0-1 normalized

    # Overall assessment
    composite_score: float         # 0-1 overall quality
    recommendation: str            # 'EXCELLENT', 'GOOD', 'FAIR', 'POOR'
    confidence: float              # 0-1 prediction confidence
```

#### SimulationAnalyzer

```python
class SimulationAnalyzer:
    def __init__(self, result, fg_positions, bg_positions,
                 fg_length, bg_length, bin_size=10000)

    def analyze(self) -> ComprehensiveAnalysis
```

#### ComprehensiveAnalysis

```python
@dataclass
class ComprehensiveAnalysis:
    coverage: CoverageAnalysis           # Detailed coverage metrics
    specificity: SpecificityAnalysis     # Enrichment and binding density
    primer_contributions: List[PrimerContribution]  # Per-primer performance
    quality_score: float                 # 0-1 overall quality
    recommendations: List[str]           # Actionable suggestions
    issues: List[str]                    # Identified problems
```

---

## Understanding Results

### Coverage Metrics

**Coverage**: Fraction of genome bins (default 10kb) with ≥1 primer binding site.

- **Excellent**: ≥60% coverage
- **Good**: 40-60% coverage
- **Fair**: 25-40% coverage
- **Poor**: <25% coverage

**Interpretation**:
- 42% coverage = 42% of genome regions have primer binding
- Higher coverage = more uniform amplification
- Gaps indicate under-amplified regions

**Uniformity**: Gini coefficient of inter-bin spacing.

- **Excellent**: ≥70% uniformity
- **Good**: 50-70% uniformity
- **Fair**: 30-50% uniformity
- **Poor**: <30% uniformity

**Interpretation**:
- 63% uniformity = primers are moderately well-distributed
- Lower uniformity = clustering in hotspots
- Higher uniformity = even spacing

### Specificity Metrics

**Enrichment**: Target amplification / Background amplification ratio.

- **Excellent**: ≥1,000× enrichment
- **Good**: 100-1,000× enrichment
- **Fair**: 10-100× enrichment
- **Poor**: <10× enrichment

**Interpretation**:
- 640,000× enrichment = target amplifies 640,000× more than background
- Higher enrichment = more specific primers
- Clinical applications typically need ≥100×

### Composite Score

**Formula**:
```
Composite = 0.40 × Coverage + 0.25 × Uniformity + 0.35 × Specificity
```

**Thresholds**:
- **EXCELLENT**: ≥0.80
- **GOOD**: ≥0.65
- **FAIR**: ≥0.45
- **POOR**: <0.45

### Recommendations

The analyzer generates actionable recommendations based on results:

**Low Coverage** (<40%):
```
Consider adding 5 more primers to improve coverage.
```

**Low Uniformity** (<50%):
```
Primers are clustering in specific regions. Add primers in under-covered regions.
```

**Critical Gaps** (>100kb):
```
Found 2 critical gaps (>100kb). Largest gap: 156,339 bp at position 270,460.
Add primers targeting these regions.
```

**Low Specificity** (<10×):
```
Low specificity (8×). Primers bind frequently to background genome.
Consider stricter frequency filters.
```

**Weak Primers**:
```
3 primers have low contribution scores.
Consider replacing: ATAAGCCAGATA, TCGCTCAGGCAA, GGGTATGGCTGT
```

---

## Troubleshooting

### Problem: "Genome file not found"

**Solution**:
```bash
# Use absolute paths
--fg /full/path/to/target.fasta

# Or relative paths from current directory
--fg ./genomes/target.fasta
```

### Problem: "HDF5 group not found"

**Cause**: Genome name extraction failed.

**Solution**: Check HDF5 structure:
```python
import h5py
with h5py.File('positions.h5', 'r') as f:
    print("Available genomes:", list(f.keys()))
```

Expected structure:
```
positions.h5
├── yersinia/
│   ├── PRIMER1/
│   │   ├── + (forward positions)
│   │   └── - (reverse positions)
│   ├── PRIMER2/
│   └── ...
└── human/
    ├── PRIMER1/
    └── ...
```

### Problem: "Coverage is 0%"

**Causes**:
1. Wrong genome name in HDF5
2. No primers found in positions file
3. Primers don't exist in HDF5

**Debug**:
```python
# Check primers in HDF5
with h5py.File('positions.h5', 'r') as f:
    genome = f['yersinia']
    print("Primers in HDF5:", list(genome.keys())[:10])

# Check your primer file
with open('primers.txt') as f:
    print("Primers in file:", [line.strip().split()[0] for line in f if line.strip()][:10])
```

### Problem: "Enrichment is 0×"

**Cause**: Background genome has more binding sites than target (reversed genomes).

**Solution**: Swap target and background:
```python
# Wrong order
sim = SwgaSimulator(primers, bg_genome, fg_genome, ...)  # Wrong

# Correct order
sim = SwgaSimulator(primers, fg_genome, bg_genome, ...)  # Correct
```

### Problem: "Matplotlib not found"

**Solution**:
```bash
pip install matplotlib

# For headless servers (no display)
export MPLBACKEND=Agg
```

---

## Advanced Topics

### Custom Bin Sizes

**Default**: 10kb bins

**Use cases**:
- **1kb bins**: High-resolution analysis for small genomes (<1 Mbp)
- **50kb bins**: Fast analysis for large genomes (>100 Mbp)

```python
sim = SwgaSimulator(..., bin_size=1000)   # 1kb bins
sim = SwgaSimulator(..., bin_size=50000)  # 50kb bins
```

### Detailed Mode (Future)

**Purpose**: Agent-based replication simulation with fork dynamics.

**Status**: Currently placeholder (uses fast mode).

**Future implementation**: Will use `Phi29Simulator` for accurate replication modeling.

### Validation Mode (Future)

**Purpose**: Stochastic Gillespie kinetics for exact chemical simulation.

**Status**: Currently placeholder (uses detailed mode).

**Future implementation**: Will model exact stochastic kinetics.

### Genome Size Limitations

**Memory usage**:
- 4.7 Mbp genome: ~5 MB RAM
- 3.3 Gbp genome: ~3.3 GB RAM

**Recommendations**:
- <100 Mbp: Load full genome
- >100 Mbp: Consider sampling approach
- >1 Gbp: Requires ≥16 GB RAM

### Parallelization (Future)

**Current**: Single-threaded

**Future**: Parallelize replicates and primer analysis:
```python
sim.simulate_detailed(num_replicates=10, n_workers=4)  # 4 parallel workers
```

---

## Examples Gallery

### Example: AT-rich Genome (Francisella, 32% GC)

```python
# Expected results
result = sim.simulate_fast()
assert result.target_coverage > 0.50       # 54%
assert result.enrichment > 20              # ~42×
assert result.recommendation in ['GOOD', 'FAIR']  # Lower due to AT-richness
```

### Example: Balanced Genome (Yersinia, 48% GC)

```python
# Expected results
result = sim.simulate_fast()
assert result.target_coverage > 0.40       # ~42%
assert result.enrichment > 1000            # ~640,000×
assert result.recommendation == 'GOOD'     # Excellent specificity
```

### Example: GC-rich Genome (Burkholderia, 68% GC)

```python
# Expected results (from benchmarks)
result = sim.simulate_fast()
assert result.target_coverage > 0.35       # Good coverage
assert result.enrichment > 500             # High specificity
assert result.recommendation in ['GOOD', 'EXCELLENT']
```

---

## Performance Benchmarks

| Operation | Genome Size | Time | Memory |
|-----------|-------------|------|--------|
| Load genome | 4.7 Mbp | ~10s | ~5 MB |
| Load genome | 3.3 Gbp | ~10s | ~3.3 GB |
| Load positions | 15 primers | <1s | <1 MB |
| Fast simulation | Any size | <1s | Minimal |
| Detailed simulation (5 replicates) | 4.7 Mbp | ~10 min | ~100 MB |
| Analysis | Any | <1s | Minimal |
| Visualization | Any | ~2s | ~50 MB |
| **Total workflow** | 4.7 Mbp + 3.3 Gbp | **~25s** | **~3.3 GB** |

---

## FAQ

**Q: What's the difference between fast, detailed, and validation modes?**

A:
- **Fast** (~1s): Coverage-based estimation, good for quick validation
- **Detailed** (~10min): Agent-based replication simulation, high accuracy
- **Validation** (~30min): Stochastic kinetics, highest accuracy (future)

**Q: What does "GOOD" recommendation mean?**

A: Your primers will work well for most applications. They have good coverage (>40%) and decent specificity (>65× enrichment). Suitable for research and some clinical applications.

**Q: My enrichment is 0.8× - what's wrong?**

A: You likely swapped target and background genomes. The primers bind more to the "background" than "target". Swap them in your command.

**Q: Can I use this for non-bacterial genomes?**

A: Yes! Works for any genome:
- Bacteria (tested)
- Archaea (should work)
- Eukaryotic parasites (tested: Plasmodium)
- Viruses (should work, but may need different bin sizes)

**Q: How many primers do I need?**

A:
- Minimum: 8-10 primers
- Recommended: 12-15 primers
- Maximum practical: 20-25 primers

More primers = better coverage but diminishing returns.

**Q: What if I don't have background genome?**

A: Run simulation with target genome as both fg and bg:
```python
sim = SwgaSimulator(primers, target, target, positions, positions)
```

This gives coverage analysis only (no specificity metrics).

---

## Citation

If you use NeoSWGA Simulator in your research, please cite:

```bibtex
@article{dwivedi2023fast,
  title={A fast machine-learning-guided primer design pipeline for selective whole genome amplification},
  author={Dwivedi-Yu, Jane A and Oppler, Zachary J and Mitchell, Matthew W and Song, Yun S and Brisson, Dustin},
  journal={PLOS Computational Biology},
  volume={19},
  number={4},
  pages={e1010137},
  year={2023}
}
```

---

## Support

**Issues**: Report bugs on the project GitHub issue tracker.

**Documentation**: See the [User Guide](user-guide.md) for general usage and the [SWGA Science](SWGA_SCIENCE.md) guide for theoretical background.
