# NeoSWGA: Selective Whole Genome Amplification Primer Design

NeoSWGA is a command-line tool for designing primer sets for selective whole-genome amplification (SWGA). It identifies primers that amplify a target genome while minimizing amplification of background genomes. The tool combines machine learning, network-based optimization, and thermodynamic modeling.

**Primary use case**: Designing primers for Phi29/EquiPhi29 polymerase-based whole genome amplification, commonly used for pathogen detection from mixed samples.

## Key Features

- **Adaptive GC filtering**: Support for extreme GC genomes (32-68% GC)
- **Network-based optimization**: 10-100x faster than greedy search
- **Thermodynamic modeling**: SantaLucia nearest-neighbor calculations with salt corrections
- **Background filtering**: Bloom filter for large background genomes (human 3 Gbp)
- **Host-free mode**: Design primers without a background genome (`--no-background`)
- **Position cache**: 1000x faster position lookups
- **Multiple optimizers**: Network, greedy, genetic algorithm, MILP, clique, and more
- **Export formats**: FASTA, vendor CSV, BED, and BedGraph for genome browser visualization

## Installation

**Requirements**: Python >= 3.11, [Jellyfish](https://www.cbcb.umd.edu/software/jellyfish/) k-mer counter

```bash
# Basic installation
pip install -e .

# With all features (Bloom filter, MILP optimizer)
pip install -e ".[improved]"

# Validate installation
neoswga validate --quick
```

## Getting Started

### Interactive Setup (Recommended for new users)

```bash
# Setup wizard - creates params.json with recommended settings
neoswga init --genome target.fasta --background host.fasta

# Or use the interactive menu to discover all features
neoswga start
```

### Validate Configuration

```bash
# Check params.json for errors before running
neoswga validate-params -j params.json

# Get reaction condition recommendations
neoswga suggest --genome-gc 0.50 --primer-length 12
```

### Interpret Results

```bash
# After pipeline completes, get quality assessment
neoswga interpret -d results/
```

## Quick Start

```bash
# Single-command pipeline:
neoswga design -j params.json

# Or run each step individually:
neoswga count-kmers -j params.json    # Generate k-mer counts
neoswga filter -j params.json         # Filter candidate primers
neoswga score -j params.json          # Score amplification potential
neoswga optimize -j params.json       # Find optimal primer sets

# Host-free optimization (no background genome required):
neoswga optimize -j params.json --no-background
```

Example `params.json`:
```json
{
  "fg_genomes": ["target_genome.fasta"],
  "bg_genomes": ["background_genome.fasta"],
  "fg_prefixes": ["target_genome"],
  "bg_prefixes": ["background_genome"],
  "data_dir": "./data/",
  "fg_seq_lengths": [3500000],
  "bg_seq_lengths": [4600000]
}
```
Or use the wizard to generate params.json automatically: `neoswga init --genome target.fasta`

## Documentation

- **[Quick Start](docs/QUICK_START.md)**: Installation and first primer design
- **[User Guide](docs/user-guide.md)**: Comprehensive usage documentation
- **[Optimization Guide](docs/optimization_guide.md)**: Choosing the right optimization method
- **[From Results to Lab](docs/FROM_RESULTS_TO_LAB.md)**: Export primers and lab workflow
- **[Report Generation](docs/user_guide_reports.md)**: Quality reports and grading
- **[Multi-Genome Guide](docs/multi-genome-guide.md)**: Pan-genome primer design
- **[SWGA Science](docs/SWGA_SCIENCE.md)**: Thermodynamics, polymerases, and reaction additives
- **[Changelog](docs/CHANGELOG.md)**: Version history

See [docs/README.md](docs/README.md) for the full documentation index.

## Based on SOAPswga

NeoSWGA extends SOAPswga, originally developed by Dwivedi-Yu et al. (2023):

```bibtex
@article{dwivedi2023fast,
  title={A fast machine-learning-guided primer design pipeline for selective whole genome amplification},
  author={Dwivedi-Yu, Jane A and Oppler, Zachary J and Mitchell, Matthew W and Song, Yun S and Brisson, Dustin},
  journal={PLOS Computational Biology},
  volume={19},
  number={4},
  pages={e1010137},
  year={2023},
  publisher={Public Library of Science San Francisco, CA USA}
}
```

## License

MIT License

## Contributing

Contributions welcome. See [CLAUDE.md](CLAUDE.md) for architecture and development guidelines.
