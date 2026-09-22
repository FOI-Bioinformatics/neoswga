# NeoSWGA: Selective Whole Genome Amplification Primer Design

NeoSWGA is a command-line tool for designing primer sets for selective whole-genome amplification (SWGA). It identifies primers that amplify a target genome while minimizing amplification of background genomes, combining graph and network optimization with thermodynamic modelling. A pre-trained random forest ships with the tool but was retired from the default path on 2026-09-05; `--amp-model` restores it.

**Primary use case**: Designing primers for Phi29/EquiPhi29 polymerase-based whole genome amplification, commonly used for pathogen detection from mixed samples.

> **Before trusting a number from a design, read [what NeoSWGA does not establish](docs/LIMITATIONS.md).**
> Per capability, [what the evidence actually is](docs/EVIDENCE.md) says whether it is implemented, connected, oracle-tested or retrospectively evaluated. Nothing here is prospectively validated.
> Coverage is a modelled geometric proxy, not predicted sequencing
> breadth, and no design from this tool has been tested in a
> laboratory as part of its development.

## Key Features

- **Adaptive GC filtering**: Support for extreme GC genomes (32-68% GC)
- **Several optimization methods**: graph set-cover, network, background-aware
  and dimer-free clique selection, with an exact solver for bounding the result
- **Thermodynamic modeling**: SantaLucia nearest-neighbor calculations with salt corrections
- **Background filtering**: Bloom filter for large background genomes (human 3 Gbp)
- **Host-free mode**: Design primers without a background genome (`--no-background`)
- **Position cache**: 1000x faster position lookups
- **Multiple optimizers**: `hybrid`, `dominating-set`, `network`, `background-aware`, plus `ensemble` (run several and keep the best by normalized score; `--ensemble-combine union` re-optimizes over the pooled primers to beat any single method)
- **Iterative design from real data**: add oligos to a validated set using in-silico and real sequencing-depth (BAM) coverage gaps (`analyze-coverage`, `expand-primers --bam`; needs the `[bam]` extra)
- **Comprehensive reports**: the technical report surfaces every in-silico result (ensemble comparison, per-target coverage, strand balance, coverage gaps, reaction conditions), each value badged MEASURED or ESTIMATED
- **Export formats**: FASTA, vendor CSV, BED, and BedGraph for genome browser visualization

## Installation

**Requirements**: Python >= 3.11, [Jellyfish](https://www.cbcb.umd.edu/software/jellyfish/) k-mer counter

```bash
pip install neoswga

# With visualization support
pip install neoswga[viz]

# With all optional features
pip install neoswga[all]

# Validate installation
neoswga validate --quick
```

### Development

```bash
# For development
pip install -e ".[dev]"
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
neoswga prepare-candidates -j params.json          # Score amplification potential
neoswga optimize -j params.json       # Find optimal primer sets

# Host-free optimization (no background genome required):
neoswga optimize -j params.json --no-background

# Not sure how many primers to ask for? `optimize` always prints a marginal
# coverage table (pp/primer) showing where extra primers stop paying, at any
# set size. For an estimate before the run, or a coverage-against-specificity
# frontier, both of which stop at 20 primers:
neoswga optimize -j params.json --auto-size       # coverage estimate, up to 20
neoswga optimize -j params.json --show-frontier   # trade-off frontier, 4 to 20
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

- **[Quick Start](docs/guides/QUICK_START.md)**: Installation and first primer design
- **[User Guide](docs/guides/user-guide.md)**: Comprehensive usage documentation
- **[Optimization Guide](docs/guides/optimization_guide.md)**: Choosing the right optimization method
- **[From Results to Lab](docs/guides/FROM_RESULTS_TO_LAB.md)**: Export primers and lab workflow
- **[Report Generation](docs/guides/user_guide_reports.md)**: Quality reports and grading
- **[Multi-Genome Guide](docs/guides/multi-genome-guide.md)**: Pan-genome primer design
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

## Deprecation Policy

NeoSWGA follows semantic versioning. When features are deprecated:

- **Deprecated features** emit a `DeprecationWarning` for at least one minor release before removal.
- **Removed features** are documented in the CHANGELOG with migration guidance.
- **params.json changes** are backwards compatible within the same major version. New parameters use sensible defaults.
- **CLI flag changes** follow the same deprecation cycle: warning first, removal in next major.

## License

AGPL-3.0-or-later. See [LICENSE](LICENSE) for details.

## Contributing

Contributions welcome. See [CLAUDE.md](CLAUDE.md) for architecture and development guidelines.
