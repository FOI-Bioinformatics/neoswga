# NeoSWGA Documentation

Documentation for NeoSWGA, a tool for designing primer sets for selective whole-genome amplification.

## Quick Links

- **[Quick Start](QUICK_START.md)** - Get started in minutes
- **[User Guide](user-guide.md)** - Comprehensive usage documentation
- **[SWGA Science](SWGA_SCIENCE.md)** - Scientific background and additive guide
- **[API Reference](API_REFERENCE.md)** - Complete API documentation
- **[Changelog](CHANGELOG.md)** - Version history and release notes

## Scientific Background

| Document | Description |
|----------|-------------|
| [SWGA Science](SWGA_SCIENCE.md) | Thermodynamics, polymerases, and additives |

## User Guides

| Guide | Description |
|-------|-------------|
| [Quick Start](QUICK_START.md) | Installation and first primer design |
| [User Guide](user-guide.md) | Complete usage documentation |
| [Report Generation](user_guide_reports.md) | Quality reports and grading |
| [Optimization Guide](optimization_guide.md) | Choosing optimization methods |
| [Multi-Genome Guide](multi-genome-guide.md) | Multiple target/background genomes |
| [Simulator Guide](simulator-guide.md) | SWGA replication simulator |
| [Adaptive QA Guide](adaptive-qa-guide.md) | AT-rich and GC-rich genomes |
| [Migration Guide](migration-guide.md) | Migrating from SOAPswga |
| [From Results to Lab](FROM_RESULTS_TO_LAB.md) | Export, ordering, and lab setup |

## Reference Documentation

| Document | Description |
|----------|-------------|
| [API Reference](API_REFERENCE.md) | CLI commands and Python API |
| [Module Reference](MODULE_REFERENCE.md) | Core module documentation |
| [Report Module](report_module.md) | Report generation system |
| [Report API Reference](api_reference_report.md) | Report module API |
| [Architecture Diagrams](ARCHITECTURE_DIAGRAMS.md) | Visual system architecture |

## Developer Documentation

| Document | Description |
|----------|-------------|
| [Developer Guide](DEVELOPER_GUIDE.md) | Contributing and extending NeoSWGA |
| [Architecture](development/architecture.md) | System architecture and design |
| [Algorithms](development/algorithms.md) | Filtering and optimization algorithms |
| [Optimizers](development/optimizers.md) | Optimization strategy details |
| [Implementation Guide](development/implementation-guide.md) | Development guidelines |
| [Deployment](development/deployment.md) | Deployment and configuration |

## Validation

Run validation tests:

```bash
neoswga validate --quick              # Installation check
neoswga validate-model                # Mechanistic model validation
pytest tests/ -v                      # Full test suite
```

## Directory Structure

```
docs/
  README.md                   # This index
  QUICK_START.md              # Quick start guide
  user-guide.md               # User manual
  user_guide_reports.md       # Report generation guide
  optimization_guide.md       # Optimization methods
  FROM_RESULTS_TO_LAB.md      # Export, ordering, and lab setup
  SWGA_SCIENCE.md             # Scientific background and additives
  API_REFERENCE.md            # API documentation
  api_reference_report.md     # Report API reference
  MODULE_REFERENCE.md         # Module documentation
  report_module.md            # Report module docs
  ARCHITECTURE_DIAGRAMS.md    # Visual architecture
  DEVELOPER_GUIDE.md          # Developer guide
  CHANGELOG.md                # Version history
  adaptive-qa-guide.md        # Adaptive QA for extreme GC genomes
  multi-genome-guide.md       # Multi-genome usage
  simulator-guide.md          # Simulator guide
  migration-guide.md          # SOAPswga migration
  development/                # Developer docs
    architecture.md
    algorithms.md
    optimizers.md
    implementation-guide.md
    deployment.md
```

## Citation

If you use NeoSWGA in your research, please cite:

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
