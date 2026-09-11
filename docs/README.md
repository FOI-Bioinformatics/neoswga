# NeoSWGA Documentation

Documentation for NeoSWGA, a tool for designing primer sets for selective whole-genome amplification.

## Quick Links

- **[Quick Start](QUICK_START.md)** - Get started in minutes
- **[User Guide](user-guide.md)** - Comprehensive usage documentation
- **[SWGA Science](SWGA_SCIENCE.md)** - Scientific background and additive guide
- **[API Reference](API_REFERENCE.md)** - Complete API documentation
- **[Troubleshooting](TROUBLESHOOTING.md)** - When a run fails or disappoints
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
| [Troubleshooting](TROUBLESHOOTING.md) | Diagnosing failed or disappointing runs |
| [Production Scenarios](production-scenarios.md) | Worked configurations for common situations |
| [Active Learning Guide](active-learning-guide.md) | Feeding wet-lab outcomes back into design |

## Reference Documentation

| Document | Description |
|----------|-------------|
| [API Reference](API_REFERENCE.md) | CLI commands and Python API |
| [Module Reference](MODULE_REFERENCE.md) | Core module documentation |
| [Report Module](report_module.md) | Report generation system |
| [Report API Reference](api_reference_report.md) | Report module API |
| [Architecture Diagrams](ARCHITECTURE_DIAGRAMS.md) | Visual system architecture |
| [Parameter Reference](params-reference.md) | Every params.json key, generated from the schema |

## Developer Documentation

| Document | Description |
|----------|-------------|
| [Developer Guide](DEVELOPER_GUIDE.md) | Contributing and extending NeoSWGA |
| [Architecture](development/architecture.md) | System architecture and design |
| [Algorithms](development/algorithms.md) | Filtering and optimization algorithms |
| [Optimizers](development/optimizers.md) | Optimization strategy details |
| [Implementation Guide](development/implementation-guide.md) | Development guidelines |
| [Deployment](development/deployment.md) | Deployment and configuration |

## Audits and validation reports

| Document | Description |
|----------|-------------|
Listed newest first. The pipeline audit is the one to read: CLAUDE.md cites its
finding identifiers (B4, B5, A6, E6), as it cites F0 and F1b from the optimize-step
audit below.

| Document | Description |
|----------|-------------|
| [Four-step pipeline audit](AUDIT_pipeline_four_steps_2026-09-06.md) | 2026-09. count-kmers through optimize, measured on three whole-genome designs. Defines findings B4, B5, A6 and E6. |
| [Optimize-step audit](AUDIT_optimize_step_2026-09-03.md) | 2026-09. Step 4 in depth. Defines findings F0 and F1b. |
| [EquiPhi29 panel audit](AUDIT_equiphi29_panel_2026-09-03.md) | 2026-09. Dimer content and Tm window of the delivered panels |
| [Alternatives and scaling audit](AUDIT_2026-08_alternatives_and_scaling.md) | 2026-08. Optimizer scaling, option reachability, and oligo-set quality against an exact optimum |

Measurement records:

| Document | Description |
|----------|-------------|
| [Additive specificity](validation/additive_specificity.md) | The additive-to-specificity lever, measured end to end |
| [Reach calibration](validation/reach_calibration.md) | Per-primer extension reach fitted to a published outcome |
| [Published primer sets](validation/published_primer_sets.md) | Scoring checked against seven wet-lab datasets |
| [Tool comparison](validation/tool_comparison.md) | NeoSWGA against swga 1.0, swga 2.0/soapswga and COATswga |
| [Optimizer cost](validation/optimizer_cost_2026-09.md) | Runtime and set quality per method across panel sizes |
| [Dimer threshold trade-off](validation/dimer_threshold_tradeoff_2026-09.md) | What panel size each pool supports at each `max_dimer_bp` |
| [Science citations](SCIENCE_CITATIONS.md) | The published results each default is grounded in |

Ten point-in-time reports from 2026-03 and 2026-04 were removed on 2026-09-11:
the four per-step audits, `AUDIT_REPORT_v3.0.md`, `OPTIMIZER_UX_AUDIT.md`, two
validation reports and two production-readiness assessments. They predated
roughly 20,000 changed lines under `neoswga/core/`, nothing linked to them, and
six carried their own "superseded in part" banner. They remain in git history.

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

AGPL-3.0-or-later. See the LICENSE file.
