# NeoSWGA Documentation

Documentation for NeoSWGA, a tool for designing primer sets for selective whole-genome amplification.

## Quick Links

- **[Quick Start](guides/QUICK_START.md)** - Get started in minutes
- **[User Guide](guides/user-guide.md)** - Comprehensive usage documentation
- **[SWGA Science](SWGA_SCIENCE.md)** - Scientific background and additive guide
- **[API Reference](reference/API_REFERENCE.md)** - Complete API documentation
- **[Troubleshooting](guides/TROUBLESHOOTING.md)** - When a run fails or disappoints
- **[Changelog](CHANGELOG.md)** - Version history and release notes

## Layout

```
docs/
  README.md                 # this index
  CHANGELOG.md               # version history
  SWGA_SCIENCE.md            # scientific background (see Scientific Background below)
  SCIENCE_CITATIONS.md       # published results each default is grounded in
  params-reference.md        # generated from params.schema.json -- see note below
  production-scenarios.md    # worked configurations for common situations
  guides/                    # task-oriented how-to documentation
  migrations/                # version-upgrade guides
  reference/                 # API, module and architecture reference
  development/               # contributing to NeoSWGA itself
  reports/                   # point-in-time audits, SWOT, tech-debt -- still cited, not superseded
  validation/                # measurement records backing specific claims/defaults
  archive/                   # historical/superseded documents, kept for the record
```

`params-reference.md` and `SCIENCE_CITATIONS.md` stay at the top level rather
than moving into `reference/`: their paths are hardcoded in
`scripts/render_schema.py` and in several test assertions
(`tests/test_schema_renderer.py`, `tests/test_set_size_flags_are_documented.py`,
`tests/test_realistic_amplicon_length.py`). `production-scenarios.md` and
`SWGA_SCIENCE.md` stay alongside them so the schema renderer's embedded links
to both don't need touching.

## Scientific Background

| Document | Description |
|----------|-------------|
| [SWGA Science](SWGA_SCIENCE.md) | Thermodynamics, polymerases, and additives |
| [Science Citations](SCIENCE_CITATIONS.md) | The published results each default is grounded in |

## Guides

| Guide | Description |
|-------|-------------|
| [Quick Start](guides/QUICK_START.md) | Installation and first primer design |
| [User Guide](guides/user-guide.md) | Complete usage documentation |
| [Report Generation](guides/user_guide_reports.md) | Quality reports and grading |
| [Optimization Guide](guides/optimization_guide.md) | Choosing optimization methods |
| [Multi-Genome Guide](guides/multi-genome-guide.md) | Multiple target/background genomes |
| [Simulator Guide](guides/simulator-guide.md) | SWGA replication simulator |
| [Adaptive QA Guide](guides/adaptive-qa-guide.md) | AT-rich and GC-rich genomes |
| [From Results to Lab](guides/FROM_RESULTS_TO_LAB.md) | Export, ordering, and lab setup |
| [Troubleshooting](guides/TROUBLESHOOTING.md) | Diagnosing failed or disappointing runs |
| [Production Scenarios](production-scenarios.md) | Worked configurations for common situations |
| [Active Learning Guide](guides/active-learning-guide.md) | Feeding wet-lab outcomes back into design |
| [Training amp_pred on real data](guides/training_amp_pred_on_real_data.md) | Retraining the retired RF scorer against lab outcomes instead of synthetic labels |

## Migrations

| Document | Description |
|----------|-------------|
| [Migration Guide](migrations/migration-guide.md) | Migrating from SOAPswga, and pre-3.5 -> 3.5 (genome-adaptive QA) |
| [3.6 -> 3.7](migrations/migration-3.6-to-3.7.md) | Polymerase-aware defaults |

## Reference Documentation

| Document | Description |
|----------|-------------|
| [API Reference](reference/API_REFERENCE.md) | CLI commands and Python API |
| [Module Reference](reference/MODULE_REFERENCE.md) | Core module documentation -- **carries an accuracy warning, see the banner at the top of the file** |
| [Report Module](reference/report_module.md) | Report generation architecture and CLI/Python usage |
| [Report API Reference](reference/api_reference_report.md) | Report module class/function signatures |
| [Architecture Diagrams](reference/ARCHITECTURE_DIAGRAMS.md) | Visual system architecture |
| [Parameter Reference](params-reference.md) | Every params.json key, generated from the schema |

## Developer Documentation

| Document | Description |
|----------|-------------|
| [Developer Guide](development/DEVELOPER_GUIDE.md) | Contributing and extending NeoSWGA |
| [Architecture](development/architecture.md) | System architecture and design -- **carries an accuracy warning, see the banner at the top of the file** |
| [Algorithms](development/algorithms.md) | Filtering and optimization algorithms |
| [Optimizers](development/optimizers.md) | **Superseded 2026-08-31** -- describes a retired 17-optimizer architecture; read [Optimization Guide](guides/optimization_guide.md) instead |
| [Deployment](development/deployment.md) | Deployment and configuration |

`MODULE_REFERENCE.md` and `development/architecture.md` both documented
modules that no longer exist in `neoswga/core/` (a 2026-09-11 survey found 12
in the former, 3 in the latter); both sets of ghost sections have now been
removed. What's still outstanding on each is coverage, not accuracy: 27 of
92 real modules have no section in `MODULE_REFERENCE.md`, and 86 of 92 have
none in `architecture.md`. See CLAUDE.md's "Core Modules" section for what
actually ships.

## Reports (audits, SWOT, technical debt)

Point-in-time analysis documents. They are not superseded by being old: the
pipeline audit is the one to read first, and CLAUDE.md cites finding
identifiers from it and from the optimize-step audit (B4, B5, A6, E6, F0,
F1b) directly.

| Document | Description |
|----------|-------------|
| [Four-step pipeline audit](reports/AUDIT_pipeline_four_steps_2026-09-06.md) | 2026-09. count-kmers through optimize, measured on three whole-genome designs. Defines findings B4, B5, A6 and E6. |
| [Optimize-step audit](reports/AUDIT_optimize_step_2026-09-03.md) | 2026-09. Step 4 in depth. Defines findings F0 and F1b. |
| [EquiPhi29 panel audit](reports/AUDIT_equiphi29_panel_2026-09-03.md) | 2026-09. Dimer content and Tm window of the delivered panels |
| [Alternatives and scaling audit](reports/AUDIT_2026-08_alternatives_and_scaling.md) | 2026-08. Optimizer scaling, option reachability, and oligo-set quality against an exact optimum |
| [SWOT: pipeline as workflow](reports/SWOT_2026-08_pipeline_as_workflow.md) | 2026-08. Strengths/weaknesses/opportunities/threats of the four-step pipeline |
| [Audit completion plan](reports/PLAN_2026-08_audit_completion.md) | 2026-08. Plan for closing out the alternatives-and-scaling audit findings |
| [Technical debt assessment](reports/TECH_DEBT_2026-09.md) | 2026-09. Structural debt (config architecture, complexity concentration, dead code), complementary to the alternatives-and-scaling audit |

Ten point-in-time reports from 2026-03 and 2026-04 were removed on 2026-09-11:
the four per-step audits, `AUDIT_REPORT_v3.0.md`, `OPTIMIZER_UX_AUDIT.md`, two
validation reports and two production-readiness assessments. They predated
roughly 20,000 changed lines under `neoswga/core/`, nothing linked to them, and
six carried their own "superseded in part" banner. They remain in git history.

## Validation

Measurement records backing specific claims and defaults elsewhere in the docs:

| Document | Description |
|----------|-------------|
| [Additive specificity](validation/additive_specificity.md) | The additive-to-specificity lever, measured end to end |
| [Reach calibration](validation/reach_calibration.md) | Per-primer extension reach fitted to a published outcome |
| [Published primer sets](validation/published_primer_sets.md) | Scoring checked against seven wet-lab datasets |
| [Tool comparison](validation/tool_comparison.md) | NeoSWGA against swga 1.0, swga 2.0/soapswga and COATswga |
| [Optimizer cost](validation/optimizer_cost_2026-09.md) | Runtime and set quality per method across panel sizes |
| [Dimer threshold trade-off](validation/dimer_threshold_tradeoff_2026-09.md) | What panel size each pool supports at each `max_dimer_bp` |

Run validation tests:

```bash
neoswga validate --quick              # Installation check
neoswga validate-model                # Mechanistic model validation
pytest tests/ -v                      # Full test suite
```

## Archive

Historical documents, kept for the record rather than as current guidance.
Each carries its own banner explaining what superseded it.

| Document | Description |
|----------|-------------|
| [Implementation Guide](archive/implementation-guide.md) | Migration notes from a pre-`core/` pipeline; references modules that no longer exist |
| [Repo cleanup design (2026-02)](archive/plans/2026-02-03-repo-cleanup-design.md) | Completed |
| [Export and lab integration plan (2026-02)](archive/plans/2026-02-04-export-and-lab-integration.md) | Completed |

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
