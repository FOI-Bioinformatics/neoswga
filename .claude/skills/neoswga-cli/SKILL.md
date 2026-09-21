---
name: neoswga-cli
description: Reference for NeoSWGA commands outside the four-step pipeline - init, start, suggest, validate, interpret, report, multi-genome, simulate, analyze-set, analyze-genome, analyze-dimers, analyze-coverage, expand-primers - plus the mechanistic-model flags on optimize, RF model retraining, and the plasmid example. Use for any neoswga subcommand other than count-kmers, filter, prepare-candidates and optimize.
---

# NeoSWGA CLI reference (beyond the four-step pipeline)

The standard pipeline (`count-kmers`, `filter`, `prepare-candidates`, `optimize`), the
optimization methods, the params.json reference and the known issues live in
the repository's `CLAUDE.md`. This file covers everything else.

## Setup commands

```bash
# Interactive workflow selector - discover all features
neoswga start

# Setup wizard - create params.json with guided configuration
neoswga init --genome target.fna [--background host.fna] [-o params.json]

# Validate params.json before running pipeline
neoswga validate params -j params.json

# Suggest optimal reaction conditions
neoswga suggest --genome-gc 0.65 --primer-length 15
neoswga suggest --genome target.fna  # Auto-calculates GC

# Interpret results after pipeline completes
neoswga interpret -d results/

# Validate mechanistic model against expected behavior
neoswga validate model               # Run all validation tests
neoswga validate model --output-json  # Output results as JSON
```

The hyphenated forms (`validate-params`, `validate-model`) still run but warn:
the subcommands are `neoswga validate {install,params,model}`.

## Quality reports

```bash
neoswga report -d results/                        # Executive summary (default)
neoswga report -d results/ --level full           # Full technical report
neoswga report -d results/ --interactive          # With interactive Plotly charts
neoswga report -d results/ --level full --interactive  # Full report with charts
neoswga report -d results/ --check                # Validate only, don't generate
```

The full technical report surfaces every in-silico result read from the
authoritative `step4_improved_df_summary.json` (preferred over CSV estimates):
the ensemble per-method comparison, per-target coverage, strand balance,
coverage gaps (in-silico +/- BAM depth), and reaction conditions. Every value
is badged MEASURED or ESTIMATED. The filtering funnel uses the real
`filter_stats.json` the filter step writes (no fabricated counts).
`--interactive` adds Plotly charts on top of the static sections; the chart
helpers in `neoswga/core/report/visualizations.py` return an empty string when
Plotly is not installed, so the report degrades rather than failing.

## Optimization with the mechanistic model

```bash
# Auto-size primer set based on application profile
neoswga optimize -j params.json --auto-size --application clinical

# Use mechanistic model for primer weighting (a non-default --mechanistic-weight
# implies --use-mechanistic-model)
neoswga optimize -j params.json --use-mechanistic-model --mechanistic-weight 0.3
```

Application profiles are `discovery` (high coverage), `clinical` (high
specificity), `enrichment` (balanced) and `metagenomics` (capture diversity);
their coverage and specificity targets are tabulated in `CLAUDE.md`.

The model in `neoswga/core/mechanistic_model.py` combines four pathways:
Tm modification (DMSO, betaine, formamide on primer-template stability),
secondary-structure accessibility (template melting, GC-dependent structure),
enzyme activity (polymerase processivity, speed, stability), and binding
kinetics (kon/koff).

## Advanced commands

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
```

## Adding oligos to an existing set (in-silico + real BAM coverage)

For iterative design: keep validated primers, exclude failed ones, and add new
primers that fill coverage gaps. Gaps can come from in-silico binding positions
and/or from real sequencing depth (a BAM mapped to the target genome). BAM
support needs the `[bam]` extra (`pip install 'neoswga[bam]'`, brings pysam).

```bash
# Inspect gaps only (read-only): writes coverage_gaps.bed + .json
neoswga analyze-coverage -j params.json --primers SEQ1 SEQ2 \
    --bam reads.bam --min-depth 5 --min-gap-size 10000 -o cov/

# Add primers, focusing the candidate pool on the merged gaps
neoswga expand-primers -j params.json --fixed-primers SEQ1 SEQ2 \
    --failed-primers SEQ3 --num-new 6 --bam reads.bam --min-depth 5 \
    --optimization-method hybrid -o expanded/
```

- BAM contigs are matched to foreground prefixes by exact/basename/`chr`-prefix
  then unique-length fallback; override with `--contig-alias FG=BAMCONTIG`.
- A base counts as a gap when its mapped depth `< --min-depth`; runs shorter
  than `--min-gap-size` bp are ignored. On circular targets (`fg_circular`) a
  gap spanning the origin is merged into one.
- Candidate selection is a HARD pre-filter (only primers binding inside a gap),
  with fallback to the full pool if too few remain to reach `--num-new`.
- Outputs: `merged_gaps.bed`, `expansion_result.json`, `expanded_primers.csv`.

## Development tasks

**Generate k-mer files for non-standard lengths**:
```bash
neoswga count-kmers -j params.json --min-k 15 --max-k 18
```

**Retrain random forest model** (for sklearn updates):
```bash
python scripts/retrain_rf_model.py --output neoswga/core/models/random_forest_filter.skops
```

**Convert a legacy pickle model to skops** (one-time migration helper):
```bash
python scripts/convert_model_to_skops.py
```

**Run pipeline on test data**:
```bash
cd tests/integration/equiphi29_baseline
neoswga count-kmers -j params.json
neoswga filter -j params.json
neoswga prepare-candidates -j params.json
neoswga optimize -j params.json
```

## Running the example

`examples/plasmid_example/` provides a quick test with two small plasmids
(pcDNA vs pLTR):

```bash
cd examples/plasmid_example
neoswga count-kmers -j params.json
neoswga filter -j params.json
neoswga prepare-candidates -j params.json
neoswga optimize -j params.json
```
