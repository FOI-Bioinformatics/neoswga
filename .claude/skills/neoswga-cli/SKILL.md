---
name: neoswga-cli
description: Reference for NeoSWGA commands outside the four-step pipeline - init, start, suggest, validate, interpret, report, multi-genome, simulate, analyze-set, analyze-genome, analyze-dimers, analyze-coverage, calibrate-reach, evaluate-set, expand-primers, swap-primer, contract-set, rescore-set, report-pool, export, doctor - plus the mechanistic-model flags on optimize, RF model retraining, and the plasmid example. Use for any neoswga subcommand other than count-kmers, filter, prepare-candidates and optimize.
---

# NeoSWGA CLI reference (beyond the four-step pipeline)

The standard pipeline (`count-kmers`, `filter`, `prepare-candidates`, `optimize`), the
optimization methods and the key parameters are summarised in the repository's
`CLAUDE.md`, with the full record in `docs/development/implementation-notes.md`,
`docs/development/design-contracts.md` and `docs/development/KNOWN_ISSUES.md`.
This file covers everything else.

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
their coverage and specificity targets are tabulated in
`docs/development/implementation-notes.md` (Application profiles). On the
default `hybrid` method the profile does not change what is selected; it
weights the ensemble winner and steers `--auto-size`.

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

- BAM contigs are bound to foreground references by NAME only: explicit alias,
  exact match, basename, then `chr`-prefix normalisation. The earlier
  unique-length fallback was removed on 2026-09-21, since equal length is not
  identity. An unmatched reference is skipped with a warning; map it with
  `--contig-alias FG=BAMCONTIG`, where the key is a FASTA record name, or the
  prefix/basename when the reference holds one record.
- A CRAM needs `--reference` (the FASTA the reads were aligned to).
- A base counts as a gap when its mapped depth `< --min-depth`; runs shorter
  than `--min-gap-size` bp are ignored. On circular targets (`fg_circular`) a
  gap spanning the origin is merged into one.
- Candidate selection is a HARD pre-filter (only primers binding inside a gap),
  with fallback to the full pool if too few remain to reach `--num-new`.
- Outputs: `merged_gaps.bed`, `expansion_result.json`, `expanded_primers.csv`.

## Working on an existing set

```bash
# Coverage, gaps and dimers for any oligo set; --genome lets -j be omitted
neoswga evaluate-set --primers SEQ1 SEQ2 --genome target.fna -o eval/

# Replace under-performing primers from a candidate pool (default: step2_df.csv)
neoswga swap-primer -j params.json --primers SEQ1 SEQ2 --max-swaps 3 -o swaps.json

# Remove redundant primers while coverage stays above a floor
neoswga contract-set -j params.json --primers SEQ1 SEQ2 --min-coverage 0.70

# Rescore under the reaction in params.json, without re-optimizing
neoswga rescore-set -j params.json --primers SEQ1 SEQ2

# Fit the per-primer coverage reach to sequencing depth (read-only)
neoswga calibrate-reach -j params.json --primers SEQ1 SEQ2 --bam reads.bam -o reach.json
```

- `rescore-set` accepts additive and salt flags (`--dmso-percent`,
  `--betaine-m`, `--mg-conc`, ...) that are not merged on this path today;
  set the reaction in params.json instead.
- `calibrate-reach` cannot describe a multi-record reference and says so. It
  has not been run against measured sequencing depth in this repository.
- The GPU flags on `evaluate-set` are accepted and change nothing.

## Pool planning, export and diagnostics

```bash
# Smallest pools meeting coverage and specificity targets; one design per grid condition
neoswga plan-pool -j params.json [--design-grid grid.json]
neoswga report-pool --input pool_plan.json -o pool_report/   # new or empty directory

# Ordering files for set 0; --set N selects an alternative set
neoswga export -d results/ -o order/ --project NAME

# Installation, params.json and optimizer capability check
neoswga doctor -j params.json [--json]
```

`export` writes nothing when the last design run failed (`design_failure.json`
is present), when the delivered set carries a recorded defect, or when
`--set` names a set other than 0, which nothing assessed.
`--allow-unqualified` overrides the last two.

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
