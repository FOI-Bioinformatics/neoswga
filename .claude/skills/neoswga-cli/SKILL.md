---
name: neoswga-cli
description: Reference for NeoSWGA commands outside the four-step pipeline - init, start, suggest, validate, interpret, report, multi-genome, simulate, analyze-set, analyze-genome, analyze-dimers, analyze-coverage, calibrate-reach, evaluate-set, improve-set, expand-primers, swap-primer, contract-set, rescore-set, report-pool, export, doctor - plus the mechanistic-model flags on optimize, RF model retraining, and the plasmid example. Use for any neoswga subcommand other than count-kmers, filter, prepare-candidates and optimize.
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

# Score a delivered set against hosts that were never in the design
neoswga evaluate-set --from-results results/ --set 0 --genome target.fna \
    --background host1.fna host2.fna -o eval/

# Report which binding sites survive in which strain, from a variant table
neoswga evaluate-set --primers SEQ1 SEQ2 --genome reference.fna \
    --variants strains.vcf -o eval/

# Diagnose a set per oligo and propose edits to it; reports, applies nothing
neoswga improve-set -j params.json --primers SEQ1 SEQ2 \
    --background host1.fna --max-edits 5 -o improvement/

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
- `evaluate-set --background FASTA...` takes hosts by path, with no params.json
  entry; `--from-results DIR` reads one delivered set (`--set N`, default 0)
  instead of `--primers`. `evaluation.json` then carries `per_target` and
  `per_host` blocks, a `target_host_pairs` table and the two reductions
  `worst_target_coverage` and `worst_host_selectivity_density`. Every field it
  carried before keeps its name and its arithmetic.
- Host sites come from k-mer counts, which is a table lookup where one exists
  and a single pass over the reference otherwise. `--scan-background` locates
  them instead, which also measures host coverage and holds the host in memory:
  332 MB on a 144 Mb reference, 2.3 GB on hg38
  (`docs/validation/query_scan_2026-09-25.md`). Without it, a counted host
  reports its coverage as unavailable WITH that reason rather than as zero.
- A reduction over a set with one unmeasured member is `None`, not the worst of
  the rest, and the reason says which member was missing. Compare
  `selectivity_density` across hosts of different size, never
  `selectivity_ratio` (Known Issue 6).
- `evaluate-set --variants FILE` states target diversity as variants against
  one reference instead of one FASTA per strain. It takes a VCF or BCF (through
  `core/variant_table.py`, the only place a variant file is opened) or a TSV
  with a header of `chrom pos ref alt` and one optional 0/1 column per strain.
  Both formats are 1-based. `evaluation.json` then carries a `per_strain` block
  and a `variant_route` block, and each per-primer record gains
  `intact_fraction_every_strain` and `intact_fraction_by_strain`. There is no
  params.json key for it.
- The model is binary and needs none: a site `[pos, pos + k)` is intact in a
  strain when no variant that strain carries falls inside it, and affected
  otherwise. An affected site is counted as lost. It is never weighted,
  scored with `occupancy.mismatch_tm`, or called tolerated: that model is a
  uniform 4.0 C per mismatch with evidence status `assumed`.
- Affected sites are split by the distance of the nearest variant from the
  primer's 3' end, on the strand the primer binds, at a stated window
  (`three_prime_window_nt`, 5 bases). The split is reported only; nothing is
  compared against it and no penalty follows from it.
- Three limits are written into `variant_route.limits` and printed: a site
  GAINED in a strain through a variant is invisible, since sites are found on
  the reference; an indel affects every site it overlaps and the coordinate
  shift it causes downstream is not modelled, so gap lengths in a carrying
  strain are approximate; and the table says nothing about sequence absent
  from the reference.
- The reader refuses rather than repairs: a contig the FASTA does not have, a
  REF allele that disagrees with the FASTA (which is what catches a table made
  against another assembly, and a 0-based one), an unsorted or interleaved
  table, a TSV genotype that is neither 0 nor 1, a TSV row with more or fewer
  columns than its header, a genotype naming an allele its row does not have
  (`ALT=.` with `1/1`), a breakend ALT, a symbolic ALT (`<DEL>`, `<DUP>`, ...)
  without `INFO/END` (declared in the header), a span past the end of its
  record, and a `.vcf.gz` compressed with gzip rather than bgzip.
- A symbolic ALT with `INFO/END` affects `[POS - 1, END)`, padding base
  included; its REF alone is only that padding base.
- VCF genotypes: any allele naming an ALT makes a carrier (`./1` is one);
  otherwise any uncalled allele makes the genotype unknown (`0/.`, `./0`,
  `.|0`, `./.`, no GT); only a fully stated reference genotype is a
  non-carrier. A strain with an unknown genotype at any row is `unavailable`
  with that row as the reason, and the `worst_strain_*` reductions are then
  `None`, not the worst of the rest.
- Circular references (the `evaluate-set` default unless `--linear`): the
  scanner wraps the whole concatenated sequence once, so a site may start in
  the last k-1 bases and end in the first ones, joining the last record to the
  first on a multi-record reference. The mask and the 3'-end split follow that
  definition. A site whose geometry the run cannot place -- a wrap position
  read as linear, an offset past the end -- is counted as `not_assessed_sites`,
  never as intact, and the strain's fraction, coverage and gaps are then
  unavailable. `intact + affected + not_assessed` always equals
  `reference_sites`, which counts a palindromic oligo's site once.
- One reference per run. With one foreground genome that is the reference.
  With more than one configured -- readable or not -- the run is refused
  unless `--variants-reference FASTA` names one of them; the refusal lists
  them. With `-j`, do not add `--genome` to pick one: it replaces
  `fg_genomes`, leaves two prefixes and one FASTA, and is refused as such.
  The chosen reference is logged, printed above the per-strain rows, and
  recorded as `variant_route.reference_fasta`. When positions come from a
  position index, its stored record starts are checked against the FASTA
  (`reference_layout.verify_layout`) and a stale index is refused; the
  configured length must equal the FASTA's. A refusal leaves no
  `design_failure.json`, for the same reason `improve-set` leaves none.
- `improve-set` takes the same set and reference options as `evaluate-set`
  (`--primers` / `--primers-file` or `--from-results DIR --set N`, `--genome`,
  `--background`, `--scan-background`, `--linear`) and writes
  `improvement_report.json`: the per-reference figures, a per-oligo attribution
  (sites on each target and host, marginal coverage per target, dimer partners
  in the set, Tm against the window) and proposed single edits in four
  sections. It writes no `step4_improved_df.csv`; applying a proposal is done
  with `swap-primer`, `expand-primers` or by editing the list.
- The sections are `drop` (drops that cost no coverage on any target), `add`,
  `swap` and `trade_off_drop`. `--max-edits N` (default 5) bounds EACH section,
  and each prints "shown M of N considered". Within a section the order is the
  gain in the WORST target's coverage, then the worst host site density, then
  the number of oligos changed. An add is listed only if it raises the worst
  target; a candidate that raises the pooled coverage and not the worst target
  is not offered.
- `trade_off_drop` is not a list of improvements. An entry says that removing
  the oligo raises the worst target-against-host selectivity density from X to
  Y and what coverage that costs on each target. The member of any set with
  the lowest target-to-host ratio always qualifies, so an ordinary set has
  entries here. No threshold decides when a host cost is too high.
- Every entry in every section was evaluated alone. Two drops that each cost
  nothing are not shown to cost nothing together.
- A proposed add is screened against every oligo that stays at `max_dimer_bp`
  (and `max_dimer_dg` when set) and against `min_tm`/`max_tm`, and its sites
  must be known on every target and host. A candidate absent from a host index
  that is not scanned is counted in `candidate_pool.unmeasured` with the reason
  and is not listed. The panel limits in params.json are advisory here: a
  proposal that misses one is listed with the limit named, and one whose limits
  could not be evaluated says so. They are judged on the evaluation's geometry
  (circular unless `--linear`), not on `fg_circular`.
- A failure of `improve-set` leaves no `design_failure.json`, and it removes
  none: it is a report on a design, not a design run, so it cannot block or
  unblock `export`.
- Adds and swaps need `-j`: the pool is `data_dir/step3_df.csv`, opened through
  `candidate_source.open_design_source`. With `--genome` alone only the
  diagnosis and the drops are produced, under the default limits.
- `fixed_oligos` in params.json marks oligos that are never proposed for
  dropping or swapping out. There is no flag for it.
- An oligo whose sites on a target could not be established is reported as
  unavailable there, and while any target is unmeasured no edit is proposed.
  The predicted figures are the evaluation run on the edited set, so they
  agree with a later `evaluate-set` run on the same params file, references
  and geometry flags by construction; that agreement is not an independent
  check.
- `evaluate-set -j` measures coverage at `coverage_reach` from params.json when
  the key is set, as `optimize`, `plan-pool`, `expand-primers` and
  `improve-set` do, and reports `extension_reach_source` beside
  `extension_reach_bp`. Until 2026-10-02 it ignored the key and reported at the
  polymerase's reach. With `--genome` alone it uses the polymerase's reach.
- `design --multi-genome` refuses: it called a pan-genome entry point that has
  never existed in this package. Several targets go in params.json as
  `fg_genomes` with one `fg_prefixes` entry each. `--min-coverage` was removed
  from `design` with it; it fed the same absent code.

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
