# Implementation notes: modules, artifacts, commands and parameters

> Moved from `CLAUDE.md` on 2026-10-01, verbatim apart from two corrections
> marked with that date, when that file was reduced to a
> working summary. Entries are dated records: a figure is what was
> measured on the date given, on the references named, and a later entry may
> correct an earlier one. `CLAUDE.md` keeps the short form and links here.

What each module, output file, command and params.json key does that its name
does not say, and why it is the way it is. For the generated key-by-key
listing see [params-reference.md](../params-reference.md); for the rules a
design run enforces see [design-contracts.md](design-contracts.md); for the
numbered defect record see [KNOWN_ISSUES.md](KNOWN_ISSUES.md).

## Architecture

### Entry Point

- `neoswga/cli_unified.py`: Main CLI entry point (all commands)
- Entry point defined in `pyproject.toml`: `neoswga = neoswga.cli_unified:main`

### Core Modules (`neoswga/core/`)

About 140 modules (143 files on 2026-10-01; the count of 75 recorded here
earlier was stale). `ls neoswga/core/` and the module docstrings are the current list;
what follows is only what the filenames do not tell you.

- **Optimizers** are dispatched through `unified_optimizer.py`:
  `hybrid_optimizer`, `dominating_set_adapter` + `dominating_set_optimizer`,
  `network_optimizer`, `background_aware_optimizer`, with
  `minimal_primer_selector` as a post-process. The registered
  `background-aware` method is `BackgroundAwareBaseOptimizer`, which delegates
  to `HybridOptimizer`. The standalone three-stage `BackgroundAwareOptimizer`
  and its module-level `optimize()` / `compare_optimizers()` were deleted on
  2026-09-10: nothing dispatched to them, and its `_prune_background` had
  diverged from the one that ships.
- **`core/exceptions.py`** holds `StepValidationResult` and
  `StepPrerequisiteError`. `core/pipeline.py` re-exports both, so existing
  importers are unaffected and the re-exported objects are identical, which is
  what keeps `except` clauses matching. They moved because `cli_unified.py` and
  `cli/pipeline.py` imported `core/pipeline.py` at module scope only to make
  the exception catchable, and that import reaches scikit-learn through
  `rf_preprocessing`. Keep this module free of dependencies beyond `typing`
  and `dataclasses`.
- **`occupancy_coverage.py`**: occupancy-weighted coverage, accumulated over
  window edges rather than over bases. Extracted from `base_optimizer` on
  2026-09-17 when the rewrite pushed that module past its size budget;
  `BaseOptimizer._compute_effective_coverage` delegates to it and supplies the
  reach and geometry. The old loop made two full passes over the target per
  primer, so it cost the same on a 1.27 Mb genome whether a primer had two
  sites or two thousand. That was 95% of one objective evaluation while the
  144 Mb host everyone blamed was 16%. 19x to 21x faster, agreeing with an
  independent float64 oracle to 1e-12; the old float32 accumulation is why the
  delivered coverage moves by up to 2.4e-8. Neither this nor `_union_coverage`
  confines a window to the record holding its site, which is Phase 6's subject
  and is deliberately unchanged here. `coverage.merged_window_intervals` is the
  interval form of `_mark_window` and is tested against it base by base.
- **`selectivity.py`**: `selectivity_from_loads`, `selectivity_density_from_loads`
  and the `MAX_SELECTIVITY` / `SELECTIVITY_REFERENCE` constants. Extracted from
  `base_optimizer` on 2026-09-30, the same way and for the same reason as
  `occupancy_coverage.py`: one new configuration field pushed that module past
  its budget, and it had been sitting at exactly its required headroom, so any
  single line would have. All four names are re-exported from `base_optimizer`
  and are the same objects, which is what leaves two modules and six test files
  untouched; the underscored spellings are kept as aliases because that is what
  those importers ask for. The two functions belong together because the
  difference between them IS Known Issue 6.
- **`lazy_dimer.py`**: owns the one decision about how to screen dimers.
  `dimer_screen(pool, max_dimer_bp)` returns the dense `dimer_matrix` below
  `LAZY_DIMER_POOL_THRESHOLD` (4,000) candidates and the pairwise
  `LazyDimerCompatibility` above it, and all three searches ask it:
  `dominating_set_optimizer`, `network_optimizer` and `refine_hybrid_stage2`.
  The last two built the dense array unconditionally, and `plan-pool` sets
  `refinement_method="swap"`, so that was on the hot path -- at the 491,836
  candidates `all_qc` retains the array is about 242 GB, a MemoryError rather
  than a slowdown (audit finding F11). A threshold the dense form cannot
  represent, 8 or above in its 4**8 code space, also takes the pairwise branch
  rather than raising, so what was configured is always enforced;
  `tests/test_one_dimer_screen_for_every_pool_size.py` holds a shrinking
  allowlist of the sites that legitimately build a dense matrix.
- **`candidate_source.py`**: where a command's candidates come from, and in
  what order. (Corrected 2026-10-01: `open_source_or_list`, described next, was
  replaced on 2026-09-21 by `open_design_source`, which refuses a mismatched
  inventory instead of falling back; see
  [design-contracts.md](design-contracts.md#which-candidate-pool-a-command-searches).)
  `open_source_or_list` is the one rule all three commands ask:
  the inventory when the directory has one, the supplied list otherwise, with
  the frontier opening at the list's own size so no delivered panel moves.
  `plan-pool`, `optimize` and `expand-primers` each read `step3_df.csv` for
  themselves before Phase 4 (audit finding F1), which made everything the
  inventory retained beyond the `max_primer` shortlist unreachable.
  `order_candidates_by_background` lives here too, because ordering the scan is
  the same concern as choosing it.
  The frontier opens at the supplied list's size and only `plan-pool` reaches
  past it: `pool_planner` is the sole caller of `source.advance()`. That is
  deliberate as of 2026-09-19. Handing `optimize` the whole 20,670-candidate
  Wolbachia inventory instead of its 2,000-primer shortlist moves the delivered
  panel a long way (Jaccard 0.500/0.263/0.171 at n=6/12/24) and trades
  specificity for coverage: density falls 12-42% while coverage rises 1-6
  points. `max_primer` cuts on `bg_count / fg_count` ascending, so the
  candidates a refill reaches bind the host twice as often and the target half
  as often, and Stage 1's greedy has no specificity term to resist them (Known
  Issue 16). Do not give `optimize` an unconditional refill; make Stage 1
  specificity-aware first, or trigger a refill only on an unmet configured
  panel limit, the way a size row does
  ([measurement](../validation/frontier_refill_on_optimize_2026-09-19.md)).
- **`position_cache.py`**: in-memory binding-position cache, about 1000x faster
  than re-reading the HDF5 files. The constructor takes a fixed primer list;
  `load` and `release` move that window afterwards, which is what a frontier
  that advances needs. A released primer is remembered as released and
  `get_positions` raises for it, because an array that is gone reads exactly
  like one that never existed. `has_entry` answers whether the cache holds an
  ANSWER for a primer on a prefix, which is not the same question as whether
  that answer is non-zero; `require_entries` is the one rule both the inventory
  provider and a `--candidates` list check a batch against.
- **`reference_panel_evaluation.py`**: one record per reference genome, plus the
  target-against-host cross table and two reductions, for `evaluate-set`. Added
  2026-10-02 (Phase 2 of the genomic-diversity plan) because every figure the
  pipeline reports over several references is pooled: a primer frequent in one
  strain and absent from another passes the foreground gate like one present in
  both, and a host that is a small share of the summed background bases is
  invisible in `selectivity_ratio` (Known Issue 6, in a second place). Nothing
  here scores, constrains or defines a threshold, and no design stage reads it.
  Each quantity is a `panel_evaluation.Measurement`, so an unmeasured reference
  is `None` with a reason rather than a 0.0 that reads as a clean host, and
  `worst_target_coverage` / `worst_host_selectivity_density` are `None` whenever
  any member is unmeasured -- `host_profile.aggregate_loads` returning 0.0 for
  an empty list is the shape deliberately not copied. Sites come from one of two
  routes and the record says which: `positions` (the index, plus a scan where
  the caller allows it) also yields coverage and gaps, `counts`
  (`kmer_tables.counts_for`, falling back to `query_scan`) yields sites alone and
  reports the positional figures as unavailable WITH that reason. The two agree
  on the site count, verified on the prepared Wolbachia design: 10 candidates
  gave 50 Drosophila sites by the design's own index and 50 by a single pass over
  the 144 Mb FASTA with no table
  (`tests/test_a_host_by_path_agrees_with_the_designs_own_host.py`). A reference
  with length 0 -- an empty or unreadable FASTA -- is unmeasured rather than
  empty, because its length is the denominator of every figure.
- **`variant_table.py`**: the one door into a variant file, added 2026-10-02
  (Phase 4 of the genomic-diversity plan). `open_variants(path,
  reference_fasta)` reads a VCF or BCF through `pysam.VariantFile` -- imported
  lazily, the pattern `bam_coverage._require_pysam` sets -- or a TSV with a
  header of `chrom pos ref alt` and one optional 0/1 column per strain, and
  returns per strain two sorted int64 arrays giving each variant's
  `[start, end)` in CONCATENATED reference coordinates. `pysam.VariantFile`
  appears nowhere else in the package
  (`tests/test_variant_table_is_the_one_door.py`); a second raw open would get
  none of the refusals below. Contig name to record offset comes from
  `reference_layout.read_layout`, not a second mapping.
  Both formats are 1-based, which is stated because they could differ: VCF
  because the specification says so, the TSV reader by choice, so that one
  convention holds for both. Every refusal is a `ReferenceDataError`: a contig
  the FASTA does not hold, a REF allele that disagrees with the FASTA, an
  unsorted or interleaved table, a TSV genotype that is neither 0 nor 1, a TSV
  row with more or fewer columns than its header, a genotype naming an allele
  its row lacks, a breakend ALT, a symbolic ALT without `INFO/END`, a span past
  the end of its record, a file pysam cannot open (including a gzip rather
  than bgzip `.vcf.gz`, which pysam reports as `NotImplementedError` during
  iteration). The REF check is the one that earns its
  keep -- it catches a table made against another assembly AND a 0-based one,
  since a shift of one base disagrees with the sequence at nearly every row.
  Alleles are checked one record at a time against a streamed FASTA, so a
  host-sized reference is never held whole. A symbolic ALT with `INFO/END`
  spans `[POS - 1, END)`; pysam folds a declared END into `record.stop`, and an
  undeclared one is not read, so it is refused rather than read as one base.
  Genotypes are read from the record's own text, because pysam reports an
  allele index the row does not have as None, the same as an uncalled allele.
  The rule: any allele naming an ALT makes a carrier; otherwise any uncalled
  allele (`0/.`, `./0`, `.|0`, `./.`, no GT) makes the strain `unavailable`
  with that row as the reason; only a fully stated reference genotype is a
  non-carrier.
- **`variant_sites.py`**: which binding sites survive in which strain. Each
  site's k bases are looked up with `numpy.searchsorted` against the variant
  starts plus a running maximum of the ends, which is exact for intervals
  longer than one base. Site geometry follows the scanner (`string_search`):
  a circular reference is ONE ring of the concatenated sequence, so a site
  starting in the last k-1 bases continues at base 0, joining the last record
  to the first on a multi-record reference, and its bases are taken modulo
  that length. A position the run's geometry cannot produce (a wrap read as
  linear, an offset past the end) is `not_assessed`, never intact, and makes
  that strain's fraction, coverage and gaps unavailable. All counts are over
  the deduplicated union of the two strand keys, so a palindromic oligo's site
  is one site in the denominator as well as the numerator. A site is intact
  when no variant falls inside it and affected
  otherwise -- binary, and chosen because it needs no model: the shipped
  mismatch model is a uniform 4.0 C per mismatch with status `assumed` in
  `core/registry/model_evidence.json`. An affected site is never weighted,
  scored with `occupancy.mismatch_tm`, or called tolerated, and neither module
  imports an occupancy or mismatch model at all (asserted, not just intended).
  Affected sites are split by the distance of the nearest variant from the
  primer's 3' end at `THREE_PRIME_WINDOW_NT` (5 bases), measured on the strand
  the primer binds: a site is stored at the forward-strand offset of the k-mer
  for either strand, so the 3' terminus is at `pos + k - 1` on the forward
  strand and at `pos` on the reverse. The split is reported and nothing is
  compared against it. Per-strain coverage and gaps come from
  `coverage.compute_per_prefix_coverage` and
  `reference_panel_evaluation.gap_statistics` -- which was renamed from
  `_gap_statistics` for this -- over an `IntactPositions` view that answers as
  a `PositionCache` does with the intact sites only, so a strain carrying no
  variants reproduces the reference figures exactly rather than nearly
  (`tests/test_variant_sites_match_a_mutated_reference.py`). Reductions over
  the strains are `None` when any strain is unavailable.
  Three limits are carried in `LIMITS` and written into the output and the
  printed report, because a reader acting on an intact fraction needs all
  three: (1) a site GAINED in a strain through a variant is invisible, since
  sites are found on the reference and a k-mer a variant creates is never
  looked for; (2) an indel affects every site it overlaps and the coordinate
  shift it causes downstream is ignored, so positions in a carrying strain are
  reference positions and gap lengths there are approximate; (3) the table says
  nothing about sequence absent from the reference.
  The oracle test is what the module rests on: SNPs generated with a seed are
  written both as a VCF and as a mutated FASTA, and the intact-site set from
  the variant route must EQUAL the reference sites still found by scanning the
  mutated FASTA at the same coordinates -- over several seeds, and on a
  two-record reference with a SNP at the last base of one record and the first
  of the next. Equality, not a tendency.
- **`gpu_acceleration.py`**: CuPy-based thermodynamics helpers. Not reached by
  any pipeline stage, and `--use-gpu` says so rather than claiming otherwise.
  `batch_binding_probability` is vectorised; `batch_calculate_tm` loops in
  Python writing element-by-element into a CuPy array, which is slower than the
  NumPy path it replaces.
- **`advanced_features.py`** and **`gc_adaptive_strategy.py`** read as optional
  but are wired into the kept paths (`rf_preprocessing.py`, and `pipeline.py` /
  `multi_genome_pipeline.py` respectively).
- **`report/`** builds the quality reports. `report/metrics.py` reads
  `effective_conditions` from the run manifest in preference to params.json.
- Pre-trained scorer: `neoswga/core/models/random_forest_filter.skops` (skops
  format, with a SHA-256 allowlist in `models/checksums.json`).

### Data Flow

```
count-kmers            filter                 prepare-candidates     optimize
     |                    |                     |                       |
     v                    v                     v                       v
 *_Xmer_all.txt  -->  step2_df.csv +    -->  step3_df.csv      -->  step4_improved_df.csv
 (k-mer counts)       positions.h5           (ordered candidates)   (final primer sets)
```

**File outputs** (in `data_dir`):
- `step2_df.csv`: Filtered primers with fg_freq, bg_freq, gini, Tm
- `filter_stats.json`: Per-stage filtering funnel counts (rendered in reports).
  The stages are `total_kmers`, `after_fg_frequency`, `after_bg_frequency`,
  `after_thermodynamic`, `after_exclusion_blacklist` (written only when an
  exclusion genome or blacklist is configured), `after_gini`,
  `after_max_primer_cut` and `final_candidates`. The frequency row is split
  because the background gate is the `bg_bool` term, not the later stage that
  used to be labelled "After background/blacklist": that one sat between two
  configuration-gated blocks and was equal to `after_thermodynamic` on every
  run without a blacklist. `after_max_primer_cut` is usually the largest single
  reduction and used to have no stage name at all. On the bundled plasmid
  example the background gate removed 13,006 of 24,809 candidates and the
  `max_primer` cut removed 4,686 of the 5,186 that reached it. Directories
  written before 2026-09-10 carry the old `after_frequency` / `after_background`
  keys and still render.
- `step3_df.csv`: The candidate pool the optimizer reads, carrying the step-2
  measurements in a deterministic order. Step 2's own ranking leads that order
  (`step2_rank`, read off step2_df.csv's row order), with Gini demoted to a
  tie-break and the primer sequence last for totality. It no longer holds an
  amplification score -- see **The `prepare-candidates` stage** below.
- `step4_improved_df.csv`: Final optimized primer sets with enrichment scores
- `step4_improved_df_summary.json`: Authoritative optimizer metrics the report reads (coverage, effective_fg_coverage, selectivity_ratio, selectivity_density, fg_total_length/bg_total_length,
  effective_fg_sites/effective_bg_sites, selectivity_mode, ensemble_comparison, per_target_coverage, strand metrics). `metrics.strand_stats` holds all five strand figures per genome, foreground and host, keyed by prefix; `metrics.primer_occupancy` holds how much of the time each delivered primer is bound, empty when no conditions were attached; `panel_regime` holds which criterion limited the panel and which had no reference.
  Also `unindexed_candidates`: how many candidates the foreground position
  index could not place. Those cover nothing and so are invisible to
  selection; the pipeline path refuses rather than reporting a coverage
  figure that describes only the rest of the pool.
- `*_positions.h5`: HDF5 files with primer binding positions, one per prefix
  and k. Since 2026-09-25 they are written as **sorted blocks**: the k-mers
  sorted as bytes, an offsets array and every position concatenated, with
  entry i at `positions[offsets[i]:offsets[i+1]]`. The earlier layout used one
  HDF5 dataset per k-mer, and on the wMel index 376 MB of its 381 MB was HDF5
  bookkeeping. Sorted blocks bring that index to 24.7 MB and a full load
  from 37 s to 5.4 s, entry for entry identical
  ([measurement](../validation/position_index_layout_2026-09-25.md)).

  `core/position_index.py` reads both layouts and writes only the new one.
  An old index is converted on its next write and read in place until then.
  Keys stay as written, never canonicalised, because the reverse strand is
  read from the reverse complement's entry. An empty entry (scanned, occurs
  nowhere) and an absent one (never scanned) stay distinct: `get` returns an
  empty array for the first and None for the second.
  `index_format_version` stays 2, because it records whether the sites are
  join-safe and a layout change alters no site; `index_layout` names the
  layout. A write rewrites the whole file beside the old one and renames it
  into place, holding the old file open read-write so two runs sharing a
  directory still collide loudly (Known Issue 21). An older release reading
  a converted index refuses at step 4. A file an older `filter` wrote into is
  refused as `MixedLayoutError` and rebuilt by the scan.
- `*_{k}mer_all.provenance.json`: A sidecar recording the genome each k-mer
  table was counted from (absolute path, content fingerprint, digest
  algorithm, k). The fingerprint is a **full SHA-256** as of 2026-09-19. It
  used to hash the size plus the first and last 1 MB, so a substitution
  anywhere in the middle of a file over 2 MB left it unchanged and a
  same-length consensus or sample-specific assembly reused the previous
  genome's counts, index and inventory silently (audit finding F6). The stated
  reason was cost and measurement does not support it: SHA-256 runs at about
  2.5 GB/s, so hg38 is about a second, cached per input per run rather than
  recomputed once per k.

  `digest_algorithm` is what makes the upgrade safe. A record written under
  the partial hash carries a value that cannot be compared with a full digest,
  so it is UNKNOWN rather than stale: step 1 recounts it once, and step 2
  SKIPS it rather than refusing. Those two must stay distinct. Making
  `_table_is_current` false for such a record without teaching
  `_tables_counted_from_another_genome` the difference made every existing
  data directory fail step 2 with "counted from a different genome", which is
  alarming and untrue; 28 tests caught it. After the one recount the records
  are comparable and the guard is stricter than it has ever been. `count-kmers`
  writes it and reuses a table only when it matches; `filter` checks the same
  record before it starts. Without it, repointing `fg_genomes` at a new assembly
  and skipping `count-kmers` built the design from the previous organism's
  counts. Tables written before the sidecar existed have none, which is treated
  as unknown rather than stale: `count-kmers` recounts them once.
- `run_manifest.json`: One appended entry per step (version, git SHA, seed, input
  checksums, CLI invocation). `resolved_params` is a copy of params.json;
  `effective_conditions` is the reaction the step actually ran under, which is
  not the same thing — `retune_for_polymerase` and the GC-adaptive strategy set
  the polymerase, temperature and additives at run time and never write back to
  the file. `export` and `report` read `effective_conditions` in preference to
  params.json, which is what makes their Tm agree with the optimizer's.

## CLI Commands

Setup, reporting, simulation, multi-genome, coverage-gap and primer-expansion
commands are documented in the `neoswga-cli` skill
(`.claude/skills/neoswga-cli/SKILL.md`), which loads on demand.

### Standard Pipeline
```bash
neoswga count-kmers -j params.json  # Step 1: Generate k-mer counts
neoswga filter -j params.json       # Step 2: Filter candidate primers
neoswga prepare-candidates -j params.json        # Step 3: Prepare the candidate pool
neoswga optimize -j params.json     # Step 4: Find optimal primer sets
```

`optimize` can refuse with a step-4 prerequisite error rather than produce a
set. It does so when `step3_df.csv` is missing or empty, when the position
files are absent, and when the position index covers only part of the
candidate pool. The last case used to return a plausible result: the
unindexed candidates cover nothing, so selection never picks one and the
coverage reported is correct for the smaller panel actually delivered. The
remediation is to re-run `neoswga filter`. A caller that passes its own
candidate list programmatically is not subject to these checks.

### Quality Assurance (`--enable-qa`)

Accepted by every pipeline step; each one routes through
`core/pipeline_qa_integration.py`:

- `filter --enable-qa` runs `apply_post_step2_qa_filter` on step2_df.csv (3'
  stability, dimer-hub degree, integrated quality score), rewrites the CSV with
  a `qa_score` column, writes `qa_report.txt`, and corrects the last stage of
  `filter_stats.json`. A QA pass that rejects every candidate fails the step
  instead of writing an empty pool.
- `prepare-candidates --enable-qa` re-orders step3_df.csv by a `composite_score`. With the
  amplification model retired there is no RF half to blend, so this is the QA
  score alone; pass `--amp-model` to get the 0.7 RF / 0.3 QA blend back. The QA
  scores come from step2_df.csv when `filter --enable-qa` produced them, and
  are computed on the spot otherwise.
- `optimize --enable-qa` drops dimer-hub primers from the candidate pool before
  optimizing. The pre-filter is pairwise, so it costs O(n^2) dimer
  calculations.
- `count-kmers --enable-qa` has no QA hook (there are no candidates yet); the
  step logs that and proceeds.

The flag is per-invocation: it is assigned to `parameter.enable_qa` on every
step, so it cannot carry over to a later step in the same process.

### The `prepare-candidates` stage

**Renamed from `score` on 2026-09-21, with no alias.** The old name described
work the stage stopped doing on 2026-09-05, and an alias would have left it
reachable and in every example someone copies. `neoswga score` now fails, but
NOT with a message naming the new command: argparse rejects it as an invalid
choice and prints all 39 subcommands, among which `prepare-candidates` has no
special standing. This entry claimed otherwise until 2026-09-24. Nothing was
left reachable, which was the point of removing it, but a user whose script
breaks is not told what to use instead -- and the repository's own Nightly E2E
workflow was one of those scripts, running red for three nights on this exact
line. `tests/test_workflows_invoke_real_commands.py` is the check that would
have caught it the day the rename landed. `--fast-score` went with it: it selected the
behaviour that had been the default since the model left the default path, so it
was a published flag that did nothing.

The stage prepares the candidate pool; it does not score it. The bundled random
forest was retired from the default path on 2026-09-05 (audit finding F0).

It was computing a prediction for every candidate and then discarding it. The
`min_amp_pred` gate removed 7 of 1222 candidates on the S. aureus panel and none
at all on E. coli (0 of 449) or M. tuberculosis (0 of 319), because the scores
cluster well above the default threshold of 10.0. And every step-4 consumer
reads only the primer column -- `unified_optimizer.py`, `dominating_set_optimizer.py`,
`background_aware_optimizer.py` and `primer_expansion.py` all call
`step3_df["primer"].tolist()`. Asked whether the score identified good primers,
taking the top half of a pool by it and optimizing over that produced the worst
of five half-pools, behind all three random halves.

The model is also fit to synthetic data generated by a hand-written rule in
`scripts/retrain_rf_model.py`, so it reproduces an opinion rather than measured
amplification, and under the default `fast_score` the delta-G features are zeroed,
making the prediction a pure function of the primer sequence -- blind to the
genome and to the reaction.

What the stage still does: it writes `step3_df.csv`, a required intermediate that
six modules read, carrying the step-2 measurements and the deterministic order
`order_step3_rows` establishes. That order is what makes an unseeded run
reproducible, and it does reach the optimizer, which is order-sensitive.

Retiring it changed no delivered panel: re-running the E. coli design returned an
identical 160-primer set. It costs 0.2 s instead of 3.1 s on 449 candidates and
writes 5 columns instead of 61.

**`--amp-model` restores the old behaviour**, score column and gate included.
`min_amp_pred` without it warns rather than silently doing nothing.

### Optimization Methods
```bash
neoswga optimize -j params.json --optimization-method=hybrid           # default
neoswga optimize -j params.json --optimization-method=dominating-set   # fast graph-based
neoswga optimize -j params.json --optimization-method=background-aware # clinical, host-aware
neoswga optimize -j params.json --optimization-method=network          # Tm-weighted, dimer-screened
neoswga optimize -j params.json --optimization-method=clique           # guaranteed dimer-free set
neoswga optimize -j params.json --optimization-method=ensemble         # run all, keep best
neoswga optimize -j params.json --optimization-method=ensemble --ensemble-combine=union  # re-optimize over pooled primers
```

**Optimization Method Comparison**:

| Method | Speed | Best For | Notes |
|--------|-------|----------|-------|
| `hybrid` | Medium | General use (default) | Combines network + set-cover approaches |
| `dominating-set` | Fast | Large primer pools | Graph-based set cover, ln(n) approximation |
| `background-aware` | Slow | Clinical applications | Three-stage. Adds a host-binding term to Stage 1.5 pruning and to the Stage 2 refinement that chooses the panel. Measured against hg38 on the three GC-tier designs at n=24 and n=36, host sites in the delivered panel fall 7-35% against `hybrid` and coverage falls 0.1-3.1 points. At n=12 on those pools it returns the same panel as `hybrid`: Stage 1 yields only 18-19 primers there, so almost every one carries coverage nothing else supplies and the coverage term decides every removal by itself |
| `clique` | Slow | Sets that must be dimer-free | Max-clique on the compatibility graph (swga 1.0's approach). The only method that GUARANTEES no dimerising pair; the others penalise dimers but can accept one. Pools of ~200 candidates; not in the default ensemble |
| `network` | Medium | Tm-weighted selection | Tm weighted, dimer-screened; stops short rather than relaxing the constraint. The `dimer_penalty` multiplier defaults to 0.0 and only ever downweighted; as of 2026-09-10 this method carries the same hard guard as `dominating-set`, and unlike that one it stops rather than admitting an unscreened primer when the pool is exhausted |
| `ensemble` | Slow | Best-of, unsure which | Runs several methods on one shared cache, keeps the best by application-weighted `normalized_score`, prints a per-method comparison table |

**Ensemble** runs a configurable set of methods (default all four) and keeps
the winner. It builds the `PositionCache` once and re-seeds before each method,
so each method's RUN is reproducible and independent of the order the methods
were listed in. Pick the subset with
`--ensemble-methods hybrid network background-aware`. Selection is by
`normalized_score` (a [0,1] value comparable across optimizers; raw `score` is
NOT comparable), weighted by `--application`, then by the smaller set, then by
method name. Those tie-breaks matter: ties are common, because
`background-aware` wraps the same `HybridOptimizer` that `hybrid` uses and the
two frequently return the identical set. Selection used to be a bare `max()`
over a dict built in `--ensemble-methods` order, so a tie was decided by flag
order while this section claimed order-independence -- reordering the same three
tied methods returned three different winners. The runner-up table is written to
`step4_improved_df_summary.json` as `ensemble_comparison`.
`--ensemble-combine union` additionally re-optimizes over the pooled primers
from all methods (can beat any single method; guarded to never worsen).

**Coverage reach (important):** optimizers SELECT for coverage at the realistic
per-primer reach (`coverage.polymerase_extension_reach('realistic')`, ~3 kb for
phi29) — the same reach the result is scored on — while amplification-network
CONNECTIVITY uses single-molecule processivity (~70 kb). Hybrid/background-aware
thread the realistic `coverage_reach` into Stage-1 set-cover so selection and
the reported `fg_coverage` agree (and ensemble comparisons are fair).

A hybrid run prints **two** coverage figures and they do not match, by
construction. `HybridOptimizer._calculate_coverage` works in bins and is used
for progress reporting and the background-pruning floor; `fg_coverage` is
computed base-by-base and is the authoritative number in
`step4_improved_df_summary.json`. Both are labelled in the output —
`(estimated, binned)` against `(measured)` — so read the measured one.

### Choosing the set size

`num_primers` is the most consequential choice in a design and three tools bear
on it. They answer different questions and two of them stop at 20 primers.

- **`--auto-size`** estimates how many primers reach the `--application`
  profile's target coverage under the configured chemistry. It inverts a
  closed-form saturation curve over genome length, primer length, processivity
  and additive effects. It never reads the candidate pool, never looks at
  background binding, and is clamped to the profile's typical range, at most 20
  primers. It does not weigh specificity, so do not read its answer as the best
  size, only as the size that reaches a coverage target.
- **`--show-frontier`** is the trade-off tool. It builds a coverage against
  fg/bg ratio frontier over the real candidate pool using the binding
  positions, and reports where the application profile lands on it. It
  evaluates 4 to 20 primers, so it cannot describe a 96- or 160-oligo panel.
- **The marginal coverage table** that `optimize` prints needs no flag and has
  no size limit. It measures cumulative foreground coverage as the delivered
  primers are added in order, at the same reach the run was scored on, and
  reports the gain per primer in percentage points:

  ```
      n   coverage   pp/primer
     32      0.627        1.10
     96      0.890        0.41
    160      0.943        0.083
  ```

  A flat `pp/primer` column means more primers buy little coverage. It is
  measured on one delivered set in its delivered order, so each row is a lower
  bound on re-optimizing at that size, and it says nothing about specificity.

The two flags work on `optimize` and on `design`. On the measured sweeps
coverage rises monotonically while selectivity density peaks near n=32 for
*M. tuberculosis* and is already falling by n=32 for *E. coli*, so the coverage
curve alone will not tell you where to stop.

```bash
neoswga optimize -j params.json --auto-size --application clinical
neoswga optimize -j params.json --show-frontier
neoswga design -j params.json --auto-size
```

### Utility Commands
```bash
neoswga plan-pool -j params.json --design-grid grid.json  # design per condition
neoswga validate --quick            # Validate installation
neoswga validate --smoke -j params.json  # Check a config: schema, unknown keys,
                                    # genome files, then all four steps against a
                                    # packaged 6 kb target under your chemistry
neoswga build-filter --genome genome.fna -o ./  # Bloom filter for a large background
neoswga show-presets                # Show reaction condition presets
```

### `evaluate-set` against several references

```bash
neoswga evaluate-set --from-results results/ --set 0 --genome target.fna \
    --background host1.fna host2.fna [--scan-background] -o eval/
```

Three additions of 2026-10-02, all after the fact: no design code changed.

- `--background FASTA...` takes hosts by path. Before it, the background came
  only from `parameter.bg_prefixes` / `bg_genomes`, so no command could score a
  delivered set against a host that was not in the design. The bookkeeping is
  `cli/_common.background_references_from_genomes`, a sibling of
  `bootstrap_params_from_genome` that deliberately does NOT write to the
  `parameter` module: a host supplied for evaluation is not a host the design
  used, and the pooled fields read those globals. Prefix and genome travel
  together to `kmer_tables.counts_for`, which is that function's stated
  requirement.
- `--from-results DIR` with `--set N` reads one delivered set through
  `delivered_set.read_delivered_set`, the reader `export`, `interpret`, `report`
  and `simulate` use. It refuses to be combined with `--primers`.
- `evaluation.json` gains `per_target`, `per_host`, `target_host_pairs`,
  `worst_target_coverage`, `worst_host_selectivity_density` and
  `per_reference_notes`. Every field it carried keeps its name and its
  arithmetic: with no `--background` the pooled figures are computed by the code
  that always computed them, which
  `tests/cli/test_evaluate_set.py::test_the_new_blocks_do_not_change_a_run_without_them`
  pins. With hosts given by path and none configured, the pooled host fields have
  a measurement where they previously had None, and they are withheld (None)
  whenever any of those hosts could not be measured.

Cost: a counted host needs no table and no memory beyond the scan;
`--scan-background` locates the sites instead, which measures host coverage and
holds the reference in memory (332 MB on a 144 Mb reference, 2.3 GB on hg38 --
`docs/validation/query_scan_2026-09-25.md`). Counting is the default for that
reason, and a counted host's coverage is unavailable with that reason rather than
zero.

`rescore-set` lost its hand-rolled copy of the per-prefix coverage loop at the
same time (`cli/iterate.py`): it caught a bare `Exception` per primer and
reported 0.0, which is the silent-zero defect, and it ignored record starts and
the circular flag. It now calls `coverage.compute_per_prefix_coverage` and emits
the same field names.

`design --multi-genome` now refuses. It called
`multi_genome_pipeline.design_pan_genome_primers`, which has never existed in
this package; run on 2026-10-02 it logged "Multi-genome mode enabled for 2
genomes" and then printed an `AttributeError` traceback and exited 1. The refusal
names the supported route (`fg_genomes` in params.json, one `fg_prefixes` entry
each). `design --min-coverage` went with it: it fed the same absent entry point
and nothing else on that path read it.

### `evaluate-set --variants`: a SNP table as the statement of diversity

```bash
neoswga evaluate-set --primers SEQ1 SEQ2 --genome reference.fna \
    --variants strains.vcf -o eval/
```

Added 2026-10-02 (Phase 4 of the genomic-diversity plan). Diversity can be
stated as one FASTA per strain, which the per-reference blocks above already
evaluate, or as variants against one reference, which is what a canonical-SNP
matrix or a VCF from a mapping pipeline holds. This is the second form. It is
an evaluation flag only: no schema key, no default, and no selection stage
reads it. A params.json key, if one is ever wanted, belongs to Phase 6.

- The reading is `core/variant_table.py` plus `core/variant_sites.py`, both
  described under Architecture, including the three limits of the route and the
  oracle test the model rests on.
- `evaluation.json` gains `per_strain` (per strain: intact sites, affected
  sites split 3'-proximal / distal with the window stated, coverage and the
  three gap figures over the intact sites only, and the per-primer intact
  fractions) and `variant_route` (the table's provenance, the reference
  figures, the two reductions, and `limits`). Each entry of the existing
  `primers` list gains `intact_fraction_every_strain` -- the lowest intact
  fraction for that primer over the strains -- and
  `intact_fraction_by_strain`. Without the flag nothing is added and no
  existing field moves.
- A fraction is `None` rather than 0.0 for a primer with no site on the
  reference, and for any primer in a run where some strain is unavailable. Zero
  of zero sites intact is not zero percent intact, and a 0.0 there reads as a
  primer whose sites the variants destroyed.
- One reference per run. A variant table is stated against one assembly. With
  one foreground genome that is the reference; with more than one CONFIGURED,
  readable or not, the run is refused unless `--variants-reference FASTA`
  names one of them. Counting configured rather than readable genomes is
  deliberate: picking the only readable one placed a table silently, and the
  REF check cannot always tell near-identical strains apart. The refusal's
  advice was tested by following it
  (`tests/cli/test_evaluate_set_variants.py::test_two_configured_genomes_are_refused_and_the_advice_works`):
  the first version told a `-j` user to pass `--genome`, which with `-j`
  replaces `fg_genomes`, leaves two prefixes and one FASTA, and was refused
  again with the same advice. That case now says to leave `--genome` off. The
  chosen reference is logged, printed above the per-strain rows, and recorded
  as `variant_route.reference_fasta`.
- When positions come from a position index, `reference_layout.verify_layout`
  checks its stored record starts against the FASTA the table is placed on,
  and the configured length against the FASTA's; either disagreement is
  refused. When this run scanned the FASTA itself the index has no record
  starts and that half of the check is recorded as not run, which is correct:
  those positions come from the same file.
- Review of 2026-10-02 found and closed, each with a test and a mutation that
  fails it: a palindromic oligo counted twice in the denominator (a strain with
  no variants reported 2 of 3 sites intact); a VCF half-call `0/.` read as the
  reference allele; a symbolic `<DEL>` read as its one padding base; a gzip
  (not bgzip) `.vcf.gz` escaping as `NotImplementedError`; wrap-around sites on
  a circular reference tested in flat coordinates, so a SNP in the wrapped part
  was missed; a TSV row longer than its header accepted; a genotype naming an
  allele its row lacks misreported as missing; and the running maximum in the
  mask, which no test exercised. `tests/variant_invariants.py` asserts
  `intact + affected + not_assessed == reference_sites` on every block a test
  builds.
- `evaluate-set` joined `cli/_failure.REPORT_ONLY_COMMANDS` with this flag. It
  reads a params.json to find the reference and now raises `ReferenceDataError`
  for a variant table it refuses, which says nothing about the design whose
  directory `-j` named; without the entry it would leave `design_failure.json`
  beside a finished design and `export` would refuse a panel no run had failed
  on. `improve-set` was on that list for the same reason.
- `run_evaluate_set` was at 199 source lines against a 200-line budget, so the
  background block moved to `_background_totals` and the new logic went into
  `_variant_blocks`, `_variant_reference` and `_attach_per_primer_intact`. No
  ratchet entry was added or raised.

### `improve-set`: from an existing set to proposed edits

```bash
neoswga improve-set -j params.json --primers SEQ1 SEQ2 \
    [--background host.fna ...] [--scan-background] [--max-edits 5] -o improvement/
neoswga improve-set --from-results results/ --set 0 --genome target.fna -o improvement/
```

Added 2026-10-02. The handler is `cli/iterate.run_improve_set` and the logic is
`core/set_improvement.py`. It reports and applies nothing: the only file it
writes is `improvement_report.json`, and no `step4_improved_df.csv` is written
or changed.

The four steps, and what each is built from:

1. **Evaluate.** `reference_panel_evaluation.evaluate_reference_panel`, on the
   `ReferenceSpec` list `cli/evaluate._reference_specs` builds, so the two
   commands agree on what a reference is. Oligos absent from the index are
   found by scanning, and mixed lengths need nothing special (Known Issue 20).
2. **Attribute, per oligo.** Sites on each target and host from the
   evaluation's `per_primer_sites`; marginal coverage per target as the
   evaluation of the set minus the evaluation of the set without the oligo;
   dimer partners in the set through `lazy_dimer.dimer_screen`; Tm from the
   resolved reaction against `min_tm`/`max_tm`. "Sole cover of some region" is
   a marginal coverage above zero, so it is unknown when the marginal is.
3. **Propose**, in four sections (`sections` in the report, in this order).
   `drop`: an oligo whose marginal coverage is exactly zero on every target.
   `add`: a candidate that raises the worst target's coverage. `swap`: a
   candidate that dimerises with exactly one oligo of the set, replacing that
   oligo. `trade_off_drop`: an oligo whose removal raises the worst
   target-against-host selectivity density; the entry carries the density
   before and after and the coverage change per target. No threshold is defined
   for any of these; each rule is a comparison of two measured figures.
4. **Check.** An add is screened against every oligo that stays, at the
   configured `max_dimer_bp` and `max_dimer_dg`, and against the Tm window, and
   its sites must be known on every reference, target and host; a candidate
   failing any of these is counted in `candidate_pool` and not listed. The
   configured panel limits are evaluated through
   `panel_acceptance.enforce_constraints` with no repair budget, against the
   references params.json names, and are advisory: a miss is named on the
   proposal. A set whose limits cannot be evaluated carries `evaluated: false`
   with the reason, not a pass.

Within a section the order is the gain in the worst target's coverage, then the
worst host site density (an unmeasured one sorts after every measured one at
the same gain), then the number of oligos changed. `--max-edits` bounds each
section separately, and each section reports `shown` and `considered`.

Four corrections made on 2026-10-02 after review, each with the case that
showed it:

- **Trade-off drops are not improvements.** The density rule is satisfied by
  the member of any set with the lowest target-to-host ratio. A set of six
  oligos with equal target sites and 5, 5, 5, 5, 5 and 6 host sites was told to
  drop the sixth "for high host load", and following the first-ranked entry
  repeatedly would empty a set. No threshold has a reference, so none was
  added: these entries moved to their own section, the label was removed, the
  wording states only what is measured, and the section carries
  `is_improvement: false` and `lowest_ratio_member_always_qualifies: true` as
  data. Every section carries `each_entry_evaluated_alone: true` and
  `jointly_applicable: false`: two oligos covering the same bases each have
  zero marginal coverage, and dropping both was not evaluated.
- **`--max-edits` is per section.** Applied to one ranked list, five helpful
  candidates filled the report and the oligo that binds nothing was never
  shown. `--max-edits 0` shows nothing and still says how many were considered.
- **A candidate with unknown host binding is not an option.** A candidate
  absent from a configured host index, with no scan, ranked first with an
  unmeasured host load. It is now counted in `candidate_pool.unmeasured`, with
  up to five examples naming the candidate, the reference and the reason. The
  same holds on a target. A host that cannot be measured at all therefore
  blocks every add, and the report says why.
- **A failure leaves no record.** `improve-set` resolves the design request to
  read `fixed_oligos`, so a refused params file reached the command boundary as
  a `DesignError`, which wrote `design_failure.json` into the design's
  `data_dir`; `export` then refused a finished design, and nothing cleared the
  record because clearing is a design step's job. `cli/_failure.py` now holds
  `REPORT_ONLY_COMMANDS`, and `write_failure_artifact` writes nothing for a
  command in it. The error is printed and the exit code is nonzero as before,
  and what the design commands record is unchanged.

The advisory limit check uses the evaluation's geometry: targets circular
unless `--linear`, not params.json's `fg_circular`. The limit evaluator has one
circular flag for every reference, where the evaluation gives hosts
`bg_circular`. When a configured host's geometry differs from the targets' and
the difference can change a figure (a host-coverage limit is set, or the
references are scanned), the limits are reported as not evaluated with that
reason. With the defaults (circular targets, linear host) this is the outcome
for `max_host_coverage` and for any limit under `--scan-background`.

With `max_dimer_dg` set and reaction conditions that carry no temperature,
`set_improvement` raises instead of screening at an assumed 37 C.

**The prediction is the evaluation.** `evaluate_reference_panel` gained a
`sources` argument (`PanelSources`), which is where one evaluation reads
positions, counts and weighted loads. `set_improvement.MemoisedSources` keeps
what was read, so the set, the set without each oligo and the set with each
candidate are all evaluated by the unchanged evaluation code while each
reference is read once. Its keys name the reference in full (prefix, genome
path and length for counts), a primer the held position cache lacks causes a
rebuild, and the weighted load is keyed on the reaction's fingerprint. It is
internal to the module. A proposal's predicted figures are that evaluation on
the edited set, and
`tests/test_improve_set_cli.py::test_the_predicted_figures_are_what_evaluate_set_then_reports`
compares them with `evaluate-set` run in a separate process. The agreement
shows the two commands share one model; it is not an independent check of the
model, and the report carries a note saying so.

**`evaluate-set` and `coverage_reach`.** Both commands resolve the reach
through `coverage.resolve_coverage_reach`. `evaluate-set -j` ignored the key
until 2026-10-02 (an inert key, Known Issue 8): with `coverage_reach: 800` it
reported 3000. It now honours the key when a params file is given and writes
`extension_reach_source` beside `extension_reach_bp`. No file it accepted is
newly refused: a value the reader rejects already fails the schema check in
`validate_params_json_file`, which was confirmed by running `evaluate-set -j`
with 0, -5, "x" and 2.5 before the change.
`tests/test_improve_set_cli.py::test_both_commands_measure_at_the_configured_coverage_reach`
pins the agreement with the key set.

What it does not do, deliberately or for now:

- It does not call an optimizer and proposes single edits only (one drop, one
  add, or one swap). Compound edits were not attempted.
- While any target is unmeasured it proposes nothing, because every ranking is
  by the worst target and that is then unknown. The attribution is still
  reported, and an oligo whose sites could not be established on a target is
  `unavailable` there with the reason.
- The candidate pool is `data_dir/step3_df.csv` through
  `candidate_source.open_design_source`, at the frontier that file names.
  Without `-j` there is no pool, and only the diagnosis and drops are produced
  under the default limits (`settings.source` in the report says which).
- `fixed_oligos` is read from the design request
  (`design_request.design_request_for_run`), where it already existed. It was
  accepted and hashed before this command and had no reader; this is its first.
  `excluded_oligos` is honoured the same way: such a candidate is not offered.
- Each evaluation of a candidate reads no reference again, but it does compute
  coverage on every target, so the cost grows with pool size times the number
  of targets. This has not been timed on a real pool here.
- A candidate's self-dimer is not rechecked here; only its pairing with the
  oligos that stay is.
- `improve-set` has no `--polymerase` flag and `evaluate-set` has one, so the
  two agree when that flag is not used to override the params file.
- Host aggregation other than per-host reporting (the plan's Phase 7
  `worst-case` mode) does not exist yet, so an add is not required to be clean
  against every host; its host binding must be KNOWN on every host, and host
  load then enters the ranking only.
- It has not been run on real data in this repository. The check planned for
  it is that the delivered Wolbachia set 0 is offered no edit that raises the
  worst target without costing density.

`--smoke` takes about 4 s against the packaged plasmid pair and exits non-zero
when the configuration would fail, so it is usable in CI. It resolves the genome
paths in params.json relative to the working directory, exactly as a real run
does: pointing it at `examples/plasmid_example/params.json` from the repository
root correctly reports both FASTAs as missing, because that file names them
relatively.

## Key Parameters (params.json)

**Primer filtering**:
- `min_k`, `max_k`: Primer length range (default: 6-12, use 12-18 for longer primers)
- `min_fg_freq`: Minimum foreground frequency (default: 1e-5)
- `max_bg_freq`: Maximum background frequency (default: 5e-6)
- `max_gini`: Maximum Gini index for binding evenness (default: 0.7,
  re-derived 2026-09-10 against delivered coverage). At 0.6 the gate removed
  primers the optimizer had selected: 8 of 160 delivered on E. coli, 1 of 200 on
  S. aureus, 5 of 36 on M. tuberculosis. The kept pools top out at 0.6877,
  0.6932 and 0.6985, and all three shipped configs already set 0.7. The Gini is
  NaN, and the primer is dropped, below `min_gini_sites` combined binding sites:
  one site gives no gap and two give a single gap, whose Gini is identically
  0.0, the best value available. Before that rule 86% of the shipped chr21 pool
  and 96.2% of the plasmid pool scored 0.0.
- `min_gini_sites`: Minimum recorded binding sites, across both strands, before
  the Gini index counts as a measurement (default: 3, the first count at which
  it can vary). Settable in params.json. The `--min-gini-sites` flag on
  `neoswga filter` is accepted and does nothing (see Known Issue 8); use the
  params.json key until that is wired. Lower
  it to 2 or 1 for a small target, where single-site primers are most of the
  pool: on the shipped plasmid example 10,158 of 10,532 indexed k-mers bind
  exactly once, so the default removes nearly all of them.
- `max_primer`: Primers to keep after filtering (default: 500). It bounds the
  working shortlist written to `step2_df.csv`, not what a design can ever reach
  -- see `candidate_retention`.
- `candidate_retention`: Which candidates a design may ever select, and which
  therefore get a background position index. `all_qc` (default) admits every
  candidate clearing the declared hard gates. `post_gini` also requires the
  evenness gate, as an ADMISSION rule rather than a ranking: a candidate that
  misses it is recorded with an explicit failed assessment naming the gate, and
  is not eligible. Both leave `max_primer` in charge of the shortlist, so the
  optimizer's runtime does not move with this setting.

  The eligible set and the indexed set are the same set in both modes, and a
  test pins that. They were not: `post_gini` indexed 20,670 candidates while
  marking all 491,836 eligible, so a design reaching one of the others would
  have scored it against an absent index and read perfect specificity. The
  mode is part of the admission-policy digest, so switching it opens a new
  generation rather than inheriting the other mode's verdicts.
  On the Wolbachia design the two index 491,836 and 20,670 candidates, costing
  443 MB and 18.9 MB
  ([benchmark](../validation/wolbachia_retention_benchmark_2026-09-16.md)).

  **Retention has not been shown to buy anything.** Measured 2026-09-17 once
  Phase 4 made the retained candidates reachable: at panel size 12 the
  shortlist (2,000), the post-Gini inventory (20,670) and all hard-QC
  candidates (491,836) all reach a selectivity density floor of 60 and all fail
  at 80, and the two larger universes deliver panels agreeing to sixteen
  significant figures on both density and coverage. The 471,166 candidates only
  `all_qc` holds changed nothing and cost 772 s against 41 s, plus 123 s of
  cache build and 785 MB of index against 20 MB. The larger universes also
  report a LOWER density on any row they cannot satisfy, which is the Stage 1
  drift recorded in `docs/validation/violation_magnitude_2026-09-17.md` rather
  than retention's doing, and the two cannot be separated until Stage 1 is
  constraint-aware. One pair, one panel size, so this is "no benefit
  demonstrated", not "no benefit exists"; the default is unchanged
  ([measurement](../validation/retention_changes_no_delivered_panel_2026-09-17.md)).

  The shortlist-only `legacy` mode was removed on 2026-09-16. It gave a
  background index to the 2,000 shortlisted candidates only, so the 489,836
  that cleared hard QC without being shortlisted -- 963,931 of their 979,672
  index entries carry real host sites -- scored against an empty background and
  read as perfectly specific. That is the silent-zero shape of Known Issues 5, 6
  and 13, reached by a fourth route. A config still naming it is refused with a
  message saying what replaced it and why.

**Thermodynamics**:
- `polymerase`: "phi29" (30C), "equiphi29" (42-45C), "bst" (60-65C), "klenow" (25-40C)
- `reaction_temp`: Reaction temperature in Celsius
- `na_conc`, `mg_conc`: Salt concentrations (mM)
- `dmso_percent`, `betaine_m`, `trehalose_m`: Common additive concentrations
- `ethanol_percent`, `urea_m`, `tmac_m`, `formamide_percent`: Advanced additives
- `min_tm`, `max_tm`: Melting temperature range

**Polymerase Presets**:

| Polymerase | Temp | Primer Length | Use Case |
|------------|------|---------------|----------|
| `phi29` | 30C | 6-12 bp | Standard SWGA, high processivity |
| `equiphi29` | 42-45C | 12-18 bp | Higher specificity, GC-rich targets |
| `bst` | 60-65C | 15-25 bp | LAMP-like applications, thermostable |
| `klenow` | 25-40C | 8-15 bp | Room temperature, lower processivity |

**Optimization**:
- `optimization_method`: read from params.json since 2026-09-05; it was inert
  before that, and Known Issue 8 records why. An explicit
  `--optimization-method` on the CLI still wins over the configured value, an
  absent flag does not. Values: 'hybrid' (default), 'dominating-set' (fast),
  'background-aware' (clinical), 'network'.
- `num_primers`, `target_set_size`: Requested primer set size (default: 6).
  **It is a request, not a guarantee** (decided 2026-09-14). The delivered panel
  is never larger, and may be smaller for two benign reasons before any pool
  deficiency: Stage 1 stops once the coverage target is met, and selection stops
  rather than admitting a pair above `max_dimer_bp`. The second is usually the
  binding one on a real pool -- measured at `max_dimer_bp` 3 the shipped pools
  support 29, 31 and 26 primers against panels of 200, 160 and 36. A short panel
  is reported with the reason; `--allow-dimer-relaxation` trades the dimer
  constraint for panel size. Guarded by
  `tests/test_delivered_panel_honours_the_dimer_limit.py`.
- `max_dimer_bp`: Longest complementary run tolerated between two different
  primers (default 3, maximum 7). The screen represents t-mers in a 4**8 code
  space, so 8 and above cannot be enforced and are refused by the schema rather
  than silently disabling the screen. A pool supports a bounded panel size at a
  given threshold: measured on the shipped pools, 3 supports 29, 31 and 26
  primers for S. aureus, E. coli and M. tuberculosis, and 4 supports 83, 72 and
  55. The shipped panels are larger than that. Selection therefore STOPS at the
  conforming size rather than growing the panel, because `num_primers` is a
  request; the 11 bp delivered heterodimer against a configured 3 came from the
  relaxation that used to be on by default.
- `allow_dimer_relaxation`: Let selection exceed `max_dimer_bp` when it stalls,
  instead of stopping (default false; `--allow-dimer-relaxation` on `optimize`).
  It trades the dimer constraint for panel size: on a 40-candidate fixture a
  request for 20 returns 12 primers with no violating pair when false, and 20
  primers with 25 violating pairs when true. Every admission is warned about by
  name. `clique` remains strict either way.
- `objective_scan_width`: How many panels the swap repair scores with the full
  objective per round (default 64; None restores an unbounded scan). The scan
  is over candidates times panel, so a 2,000-candidate shortlist against a
  12-primer panel is 24,000 pairs. The cheap bin gain ranks them and only the
  leaders are scored, and the prescreen is not a new criterion -- it is what
  `refine_by_swaps` has always used when given no objective. Given a budget it
  cannot exhaust, widths 16 and 64 and an unbounded scan converge to the
  IDENTICAL panel on the Wolbachia pool at Jaccard 1.000, costing 64, 320 and
  17,913 objective evaluations. The stronger reason to ship a width is not the
  speed: without one the default budget truncated every size measured, landing
  at Jaccard 0.500 against that optimum, so the answer was wherever the budget
  ran out. The cheap pass is deliberately NOT charged against
  `swap_max_evaluations`, since charging it would rank a prefix of the pool and
  reintroduce the blindness the bound removes. Where no constraint binds the
  repair never runs and every width returns the same panel.
  ([measurement](../validation/scan_width_2026-09-17.md))
- `max_frontier_refills`: How many times a size row may widen the candidate
  frontier when it cannot satisfy its constraints (default 4; 0 restores the
  single-frontier behaviour). The inventory holds every candidate that cleared
  hard QC, 20,670 on the Wolbachia design against a 2,000 shortlist, and
  `advance()` returned False from the day it was written, so the rest could not
  affect any panel. Each refill doubles the frontier, so four reach that whole
  universe, and a row that already qualifies never refills: at floors of 40 and
  60 the delivered panel and the runtime are unchanged. At a floor of 100, which
  the shortlist cannot reach, the run examines all 20,670 and reports
  `inventory_exhausted` rather than `frontier_exhausted` after looking at under
  a tenth of what it was allowed to reach. The row carries `frontier_refills`
  and `candidates_exhausted`; the widened frontier is vetted through increment
  3's position check rather than assumed.
  ([measurement](../validation/frontier_refill_2026-09-17.md), which also
  records a pre-existing objective defect this makes reachable: two panels
  failing the same single constraint tie on violation COUNT, so coverage breaks
  the tie and the deciding metric drifts the wrong way.)
- **How failing panels are ranked**: `PoolObjective.shortfall`, not the NUMBER
  of violated constraints. Both objective-scored searches used
  `len(violations)`, so two panels failing the same single limit tied and
  coverage broke the tie, letting the deciding metric drift away from the limit
  it was chasing. Each shortfall term is relative to its own limit so a density
  floor and a site ceiling are comparable, terms sum, and it is zero exactly
  when `violations` is empty -- which is what keeps every feasible panel ahead
  of every infeasible one. A repair that does NOT succeed now returns the panel
  it was given, which is what makes the ordering safe: on a limit no panel can
  meet, chasing it would otherwise trade real coverage for a step toward a floor
  it never reaches. Measured on the Wolbachia pool at an unreachable floor,
  delivered density rose 20.9 to 28.8 and 14.6 to 19.2 for 2 points of coverage
  ([measurement](../validation/violation_magnitude_2026-09-17.md)). Not fixed:
  density still falls as the frontier refills, and that drift is the optimizer's
  own selection rather than the repair's.
- `max_dimer_dg`: Optional ADDITIONAL dimer floor in kcal/mol on the free
  energy of the longest complementary region between two primers, evaluated at
  the reaction temperature. Unset by default. Applied only to a pair
  `max_dimer_bp` has already passed, so it can make the screen stricter and
  never looser, and a configured floor forces the pairwise screen because the
  dense matrix codes t-mers and cannot express free energy.

  **Do not read it as a way to relax `max_dimer_bp`.** A -6 floor with no
  length cap admits 8 bp complementary runs, and the 11 bp delivered
  heterodimer this project recorded is what that looks like. Its use is the
  opposite: raise `max_dimer_bp` for a larger panel and keep a stability bound.
  Measured on 200-primer pools, `run <= 3` supports a greedy panel of 17 to 21
  while `run <= 5` with a -4 floor supports 50 to 77. At the shipped default a
  floor decides nothing at all, because every pair it rejects the run screen
  already rejects. Cost is 1.1x the run screen, not the O(n^2) problem the
  audit guessed. -6.0 follows Rychlik (1995); nothing validates it against a
  reaction
  ([measurement](../validation/dimer_stability_floor_2026-09-18.md)).

- `min_per_target_coverage`: Multi-genome runs only. Minimum coverage required
  on EVERY individual target, unset by default; 0.0 also means disabled. Set it
  and `optimize` prints a per-target table naming the starved targets.
  Aggregate coverage hides them: a panel covering one target 0.9 and another
  0.1 beats a balanced 0.5/0.5 panel on the mean, and nothing in selection
  balances across targets.

  **Checked and reported, deliberately not repaired.** The repair scores
  candidate panels through `compute_metrics`, which does not populate
  `per_target_coverage` -- that is filled in by the caller so all methods get
  it uniformly -- so a floor chased through the repair would score every
  candidate against an empty dict. Same reason `pool_planner.repair_panel`
  leaves a dimer violation alone. `--min-per-target-coverage` previously
  carried an argparse default of 0.0 and now uses the `None` sentinel, so a
  configured value is not beaten on every run.

- `max_sets`: How many distinct primer sets to offer, best first (default: 5).
  Alternatives are found by excluding the primers already chosen and selecting
  again, so each is a different set rather than a reordering. They are numbered
  in the `set_index` column of `step4_improved_df.csv`; set 0 is the one the
  metrics and the summary describe. Fewer than `max_sets` is normal on a small
  candidate pool.
- `iterations`: How many attempts to make when searching for those alternatives
  (default: 8). It deliberately does NOT bound the primary selection — doing so
  would cap how many primers a run can choose, so `iterations: 8` would quietly
  truncate a 96-oligo panel.

**Panel limits** (params.json only; every one unset by default):

| Key | Holds | Needs a background |
|---|---|---|
| `min_selectivity_density` | occupancy-weighted fg load per base over bg load per base, at least | yes |
| `max_background_sites` | total host binding sites, at most | yes |
| `max_worst_hole` | largest foreground gap in bp (`max_gap`), at most | no |
| `max_mean_gap` | mean foreground gap in bp, at most | no |
| `max_evenness` | Gini of the PANEL's foreground gaps, at most (distinct from `max_gini`, which gates candidates) | no |
| `max_host_coverage` | fraction of the host within reach of a panel site (`bg_coverage`), at most | yes |

`core/panel_acceptance.py` reads them. **Set none and nothing changes**:
`constraints_from_parameter` returns None, no objective is built, and the
delivered panel is byte-identical to what it was. That is deliberate rather
than cautious -- no spacing threshold derived from the polymerase reach
separates the 18 published sets with wet-lab outcomes, the winners included, so
NeoSWGA must not pick one, and a fitted weight is wrong for one of the two
benchmarks either way. A user drawing a line is a different claim, and the
"What limits this panel" report is what tells them which properties had no
reference at all.

Set one and `optimize` prints a "Configured limits" table, attempts ONE bounded
repair through `pool_planner.repair_panel` (the same repair `plan-pool` uses, so
there is one in the codebase rather than two that can disagree), and reports
whether it succeeded. A repair that does not resolve the violation returns the
panel it was given: on a limit no panel can meet, chasing it trades real
coverage for a step toward a limit it never reaches. A background-measured limit
set without a background genome is refused rather than reported as satisfied.

`strand_coverage_ratio` and `strand_alternation_score` are deliberately NOT
constrainable. Both read 0.0 when measured zero and when the position cache
could not supply them, and nothing distinguishes the two, so a limit would
reject a panel for a missing measurement while reporting a violated constraint.
The dimer limit is outside for a different reason: it is a hard constraint on
the delivered panel, not a tradeable term.

**Application profiles** (`--application`). **On the default `hybrid` method
this changes nothing about what is selected** -- measured 2026-09-21, all four
profiles deliver an identical panel. The profile sets `tm_weight` and
`uniformity_weight`, both of which `HybridOptimizer` hands to a
`NetworkOptimizer` that nothing ever reads back; Stage 2 is
`_network_refine`, a method on the class itself, with no Tm or uniformity
term. The comment beside that construction claimed it was "the object that
performs refinement" and has been corrected. Known Issue 8's class again, in
`attach_search_config`'s shape: both ends exist and the path does not.

What the profile DOES still do: weight `normalized_score` when picking an
ensemble winner, and steer `--auto-size`. `network` reads `tm_weight`
properly. Setting either weight now warns
([measurement](../validation/selection_weights_are_inert_2026-09-21.md)).

Wiring them into `_network_refine` would move every delivered panel and is a
decision, not a repair -- more so because the Tm term is a Gaussian peaked at
`reaction_temp + 5` while `occupancy.site_occupancy` is monotone increasing in
Tm, so connecting it silently picks one of two unreconciled models. Its span
across oligo lengths is large: median `tm_score` on the plasmid pool runs
0.000014 at k=7 to 0.697 at k=9, 48,488-fold.

The table below is the weighting used to pick an ensemble winner:

| Application | Coverage Target | Specificity | Typical Size | Use Case |
|-------------|-----------------|-------------|--------------|----------|
| `discovery` | 90% | 60% | 10-15 | Pathogen discovery, maximize sensitivity |
| `clinical` | 70% | 90% | 6-10 | Diagnostics, minimize false positives |
| `enrichment` | 80% | 75% | 8-12 | Sequencing enrichment, balanced |
| `metagenomics` | 95% | 50% | 15-20 | Capture diversity |
