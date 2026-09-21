# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

NeoSWGA is a command-line tool for selecting primer sets for selective whole-genome amplification (SWGA). See [README.md](README.md) for user-facing documentation and quick start.

**External dependency**: Jellyfish k-mer counter must be in PATH.

## Architecture

### Entry Point

- `neoswga/cli_unified.py`: Main CLI entry point (all commands)
- Entry point defined in `pyproject.toml`: `neoswga = neoswga.cli_unified:main`

### Core Modules (`neoswga/core/`)

75 modules. `ls neoswga/core/` and the module docstrings are the current list;
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
  what order. `open_source_or_list` is the one rule all three commands ask:
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
  ([measurement](docs/validation/frontier_refill_on_optimize_2026-09-19.md)).
- **`position_cache.py`**: in-memory binding-position cache, about 1000x faster
  than re-reading the HDF5 files. The constructor takes a fixed primer list;
  `load` and `release` move that window afterwards, which is what a frontier
  that advances needs. A released primer is remembered as released and
  `get_positions` raises for it, because an array that is gone reads exactly
  like one that never existed. `has_entry` answers whether the cache holds an
  ANSWER for a primer on a prefix, which is not the same question as whether
  that answer is non-zero; `require_entries` is the one rule both the inventory
  provider and a `--candidates` list check a batch against.
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
- `*_positions.h5`: HDF5 files with primer binding positions
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
reachable and in every example someone copies. `neoswga score` now fails with a
message naming the new command. `--fast-score` went with it: it selected the
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
neoswga build-filter genome.fna ./  # Build Bloom filter for large background
neoswga show-presets                # Show reaction condition presets
```

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
  ([benchmark](docs/validation/wolbachia_retention_benchmark_2026-09-16.md)).

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
  ([measurement](docs/validation/retention_changes_no_delivered_panel_2026-09-17.md)).

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
  ([measurement](docs/validation/scan_width_2026-09-17.md))
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
  ([measurement](docs/validation/frontier_refill_2026-09-17.md), which also
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
  ([measurement](docs/validation/violation_magnitude_2026-09-17.md)). Not fixed:
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
  ([measurement](docs/validation/dimer_stability_floor_2026-09-18.md)).

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

**Application profiles** (`--application`, and the weighting used to pick an
ensemble winner):

| Application | Coverage Target | Specificity | Typical Size | Use Case |
|-------------|-----------------|-------------|--------------|----------|
| `discovery` | 90% | 60% | 10-15 | Pathogen discovery, maximize sensitivity |
| `clinical` | 70% | 90% | 6-10 | Diagnostics, minimize false positives |
| `enrichment` | 80% | 75% | 8-12 | Sequencing enrichment, balanced |
| `metagenomics` | 95% | 50% | 15-20 | Capture diversity |

## The design-failure contract (2026-09-21)

A required calculation that fails now fails the run. It does not return a
substitute value. Four errors in `core/exceptions.py` carry this, all under
`DesignError`: `InvalidDesignRequest`, `ReferenceDataError`,
`UnsupportedModelError` and `ModelEvaluationError`. `SearchBudgetExhausted` is
deliberately outside the family, because spending an allowance is a recorded
stopping point rather than a failure.

The distinction the family draws is between a measurement and its absence. A
candidate that misses a Tm window has been measured and rejected; a candidate
whose Tm raised has not been measured at all. `DesignError.qc_reason` is always
None, so code asking "was this a QC rejection" gets a definite no.

What changed, and what each substitution used to cost:

| Site | Was | Now |
|---|---|---|
| `thermodynamics.calculate_tm_batch` | NaN for any failure | raises; a non-ACGT base is a named `InvalidSequenceError` with `qc_reason`, which is a QC rejection and stays one |
| `thermodynamic_filter._check_heterodimer_pair` | 0.0 free energy, the most permissive answer the screen has | raises; a positive or infinite duplex energy still returns 0.0, because that is a measurement |
| `coverage.polymerase_extension_reach`, `product_reach` | a default reach for an unknown polymerase | raises `UnsupportedModelError` naming the supported set |
| `coverage._record_starts_for` | None when the getter raised | raises; a cache with NO getter still returns None, which is absence rather than failure |
| `occupancy.discrimination_profile` | skipped a failed primer | raises; a mean over an unknown subset was reported with the authority of a mean over the pool |
| `unified_optimizer` per-target coverage | empty dict | raises; `base_optimizer` gates the floor on a non-empty dict, so a requested `min_per_target_coverage` passed vacuously |
| `unified_optimizer` application weights | debug line, defaults applied | raises; `--application clinical` silently had no effect |
| `unified_optimizer` ensemble winner evaluator | debug line, `optimizer=None` | raises; with None every configured panel limit went unenforced |
| `unified_optimizer` post-optimization validator | skipped | raises; skipping disarms the duplicate, size-drift, zero-coverage, blacklist and delivered-dimer checks at once, and writes no validation file, so `export` prints "ready for ordering" |
| `_reseed` | `pass` | raises; the caller logged "set for reproducibility" either way |
| ensemble member failure | any exception became an `status: "error"` row | a `DesignError` propagates, because the next member computes the same quantity from the same data; an algorithm that cannot run on this pool is still a visible row |

`core/design_result.py` separates three things that were one. **Run state** is
what happened to the process (`finished`, `failed`, `interrupted`). **Termination
reason** is why the search stopped (`qualified`, `budget_exhausted`,
`candidates_exhausted`, `refill_exhausted`, `error`). **Qualification** is a
property of the panel. `recommendation_allowed(run_state, qualified)` needs both,
and refuses an unknown state rather than defaulting either way.

At the command boundary a `DesignError` prints the stage, field, artifact, model
and input, then writes `design_failure.json` into the run directory and exits
nonzero. The record exists because an output directory holding last week's
`step4_improved_df.csv` reads exactly like one holding this morning's. Each
pipeline step re-raises `DesignError` rather than reducing it to "step N
failed"; `cli/_failure.py` owns the record.

## One resolved design request

`core/design_request.py` resolves a params mapping into a frozen `DesignRequest`
carrying references, chemistry, candidate-source identity, fixed and excluded
oligos, panel limits, size policy, search budgets, seed, model identifiers and
the concentration policy. Nested content is tuples, so a stage cannot append to
a list it was handed.

`default_sources` records, per setting, whether the request supplied it or which
default did. `request_hash` is a SHA-256 over a canonical JSON form, so it is
stable across processes and independent of key order; it is recorded in the run
manifest for `optimize`.

`optimize` resolves the request from the params FILE before the search starts,
not from the `parameter` module: `get_params` runs inside `optimize_step4`, so
at that point every reaction global still holds its default. This is the same
ordering trap `warn_on_condition_drift` documents.

Refusals it makes that used to be silent: unknown and retired keys (a leading
underscore marks a comment and is accepted), non-finite values, `coverage_reach`
of 0, negative budgets, an oligo that is both fixed and excluded, a
background-measured panel limit with no background genome, and an unsupported
polymerase. The explicit zero matters on its own: `design_context_from_params`
used `override or params.get("coverage_reach")`, and 0 is falsy, so it silently
became 3 kb and every coverage figure was reported at a reach the request did
not ask for.

All three design commands resolve it from the same file: `optimize`,
`plan-pool` and `expand-primers`. A command that skipped the gate would accept
what the others refuse, which is how one params file came to mean different
chemistry depending on which command was run.
`tests/test_resolved_design_request.py` walks the call path from each handler
rather than checking that a call appears somewhere in the module.

**Not yet done from the plan's Task 2**: evaluator code still reads `parameter`
globals at run time, and `OptimizationRequest.optimizer` still owns the
scientific settings. The request is a validation gate and a provenance record,
not yet the single channel those settings travel through.

One instance of the mutable-global read was found and removed the hard way. An
index-identity check placed inside `run_optimization` read
`parameter.fg_genomes` and paired it with the prefixes the CALL was given;
under `pytest -n 8` that paired a test's own prefix with another test's FASTA
and the design refused its own index. Reference identity now travels on the
request, where a prefix and its genome are named together.

## Which candidate pool a command searches

`open_source_or_list` is gone, replaced by three functions in
`candidate_source.py` that keep absence and failure apart.

- `open_explicit_source(candidates)`: the pool the user named. It always wins.
- `open_inventory_source(...)`: the inventory, or **None** when the directory has
  none, which is a fact about the directory. It raises `ReferenceDataError` when
  the directory HAS an inventory that holds nothing under this reaction.
- `open_design_source(...)`: the rule every command asks, built from those two.

The old function wrapped the inventory open in `except ValueError` and fell back
to the caller's CSV for both cases. A reaction fingerprint mismatch therefore
became a quiet run over the `max_primer` shortlist, at the shortlist's frontier,
with everything the inventory held unreachable and one `logger.info` line to say
so. That is how the occupancy-gate measurement in Known Issue 17 produced an
apparent density improvement that was not real.

So `filter --preset enhanced_equiphi29` followed by a plain `optimize` now
refuses, naming the remedy, where it used to warn and proceed.
`tests/test_optimize_warns_on_condition_drift.py` was inverted to match.

A self-dimer screen that empties a non-empty frontier now raises
`NoCandidatesError` naming the screen and its threshold, instead of handing an
empty list to an optimizer that answered "candidates list cannot be empty".

## Two scanners, one quantity

Fixed 2026-09-21. `string_search` has two position scanners: the
Aho-Corasick `get_all_positions_multi_k`, and the sliding-window
`get_all_positions_per_k` used when that package is absent or one k is being
scanned. Records are concatenated with no separator, so the last k-1 bases of
one record and the first bases of the next form k-mers that occur in neither.
The Aho-Corasick path rejected those matches. The sliding-window path did not.

So the same reference produced different site sets depending on which path
ran, and the sliding-window path stored up to k-1 fabricated sites per record
join. On a two-record fixture `ACGGTA` is absent from both records and was
stored at offset 4. A fabricated foreground site inflates coverage, a
fabricated background site deflates specificity, and nothing downstream can
tell either from a real one.

Single-record references are unaffected: with no joins the two scanners always
agreed, which is why every complete bacterial genome and every plasmid in this
repository reads the same before and after. Draft assemblies and hg38 are where
it bit.

`spans_a_record_join` now holds the rule once and both scanners call it. A
match beginning exactly ON a boundary starts a record and is kept; off by one
here would delete the first k-1 sites of every contig.

`INDEX_FORMAT_VERSION` is 2. A version 1 index is refused for a MULTI-RECORD
reference only, because a single-record one cannot carry the defect and
refusing it would force a recount for something that never applied to it.

`tests/test_positions_agree_with_an_independent_count.py` checks the scan
against a brute-force sliding window written in that file, which calls neither
scanner nor any helper they call. It covers overlapping occurrences,
palindromes, a reverse complement that also occurs forward, ambiguous bases, a
circular origin, record joins, and an exhaustive pass over every window of a
small reference. It also asserts the two scanners agree, which is the check
that would have caught this.

## What the chemistry model supports

`neoswga/core/registry/model_evidence.json` records, per constant, what it is
and over what domain the model supports it. `model_evidence.py` loads it and
`require_model_support(request)` runs inside `resolve_design_request`, so a
computation outside a recorded domain is refused before any index is opened.

Five statuses, and the line that carries the weight is between the first two.

| Status | Meaning |
|---|---|
| `measured` | the cited work reports this value for a case the model applies it to |
| `estimated` | extrapolated from data at another temperature, on longer DNA, or in another buffer |
| `empirical` | chosen so the model behaves plausibly; no source reports it |
| `assumed` | a modelling decision with a stated reason and no measurement |
| `absent` | nothing computes this effect, and the code must not report zero for it |

Most additive coefficients are 37 C figures for PCR-length duplexes applied to
12-mers at 30 C. That may well be fine; it is not a measurement of it. Only
`tm_urea` was chosen because its source concerns short oligos.

**The registry does not claim the literature was re-read.** Every `source` is
the attribution the repository already carried. What the registry adds is the
second judgement the prose ledger never made: whether the cited work covers
this case. A registry that implied verification would break, in the act of
recording it, the rule it exists to enforce.

**What is refused**: an unknown polymerase; an oligo length outside the
enzyme's modelled range, which is the defect where a Bst design was filtered
through phi29's 6-12 bp window; and an additive whose duplex effect nothing
computes while the literature expects one. Estimates are NOT refused. A model
that declined to run on an extrapolated coefficient would decline to run.

Two findings from compiling it, both in
[docs/validation/chemistry_model_evidence.md](docs/validation/chemistry_model_evidence.md):

- **Glycerol is accepted, range-validated, printed, and changes no Tm.**
  Measured here: a 12-mer at 10% glycerol returns the same effective Tm as at
  0%, to the last digit. The literature expects a real destabilisation, and the
  shipped `q_solution` preset sets 10%. No coefficient is invented to close
  this, because inventing one is the promotion of an assumption the registry
  exists to prevent; a design that sets glycerol is refused instead. BSA and
  PEG also have no Tm term and are `assumed` rather than `absent`: they act on
  the enzyme and on crowding, not on duplex stability.
- **`neoswga.core.registry` was not installed.** It was missing from
  `pyproject.toml`'s explicit `packages` list, so a built wheel contained none
  of it, while `core/parameter.py` imports `registry.views` at module scope.
  Verified by building a wheel. Package data could not have helped:
  `include-package-data` applies to packages that are being installed. A test
  now compares the declared list against the packages on disk.

`docs/SCIENCE_CITATIONS.md` had drifted: it states Klenow processivity as
10,000 bp citing Bambara (1978) while the shipped registry says 40 bp. The
prose was right when written and the code moved. The registry is checked
against `registry/views.as_characteristics()` by a test, so the two cannot
disagree silently; read the prose document as commentary rather than as the
record.

## One panel assessment, and a coverage oracle that is not the code

`core/panel_evaluation.py` holds `evaluate_panel(request, oligos, metrics) ->
PanelAssessment`: one immutable record carrying the panel, the request hash,
every metric WITH its units and the reach and denominator it was computed at,
per-target results, every hard-constraint violation, and qualification as a
boolean that is true exactly when there are none.

Three rules it enforces, each for a failure this project has seen the shape of:

- **A non-finite required quantity fails the run.** NaN compares False against
  every threshold, so a panel carrying one passes no limit and fails no limit.
- **Unavailable and zero have different representations.** A panel that binds
  the host nowhere and a panel whose host index was never opened both read zero
  otherwise. `Measurement` refuses to hold both a value and a reason for not
  having one.
- **A verified zero background is a zero denominator, not a ratio.**
  `base_optimizer` reports `MAX_SELECTIVITY` (1e6) there, deliberately, because
  it is finite and JSON carries it; its own docstring says that means "no
  background binding was detected", not "measured this well", which concedes a
  reader cannot tell them apart from the number. The assessment says undefined
  and why, and keeps the site count. The sentinel is left in place because
  changing it moves every saved summary.

It is NOT yet routed through the optimizer, the acceptance path or the report.
Those still assemble their own answers, which is the rest of Task 5.

`tests/test_coverage_independent_oracle.py` checks coverage against a
base-by-base oracle written in that file. It calls neither
`merged_window_intervals` nor `_mark_window` nor anything they call, so
agreement is evidence rather than a restatement. It is deliberately the slow
implementation production replaced: the fast one accumulates log(1 - theta) at
window edges, and an edge-accounting error is invisible from inside that
formulation. Verified load-bearing by mutating the production grouping from
per-primer to per-site, which fails four of its cases.

The window convention, measured rather than assumed: a site at `pos` with reach
`r` covers `[pos - r, pos + r)`, so a window is `2r` wide.

**The two production coverage paths disagree across a record join, measured and
not fixed.** `compute_per_prefix_coverage` marks through `_mark_window` WITH
record starts, so a window stops at a contig edge. `_union_coverage` and the
occupancy path go through `merged_window_intervals`, which takes no record
starts by design. On the fixture in that file a site 2 bases before a join
covers 12 bases confined and 20 unconfined. So a multi-record reference has two
coverage figures and which one a reader sees depends on the code path. The test
asserts both numbers so the gap cannot grow unnoticed.

## The report and the saved result must agree

`tests/test_report_agrees_with_the_saved_result.py` asserts that every quantity
the report renders equals the one the summary holds, for the exact exported
panel. The report computes its own estimates from the results CSV and overrides
them with the optimizer summary where the summary has an opinion, which is the
right order; `from_optimizer` says which a reader is looking at.

**Writing that check found a favourable default.** `mean_gap`, `max_gap`,
`gap_gini` and `gap_entropy` were read with `.get(key, 0.0)`, and zero is the
BEST value for every one of them. A `max_gap` of 0.0 says the panel leaves no
coverage hole anywhere. A summary that did not carry the key therefore rendered
as the best possible measurement rather than as no measurement, and directories
written before these keys existed are explicitly supported, so this was
reachable rather than theoretical.

The four fields are now `Optional[float] = None`, read without a literal
default, and both render sites skip the gap section unless every figure in it
was measured. Rendering it with one missing is what put a "no coverage hole"
verdict beside three real numbers.

This is the silent-zero family again, in its fourth shape: not a scan that
found nothing, an integer that saturated, a cache asked for what it does not
hold, or a guard with the wrong predicate, but a dictionary lookup whose
fallback happens to be the answer everyone wants.

## Two counters, and only one of them bounds the run

`SearchBudget` is the SHARED ledger. It counts uncached evaluations of the
shared objective across every stage and raises when spent, so a later stage
cannot get a fresh allowance. `swap_max_evaluations` is a PER-STAGE allowance
inside the deletion and swap loops, and it starts at zero every time one of
them is entered.

The per-stage one is not a defect; bounding one loop is reasonable. What would
be a defect is believing it bounds the run. It defaults to 10,000 and looks
like a total. The only setting that is a total is `total_search_evaluations`,
and it is None by default, so **by default there is no total bound at all.**

`describe()` now carries `uncounted_scopes`, naming the two kinds of work the
ledger does not see, rather than leaving a reader to infer them from a count
lower than they expected:

- **proposal generation**: an optimizer scoring candidate panels through
  `compute_metrics` directly rather than through the shared objective. Only
  `clique` does this, in a loop bounded by its own `max_scored_sets`.
- **final assessment**: one `compute_metrics` per stage once the panel is
  decided, deliberately uncharged so reporting cannot consume a search's
  allowance.

`tests/test_search_budget_contract.py` holds the ratchet.
`UNCOUNTED_SEARCH_LOOPS` lists every function that evaluates panels in a loop
outside the objective, with its reason and the bound that does apply, and the
list can only shrink. A call made ONCE per stage is final assessment and is
not flagged; one inside a `for` or `while` is search work, and search work the
ledger cannot see is what the check is for. Verified load-bearing by wrapping
an existing single call in a loop, which fails it.

## The smallest pool, and what deletion cannot reach

`--minimize-primers` reduces a panel by removing one primer at a time and
keeping it when no single removal still qualifies. That is a local optimum, not
the smallest pool, and the two differ whenever one candidate covers what two
others cover between them.

`tests/test_smallest_pool_search.py` enumerates every subset of a four-candidate
set-cover fixture and compares the search against the enumerated answer. The
fixture is constructed so the structure is visible: one candidate covers the
union of two others, so `{P1,P2,P3}` qualifies at size 3, no single deletion
from it qualifies, and `{P4,P3}` qualifies at size 2. Deletion stops at 3;
`panel_beam.beam_search` asked for 2 finds `{P4,P3}`.

**The capability exists and `optimize` cannot reach it.** `beam_search` is
called only from `pool_planner`, so `plan-pool` can escape this local optimum
and `optimize --minimize-primers` cannot.

**It is not wired, because the gap did not reproduce on either real instance
available** (measured 2026-09-21). On `examples/plasmid_example` one primer
covers the target completely at 3 kb, so there is nothing to reduce. On a
300 kb random sequence with 60 candidates, deletion stopped at the requested 20
and a beam over the same pool found nothing smaller at the same coverage.
Random sequence rarely produces the dominance structure the fixture has.

So: the local optimum is real, the beam escapes it, and whether it costs
anything on a pool anyone would design from is unmeasured. Same resolution as
Known Issues 11 and 16 -- a mechanism that fires without a demonstrated benefit
does not ship on by default. Measuring it needs a real multi-kb target with a
candidate pool that is not saturated, which this repository does not contain.

The file also states the claim the plan forbids, as arithmetic: a search that
examined every candidate examined `n` things, while the subsets number `2**n`.
Exhausting candidates is not exhausting candidate subsets, and only a toy case
like this one can produce a minimum certificate.

## A failed run leaves nothing exportable

An output directory is the only thing a later command sees, and nothing in it
carries a timestamp anyone compares. A directory whose most recent run FAILED
therefore looked exactly like one whose run succeeded: the previous run's
`step4_improved_df.csv` was still sitting there, real and stale, and `export`
turned it into an oligo order.

Task 1 wrote `design_failure.json` for exactly this, **and nothing read it.**
That is the Known Issue 8 class in artifact form: the evidence exists and the
check does not.

`export.export_is_blocked(results_dir)` is that check, and `neoswga export`
now consults it before loading anything, exiting nonzero with the recorded
stage and reason. It refuses a `failed` or `interrupted` run, and also a
`finished` one that did not qualify, since finishing is not finding something.

A record that cannot be parsed blocks rather than passes. Unknown is not
success, and defaulting the other way would make a corrupted artifact the most
permissive state available, which is every silent-zero in this file.

`cli/_failure.clear_failure_artifact` removes the record when step 4 finishes.
A record that is never cleared is as wrong as one that is never read: the user
fixes the problem, the next run succeeds, and the export refuses on evidence
that no longer describes anything.

Verified end to end: a directory carrying a failure record exits 1 and writes
no FASTA; the same directory with the record cleared exits 0 and writes one.
`tests/test_design_report_provenance.py` also pins that the CSV, the summary
JSON, the rendered report and the exported FASTA all name the same panel.

## Testing

```bash
pytest tests/                         # All unit tests
pytest tests/test_hybrid_optimizer.py  # Specific test
neoswga validate --quick              # Quick validation
```

**Integration tests** (`tests/integration/`):
- `phi29_baseline/`, `phi29_with_bg/`: Phi29 polymerase scenarios (no background / with background)
- `equiphi29_baseline/`: EquiPhi29 scenario
- End-to-end tests: `test_pipeline_e2e.py`, `test_integration.py`, `test_optimizer_method_coverage.py`, etc.

Two things a contributor should know:

- A full suite run leaves `git status --porcelain` byte-identical.
  `tests/test_the_suite_leaves_no_files_behind.py` enforces that, and also
  fails if a root-anchored `.gitignore` entry is added for a pipeline artifact
  instead of stopping the write.
- Tests needing `examples/plasmid_example`'s generated files guard on
  `tests.conftest.plasmid_example_ready()`, not on the directory existing --
  the directory is committed, so its presence proves nothing. Without jellyfish
  those tests skip with a reason naming it. A new test that pins a quantity
  should follow `tests/test_hybrid_optimizer_run.py`'s fixtures, which write
  HDF5 directly and need no external tool.

## Code Patterns

**Parameter handling**:
```python
from neoswga.core.parameter import get_params
params = get_params('params.json')
```

**Multiprocessing**:
```python
from neoswga.core.utility import create_pool
with create_pool(cpus) as pool:
    results = pool.map(process_func, items)
```

**Position data** (HDF5 format):
```python
import h5py
with h5py.File('positions.h5', 'r') as f:
    positions = f[primer_sequence][:]
```

## Known Issues

1. **Large background genomes**: Use `neoswga build-filter` to pre-build a Bloom filter for human genome. Note this is not always needed: at k=12 exact jellyfish counting of the whole human genome costs about 7 minutes and a 138 MB table (8,368,418 canonical 12-mers), which is well within reach. The Bloom path matters at longer k, where the count table stops being small. Measure before reaching for it -- the sampled-index path has its own resolution trap (see `_warn_if_sample_too_sparse`).

2. **sklearn compatibility**: The RF model ships in skops format (version-tolerant, no arbitrary-code deserialization), so minor sklearn upgrades no longer require retraining. A major sklearn upgrade may still warrant re-validating the model.

   **The format is version-tolerant; its default trust list is not**, and the
   two are easy to conflate. skops 0.15.0 stopped implicitly trusting
   `sklearn.tree._tree.Tree`, and every model-loading test went red on CI while
   passing locally on 0.14.0 -- `pip install -e ".[dev]"` resolves
   `skops>=0.11,<1` to whatever is newest, and `requirements-dev.lock` is not
   used by the test job.

   `rf_preprocessing._TRUSTED_MODEL_TYPES` now names the types a forest
   legitimately needs and `unexpected_model_types` refuses anything else,
   naming it. That is narrower than trusting the archive wholesale and does not
   depend on a future release keeping today's defaults; pinning skops would
   have worked until the next release did the same thing. The digest check
   against `models/checksums.json` still runs first, so this narrows what an
   already-vouched-for file may reconstruct rather than replacing provenance.

   The refusal rule is a pure function because which types skops reports as
   untrusted depends on the installed version: a test driving the loader would
   exercise it on 0.15 and skip straight past it on 0.14.

3. **Memory usage**: The filter command loads all background k-mers into memory. Use Bloom filter for large backgrounds.

4. **PositionCache strand parameter**: Uses 'forward', 'reverse', 'both' (not '+' or '-').

5. **pyahocorasick silently finds nothing past 2 Gb** (upstream, unfixed): `Automaton.iter()`
   indexes with a 32-bit int, so for a string longer than `2**31 - 1` its scan loop never runs.
   It yields nothing, raises nothing and warns nothing. Verified on **2.3.1, the latest release
   as of 2026-08**, with the needle planted at offset 1000: length `2**31 - 10` finds it,
   `2**31 + 10` does not. Not found in the upstream tracker.

   The re-check is automated rather than left to memory. `tests/test_pyahocorasick_limit_canary.py`
   fails as soon as the installed version leaves the set that has actually been measured, and
   names the command that re-measures it:

   ```bash
   NEOSWGA_VERIFY_AHOCORASICK_LIMIT=1 pytest tests/test_pyahocorasick_limit_canary.py -k still_present
   ```

   That run allocates over 2 GB, so it is opt-in. If the limit is gone, the chunking is still
   correct but no longer load-bearing, and keeping or dropping it is a deliberate choice.

   `string_search.MAX_SCAN_CHUNK` (2**30) works around this by scanning in overlapping windows,
   so this is **load-bearing, not defensive** -- removing the chunking silently breaks any
   background above 2.147 Gb, which includes human (3.1 Gb) and mouse (2.7 Gb).

   Symptom before the fix: `total_bg_sites` and `bg_coverage` read 0 for a whole-genome host,
   which is indistinguishable downstream from a perfectly specific primer set. A 27-primer panel
   scored 0 against hg38 where the jellyfish counts put the true figure at 860.

6. **A partial background is not a specific design**: `selectivity_ratio` is a ratio of counts
   with no genome length in it, so it moves with how much background sequence you supply --
   about 66x between human chr21 and whole hg38, with nothing about the primers changed. Read
   `selectivity_density` (added beside it) when comparing designs scored against different
   backgrounds. At k=12 note that 99.7% of all canonical 12-mers occur in the human genome, so
   a gate demanding zero host sites returns a single-site candidate pool rather than failing.
   See [docs/validation/additive_specificity.md](docs/validation/additive_specificity.md#the-background-was-one-chromosome-and-that-mattered).

7. **Genome coordinates are int64** (`position_cache.POSITION_DTYPE`), and must stay so.
   Positions were cast to `np.int32` on the way out of HDF5. int32 tops out at 2,147,483,647,
   so for human (3.1 Gb) and mouse (2.7 Gb) every site past that offset **saturated at the
   ceiling**; the `strand="both"` path then calls `np.unique`, which collapsed all of them into
   a single site.

   This is a second, independent instance of the 2**31 failure in Known Issue #5 -- different
   place, same symptom, and neither one is visible on a single chromosome. Measured on an
   *M. tuberculosis*-vs-hg38 run: `CACCGACGACGA` occurs 48 times in hg38 (jellyfish and a
   direct string count agree), the scan stored all 48 correctly, and the cache returned 5.
   Across the twelve-primer set `total_bg_sites` read 48 against a true 114.

   It changed the design, not just the report: with the true background visible the optimizer
   **drops** `CACCGACGACGA`, the primer whose host load int32 had been hiding. After the fix
   `total_bg_sites` matches the jellyfish count exactly (77 = 77).

   The lesson both issues share: **test against a whole genome, not a chromosome.** chr21 is
   46 Mb and cannot reach either limit, so both bugs sat behind a passing test suite.
   `tests/test_position_cache.py::TestPositionsPastTheInt32Ceiling` pins this one.

**`--design-grid` on `plan-pool`** designs once per condition and length in a
JSON grid and writes `design_sweep.json` beside the usual report, rather than
one `pool_plan`. The grid names `lengths` and `conditions`, where each condition
names only the fields it changes: the baseline is the reaction this run
resolved, so a grid varying DMSO alone keeps the buffer, salts and oligo
concentration, and the comparison is between chemistries rather than against
library defaults. A cache and optimizer are rebuilt per condition and length,
since the index is per length and the chemistry is what varies.

It needs the candidate inventory, and it looks each condition up by reaction
fingerprint, so a condition the filter never recorded is reported as having no
eligible candidate rather than designed with an empty pool. Wired on 2026-09-17
in Phase 4 increment 6; it was audit finding F4, parsed and documented and read
by nothing, and its entries are now gone from both the inert-option and
unreachable-capability allowlists.

**`--min-fg-bg-ratio` was read and then overruled** -- FIXED 2026-09-17.
`optimize`'s background prefilter kept every candidate at or above the ratio,
then, if that removed more than `max_removal_fraction` of them, discarded the
threshold and kept the top 80% by ratio instead. On the 2,000-candidate
Wolbachia shortlist the threshold removes 64.8% at its default of 1.0, so the
clause fired at 1.0, 2.0, 5.0 and 20.0 and removed exactly 400 every time. The
flag changed nothing above about 1.0 and the rule in force was "drop the worst
20%".

That is the Known Issue 8 class in a shape none of its ratchets look for: not a
flag nobody reads, but a flag that is read and then overruled by a second rule
on the same decision. `max_removal_fraction` was also a bound on the fraction
of a BATCH, so which candidates survived depended on how many others were below
the threshold alongside them.

`order_candidates_by_background` replaces it. Candidates at or above the ratio
are searched first and the rest are searched last; nothing is deleted, so the
400 the old path made unreachable at every setting are reachable again, which
matters because increment 5's refill can now reach them. The partition is
stable, preserving the inventory's `search_rank` traversal. Delivered panel on
the measured design: 11 of 12 primers shared, Jaccard 0.846
([measurement](docs/validation/background_ordering_2026-09-17.md)).
`bg_max_removal` is retired with the clause.

**The objective never reached the stage that refines** -- FIXED 2026-09-18.
`plan_pool` attached `pool_objective` to the optimizer it was handed, which on
every command-line path is a wrapper (`HybridBaseOptimizer` or
`BackgroundAwareBaseOptimizer`) that delegates the search to an inner
`HybridOptimizer`. `_swap_refine` is a method of the INNER one, so
`refine_hybrid_stage2` read the attribute off an object nobody had set. Measured
through a real design, the refinement ran once and received None: Stage 2
refined on raw covered bases while the row was accepted on occupancy-weighted
coverage under a specificity floor. Delivered density on a failing row went from
28.78 to 42.62 once connected.

Two tests covered it and neither could see it -- one asserted by AST that
`plan_pool` assigns an attribute of that name, the other by source text that the
refinement reads one. Both ends existed and the path did not. Use
`swap_refinement.attach_search_config` for anything a delegate must read, and
assert the PATH: `tests/test_the_objective_reaches_the_stage_that_refines.py`
drives a real factory-built optimizer under both methods.

**Stage 1 is deliberately NOT constraint-aware**, and there is no demonstrated
reason for it to be. The specificity density floor is exactly additive over
primers, so a density-only ceiling is computable and comes to 79.807 for a
12-primer panel on the shipped pool. **That figure bounds nothing deliverable**:
the panel achieving it has coverage 0.4042 against a 0.5 target, and the
constructions that appear to beat the search carry 30 dimerising pairs out of 66
at the configured `max_dimer_bp` of 3. Only 6 of the 16 most selective
candidates are mutually compatible, so a 12-primer panel cannot be built from
them at all.

With dimers and coverage both in force, no deterministic construction beats the
search: the best reach density 50.6 at coverage 0.527, or 64.0 at coverage 0.385
which fails the target, against the search's **60.112 at 0.6535**. So its
failure at a floor of 65 is probably correct. Three Stage 1 rules built on the
accounting each failed to improve a delivered panel, which is best explained by
there being nothing to find.

Settled with controls on 2026-09-18: **no construction respecting the dimer
screen beats the search.** Nine of them, over three reference densities and
three pool sizes, all land between 0.375 and 0.387 coverage against a 0.5 target
and none exceeds 64.0 density. Selective candidates are GC-richer (0.44-0.48
against 0.39) and pairwise less compatible (58-63% against 69-71%), but the
largest mutually compatible subset is NOT smaller -- 18 among the top 64 by
slack, more than a 12-primer panel needs -- so compatibility is not the barrier.
The specificity against coverage trade-off is
([measurement](docs/validation/no_search_headroom_on_this_pool_2026-09-18.md)).

The lesson is the reusable part: **an achievability figure that omits a
constraint bounds nothing, and a striking ratio without a control is not a
finding.** Three claims of mine died in sequence here -- a 79.807 ceiling that
ignored coverage and dimers, an existence proof carrying 30 dimerising pairs out
of 66, and a "6 of 16" compatibility barrier that sits inside the random range
of 6 to 8. Three Stage 1 search rules were built on the first two before the
third was tested. The check that would have killed all three at the outset is
the same one: evaluate a candidate panel through the acceptance path a delivered
panel takes, and compare it against a control. Quote 60.112 at coverage 0.6535
as the reference for this pool. The accounting lives in
`scripts/benchmarking/selectivity_budget.py` as a diagnostic, not in the
package, because nothing in the search uses it.

8. **`optimization_method` in params.json did nothing** — FIXED 2026-09-05
   (audit finding F1b). The key was declared in `params.schema.json`,
   documented above, accepted by the validator, and read by nothing:
   `get_params` assigned no module global, and `run_step4` passed
   `args.optimization_method` straight through with an argparse default of
   `'hybrid'`, so the flag's default beat the config every time.

   ```
   params.json "optimization_method": "dominating-set"
     -> parameter module global: <UNSET>,  optimizer actually run: hybrid
   ```

   It cost more than provenance: `hybrid` returns a set **identical** to
   `dominating-set` (Jaccard 1.000) at 7.8x the cost at 32 primers and 260x at
   128 (2239 s against 8.6 s), so every params.json user ran the slowest method
   for the same answer. `design` pinned it a second way — that subparser has no
   `--optimization-method` and `run_design` hardcoded `"hybrid"`.

   Fixed in three places, because one alone was not enough. The global is
   assigned in `_apply_params_only_keys`. The flag's argparse default is now
   `None`, the sentinel that distinguishes an explicit `--optimization-method
   hybrid` — which must beat a configured `dominating-set` — from an absent
   flag, which must not; do not give it a real default again.
   And the lookup goes through `optimization_method_from_params`, beside the
   other pre-read resolvers, because `run_step4` builds its argument list
   before `optimize_step4` triggers `get_params`, so reading the global at that
   point sees nothing. Tests:
   `tests/test_optimization_method_routes_from_params.py`.

   It was NOT the last instance of its class -- a config key or flag that is
   documented, accepted, and read by nothing. An audit on 2026-09-14 found ten
   more. The class is now closed, and held closed by a ratchet.

   Closed on 2026-09-14 in the two ways available. Wired, because each already
   had a reader taking its fallback: `mismatch_penalty` (whose consumer
   `occupancy.default_mismatch_penalty` was written for it and received None on
   every call), `max_homopolymer_run`, `gc_clamp_window` and `max_gc_in_clamp`.
   Retired from the schema, because nothing implemented what they named:
   `retries`, `drop_iterations`, `top_set_count`, `selection_metric` and
   `bl_penalty`. The first four appeared only in a module-level `defaults` dict
   in `core/pipeline.py` that itself had no reader; the dict is gone. Setting
   any of the five now produces the unknown-key warning rather than silence.

   `tests/test_no_schema_key_is_inert.py` is the ratchet: every schema key must
   bind a `parameter` global or appear on a short list of keys consumed during
   loading, each with its reason.

   Fixed on 2026-09-14: `filter --gc-tolerance`, `filter --excl-threshold`,
   `expand-primers --optimization-method` and
   `plan-pool --swap-max-evaluations` all carried a real argparse default and so
   beat params.json on every run. They now use the `None` sentinel, and
   `expand-primers` routes through `resolve_optimization_method` rather than
   reading the attribute. `--gc-tolerance` was the costly one: the block it fed
   also computed its own GC window, clamping the lower bound at 0.20 where
   `adaptive_gc_window` releases it to zero below the extreme-AT threshold, so
   on a 19% GC target it excluded exactly the zero-GC primers published AT-rich
   designs are built from. It now routes through `adaptive_gc_window`.

   `tests/test_design_options_have_effect.py`,
   `tests/test_params_json_routes_optional_keys.py` and
   `tests/test_optimizer_config_reaches_optimizers.py` are named as the tests
   that hold the line, but between them they cover about fifty keys and none of
   the ten above -- which is why the class survived being declared closed.
   `tests/test_cli_defaults_do_not_beat_params_json.py` now pins the four flag
   defaults. Extend all four when adding an option.

   **Declared closed twice, and closed neither time.** The audit of 2026-09-16
   found `--design-grid` on `plan-pool` parsed, documented in the help text, and
   never read -- added *after* the second closure. Checking for more found
   `--data-dir` inert on `count-kmers`, `filter`, `prepare-candidates`,
   `optimize`, `design`
   and `evaluate-set`, and `--min-gini-sites` inert on `filter`, which the Key
   Parameters section above documented as working. Both were verified by
   resolving a config whose flag value differed from the file value: the file
   won each time.

   The reason none of the four ratchets caught them is structural.
   `test_no_schema_key_is_inert.py` and `test_params_json_routes_optional_keys.py`
   iterate params.json **schema keys**, so a CLI flag is invisible to them.
   `test_cli_defaults_do_not_beat_params_json.py` checks four argparse
   **defaults** and asserts nothing about whether a flag is read.
   `test_design_options_have_effect.py` calls `run_optimization` directly, so it
   covers the **optimize path only**. None of them asks "does this flag do
   anything".

   `tests/test_every_cli_option_has_an_effect.py` now does. It reads the dispatch
   table out of `main()`, walks from each handler through the functions it calls,
   collects every attribute read off the argparse namespace -- including the
   `merge_args_to_parameter` and `@params_command(merge=...)` routes, which are
   reads by another name -- and fails on any declared option it cannot account
   for. Currently-inert options are listed in `KNOWN_INERT` with a reason, and a
   second test fails on an entry that has since been wired, so the list can only
   shrink.

   The same defect exists one layer down, where a capability is built and tested
   and no command can reach it. The audit found six at once, all from the
   condition-aware pool design work: `CandidateProvider`, `design_sweep`,
   `load_design_grid`, `load_grid_file`, `ensure_positions` and, in practice,
   `beam_search`. Unit tests cannot see this, because a test that constructs the
   thing directly and asserts it behaves passes whether or not anything calls it.
   `tests/test_no_capability_is_unreachable.py` walks transitive reach from the
   dispatch table -- not bare references, since `design_sweep` calling
   `provider.expand` must not make `expand` count -- and holds the same kind of
   shrinking allowlist.

   Do not record this class as closed again. Record what the ratchets cover.

   A third variant surfaced on 2026-09-17, and none of the five ratchets looks
   for it: not a CLI flag and not a schema key nobody reads, but a schema key
   read in some places and not in the one that would spend it. `cpus` reached
   `create_pool` and did not reach step 1. `kmer_counter.run_jellyfish` declares
   `cpus: int = 4` and computes
   `max_workers = min(num_k, cpu_count // max(cpus, 1))`, so that default set
   both the threads per jellyfish process AND how many k values ran at once,
   while the configured value set neither. All four of step 1's call sites
   omitted it. On a 64-core machine with `cpus: 16` the run used 4 threads per
   process and 16 concurrent k values, the transpose of what was asked for.
   Fixed by passing `parameter.cpus`; pinned by
   `tests/test_the_configured_cpu_count_reaches_jellyfish.py`, which walks the
   AST of step 1 rather than asserting a thread count.

9. **Evenness is not measurable from one or two sites** -- FIXED 2026-09-10
   (audit finding B4). `filter.get_gini` keeps a primer when
   `gini.notna() & (gini < max_gini)`. The `.notna()` half was written for
   exactly the case where evenness cannot be measured, but a single-site primer
   produced 0.0, the best value available, so the guard never fired and the
   gate ranked that primer first. 86% of the shipped Prevotella-against-chr21
   pool and 96.2% of the plasmid example sat at 0.0, all with two or fewer
   foreground sites; the three whole-genome GC-tier pools have none at 0.0,
   which is why this was invisible in the runs this project usually inspects.
   `min_gini_sites` is the threshold, default 3, settable in params.json and as
   `--min-gini-sites` on `filter`; `primer_attributes.DEFAULT_MIN_GINI_SITES`
   holds the default. It is threaded into `get_gini_from_txt_for_one_k` as an
   argument rather than read from a module global, because that function runs in
   a spawned multiprocessing worker which would otherwise see the default
   instead of the configured value. `pipeline.check_gini_stage_kept_something`
   refuses to write an empty pool and names the threshold in force.
   Tests: `tests/test_gini_needs_enough_sites.py`,
   `tests/test_min_gini_sites_is_configurable.py`.

10. **The candidate loader filtered on a Tm known to be wrong** -- FIXED
    2026-09-10 (audit finding B5). `melting_temp.py:153` keeps the original melt
    package's GC-fraction bug for compatibility with the random forest retired
    on 2026-09-05. It reads 10.12 C high at k=12 (measured over 20,000 random
    12-mers, sd 0.30 C), so with the old symmetric 15 C margin the loader's
    window was `[min_tm - 25, max_tm + 5]` in true-Tm terms. Nothing was lost on
    plain phi29 and 9.6% of k=12 candidates were lost under DMSO 10% plus
    betaine 1.5 M, silently. `kmer_counter.get_primer_list_from_kmers` now calls
    the same `ReactionConditions.calculate_effective_tm` the gate calls, on the
    window `filter._resolve_tm_window` resolves, with 2 C of stated headroom.
    The shim itself is still used by `rf_preprocessing` and is correct to leave
    there: it is what the bundled model was fitted against.

11. **Near-duplicate primers reached the delivered panel** -- measured, and the
    remedy is OFF by default (audit finding A6). Delivered E. coli set 0 held 9
    pairs at Hamming distance 1 or less and 5 primers sharing the 3' hexamer
    GCGAAA. No optimizer had a similarity rejection test: the greedy picked the
    largest ABSOLUTE new coverage, so a primer with fifty sites and forty-eight
    already covered beat one with five sites all new.
    `dominating_set_optimizer.DEFAULT_REDUNDANCY_THRESHOLD` can skip a candidate
    whose covered bins are already covered above the threshold, on site sets
    rather than on sequences.

    It defaults to **1.0, which disables it**, because measurement did not
    support switching it on. Across three real pools at ten combinations of tier
    and panel size, a 0.9 threshold fired 328,846 times and changed nothing:
    coverage identical in five of six cases and 0.03 points lower in the sixth,
    with the Hamming-1 pair count and the duplicate 3' hexamer count -- the two
    things it was built to reduce -- identical in all six. The reason is
    structural: a candidate more than 90% already covered has a small marginal
    gain by construction, so the greedy's argmax was never going to pick it. The
    criterion and the objective are nearly the same signal. No threshold beats
    disabled on average, and gains sit beside large losses in the same tier:
    M. tuberculosis at n=36 gains 5.53 coverage points at 0.05 and loses 10.72
    at 0.00. The mechanism and its tests are kept; pass an explicit threshold to
    use it. It affects `hybrid` and `background-aware` too, which both call
    `optimize_greedy` for their Stage-1 set cover.

    Taken together, entries 9 to 11 moved the three shipped whole-genome designs
    by almost nothing. Coverage changed by at most 0.17 percentage points, the
    delivered panels have Jaccard 0.993, 0.976 and 1.000 against their
    baselines, and the candidate pool size is identical on all three. The
    evenness rule bites on small targets, which is where the defect was
    measurable in the first place.

12. **CLI startup used to import scikit-learn** -- FIXED 2026-09-10 (audit
    finding E6). `cli_unified.py` imported `core.pipeline` at module scope only
    to make `StepPrerequisiteError` catchable, and `core/pipeline.py` imports
    `rf_preprocessing`, which imported sklearn at module scope. Every invocation
    paid it, `--help` and `show-presets` included, for a model retired from the
    default path on 2026-09-05.

    `StepPrerequisiteError` and `StepValidationResult` now live in
    `core/exceptions.py` (re-exported from `core/pipeline.py`, the same objects,
    so `except` clauses and the two importing tests are unaffected), and the
    sklearn alias fix runs inside `load_model_safely` instead of at import.

    `neoswga --help` now costs **about a third of what it did**. Quote the ratio
    rather than an absolute pair: measured twice hours apart the saving was
    about 3.3x both times, while the before figure itself moved from 1.21 s to
    1.40 s between sessions with no code change, purely with machine load.
    `python -X importtime -c "import neoswga.cli_unified"` reports no sklearn
    entry at all. `tests/test_cli_import_is_light.py` fails if either import
    comes back. What remains is not sklearn: it is pandas, reached through
    `cli/_common.py` -> `reaction_conditions` -> `thermodynamics` -> `utility`.
    That chain is pre-existing and is the obvious next target if CLI startup is
    worth more work.

13. **Five commands measured the host genome they were told to ignore** --
    FIXED 2026-09-10. Each read `bg_prefixes` from params, built a
    `PositionCache` over `fg_prefixes` alone, then handed the background
    prefixes to something that queries that cache by prefix.
    `PositionCache.get_positions` answered an unindexed prefix with an empty
    array, silently, so every background lookup read zero.

    The manifestation is `NetworkOptimizer._evaluate_primer_addition`, whose
    score is `fg_improvement / (1.0 + bg_added)`. With an fg-only cache
    `bg_added` is always 0.0, so every candidate scored as perfectly selective.
    This is the same symptom as Known Issues 5 and 6, reached by a third route:
    not a scan that found nothing and not an integer that saturated, but a cache
    asked for something it does not hold.

    It mattered most in `expand-primers`, which exists to add primers to an
    existing panel, so specificity is the property the user is asking it to
    preserve.

    `get_positions` now raises `MissingPositionsError` for a prefix the cache
    was not built over. That uses a separate `on_unindexed_prefix` knob, not the
    existing `on_missing`: a primer with no hits on an INDEXED prefix is a
    plausible measurement of zero and warns, while a prefix nobody indexed is a
    caller error and raises.

    `tests/test_expansion_counts_background.py` walks the AST for any function
    that forwards `bg_prefixes` while building a cache without them. That check
    found the fifth site after a manual review had settled on four.

14. **Reading the background is not acting on it.** A background-aware stage
    that does not choose the panel changes nothing useful. `expand-primers` was
    fixed on 2026-09-10 to build its `PositionCache` over
    `fg_prefixes + bg_prefixes`, and a real run afterwards still queried the
    host prefix zero times. Three further seams had to be closed before the data
    was read at all: `background_pruning` defaulted to False on the expansion
    path, `PrimerExpander.expand` silently substituted `hybrid` for every method
    it did not recognise including `background-aware`, and `_prune_background`
    would have removed primers from the very panel the user asked to extend.

    Even then it read the host without acting on it. Stage 1.5 pruning is not
    the stage that picks the panel; Stage 2 `_network_refine` is, and it ranked
    on amplification connectivity and unique coverage bins alone. Enabling
    pruning therefore only shrank the pool Stage 2 drew from, and on a
    40-candidate expansion over a 300 kb synthetic pair it moved delivered host
    binding the wrong way, 32 sites to 45. `_STAGE2_BACKGROUND_WEIGHT` adds the
    host as a third normalised axis in that stage, gated on `background_pruning`
    so `hybrid` panels are unchanged (verified identical on all three GC tiers
    at n=12/24/36).

    The general lesson: check which stage produces the delivered result before
    concluding that a measurement reaching the code means it reached the user. A
    query count answers "was it read", not "did it matter".

    **The two Stage 2s carry different things and neither carries both** --
    found 2026-09-19 while wiring Phase 6's deficit objective. The host term
    above lives in `_network_refine`. The objective a search can be steered by,
    `pool_objective`, is read only by `_swap_refine`. So a host-aware expansion
    cannot rank by recovered deficit, and a deficit-targeted one is not
    host-aware. `PrimerExpander._expand_hybrid` chooses between them on
    `background_pruning` and WARNS when target gaps are present but cannot
    steer selection, rather than narrowing the pool to the gaps and then
    ranking by something else. Switching expansion to `swap` wholesale was the
    first attempt and `tests/test_expansion_uses_the_background.py` caught it
    immediately: background-aware and hybrid returned the same panel, because
    the host term had been left behind. Combining them means putting the host
    term into the swap score as a weighted axis rather than its current
    lexicographic tie-break, which is unmeasured.

    `examples/plasmid_example` cannot demonstrate any of this. Six primers
    already cover its 5.4 kb target completely at 3 kb reach, so expansion adds
    nothing and `optimize` early-returns before Stage 1.5. Pin this behaviour on
    a target large enough that Stage 1 over-selects;
    `tests/test_expansion_uses_the_background.py` builds one at 300 kb with no
    external tool.

15. **The guard against the silent zero was itself silent** -- FIXED 2026-09-17
    (Phase 4 increment 3 of the 2026-09-16 pipeline audit).
    `CandidateProvider.ensure_positions` exists to refuse a candidate whose
    binding data is absent, so an unmeasured primer cannot be scored as though
    it bound nothing. It had two defects and each one alone made it useless.

    It returned quietly when no position cache was attached, and nothing in
    production attached one. So the single configuration it was written to
    catch was the configuration in which it did not run.

    And its predicate asked whether the candidate had a hit on ANY prefix:

    ```
    not any(len(cache.get_positions(prefix, sequence, "both"))
            for prefix in cache.fname_prefixes)
    ```

    A candidate with fifty foreground sites and no background entry at all
    therefore passed, which is unknown specificity reported as perfect
    specificity. A candidate indexed against a host it binds nowhere failed,
    though that zero is a measurement and a good one. The question is whether
    there is an ENTRY on EVERY prefix the design scores against, which is the
    distinction `_resolve_missing` already drew for the constructor's primer
    list and `PositionCache.has_entry` now exposes.

    `PositionCache.require_entries` holds the rule once, for both the inventory
    provider and a `--candidates` list, and `pool_planner._prepare_candidate_pool`
    is where `plan-pool` attaches the cache and runs the check, before any panel
    is evaluated. Tests:
    `tests/test_positions_arrive_on_demand.py` for the behaviour and
    `tests/test_the_frontier_is_vouched_for_before_it_is_scored.py` for the
    wiring, the second because the first would have passed throughout the years
    the check was inert.

    `load` and `release` arrive with it. The cache took a fixed primer list at
    construction, which was sufficient only while a design never looked past
    the `max_primer` shortlist. `load` also drops the memoized `both` key for
    the primers it admits: `get_positions` writes one for any primer it is
    asked about, including the empty one it returns for a primer the cache does
    not hold, so without that invalidation a candidate would keep answering
    with the zero it gave before its positions arrived.

16. **The stage that picks the panel is the least informed one** -- WIRED and
    measured 2026-09-19, and deliberately OFF by default (audit
    [pool_selection_audit_2026-09-18.md](docs/validation/pool_selection_audit_2026-09-18.md)).
    `optimize_greedy` takes an `objective`, and supplying it makes Stage 1
    select on occupancy-weighted coverage with a background tie-break instead
    of on unweighted coverage bins. Commit `59a4ee3` added it under the heading
    "The greedy now chooses on the quantity the design is judged on". **No
    production caller passes it.** `hybrid_optimizer.py:711`,
    `dominating_set_adapter.py:164` and `primer_expansion.py:652` all omit it;
    the only caller that supplies it is
    `tests/test_partial_panel_pruning.py:214`.

    It matters exactly where additives matter. Occupancy depends only on the
    primer, so an unweighted bin count misranks two candidates by the ratio of
    their occupancies, and across the pool the Tm gate admits that ratio is 1.8
    on phi29 at 30 C, 7.8 on equiphi29 at 42 C and 8.3 under DMSO 5% plus
    betaine 1 M. On phi29 occupancy is saturated and the unweighted count is
    nearly right, which is the same reason phi29 offers no discrimination.

    The ratchets cannot see this class.
    `tests/test_no_capability_is_unreachable.py` walks reach to FUNCTIONS and
    `optimize_greedy` is reachable; a PARAMETER no caller supplies is invisible
    to it. This is a fifth route into Known Issue 8's class, and the list there
    should be read as covering options and capabilities but not arguments.

    **Now reachable, as `stage1_objective_width` in params.json, defaulting to
    None.** An integer turns the objective on and bounds its cost: the cheap
    bin gain ranks every candidate and only that many leaders are scored. The
    bound is not optional -- one `compute_metrics` call costs 36 ms on the
    Wolbachia design, so a full scan is 14.4 minutes against 39 s for the whole
    run, which is the likeliest reason this stayed unwired.

    **It is off by default because measurement does not support switching it
    on.** At n=6/12/24 it improves the metric it now selects on (effective
    coverage +0.0073, +0.0139, +0.0583) and costs specificity every time
    (density -1.64, -6.50, -7.60; host sites 261 to 456 at n=24) for 3.5x to
    8.4x the runtime. With no width set the delivered panel is identical to
    before, verified at n=12 to the last digit. Same resolution as Known Issue
    11, for the same reason
    ([measurement](docs/validation/stage_one_objective_2026-09-19.md)).

    Two things the wiring clarified. `PoolObjective.coverage()` is coverage and
    nothing else, and `_objective_gain` is a coverage delta, so constraints
    reach Stage 1 only as a STOP rule and never as a selection criterion --
    "select on the quantity the design is judged on" changes what COVERAGE
    means, not whether specificity is weighed. And occupancy weighting favours
    primers whose Tm sits near the reaction temperature, which is the same
    property that makes them bind the host, so the unweighted bin count was
    accidentally the more specific rule. That is Known Issue 17's axis.

    Stage 1 uses `stage1_pool_objective`, NOT `pool_objective`. Writing it into
    the latter overwrote the objective `plan_pool` attaches for Stage 2 -- with
    None on every default run -- silently undoing the Phase 6 fix recorded in
    `attach_search_config`. Four tests caught it; keep the two names apart.

17. **The pool cannot discriminate, and the candidate filter is not the fix**
    -- measured 2026-09-19. Half of this entry's original diagnosis does not
    survive measurement, and the remedy it implied makes panels worse.

    **The floor is not padding the pool.** phi29's default floor does sit 10 C
    below its reaction temperature, but at k = 12 there is nothing down there:
    of 40,000 random 12-mers, 8 fall in the Tm 20-25 band (0.02%). Moving the
    floor changes essentially nothing.

    **Saturation is real and severe.** 86% of random 12-mers sit at or above
    0.998 occupancy at phi29 30 C, where a 4 C mismatch penalty leaves
    discrimination -- matched over single-mismatch occupancy -- at 1.01 or
    less. On the real Wolbachia shortlist, 65% of 2,000 candidates are above
    0.99 occupancy and mean discrimination is 1.11. Specificity in such a pool
    is a property of where sites fall, not of binding.

    **But gating on occupancy delivers a worse panel.** Measured at n=12 with
    the candidate list authoritative: capping occupancy at 0.95 raises the
    delivered panel's discrimination 1.098 to 1.400 and costs coverage 0.7334
    to 0.4397, selectivity density 25.62 to 6.60, with host sites RISING 149 to
    237. Discrimination lives in a tail too small to build a panel from --
    candidates above 2 are 1.6% of the space at k = 12 and 30 C.

    **The lever that works is the reaction.** Same 20,000 12-mers: phi29 30 C
    gives mean discrimination 1.065 with 1.6% above 2; DMSO 10% plus betaine
    1.5 M gives 1.374 and 11.3%; equiphi29 at 42 C gives 2.190 and 36.0%.
    Nothing about the candidates changes in any row. Saturation is a
    phi29-at-30-C problem, not a filtering problem.

    So what ships is a measurement, not a gate. `occupancy.discrimination_profile`
    computes the regime and `log_discrimination_profile` reports it at the end
    of `filter`, warning below `DISCRIMINATION_FLOOR` (1.5, between the two
    measured regimes) and naming the lever that works and the one that does
    not. No candidate is filtered and no delivered panel moves
    ([measurement](docs/validation/occupancy_and_discrimination_2026-09-19.md)).

    There is still no occupancy gate anywhere, and that is now a decision
    rather than an omission; `occupancy_ranking` remains the nearest thing and
    ranks on background load rather than on whether the candidate binds the
    target. Untested: whether a discrimination TERM in selection, as opposed to
    a gate on the pool, would help.

    A methodological note worth keeping. The first run of the gate experiment
    appeared to show density IMPROVING to 44.28. It did not:
    `open_source_or_list` prefers the inventory over the supplied list, so
    swapping `step3_df.csv` only set the frontier SIZE and the run searched the
    inventory as usual. The apparent gain was the smaller frontier. It was
    caught by checking that the delivered primers were actually in the capped
    pool -- none of them were.

    The consequence for the additive lever is measured in the audit. Occupancy
    and mismatch discrimination move in opposite directions along the Tm axis,
    so an additive improves every GC class at or above 6 of 12 and degrades
    every class below it, moving the best class up one step. The two routes to
    specificity conflict: compositional rarity favours GC-rich against an AT-rich
    host, thermodynamic discrimination favours AT-rich, and their correlation at
    k = 12 is about -0.89. An additive is the only lever that moves a candidate
    along the thermodynamic axis without changing its composition, which is why
    the best design measured in `docs/validation/additive_specificity.md` is an
    additive design at k = 12 rather than a longer-primer one.

    Do not conclude from the pool-size table that longer primers help. Above
    k = 15 an additive admits more candidates rather than fewer, and that regime
    is saturated: mean discrimination is 1.09 at k = 18 against 2.99 at k = 12,
    and occupancy spread across the admitted pool collapses to 1.0. A draft of
    the audit recommended k >= 15 before the discrimination column was measured.

18. **`max_gap` and `bg_coverage` are computed and read by nothing that
    selects** -- found 2026-09-18, partly acted on. Both reach
    `step4_improved_df_summary.json` and the reports. Neither appears in
    `normalized_score` or in any optimizer's scoring, and deliberately still
    does not.

    What changed the same day: both are now **constrainable** via
    `max_worst_hole` and `max_host_coverage` (see **Panel limits** above) and
    both are **reported** by the "What limits this panel" table, which also
    names them as having no reference. So a user can hold a panel to either,
    and neither has acquired a default -- no threshold derived from the reach
    separates the published wet-lab winners, so picking one would be the
    scoring change that evidence refuses.

    Also fixed on 2026-09-18: three of the five strand quantities
    `PositionCache.compute_strand_alternation_stats` returns were computed and
    discarded at the call site, and the loop stopped after the first foreground
    prefix so the host was never measured at all.
    `core/strand_metrics.py` collects all five for every foreground genome AND
    the background onto `PrimerSetMetrics.strand_stats`, keyed by prefix.
    `strand_alternation_gap_max` is the one that mattered: exponential
    amplification needs two sites in convergent orientation within the
    polymerase's reach, so the widest gap between opposite-strand sites is the
    closest quantity here to the mechanism, and on the host it is what swga 2.0
    approximates with `within_mean_gap_ratio` and fits against measured
    sequencing breadth. It reaches the report as `convergent_gap` and
    `host_convergent_gap`.

    A prefix the cache cannot answer for is now ABSENT from `strand_stats`
    rather than zero, and the two headline scalars are `None` rather than 0.0.
    They were initialised to 0.0 and left there, so a zero meant either
    "measured zero" or "never asked".

    **The source was fixed on 2026-09-19 and the consumers audited.** A
    one-site panel does NOT genuinely score 0.0 for alternation, which an
    earlier note here got wrong: alternation is the fraction of ADJACENT site
    pairs on opposite strands, so below two sites there is no pair and the
    fraction is 0/0. `compute_strand_alternation_stats` now returns None there,
    while two same-strand sites still return a measured 0.0.
    `strand_coverage_ratio` is min/max over the two strand counts and needs
    only one site, so it is None only with no sites at all; a lone forward site
    really is maximally unbalanced. The gap figures keep the genome length,
    which encodes "no convergent pair anywhere" rather than a missing
    measurement, and `worst_convergent_gap` reads it that way.

    Every consumer was checked and none needed changing: `panel_regime._as_float`
    and `report/metrics._safe_float` preserve None, `headline_strand_scalars`
    already returned `(None, None)`, and the technical report skips a row whose
    value is None. Pinned by
    `tests/test_strand_scores_say_when_they_are_unmeasurable.py`.

    They are still NOT constrainable, but the reason is now different: no
    threshold has a reference, which is why `max_worst_hole` and
    `max_host_coverage` ship unset rather than defaulted. Adding a strand limit
    is a decision about evidence, not a blocked repair.

    `bg_coverage` is the only computed quantity that sees background site
    POSITION. `selectivity_density` and `total_bg_sites` are additive in
    per-primer counts -- `occupancy.weighted_site_load` sums `count * theta` per
    mismatch class and no position enters -- so two backgrounds with identical
    per-primer counts score identically whether their sites are clustered or
    dispersed. That distinction is most of off-target amplification, since SWGA
    needs two convergent sites within the polymerase's reach.

    swga 1.0 made both criteria hard in 2017, and they are the only two things
    that constrain its selection: the clique search runs `--unweighted --all`
    and stores every clique passing the `max_fg_bind_dist` gap cut, so its score
    expression ranks the output and never steers the search. The background side
    is a pruning budget inside the recursion -- each vertex weight is the
    primer's raw background site count (`weight = primer.bg_freq` in
    `graph.py`), and a partial clique is pruned once the summed weight exceeds
    `bg_length / min_bg_bind_dist`. It is NOT a ranking by mean background
    binding distance; that quantity is what the search emits, as
    `bg_len / graph_subgraph_weight`, and an earlier draft of this entry
    conflated the two.

    **The field does not agree that gap statistics belong in the objective.**
    swga 2.0 fits both as `on_gap_gini` and `off_gap_gini`, where `off_gap_gini`
    carries the second largest recorded weight in the only set-level model
    fitted against measured sequencing breadth -- a value nobody has been able
    to verify from a source that opens, since it sits in a CAPTCHA-gated table.
    COATswga (2025) computes no Gini and no gap statistic at all, on the stated
    ground that a per-primer Gini cannot speak for a whole set, which is an
    argument against `max_gini` as much as for the interval-union objective this
    project already uses. So treat background evenness as a measurement to make
    before it is a term to add. No background amplification network is built
    here, in contrast to the foreground network the hybrid and network methods
    build at about 70 kb.

    Worth knowing about the ancestry: swga 2.0 is this project's direct
    ancestor, and the Known Issue 8 class is partly inherited. In its shipped
    master `filter.filter_extra` implements the GC, homopolymer, GC-clamp and
    self-dimer rules the paper describes and nothing calls it, and it would
    raise if called, reading a `default_max_self_dimer_bp` that `parameter.py`
    never assigns. Its step 2 also computes `ratio = bg_count / fg_count`, where
    lower is more specific, then keeps `sort_values(by=["ratio"],
    ascending=False)[:max_primer]`, retaining the LEAST specific survivors.
    NeoSWGA sorts that ascending. All three were reported from that repository's
    source and were not re-verified here.
