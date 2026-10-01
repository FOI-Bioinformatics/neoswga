# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

NeoSWGA is a command-line tool for selecting primer sets for selective whole-genome amplification (SWGA). See [README.md](README.md) for user-facing documentation and quick start.

**External dependency**: a k-mer counter in PATH. KMC3 is used when installed and jellyfish otherwise; setting `"kmer_counter"` in params.json requires the named one (Known Issue 24).

Python 3.13 only (`requires-python = ">=3.13"`).

## Where the detail lives

This file is the working summary. The dated record behind each statement, with
the measurements, was moved out on 2026-10-01 and is read on demand:

| Document | Holds |
|---|---|
| [docs/development/implementation-notes.md](docs/development/implementation-notes.md) | Per module, output file, command and params.json key: what it does and why |
| [docs/development/design-contracts.md](docs/development/design-contracts.md) | The rules a design run enforces, each with the failure it was written for |
| [docs/development/KNOWN_ISSUES.md](docs/development/KNOWN_ISSUES.md) | Known Issues 1-24 in full. The numbering is stable; code and tests cite it |
| [docs/params-reference.md](docs/params-reference.md) | Every params.json key, generated from the schema by `scripts/render_schema.py` |
| `docs/validation/` | One measurement record per claim or default |
| `.claude/skills/neoswga-cli/SKILL.md` | Commands outside the four-step pipeline; loads on demand |

Read the relevant section there before changing the behaviour it describes.
When adding a record, add it there and keep this file to a line or two.

## Architecture

- Entry point: `neoswga/cli_unified.py` (`neoswga = neoswga.cli_unified:main`
  in `pyproject.toml`). Handlers live in `neoswga/cli/`.
- `neoswga/core/` holds about 140 modules; `ls` and the module docstrings are
  the current list. What the filenames do not say:

| Module | Note |
|---|---|
| `unified_optimizer.py` | Dispatches the optimizers: `hybrid_optimizer`, `dominating_set_adapter` + `dominating_set_optimizer`, `network_optimizer`, `background_aware_optimizer`, with `minimal_primer_selector` as a post-process. `background-aware` is `BackgroundAwareBaseOptimizer`, which delegates to `HybridOptimizer` |
| `exceptions.py` | `StepPrerequisiteError`, `StepValidationResult` and the `DesignError` family. Keep it free of dependencies beyond `typing` and `dataclasses`; `core/pipeline.py` re-exports the same objects |
| `occupancy_coverage.py`, `selectivity.py` | Extracted from `base_optimizer` when it reached its size budget; re-exported from it as the same objects |
| `lazy_dimer.py` | `dimer_screen(pool, max_dimer_bp)` is the one decision about dense matrix against pairwise screen (threshold 4,000 candidates) |
| `candidate_source.py` | `open_design_source` is the rule every design command asks. An explicit list wins; an absent inventory is None; a mismatched inventory raises. Only `plan-pool` advances the frontier, deliberately |
| `position_index.py` | The only reader and writer of `*_positions.h5` (sorted-blocks layout) |
| `position_cache.py` | In-memory positions. `load`/`release` move the window; `has_entry` and `require_entries` ask whether an ANSWER exists, which differs from a non-zero answer |
| `kmer_backend.py`, `kmer_tables.py`, `query_scan.py` | Counter invocation; the one place anything asks about a k-mer table; a direct scan for a host with no table |
| `design_request.py`, `design_result.py`, `panel_evaluation.py` | The resolved frozen request and its hash; run state against termination reason against qualification; one `PanelAssessment` |
| `pool_planner.py`, `pool_objective.py`, `panel_acceptance.py` | `plan-pool`, the shared objective and `shortfall`, the configured panel limits. `repair_panel` is the one bounded repair |
| `bam_coverage.py`, `depth_policy.py` | `open_alignment` is the only door to a BAM or CRAM |
| `delivered_set.py` | Which `set_index` a command reads (default 0) |
| `concurrent_runs.py` | Translates the HDF5 lock error from two runs sharing a directory |
| `registry/model_evidence.json` | Per chemistry constant, its evidence status and supported domain; checked in `resolve_design_request` |
| `gpu_acceleration.py` | Not reached by any pipeline stage; `--use-gpu` says so |
| `advanced_features.py`, `gc_adaptive_strategy.py` | Read as optional but are wired into kept paths |
| `report/` | Quality reports. Reads the saved panel assessment last, and `effective_conditions` from the run manifest in preference to params.json |
| `models/random_forest_filter.skops` | Retired from the default path; reachable through `--amp-model` |

### Data flow

```
count-kmers            filter                 prepare-candidates     optimize
     |                    |                     |                       |
     v                    v                     v                       v
 k-mer tables    -->  step2_df.csv +    -->  step3_df.csv      -->  step4_improved_df.csv
 (KMC db or text)     positions.h5           (ordered candidates)   (final primer sets)
```

Files in `data_dir`:

- `step2_df.csv`, `filter_stats.json`: the shortlist and the per-stage funnel.
- `step3_df.csv`: the candidate pool in a deterministic order. It carries no
  amplification score since 2026-09-05.
- `step4_improved_df.csv`: up to `max_sets` ALTERNATIVE sets, numbered by
  `set_index`. Set 0 is the one the summary describes (Known Issue 19).
- `step4_improved_df_summary.json`: the authoritative optimizer metrics the
  report reads.
- `*_positions.h5`: one per prefix and k. Keys are never canonicalised. An
  empty entry (scanned, binds nowhere) and an absent one (never scanned) are
  distinct.
- `*_{k}mer_all.provenance.json`: which genome a table was counted from, by
  full SHA-256. A record under the older partial hash is unknown, not stale.
- `run_manifest.json`: one entry per step. `effective_conditions` is the
  reaction the step ran under, which can differ from params.json.
- `design_failure.json`: written when a design fails; `export` refuses while
  it is present.

## CLI Commands

```bash
neoswga count-kmers -j params.json         # Step 1: k-mer counts
neoswga filter -j params.json              # Step 2: filter candidates, write position indexes
neoswga prepare-candidates -j params.json  # Step 3: write the ordered candidate pool
neoswga optimize -j params.json            # Step 4: select primer sets
neoswga design -j params.json              # all four
neoswga plan-pool -j params.json [--design-grid grid.json]
neoswga validate --quick                   # installation check
neoswga validate --smoke -j params.json    # config check against a packaged 6 kb target, ~4 s
neoswga build-filter --genome genome.fna -o ./
neoswga show-presets
```

- `prepare-candidates` was named `score` until 2026-09-21; there is no alias.
  It prepares the pool and does not score it. `--amp-model` restores the
  retired random-forest score and its `min_amp_pred` gate.
- `optimize` refuses with a step-4 prerequisite error when `step3_df.csv` is
  missing or empty, the position files are absent, or the index covers only
  part of the pool. The remedy is to re-run `neoswga filter`.
- `--enable-qa` is accepted by every step and is per-invocation.
- `export`, `interpret`, `report` and `simulate` read one set, set 0 by
  default; `--set N` on `export` and `interpret` selects another. `export`
  refuses a set other than 0 because nothing assessed it;
  `--allow-unqualified` overrides.

### Optimization methods

`--optimization-method` on the CLI, or `optimization_method` in params.json.
An explicit flag wins; an absent flag does not.

| Method | Speed | Notes |
|---|---|---|
| `hybrid` (default) | Medium | Set cover, then network refinement |
| `dominating-set` | Fast | Graph-based set cover |
| `background-aware` | Slow | `hybrid` with a host-binding term in pruning and in Stage 2 |
| `network` | Medium | Tm-weighted, dimer-screened; stops short instead of relaxing the screen |
| `clique` | Slow | The only method that guarantees no dimerising pair; pools of about 200 |
| `ensemble` | Slow | Runs several on one shared cache and keeps the best by `normalized_score`, then smaller set, then method name. `--ensemble-combine union` re-optimizes over the pooled primers |

Raw `score` is not comparable across methods; `normalized_score` is.

**Coverage reach**: optimizers select and report coverage at the realistic
per-primer reach (about 3 kb for phi29), while network connectivity uses
processivity (about 70 kb). A hybrid run prints two coverage figures, labelled
`(estimated, binned)` and `(measured)`; the measured one is authoritative.

### Choosing the set size

`num_primers` is a request, not a guarantee. The delivered panel is never
larger and may be smaller: Stage 1 stops at the coverage target, and selection
stops instead of admitting a pair above `max_dimer_bp`.

- `--auto-size` inverts a closed-form coverage curve. It does not read the
  candidate pool or the background and is clamped to at most 20 primers.
- `--show-frontier` builds a coverage against fg/bg frontier over the real
  pool, for 4 to 20 primers.
- The marginal coverage table `optimize` prints has no size limit. Each row is
  a lower bound on re-optimizing at that size and says nothing about
  specificity.
- `--minimize-primers` removes one primer at a time, which reaches a local
  optimum. It logs what it removed and the figure the target was compared
  against. `--coverage-metric effective|raw` chooses that figure; the values
  `realistic`/`processivity` belong to a different `coverage_metric` in
  `coverage.polymerase_extension_reach` and are refused here.

## Key Parameters (params.json)

Full list: [docs/params-reference.md](docs/params-reference.md). Rationale and
measurements: implementation-notes.md.

- **Filtering**: `min_k`/`max_k` (6-12; 12-18 for longer primers),
  `min_fg_freq` (1e-5), `max_bg_freq` (5e-6), `max_gini` (0.7),
  `min_gini_sites` (3; the `--min-gini-sites` flag is inert, use the key),
  `max_primer` (500; bounds the shortlist, not what a design can reach),
  `candidate_retention` (`all_qc` default, or `post_gini`; `legacy` is refused).
- **Thermodynamics**: `polymerase` (`phi29` 30 C, `equiphi29` 42-45 C, `bst`
  60-65 C, `klenow` 25-40 C), `reaction_temp`, `na_conc`, `mg_conc`,
  `min_tm`/`max_tm`, and the additives `dmso_percent`, `betaine_m`,
  `trehalose_m`, `ethanol_percent`, `urea_m`, `tmac_m`, `formamide_percent`.
  A design that sets glycerol is refused: nothing computes its Tm effect.
- **Selection**: `num_primers`/`target_set_size` (6), `max_dimer_bp` (3,
  maximum 7), `allow_dimer_relaxation` (false), `max_dimer_dg` (unset; can only
  tighten the screen), `max_sets` (5), `iterations` (8; bounds the search for
  alternatives, not the primary selection), `coverage_metric`,
  `min_per_target_coverage` (checked and reported, not repaired).
- **Search control**: `objective_scan_width` (64), `max_frontier_refills` (4),
  `swap_max_evaluations` (per stage, not a total), `total_search_evaluations`
  (the only total; None by default, so by default no total bound exists),
  `stage1_objective_width` (None; off by default after measurement).
- **Panel limits**, all unset by default and params.json only:
  `min_selectivity_density`, `max_background_sites`, `max_worst_hole`,
  `max_mean_gap`, `max_evenness`, `max_host_coverage`. Set none and the
  delivered panel is unchanged. Set one and `optimize` attempts one bounded
  repair and reports the outcome. A background-measured limit with no
  background genome is refused.
- **`--application`** (`discovery`, `clinical`, `enrichment`, `metagenomics`)
  weights the ensemble winner and steers `--auto-size`. On `hybrid` it does
  not change what is selected.

## Rules a change must keep

Each of these was learned from a defect; the record is in the linked documents.

**Unknown is not zero, and not success.**
- A quantity that could not be measured is `None`, absent, or an error. It is
  never 0.0, an empty array, or a default that happens to be the best value.
  `.get(key, 0.0)` on a gap, a background count or a dimer energy is this
  defect.
- A required calculation that fails raises a `DesignError` subclass
  (`InvalidDesignRequest`, `ReferenceDataError`, `UnsupportedModelError`,
  `ModelEvaluationError`). Do not catch it and substitute. A measured QC
  rejection is a different thing and stays a rejection.
  `SearchBudgetExhausted` is a recorded stopping point, not a failure.
- An artifact that cannot be parsed blocks; it does not pass.
- Genome coordinates are int64 (`position_cache.POSITION_DTYPE`). The scan is
  chunked at `string_search.MAX_SCAN_CHUNK` because pyahocorasick finds nothing
  past 2**31 characters; that chunking is load-bearing.

**One door per resource.**
- Position indexes: `core/position_index.open_index`. A raw `h5py` lookup finds
  nothing in the sorted-blocks layout.
- Alignment files: `bam_coverage.open_alignment`.
- K-mer tables: `kmer_tables` (`table_exists`, `discover_prefixes`,
  `counts_for`), never a filename glob. A table is a KMC database or a text
  file. A count lookup has a fixed cost per call, so batch by k.
- Dimer screening: `lazy_dimer.dimer_screen`. Candidate pool:
  `candidate_source.open_design_source`. Delivered set: `delivered_set`.
  Bounded repair: `pool_planner.repair_panel`.

**Options must reach the code that decides (Known Issue 8).**
- A CLI flag that also has a params.json key takes `None` as its argparse
  default, so an absent flag does not beat the file.
- A new params.json key needs the schema entry, a module-level default in
  `core/parameter.py`, a reader, and a regenerated `docs/params-reference.md`.
- Anything a delegate optimizer must read goes through
  `swap_refinement.attach_search_config`. Stage 1 reads `stage1_pool_objective`;
  Stage 2 reads `pool_objective`; keep them apart.
- Assert the path, not the two ends. Extend the ratchets when adding an option:
  `test_every_cli_option_has_an_effect.py`,
  `test_no_capability_is_unreachable.py`, `test_no_schema_key_is_inert.py`,
  `test_cli_defaults_do_not_beat_params_json.py`,
  `test_params_json_routes_optional_keys.py`,
  `test_design_options_have_effect.py`. Their allowlists may only shrink.
- Do not record this class as closed.

**Ordering and globals.**
- `get_params` runs inside the step, so code that runs before it must resolve
  from the params FILE. `resolve_design_request` does; all three design
  commands call it.
- Do not pair a `parameter` global with an argument the call was given. Under
  `pytest -n 8` that paired one test's prefix with another's FASTA. Reference
  identity travels on the request.
- Step 2's sort ends on the primer sequence, so the counter's emission order
  cannot change a result.

**Structure.**
- CLI startup must not import scikit-learn (`tests/test_cli_import_is_light.py`).
- `base_optimizer.py` and `pipeline.py` are at their size budgets. Extract to a
  new module and re-export the same objects.
- New packages must be listed in `pyproject.toml`'s `packages`; a test compares
  the list against disk.
- ILP models use HiGHS. CBC kills the interpreter on Python 3.13, and
  `core/ilp_solver.py` refuses instead of falling back.
- KMC is invoked with `-ci1 -cs1000000`, and intersections with `-ocleft`.
  All three are load-bearing and silent when wrong.

**Defaults follow measurement.**
- A capability with no demonstrated benefit ships off by default (redundancy
  threshold, Stage 1 objective, beam search on `optimize`, frontier refill on
  `optimize`, an occupancy gate). Do not switch one on without a new
  measurement.
- Panel limits and strand limits have no defaults because no threshold has a
  reference.
- The profile weights `tm_weight` and `uniformity_weight` are inert on
  `hybrid`. Wiring them moves every delivered panel and is a decision, not a
  repair.

## Known Issues (index)

Full text: [docs/development/KNOWN_ISSUES.md](docs/development/KNOWN_ISSUES.md).

1. **Large backgrounds**: `build-filter` exists, but exact counting of hg38 at
   k=12 is minutes and a 138 MB table. Measure before reaching for it.
2. **skops**: the format is version-tolerant, its default trust list is not;
   `rf_preprocessing._TRUSTED_MODEL_TYPES` names what a model may reconstruct.
3. **Memory**: `filter` loads all background k-mers.
4. **PositionCache strand** is `'forward'`, `'reverse'` or `'both'`, not `+`/`-`.
5. **pyahocorasick past 2 Gb** finds nothing, silently. The chunked scan works
   around it; a canary test fails when the installed version changes.
6. **A partial background is not a specific design**: `selectivity_ratio`
   moves with background size. Compare `selectivity_density`.
7. **int64 coordinates**. int32 saturated on hg38 and hid host sites. Test
   against a whole genome, not a chromosome.
   - Unnumbered, between 7 and 8: `--design-grid` on `plan-pool`;
     `--min-fg-bg-ratio` orders candidates and no longer deletes them; the
     objective must reach the stage that refines; Stage 1 is deliberately not
     constraint-aware (reference figure for the Wolbachia pool: density 60.112
     at coverage 0.6535).
8. **Inert options**: declared, documented, read by nothing. Found repeatedly,
   in flags, schema keys, capabilities and arguments.
9. **Evenness needs sites**: Gini is NaN below `min_gini_sites`.
10. **The candidate loader** uses the same Tm as the gate.
11. **Near-duplicate primers**: the redundancy threshold exists and is disabled
    (1.0); measurement did not support enabling it.
12. **CLI startup** no longer imports scikit-learn.
13. **Unindexed prefix**: `PositionCache.get_positions` raises
    `MissingPositionsError` instead of answering with an empty array.
14. **Reading the background is not acting on it**. `_network_refine` carries
    the host term and `_swap_refine` reads the objective; neither carries both.
15. **`require_entries`** asks for an entry on every prefix, not a hit on any.
16. **Stage 1 objective**: wired as `stage1_objective_width`, off by default;
    it gains coverage and costs specificity and runtime.
17. **The phi29 30 C pool cannot discriminate**, and gating candidates on
    occupancy gives a worse panel. The lever is the reaction. `filter` reports
    the regime and does not enforce it.
18. **`max_gap` and `bg_coverage`** are constrainable and enter no score.
    Strand figures are `None` when unmeasurable.
19. **One set, not every alternative**: see `delivered_set.py`.
20. **Mixed oligo lengths** work; on the one pair measured they bought nothing.
21. **HDF5 lock error (errno 35)** is two processes sharing a data directory.
    Do not suggest `HDF5_USE_FILE_LOCKING=FALSE`.
22. **Bloom path**: capacity bounds distinct k-mers; filters record their k
    range. Use `--from-kmers` at host scale. Not yet built against a host
    genome here.
23. **Python 3.13 only**; HiGHS, not CBC; bioconda has no 3.13 build of
    `kmer-jellyfish`; CI installs the optional extras so those tests run.
24. **KMC3 preferred**, tables read as databases. About 7x faster and about
    18x the memory of jellyfish on hg38 at k=12. A host with no table is
    scanned by `query_scan`.

## Reference genomes available locally

All gitignored, so `git ls-files` shows none. Check here before concluding
that something cannot be measured.

| Reference | Size | Records | Location |
|---|---|---|---|
| hg38 | 3.30 Gb | 705 | `tests/validation/genomes/human_full.fna` |
| *Drosophila* | 144 Mb | 1,870 | `examples/wolbachia_pool_design/input/` |
| human chr21 | 46.7 Mb | 1 | `tests/validation/genomes/` |
| E. coli | 4.64 Mb | 1 | `tests/validation/genomes/` |
| M. tuberculosis | 4.41 Mb | 1 | `tests/validation/genomes/` |
| Prevotella | 3.17 Mb | 2 | `tests/validation/genomes/` |
| S. aureus | 2.82 Mb | 1 | `tests/validation/genomes/` |
| Wolbachia wMel | 1.27 Mb | 1 | both locations |
| two plasmids | ~6 kb | 1 | `examples/plasmid_example/`, packaged in `core/smoke/` |

`examples/wolbachia_pool_design/work` is prepared (12-mer indexes for wMel and
*Drosophila*, 2,000-candidate shortlist; one `compute_metrics` call costs about
40 ms) and is the one to reach for. Prevotella and chr21 also carry tables and
indexes. The multi-record references (*Drosophila*, hg38, Prevotella) are the
ones that exercise record-join geometry. There is no BAM or CRAM anywhere, so
no coverage reach has been measured against sequencing depth.

## Testing

```bash
pytest tests/ -n 8                     # full suite, foreground
pytest tests/test_hybrid_optimizer.py  # one file
pytest tests/ -rs                      # show skip reasons; compare the skip count, not only passes
neoswga validate --quick
```

- A full run leaves `git status --porcelain` unchanged; a test enforces it.
- Tests needing the generated plasmid example request the
  `primed_plasmid_example` fixture. A `skipif` on `plasmid_example_ready()` is
  evaluated at collection, before the fixture primes the example. Fixtures
  that write HDF5 directly (`tests/test_hybrid_optimizer_run.py`) need no
  external tool.
- Import-time code in a test module runs at collection for the whole session.
- CI runs serially on Python 3.13 with
  `.[dev,improved,bam,viz,interactive]` and both k-mer counters.
- Integration tests are in `tests/integration/`;
  `test_strict_design_pipeline.py` drives each refusal end to end.

## Code Patterns

```python
from neoswga.core.parameter import get_params
params = get_params('params.json')

from neoswga.core.utility import create_pool
with create_pool(cpus) as pool:
    results = pool.map(process_func, items)

from neoswga.core.position_index import open_index
with open_index('g_12mer_positions.h5') as index:
    positions = index.get(primer_sequence)  # None: never scanned; empty: binds nowhere
```

## Working habits this project has paid for

- Check which stage produces the delivered result before concluding that data
  reaching the code changed the answer. A query count shows it was read.
- Compare a baseline against the project's recorded figure for the same
  quantity before interpreting a difference.
- An achievability figure that omits a constraint bounds nothing. Evaluate a
  candidate panel through the acceptance path and against a control.
- Keep the script behind a measurement. An inline measurement cannot be
  re-run.
- Test against a whole genome and a multi-record reference; a single
  chromosome reaches neither the 2**31 limits nor a record join.
- State the size of a claim. One pair, one panel size or one seed is "no
  benefit demonstrated", and "not determined" is for cases where more evidence
  would settle it.
