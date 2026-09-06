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
  `minimal_primer_selector` as a post-process.
- **`position_cache.py`**: in-memory binding-position cache, about 1000x faster
  than re-reading the HDF5 files.
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
count-kmers            filter                 prepare (`score`)      optimize
     |                    |                     |                       |
     v                    v                     v                       v
 *_Xmer_all.txt  -->  step2_df.csv +    -->  step3_df.csv      -->  step4_improved_df.csv
 (k-mer counts)       positions.h5           (ordered candidates)   (final primer sets)
```

**File outputs** (in `data_dir`):
- `step2_df.csv`: Filtered primers with fg_freq, bg_freq, gini, Tm
- `filter_stats.json`: Real per-stage filtering funnel counts (rendered in reports)
- `step3_df.csv`: The candidate pool the optimizer reads, carrying the step-2
  measurements in a deterministic order (gini, then the primer sequence). It no
  longer holds an amplification score -- see **The `score` stage** below.
- `step4_improved_df.csv`: Final optimized primer sets with enrichment scores
- `step4_improved_df_summary.json`: Authoritative optimizer metrics the report reads (coverage, effective_fg_coverage, selectivity_ratio, selectivity_density, fg_total_length/bg_total_length,
  effective_fg_sites/effective_bg_sites, selectivity_mode, ensemble_comparison, per_target_coverage, strand metrics).
  Also `unindexed_candidates`: how many candidates the foreground position
  index could not place. Those cover nothing and so are invisible to
  selection; the pipeline path refuses rather than reporting a coverage
  figure that describes only the rest of the pool.
- `*_positions.h5`: HDF5 files with primer binding positions
- `*_{k}mer_all.provenance.json`: A sidecar recording the genome each k-mer
  table was counted from (absolute path, content fingerprint, k). `count-kmers`
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
neoswga score -j params.json        # Step 3: Prepare the candidate pool
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
- `score --enable-qa` re-orders step3_df.csv by a `composite_score`. With the
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

### The `score` stage

`score` prepares the candidate pool; it does not score it. The bundled random
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
neoswga optimize -j params.json --optimization-method=background-aware # clinical, 10-20x bg reduction
neoswga optimize -j params.json --optimization-method=network          # Tm-weighted, dimer-aware
neoswga optimize -j params.json --optimization-method=clique           # guaranteed dimer-free set
neoswga optimize -j params.json --optimization-method=ensemble         # run all, keep best
neoswga optimize -j params.json --optimization-method=ensemble --ensemble-combine=union  # re-optimize over pooled primers
```

**Optimization Method Comparison**:

| Method | Speed | Best For | Notes |
|--------|-------|----------|-------|
| `hybrid` | Medium | General use (default) | Combines network + set-cover approaches |
| `dominating-set` | Fast | Large primer pools | Graph-based set cover, ln(n) approximation |
| `background-aware` | Slow | Clinical applications | 10-20x background reduction, three-stage |
| `clique` | Slow | Sets that must be dimer-free | Max-clique on the compatibility graph (swga 1.0's approach). The only method that GUARANTEES no dimerising pair; the others penalise dimers but can accept one. Pools of ~200 candidates; not in the default ensemble |
| `network` | Medium | Tm-weighted selection | Dimer penalty aware |
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

### Utility Commands
```bash
neoswga validate --quick            # Validate installation
neoswga build-filter genome.fna ./  # Build Bloom filter for large background
neoswga show-presets                # Show reaction condition presets
```

## Key Parameters (params.json)

**Primer filtering**:
- `min_k`, `max_k`: Primer length range (default: 6-12, use 12-18 for longer primers)
- `min_fg_freq`: Minimum foreground frequency (default: 1e-5)
- `max_bg_freq`: Maximum background frequency (default: 5e-6)
- `max_gini`: Maximum Gini index for binding evenness (default: 0.6)
- `max_primer`: Primers to keep after filtering (default: 500)

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
- `optimization_method`: **currently inert in params.json — see Known Issue #8.**
  Use `--optimization-method` on the CLI. Values: 'hybrid' (default),
  'dominating-set' (fast), 'background-aware' (clinical), 'network'
- `num_primers`, `target_set_size`: Desired primer set size (default: 6)
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

**Application profiles** (`--application`, and the weighting used to pick an
ensemble winner):

| Application | Coverage Target | Specificity | Typical Size | Use Case |
|-------------|-----------------|-------------|--------------|----------|
| `discovery` | 90% | 60% | 10-15 | Pathogen discovery, maximize sensitivity |
| `clinical` | 70% | 90% | 6-10 | Diagnostics, minimize false positives |
| `enrichment` | 80% | 75% | 8-12 | Sequencing enrichment, balanced |
| `metagenomics` | 95% | 50% | 15-20 | Capture diversity |

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

   This was the last known instance of one class: a config key or flag that is
   documented, accepted, and read by nothing. `additionalProperties: true`
   means none of them warn. `tests/test_design_options_have_effect.py`,
   `tests/test_params_json_routes_optional_keys.py` and
   `tests/test_optimizer_config_reaches_optimizers.py` are the tests that hold
   the line; extend them when adding an option.
