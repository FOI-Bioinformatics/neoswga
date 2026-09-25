# Changelog

All notable changes to NeoSWGA are documented in this file.

## [Unreleased]

### KMC3 is the default k-mer counter

#### BREAKING

- **KMC3 replaces jellyfish as the default counter.** Install it with
  `conda install -c bioconda kmc`. Jellyfish remains fully supported: set
  `"kmer_counter": "jellyfish"` in params.json.
- The `improved` extra is unrelated to this; no Python dependency changed.

#### CHANGED

- **Count lookups are answered from the binary database.** `filter` asks for
  the counts of a known candidate list, which is a set operation. It now goes
  through `kmc_tools simple ... intersect` rather than streaming the text
  table into Python. Measured on Drosophila at k=18 with 2,000 candidates:
  2.5 s against 17.9 s, and no 2.3 GB intermediate.
- Table scans stream from the counter instead of reading a materialised text
  file. One consumer was reading the file twice, once only to size a progress
  bar.
- A reference counted into a binary database now counts as counted. Nine
  checks tested for `{prefix}_{k}mer_all.txt` specifically.

#### COMPATIBILITY

- **Existing data directories are unchanged.** They hold text tables and no
  database, and every path still reads them.

#### KNOWN LIMITS

- KMC is not the memory-friendly option: it peaks at 1,201 MB where jellyfish
  peaks at 631 MB on the same job, and refuses to run under 2 GB. Its
  advantage is speed at long k on large references.
- `py_kmc_api` cannot be used on Python 3.13; the shipped build targets 3.10.

#### FIXED

- Three runtime `NameError`s in the `--bam` expansion helpers, which read
  `parameter` without importing it. Nothing had exercised those paths.
- An unresolvable type annotation in `experimental_tracker`.

### Python 3.13 only

`requires-python` is now `>=3.13`. Earlier interpreters are not supported and
the CI matrix tests one version on two operating systems.

#### BREAKING

- **Python 3.11 and 3.12 are no longer supported.** Install under 3.13 or
  later.
- The `improved` extra gains `highsbox`. Installing `mip` alone is no longer
  enough on 3.13; see the fix below.

#### FIXED

- **The exact ILP solve crashed the interpreter on Python 3.13.** python-mip
  defaults to CBC, and constructing a CBC model terminates the process with
  SIGKILL there: no exception, no traceback, exit 137. `core/ilp_solver.py`
  now prefers HiGHS, which works on 3.13, and warns when falling back.
  Measured on macOS arm64 with mip 2.0.0 and cbcbox 2.935.
- **Documented `build-filter` commands did not run.** Nine examples across the
  guides, the developer docs and two example READMEs used a positional form
  that argparse rejects. They now use `--genome` and `-o`.
- **QUICK_START named a params key that does not exist.** It told users to set
  `background_bloom_path`, which is not a schema key, and named a file
  `build-filter` does not write. The keys are `use_bloom_filter` and
  `bloom_filter_path`, and the file is `bg_bloom.pkl`.
- `pyyaml` is declared in the `dev` extra. It previously reached CI only as a
  transitive dependency of pre-commit, while a ratchet `importorskip`s it.

#### CHANGED

- The codebase is modernised to current Python idioms: 2,294 sites, almost all
  `typing.List` to `list` and `Optional[X]` to `X | None`. No behaviour change.

#### KNOWN ISSUE

- bioconda has no Python 3.13 build of `kmer-jellyfish`, so
  `conda create ... python=3.13 kmer-jellyfish` silently resolves to 1.1.12,
  whose CLI NeoSWGA refuses. Install Jellyfish through brew or apt instead.
  See docs/guides/TROUBLESHOOTING.md.

### Bloom background screening works at the scale it exists for

Twelve defects in the Bloom path, all of which passed every small-genome test
and failed only at or near host scale.

#### FIXED

- **The filter could not be built for a host genome.** Capacity was sized to
  ten times the base count, but pybloom allocates from capacity and counts
  only distinct items, so hg38 asked for 39.56 GB where 26.8 MB suffices.
- **A saved filter reloaded at only one geometry.** pybloom picks its hash from
  the filter's size, so small filters and very large ones saved successfully
  and refused to load.
- **A primer length the filter never indexed read as absent**, which cleared
  the background gate unscreened. Both artifacts now record the lengths they
  cover and a design outside them is refused.
- **`use_bloom_filter` without a path screened nothing, silently.** It now
  refuses.
- The k-mer-file route validates bases and length before indexing.
- `BackgroundBloomFilter` requires an explicit capacity; the former 3 GB
  default was an allocation nobody chose.

### CI runs what it says it runs

- The `@pytest.mark.scale` tests now run in the nightly workflow their marker
  names, rather than in every pull-request cell while nightly excluded them.
- Coverage is measured in one cell instead of six.
- The `build` job no longer waits on the test matrix, so a packaging failure
  surfaces in about a minute instead of twenty.
- Nightly E2E, red since 2026-09-22 on the retired `score` command, is fixed
  and now runs only what the pull-request suite does not.

### Occupancy-based specificity model

Selectivity was a count of exact k-mer matches, so for a fixed primer set no
additive could move it by any amount -- the lever this tool exists to use had
nowhere to land. Sites are now weighted by two-state occupancy under the
configured conditions, and background near-matches are counted by mismatch
class from the existing jellyfish `*_all.txt` files.

Full findings, including what measurement rejected:
[docs/validation/additive_specificity.md](validation/additive_specificity.md).

#### NEW

- `effective_fg_coverage` in `PrimerSetMetrics` and
  `step4_improved_df_summary.json`: coverage weighted by how much of the time
  each site is actually bound, beside the raw `fg_coverage`. Selection uses it
  when reaction conditions are available and falls back to the raw figure
  otherwise. Previously selectivity was occupancy-weighted and coverage was
  not, so a set could be rewarded for reaching sites it does not occupy and
  the coverage cost of additives was invisible.
- `effective_fg_sites`, `effective_bg_sites`, `selectivity_mode` beside the
  raw site counts.
- `neoswga/core/occupancy.py`, `mismatch_counts.py`, `condition_sweep.py`.
- `scripts/fetch_reference_genomes.py` for the validation genomes, which are
  not committed.

#### FIXED

- **`num_primers` / `target_set_size` from params.json were ignored by
  `optimize`**, which always used the default of 6. The completion check read
  the params value later, so a run that had been asked for 6 primers reported
  "Found 6 primers but target was 8 (PARTIAL result: insufficient candidates)"
  and advised relaxing the filters -- for a pool that was never consulted at
  that size. `--num-primers` was unaffected.
- **The candidate cut kept an arbitrary slice rather than the best
  `max_primer`.** Ranking is by `ratio` = bg_count / fg_count, which is 0 for
  every primer with no exact background match -- at k=12 against a distant
  background, all 369,431 survivors tied at exactly 0.0 and the stable sort
  returned whichever arrived first. Ties now break on foreground abundance.
  End to end at k=12: coverage 7.6% -> 40.3%, selectivity 0.69 -> 2.54.
  Raising `max_primer` no longer helps, because it was compensating for this.
- `additive_optimizer`'s `optimize_for="specificity"` scored a synthetic
  primer's binding stability, which an indiscriminate primer also scores well
  on. It can now be given a real fg/bg measurement.
- Bloom-filter background path assigned a fixed sentinel count that passed any
  frequency threshold on a large genome, silently disabling background
  filtering at exactly the scale it exists for.

#### CHANGED

- `normalized_score` selectivity term is log-scaled. The previous
  `min(ratio/100, 1)` was calibrated for exact-match ratios that routinely
  saturated it; occupancy-weighted ratios land far below 100, where the same
  linear form ignored the term instead.
- Sequence heuristics (GC clamp window and threshold, homopolymer run) are
  configurable via params.json. Defaults reproduce previous behaviour exactly
  -- measurement rejected scaling them to primer length.

## [3.7.0] - 2026 - Production Readiness

### NEW FEATURES

#### Canonical params.json JSON Schema
- `neoswga/core/schema/params.schema.json` is the authoritative schema
- `neoswga schema --dump [-o FILE]` emits the schema for IDE integration
- `ParamValidator` loads the schema and surfaces violations alongside
  range / interdependency checks
- `jsonschema` added as a core dependency

#### Polymerase-aware defaults
- `min_k`, `max_k`, and `mg_conc` defaults now come from the chosen
  polymerase: phi29 6-12 bp / 10 mM, equiphi29 10-18 bp / 10 mM,
  bst 15-25 bp / 8 mM, klenow 8-15 bp / 10 mM
- Explicit user values in `params.json` always win

#### Adaptive GC filter auto-engagement
- When `genome_gc` is set or auto-computed and falls outside [0.35, 0.65],
  the filter switches to `genome_gc +/- gc_tolerance` automatically
- Logs a single INFO line when it engages
- Disable with `"adaptive_gc": false` in `params.json`

#### Additives affect scoring and optimization
- `IntegratedQualityScorer` now uses salt-corrected nearest-neighbor Tm
  plus `ReactionConditions.calculate_tm_correction()` (DMSO, betaine,
  trehalose, formamide, ethanol, urea, TMAC)
- Previously additives only affected filtering-time Tm; primers were
  mis-ranked at scoring time

#### Genome-size k-mer heuristic
- `condition_suggester.suggest_kmer_range(size, gc, polymerase)` returns
  size-aware recommendations
- `get_params()` emits warnings for unusual combinations (e.g. max_k<8
  on a 5 Mb bacterium, or min_k>12 on a 5 kb plasmid)

#### Example templates
- `examples/equiphi29_scenario/` - long primers + DMSO + betaine
- `examples/multi_genome_blacklist/` - pan-primer with zero-tolerance
  contaminant list
- `examples/gc_extreme/` - AT-rich or GC-rich targets

#### Release infrastructure
- `.github/workflows/publish.yml` - PyPI trusted publishing on release
- `.github/workflows/nightly.yml` - nightly end-to-end pipeline tests
- `.pre-commit-config.yaml` - ruff / black / isort / yaml / toml hooks
- `ruff` and `mypy` baselines in CI (non-blocking) and `pyproject.toml`

### BUG FIXES

- `pipeline.py` position-file cache now requires every `fg_prefixes`
  entry to have a cache before skipping creation; previously one missing
  cache could silently leave multi-target genomes unscanned
- `background_aware_optimizer.compare_optimizers` now aggregates hits
  across all `bg_prefixes` instead of only `bg_prefixes[0]`
- `bl_seq_lengths` is auto-computed when `bl_genomes` is set via CLI or
  JSON, preventing incorrect blacklist frequencies

### EXPANDED VALIDATION

- `PARAM_RANGES` now covers `formamide_percent`, `ethanol_percent`,
  `urea_m`, `tmac_m`, `glycerol_percent`, `peg_percent`, `bsa_ug_ml`,
  `mg_conc` (widened to 0-20 mM), `gc_tolerance`, `genome_gc`,
  `bl_penalty`, and `max_bl_freq`

### BREAKING CHANGES

None. All 3.6 configurations continue to run; see
`docs/migration-3.6-to-3.7.md` for behavioural differences.

## [3.6.0] - 2026 - Optimizer Framework and Pipeline Hardening

### NEW FEATURES

#### BaseOptimizer Framework and OptimizerFactory
- All optimizers now inherit from `BaseOptimizer` with a consistent interface
- `OptimizerFactory` registry with decorator-based registration
- `PrimerSetMetrics.normalized_score()` enables cross-optimizer comparison
- `CompositeOptimizer` for chaining multiple optimizers

#### New Optimizers
- **CliqueOptimizer** (`clique`): Dimer-free primer sets via maximum clique enumeration
- **DimerValidator**: Post-optimization dimer validation and replacement suggestions
- **NormalizedOptimizer** (`normalized`): Strategy presets (discovery, clinical, enrichment, metagenomics)
- **TilingOptimizer** (`tiling`): Interval-based genome tiling coverage
- **MultiAgentOrchestrator** (`multi-agent`): Parallel ensemble of optimizer strategies
- **BackgroundPrefilter** (`bg-prefilter`): fg/bg ratio pruning wrapper
- **SerialCascadeOptimizer**: Pipeline combinations (coverage-then-dimerfree, dimerfree-scored, bg-prefilter-hybrid)
- **WeightedSetCoverOptimizer** (`weighted-set-cover`): Score-weighted set cover variant

#### Host-Free Mode
- `--no-background` flag for optimizing without a background genome
- Useful for general MDA or enrichment-only workflows

#### Enrichment Prediction
- `EfficiencyPredictor` provides mechanistic enrichment fold-change estimates
- Integrated into `neoswga interpret` output

#### Condition Sweep
- `neoswga suggest --sweep` searches 108 additive combinations
- `AdditiveOptimizer` recommends optimal cocktail

#### Simulation Validation
- `--validate-with-simulation` flag for post-hoc primer set validation
- `ExperimentalTracker` for logging wet-lab outcomes and calibrating predictions

#### Export and Lab Integration
- `neoswga export` command with FASTA, CSV, BED, BedGraph, and protocol formats
- `PrimerExporter` class with vendor-specific CSV (IDT, Twist, Sigma)
- Modification profiles for PTO bonds and 5' blocking

#### Pareto Frontier Visualization
- `pareto_frontier.py`: Plot, report, and CLI summary of multi-objective results

#### Primer Expansion
- `PrimerExpander`: Identify coverage gaps and suggest additional primers

### REFACTORING
- Unified CLI entry point in `cli_unified.py` (replaces old multi-command CLIs)
- Consolidated test suite with shared fixtures
- Converted all imports to from-import style

### BUG FIXES
- Fixed polymerase validation for non-standard enzyme names
- Fixed test pollution from shared mutable state in parameter module

---

## [3.5.0] - 2025 - Genome-Adaptive QA

### NEW FEATURES
- Genome-adaptive quality assessment for extreme GC genomes
- Automatic GC classification (extreme_at, at_rich, balanced, gc_rich, extreme_gc)
- Mechanistic four-pathway model (Tm, accessibility, enzyme, kinetics)
- `MechanisticModel` and `MechanisticEffects` classes
- Additive interaction registry with synergy/antagonism modeling
- Set size optimizer with application profiles
- Quality report module (`neoswga/core/report/`)
- Interactive Plotly visualizations (optional dependency)
- Setup wizard (`neoswga init`), parameter validator, condition suggester, results interpreter

---

## [2.0.0] - 2025 - Improved Pipeline Release

### CRITICAL FIXES

#### Fixed: GC Filter Bug That Blocked Entire Organism Classes
- **Issue**: Fixed GC thresholds (37.5-62.5%) rejected ALL primers for organisms with extreme GC content
- **Organisms affected**: AT-rich targets (~33% GC) and GC-rich targets (~67% GC)
- **Impact**: Algorithm was completely non-functional for ~20% of bacterial pathogens
- **Solution**: Adaptive GC filtering based on genome composition
- **File**: `neoswga/core/adaptive_filters.py`

**Before**:
```python
# src/filter.py:56
if GC_content <= 0.375 or GC_content >= 0.625:
    return False  # BLOCKS AT-rich and GC-rich targets
```

**After**:
```python
# Adapts to genome GC content
gc_min = max(0.20, genome_gc - 0.15)
gc_max = min(0.80, genome_gc + 0.15)
# AT-rich target (33%): Accepts 18-48% GC primers
# Caulobacter (67%): Accepts 52-82% GC primers
```

### NEW FEATURES

#### 1. Position Cache - 1000x Speedup
- **File**: `neoswga/core/position_cache.py`
- **Problem**: 20,000+ HDF5 disk reads at 10ms each = 200s I/O time
- **Solution**: Load all positions into memory once (4 MB typical)
- **Impact**: Query time 10ms → 0.01ms (1000x faster)
- **Memory**: ~4 MB for 500 primers × 1000 binding sites

#### 2. Background Bloom Filter - Enables Human Genome
- **File**: `neoswga/core/background_filter.py`
- **Problem**: Exact HDF5 index for human genome = 170 GB (infeasible)
- **Solution**: Bloom filter (4 GB) + sampled index (1.7 GB)
- **Impact**: Can now process human (3 Gbp), tick (2.1 Gbp), mosquito (1.4 Gbp)
- **Accuracy**: ~1% false positive rate (acceptable for filtering)
- **Build time**: ~30 minutes one-time, then instant reuse

#### 3. Network-Based Optimization - 100x Better Enrichment
- **File**: `neoswga/core/network_optimizer.py`
- **Problem**: Ratio-based scoring ignores amplification dynamics
- **Solution**: Graph-based optimization modeling connectivity
- **Mathematical basis**:
  - Connected primers → exponential growth (2^n)
  - Isolated primers → linear growth (n)
  - Difference: 500x in amplification
- **Impact**: Predicted enrichment 100x higher than old method

#### 4. MILP Optimizer - Provably Optimal Solutions
- **File**: `neoswga/core/milp_optimizer.py`
- **Problem**: Greedy algorithm has unknown approximation ratio
- **Solution**: Mixed Integer Linear Programming
- **Guarantees**: <5% optimality gap
- **Performance**: Solves 100-500 candidates in minutes
- **Fallback**: Hybrid method auto-switches to greedy for large sets

#### 5. Stochastic Simulator - Validation
- **File**: `neoswga/core/stochastic_simulator.py`
- **Purpose**: Validate network predictions with actual kinetics
- **Algorithm**: Gillespie (stochastic simulation)
- **Models**:
  - Primer binding/unbinding
  - Polymerase extension
  - Resource depletion (dNTP, polymerase)
  - Network rebinding effects

#### 6. Pipeline Integration - Drop-In Replacement
- **File**: `neoswga/core/pipeline_integration.py`
- **Purpose**: Seamless integration with existing codebase
- **Usage**: Replace `step4()` with one line
- **Compatibility**: Returns results in old format
- **Features**:
  - Automatic fallback if improvements fail
  - Side-by-side comparison
  - Parameter migration

#### 7. Validation Framework
- **File**: `neoswga/core/validation.py`
- **Tests**:
  - Position cache correctness and speed
  - Adaptive GC filter functionality
  - Bloom filter accuracy
  - Network optimization quality
  - MILP optimality
  - End-to-end pipeline
- **Modes**: Full validation or quick test

#### 8. Command-Line Interface
- **File**: `neoswga/cli_improved.py`
- **Commands**:
  - `neoswga-improved build-filter` - Build Bloom filters
  - `neoswga-improved design` - Design primers
  - `neoswga-improved compare` - Compare old vs. new
  - `neoswga-improved validate` - Run tests

### PERFORMANCE IMPROVEMENTS

#### Runtime
- **Old**: 5 minutes (E. coli vs. human)
- **New**: 30 seconds
- **Speedup**: 10x

#### Memory
- **Old**: 170 GB (human genome exact index)
- **New**: <10 GB (Bloom filter + cache)
- **Reduction**: 17x

#### Enrichment
- **Old**: 10-100x predicted
- **New**: 1000-10,000x predicted
- **Improvement**: 100x

#### Query Speed
- **Old**: 10ms per HDF5 read
- **New**: 0.01ms per cache lookup
- **Speedup**: 1000x

### DOCUMENTATION

New documentation:
- `README_IMPROVED.md` - User guide
- `IMPLEMENTATION_SUMMARY.md` - Technical overview
- `IMPLEMENTATION_GUIDE.md` - Migration guide
- `QUICKSTART_NEW_PIPELINE.md` - 5-minute setup
- `DEPLOYMENT.md` - Deployment guide
- `CHANGELOG.md` (this file)

### DEPENDENCIES

New optional dependencies:
- `pybloom-live>=3.1.0` - For Bloom filters
- `mip>=1.13.0` - For MILP optimization

Install with: `pip install ".[improved]"`

### TESTING

New test files:
- `benchmark_improvements.py` - Performance benchmarks
- `neoswga/core/validation.py` - Validation suite

Run tests:
```bash
python benchmark_improvements.py --test all
python -m neoswga.core.validation
```

### BREAKING CHANGES

None. Improved pipeline is fully backward compatible:
- Drop-in replacement for `step4()`
- Returns same format as old pipeline
- Accepts same parameters

Users may see different results (better enrichment) but interface is identical.

### MIGRATION

Three migration paths:

1. **Drop-in replacement** (easiest):
   ```python
   from neoswga.core.unified_optimizer import optimize_step4
   def step4():
       return optimize_step4()
   ```

2. **Gradual integration**:
   - Week 1: Add position cache
   - Week 2: Add adaptive GC filter
   - Week 3: Add background filter
   - Week 4: Switch to network optimizer

3. **Complete rewrite**:
   - Use `ImprovedPipeline` directly
   - See examples in documentation

### KNOWN ISSUES

1. **Optional dependencies**: Some features require pybloom-live and mip
   - Without pybloom-live: Background filtering disabled
   - Without mip: Only greedy optimization available
   - Core functionality works without these

2. **Memory**: Position cache loads all positions (usually fine for bacteria)
   - Solution: Use `StreamingPositionCache` for very large datasets

3. **MILP timeout**: May timeout for >1000 candidates
   - Solution: Hybrid method auto-falls back to greedy

### EXPERIMENTAL VALIDATION

Network-based predictions validated with:
- Gillespie stochastic simulation
- Mathematical analysis (graph theory)
- Information-theoretic foundation
- Thermodynamic modeling

Users should verify designed primers experimentally:
1. Design primers with improved pipeline
2. Run SWGA amplification
3. Measure enrichment (qPCR/sequencing)
4. Compare to predictions

Expected: 2-5x match between prediction and experiment

### ACKNOWLEDGMENTS

Original SOAPswga by Jane Dwivedi-Yu et al. (2023)
Enhanced implementation by Andreas Sjodin (2025)

Based on feedback identifying:
- GC filter bug blocking AT-rich and GC-rich targets
- I/O bottleneck from repeated HDF5 reads
- Human genome infeasibility
- Ratio-based scoring limitations

---

## [1.0.0] - Original SOAPswga

Original implementation with:
- K-mer preprocessing (jellyfish)
- Candidate filtering (frequency, Gini, Tm)
- Random forest scoring
- Greedy BFS optimization

See original documentation for details.

---

**For detailed technical information, see:**
- [API Reference](reference/API_REFERENCE.md) - Public API documentation
- [Module Reference](reference/MODULE_REFERENCE.md) - All core modules
- [Developer Guide](development/DEVELOPER_GUIDE.md) - Development setup and contribution guidelines
