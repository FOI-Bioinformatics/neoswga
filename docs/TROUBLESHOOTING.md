# NeoSWGA Troubleshooting Guide

This guide covers common issues encountered when installing and running NeoSWGA,
with symptoms, causes, and step-by-step solutions.

---

## Installation Issues

### Jellyfish not found

**Symptom:**

```
ERROR: Jellyfish is required but not found in PATH.
Jellyfish is an external k-mer counting tool that must be installed separately.
```

**Cause:** Jellyfish is an external C++ tool for k-mer counting. It is not a
Python package and must be installed separately from NeoSWGA.

**Solution:**

1. Install Jellyfish using one of the following methods:
   ```bash
   # Via conda (recommended)
   conda install -c bioconda jellyfish

   # Via Homebrew (macOS)
   brew install jellyfish

   # From source
   # https://github.com/gmarcais/Jellyfish
   ```
2. Verify installation: `jellyfish --version`
3. Ensure the binary is on your PATH. If installed via conda, make sure the
   correct conda environment is activated.

### Jellyfish version 1.x detected

**Symptom:**

```
ERROR: Jellyfish version 1.x.x detected. NeoSWGA requires Jellyfish 2.x.
```

**Cause:** NeoSWGA requires Jellyfish 2.x because version 1.x has an
incompatible command-line interface.

**Solution:**

1. Remove the old version: `conda remove jellyfish` or `brew uninstall jellyfish`
2. Install version 2.x: `conda install -c bioconda 'jellyfish>=2'`
3. Confirm: `jellyfish --version` should report 2.x.x

### Python version incompatibility

**Symptom:**

```
ERROR: Requires-Python >=3.11
```

Or import errors referencing missing standard library features.

**Cause:** NeoSWGA requires Python 3.11 or later. Earlier versions are not
supported.

**Solution:**

1. Check your Python version: `python --version`
2. If below 3.11, install a newer version via conda, pyenv, or your system
   package manager.
3. Create a dedicated environment:
   ```bash
   conda create -n neoswga python=3.12
   conda activate neoswga
   pip install -e .
   ```

### scikit-learn model loading fails

**Symptom:**

```
ModuleNotFoundError: No module named 'sklearn'
```

Or unpickling errors when the score command loads the random forest model:

```
ValueError: node array from the pickle has an incompatible dtype
```

**Cause:** The pre-trained random forest model in
`neoswga/core/models/random_forest_filter.p` was serialized with a specific
scikit-learn version. Major version changes can break deserialization.

**Solution:**

1. Ensure scikit-learn is installed within the supported range: `pip install 'scikit-learn>=1.0,<2'`
2. If the model still fails to load, retrain it:
   ```bash
   python scripts/retrain_rf_model.py --output neoswga/core/models/random_forest_filter.p
   ```
3. If the retrain script is not available, reinstall NeoSWGA from source:
   ```bash
   pip install -e . --force-reinstall
   ```

### Optional dependencies not installed

**Symptom (Bloom filter):**

```
WARNING: pybloom_live not installed. Install with: pip install pybloom-live
```

**Symptom (MILP optimizer):**

```
ImportError: python-mip required. Install: pip install mip
```

**Symptom (interactive charts):**

Report generation skips interactive charts silently, or:

```
Report module not available
```

**Cause:** NeoSWGA has several optional dependency groups that are not installed
by default.

**Solution:**

Install the optional groups you need:

```bash
# Bloom filter for large backgrounds + MILP optimizer
pip install -e ".[improved]"

# Interactive Plotly charts in reports
pip install -e ".[interactive]"

# GPU acceleration (CUDA 11.x)
pip install -e ".[gpu]"

# GPU acceleration (CUDA 12.x)
pip install -e ".[gpu-cuda12]"

# All common extras
pip install -e ".[improved,interactive]"
```

---

## Pipeline Errors

### "No candidates remaining" after filtering

**Symptom:**

```
NoCandidatesError: No candidates remaining from <N> primers after <stage>
```

Or the filter step produces a `step2_df.csv` with zero rows.

**Cause:** Filtering thresholds are too stringent for the genome pair. This
commonly occurs with:

- Very similar target and background genomes (few discriminating k-mers)
- Extreme GC content genomes
- Very short or very long primer length ranges

**Solution:**

1. Relax filtering parameters in `params.json`:
   ```json
   {
     "min_fg_freq": 5e-6,
     "max_bg_freq": 1e-5,
     "max_gini": 0.8,
     "max_primer": 1000
   }
   ```
2. Widen the primer length range (e.g., `"min_k": 6, "max_k": 14`).
3. For extreme GC genomes, NeoSWGA applies adaptive GC filtering automatically.
   Check the log output for adaptive filter messages.
4. Validate your parameters before running:
   ```bash
   neoswga validate-params -j params.json
   ```
5. As a diagnostic step, run filter with `--verbose` to see how many candidates
   are removed at each stage.

### Running pipeline steps out of order

**Symptom (filter without count-kmers):**

```
ERROR: Required file not found: <path>
Ensure Step 1 (count-kmers) has completed successfully.
Run: neoswga count-kmers -j params.json
```

**Symptom (score without filter):**

```
ERROR: Required file not found: <path>
Ensure Step 2 (filter) has completed successfully.
Run: neoswga filter -j params.json
```

**Symptom (optimize without score):**

```
ERROR: Required file not found: <path>
Ensure Step 3 (score) has completed successfully.
Run: neoswga score -j params.json
```

**Cause:** Each pipeline step depends on output files from the previous step.
The steps must run in order: count-kmers, filter, score, optimize.

**Solution:**

Run the full pipeline in sequence:

```bash
neoswga count-kmers -j params.json
neoswga filter -j params.json
neoswga score -j params.json
neoswga optimize -j params.json
```

Or use the single-command pipeline:

```bash
neoswga design -j params.json
```

### Missing or invalid params.json

**Symptom (file not found):**

```
Error: parameter file 'params.json' not found.
```

**Symptom (invalid JSON):**

```
Error: 'params.json' is not valid JSON: Expecting property name ...
```

**Symptom (missing required fields):**

```
Error: 'params.json' is missing required field 'data_dir'.
```

```
Error: 'params.json' is missing genome specification.
```

**Cause:** The parameter file is missing, contains syntax errors, or lacks
required fields (`data_dir` and a genome specification).

**Solution:**

1. Use the setup wizard to create a valid configuration:
   ```bash
   neoswga init --genome target.fasta --background host.fasta
   ```
2. At minimum, `params.json` must contain:
   ```json
   {
     "fg_genomes": ["target_genome.fasta"],
     "bg_genomes": ["background_genome.fasta"],
     "data_dir": "./data/"
   }
   ```
3. Validate your file before running: `neoswga validate-params -j params.json`

### Invalid genome FASTA file

**Symptom:**

```
GenomeFileError: Cannot read genome file: <path>
```

Or unexpected characters in sequences during filtering.

**Cause:** The FASTA file is malformed, empty, or contains non-standard
characters.

**Solution:**

1. Verify the file is a valid FASTA: it should start with `>` followed by a
   header line, then sequence lines containing only A, C, G, T, and N.
2. Check for common issues:
   - Windows line endings (use `dos2unix` to convert)
   - Compressed files (decompress `.gz` files first)
   - Empty files or files with zero-length sequences
3. Confirm the file path in `params.json` is correct and the file is readable.

### Memory errors on large genomes

**Symptom:**

```
MemoryError
```

Or the process is killed by the OS, or Jellyfish reports insufficient memory.

**Cause:** Loading all k-mers from a large background genome (e.g., human at
3 Gbp) into memory exceeds available RAM. The filter step is particularly
memory-intensive.

**Solution:**

1. Use the Bloom filter for large backgrounds:
   ```bash
   neoswga build-filter background_genome.fna ./
   ```
   Then reference the filter in your pipeline run. The Bloom filter uses
   approximately 1% of the memory of a full k-mer table.
2. Reduce the k-mer range to lower memory usage (e.g., `"min_k": 8, "max_k": 12`).
3. Increase available swap space on the system.
4. See the System Requirements section below for RAM estimates.

### Jellyfish timeout or failure

**Symptom:**

```
jellyfish count failed for k=<N>: <error message>
```

Or the pipeline hangs during count-kmers.

**Cause:** Jellyfish has a default timeout of 1 hour per k-value. Very large
genomes or insufficient system resources can cause timeout. A small default hash
size can also force repeated resizing, which slows execution.

**Solution:**

1. Ensure sufficient disk space for temporary Jellyfish files (approximately
   2-4x the genome size).
2. Run with fewer k-values to reduce total time:
   ```bash
   neoswga count-kmers -j params.json --min-k 8 --max-k 12
   ```
3. Verify Jellyfish works independently:
   ```bash
   jellyfish count -m 10 -s 100M -t 4 -C genome.fasta -o test.jf
   jellyfish dump -c test.jf | head
   ```
4. If the hash size is too small, NeoSWGA auto-estimates it from the genome file
   size. For very large genomes, ensure at least 8 GB of RAM is available.

---

## Optimization Issues

### MILP optimizer not available

**Symptom:**

```
ImportError: MILP optimizer requires python-mip package. Install: pip install mip
```

**Cause:** The `mip` package is an optional dependency not included in the
base installation.

**Solution:**

```bash
pip install mip
# Or install the improved extras group:
pip install -e ".[improved]"
```

If you do not need exact solutions, use the default `hybrid` optimizer or the
faster `dominating-set` method instead.

### Optimizer returns fewer primers than requested

**Symptom:** The optimized set contains fewer primers than the
`num_primers` / `target_set_size` parameter.

**Cause:** The candidate pool after scoring may be smaller than the requested
set size, or the optimizer could not find enough primers that meet the quality
constraints (coverage, dimer avoidance).

**Solution:**

1. Increase the candidate pool by relaxing filter thresholds (see "No candidates
   remaining" above).
2. Increase `max_primer` in `params.json` to retain more candidates through
   filtering.
3. Lower `target_set_size` to a value achievable with the available candidates.
4. Check the filter and score output files (`step2_df.csv`, `step3_df.csv`) to
   confirm enough candidates are available.

### Non-reproducible results

**Symptom:** Running the same pipeline twice produces different primer sets.

**Cause:** Several optimizers (genetic, hybrid, moea) use stochastic search.
Without a fixed random seed, results vary between runs.

**Solution:**

Set the `--seed` flag during optimization:

```bash
neoswga optimize -j params.json --seed 42
```

Deterministic optimizers (`greedy`, `dominating-set`, `milp`) produce
reproducible results without a seed.

### Optimization takes too long

**Symptom:** The optimize step runs for hours without completing.

**Cause:** Optimization time depends on the number of candidates, the requested
set size, and the optimization method. Methods such as `background-aware` and
`moea` are computationally expensive.

**Solution:**

1. Use a faster optimization method:
   ```bash
   # Fast graph-based (8x faster than greedy)
   neoswga optimize -j params.json --optimization-method=dominating-set

   # Standard balanced approach
   neoswga optimize -j params.json --optimization-method=hybrid
   ```
2. Reduce the candidate pool by lowering `max_primer` in `params.json`.
3. Reduce `iterations` (default: 8) for faster but potentially less optimal
   results.
4. Reduce `max_sets` (default: 5) to build fewer parallel primer sets.

---

## Result Quality Issues

### Low coverage score

**Symptom:** The report or interpreter indicates low genome coverage (below 50%),
or the quality grade is penalized for coverage.

**Cause:** The selected primer set does not bind frequently enough across the
target genome. This can result from a small primer set, a large target genome,
or primers with uneven binding distributions.

**Solution:**

1. Increase the primer set size:
   ```bash
   neoswga optimize -j params.json --target-set-size 12
   ```
2. Use `--auto-size` to let NeoSWGA recommend a set size for your application:
   ```bash
   neoswga optimize -j params.json --auto-size --application discovery
   ```
3. Lower the `max_gini` threshold in `params.json` to favor more evenly
   distributed primers (e.g., `"max_gini": 0.5`).
4. Widen the primer length range to increase the candidate pool.

### Poor enrichment ratio

**Symptom:** The enrichment ratio (target vs. background amplification) is below
5x, or the report flags poor specificity.

**Cause:** The primers bind too frequently to the background genome relative to
the target. This is common when the target and background share substantial
sequence homology.

**Solution:**

1. Use the `background-aware` optimizer, which explicitly minimizes background
   binding:
   ```bash
   neoswga optimize -j params.json --optimization-method=background-aware
   ```
2. Tighten the background frequency filter: set `max_bg_freq` to a lower value
   (e.g., `1e-6`).
3. Build a Bloom filter for the background genome to improve background
   frequency estimates:
   ```bash
   neoswga build-filter background_genome.fna ./
   ```
4. Consider using a longer primer length range (12-18 bp with EquiPhi29) for
   greater specificity.

### High dimer risk

**Symptom:** The interpreter or report warns about high dimer risk among the
selected primers.

**Cause:** Primers in the set have complementary 3' ends that can form
primer-dimer products, competing with genome amplification.

**Solution:**

1. Use the `network` optimizer, which applies a dimer penalty during selection:
   ```bash
   neoswga optimize -j params.json --optimization-method=network
   ```
2. Run dimer analysis on your primer set to identify problematic pairs:
   ```bash
   neoswga analyze-dimers --primers SEQ1 SEQ2 SEQ3 --output dimers/
   ```
3. Manually remove the most problematic primer and re-optimize with a larger set
   to compensate.
4. Increase the candidate pool so the optimizer has more dimer-free options.

### Low report grade (D or F)

**Symptom:**

```
Quality Grade: D
```

or

```
Quality Grade: F
```

**Cause:** The overall quality grade is a weighted composite of coverage,
enrichment, primer count, Tm uniformity, and dimer risk. A D or F grade
indicates that multiple quality components scored poorly.

**Solution:**

1. Generate a full technical report to identify which components are weak:
   ```bash
   neoswga report -d results/ --level full
   ```
2. Run the interpreter for targeted recommendations:
   ```bash
   neoswga interpret -d results/
   ```
3. Address the weakest components first (see the specific sections above for
   coverage, enrichment, and dimer issues).
4. Common improvements:
   - Increase primer set size for better coverage
   - Use `background-aware` optimization for better enrichment
   - Narrow the Tm range in `params.json` for better Tm uniformity
   - Use `network` optimization for lower dimer risk

---

## System Requirements

### Supported platforms

NeoSWGA runs on Linux and macOS. Windows is not officially supported but may
work under WSL2. Python 3.11 or later is required.

### Resource estimates by genome size

The following are approximate resource requirements. Actual usage depends on
k-mer range, number of candidates, and optimization method.

| Genome size | RAM (filter step) | Disk space | Time (full pipeline) |
|-------------|-------------------|------------|----------------------|
| Plasmid (5-50 kb) | < 1 GB | < 100 MB | < 1 minute |
| Bacterial (1-10 Mb) | 1-4 GB | 100 MB - 1 GB | 5-30 minutes |
| Parasite (10-100 Mb) | 4-8 GB | 1-5 GB | 30 min - 2 hours |
| Human background (3 Gb) | 16-32 GB* | 10-30 GB | 1-6 hours |

*Use a Bloom filter (`neoswga build-filter`) to reduce memory to 2-4 GB for
large background genomes.

### Minimum requirements

- **CPU:** 2 cores (4+ recommended for Jellyfish parallelism)
- **RAM:** 2 GB minimum; 8 GB recommended for bacterial genomes
- **Disk:** 1 GB free space plus 2-4x the genome file sizes for intermediate files

### Quick installation check

```bash
neoswga validate --quick
```

This verifies that NeoSWGA is installed correctly, Jellyfish is available, and
core dependencies are functional.
