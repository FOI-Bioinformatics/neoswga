# E. coli Environmental Example: Long Primers, Broad Coverage

Designs 15 bp SWGA primers for **E. coli K-12 MG1655** (NC_000913.3, 4.64 Mb),
targeting broad recovery of the chromosome. Intended as a starting point for
amplifying E. coli DNA in environmental samples where the exact non-target
community is not specified in advance.

## Design choices

- **Polymerase: equiphi29 at 43 C.** Required for 15 bp primers; phi29 at 30 C is
  tuned for 6-12 bp primers.
- **Primer length: 15 bp (fixed).** A random 15-mer appears roughly once per 1 Gb
  of random DNA (4^15 approx 1e9 possible sequences), so 15-mers are
  intrinsically rare in foreign genomes. This is the only specificity mechanism
  used here. A wider range (e.g. 15-18 bp) is possible in principle but becomes
  compute-heavy because the pipeline auto-scales `min_fg_freq` by primer length
  (see `neoswga/core/filter.py:90`), admitting nearly every unique long k-mer
  as a candidate and inflating position-caching cost.
- **No background or blacklist.** `bg_genomes` is empty. The pipeline does not
  evaluate or filter against any non-target sequence.
- **Optimization: dominating-set (CLI flag).** A graph-based set-cover
  formulation with an ln(n) approximation. Well matched to a pure-coverage
  objective and 1000x faster than the default `hybrid` for this pool size.
- **Tight `min_fg_freq: 2e-3`.** Scales down to a lower bound of roughly 9
  occurrences per primer in the genome after the length-scaling in filter.py,
  yielding a manageable candidate pool (a few hundred primers) instead of the
  ~4 million unique 15-mers present in E. coli.
- **`min_amp_pred: 0.0`.** The pretrained random-forest amplification predictor
  is calibrated on 6-12 bp primers; 15-mers systematically score low. Setting
  the threshold to 0 retains every scored primer and lets the optimizer handle
  selection.

## Specificity caveat

This setup does not measure cross-amplification against any defined community.
Primer length reduces random-match probability but does not rule out hits to
abundant non-target k-mers in real environmental samples. Empirical validation
is the user's responsibility. If you have a relevant host or contaminant
genome (human, common soil/water taxa), add it under `bg_genomes` and set
`max_bg_freq` to enable background filtering.

## Setup

Fetch the reference chromosome once:

```bash
cd examples/ecoli_environmental
sh download_genome.sh
```

This writes `ecoli_K12_MG1655.fasta` (~4.5 MB). The file is gitignored.

## Running the pipeline

```bash
neoswga count-kmers -j params.json
neoswga filter      -j params.json
neoswga score       -j params.json
neoswga optimize    -j params.json --optimization-method=dominating-set
```

Note the explicit `--optimization-method` flag: it is a CLI option, not a
`params.json` key, and defaults to `hybrid` if omitted.

## Expected runtime

| Step | Time (4 cores) |
|------|---------------|
| count-kmers | a few seconds |
| filter | about 90 seconds |
| score | about 15-20 minutes (random-forest feature computation dominates) |
| optimize | under 1 second |

## Expected output

| Step | Key output file | Notes |
|------|-----------------|-------|
| 1 | `ecoli_K12_MG1655_15mer_all.txt` | k-mer counts from jellyfish |
| 2 | `step2_df.csv` | filtered candidates (~600 primers expected) |
| 3 | `step3_df.csv` | scored candidates with `amp_pred` column |
| 4 | `step4_improved_df.csv`, `step4_improved_df_summary.json` | final set |

A representative run with the committed parameters selects 6 primers with:

- `fg_coverage` approximately 0.997 (optimizer-reported coverage of the
  genome using the effective per-primer reach window)
- `per_target_coverage` approximately 0.24 (a stricter bin-based metric)
- Mean Tm approximately 45 C (range 42-48 C)
- Max gap between adjacent binding sites approximately 150 kb
- Coverage of 463/465 genomic bins at the default bin size

The two coverage numbers reflect different assumptions about per-primer reach
and are both reported for transparency. The optimizer uses `fg_coverage` as
its objective.

## Reviewing results

```bash
neoswga interpret -d ./
neoswga report    -d ./ --interactive
```

## Alternative: automatic set sizing

Instead of fixing `num_primers` / `target_set_size` in `params.json`, you can
let the pipeline derive a size from the `metagenomics` application profile
(target coverage 0.95):

```bash
neoswga optimize -j params.json --optimization-method=dominating-set \
                 --auto-size --application metagenomics
```
