# E. coli Environmental Example: Long Primers, Broad Coverage

Designs 17 bp SWGA primers for **E. coli K-12 MG1655** (NC_000913.3, 4.64 Mb),
targeting broad recovery of the chromosome. Intended as a starting point for
amplifying E. coli DNA in environmental samples where the exact non-target
community is not specified in advance.

## Design choices

- **Polymerase: equiphi29 at 43 C.** Required for 17 bp primers; phi29 at 30 C
  is tuned for 6-12 bp primers.
- **Primer length: 17 bp (fixed).** A random 17-mer appears roughly once per
  16 Gb of random DNA (4^17 approx 1.7e10 possible sequences), so 17-mers are
  very unlikely to occur in foreign genomes by chance. This is the only
  specificity mechanism used here.
- **No background or blacklist.** `bg_genomes` is empty. The pipeline does not
  evaluate or filter against any non-target sequence.
- **Optimization: dominating-set with auto-size.** A graph-based set-cover
  formulation with an ln(n) approximation, sized by the `metagenomics`
  application profile (target coverage 0.95, set size 15-20 primers). The
  dominating-set greedy left on its own stops when marginal coverage falls
  to zero, which under-sizes the set on genomes with patchy candidate
  distributions; auto-size forces the optimizer to keep picking primers
  until the profile's target size is reached.
- **`min_fg_freq: 0.015`.** The frequency threshold in `params.json` is
  auto-scaled by primer length as `base * 4^(10 - k)`, so at 17 bp the
  effective threshold is about 9e-7 and requires a primer to appear roughly
  5 or more times in the chromosome. This keeps the candidate pool at
  around 1000 primers (after the `max_primer` cap).
- **`min_amp_pred: 0.0`.** The pretrained random-forest amplification predictor
  is calibrated on 6-12 bp primers and systematically under-scores 17-mers.
  Setting the threshold to 0 retains every scored primer and lets the
  optimizer handle selection. The RF score is advisory and not read by any
  optimizer, so this change has no effect on the selected set.

## Specificity caveat

This setup does not measure cross-amplification against any defined community.
Primer length reduces random-match probability but does not rule out hits to
abundant non-target k-mers in real environmental samples. Empirical validation
is the user's responsibility. If you have a relevant host or contaminant
genome (human, common soil/water taxa), add it under `bg_genomes` and set
`max_bg_freq` to enable background filtering.

## Tm caveat

17-mers at 40-60% GC have a mean Tm around 52-56 C, which is above the 43 C
reaction temperature. In practice this means strong binding and slower primer
turnover on the template. Equiphi29 tolerates this range, but the setup is at
the upper edge of the polymerase's comfortable window. If this matters for
your application, you can either lower `max_tm` in `params.json` to prefer
cooler primers, or drop to 15 bp and re-tune `min_fg_freq` (see the
auto-scaling note above).

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
neoswga optimize    -j params.json --optimization-method=dominating-set \
                    --auto-size --application metagenomics
```

Note the explicit flags on `optimize`:

- `--optimization-method=dominating-set` is a CLI option, not a `params.json`
  key, and defaults to `hybrid` if omitted.
- `--auto-size --application metagenomics` tells the optimizer to size the
  primer set from the `metagenomics` profile (target coverage 95%, typical
  size 15-20 primers) rather than stopping when marginal coverage is zero.

## Expected runtime

| Step | Time (4 cores) |
|------|---------------|
| count-kmers | a few seconds |
| filter | about 90 seconds |
| score | a few seconds (thermodynamic histogram features are skipped by default; add `--full-score` to include them, which takes 15-20 minutes) |
| optimize | under 1 second |

## Expected output

| Step | Key output file | Notes |
|------|-----------------|-------|
| 1 | `ecoli_K12_MG1655_17mer_all.txt` | k-mer counts from jellyfish |
| 2 | `step2_df.csv` | filtered candidates (around 1000 primers after the `max_primer` cap) |
| 3 | `step3_df.csv` | scored candidates with `amp_pred` column |
| 4 | `step4_improved_df.csv`, `step4_improved_df_summary.json` | final set |

A representative run with the committed parameters selects 20 primers with:

- `fg_coverage` approximately 0.25 (union of 3 kb reach windows around each
  primer binding site, summed across both strands, as a fraction of the
  4.64 Mb chromosome)
- `per_target_coverage` approximately 0.32 (same idea but computed per target
  genome)
- Graph-level coverage 250/273 of 10 kb bins (91.6%) — the proportion of
  10 kb chromosomal windows with at least one primer binding site within
  3 kb of that window's start
- Mean Tm approximately 52 C (range 48-56 C)
- Max gap between adjacent binding sites approximately 145 kb
- Dimer risk score approximately 0.15

## Understanding the coverage numbers

The setup reports three coverage numbers that measure different things, and
they are all legitimately different:

- **Graph-level bin coverage (91.6%, 250/273 regions).** Did every 10 kb
  stretch of the chromosome get a primer binding site somewhere within 3 kb
  of it? This is the optimizer's own objective.
- **fg_coverage (0.25).** Take the union of 3 kb-radius windows around every
  selected-primer binding site and divide by genome length. This is the
  realistic "what fraction of the chromosome is within per-primer reach" in a
  dense SWGA reaction, per Clarke et al. (2017) and Dwivedi-Yu et al. (2023).
- **per_target_coverage (0.32).** Same union-of-windows idea, computed per
  target genome, reported for multi-genome runs.

Earlier versions of neoswga reported `fg_coverage` using phi29's 70 kb
single-molecule processivity as the reach, which inflated the number by
5-20x over what the set can actually amplify (see the Phase 16 note in
`neoswga/core/coverage.py`). The current build uses the realistic 3 kb
(phi29) / 4 kb (equiphi29) reach by default; the ~0.25 above is the honest
number. If you need the legacy behavior for comparison, it is still
available via `coverage.polymerase_extension_reach(..., coverage_metric='processivity')`.

## Reviewing results

```bash
neoswga interpret -d ./
neoswga report    -d ./ --interactive
```

## Alternative: fixed set size

Instead of `--auto-size`, you can fix the set size directly. The committed
`params.json` sets `target_set_size: 20`, which will be honored by optimizers
that respect it (`genetic`, `moea`). Dominating-set uses it as an upper bound
only and terminates early when marginal coverage is zero.
