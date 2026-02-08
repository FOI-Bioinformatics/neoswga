# Plasmid Example: pcDNA vs pLTR

A minimal working example that designs SWGA primers to selectively amplify
**pcDNA** (target, 6157 bp) over **pLTR** (background, 6258 bp).

These are small plasmids, so the pipeline completes in under a minute and is
suitable for verifying a new installation.

## Running the example

```bash
cd examples/plasmid_example
neoswga count-kmers -j params.json   # Step 1 - count 6-12mer frequencies
neoswga filter -j params.json        # Step 2 - apply sequence filters
neoswga score -j params.json         # Step 3 - predict amplification efficacy
neoswga optimize -j params.json      # Step 4 - select optimal primer set
```

## Expected output

| Step | Approximate primer count | Key output file |
|------|--------------------------|-----------------|
| 1 | all 6-12mers | `*_Xmer_all.txt` |
| 2 | ~100-300 | `step2_df.csv` |
| 3 | ~20-80 | `step3_df.csv` |
| 4 | 6 (default set size) | `step4_improved_df.csv` |

Exact counts vary slightly between runs due to stochastic elements in the
scoring and optimization steps.

## Interpreting results

After the pipeline completes, the final primer set is in `step4_improved_df.csv`.
Each row is a primer with columns for sequence, coverage, enrichment score, and
Gini index. Lower Gini values indicate more uniform binding across the target.

For a more detailed quality assessment:

```bash
neoswga interpret -d .
```
