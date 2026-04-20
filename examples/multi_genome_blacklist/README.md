# Multi-genome + blacklist scenario

A template with two target genomes (pan-primer design), one background,
and one blacklist. Primers must bind at least one target frequently,
avoid the background to within `max_bg_freq`, and have zero tolerance
for binding the blacklist (`max_bl_freq: 0.0`).

## Use cases

- Designing primers that must amplify multiple strains simultaneously
  (pan-pathogen / pan-serotype)
- Excluding contaminant sequences that could produce false positives
  (e.g., nearby commensal species)

## Running

Replace the FASTA paths with your real files:

```bash
neoswga count-kmers -j params.json
neoswga filter -j params.json
neoswga score -j params.json
neoswga optimize -j params.json
```

## Blacklist vs background vs exclusion

| Category   | Parameter          | Behaviour                                 |
|-----------|--------------------|-------------------------------------------|
| Background | `bg_genomes`       | Tolerated below `max_bg_freq`.           |
| Blacklist  | `bl_genomes`       | Weighted penalty, rejected if `freq > max_bl_freq`. |
| Exclusion  | `excl_genomes`     | Hard-reject any hit (`excl_threshold=0`). |

Use blacklist for "avoid unless unavoidable" and exclusion for strict
zero-tolerance contaminants such as mtDNA or chloroplast when you know
you want no amplification of those sequences at all.

## Scaling

For Gb-scale host backgrounds (e.g. human genome), pre-build a Bloom
filter once and reuse it:

```bash
neoswga build-filter GRCh38.fna ./bloom
```

Then set `use_bloom_filter: true` and `bloom_filter_path: ./bloom/background.pkl`
in params.json.
