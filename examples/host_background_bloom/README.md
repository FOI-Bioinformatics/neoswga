# Gb-scale host background with Bloom filter

Template for SWGA against a pathogen target where the host background
(human, mouse, plant) is too large to load as raw k-mer counts. The
Bloom filter stores a probabilistic representation of host k-mers,
trading a small false-positive rate for a ~100x memory reduction.

## Workflow

1. Build the Bloom filter once from the host FASTA:

   ```bash
   neoswga build-filter path/to/GRCh38.fna ./bloom
   ```

   This writes `./bloom/background.pkl` and `./bloom/sampled.pkl`.

2. Point `params.json` at the pre-built filter (this template already
   does):

   ```json
   {
     "use_bloom_filter": true,
     "bloom_filter_path": "./bloom/background.pkl",
     "sampled_index_path": "./bloom/sampled.pkl"
   }
   ```

3. Run the pipeline as usual:

   ```bash
   neoswga count-kmers -j params.json   # only counts target k-mers
   neoswga filter -j params.json        # queries the Bloom filter for bg freq
   neoswga score -j params.json
   neoswga optimize -j params.json
   ```

## When to use this template

- Host backgrounds >= 100 Mb where the all-k-mer text file would exceed
  available RAM.
- Multiple pipeline runs against the same host (the filter is built once
  and reused across targets / conditions).
- Clinical or metagenomic scenarios where the host genome is fixed.

## Limitations

- The Bloom filter has a ~0.1% false-positive rate by default. Primers
  that "just barely" pass the `max_bg_freq` threshold may occasionally
  be rejected or admitted incorrectly. Tighten `max_bg_freq` if this
  matters for your use case.
- The `sampled_index_path` is used to estimate `bg_freq` values for
  reporting; it is not strictly required but improves the report's
  accuracy.

## Scaling pointer

For host genomes above 10 Gb, also consider:
- `min_k: 12` or higher for better specificity on long genomes.
- Running the filter on a machine with enough RAM to hold the Bloom
  bitmap (default 2 GB for human, see `neoswga build-filter --help`).
