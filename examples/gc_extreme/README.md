# Extreme GC scenario

Template for genomes with strongly skewed GC content (e.g. *Plasmodium
falciparum* ~20% GC, *Mycobacterium tuberculosis* ~65% GC, *Burkholderia*
~67% GC, *Streptomyces* ~72% GC).

## Why this scenario exists

The default primer GC window (37.5-62.5%) excludes most candidate primers
from extreme-GC targets. When `genome_gc` is known or auto-computed, the
pipeline switches to `genome_gc +/- gc_tolerance`, which for a 25% GC
genome yields a 10-40% window instead.

Additives help further:
- Betaine 1-2 M equalises AT / GC binding.
- TMAC tightens the Tm distribution for GC extremes.
- DMSO reduces secondary structure in GC-rich targets.

## Running

```bash
neoswga count-kmers -j params.json
neoswga filter -j params.json    # adaptive GC engages automatically
neoswga score -j params.json
neoswga optimize -j params.json
```

To verify the adaptive path engaged, look for this line in the pipeline
output:

```
Extreme GC genome detected (GC=...%). Adaptive GC filter engaged: primer GC window 0.XX-0.XX.
```

## Pinning your own GC window

If you want the old fixed window regardless of genome GC, set:

```json
"adaptive_gc": false,
"gc_min": 0.375,
"gc_max": 0.625
```

or explicitly pin `gc_min` / `gc_max` to your chosen values.
