# equiphi29 scenario

A template configuration for selective whole-genome amplification with
EquiPhi29 polymerase at 42-45 C using longer (10-18 bp) primers plus
5% DMSO and 1 M betaine.

## When to use this template

- Target is GC-rich or has strong secondary structure where higher
  amplification temperature helps.
- You want longer primers for improved specificity, which EquiPhi29
  tolerates better than phi29.
- Betaine equalises AT / GC binding kinetics; DMSO reduces secondary
  structure.

## Running

Replace `target.fasta` / `host.fasta` with your real files (or remove
`host.fasta` entirely if you have no off-target genome), then:

```bash
neoswga count-kmers -j params.json
neoswga filter -j params.json
neoswga score -j params.json
neoswga optimize -j params.json
neoswga report -d .
```

## Tuning pointers

- `reaction_temp`: anywhere in 42-45 C is supported.
- Drop `dmso_percent` and `betaine_m` to 0 if your target is not
  structurally difficult.
- Bump `max_k` down to 15 if amplification specificity is already
  sufficient; shorter primers are cheaper in oligo synthesis.
- For clinical applications, add `optimization_method = "background-aware"`
  to minimise host-background amplification further.
