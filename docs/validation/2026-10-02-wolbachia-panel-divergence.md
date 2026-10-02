# The Wolbachia panel resolves three divergence levels, and the ordering agrees with the published supergroups

2026-10-02. Phase 1 task 4 of
`docs/superpowers/plans/2026-10-02-genomic-diversity-and-multi-host.md`: record
pairwise divergence per target pair so later per-strain results can be read
against it. Measured on the five target genomes of
`tests/validation/genomes/diversity_panel.json` -- wMel, wRi, wPip, wAlbB and
wBm -- counted at k=12 with jellyfish 2.3.1. KMC3 is the project default and is
not installed on this machine; both counters emit canonical k-mers with the same
content, so the tables are interchangeable here. Counting the four new
*Wolbachia* genomes plus the *Lactobacillus* host took 2.2 s of step time
(2.8 s wall clock) and wMel a further 1.0 s; the measurement itself runs in
2.5 s, of which 1.25 s is reading the five tables, at 599 MiB peak RSS.

Every figure below is a k-mer measure. **It is not ANI, and it is not a
phylogenetic distance.** A single substitution removes up to k shared k-mers, an
indel or a rearrangement removes a different number again, and nothing here
aligns anything.

## The panel

| key | strain | accession | length | records | distinct 12-mers | per Mb |
|---|---|---|---|---|---|---|
| `wolbachia` | wMel (reference) | GCF_000008025.1 | 1,267,782 | 1 | 921,553 | 726,902 |
| `wolbachia_wri` | wRi | GCF_000022285.1 | 1,445,873 | 1 | 944,512 | 653,247 |
| `wolbachia_wpip` | wPip Pel | GCF_000073005.1 | 1,482,455 | 1 | 993,400 | 670,105 |
| `wolbachia_walbb` | wAlbB | GCF_004171285.1 | 1,484,007 | 1 | 912,773 | 615,073 |
| `wolbachia_wbm` | wBm (TRS) | GCF_000008385.1 | 1,080,084 | 1 | 849,923 | 786,905 |

The canonical 12-mer space holds 8,390,656 k-mers, so each genome occupies
between 0.1013 (wBm) and 0.1184 (wPip) of it.

## Pairwise shared canonical 12-mers

| pair | shared | union | Jaccard | A in B | B in A |
|---|---|---|---|---|---|
| wMel vs wRi | 745,422 | 1,120,643 | **0.6652** | 0.8089 | 0.7892 |
| wAlbB vs wPip | 667,905 | 1,238,268 | **0.5394** | 0.7317 | 0.6723 |
| wMel vs wPip | 391,624 | 1,523,329 | 0.2571 | 0.4250 | 0.3942 |
| wPip vs wRi | 395,976 | 1,541,936 | 0.2568 | 0.3986 | 0.4192 |
| wAlbB vs wRi | 377,251 | 1,480,034 | 0.2549 | 0.4133 | 0.3994 |
| wMel vs wAlbB | 370,325 | 1,464,001 | 0.2530 | 0.4018 | 0.4057 |
| wBm vs wRi | 314,720 | 1,479,715 | 0.2127 | 0.3703 | 0.3332 |
| wMel vs wBm | 309,175 | 1,462,301 | 0.2114 | 0.3355 | 0.3638 |
| wBm vs wPip | 318,594 | 1,524,729 | 0.2090 | 0.3749 | 0.3207 |
| wAlbB vs wBm | 300,667 | 1,462,029 | 0.2057 | 0.3294 | 0.3538 |

Jaccard is the symmetric figure and is the one to quote. The two one-sided
fractions are reported because the genomes differ in size by up to 404 kb, so a
"shared fraction" with no direction stated is ambiguous.

## Three levels, and the ordering agrees with the supergroups

**The supergroup assignment is from the literature, not from this
measurement.** A: wMel and wRi. B: wPip and wAlbB. D: wBm, from the nematode
*Brugia malayi*. The assignments are recorded in the catalogue entries in
`scripts/fetch_reference_genomes.py` and were used to choose the panel; they
were not derived here. They were checked on 2026-10-02 against Klasson et al.
2009 (the wRi genome, PNAS 106:5725, which places wMel and wRi in A and wPip
in B) and Sinha et al. 2019 (the wAlbB genome, Genome Biol Evol 11:706, which
places wAlbB in B); the wBm supergroup D assignment appears in both as the
nematode outgroup.

| level | pairs | Jaccard range | spread within the level |
|---|---|---|---|
| within supergroup A | 1 | 0.6652 | - |
| within supergroup B | 1 | 0.5394 | - |
| A against B | 4 | 0.2530 - 0.2571 | 0.0041 |
| D against A or B | 4 | 0.2057 - 0.2127 | 0.0070 |

The four bands do not overlap. Both within-supergroup pairs sit above every
cross-supergroup pair, and every pair involving wBm sits below every A-against-B
pair. The gap between the A-against-B band and the wBm band is 0.0403, which is
5.8 times the larger of the two within-band spreads; the gap between within-B
and A-against-B is 0.2823.

So the measured ordering agrees with the published grouping, on all ten pairs.
It is worth stating what that agreement is worth: the panel was chosen to span
supergroups, so agreement here confirms that an exact-12-mer count sees the
divergence that was intended, and it establishes a usable gradient. It is not
independent evidence for the phylogeny.

## The chance floor is not zero, and it compresses the distant end

Two unrelated genomes that each occupy a ninth of the canonical 12-mer space
share a measurable fraction of it by coincidence. The expected figures assume
the two sets are drawn independently from that space.

| pair | Jaccard | Jaccard if independent | ratio |
|---|---|---|---|
| wMel vs wRi | 0.6652 | 0.0589 | 11.30 |
| wAlbB vs wPip | 0.5394 | 0.0601 | 8.97 |
| wMel vs wAlbB | 0.2530 | 0.0578 | 4.38 |
| wAlbB vs wRi | 0.2549 | 0.0586 | 4.35 |
| wMel vs wPip | 0.2571 | 0.0604 | 4.26 |
| wPip vs wRi | 0.2568 | 0.0612 | 4.19 |
| wMel vs wBm | 0.2114 | 0.0556 | 3.80 |
| wBm vs wRi | 0.2127 | 0.0563 | 3.78 |
| wAlbB vs wBm | 0.2057 | 0.0554 | 3.72 |
| wBm vs wPip | 0.2090 | 0.0577 | 3.62 |

No pair is at the floor; the least similar pair is still 3.6 times above it.
But the eight cross-supergroup pairs span 0.2057 to 0.2571 against a floor of
about 0.058, so the whole distant end of the panel occupies 0.05 of Jaccard
while the near end spans 0.28. **The measure discriminates poorly at distance
at this k.** A finer separation between the A-against-B and the wBm levels would
need a longer k, not a different index. No longer k has been counted for these
genomes.

## The asymmetry is real and tracks set size

| pair | A in B | B in A | difference |
|---|---|---|---|
| wAlbB vs wPip | 0.7317 | 0.6723 | 0.0594 |
| wBm vs wPip | 0.3749 | 0.3207 | 0.0541 |
| wBm vs wRi | 0.3703 | 0.3332 | 0.0371 |
| wMel vs wPip | 0.4250 | 0.3942 | 0.0307 |
| wMel vs wBm | 0.3355 | 0.3638 | 0.0283 |
| wAlbB vs wBm | 0.3294 | 0.3538 | 0.0244 |
| wPip vs wRi | 0.3986 | 0.4192 | 0.0206 |
| wMel vs wRi | 0.8089 | 0.7892 | 0.0197 |
| wAlbB vs wRi | 0.4133 | 0.3994 | 0.0139 |
| wMel vs wAlbB | 0.4018 | 0.4057 | 0.0039 |

In every one of the ten pairs the higher fraction belongs to the genome with
fewer distinct 12-mers, which is what a containment measure does. The largest
difference, 0.0594 on wAlbB against wPip, is larger than the entire spread of
the A-against-B level, so a one-sided figure quoted without its direction could
reorder the panel.

## Distinct 12-mers do not follow genome length

wAlbB is the longest genome in the panel at 1,484,007 bp and holds the fewest
distinct 12-mers per Mb, 615,073, while wBm is the shortest at 1,080,084 bp and
holds the most, 786,905. wAlbB holds fewer distinct 12-mers in absolute terms
(912,773) than wMel does (921,553) in 216 kb less sequence.

This was not expected and matters beyond this record: the candidate shortlist is
drawn from distinct k-mers, so the genome offering the largest pool is not the
longest one. Repeat content is the obvious explanation and is not measured here.

## What this does not establish

- **Not ANI and not a phylogenetic distance.** Nothing was aligned. The figures
  order the panel; they do not quantify identity, and they must not be reported
  as a percentage identity.
- **The supergroup assignment is not a result.** It comes from the literature,
  cited above. The agreement reported here is agreement with an assumption the
  panel was built on, not a test of it. Nothing here could overturn a
  supergroup assignment, and a panel chosen to span supergroups was always
  going to separate into bands; what the measurement adds is the size of the
  separation and the fact that an exact-12-mer count sees it.
- **One k and one measure.** Everything is at k=12 with exact canonical
  matching. No other k was counted, so the claim that a longer k would separate
  the distant pairs better is an expectation and not a measurement.
- **Nothing about primer coverage.** A shared 12-mer is not a usable primer: no
  GC window, Tm, dimer or occupancy criterion was applied, and the SWGA
  candidate pool is a small filtered subset. How far per-strain coverage falls
  across this panel is Phase 3's question and is open.
- **Nothing about the host panel.** *Lactobacillus* was counted at k=12 but is
  a host, so it is not in the pairwise measurement. No target-against-host
  figure is reported here.
- **Nothing at SNP resolution.** The measure cannot distinguish a substitution
  from an indel or a rearrangement, so it cannot say whether a variant table is
  a reasonable description of any pair. That is Phase 5's gate B.
- **No claim about mismatch tolerance.** This is exact matching throughout,
  consistent with the rest of the package.

## Reproducing

```bash
# Tables, if tests/validation/genomes/ does not already hold them. A params.json
# naming the five genomes as fg_genomes with fg_prefixes alongside them, then:
neoswga count-kmers -j params.json       # min_k 12, max_k 12

# The measurement.
python scripts/benchmarking/wolbachia_panel_divergence.py \
    --manifest tests/validation/genomes/diversity_panel.json \
    --k 12 \
    --output divergence_k12.json
```

The script is `scripts/benchmarking/wolbachia_panel_divergence.py`. It reads
the panel manifest, takes the entries with `role: target`, and reads counts only
through `kmer_tables.iter_table`, which raises on a prefix that was never
counted rather than reporting an empty set. Output ordering is by key and by
sorted key pair, so a re-run is comparable line by line. The genome files are
pinned by SHA-256 in the manifest and in
`scripts/fetch_reference_genomes.py`.

The exact sets are held in memory, which is affordable only because a
*Wolbachia* genome has about a million distinct 12-mers. This does not transfer
to a host-sized genome; the script's docstring says so.
