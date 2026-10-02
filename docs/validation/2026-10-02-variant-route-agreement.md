# A variant table against one reference does not substitute for a strain genome across supergroups, and it fails in the optimistic direction

2026-10-02. Phase 5 of
`docs/superpowers/plans/2026-10-02-genomic-diversity-and-multi-host.md`:
decision gate B. Measured on the five *Wolbachia* genomes of
`tests/validation/genomes/diversity_panel.json` with wMel (`wolbachia`,
1,267,782 bp, one record, treated as circular) as the reference, and on the
twelve panels Phase 3 delivered in its D1 and D2 groups (6, 12 and 24 requested
primers, all 12-mers, each delivered at full size). Variant tables were derived
here with minimap2 2.28-r1209 and `paftools.js call` 2.28-r1209; neither is a
dependency of this package and neither is on PATH, so both were taken as script
arguments from a conda environment. Nothing was installed. Per-primer reach
3,000 bp, as Phase 3 ran it. The whole measurement takes about 25 s and peaks
well under 1 GB.

**The headline.** The variant route reports 77 to 90 percent of the panel's
reference binding sites intact in every strain, including the most divergent
one, while the alignment the table was derived from says only 0.7 to 54 percent
of those sites exist at their homologous locus in that strain. The route's
coverage figure exceeds the like-for-like figure by 0.17 to 0.26 for wRi (within
supergroup A, k-mer Jaccard 0.6652 against wMel) and by 0.41 to 0.75 for the
three cross-supergroup strains. In all twelve panels the variant route places
wBm -- the most divergent strain in the panel -- second best of the four; the
genome route places it worst in ten of twelve and the alignment route in eleven
of twelve. The two routes name the same worst strain in 2 of 12 panels. The
error is not noise: it is one-directional, and it is optimistic.

The cause is measured, not inferred. Essentially the whole disagreement is
reference sequence the aligner never aligned: a site there carries no variant,
so the route calls it intact. Of the sites where the variant route says intact
and the alignment says the site does not survive, the number that lie outside
the region the variant caller was allowed to work in is 1,827 of 1,827 for wRi,
4,327 of 4,342 for wPip, 4,258 of 4,258 for wAlbB and 4,614 of 4,614 for wBm,
summed over the twelve panels -- so 15 such sites in 5,441 x 4 verdicts are not
explained by the blind region. Inside that region the two routes agree on 98.7
to 100 percent of sites.

## Three routes, three denominators

Mixing them up is the easiest way to misread every table below.

| route | what it measures | denominator | code |
|---|---|---|---|
| variant | reference sites with no variant inside them | the REFERENCE, wMel, 1,267,782 bp | `variant_table.open_variants` + `variant_sites.evaluate_strain_panel`, as `evaluate-set --variants` drives them |
| alignment | reference sites whose k bases are all aligned to an identical base with no insertion inside | the REFERENCE, wMel | derived in the script from minimap2's `cs` tag, independently of the variant caller |
| genome | every site found in the strain's own sequence, gained sites included | the STRAIN's own genome | `reference_panel_evaluation.evaluate_reference_panel` over the strain FASTA |

The **variant and alignment routes are the like-for-like pair**: the same sites,
the same reference denominator, the same coverage function
(`coverage.compute_per_prefix_coverage` over a `variant_sites.IntactPositions`
view). Their difference is the main result and the only difference reported as a
subtraction.

The **genome route has a different denominator** and its figures are reported
beside, never subtracted from, the other two. It is the route Phase 3 used and
it is what a user would actually run; its disagreement with the variant route is
shown as a ranking comparison rather than a difference.

The genome-route figures quoted here were recomputed through
`evaluate_reference_panel` for two of the twelve panels and reproduced Phase 3's
`results.json` exactly -- all 10 site counts identical, all 10 coverage figures
agreeing to 0.000000 -- so the figures taken from that file for the remaining
ten panels are being compared like with like.

## What aligns, and under which preset

This is a primary result, not a diagnostic: reference sequence with no alignment
is exactly what the variant route cannot see.

minimap2's assembly presets are tied to divergence -- asm5 to about 0.1 percent,
asm10 to about 1 percent, asm20 to about 5 percent. The panel spans
within-supergroup to cross-supergroup, so no single preset was assumed to suit
every pair; all three were run and all three are reported. "callable" is the
union of blocks that pass the variant caller's own filters, below.

| strain | Jaccard vs wMel | preset | primary blocks | median block | longest block | reference aligned | callable |
|---|---|---|---|---|---|---|---|
| wRi | 0.6652 | asm5 | 224 | 2,384 | 57,962 | 0.6329 | 0.6114 |
| wRi | 0.6652 | asm10 | 115 | 6,059 | 70,249 | 0.6743 | 0.6657 |
| wRi | 0.6652 | asm20 | 80 | 9,052 | 97,314 | **0.6983** | 0.6932 |
| wPip | 0.2571 | asm5 | 46 | 766 | 4,781 | 0.0319 | 0.0186 |
| wPip | 0.2571 | asm10 | 119 | 1,420 | 6,152 | 0.1594 | 0.1485 |
| wPip | 0.2571 | asm20 | 289 | 1,848 | 22,146 | **0.5744** | 0.5561 |
| wAlbB | 0.2530 | asm5 | 34 | 968 | 5,183 | 0.0346 | 0.0273 |
| wAlbB | 0.2530 | asm10 | 91 | 1,547 | 8,786 | 0.1523 | 0.1429 |
| wAlbB | 0.2530 | asm20 | 247 | 1,880 | 22,141 | **0.5715** | 0.5543 |
| wBm | 0.2114 | asm5 | 2 | 2,258 | 2,910 | 0.0035 | 0.0035 |
| wBm | 0.2114 | asm10 | 6 | 2,084 | 3,111 | 0.0101 | 0.0094 |
| wBm | 0.2114 | asm20 | 237 | 1,872 | 22,448 | **0.5013** | 0.4843 |

Jaccard over canonical 12-mers is quoted from
`docs/validation/2026-10-02-wolbachia-panel-divergence.md`; it is a k-mer
measure, not ANI.

**The preset is not a detail.** At asm5 the three cross-supergroup strains align
3 percent of the reference or less, and wBm aligns 0.35 percent in two blocks. A
variant table built that way would describe the aligned 0.35 percent and say
nothing at all about the rest, while the route would still report almost every
site intact. asm20 is the preset used for the variant tables below, chosen
because it aligns the most of the reference for every pair. Even at asm20 the
reference is only half aligned for the three distant strains. No preset was
found that aligns most of a cross-supergroup *Wolbachia* pair, and minimap2
ships no assembly preset more permissive than asm20.

## The variant tables

`paftools.js call` was given `-l 1000 -L 1000 -q 5`. Its defaults are `-l 10000
-L 50000`, written for large eukaryotic assemblies; at `-L 50000` almost no
block in a cross-supergroup pair is eligible and almost no variant is called,
which would read as "these genomes are nearly identical". The callable region in
every table here is computed under the same two filters the caller was given, so
the region the route can see is the region the caller actually used.

| strain | Jaccard | reference aligned | reference identical | callable | VCF records | SNPs | insertions | deletions |
|---|---|---|---|---|---|---|---|---|
| wRi | 0.6652 | 0.6983 | 0.6712 | 0.6932 | 18,751 | 17,718 | 539 | 494 |
| wPip | 0.2571 | 0.5744 | 0.4980 | 0.5561 | 80,907 | 79,163 | 674 | 1,070 |
| wAlbB | 0.2530 | 0.5715 | 0.4958 | 0.5543 | 81,116 | 79,417 | 637 | 1,062 |
| wBm | 0.2114 | 0.5013 | 0.4278 | 0.4843 | 85,153 | 83,672 | 764 | 717 |

`open_variants` accepted all four VCFs as written by `paftools.js call`, with no
conversion and no edit. **Nothing was refused**: no REF mismatch, no symbolic
allele, no unsorted record, no unknown contig. Each table resolved to one
measured strain, and each carried the note the door writes about multi-base REF
alleles. `evaluate-set --variants` was also run once end to end on the wRi table
and the D1 n=6 panel through the CLI, and reported intact 336, affected 77,
intact fraction 81.4 percent, coverage 61.5 percent -- the same figures the
script computes through the underlying functions.

## The main table: the like-for-like pair, D1 n=6

Both columns are the same 413 reference sites of the same panel over the same
1,267,782 bp denominator. The genome-route columns are shown beside them with
their own, different denominator.

`D1__wolbachia__vs__lactobacillus`, 6 primers, 413 reference sites, reference
coverage 0.6675:

| strain | Jaccard | intact, variant | surviving, alignment | cov. of reference, variant | cov. of reference, alignment | difference | sites in strain, genome | cov. of strain, genome |
|---|---|---|---|---|---|---|---|---|
| wRi | 0.6652 | 336 | 196 | 0.6153 | 0.4131 | **+0.2022** | 398 | 0.5971 |
| wPip | 0.2571 | 346 | 7 | 0.5749 | 0.0331 | **+0.5418** | 120 | 0.2728 |
| wAlbB | 0.2530 | 342 | 10 | 0.5743 | 0.0473 | **+0.5269** | 130 | 0.2875 |
| wBm | 0.2114 | 369 | 8 | 0.6129 | 0.0358 | **+0.5771** | 55 | 0.2494 |

The difference is positive in every row of every panel: the variant route is
never pessimistic in aggregate.

### It scales the same way across all twelve panels

| strain | Jaccard | panels | difference, min | median | max | panels above the tolerance |
|---|---|---|---|---|---|---|
| wRi | 0.6652 | 12 | 0.1660 | 0.2075 | 0.2623 | 12 of 12 |
| wPip | 0.2571 | 12 | 0.4257 | 0.5380 | 0.6732 | 12 of 12 |
| wAlbB | 0.2530 | 12 | 0.4071 | 0.5230 | 0.6423 | 12 of 12 |
| wBm | 0.2114 | 12 | 0.4556 | 0.5836 | 0.7457 | 12 of 12 |

The difference grows with panel size, because a larger panel has more sites and
more of them fall outside the aligned core. It does not shrink: the smallest
difference anywhere in the measurement is 0.1660, on the closest pair at the
smallest panel.

## The central weakness, with a number on it

A reference site in sequence the aligner did not align carries no variant,
because none was looked for there, and the route reports it intact. Counted per
strain over the twelve panels, which hold 5,441 reference sites in total:

| strain | reference sites outside the callable region | share of the panel's sites | of those, surviving by alignment | of those, called intact by the variant route |
|---|---|---|---|---|
| wRi | 1,913 | 0.335 - 0.365 | 86 | **1,913** |
| wPip | 4,329 | 0.707 - 0.848 | 2 | **4,329** |
| wAlbB | 4,260 | 0.695 - 0.842 | 2 | **4,260** |
| wBm | 4,620 | 0.766 - 0.886 | 6 | **4,620** |

Every one of them is reported intact. For the three cross-supergroup strains
that is 70 to 89 percent of the panel's sites being reported on without
evidence.

Inside the callable region the two routes very nearly agree, which is what
isolates the blind region as the cause:

| strain | inside-callable sites, 12 panels | variant route intact where the alignment says no | alignment keeps it, variant route calls it affected |
|---|---|---|---|
| wRi | 3,528 | 0 | 36 (1.0%) |
| wPip | 1,112 | 15 (1.3%) | 13 (1.2%) |
| wAlbB | 1,181 | 0 | 12 (1.0%) |
| wBm | 821 | 0 | 0 |

The residual is small and runs both ways: 15 optimistic sites in 6,642
inside-callable verdicts, and 61 pessimistic ones. The 15 optimistic sites for
wPip and the 0 to 36 pessimistic sites per strain were not traced to a cause;
block edges
and overlapping primary alignments are the obvious candidates and neither was
checked.

### The panel's sites are blinder than the genome is

The share of the reference that is callable is not the share of the panel's
sites that are callable. The control below applies the same two verdicts to
**every** reference 12-mer position rather than to the panel's sites:

| strain | all 12-mer windows inside the callable region | all 12-mer windows surviving by alignment | the panel's sites inside the callable region (D1 n=6) | the panel's sites surviving (D1 n=6) |
|---|---|---|---|---|
| wRi | 0.6926 | 0.5766 | 0.637 | 0.475 |
| wPip | 0.5542 | 0.1539 | 0.177 | 0.017 |
| wAlbB | 0.5525 | 0.1546 | 0.194 | 0.024 |
| wBm | 0.4825 | 0.0918 | 0.126 | 0.019 |

For wRi the panel's sites behave roughly like a random window. For the three
distant strains they do not: a random reference 12-mer window is inside the
callable region 48 to 55 percent of the time, while the panel's sites are 13 to
19 percent of the time. The panel's binding sites are concentrated in reference
sequence that does not align across supergroups. The candidate pool is drawn on
foreground frequency, and the delivered 12-mers at the low-Tm end of these
panels are strongly AT-rich (`TAGTAGAAGAAA`, `AAAGAAGCAAAA`), so composition is
the obvious explanation; **it was not measured here** and is stated as a
candidate explanation only.

The consequence is that reasoning from "the aligner aligned half the reference,
so the route sees half the sites" understates the blindness by a factor of three
for these panels.

## The ranking inversion

Coverage figures are not the only thing a user reads; which strain is worst is.

| design | n | worst by variant route | worst by alignment | worst by genome route | wBm's rank, variant | alignment | genome |
|---|---|---|---|---|---|---|---|
| D1 vs lactobacillus | 6 | wAlbB | wPip | wBm | 2 | 3 | 4 |
| D1 vs lactobacillus | 12 | wAlbB | wBm | wAlbB | 2 | 4 | 2 |
| D1 vs lactobacillus | 24 | wAlbB | wBm | wAlbB | 2 | 4 | 2 |
| D1 vs drosophila | 6 | wAlbB | wBm | wBm | 2 | 4 | 4 |
| D1 vs drosophila | 12 | wPip | wBm | wBm | 2 | 4 | 4 |
| D1 vs drosophila | 24 | wPip | wBm | wBm | 2 | 4 | 4 |
| D2 vs lactobacillus | 6 | wAlbB | wBm | wBm | 2 | 4 | 4 |
| D2 vs lactobacillus | 12 | wPip | wBm | wBm | 2 | 4 | 4 |
| D2 vs lactobacillus | 24 | wAlbB | wBm | wBm | 2 | 4 | 4 |
| D2 vs drosophila | 6 | wPip | wBm | wBm | 2 | 4 | 4 |
| D2 vs drosophila | 12 | wPip | wBm | wBm | 2 | 4 | 4 |
| D2 vs drosophila | 24 | wAlbB | wBm | wBm | 2 | 4 | 4 |

The variant route never names wBm as the worst strain and ranks it second of
four in all twelve panels. The genome route names wBm worst in ten of twelve,
the alignment route in eleven of twelve. Agreement between the variant and
genome routes on the identity of the worst strain: 2 of 12.

`evaluate-set --variants` reports `worst_strain_coverage` and
`worst_strain_intact_fraction` as reductions. On this panel those two reductions
pick the wrong strain most of the time, and they are optimistic when they do.

## Gate B, read against the figures

The tolerance is **0.05 absolute difference in coverage fraction on the
reference denominator**. It is a stated choice, not a derived one, and no code
reads it: Phase 3 measured per-strain coverage falling by 0.25 to 0.40 between
wMel and a cross-supergroup strain, so 0.05 is about an eighth of the effect the
two routes are being used to measure. Below it the routes would rank the strains
alike and a reader would draw the same conclusion from either; above it they
would not. The raw differences are in the results file, so another tolerance can
be applied without re-running anything.

**Divergent regime.** The plan expected the variant route to diverge from the
genome route here, and it does, on all twelve panels and all four strains. The
difference exceeds the tolerance at **every** divergence level this panel
offers, including the closest pair: 0.1660 to 0.2623 at Jaccard 0.6652 (within
supergroup A), and 0.4071 to 0.7457 at Jaccard 0.2114 to 0.2571 (across
supergroups). So the level at which the difference exceeds the tolerance is
**below the closest pair on disk**; this panel cannot locate it, and it is not
located here. Above Jaccard 0.6652 nothing was measured.

**Clonal regime: not measured.** The plan expected close agreement in a clonal
regime. **There is no clonal pair in this panel.** The closest pair is wMel
against wRi, two strains within supergroup A at k-mer Jaccard 0.6652, 30 percent
of the reference unaligned at the best preset, and 18,751 called variants over
1.27 Mb. That is not a clonal relationship, and the clonal half of gate B is
therefore **not measured** -- not "measured and found adequate", and not
"measured and found inadequate". What the measurement bounds is one end only:
agreement is already outside the tolerance at this divergence, so the regime in
which the route is adequate, if it exists, lies strictly closer than wMel
against wRi. Testing it needs two genomes that differ by SNPs over nearly all
their length; none is on disk, and Phase 4's oracle test (a reference and a
seeded mutant of it) is the existing in-repo case for that end, which the route
passes exactly by construction.

**What follows for Phase 6.** Phase 6b -- accepting a variant table at design
time -- is **not supported by this measurement** for any pair this panel
contains. On these twelve panels and these five genomes a variant table derived
by this aligner and caller reports coverage too high by 0.17 to 0.75 and names
the wrong worst strain in ten of twelve panels. Strain genomes are the supported
input across supergroups, and the plan's statement that the divergence level
should be documented resolves to: this panel cannot bracket it, and no level at
which the route is adequate was demonstrated.

Three things would change the reading and none was done: a clonal pair, a
different aligner or preset (the route inherits the aligner's blind spots
entirely), and an aligned-region-only variant route that refused to report on
sequence it could not see. The third is the one this measurement argues for: the
failure is not that the variant model is wrong -- inside the callable region it
agrees to within about 1 percent of sites -- but that it reports on regions
where it has no evidence, which is the "unknown is not zero, and not success"
rule in a new form. A route that carried the callable region alongside the
variants and reported sites outside it as not assessed rather than intact would
have been honest on all twelve panels. `variant_sites.StrainSites` already
carries a `not_assessed_sites` field for a related case, so the shape exists.
That is a design proposal, not a measurement, and it is for Phase 6 to accept or
reject.

## What this does not establish

- **Five genomes of one genus, one reference, twelve panels, one panel size
  family.** Everything here is *Wolbachia* against wMel at k=12 with 12-mer
  oligos. Nothing about another genus, another k, mixed oligo lengths, or a
  reference other than wMel.
- **One aligner, one caller, one preset for the tables.** minimap2 2.28 with
  `-x asm20 --cs=long` and `paftools.js call -l 1000 -L 1000 -q 5`. The three
  presets measured show the aligned fraction moves by a factor of 140 for wBm
  between asm5 and asm20, so the preset is a first-order input and a different
  aligner (`nucmer`/`show-snps`, or a whole-genome aligner tolerant of
  rearrangement) was not tried. **A result measured under one aligner says
  nothing about another.**
- **The alignment route is independent of the caller, not of the aligner.** It
  is derived here from the same minimap2 `cs` tag the variant table was called
  from. It is therefore a check on the caller and on the route's blindness, and
  not an independent ground truth for what a strain contains. The genome route
  is the independent one, and its denominator differs.
- **No clonal pair.** See gate B above. The clonal half of the gate is not
  measured.
- **The residual disagreement inside the callable region was not traced.** 0 to
  1.4 percent of inside-callable sites disagree in each direction; the cause was
  not established.
- **Why the panel's sites avoid aligned sequence was not measured.** AT-rich
  composition is a candidate explanation and no composition figure was computed.
- **Nothing about mismatch tolerance.** Both routes are exact-match throughout,
  as the rest of the package is. An "affected" site is counted as lost and is
  never weighted or called tolerated.
- **Nothing about amplification.** These are binding-site and coverage figures.
  No reaction was simulated and no yield was measured.
- **No host and no specificity figure.** Only the five target genomes appear.
  hg38 was deliberately not touched.
- **The genome route's excess is not "gained sites" exactly.** The results file
  reports, per panel and strain, the genome route's sites minus the alignment
  route's surviving sites (138 to 330 for wRi, 43 to 139 for wBm, 111 to 409
  across all four). That count mixes sites genuinely gained in the strain with
  sites at loci the aligner did not align and with duplicated loci, and it is
  named `sites_not_explained_by_a_surviving_reference_site` for that reason. It
  is not a measurement of gain.

## Reproducing

```bash
# Tool paths are arguments; minimap2 and paftools.js are not dependencies of
# this package and are not on PATH. Any conda environment carrying minimap2
# 2.28 works; paftools.js needs the k8 interpreter beside it.
MM2=/path/to/minimap2
PAF=/path/to/paftools.js

# Verify the cs-tag reader the alignment route rests on, and stop. It checks a
# hand-written cs tag against a hand-computed expectation, checks site_survives
# on five hand-chosen offsets, and then plants one substitution and one 3 bp
# deletion in a 5 kb slice of wMel, aligns the mutant back, and asserts the
# reference bases reported as not identically aligned are exactly those four.
python scripts/benchmarking/variant_route_agreement.py --self-check \
    --minimap2 $MM2 --paftools $PAF

# The measurement. About 25 s.
python scripts/benchmarking/variant_route_agreement.py \
    --minimap2 $MM2 --paftools $PAF \
    --results tests/validation/genomes/diversity_baseline/results.json \
    --workdir tests/validation/genomes/variant_route_agreement

# Every table above, from the finished results file.
python scripts/benchmarking/variant_route_agreement.py \
    --tables tests/validation/genomes/variant_route_agreement/results.json
```

The script is `scripts/benchmarking/variant_route_agreement.py`. It writes
alignments, VCFs and `results.json` into
`tests/validation/genomes/variant_route_agreement/`, which is gitignored. It
reads Phase 3's `results.json` only, never Phase 3's position indexes: the
reference is scanned once from its FASTA through `PositionCache` with
`on_missing='scan'` against a prefix in its own work directory, so no HDF5 file
of another run is opened. The variant route is driven through
`variant_table.open_variants` and `variant_sites.evaluate_strain_panel`, the
genome route through `reference_panel_evaluation.evaluate_reference_panel`, and
both coverage figures on the reference denominator through the one
`coverage.compute_per_prefix_coverage`; the alignment route's site subsets are
presented to it through `variant_sites.IntactPositions`, so all three coverage
numbers come out of one function.

The one conversion that lives in the script is the `cs`-tag reader, whose output
every alignment-route figure rests on. It is checked by `--self-check`, above,
and by a consistency test that runs in the measurement itself: a reference site
that survives at its homologous locus means the strain holds that exact 12-mer,
so the genome route must find at least as many sites of that primer in the
strain as the alignment says survive. Over all 48 panel-strain rows and every
primer in them, **no violation was found**. The script also asserts for every
alignment that walking the `cs` tag ends on the reference coordinate the PAF
records as the alignment's end, and errors rather than warns if it does not.

The genomes are pinned by SHA-256 in
`tests/validation/genomes/diversity_panel.json` and in
`scripts/fetch_reference_genomes.py`.
