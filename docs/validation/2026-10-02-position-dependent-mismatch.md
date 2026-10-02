# Does a position-dependent mismatch model predict the measured outcomes better?

Measured 2026-10-02. Phase 4b, validation items 1-4, of
`docs/superpowers/plans/2026-10-02-genomic-diversity-and-multi-host.md`, which
asks for the model to be scored against the wet-lab outcomes already in the
repository, with the controls of
[heldout_panel_discrimination_2026-09-16.md](heldout_panel_discrimination_2026-09-16.md).

Reproduce with `scripts/benchmarking/mismatch_model_ranking.py`.

**No benefit is demonstrated on these three studies.** Within a study the
position-dependent model ranks the measured outcomes slightly WORSE than the
uniform one (mean area under the ROC curve 0.39 against 0.43, with or without
the 3'-terminal rule). Held out by study, thresholded accuracy is 0.57 for
either variant against the uniform model's 0.48, under a majority baseline of
0.62 that nothing but the published control clears. Twenty-one scorable panels,
eight of them effective, three studies, two of which are scored against a
stand-in host.

**Revised after review.** A first pass reported the model AHEAD on held-out
accuracy (0.67 and 0.76). That advantage was an artefact: the internal-mismatch
table was being applied at the first and last base of the duplex, where it has
no domain, and a clamp then substituted exactly zero penalty whenever the
arithmetic came out net stabilising -- weighting a mismatched host site as well
as a perfect one. With terminal positions taken out of the table's domain
(section 1b) the advantage disappears and the uniform model's own figures are
unchanged. The earlier numbers are superseded, not averaged in.

## What was built

| Term | Where | Switched by | Evidence record |
|---|---|---|---|
| Duplex stability, graded by where the mismatch sits | `core/mismatch_model.py` | `mismatch_model: position-dependent` | `mismatch_duplex_delta_g`, status `estimated` |
| 3'-terminal extension block, near binary | same | `mismatch_model: position-dependent-3prime` | `three_prime_mismatch_extension`, status `assumed` |
| Uniform Tm penalty per mismatch | `core/occupancy.py` | `mismatch_model: uniform` (the default) | `mismatch_penalty`, status `assumed` |

The duplex term is an extrapolation and the record says so: the SantaLucia
(1998) internal-mismatch nearest-neighbour set was measured at 37 C in 1 M Na+
on duplex DNA, and it is applied here near 30 C in a magnesium buffer on
8-12mers, converted to a melting-temperature shift by a first-order
perturbation of `Tm = dH / (dS + R ln C)`. The table is a measurement; this use
of it is not. The 3'-terminal rule rests on a docstring citation to Innis
(1988) for Taq in PCR and no measurement on a strand-displacing polymerase was
found, which is why it is a separate value of the key.

No threshold refuses a design for being far from the table's reference,
because no threshold has a reference. The distance is reported instead: the
panel assessment's evidence line, the run manifest's `effective_conditions`
and one warning per run all state how many degrees separate the reaction from
37 C and that the buffer differs. The model also enters `request_hash`, so two
designs that deliver different panels no longer record one identity.

Nothing about a default run moves. With the key unset or set to `uniform`,
`tests/test_design_options_have_effect.py::test_the_mismatch_model_default_moves_nothing`
compares the delivered metrics field by field.

## 1b. The table's domain is internal positions

`DELTA_G_MISMATCH` is an INTERNAL single-mismatch set. A terminal position has
one flanking stack instead of two, and the applicable parameters there are
dangling-end and terminal-mismatch sets (Bommarito et al. 2000) which this
repository does not carry. The first implementation applied the internal set at
both termini anyway, and additionally dropped that end's initiation correction
(+0.98 / +1.03) because the terminus was no longer Watson-Crick paired.
Dropping a positive correction is stabilising, and a stabilising mismatch
doublet can out-stabilise a weak Watson-Crick stack, so the two could net out
positive -- a mismatch making the duplex look MORE stable than the perfect
match. A clamp at zero then turned that into exactly zero penalty, the most
favourable value available, for a quantity the source does not parameterise.

Measured over 1,400 random primers at three polymerases and seven oligo
lengths, before the fix:

| where the mismatch sits | neighbours | would have been clamped to zero penalty |
|---|---|---|
| at either terminus | 8,400 | **771 (9.2%)** |
| at an interior position | 43,199 | **0** |

So the defect was entirely a terminal-position effect, and at interior
positions the arithmetic never produces a stabilising mismatch at all.

**The fix**: a neighbour with a mismatch at either terminus does not use the
table. It takes the uniform assumed penalty, `- distance * mismatch_penalty`,
which is what the shipped model charges it, and the 3'-terminal extension
factor of the `-3prime` variant still applies on top because that is an
extension effect rather than a duplex one. The clamp stays as a guard, and
`test_the_clamp_never_binds_at_an_internal_position` is what makes "it never
binds" a measured statement rather than a hope.

The domain limit is stated in the module docstring, in the
`mismatch_duplex_delta_g` record's `sequence_domain` and `notes`, in the schema
description, in `docs/params-reference.md`, in the implementation notes and in
"Constants this rests on" below.

One consequence worth stating: with terminal positions on the uniform penalty,
the duplex term penalises a 3'-terminal mismatch LESS than most interior ones
(4.0 C against 10 to 19 C on a 12-mer). That is not a claim about extension;
it is the reason the 3'-terminal extension rule has to be a separate term
rather than something the duplex term already expresses.

## 1. The oracle

One 5-mer duplex with one internal mismatch, against a value computed by hand
from the published table. Primer `5'-GCATG-3'` against the genomic site
`5'-GCGTG-3'`, so the duplex is

```
x (primer,   5'->3'):  G  C  A  T  G
y (template, 3'->5'):  C  G  C  A  C
```

| contribution | key | value |
|---|---|---|
| 5' terminal correction | `G/C` | +0.98 |
| 3' terminal correction | `G/C` | +0.98 |
| doublet 1 | `GC/CG` | -2.24 |
| doublet 2 | `CA/GC` | +0.75 |
| doublet 3 | `AT/CA`, absent; the same stack rotated is `AC/TA` | +0.77 |
| doublet 4 | `TG/AC` | -1.45 |
| total | | **-0.21 kcal/mol** |

The code returns -0.21 kcal/mol. A perfect duplex returns exactly what
`thermodynamics.compute_free_energy_for_two_strings(primer,
complement(primer))` returns, for every primer tested.

**A finding on the way.** `DELTA_G_MISMATCH` holds only the 64 doublets whose
5' position is Watson-Crick paired, so the doublet on the 3' side of any
mismatch is absent as written and the shipped walk charges it the flat 4.0
kcal/mol penalty. A nearest-neighbour doublet is unchanged by rotating the
duplex 180 degrees, and the rotated key is tabulated, so the new walk looks it
up. Without the rotation the mismatch above would cost 3.23 kcal/mol more than
it does. The existing function is unchanged; the new module walks the duplex
itself, which is also what removes the early stop at `penalty * 10` that would
truncate a whole-duplex sum.

## 2. The reduction

With the position-dependent path configured and every mismatch weighted alike,
the site load equals `occupancy.weighted_site_load` to a relative tolerance of
1e-9. The two then differ only in the order the sum is taken: the uniform path
sums the counts in a mismatch class and multiplies the total by one occupancy,
and the per-neighbour path multiplies each group's count by the same occupancy
and sums. **The tolerance is floating-point reassociation and nothing else**,
which is some six orders of magnitude tighter than any difference the model
produces, so a failure at it would mean the two paths disagree about the
counting rather than about the position.

## 3. The real test: ranking the measured outcomes

Three studies carry a published per-panel geometry alongside a binary outcome,
as in the 2026-09-16 record. Unlike that record, the predictors here are
**recomputed from k-mer tables rather than read from the papers**, which makes
this a check on this implementation and not on the metric definitions, and
which is why it needs the genomes.

| study | panels scored | effective | foreground | background |
|---|---|---|---|---|
| Clarke 2017, *M. tuberculosis* | 12 (10 without controls) | 4 | `mtb.fna` | human chr21, a stand-in |
| Clarke 2017, *Wolbachia* | 5 | 1 | wMel | *Drosophila* 144 Mb, as published |
| Dwivedi-Yu 2023, *Prevotella* | 6 | 3 | `prevotella.fna` | human chr21, a stand-in |

Negative controls excluded unless stated. Loads at 30 C, phi29, one mismatch
class, through `occupancy.weighted_site_load`.

### Within-study ranking

Area under the ROC curve; 0.5 is no information and 0.0 a perfect inversion.
Negative controls excluded.

| predictor | Mtb | Wolbachia | Prevotella | mean |
|---|---|---|---|---|
| selectivity density, uniform | 0.96 | 0.00 | 0.33 | 0.43 |
| selectivity density, position-dependent | 0.96 | 0.00 | 0.22 | 0.39 |
| selectivity density, position-dependent-3prime | 0.83 | 0.00 | 0.33 | 0.39 |
| published fg/bg ratio (control) | 0.00 | 1.00 | 0.89 | 0.63 |
| panel size (control) | 0.58 | 0.75 | 0.72 | 0.69 |

Including the two deliberate *M. tuberculosis* negative controls: 0.97 / 0.97 /
0.88 in the Mtb column, means 0.43 / 0.40 / 0.40, panel size 0.72. No ordering
changes. `selectivity` and `selectivity_density` rank identically in every
cell, because within one study the two differ by a constant factor.

The model changes nothing on *M. tuberculosis*, changes the *Prevotella*
ranking for the worse, and cannot change *Wolbachia*, where every quantity
computed here is a perfect inversion: the one panel the authors called
effective (`TmL/Even`) carries the **highest** host load of the five, 18,598
sites against 552 for the panel they called selective. That is the same shape
the 2026-09-16 record reported for two quantities from this same paper, and it
is the finding that bounds everything else here.

### Leave-one-study-out thresholded accuracy

The threshold is chosen on two studies and applied to the third, for each study
in turn. Negative controls excluded.

| predictor | held-out accuracy | majority baseline |
|---|---|---|
| published fg/bg ratio (control) | 0.67 | 0.62 |
| selectivity density, position-dependent | 0.57 | 0.62 |
| selectivity density, position-dependent-3prime | 0.57 | 0.62 |
| selectivity, any model | 0.52 | 0.62 |
| panel size (control) | 0.52 | 0.62 |
| selectivity density, uniform | 0.48 | 0.62 |

Including the negative controls: 0.70, 0.61, 0.61, 0.48, 0.52, 0.48, against a
baseline of 0.65.

**Nothing computed here clears the baseline.** The position-dependent model is
9 points above the uniform model and 5 points below the majority baseline;
on 21 panels those are two panels and one panel respectively. The published
fg/bg ratio, which is the papers' own number and not this tool's, is the only
predictor above the baseline, by 5 points.

The count-ratio form sits at or below the baseline under every model, which is
Known Issue 6 in its usual shape, and two of the three studies are scored
against a host sixty-six times smaller than the published one.

### Per-panel loads

Foreground over background site load, by model.

| panel | effective | uniform | position-dependent | +3prime |
|---|---|---|---|---|
| Mtb4 | yes | 8770/3384 | 3068/720 | 2057/445 |
| Mtb6 | yes | 8485/3179 | 2748/824 | 1866/446 |
| Mtb8 | yes | 7251/3343 | 3027/714 | 1987/472 |
| Mtb9 | yes | 6864/2823 | 2820/664 | 1791/454 |
| Mtb1 | no | 3260/2871 | 1052/405 | 791/308 |
| Mtb2 | no | 3405/2244 | 1197/422 | 922/258 |
| Mtb3 | no | 3527/2442 | 1080/391 | 832/284 |
| Mtb5 | no | 4263/2710 | 1109/616 | 846/296 |
| Mtb7 | no | 7737/3423 | 2677/722 | 1820/388 |
| Mtb10 | no | 4198/3230 | 1182/975 | 863/783 |
| MtbSparse | control | 612/1114 | 259/171 | 186/121 |
| MtbUneven | control | 1798/886 | 404/322 | 334/217 |
| TmL/Even | yes | 465/18598 | 241/4661 | 213/3099 |
| TmL/Selective | no | 507/14656 | 312/3538 | 295/2337 |
| TmH/Even | no | 352/2975 | 220/1557 | 210/1259 |
| TmH/Selective | no | 284/552 | 201/334 | 198/309 |
| Leichty and Brisson 2014 | no | 213/4952 | 150/1458 | 145/934 |
| Prev03 | yes | 584/1949 | 201/230 | 153/168 |
| Prev04 | yes | 1402/5217 | 642/766 | 477/542 |
| Prev06 | yes | 878/2188 | 319/465 | 219/317 |
| Prev01 | no | 79/241 | 31/30 | 22/18 |
| Prev02 | no | 435/1395 | 186/202 | 133/100 |
| Prev05 | no | 403/1230 | 174/234 | 122/184 |

The position-dependent model cuts background load harder than foreground load
everywhere, which is the direction it was built to have: a background site
carries a mismatch and a foreground site does not. On the human-chr21 studies
background load falls by a factor of 4 to 7 while foreground falls by a factor
of 2 to 3. That it nonetheless does not rank better is the result.

### What could not be scored, and why

| dataset | why not |
|---|---|
| Clarke 2017 Mtb, Dwivedi-Yu Prevotella, against the published host | hg38 is on disk but no hg38 k-mer table exists at k=7, 8, 9 or 11, and counting one was out of scope for this run. **Both are scored against human chr21 instead, and every figure above carries that substitution.** Known Issue 6: a partial background is not a specific design, and `selectivity_density` is the figure to compare for that reason |
| `dwivedi_yu_2023_leishmania` | 4 panels, none with the published per-panel geometry the controls need. No `published` block in the dataset |
| `leichty_brisson_2014` | 3 panels, no `published` block, and the outcome is recorded as a fold range (`>1e5`) rather than a label |
| `oyola_2016_pfalciparum` | 1 panel, all of it effective. A single-class sample has no ranking statistic |
| `dwivedi_yu_2023_plasmid_amplification` | measurements of plasmid amplification, not primer sets. No `sets` block |

So four of the seven datasets in `tests/validation/data/` carry no label-plus-
geometry pair and are not scorable by this methodology at all, which is the same
subset the 2026-09-16 record used.

## 4. Saturation of the host load

Known Issue 6 records that 99.7% of canonical 12-mers occur in the human
genome, so allowing one mismatch could make nearly every 12-mer's neighbourhood
present and the host load could stop separating candidates.

**Measured against human chr21, it does not.** 1,000 12-mers from the
checked-in candidate pool (`tests/validation/genomes/step3_df.csv`), host load
against `human_chr21`:

| model | median | p10 | p90 | mean | sd | coefficient of variation |
|---|---|---|---|---|---|---|
| uniform | 7.98 | 3.99 | 14.94 | 9.00 | 5.10 | 0.566 |
| position-dependent | 5.40 | 2.43 | 10.24 | 5.99 | 3.41 | 0.570 |
| position-dependent-3prime | 5.00 | 2.11 | 9.51 | 5.53 | 3.23 | 0.583 |

Re-measured after the terminal-domain fix of section 1b. The figures moved in
the third decimal place (position-dependent mean 5.983 to 5.990, coefficient of
variation 0.570 both times), so nothing here turned on the defect.

The spread is essentially unchanged: the position-dependent model shifts the
whole distribution down by about a third and separates candidates neither
better nor worse. Some oligos have a host load of exactly 0.0, so the max/min
ratio is undefined and the coefficient of variation is the figure to read.

**This is chr21, not hg38, and the question the plan asks is about hg38.** A
3.30 Gb host is 66 times more sequence and the saturation argument is about
exactly that difference, so this measurement does not answer it. It was NOT run
against hg38, deliberately: no hg38 12-mer table exists on this machine and
counting one is minutes of CPU and a 138 MB table (Known Issue 1). Once such a
table exists at prefix `PREFIX`:

```
python scripts/benchmarking/mismatch_model_ranking.py --saturation \
    --pool tests/validation/genomes/step3_df.csv \
    --host-prefix PREFIX --k 12 --limit 2000
```

Expect the peak resident set to be roughly 1.2 GB larger than the chr21 run's
940 MB, because `mismatch_counts.load_kmer_counts` materialises the whole table
as a dict and an hg38 12-mer table holds close to all 8.39 million canonical
12-mers.

## Cost

Measured on 20 random 12-mers against a 200 kb synthetic table, median of five
repeats, table already cached.

| depth | model | per primer | site groups per primer |
|---|---|---|---|
| 1 | uniform | 0.046 ms | 37 |
| 1 | position-dependent | 0.126 ms | 37 |
| 1 | position-dependent-3prime | 0.116 ms | 37 |
| 2 | uniform | 0.867 ms | 631 |
| 2 | position-dependent | 2.27 ms | 631 |
| 2 | position-dependent-3prime | 2.27 ms | 631 |

The position-dependent model costs about 2.7x the uniform one at the shipped
depth of one mismatch; the neighbour enumeration and the count lookups are
shared, and the extra is the nearest-neighbour walk per neighbour. Skipping the
walk entirely for a terminal mismatch (section 1b) made both models faster in
absolute terms than the figures first recorded, and left the ratio where it was.

**Distance 2 costs about 19x distance 1**, in both models, and the group count
goes from 37 to 631 as the plan predicted (1 exact + 36 at distance 1 + 594 at
distance 2 for a 12-mer). It stays behind the existing `max_mismatches` key,
default 1, and no CLI flag was added. The loads also diverge far more at depth
2 (uniform 270 against position-dependent 108 on the synthetic table), because
the uniform model charges 8 C for two mismatches where the per-neighbour walk
charges 20 to 35 C for two interior ones.

## What this agrees with

[published_primer_sets.md](published_primer_sets.md) already records that the
*M. tuberculosis* and *Prevotella* benchmarks contradict each other about which
statistic separates the effective panels: mean binding distance separates on
Mtb and not on Prevotella, while worst gap and Gini do the opposite. The
columns above are the same disagreement reached from a different quantity, and
they add the *Wolbachia* study as a third direction. That record's conclusion,
that four datasets of 6, 5, 12 and 4 sets cannot support a fitted weighting,
applies to this measurement unchanged: 21 panels cannot establish a mismatch
model either.

## Constants this rests on

- `thermodynamics.DELTA_G_MISMATCH`, 64 doublets attributed to SantaLucia
  (1998). **Not independently verified**, and it sits in a block of
  `thermodynamics.py` headed "Legacy API compatibility for rf_preprocessing" --
  that is, a table kept for the retired random-forest feature path. The
  registry's own `provenance_note` says the primary literature was not re-read
  when it was compiled, and the ten CANONICAL stacks are the part the repository
  has checked; this 64-entry mismatch set is covered by neither a check nor a
  recorded correction. If those values move, every figure in this record moves
  and it should be re-measured rather than edited.
  The set is INTERNAL-mismatch only, and section 1b is what that cost: a
  terminal position is outside its domain and now falls back to the uniform
  penalty rather than being scored from a table that does not describe it.
- `NN_INIT_CORRECTIONS`, reproduced exactly from the existing walk so that the
  perfect-duplex equality holds. Also unverified.
- **No dangling-end or terminal-mismatch parameter set** (Bommarito et al. 2000
  would be the applicable one) exists in this repository. That absence, not a
  choice, is why a terminal mismatch takes the uniform penalty. Adding that set
  is what would let the model grade a terminal position, and it would change
  every figure here.
- `conditions.calculate_effective_tm` and `calculate_enthalpy_entropy`, both
  unchanged and used exactly as the uniform path uses them. The model therefore
  inherits whatever those paths do with additives, including the DMSO
  coefficient recorded as sitting below its own citation. The figures here were
  measured with no additives, so none of that enters them.
- `THREE_PRIME_EXTENSION_FACTOR = 0.05` is not a measurement. It states
  "largely a non-site" as a number so the claim is visible and can be
  overturned; see its registry record.

## Limitations

- Twenty-one panels and eight positives, three studies. Every figure here has
  an interval wide enough to contain the baseline, and the largest gap between
  the two models on either statistic is 9 points, which is two panels.
- An earlier pass of this measurement reported the model ahead on held-out
  accuracy. That was an artefact of a defect found in review (section 1b), and
  it is worth recording as a limitation of the method and not only as a fixed
  bug: a 9-point difference on 21 panels is the size of thing a single
  arithmetic error in one position of the oligo can produce.
- Two of the three studies are scored against human chr21 in place of the whole
  human genome. That is a declared substitution, not a measurement of the
  published design, and `selectivity_ratio` in particular moves with background
  size.
- One reaction (30 C, phi29, no additives) is used for every panel, because the
  papers do not state conditions fully enough to reconstruct each one. Three of
  the studies did use a phi29-class isothermal amplification.
- The *Wolbachia* study inverts every quantity computed here, and it has one
  positive panel, so its area under the curve is the rank of one observation.
- The outcome labels are not one endpoint: two are the authors' judgement of
  which sets were most effective and one is reads on target against a control.
- The duplex term is an extrapolation of a 37 C 1 M Na+ table, and the
  3'-terminal factor of 0.05 is a stated claim rather than a measured number.
  Neither is a measurement of mismatch discrimination, and a load computed
  under either is reported with its own `selectivity_mode` so it cannot be read
  as the uniform one.
- The saturation check is chr21, and the question it was written for is hg38.

## Decision gate C

The plan's criterion, read against the figures above and nothing else.

**The model is not shown to rank the measured outcomes better.** Within study
it is lower than the uniform model in two of the three studies' terms (mean
area under the curve 0.39 against 0.43; equal on *M. tuberculosis* for the
duplex term, lower on *Prevotella*, and *Wolbachia* is a perfect inversion
under every quantity computed here). Held out by study it is higher, 0.57
against 0.48, which on twenty-one panels is a margin of about two panels, and
both figures sit below the 0.62 majority baseline.

Both branches of the gate reach the same disposition, which is what the plan
says: it stays reachable and off by default. The branch taken is the second
one -- "no benefit demonstrated on these datasets" -- and that is a result
worth having, because it bounds what the uniform simplification costs. On this
evidence, not more than a couple of panels either way.

**A second panel would be needed before any change of default.** Nothing here
would support switching it on, and the plan's first branch would not have
either: it says explicitly that a better ranking still leaves the model off
until a second panel confirms it.

**The size of the claim.** Three of the seven datasets in
`tests/validation/data/` are scorable by this methodology; the other four carry
no label-plus-geometry pair. Twenty-one panels, eight of them effective, one
reaction (30 C, phi29, no additives) for every panel, one polymerase class
(phi29-type isothermal amplification), and two of the three studies scored
against human chr21 standing in for the hg38 their papers used. The outcome
labels are three different endpoints pooled. No figure here has an interval
that excludes the baseline.

## Position

The model is built, reachable, separately switchable in its two terms, and off
by default. Phase 4's affected-site classification may now report a graded
figure beside the binary one; it does not replace it, because the binary
reading needs no model.
