# Using additives to buy specificity: what the model says, and what it took to see it

NeoSWGA exists to design phi29/equiphi29 SWGA sets, and the lever of interest is
using reaction additives to obtain better specificity. This document records
what happened when that lever was built and then run end to end against real
genomes, including two defects that had to be fixed before any of it was
measurable and three ideas that measurement rejected.

Everything below is measured on one target/background pair unless stated
otherwise. Reproduce with `scripts/fetch_reference_genomes.py prevotella
human_chr21`, then the `tests/validation/genomes/equiphi29_{k}mer.json` and
`phi29_{k}mer.json` parameter files. The tests that pin these numbers are in
`tests/validation/test_prevotella_specificity.py` and skip without the genomes.

| | |
|---|---|
| target | *Prevotella melaninogenica* ATCC 25845, 3,168,282 bp, 41% GC |
| background | human chr21 (hg38), 46,709,983 bp |
| platforms | equiphi29 at 42 C, phi29 at 30 C |
| set size | 8 primers unless stated |

A caution that applies to the whole document: **none of this is validated
against a measured wet-lab outcome.** No dataset in the repository carries
per-set yield at these primer lengths, and the one benchmark with a continuous
specificity outcome cannot serve, for the reasons in
[published_primer_sets.md](published_primer_sets.md#specificity-validation-against-prevotella-a-negative-result-and-why).
What follows is internally consistent and mechanistically grounded, and it is
not evidence about reactions.

## Why the lever did not exist before

Selectivity was a count of exact k-mer matches from the HDF5 position files:

```python
total_fg = len(fg_positions)
total_bg = len(bg_positions)
selectivity = total_fg / max(total_bg, 1)
```

No temperature, no free energy, no additive term. For a fixed primer set, DMSO
and betaine could not move it by any amount. Two absences, which are one absence
seen twice: background near-matches were not counted, and site occupancy was not
modelled, so every site contributed 1.0 regardless of whether it was bound.

The replacement weights each site by two-state occupancy,

```
theta(s) = 1 / (1 + exp( (dH/R) * (1/T - 1/Tm_eff(s)) ))
selectivity = SUM_j n_fg,j * theta_j  /  SUM_j n_bg,j * theta_j
```

with sites grouped by mismatch class `j` and `theta_j` evaluated at
`Tm_eff - j * mismatch_penalty`. Counting the classes is cheap because the
jellyfish `*_all.txt` files hold every k-mer: one-mismatch load is a sum over
the primer's `3k` neighbours.

## The regime where additives work at all

Additives lower effective Tm, which discriminates only where a duplex is near
its melting transition. Far below it, perfect and mismatched sites are both
saturated and the occupancy terms cancel in the ratio. That makes the lever a
property of primer length, and the effect is sharp.

Occupancy selectivity across 38-50 C, equiphi29:

| set | length | 38 C | 42 C | 46 C | 50 C | fold |
|---|---|---|---|---|---|---|
| Prev03 (published) | 8 | 0.30 | 0.30 | 0.30 | 0.30 | **1.00x** |
| Prev06 (published) | 7, 8 | 0.41 | 0.41 | 0.41 | 0.41 | **1.00x** |
| pipeline | 10 | 0.53 | 0.61 | 0.68 | 0.71 | 1.36x |
| pipeline | 11 | 0.75 | 0.89 | 1.14 | 1.31 | 1.74x |
| pipeline | 12 | 2.00 | 2.35 | 3.54 | 6.05 | **3.03x** |

The published 7-8mer sets are flat to two decimal places at every temperature.
The structural reason: all 8,192 canonical 7-mers and all 32,896 canonical
8-mers occur in *both* a 3 Mb bacterium and a 46 Mb human chromosome, at means
of ~97 and ~1219 occurrences. At that length a primer's background load is set
by how much sequence it is counted against, not by its sequence.

**Additives cannot buy specificity for a set of 7-8mers.** There is no
mismatched population to discriminate against when the perfect matches already
number in the thousands. The lever needs primers long enough for exact matches
to be rare, which is the equiphi29 regime rather than phi29's 6-12mers.

The fold value is set-dependent (a six-primer 12-mer set measured 2.05x where an
eight-primer one measured 1.51x), so no particular fold is worth relying on.
What is stable is the gap: everything at 10-12 responds, everything at 7-8 is
flat, nothing in between.

## Why the platform is equiphi29 and not phi29

The same argument applies to temperature, and isolates cleanly. One 12-mer set,
both platforms, 5% DMSO + 1 M betaine:

| platform | selectivity | with additive | | effective coverage cost |
|---|---|---|---|---|
| phi29 30 C | 1.82 | 1.84 | **1.01x** | -0% |
| equiphi29 42 C | 2.35 | 3.54 | **1.51x** | -9% |

At 30 C a 12-mer melting near 42 C is far below its Tm; lowering Tm by a few
degrees moves neither the perfect nor the mismatched population. The additive
lever is a property of operating near Tm, not of the additive, and phi29's
temperature puts it out of reach.

Additives and temperature are near-interchangeable where both are available:
46 C with no additive gives the same selectivity (3.54) and effective coverage
(0.348) as 42 C with DMSO + betaine (0.352). Additives are the lever when the
enzyme's ceiling stops you raising temperature, which for equiphi29 is 45 C.

## Both halves of the objective have to be on one footing

Selectivity became occupancy-weighted while coverage stayed a count of binding
sites, treating every exact match as fully occupied. That rewarded a set for
reaching sites it does not occupy at the reaction temperature, and it made the
coverage cost of additives invisible.

`effective_fg_coverage` unions the same binding windows, with each site covering
its window with probability theta:

```
P(covered at x) = 1 - PRODUCT over sites s reaching x of (1 - theta_s)
```

Sites of one primer are marked together before being applied, because
overlapping windows of the *same* primer are one opportunity to prime, not
several. Independence across primers is an approximation: sites on one template
molecule compete for polymerase, so this reads as an upper bound per molecule.

The bias is large and not uniform, so it does not cancel in a comparison. At
equiphi29 42 C:

| set | raw coverage | effective | mean occupancy |
|---|---|---|---|
| k=10 | 0.520 | 0.245 | 0.36 |
| k=11 | 0.324 | 0.220 | 0.56 |
| k=12 | 0.403 | **0.387** | 0.86 |

Raw coverage ranks k=10 first, effective coverage ranks k=12 first. Selection
now uses the effective figure when reaction conditions are available, falling
back to raw otherwise (0.0 means "not computed", not "nothing covered").

With the cost visible, the design rule is two lines:

| additive: 5% DMSO + 1 M betaine, 42 C | effective coverage | selectivity |
|---|---|---|
| k=12 | -9% | +51% |
| k=10 | -64% | +13% |

**Stringency is affordable when the primers melt above the reaction temperature
and ruinous when they do not.**

## Two defects that had to be fixed first

### params.json's set size never reached the optimizer

`optimize` resolved its target primer count from `parameter.num_primers` before
anything had read params.json. That global defaults to 6 and is only overwritten
when `get_params()` runs, which happens lazily inside `optimize_step4` -- after
the target had been resolved and passed down. A params.json asking for 8 primers
silently got 6.

The failure was quiet in the worst way. The completion check reads
`parameter.target_set_size`, and by then `get_params()` *has* run, so it compared
the 6 produced against the 8 requested and advised relaxing the filters and
enlarging the candidate pool. Neither was the problem. Acting on that advice
would have degraded the design to find primers that were never missing.

### the candidate cut was a slice, not a selection

`step2` ranks survivors by `ratio = bg_count / fg_count` and keeps the best
`max_primer`. That key is correct but it is not a total order, and in this
tool's own regime it is barely an order at all: a primer long enough to have no
exact background match has ratio `0/fg_count = 0`.

At k=12 **all 369,431 survivors tied at exactly 0.0**. `sort_values` is stable,
so the cut returned whichever 10,000 arrived first: mean foreground count 1.29,
against 9.5 for the 10,000 most abundant, with 10 primers in common between the
two. Coverage is set by how many sites a design binds, so the optimizer was
handed single-site candidates and coverage was capped before it began.

Breaking the tie on foreground abundance, end to end at k=12:

| | before | after |
|---|---|---|
| coverage | 7.6% | **40.3%** |
| foreground sites | 43 | 528 |
| selectivity | 0.69 | **2.54** |
| max gap | 326 kb | 108 kb |

Both halves improve, which is the signature of a defect rather than a frontier:
the extra sites are foreground sites and the background count stayed at zero.
k=10 is unchanged, as it should be -- its background counts are nonzero, so its
ratios were never all-tied.

A correction this forced: raising `max_primer` from 2,000 to 10,000 had appeared
to buy +43% coverage. With the tie broken, the 2,000 pool matches the 10,000
pool exactly. That earlier gain was extra draws from an arbitrary slice.

## The coverage ceiling, and what it is not

After both fixes, effective coverage at k=12 with 8 primers is 0.39. Four checks
say this is not an optimizer defect:

- hybrid and dominating-set return **identical** sets, so the two-stage
  narrowing is not discarding coverage
- an independent greedy union-coverage selection over the same pool reaches
  0.418 against the optimizer's 0.403
- marginal gain collapses after the *first* primer: 0.290, then 0.022, 0.021,
  0.018, ...
- among the 250 most abundant candidates the median site count is 8, and spread
  efficiency (actual coverage / if-evenly-spaced) has median **1.00**

That last one refuted the working hypothesis. Three of the eight selected
primers are one repeat family (`AAAAGGGCGTTA`, `AAGGGCGTTAGT`, `CGCCCTTTTGAA`
share `GGGCGTTA`), so redundant abundance looked like the culprit -- but the
sites are well spread. There simply are not many of them: only 2,286 of 2.18M
canonical 12-mers occur 10 or more times in 3.2 Mb.

**The ceiling is the 12-mer abundance distribution.** Set size is the lever that
works:

| primers | effective coverage | selectivity | max gap |
|---|---|---|---|
| 8 | 0.382 | 2.54 | 108 kb |
| 16 | 0.435 | 2.09 | 79 kb |
| 24 | 0.488 | 1.63 | 62 kb |
| 32 | 0.539 | 1.39 | 46 kb |

Set size and additives trade in opposite directions, so they compose: n=32 with
DMSO + betaine gives effective coverage 0.474 at selectivity 1.82.

## Three ideas measurement rejected

**Mixed-length pools do not help at equiphi29.** Length and the Tm window are
not independent knobs. With `min_k=10, max_k=12` and a Tm floor at the reaction
temperature, the filter loaded 3.6M candidates across all three lengths and
returned an all-12-mer pool, because 10-mers melt at 39-42 C and cannot clear a
42 C floor. You cannot buy the abundance of short primers at a temperature where
they are not bound.

**Raising salt helps occupancy, not selectivity.** `mg 20 mM + na 200 mM` raises
the 10-mer set's mean Tm 40.4 -> 42.9 C and effective coverage 0.245 -> 0.361
(+47%), but its selectivity stays at 0.55. At k=12 it buys almost nothing:
0.387 -> 0.396, with selectivity falling 2.35 -> 2.08. Note that 20 mM Mg is
above the usual phi29 working range and the model does not capture enzyme
inhibition there, so this needs wet-lab confirmation before anyone acts on it.

**Enlarging the candidate pool stopped helping** once the cut was fixed, as
above. `max_primer` is schema-capped at 10,000; the filter refuses to run above
it rather than silently truncating.

## Platform summary

8 primers, same target and background:

| k | phi29 30 C eff cov | sel | equiphi29 42 C eff cov | sel |
|---|---|---|---|---|
| 10 | **0.498** | 0.40 | 0.245 | 0.61 |
| 11 | 0.338 | 0.73 | 0.220 | 0.89 |
| 12 | 0.403 | 1.82 | 0.387 | **2.35** |

phi29 delivers the highest effective coverage measured anywhere here -- its
lower temperature keeps abundant 10-mers occupied -- and cannot make any of it
specific. Its k=12 effective coverage equals its raw coverage exactly, because
at 30 C occupancy is ~1.0, which is precisely why there is no discrimination to
be had. equiphi29 gives up a little coverage at k=12 and is the only platform
where the additive lever exists.

The best design measured: **equiphi29 42 C, k=12, 8 primers, 5% DMSO + 1 M
betaine** -- effective coverage 0.352, selectivity 3.54.

## What is not established

- No wet-lab validation, as above. This is a model that is now internally
  consistent, not a predicted yield.
- The amplification figures that set the additive ceiling come from the
  mechanistic model's enzyme pathway, whose magnitudes are marked EMPIRICAL in
  [SCIENCE_CITATIONS.md](../SCIENCE_CITATIONS.md). Treat the coverage and
  selectivity columns as far firmer than any claim about reaction yield.
- One target/background pair. The all-zero `ratio` tie is guaranteed wherever
  background counts vanish, but how much the abundance ordering helps will
  depend on the genome.
- Mismatch classes are capped at 1 by default. Two-mismatch load is 405 lookups
  for a 10-mer rather than 30, and has not been swept.
