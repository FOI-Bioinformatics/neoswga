# Fitting the per-primer reach to a measured outcome

**Date:** 2026-08-19
**Data:** Clarke EL, Sundararaman SA, Seifert SN, Bushman FD, Hahn BH, Brisson D
(2017) *Bioinformatics* 33(14):2071-2077, Table 1 / Supplementary Table 1
(*Wolbachia pipientis* sets) with outcomes from Results 3.1.
**Genome:** *W. pipientis* wMel, 1,267,782 bp, 35.2% GC, via
`scripts/fetch_reference_genomes.py wolbachia`.
**Fixture:** `tests/validation/data/clarke_2017_wolbachia.json`

Related: [tool_comparison.md](tool_comparison.md) records why NeoSWGA's reach
convention differs from swga 2.0's; this fits the number rather than arguing it.

## Why this was the open question

Coverage is meaningless without the reach it is computed at. The same primer set
measures 0.418 at 3 kb, 0.836 at 10 kb and ~1.0 at 70 kb, so a coverage figure
is a statement about an assumption as much as about a design. NeoSWGA selects at
a 3 kb realistic per-primer reach and reports `headline_coverage` at a 10 kb
measured product length, and **neither had been measured**. `calibrate-reach
--bam` exists to fit it from sequencing depth and has never been run: there is no
BAM in the repository.

Clarke's *Wolbachia* experiment is the substitute. It publishes primer sequences
alongside one quantitative breadth outcome, which is enough to ask: what reach
makes the model reproduce what the reaction did?

## The outcome being fitted

Five sets, one with a number attached:

| set | n | published mean gap | outcome |
|---|---|---|---|
| `TmL/Even` | 7 | 6853 bp | **10x depth over 60-75% of the genome** |
| `TmL/Selective` | 9 | 5327 bp | partial: improved coverage, not efficiency; uneven |
| Leichty & Brisson 2014 | 2 | 5305 bp | partial: same |
| `TmH/Even` | 2 | 13070 bp | ineffective |
| `TmH/Selective` | 2 | 12074 bp | ineffective |

Unamplified control: 10x over 2.8%.

**The assumption, stated plainly:** breadth at >=10x is treated as a proxy for
"was amplified". It is not the same quantity — whether a region within reach
reaches 10x also depends on total sequencing and on evenness — but the control
at 2.8% says the amplification, not the sequencing, is what moved it.

## Validating the site-finding first

Binding sites were found by direct string search for each primer and its reverse
complement on the plus strand, treating the chromosome as circular. Comparing
the implied site count (genome length / published mean gap) with what that
finds:

| set | sites found | implied by the paper | ratio |
|---|---|---|---|
| Leichty & Brisson 2014 | 238 | 239 | 1.00 |
| `TmH/Even` | 96 | 97 | 0.99 |
| `TmH/Selective` | 104 | 105 | 0.99 |
| `TmL/Selective` | 441 | 238 | 1.85 |
| `TmL/Even` | 324 | 185 | 1.75 |

Three of five agree to within one site, which is a strong check on the method
and on the genome build. The two that do not are both the multi-primer TmL sets,
and no convention tried reconciles them: pooling one strand instead of two
overshoots (7592 against 6853 for `TmL/Even`), and averaging per-primer means
is an order of magnitude out. There are no reverse-complement pairs within
either set, and every primer is found with a plausible per-primer count.

The likeliest reading is that the Table 1 statistics for those two rows were
computed over a smaller primer list than Supplementary Table 1 supplies. That is
not established, and the calibration below carries it as an uncertainty rather
than resolving it.

## The fit

Union coverage of the symmetric extension windows, the convention the tool uses:

| set | 1 kb | 2 kb | 3 kb | 4 kb | 6 kb | 8 kb | 10 kb | outcome |
|---|---|---|---|---|---|---|---|---|
| `TmL/Even` | 0.28 | 0.48 | 0.62 | 0.72 | 0.84 | 0.92 | 0.96 | **effective** |
| `TmL/Selective` | 0.26 | 0.45 | 0.60 | 0.71 | 0.85 | 0.92 | 0.96 | partial |
| Leichty 2014 | 0.20 | 0.36 | 0.48 | 0.58 | 0.71 | 0.80 | 0.85 | partial |
| `TmH/Even` | 0.12 | 0.24 | 0.34 | 0.43 | 0.57 | 0.68 | 0.76 | ineffective |
| `TmH/Selective` | 0.11 | 0.20 | 0.28 | 0.36 | 0.49 | 0.59 | 0.68 | ineffective |

`TmL/Even` enters the measured 60-75% band at **3000 bp** and leaves it at
**4500 bp**. Under the site count the paper's own statistic implies (185 rather
than 324, sampled down and averaged over five draws) the band moves to
**4100-6200 bp**.

**Calibrated range: 3.0-6.2 kb**, centred near 4-5 kb.

Two things follow.

**The 3 kb selection reach is corroborated**, sitting at the bottom of the range.
It was derived from the inter-primer spacings published tools designed for — a
designer's choice, not an observation — and it survives contact with an outcome.

**The 10 kb reporting reach is not.** `headline_coverage` uses the measured MDA
product length, and that is a defensible quantity in itself, but it is outside
the band this outcome supports. A set reported at 10 kb reads about 1.5x higher
than the same set at the calibrated reach.

## It discriminates, which is the real test

An absolute fit could be a coincidence. The ordering is harder to get right by
accident. At any reach inside the calibrated band the two ineffective sets sit
at 0.28-0.43 while the effective one sits at 0.62-0.72 — they need 6.8-8.3 kb
just to reach 60%. The model separates effective from ineffective at the
calibrated reach without being asked to.

It does **not** separate `TmL/Even` from `TmL/Selective`: 0.62 against 0.60 at
3 kb, and the second was only partially effective. That is expected and is the
paper's own finding — those two were built as a controlled contrast between
evenness and selectivity, and evenness is what separated them. Coverage-within-
reach is not an evenness metric, and `max_gap` and `gap_gini` exist for that.

## What this does to the completeness question

Re-running the *Prevotella melaninogenica* sweep (3,168,282 bp against human
chr21, 12-mers, background-clean candidates deduplicated by reverse-complement
pair, greedy maximum coverage) inside the calibrated band:

| reach | >=90% | >=95% | >=99% | coverage at n=8 / 16 / 32 / 64 |
|---|---|---|---|---|
| 3.0 kb | 76 | **96** | 126 | 0.44 / 0.54 / 0.67 / 0.86 |
| 4.5 kb | 40 | **51** | 70 | 0.60 / 0.70 / 0.85 / 0.98 |
| 6.2 kb | 24 | **31** | 43 | 0.72 / 0.83 / 0.96 / — |

**A >95% design needs 31-96 primers.** The earlier figure of ~15, taken at a
10 kb reach, is not supported by this calibration.

The specificity half is not the constraint. The candidate pool used here is
219,730 distinct 12-mers with **zero** exact matches anywhere in chr21, so
taking 96 of them forces no background binding by construction. Coverage and
specificity are not in tension at these sizes on this target/background pair;
what binds is how many oligos you are willing to buy.

Which makes the set-size ceilings the operative limit. `num_primers` and
`target_set_size` were schema-capped at 50 — below most of the measured range —
and are now capped at 200. The `--auto-size` ceiling of 20 is **not** raised,
for a reason recorded in the next section.

## A defect found while doing this

`estimate_optimal_set_size` is genome-independent. `genome_length` appears in
`expected_sites` and again in the denominator of `coverage_per_primer`, so it
cancels, and the estimate depends only on primer length:

| primer length | estimated N, any genome from 6 kb to 50 Mb |
|---|---|
| 8 | 53 |
| 10 | 851 |
| 12 | 13,612 |

At every SWGA primer length that is above the clamp, so the clamp has always
been the entire answer and the model has contributed nothing.

The cause is the random-sequence site model: `expected_sites` assumes a k-mer
occurs 4^-k times per base, while SWGA candidates are selected precisely for
being abundant in the target. The Prevotella 12-mer pool has a median of 8 sites
per primer against the 0.38 that model predicts — a factor of about 20. Primers
chosen for non-randomness cannot be modelled as random.

The ceiling stays at 20 rather than being raised, because `--auto-size` output
is spent on oligos and a wrong number that is small does less damage than a
wrong number that is large. Recalibrating the model against the table above is
the fix; `num_primers` in params.json is the way to ask for more meanwhile.

## What is not established

- **One target, one outcome.** The reach is fitted to a single set with a single
  published breadth figure. Two of the four other sets corroborate the ordering
  but carry no number.
- **Breadth at 10x is not "within reach".** The proxy is argued above and is not
  proven.
- **The site-count discrepancy on the two TmL sets is unexplained**, and it is
  the reason the range is 3.0-6.2 kb rather than 3.0-4.5.
- **The symmetric window is still an approximation.** Extension from a binding
  site is one-sided and strand-determined; measured on Prevotella, symmetric
  overstates coverage by 4-14% (growing with set size, shrinking with reach).
  Correcting it would move the fitted reach up slightly.
- **A BAM would settle this properly.** `calibrate-reach --bam` fits the reach
  from sequencing depth directly, without the breadth-at-depth proxy or the
  single-outcome dependence.
