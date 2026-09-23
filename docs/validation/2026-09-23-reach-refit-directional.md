# The coverage reach, refitted under directional windows

Measured 23 September 2026.

[reach_calibration.md](reach_calibration.md) fitted the per-primer reach to
**3.0-4.5 kb** by asking which value reproduces the 60-75% breadth Clarke et al.
(2017) measured at 10x depth for their `TmL/Even` Wolbachia set. It used
symmetric windows, and records in its own limitations that extension is
one-sided and that "correcting it would move the fitted reach up slightly".

The geometry and the value were fitted together, so correcting one without
refitting the other is unsound in either direction. This refits it, using the
same five published sets, the same genome and the same outcome. No new data.

```bash
python scripts/benchmarking/refit_reach_directional.py
```

## Step 1: the existing fit rebuilds

Before changing anything. A table that cannot be rebuilt cannot be refitted.

| set | 1 kb | 2 kb | 3 kb | 4 kb | 6 kb | 8 kb | 10 kb |
|---|---|---|---|---|---|---|---|
| `TmL/Even` | 0.28 | 0.48 | 0.62 | 0.72 | 0.84 | 0.92 | 0.96 |
| `TmL/Selective` | 0.26 | 0.45 | 0.60 | 0.71 | 0.85 | 0.92 | 0.96 |
| Leichty 2014 | 0.20 | 0.36 | 0.48 | 0.58 | 0.71 | 0.80 | 0.85 |
| `TmH/Even` | 0.12 | 0.24 | 0.34 | 0.43 | 0.57 | 0.68 | 0.76 |
| `TmH/Selective` | 0.11 | 0.20 | 0.28 | 0.36 | 0.49 | 0.59 | 0.68 |

**Largest disagreement with the published table: 0.005.** The symmetric fit
reproduces as well: the band is entered at 2,866 bp and left at 4,551 bp, which
is the recorded 3.0-4.5 kb. Everything below is therefore comparable to the
record rather than to a re-derivation of it.

## Step 2: the same sets under directional windows

| set | 1 kb | 2 kb | 3 kb | 4 kb | 6 kb | 8 kb | 10 kb |
|---|---|---|---|---|---|---|---|
| `TmL/Even` | 0.20 | 0.35 | 0.48 | 0.58 | 0.72 | 0.81 | 0.88 |
| `TmL/Selective` | 0.18 | 0.33 | 0.45 | 0.55 | 0.70 | 0.79 | 0.86 |
| Leichty 2014 | 0.12 | 0.23 | 0.32 | 0.40 | 0.54 | 0.65 | 0.73 |
| `TmH/Even` | 0.07 | 0.13 | 0.19 | 0.25 | 0.36 | 0.45 | 0.52 |
| `TmH/Selective` | 0.06 | 0.11 | 0.16 | 0.20 | 0.27 | 0.34 | 0.39 |

## Step 3: the refit

| geometry | enters the 60-75% band | leaves it |
|---|---|---|
| symmetric | 2,866 bp | 4,551 bp |
| **directional** | **4,354 bp** | **6,710 bp** |

**Calibrated directional range: 4.4-6.7 kb**, against 2.9-4.6 kb symmetric.

### The prediction was 2x and the answer is 1.52x

The plan predicted the directional band would land near twice the symmetric
one, on the reasoning that a one-sided window credits half the bases of a
symmetric window at the same reach. It lands at **1.52x**.

The reason is that coverage is concave in reach, not linear. Windows overlap,
so restoring a lost fraction of coverage costs less than doubling the reach.
The prediction assumed the per-site arithmetic carries through to the union
over a set, and it does not.

This matters for what the correction IS. At 2x, directional-at-the-refit and
symmetric-at-today's-value would have been the same total width per site, and
the change would have been pure placement. At 1.52x they are not: the
directional fit credits **4,354 bp per site against the symmetric fit's 5,732**,
about 24% less total width. The directional model reaches the same measured
breadth with less width because it stops covering, from one side, ground
already covered from the other.

The oligo's own footprint is neglected in both geometries, matching what
`dominating_set_optimizer.add_primer_coverage:205-224` already does. Including
it moves the fitted edges by 1 bp -- 4,355 against 4,354 -- which is why one
shared approximation was preferred to two models differing by less than either
one's uncertainty.

## What it does to a design at the shipped default

The shipped `reach_phi29` is 3000, which sits just above the symmetric band's
lower edge of 2,866. The directional equivalent of that position is about 4,400.
Pairing them:

| set | outcome | symmetric @3000 | directional @4400 | delta |
|---|---|---|---|---|
| `TmL/Even` | **effective** | 0.619 | 0.609 | -1.6% |
| `TmL/Selective` | partial | 0.598 | 0.589 | -1.5% |
| Leichty 2014 | partial | 0.484 | 0.434 | -10.3% |
| `TmH/Even` | ineffective | 0.341 | 0.277 | -18.8% |
| `TmH/Selective` | ineffective | 0.283 | 0.214 | -24.3% |

**The geometry and the reach have to move together.** Switching to directional
windows while leaving the reach at 3000 puts the default below the calibrated
band: `TmL/Even` reads 0.48 there, outside the 60-75% its wet-lab outcome
recorded. That would not be a more correct model, it would be a model
mis-calibrated by construction.

### Separation between outcomes widens, and the anchor makes that hard to read

The effective-to-ineffective ratio goes from 1.82 (0.619 / 0.341) to 2.20
(0.609 / 0.277), and the gap between the weakest partial set and the strongest
ineffective one goes from 0.143 to 0.157. On its face the directional model
distinguishes the measured outcomes more sharply.

**That reading is not safe from this experiment.** The directional reach was
fitted ON `TmL/Even`, so that set's value is pinned by construction and the
widening is entirely the other four moving down. Whether that is better
discrimination or an artefact of the anchor cannot be told with one anchor.
Fitting instead to a different set, or holding out `TmL/Even` and fitting to
`Leichty`, would separate the two; neither was done here.

The ordering of the five sets is unchanged under both geometries, and coverage
still does not separate `TmL/Even` from `TmL/Selective` -- 0.609 against 0.589.
`reach_calibration.md` already records that as expected, since those two were
built as a controlled contrast on evenness and coverage-within-reach is not an
evenness metric.

## What this does not establish

- **One target, one published outcome.** Every limitation
  `reach_calibration.md` records under **What is not established** carries over
  unchanged: breadth-at-10x is an argued proxy rather than a proven one, and
  the two `TmL` rows carry an unexplained site-count discrepancy that is why
  its range is 3.0-6.2 rather than 3.0-4.5. The same discrepancy is present
  here, reproduced rather than resolved.
- **No BAM.** `calibrate-reach --bam` fits the reach from sequencing depth
  directly and has never been run; there is none in this repository. Both
  numbers above are fitted to a breadth proxy.
- **Nothing about delivered panels.** This fits a reach; it does not say
  whether a design run under directional windows selects a different set of
  oligos. That is the next measurement and it decides whether this is a
  reporting change or a design change.
- **The binned Stage-1 model is untouched.** It already selects directionally,
  at one-sided `r`, while this fit concerns the base-resolution scoring path.
  Aligning them is a separate increment.
