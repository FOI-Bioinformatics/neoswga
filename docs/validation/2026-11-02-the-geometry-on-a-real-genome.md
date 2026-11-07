# The coverage geometry on a real genome: the panel moves at every size, and the inflation is 1.5x

Measured 2 November 2026, on `main` at `c9649b6`. Phase 3 of
[the plan](../reports/PLAN_2026-11-02_known_issues.md), against Known Issue 25
and step 1 of "What remains for whoever decides to switch" in
[the geometry decision](2026-10-29-the-coverage-geometry-decision.md).
Reproduce with:

```
PYTHONPATH=. python docs/validation/geometry_second_pool_2026-11-02/measure.py \
    --work <dir with params.json and a prepared pool> --out results.json
```

Results and the pool's configuration in
[`geometry_second_pool_2026-11-02/`](geometry_second_pool_2026-11-02/).
**No code changed.** Nothing here moves a default; it supplies the figures the
three decisions in Phase 3 were waiting on.

## Why this run exists

`2026-09-23-does-the-panel-move.md` measured the geometry's effect on the
delivered panel on ONE pool and recorded two limitations that bound what it
could conclude:

- **One target, one background, one candidate pool** -- Prevotella against
  human chr21.
- **"The background axis was not exercised."** `bg_coverage` read exactly
  0.00000000 in every run, because that candidate pool has zero exact matches
  anywhere in chr21 by construction. `bg_coverage` is the one quantity that
  sees background site POSITION, so the geometry's effect on host binding was
  untested.

And the geometry decision itself recorded a third: **"It does not establish the
inflation on a real genome."** Every figure behind the 1.3-1.9x in Known
Issue 25 came from synthetic sequence, where sites are Poisson-distributed and
real genomes cluster.

This run addresses all three. It is on real sequence, on a second pair, with a
background that answers.

## Setup

*E. coli* K-12 MG1655 (4,639,675 bp, 1 record, circular, GC 50.8%) as the
target against human **chr1** (hg38, 248,956,422 bp, 7.4% N) as the host,
12-mers, equiphi29 at 43 C, `max_primer` 1000, `hybrid`. One
`neoswga optimize` per arm against one prepared `step3_df.csv`, four panel
sizes, four arms:

| arm | what it isolates |
|---|---|
| `symmetric` @4000 | the shipped default at equiphi29's realistic reach |
| `directional` @4000 | the geometry alone, at equal reach |
| `directional` @6080 | the refitted pairing (1.52x, from `2026-09-23-reach-refit-directional.md`) |
| `symmetric` @6080 | **new**: the reach alone, which the first pool's design could not separate from the geometry |

**Why chr1 and not chr21.** A 12-mer occurs about 15 times in 249 Mb, so
essentially no candidate avoids chr1; the shortlist's best occupancy-weighted
background load is 0.0337 rather than 0. Against chr21 the same target's
shortlist ranks 1,000 candidates with zero host sites to the top and
`bg_coverage` returns to exactly 0 -- the first pool's dead axis, reproduced.
**The host has to be large enough that avoidance is impossible before
`bg_coverage` measures anything.**

Both genomes were fetched from the AWS iGenomes mirror, which this container's
egress policy allows where NCBI, UCSC and EBI are refused, and verified by
length (both match the published assembly exactly), record count, GC and
SHA-256; `scripts/fetch_reference_genomes.py` carries the entries and the
checksums, and the pipeline's own `*_12mer_all.provenance.json` recomputed the
same digests.

## The measurement

| arm | n | fg_coverage | bg_coverage | bg sites | sel_density | Jaccard vs `symmetric`@4000 |
|---|---|---|---|---|---|---|
| `symmetric`@4000 | 4 | 0.26911 | 6.43e-05 | 2 | 139.33 | — |
| `directional`@4000 | 4 | 0.17579 | 3.21e-05 | 2 | 131.87 | **0.600** |
| `directional`@6080 | 4 | 0.21841 | 4.88e-05 | 2 | 139.33 | 1.000 |
| `symmetric`@6080 | 4 | 0.37855 | 9.77e-05 | 2 | 139.33 | 1.000 |
| `symmetric`@4000 | 8 | 0.31464 | 6.43e-05 | 2 | 133.35 | — |
| `directional`@4000 | 8 | 0.20357 | 3.21e-05 | 2 | 137.74 | **0.778** |
| `directional`@6080 | 8 | 0.29226 | 4.88e-05 | 2 | 137.74 | **0.778** |
| `symmetric`@6080 | 8 | 0.44049 | 9.77e-05 | 2 | 133.35 | 1.000 |
| `symmetric`@4000 | 16 | 0.35800 | 9.64e-05 | 3 | 137.61 | — |
| `directional`@4000 | 16 | 0.22965 | 3.21e-05 | 2 | 137.67 | **0.778** |
| `directional`@6080 | 16 | 0.33148 | 7.33e-05 | 3 | 135.93 | **0.600** |
| `symmetric`@6080 | 16 | 0.49615 | 1.95e-04 | 4 | 134.42 | **0.600** |
| `symmetric`@4000 | 24 | 0.38617 | 1.29e-04 | 4 | 102.79 | — |
| `directional`@4000 | 24 | 0.24522 | 6.43e-05 | 4 | 103.74 | **0.920** |
| `directional`@6080 | 24 | 0.35426 | 9.77e-05 | 4 | 91.64 | **0.548** |
| `symmetric`@6080 | 24 | 0.52831 | 1.95e-04 | 4 | 91.64 | **0.548** |

## 1. The panel moves at every size, not at half of them

At equal reach the delivered panel differs from the symmetric one at **four
sizes of four**: Jaccard 0.600, 0.778, 0.778, 0.920. The first pool saw two of
four. So the geometry being a design change rather than a reporting change now
holds on two pools, and on this one it holds everywhere it was asked.

## 2. The inflation on a real target is 1.49x to 1.80x

The delivered panel's own `coverage_geometries` block, across all sixteen runs:

| | min | median | max |
|---|---|---|---|
| inflation (symmetric / directional) | **1.491** | **1.546** | **1.797** |

Known Issue 25 claims 1.3x to 1.9x from synthetic sequence. The real-genome
range sits **inside** that interval, narrower and higher-centred. The entry's
headline survives the move to real sequence; its stated range is wider than
what a real target produced.

## 3. Clustering buffers the target and not the host, so the fg/bg ratio IMPROVES

This is the finding worth keeping, and it was not visible on synthetic
sequence.

Halving the window does not halve target coverage on a real genome:

| n | directional/symmetric, fg | directional/symmetric, bg |
|---|---|---|
| 4 | 0.653 | 0.500 |
| 8 | 0.647 | 0.500 |
| 16 | 0.641 | 0.333 |
| 24 | 0.635 | 0.500 |

On the synthetic pool the fg ratio was close to 0.5 throughout -- "what a
half-width window gives when the spans mostly do not overlap". On *E. coli*
it is **0.64**, because real binding sites cluster: neighbouring spans already
overlapped under the symmetric window, so removing one side of each span
destroys less coverage than its width suggests.

The host does not get that protection. The panel's two to four chr1 sites are
isolated in 249 Mb, nothing overlaps, and `bg_coverage` halves **exactly** --
3.21e-05 is 8,000 bp of 248,956,422 against the symmetric arm's 16,000.

So the two figures do not scale together, and the coverage-based specificity
ratio moves in the directional model's favour:

| n | fg/bg coverage, symmetric | fg/bg coverage, directional | change |
|---|---|---|---|
| 4 | 4,187 | 5,470 | +31% |
| 8 | 4,896 | 6,335 | +29% |
| 16 | 3,714 | 7,147 | +92% |
| 24 | 3,004 | 3,816 | +27% |

**A switch to directional does not scale every coverage figure by one factor.**
It lowers target coverage by about a third and host coverage by a half, on this
pair, and the gap between those two is a property of how clustered each genome's
sites are -- not a constant to correct for. `selectivity_density`, which counts
sites rather than covered bases, is unmoved by the geometry (133-139 at n=4-16
in every arm), as it should be.

## 4. The reach and the geometry move the panel by turns, and the fourth arm was needed

| n | `symmetric`@4000 vs `symmetric`@6080 (reach alone) | `symmetric`@4000 vs `directional`@4000 (geometry alone) | `directional`@6080 vs `symmetric`@6080 |
|---|---|---|---|
| 4 | 1.000 | 0.600 | 1.000 |
| 8 | 1.000 | 0.778 | 0.778 |
| 16 | 0.600 | 0.778 | 0.882 |
| 24 | 0.548 | 0.920 | 1.000 |

At the small sizes the geometry moves the panel and the reach does not; at the
large sizes it reverses. And at 6,080 the two geometries deliver the same panel
at n=4 and n=24. **So "the refitted pairing changes the panel" cannot be
attributed to the geometry**: at three of four sizes the panel delivered by
`directional`@6080 is the one `symmetric`@6080 delivers, and what moved it away
from `symmetric`@4000 was the reach. The first pool's design -- symmetric@4000
against directional@6080, with no symmetric@6080 arm -- could not see that, and
this arm was added for exactly that reason.

## 5. A third data point on the 1.52x factor not transferring

`directional`@6080 against `symmetric`@4000, the refitted pairing:

| n | symmetric@4000 | directional@6080 | delta |
|---|---|---|---|
| 4 | 0.26911 | 0.21841 | **-18.8%** |
| 8 | 0.31464 | 0.29226 | -7.1% |
| 16 | 0.35800 | 0.33148 | -7.4% |
| 24 | 0.38617 | 0.35426 | -8.3% |

Three targets now: Wolbachia -1.6% (where the factor was fitted), Prevotella
-18 to -25%, *E. coli* -7 to -19%. The factor is coverage-neutral only on the
target it was fitted on, which
`2026-09-23-reach-refit-directional.md` stated and this confirms on a third
genome.

## 6. The control: the same target against chr21, where the axis dies

Sixteen more runs, identical but for the host -- human chr21 (46,709,983 bp)
instead of chr1, which means a different `max_primer` shortlist because the
ranking is by background load.
[`results_ecoli_vs_chr21.json`](geometry_second_pool_2026-11-02/results_ecoli_vs_chr21.json).

| | chr1 | chr21 |
|---|---|---|
| `bg_coverage`, every arm and size | 3.2e-05 to 2.0e-04 | **exactly 0.0, 0 host sites** |
| inflation | 1.491-1.797 | 1.518-1.637 |
| fg ratio, directional/symmetric at equal reach | 0.635-0.653 | 0.606-0.645 |
| panel moves at equal reach | 4 sizes of 4 | 3 sizes of 4 |

Two things this settles:

- **The dead background axis is the host's size, not the first pool's bad
  luck.** Against chr21 the shortlist again fills with candidates that bind
  the host nowhere, and `bg_coverage` reads exactly 0 in all sixteen runs --
  the first pool's result reproduced on a different target. A geometry
  measurement that wants the background axis has to choose a host a k-mer
  cannot avoid, and this is the second pair where that is the whole
  difference.
- **The 0.64 fg ratio is the target's property.** Two different candidate
  pools over the same genome give 0.606-0.653, against about 0.50 on synthetic
  Poisson sites. So the buffering is site clustering in *E. coli*, not an
  artifact of one shortlist.

Across both hosts the panel moves at **seven of eight** sizes at equal reach.

## What this changes for the three decisions

Phase 3 lists three decisions that are the maintainer's. None is settled here;
each now has a figure it lacked.

1. **Does the default move at all?** The re-baselining cost is now bounded: a
   recorded `fg_coverage` would fall by a factor of **1.49 to 1.80** on real
   sequence at equal reach, or by **7 to 19%** under the refitted pairing. It
   is not a uniform rescale, so a record quoting a coverage figure cannot be
   corrected arithmetically; each has to be re-run or re-stated.
2. **Scope by polymerase?** Unchanged by this run, and still the blocker: the
   directional band was fitted for phi29 against one published outcome, and
   this pool is equiphi29 at the 1.52x-scaled reach precisely because no
   equiphi29 band exists. The non-transfer in §5 now has three targets behind
   it, which argues against scaling a band to a new polymerase at all.
3. **What does `--target-coverage 0.70` mean across the change?** A 0.70
   request is a different request: on this pool a panel reading 0.386
   symmetric reads 0.245 directional. And §3 says the panel LIMITS do not all
   move together -- a `max_host_coverage` limit tightens by 2x while a
   coverage target loosens by about 1.5x, so the two cannot be rescaled by one
   factor either.

## Scope of the claim

One real target against two real hosts, two candidate pools, four panel sizes,
one run per arm, 12-mers, equiphi29 at 43 C, `hybrid` only -- thirty-two runs.
The two pools share their target, so they are not two independent pairs.
Specific caveats:

- **The background axis answers, but quietly.** The delivered panels hold two
  to four chr1 sites, so `bg_coverage` is 1e-4 at most and the halving is the
  arithmetic of a handful of non-overlapping spans. That is a real measurement
  where the first pool had none, and it is not a crowded host: a design whose
  panel binds the host hundreds of times could behave differently, and nothing
  here says it would not.
- **This pool cannot discriminate.** `filter` warned that mean
  matched/mismatched occupancy is 1.23 with 27% of candidates saturated, the
  regime Known Issue 17 describes. Specificity is weakly selected for whatever
  the optimizer does, which is one reason `selectivity_density` sits flat
  across arms.
- **chr1 is 7.4% N.** `count-kmers` warns that sites in those regions are not
  real and coverage is overstated; the denominator includes them in every arm,
  so the comparison is unaffected and the absolute `bg_coverage` is a little
  low.
- **Nothing here says which panel amplifies better.** Both arms are scored
  under their own geometry. `calibrate-reach --bam` fits the reach from
  measured depth and still has never been run: there is no BAM or CRAM in this
  repository. No wet-lab outcome is involved in any figure above.
