# Does the delivered panel move under the directional geometry?

Measured 23 September 2026, once the geometry reached every path a run uses:
Stage-1 set cover (which already selected directionally), the Stage-2 swap
through `effective_fg_coverage`, and the reported `fg_coverage` and
`bg_coverage`.

**Yes, at some panel sizes, under both reach settings tried.** So the geometry
is a design change on this pool and not only a reporting one.

## Setup

*Prevotella melaninogenica* (3,168,282 bp) against human chr21, 12-mers, the
shipped `equiphi29_12mer` configuration and its `step3_df.csv`, at four panel
sizes. Two runs per size per configuration through `neoswga optimize`, differing
only in `coverage_geometry` and `coverage_reach`.

## At equal reach: coverage halves, and the panel moves at two sizes of four

Both at `coverage_reach` 4000, the configuration's own value.

| n | symmetric | directional | ratio | Jaccard |
|---|---|---|---|---|
| 4 | 0.05303 | 0.02646 | 0.499 | **0.600** |
| 8 | 0.09012 | 0.04545 | 0.504 | 1.000 |
| 16 | 0.14950 | 0.07511 | 0.502 | **0.882** |
| 24 | 0.18558 | 0.10027 | 0.540 | 1.000 |

The ratio is close to 0.5 throughout, which is what a half-width window gives
when the spans mostly do not overlap. The panel is identical at n=8 and n=24
and differs at n=4 and n=16.

## At the refitted reach: coverage is 18-25% lower and the panel moves at three sizes of four

`docs/validation/2026-09-23-reach-refit-directional.md` fits the directional
band at 1.52x the symmetric one on Wolbachia. Applying that factor here gives
6080 against 4000.

| n | symmetric @4000 | directional @6080 | delta | Jaccard |
|---|---|---|---|---|
| 4 | 0.05303 | 0.04030 | -24.0% | 1.000 |
| 8 | 0.09012 | 0.07275 | -19.3% | **0.600** |
| 16 | 0.14950 | 0.11229 | -24.9% | **0.600** |
| 24 | 0.18558 | 0.15117 | -18.5% | **0.548** |

### The refit factor does not transfer between targets

This is the finding worth keeping. On Wolbachia the refitted pairing is nearly
coverage-neutral -- `TmL/Even` reads 0.619 symmetric at 3000 and 0.609
directional at 4400, a difference of 1.6%. On Prevotella the same 1.52x factor
leaves coverage **18 to 25% lower**.

The factor was fitted by asking which reach puts one Wolbachia set inside its
measured breadth band. Coverage is concave in reach and the concavity depends on
site density, so a ratio of band edges on one target does not preserve coverage
on another. Nothing in
`docs/validation/2026-09-23-reach-refit-directional.md` claimed it would, and
this says explicitly that it does not.

A reach refitted on Prevotella would land somewhere else again. That fit is not
possible here: there is no published wet-lab outcome for a Prevotella set in
this repository to fit against, which is the whole reason the Wolbachia sets
were used.

## What this means for the default

The plan for the next increment fixed its decision rule before the measurement:

> **Panels differ:** the default stays `symmetric`. Record the Jaccard and the
> direction, and state plainly that the geometry is a design change on this
> pool.

Panels differ, at two of four sizes at equal reach and three of four at the
refitted one. **The default stays `symmetric`.**

That is the same resolution Known Issues 11, 16 and 17 reached: the capability
ships reachable and off, with the measurement recorded, because no delivered
panel has been shown to improve. What would change it is an outcome to fit
against -- a measured sequencing breadth for a design from this pool -- and
`calibrate-reach --bam` exists for exactly that and has never been run, there
being no BAM in this repository.

## What this does not establish

- **One target, one background, one candidate pool.** Four panel sizes on
  Prevotella against chr21.
- **The background axis was not exercised.** `bg_coverage` reads exactly
  0.00000000 in every run above and `selectivity_density` is the
  `MAX_SELECTIVITY` sentinel, because this candidate pool has zero exact
  matches anywhere in chr21 by construction -- `reach_calibration.md` records
  it as 219,730 12-mers chosen that way. So the geometry's effect on host
  binding is untested here, and `bg_coverage` is the one quantity that sees
  background site POSITION.
- **Which panel is better.** Both are scored under their own geometry, and
  nothing here says which figure predicts amplification. The panels differ; no
  claim is made that one is preferable.
- **Whether the differing panels differ materially.** Jaccard 0.548 at n=24
  means eleven oligos of twenty-four are shared. Whether the two sets would
  behave differently in a reaction is exactly the question no measurement in
  this repository can answer.
