# The coverage geometry: the default stays, and both figures are now reported

Decided and changed 29 October 2026, on `main` at `366df9b` plus the change
this document records. Phase 5.2 of
[the plan](../reports/PLAN_2026-10-28_pool_design_trustworthiness.md), which
required a measurement before any fix and **did not pre-approve one**.

The measurement is
[`2026-10-29-directional-coverage-measured.md`](2026-10-29-directional-coverage-measured.md):
the shipped symmetric window overstates the directional figure by **1.27x to
1.85x** on three real panels at the shipped 3 kb reach, 19.5 to 25.8 percentage
points, with the worst absolute error at the intermediate site density where a
real SWGA design sits.

The regression is `tests/test_both_coverage_geometries_are_reported.py`
(18 tests).

## The decision

**`symmetric` stays the default.** Three reasons, in order of weight:

1. **The reach was fitted against it.** `reach_calibration.md` enters the
   published 60-75% breadth band at 2,866 bp and leaves it at 4,551 bp under
   symmetric windows; under directional windows the same five sets enter at
   4,354 and leave at 6,710
   ([`2026-09-23-reach-refit-directional.md`](2026-09-23-reach-refit-directional.md)).
   The shipped phi29 3,000 sits just inside the symmetric band and well below
   the directional one, where `TmL/Even` reads 0.48 against the 60-75% its
   wet-lab outcome recorded. Switching the geometry alone does not produce a
   more correct model; it produces a model mis-calibrated by construction.
2. **Every recorded figure in this repository was produced under it.**
   `fg_coverage` is authoritative in every `step4_improved_df_summary.json` and
   quoted throughout `docs/validation/`. Changing the default moves all of
   them, and re-baselining a research record is a maintainer's decision, not a
   repair.
3. **A correct fix changes selection, not only reporting.** Stage 1's set cover
   would have to learn that a site reaches `r` and not `2r`;
   [`2026-09-23-does-the-panel-move.md`](2026-09-23-does-the-panel-move.md)
   measured that the delivered panel then differs at two of four sizes at equal
   reach and three of four at the refitted one. That is a design change with
   its own measurement to make, and the directional column of the 5.1 table is
   explicitly NOT what such a run would deliver, because those panels were
   selected under the symmetric model.

What the plan listed as prerequisites for acting already existed: a geometry
key (`coverage_geometry`, both values implemented and schema'd), a refitted
reach, and a measurement of whether the panel moves. The missing piece was not
a capability. It was that a reader of a summary could not tell which convention
produced the number in front of them.

## What changed instead

### 1. Both figures are reported for the delivered panel, labelled

`core/coverage_geometries.py` computes the panel's coverage under both windows
once, after the search, and `step4_improved_df_summary.json` carries it:

```json
"coverage_geometries": {
  "reach": 4000,
  "reported_geometry": "symmetric",
  "symmetric": 0.4749259842519685,
  "directional": 0.26096929133858265,
  "inflation": 1.8198539062429286,
  "note": "The delivered panel was selected and scored under 'symmetric'. ..."
}
```

and `optimize` logs one line, as a hybrid run already prints two coverage
figures:

```
Coverage geometry at 4,000 bp: 0.4749 (symmetric) against 0.2610 (directional),
1.82x. Scored on 'symmetric'.
```

That is a real run: equiphi29 at 43 C on the synthetic 1.27 Mb target, twelve
primers, at equiphi29's 4,000 bp realistic reach. `symmetric` equals
`fg_coverage` to the last digit, which is the point --
`test_the_symmetric_column_equals_the_optimizers_own_coverage` asserts it
against the optimizer's own `compute_metrics` rather than against a
re-derivation, because three implementations of one quantity is how this
codebase has produced disagreeing coverage numbers before.

**Once, after the search, not per candidate.** `_compute_metrics` runs for every
candidate panel and `_oriented_positions` fetches the second orientation only
when the configured geometry asks for it, because `get_positions` is the
hottest method in the package. Paying for both geometries on every evaluation,
to produce a figure nothing selects on, is not a trade worth making.

**The second column bounds the overstatement, not the loss.** The panel was
selected under the geometry it was scored on, so the other figure is the same
sites under the other window. A set cover taught the directional geometry would
pick differently and do better than that column. The `note` field says so, in
the artifact, where a reader of the number will be.

### 2. A directional run must name its own reach

`coverage_geometry: "directional"` with no `coverage_reach` was reachable from
params.json with one key and silently used the polymerase default -- a reach
fitted under the other geometry. `design_request._resolve_reach` now refuses
it, naming the calibrated bands and the record.

**Nothing is substituted, deliberately.** Substituting 4,400 is the tempting
move and it cannot be made honestly: the directional band was fitted for phi29
against one published outcome, and no directional band exists for `bst`,
`equiphi29` or `klenow`, so a default would invent three numbers from one
measurement. The same reasoning refuses a design that sets glycerol because
nothing computes its Tm effect.

`_resolve_geometry` also refuses an unrecognised geometry up front.
`coverage.site_spans` already raised on one, but mid-search, after the counting
and the filtering had been paid for; the schema's enum only reaches a caller
that validates against the schema.

### 3. The geometry is on the recorded request

`DesignRequest.coverage_geometry`, so the request hash moves with it and a run
under one geometry can be told from a run under the other.
`default_sources["coverage_geometry"]` distinguishes a request that named it
from one that took the default.

## What this does not do

- **It does not fix the defect.** The symmetric window still credits every
  occurrence in the direction it cannot extend in, and `fg_coverage` is still
  that figure. This makes the overstatement visible and quantified in every
  artifact; it does not remove it.
- **It does not refit anything.** No reach, no threshold and no weight moves.
- **It does not change any delivered panel.** Selection, scoring and every
  panel limit are untouched; `normalized_score` is unchanged. A re-run of any
  recorded design produces the same primers and the same `fg_coverage`.
- **It does not establish the inflation on a real genome.** The 5.1 measurement
  is on synthetic sequence, where sites are Poisson-distributed and real
  genomes cluster. The direction and the mechanism generalise; the percentages
  are for that site geometry.

## What remains for whoever decides to switch

In the order the switch needs them, all of them measurements rather than
opinions:

1. Teach Stage 1's set cover the directional geometry, and measure the panel it
   then delivers against the symmetric one at each geometry's own calibrated
   reach. `2026-09-23-does-the-panel-move.md` is the start of this and covers
   one pool.
2. Re-baseline every `fg_coverage` quoted in `docs/validation/` and say, per
   record, whether its conclusion survives the new figure.
3. Fit a directional band for the other three polymerases, or scope the default
   switch to phi29 and refuse the rest, which is what the reach refusal does
   today.
4. Decide what `--target-coverage` and the panel limits mean across the change,
   since a 0.70 target is a different request under a figure that is 1.3x to
   1.9x smaller.

Nothing above is started here, and the plan did not ask for it.
