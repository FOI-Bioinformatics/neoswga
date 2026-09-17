# Positions for a frontier that moves, and the panel it does not change

Measured 2026-09-17. Phase 4 increment 3 of the plan for
`docs/validation/pipeline_audit_2026-09-16/`.

The increment does two things: it lets the position cache admit and drop
candidates after construction, and it makes the check that refuses an
unmeasured candidate actually run. Neither is supposed to move a delivered
panel, and this note records the run that confirms it.

## The acceptance check

wMel `NC_002978.6` against the Drosophila `GCF_000001215.4` background, counted
and filtered once under `candidate_retention="post_gini"`, then `score`, then
`plan-pool` at panel sizes 8 to 16, occupancy-weighted coverage, selectivity
density floor 1.0, `--swap-max-evaluations 10000`. The same design was run on
this branch and on `main` over the same index.

| Panel size | Delivered panels identical | Jaccard | Coverage |
|---|---|---|---|
| 8 to 16, every size | yes | 1.000 | identical to 8 decimal places |

The recommendations are identical too, and both runs read the inventory rather
than falling back to the CSV: universe 20,670, frontier 2,000, examined 2,000,
unexamined 18,670.

The rebuilt pool reproduces the funnel this design is documented with, which is
worth stating because the acceptance would otherwise only compare two runs of
the same possibly-wrong thing.

| Stage | Candidates |
|---|---|
| Total 12-mers | 874,596 |
| After background frequency | 706,756 |
| After thermodynamics (hard QC) | 491,836 |
| After the evenness gate | 20,670 |
| After the `max_primer` cut | 2,000 |

Coverage also matches
[the search-budget table](wolbachia_search_budget_2026-09-16.md) at the same
budget: 0.6894 at size 8, 0.7031 at size 10, 0.7200 at size 12.

## What the check refuses now, and did not before

`CandidateProvider.ensure_positions` asked whether a candidate had a hit on any
prefix. A candidate with foreground sites and no background entry passed that
question, and a design would then read its absent host sites as perfect
specificity. A candidate indexed against a host it binds nowhere failed it,
although that zero is a measurement.

The question is per prefix and about entry presence, not hit count. On the run
above nothing was refused, which is the expected result under `post_gini`:
every eligible candidate is indexed against both references, so the 2,000 in
the frontier all have entries. The refusal is pinned on fixtures instead, in
`tests/test_positions_arrive_on_demand.py`.

It also raises when no cache is attached, where it used to return quietly.
That mattered more than the predicate did, because no production caller
attached one. `pool_planner._prepare_candidate_pool` now attaches the
optimizer's cache and runs the check before any panel is evaluated, and
`tests/test_the_frontier_is_vouched_for_before_it_is_scored.py` pins the wiring
separately from the behaviour.

## What this does not establish

- It does not exercise `load` or `release` at scale. The frontier does not move
  yet; `advance()` still returns False, and expansion arrives in increment 5.
  Both methods are pinned on fixtures only.
- It does not measure memory. The claim that a moving window needs `release`
  rests on the 443 MB index figure from the
  [retention benchmark](wolbachia_retention_benchmark_2026-09-16.md), not on a
  measurement of a window that moves, because there is not one to measure.
- One machine, one pair, one retention mode. An `all_qc` design admits 491,836
  candidates rather than 20,670, and the check has not been run against that
  index.
- Identical panels mean this increment is behaviour-preserving. They say
  nothing about whether the panels are good.
