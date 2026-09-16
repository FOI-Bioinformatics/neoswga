# Where a plan-pool run spends its time

Measured 2026-09-16. Plan step 204 of
`docs/superpowers/plans/2026-09-15-condition-aware-pool-design.md`, which asks
for incremental interval and occupancy updates only where profiling warrants
them. This is the profile, and what it warranted.

Reproduce with `scripts/benchmarking/profile_pool_plan.py <tmp> <density floor>`.

## The workload

A 2 Mb synthetic target and a background derived from it, 300 candidate 12-mers
with 3 to 30 planted sites each, panel sizes 4, 8, 12 and 16, and a selectivity
density floor of 400. The floor matters: at a floor every panel clears, no
repair runs and the profile describes the optimizer rather than the search. At
400 every size needs a repair, which is the path this measurement is about.

Single runs on one machine. The absolute seconds carry no spread; the ratio
between the two configurations is the result.

## Before

`plan_pool` took 29.7 s. `_repair` accounted for 27.3 s of it, and inside that
`compute_metrics` took 27.6 s across 4,040 calls, about 6.8 ms each.

`BaseOptimizer.compute_metrics` answers every question anything has ever asked
of a primer set: both coverage figures, a coverage sweep across five reaches,
background coverage, gap mean, maximum, Gini and Shannon entropy, melting
temperatures, dimer risk and strand alternation. That is the right shape for a
report, which is computed once.

`PoolObjective` and the plan row between them read five fields: geometric
coverage, occupancy-weighted coverage, selectivity density, total background
sites and the maximum gap. Everything else in each of those 4,040 calls was
discarded. The coverage sweep alone ran `_compute_coverage_at` 24,240 times,
six per call, for a figure nothing on this path reads.

## What changed

`neoswga/core/pool_metrics.py` computes those five fields, through the same
helpers, and nothing else. `BaseOptimizer.compute_pool_metrics` exposes it, and
`plan_pool` uses it when the optimizer has it. A caller passing its own
optimizer-shaped object keeps the full evaluation.

`tests/test_pool_metrics_agree_with_the_full_evaluation.py` pins the two to the
same numbers on randomised panels, with and without reaction conditions, so the
occupancy path is covered rather than assumed. Two definitions of one quantity
drift while each stays self-consistent, and that test is what stops it.

## After

| | Before | After |
|---|---|---|
| `plan_pool` wall clock, seconds | 29.7 | 5.8 |
| Largest item in the profile | `compute_metrics`, 27.6 s | `_union_coverage`, 2.5 s |

The delivered panels are the same. Each size needed the same number of repair
swaps, 2, 3, 4 and 2, and reported the same selectivity density.

One caveat on "the same". The swap loop carries a wall-clock budget as well as
an evaluation budget, so where the time limit binds a faster evaluator fits more
evaluations into it and can reach a different panel. On this workload one size
moved from 7,646 evaluations to 8,352 without changing its result.

## Is an incremental interval update warranted now

Not yet, on this evidence.

`_union_coverage` is now the largest single item at 2.5 s of 5.8 s. An
incremental version is the obvious next step and it is also the riskiest: a swap
removes a primer as well as adding one, so removing its contribution from a
union needs per-base multiplicity rather than the union alone, and an off-by-one
there produces a coverage figure that is wrong without being obviously wrong.

The ceiling on the gain is 2.5 s out of 5.8 s on a path that now completes in
seconds, against a correctness risk in the quantity the whole design is accepted
on. Removing work that was being discarded bought 5.1x for no change in what is
computed. Taking the remaining 1.8x by changing how it is computed is a
different trade, and the measurement does not yet support it.

If it is taken later, step 204's own instruction is the right one: compare each
incremental marginal and swap against a full evaluation on randomised small
cases, the way `test_pool_metrics_agree_with_the_full_evaluation.py` does here.

## Limitations

- One synthetic organism pair, one length, one density floor, one machine.
  Single runs.
- Synthetic sequence has a uniform base composition and planted sites, so site
  clustering is unlike a real genome. It exercises the code paths, not the
  distribution.
- The profile measures `plan_pool`. It says nothing about `optimize`, which
  reaches `compute_metrics` by a different route and still needs the full one.
