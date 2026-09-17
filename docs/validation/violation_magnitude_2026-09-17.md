# Ranking failing panels by how far they missed

Measured 2026-09-17 on the real Wolbachia pair. Follows the defect
[increment 5's measurement](frontier_refill_2026-09-17.md) exposed.

Both objective-scored searches ranked on `len(violations)`. Two panels failing
the same single constraint therefore tied, coverage broke the tie, and the
deciding metric was free to drift. `refine_by_swaps` already stated the intended
rule -- "before a panel is feasible the useful direction is out of violation
rather than up the coverage curve" -- and counting implements that only across
different numbers of violated constraints, never within one.

`PoolObjective.shortfall` measures the distance instead. Each term is relative
to its own limit, so a density floor and a background-site ceiling are
comparable and neither dominates by being measured on a larger scale, and the
terms are summed so failing two limits is worse than failing one. It is zero
exactly when `violations` is empty, which is what keeps every feasible panel
ahead of every infeasible one.

## What changed on the real pool

Panel size 12 at a selectivity density floor of 100, which no panel from this
inventory reaches, so every row is rejected and the ranking among failing panels
is all that moves.

| Refills | Density before | Density after | Coverage before | Coverage after |
|---|---|---|---|---|
| 0 | 20.852 | 28.776 | 0.739548 | 0.720000 |
| 4 | 14.589 | 19.162 | 0.765454 | 0.741632 |

The deciding metric improved at both refill levels, by 7.9 and 4.6 points of
density, at a cost of 2.0 and 2.4 points of coverage. That is the trade the rule
asks for: on a rejected row the coverage figure is informational and the
distance from the constraint is the thing a reader wants.

## What makes the ordering safe

Ranking on distance is right while feasibility is reachable. On a constraint no
panel can meet, a search that chases it will trade real coverage for a step
toward a floor it never reaches. Two existing guards in
`tests/test_pool_plan_repair.py` caught exactly that, and one of them is named
for it: a failed repair must not deliver something worse than it started with.

So a repair that does not succeed now returns the panel it was given rather than
the one the search wandered to. The row is rejected either way, so a partial move
buys nothing and can cost. Both guards then pass unchanged, which is the sign
the two changes belong together: the ordering pursues feasibility and failing to
reach it costs nothing.

On the real pool that accounts for most of the improvement at zero refills,
where the reported panel is now the optimizer's rather than the swap loop's.

## What it did not fix

Density still falls as the frontier is refilled, 28.8 to 19.2. That drift is no
longer the repair's: with a failed repair returning its input, the panel a
rejected row reports is the optimizer's own selection, and the optimizer finds
higher coverage when handed more candidates. Stage 1 does not see the objective
at all.

So the defect has moved rather than gone, and it has moved somewhere narrower
and better understood. Fixing it means making the optimizer's own selection
constraint-aware under an unreachable limit, which is a larger change than this
one and wants its own measurement.

## What this does not establish

- One pool, one pair, one panel size, single runs. The floor of 100 was chosen
  because nothing reaches it, which is what isolates the ranking among failing
  panels; it is not a sensible design constraint.
- It does not show the change helps when a constraint IS reachable. There the
  rule should let a search find a qualifying panel it previously missed, and no
  such case was constructed on real data. The unit tests cover the ordering;
  the benefit is inferred.
- The relative scaling of the two constraint terms is a choice. Missing a floor
  by half scores the same as exceeding a ceiling by half, which makes them
  comparable but does not make them equally important, and no measurement
  supports treating them as equal.
- An unmeasurable coverage scores infinite, which is a judgement that such a
  panel is not nearly-feasible rather than a measured distance.
- Coverage and selectivity remain modelled site geometry, not measured
  amplification. See [evidence_matrix_2026-09-15.md](evidence_matrix_2026-09-15.md).
