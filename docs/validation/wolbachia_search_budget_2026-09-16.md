# The search budget changes the panel

Measured 2026-09-16. The second half of plan step 249 of
`docs/superpowers/plans/2026-09-15-condition-aware-pool-design.md`, which asks
for designs compared at equal documented search budgets and for runs allowed to
spend more time.

Reproduce with `scripts/benchmarking/wolbachia_search_budget.py`. Raw numbers in
`examples/wolbachia_pool_design/search_budget_2026-09-16/results.json`.

**A larger budget buys coverage and spends selectivity.** It is not an accuracy
knob. It moves the delivered panel along the trade-off in the direction the
objective points, and on this design it moves it a long way.

## What was run

wMel `NC_002978.6` against the Drosophila `GCF_000001215.4` background, counted
and filtered once under `candidate_retention="post_gini"`, then `score`, then
`plan-pool` three times over the same 2,000-candidate shortlist. Panel sizes 8
to 16, coverage metric occupancy-weighted, selectivity density floor 1.0, which
nothing violates, so every panel is feasible and the objective's deciding term
is coverage.

Only `--swap-max-evaluations` differs between the three runs.

## Result

| Panel size | Coverage at 1,000 | at 10,000 | at 100,000 | Density at 1,000 | at 10,000 | at 100,000 |
|---|---|---|---|---|---|---|
| 8 | 0.6837 | 0.6894 | 0.6912 | 32.38 | 30.25 | 28.79 |
| 10 | 0.7000 | 0.7031 | 0.7138 | 34.34 | 32.62 | 22.65 |
| 12 | 0.7180 | 0.7200 | 0.7358 | 32.25 | 28.78 | 20.36 |
| 14 | 0.7216 | 0.7227 | 0.7468 | 31.68 | 28.35 | 24.21 |
| 16 | 0.7438 | 0.7449 | 0.7642 | 26.93 | 24.45 | 21.86 |

Coverage rises with the budget at every one of the nine sizes. Selectivity
density falls at every one of them. Sizes 9, 11, 13 and 15 are omitted from the
table only for width and behave the same way.

The panels are substantially different, not reorderings of one answer. Jaccard
against the smallest budget:

| Panel size | 10,000 | 100,000 |
|---|---|---|
| 8 | 0.78 | 0.60 |
| 10 | 0.67 | 0.43 |
| 12 | 0.85 | 0.41 |
| 14 | 0.87 | 0.65 |
| 16 | 0.88 | 0.68 |

At size 11 the largest budget shares 0.38 of its panel with the smallest.

## The budget is nearly free in time

| Budget | Wall clock, seconds |
|---|---|
| 1,000 | 19.9 |
| 10,000 | 19.4 |
| 100,000 | 23.0 |

A hundredfold budget costs 15 percent more wall clock. That is a consequence of
the focused evaluator: before that change each panel evaluation cost about 6.8
ms and the same sweep would have been several minutes
(`docs/validation/pool_plan_profile_2026-09-16.md`). The 1,000 and 10,000 rows
are within each other's noise, so the ordering between those two means nothing;
only the 100,000 row is separated.

## What follows

**Comparisons must fix the budget.** Two designs run at different budgets differ
by more than half their primers here, which is larger than most effects this
tool is used to measure. A benchmark that varies anything else while the budget
also moves is not measuring the thing it names.

**A loose specificity floor lets the budget spend specificity.** The floor here
is 1.0 and no panel is near it, so the objective's deciding term is coverage all
the way and more search means more coverage at whatever selectivity cost. A run
that cares about specificity should say so with
`--min-selectivity-density`, not by keeping the search short. Stopping the search
early is a way of not finding the panel the objective was asked for, and it is
not the same as asking for a different panel.

**More search is not a better design.** Every one of these panels is a correct
answer to the question as posed. Whether 0.7642 coverage at density 21.9 beats
0.7438 at 26.9 is a question about the experiment, and nothing in this tool
answers it.

## Limitations

- One organism pair, one length, one chemistry, one machine, single runs. The
  wall-clock figures carry no spread and two of the three are within noise.
- Three budget points, not a curve.
- The density floor of 1.0 is deliberately slack. At a floor that bites, the
  objective's first term is the constraint and the direction of these results
  could differ. That case is not measured here.
- The repair path did not fire in any of these runs, because no panel violated
  the floor. This measures the search, not the repair.
- Coverage and selectivity density are modelled site geometry, not measured
  amplification. See `docs/validation/evidence_matrix_2026-09-15.md`.
