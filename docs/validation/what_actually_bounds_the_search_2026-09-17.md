# What actually bounds the objective-scored search

Measured 2026-09-17 on the real Wolbachia pair, while checking two premises of
Phase 4 increment 4 of the plan for `docs/validation/pipeline_audit_2026-09-16/`.

Both premises were wrong, in the same direction: the objective is used far less
than the cost projections assume, and the knob that appears to control how much
it is used does not.

Reproduce with `scripts/benchmarking/count_objective_calls.py`.

## The greedy never scores with the objective

`docs/validation/objective_evaluation_cost_2026-09-17.md` prices scoring every
candidate with the full objective at every greedy step, and concludes 48.9 h for
a 24-primer panel over 491,836 candidates. The arithmetic is right and the
architectural conclusion follows from it, but no user can run that scan.

`optimize_greedy` takes `objective=None`
(`dominating_set_optimizer.py`) and threads it inward. Every `objective=`
in the package is either that default, the internal forwarding of it, or one of
two real call sites: `pool_planner` into `refine_by_swaps`, and
`swap_refinement` into its own scorer. Nothing supplies one to the greedy. The
greedy scores on coverage bins.

So the objective-scored searches that exist are the swap repair and the beam,
and both are budgeted already.

## The objective is barely used when no constraint binds

One `plan-pool` size row at panel size 12, occupancy-weighted coverage,
selectivity density floor 1.0:

| Selectivity floor | Objective evaluations | Wall clock | Repair |
|---|---|---|---|
| 1.0 | 1 | 2 to 6 s | not attempted |
| 40.0 | about 3,200 | 22 to 25 s | fires, stops on time |

Counted as cache misses on the objective's evaluator, which is the quantity that
costs anything.

The floor of 1.0 is the one
[the published search-budget sweep](wolbachia_search_budget_2026-09-16.md) used,
and that note says so: nothing violates it. With nothing to repair, the
objective is called once, for the final assessment of the delivered row. The
23 s that note measures for nine sizes is therefore almost entirely the
optimizer, not the objective.

## The evaluation budget is not the bound

With a floor that does bind, the repair runs. Changing only
`swap_max_evaluations`:

| Evaluation budget | Objective evaluations | Repair evaluations | Stop reason |
|---|---|---|---|
| 10,000 | 3,292 | 8,944 | `time_limit` |
| 100,000 | 3,212 | 7,391 | `time_limit` |

Both stop on the clock. `OptimizerConfig.swap_max_seconds` defaults to 10.0 s
and has no CLI flag, while `--swap-max-evaluations` is exposed on `plan-pool`.
So on this pair the exposed knob is inert above roughly 8,000 and the real bound
is a default the user cannot set.

That is the class Known Issue 8 tracks, in yet another shape: not a flag nobody
reads, but a flag that is read and then overruled by a second limit on the same
loop. It also means the 10,000 against 100,000 columns of the search-budget
table are not a budget comparison. Whatever separates them, it is not the
evaluation budget, because neither run reached it.

## The dimer guard rejects most pairs before they are scored

`refine_by_swaps` increments `evaluations` before the dimer guard and scores
only what survives it. The repair counted 7,391 to 8,944 evaluations against
about 3,200 objective calls in the same run, so roughly two thirds of candidate
pairs are rejected without being scored.

This answers an open question in
[the parallelism audit](parallelism_opportunities_2026-09-17.md), which offered
two explanations for the swap loop's wall clock and could not choose between
them. It is the guard, not early convergence.

## What this means for increment 4

The increment is still worth doing and its shape is unchanged: a cheap prescreen
over the frontier, then the full objective over a bounded number of leaders. Two
things change.

The scan it bounds is the swap and beam scan, not a greedy scan, because the
greedy does not score with the objective. Bounding a greedy scan that does not
exist would be building the thing the audit is about.

And the budget the increment introduces has to be one budget. Today two limits
govern the same loop, one exposed and one not, and the exposed one loses. The
plan already asks for `stop_detail` to distinguish `panel_evaluations` from
`seconds`, which is exactly the distinction that is missing here; the point of
this note is that the distinction is not hypothetical and that the current
reporting hides it behind a single `budget_exhausted`.

## What this does not establish

- One pair, one panel size, one machine, single runs. The 3,200 figure will move
  with the target, the panel size and the floor chosen.
- It does not measure the beam separately from the swaps. Both run inside the
  repair and both were counted together.
- It does not explain what separates the 10,000 and 100,000 columns of the
  search-budget table. It establishes only that the evaluation budget is not the
  cause, since neither run reached it. Timing jitter under a 10 s deadline is a
  candidate and is not measured.
- A floor of 40.0 was chosen because it binds on this pair, not because it is a
  reasonable design constraint.
- Counting cache misses understates the loop's work and is meant to: a cache hit
  costs a dictionary lookup, and the question here is cost.
