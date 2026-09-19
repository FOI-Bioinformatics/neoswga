# What Stage 1 selecting on the accepted metric actually buys

Measured 2026-09-19 on the Wolbachia/Drosophila pair, closing Known Issue 16.

`optimize_greedy` has always accepted an `objective`, and supplying it makes
the set cover rank candidates by occupancy-weighted coverage gain instead of by
unweighted coverage bins. Commit `59a4ee3` added it under the heading "The
greedy now chooses on the quantity the design is judged on". **No production
caller passed it.** `hybrid_optimizer.py`, `dominating_set_adapter.py` and
`primer_expansion.py` all omitted it; the only caller that supplied one was a
test.

It is now wired, reachable from `params.json` as `stage1_objective_width`, and
**off by default** -- on measurement, not caution.

## Why it was never wired

The scan recomputes the objective for every candidate at every pick. One
`compute_metrics` call on this design costs **36 ms**, so scoring all 2,000
shortlisted candidates at each of Stage 1's picks is 14.4 minutes of objective
evaluation, and the 20,670-candidate inventory is about 2.5 hours -- against
39 s for the whole `optimize` run. A correct rule nobody can afford to run is
not a fix.

So the objective arrives with a prescreen, the shape this codebase already uses
for the swap repair: the cheap bin gain ranks every candidate, and only the top
`stage1_objective_width` are scored. The prescreen is not a new criterion. It
is exactly what the greedy used when it had no objective at all.

## What turning it on does

`neoswga optimize`, phi29, `max_primer` 2000, seed 20260919, width 64. Stage 1
selects 20 primers for a 12-primer panel, so the width bounds it to 1,280
objective evaluations per run rather than 40,000.

| n | Jaccard | effective coverage | selectivity density | host sites | seconds |
|---|---|---|---|---|---|
| 6 | 0.500 | 0.6452 -> 0.6525 (+0.0073) | 33.59 -> 31.96 (-1.64) | 92 -> 110 | 22 -> 186 |
| 12 | 0.200 | 0.7283 -> 0.7422 (+0.0139) | 25.62 -> 19.11 (-6.50) | 149 -> 193 | 39 -> 278 |
| 24 | 0.116 | 0.7748 -> 0.8331 (+0.0583) | 19.46 -> 11.87 (-7.60) | 261 -> 456 | 213 -> 742 |

It does exactly what it claims: the metric it now selects on improves at every
size, and by more the larger the panel. It also costs specificity at every
size, and host binding nearly doubles at n=24. Runtime is 3.5x to 8.4x.

**That is a trade, not an improvement**, so the shipped default is unchanged
and a user who wants coverage asks for it. This is the resolution Known Issue
11 reached for `DEFAULT_REDUNDANCY_THRESHOLD`, for the same reason.

With no width configured the delivered panel is identical to the previous
release -- verified at n=12, same primers, same coverage and density to the
last digit, 43 s against 39 s.

## Why it costs specificity

Two things, and the second is the interesting one.

**The objective is coverage, not the acceptance criterion.**
`PoolObjective.coverage()` returns `effective_fg_coverage` and nothing else.
`_objective_gain` is a coverage delta. Constraints reach Stage 1 only as a
STOP rule (`_should_stop_extending`, which halts growth once a panel is past a
background-site cap), never as a selection criterion. So "select on the
quantity the design is judged on" is half true: it changes what COVERAGE means,
from bins to occupancy-weighted, and leaves specificity out of the ranking
exactly as before. Specificity enters only through the tie-break on
`total_bg_sites`, which fires on exact ties.

**And occupancy weighting actively favours host binders.** Occupancy is high
for a primer whose melting temperature sits near the reaction temperature, and
that is the same property that makes it bind the host well. It is the axis
Known Issue 17 describes: occupancy and mismatch discrimination move in
opposite directions. Weighting coverage by occupancy therefore promotes exactly
the candidates a specificity-led design would demote. The unweighted bin count
was accidentally more specific.

This is the same shape as
[frontier_refill_on_optimize_2026-09-19.md](frontier_refill_on_optimize_2026-09-19.md):
giving a coverage-only objective more power -- a wider candidate universe
there, a better coverage metric here -- buys coverage and spends specificity,
because nothing in the selection resists it.

## When to turn it on

Set `stage1_objective_width` in `params.json` when coverage is the binding goal
and specificity is held by a configured panel limit rather than by the search.
64 is the tested value; None is the default and keeps the bin count. An
unbounded scan is correct and unaffordable on a real pool.

## What this does not establish

One organism pair, one polymerase, three panel sizes, one width. It does not
show what the option is worth under equiphi29 or an additive, where the
occupancy spread across the admitted pool is 7.8 and 8.3 rather than phi29's
1.8 -- so the effect should be LARGER there, in both directions, and that is an
expectation rather than a measurement.

It also does not test what a genuinely specificity-aware Stage 1 would do. That
remains the open question, and CLAUDE.md's "Stage 1 is deliberately NOT
constraint-aware" records three attempts that failed to improve a delivered
panel.
