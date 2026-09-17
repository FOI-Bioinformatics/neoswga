# What refilling the frontier buys, and what it exposes

Measured 2026-09-17 on the real Wolbachia pair. Phase 4 increment 5 of the plan
for `docs/validation/pipeline_audit_2026-09-16/`.

`advance()` returned False from the day it was written, so the retained
inventory was decorative: 20,670 candidates eligible, 2,000 searched, 18,670
unable to affect any panel. It now widens the frontier, `plan_pool` refills when
a size row cannot satisfy its constraints, and each refill doubles, so four
reach the whole universe.

## Where the shortlist already suffices, nothing happens

Panel size 12, occupancy-weighted coverage, `max_frontier_refills=4`:

| Selectivity floor | Refills used | Feasible | Coverage | Density | Examined |
|---|---|---|---|---|---|
| 40 | 0 | yes | 0.711390 | 41.615 | 2,000 of 20,670 |
| 60 | 0 | yes | 0.653493 | 60.112 | 2,000 of 20,670 |

A row that qualifies never refills, so the delivered panel and the runtime are
unchanged wherever the shortlist is enough. Both of these are identical to the
same run with refilling disabled.

## Where it does not, the refill turns a shrug into an answer

At a floor of 100, which the shortlist cannot reach:

| Refills allowed | Refills used | Ended | Examined | Feasible | Coverage | Density | Host sites | Seconds |
|---|---|---|---|---|---|---|---|---|
| 0 | 0 | frontier exhausted | 2,000 | no | 0.739548 | 20.852 | 151 | 13 |
| 4 | 4 | inventory exhausted | 20,670 | no | 0.765454 | 14.589 | 211 | 64 |

The panel still does not qualify, and that is the point. Before, the run
reported `frontier_exhausted` having looked at under 10% of what it was allowed
to reach, which tells a user nothing about whether the constraint is achievable.
Now it reports `inventory_exhausted` after examining all 20,670, which is a
meaningful negative: no panel of twelve from this inventory reaches a density of
100.

The cost is proportional. Searching ten times the candidates took five times as
long, and building the position cache over the whole eligible set took 3 s
rather than being avoided.

## What it exposes, which is not its own fault

Refilling moved the panel **away** from the floor it was trying to reach: the
density fell from 20.9 to 14.6 while coverage rose from 0.740 to 0.765.

Both panels violate exactly one constraint. `refine_by_swaps` scores
lexicographically on `(-violations, coverage, -background)`, and `violations` is
a COUNT. Once two panels fail the same single constraint the count ties,
coverage breaks the tie, and the deciding metric is free to drift. More
candidates give the coverage term more to work with, so the drift gets worse
rather than better.

The swap loop's own docstring states the intended rule: "before a panel is
feasible the useful direction is out of violation rather than up the coverage
curve." Comparing counts implements that only across different numbers of
violated constraints, not within one. Ordering by violation magnitude would
implement it properly, and that is a change to acceptance semantics with its own
measurement rather than something to slip into an increment about expansion.

**Fixed the same day**, once it was measured. `PoolObjective.shortfall` ranks
failing panels by distance from feasibility instead of by how many constraints
they fail, and a repair that does not succeed now returns the panel it was
given. On this configuration the density rose from 20.9 to 28.8 at zero refills
and from 14.6 to 19.2 at four. The numbers in the table above are from before
that change. See
[violation_magnitude_2026-09-17.md](violation_magnitude_2026-09-17.md), which
also records what the fix did NOT address: density still falls across refills,
and that drift belongs to the optimizer's own selection rather than to the
repair.

## What this does not establish

- One pool, one pair, one panel size, single runs, one machine. The floor of 100
  was chosen because the shortlist misses it, not because it is a sensible
  design constraint.
- It does not show that refilling ever produces a feasible panel that the
  shortlist could not. On this pool at these floors it either was not needed or
  was not sufficient. Finding a floor where the extra candidates decide
  feasibility is Phase 5's benchmark, not this note.
- The doubling growth rule is a choice, not a measurement: four refills reach
  this universe where fixed 2,000-candidate batches would take nine. No
  comparison of delivered panels between growth rules was made.
- Nothing here measures memory. The position cache over all 20,670 eligible
  candidates was built without incident on this pair; an `all_qc` design retains
  491,836 and the
  [retention benchmark](wolbachia_retention_benchmark_2026-09-16.md) puts that
  index at 443 MB for the background alone.
- Coverage and selectivity remain modelled site geometry, not measured
  amplification. See [evidence_matrix_2026-09-15.md](evidence_matrix_2026-09-15.md).
