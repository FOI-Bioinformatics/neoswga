# How narrow the objective-scored scan can be

Measured 2026-09-17 on the real Wolbachia pair. Phase 4 increment 4's
acceptance measurement, from the plan for
`docs/validation/pipeline_audit_2026-09-16/`.

Reproduce with `scripts/benchmarking/scan_width_sweep.py`.

**A width of 16 reaches the same panel an unbounded scan reaches, using 280
times fewer objective evaluations.** The shipped default is 64, a four-fold
margin on that.

## What was bounded, and why it is not what the plan said

The plan describes bounding the per-step scan of a greedy that scores every
candidate with the objective. That greedy does not exist: `optimize_greedy`
takes `objective=None` and nothing in the package supplies one, so the only
objective-scored searches are the swap repair and the beam. See
[what_actually_bounds_the_search_2026-09-17.md](what_actually_bounds_the_search_2026-09-17.md).

The scan that needed bounding is the repair's, over `candidates x panel`: two
thousand candidates against a twelve-primer panel is twenty-four thousand pairs
per round, and every pair surviving the dimer guard was scored in full.

The prescreen was already present, as the other branch of the same function.
`refine_by_swaps` without an objective ranks swaps by bin gain minus bin loss.
Bounded mode uses exactly that to rank the pairs and scores only the leaders, so
this introduces no new criterion.

## The prescreen is a filter, not a decider

The plan's acceptance was Jaccard against an unbounded scan, with 0.8 at width
64 as the line below which the prescreen would be deciding.

Taken at the default budget the answer looked alarming: Jaccard 0.429 to 0.750.
It was not the prescreen. The unbounded scan stopped at `evaluation_limit` on
every size measured, so it was not a gold standard but a truncated search.

Given a budget it cannot exhaust, at panel size 12:

| Configuration | Objective evaluations | Wall clock | Coverage | Density | Stop | Jaccard vs optimum |
|---|---|---|---|---|---|---|
| unbounded, default budget | 1,198 | 12.6 s | 0.712284 | 40.009 | evaluation limit | 0.500 |
| unbounded, large budget | 17,913 | 85.0 s | 0.719019 | 41.260 | local optimum | 1.000 |
| width 64, large budget | 320 | 45.9 s | 0.719019 | 41.260 | local optimum | 1.000 |
| width 16, large budget | 64 | 44.8 s | 0.719019 | 41.260 | local optimum | 1.000 |

All three converge to the identical panel. Width 16 gets there on 64 objective
evaluations against 17,913, and the truncated run was the odd one out at
Jaccard 0.500 against the optimum it never reached.

A single-round probe on the real pool agreed in advance: ranking 9,456 swap
pairs by bin gain and by the objective put the objective's ten best pairs at bin
positions 0 through 9. That probe covered only the feasible regime, since all
9,456 pairs satisfied the constraint, which is why the end-to-end comparison
above is the result and the probe is only corroboration.

## A width also fixes what the default budget did

Without a width, the default evaluation budget truncated every size measured.
With one it converges, so the answer is a local optimum rather than wherever the
budget ran out. At width 64 the repair spent 512 to 768 objective evaluations
against a budget of 10,000.

That is the stronger reason to ship a width than the speed is.

## Where nothing binds, nothing changes

At a selectivity floor of 1.0, which
[the published sweep](wolbachia_search_budget_2026-09-16.md) used and which
that note says nothing violates, the repair is never attempted. Every width
from 16 to unbounded returns Jaccard 1.000 with identical coverage at sizes 10,
12 and 14. The shipped configuration is unaffected.

## Where a constraint binds, at the default budget

| Width | Size | Jaccard vs untruncated | Coverage | Density | Objective evaluations | Pairs seen |
|---|---|---|---|---|---|---|
| 16 | 12 | see note | 0.711390 | 41.615 | 128 | 190,848 |
| 64 | 12 | see note | 0.711390 | 41.615 | 512 | 190,848 |
| 256 | 12 | see note | 0.711390 | 41.615 | 2,048 | 190,848 |
| unbounded | 12 | 0.500 | 0.712284 | 40.009 | 1,198 | 10,000 |

Every bounded width reached `local_optimum`; the unbounded run reached
`evaluation_limit`. The pairs column is the point: a bounded run examines 190,848
pairs of the pool where the unbounded one examined 10,000, because the cheap pass
is not charged against the evaluation budget. Charging it would rank a prefix of
the pool and reintroduce the blindness the bound removes.

## The cost that remains

Wall clock did not fall by the factor the objective evaluations did, and at the
default budget a bounded run is slower: 26.6 s at width 16 against 22.3 s
unbounded over three sizes. The unbounded run is quick because it gives up.

With the objective no longer dominant, the remaining cost is the prescreen's own
pass: for each round it walks every admissible pair, and the per-pair work is a
dimer lookup and two set operations. At 190,000 to 330,000 pairs per round over
several rounds that is now the largest term. Reducing it is a separate change
and is not attempted here.

## What this does not establish

- One pool, one pair, one panel size for the untruncated comparison, one
  machine, single runs. The width that suffices will move with the pool and the
  constraint.
- The floor of 40.0 was chosen because it binds on this pair, not because it is
  a sensible design constraint.
- Width 16 sufficed here. 64 is a margin chosen because its extra cost is
  negligible beside the prescreen's pass, not because any measurement showed 16
  to be insufficient.
- Nothing here bounds the beam, which is the other objective-scored search. It
  carries its own evaluation budget and was not changed.
- The prescreen's agreement with the objective was measured in the feasible
  regime only. A pool where the constraint separates candidates sharply could
  rank differently, and the bin gain knows nothing about background load.
- Converging on a local optimum is not finding a global one. Every figure here
  is a heuristic search result, and modelled site geometry rather than measured
  amplification. See
  [evidence_matrix_2026-09-15.md](evidence_matrix_2026-09-15.md).
