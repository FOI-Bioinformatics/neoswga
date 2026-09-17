# Occupancy coverage over intervals: what it cost and what it changed

Measured 2026-09-17, on the real Wolbachia pair. Done before Phase 4 increment 4
of the plan for `docs/validation/pipeline_audit_2026-09-16/`, because the
measurements in
[objective_evaluation_cost_2026-09-17.md](objective_evaluation_cost_2026-09-17.md)
and [parallelism_opportunities_2026-09-17.md](parallelism_opportunities_2026-09-17.md)
put this ahead of bounding the scan.

`_compute_effective_coverage` allocated a boolean window and a float32
accumulator the length of the target and made two full passes over both per
primer. Its cost was therefore linear in the genome and in the panel and
independent of how many sites there were, which is the wrong way round: a
1.27 Mb target with twenty sites per primer paid for 1.27 million bases either
way.

Taking logs turns the product over primers into a sum, and a sum has a
difference array. Each primer's merged window contributes its log(1 - theta) at
one edge and removes it at another, so the answer is a walk over the edges.
Those number twice the merged windows rather than once the genome.

## The function

Against the real implementation on the real index, from
`scripts/benchmarking/check_interval_sweep.py`:

| Geometry | Panel | Loop | Sweep | Speedup | Absolute difference |
|---|---|---|---|---|---|
| linear | 12 | 8.69 ms | 0.46 ms | 18.8x | 4.0e-9 |
| linear | 24 | 14.63 ms | 0.68 ms | 21.5x | 6.2e-9 |
| circular | 12 | 8.19 ms | 0.44 ms | 18.6x | 4.0e-9 |
| circular | 24 | 14.44 ms | 0.78 ms | 18.5x | 6.2e-9 |

The residual difference is the old float32 accumulation, not the new
arithmetic. The sweep sums in float64 and exponentiates once, where the loop
multiplied in float32 once per primer, so the new value is the more accurate
one. `tests/test_effective_coverage_is_computed_over_intervals.py` pins
agreement against an independent float64 oracle at 1e-12.

## The design, where no constraint binds

Panel sizes 8 to 16, occupancy-weighted coverage, selectivity density floor
1.0. This is the configuration
[the published sweep](wolbachia_search_budget_2026-09-16.md) used, and nothing
violates that floor.

| Quantity | Result |
|---|---|
| Delivered panels, all nine sizes | identical |
| Largest coverage difference | 2.4e-8 |
| Recommendations | identical |
| CPU time, new | 16.30 s, 16.61 s |
| CPU time, old | 17.79 s, 16.71 s |

Within noise, and that is expected rather than disappointing. With nothing to
repair the objective is evaluated once per size row, so a 19x saving on it has
almost nothing to act on. The point of this configuration is that the panels do
not move.

## The design, where a constraint binds

One size row at panel 12 with the density floor raised to 40.0, which the
optimizer's first answer misses, so the repair runs. Two runs of each, taken
alternately.

| Run | Repair | Stop reason | Beam | Density | Meets 40.0 | Coverage |
|---|---|---|---|---|---|---|
| new, rep 1 | succeeded | evaluation limit | target met | 40.0088 | yes | 0.712284 |
| new, rep 2 | succeeded | evaluation limit | target met | 40.0088 | yes | 0.712284 |
| old, rep 1 | failed | time limit | budget exhausted | 28.6882 | no | 0.724030 |
| old, rep 2 | failed | time limit | budget exhausted | 28.6882 | no | 0.724030 |

CPU time was 12.99 s and 12.85 s for the new implementation against 24.03 s and
24.07 s for the old, a factor of 1.87 end to end.

**The speed is not the result. The result is that the design now completes.**
The old implementation could not finish the repair inside its ten second
deadline, so it returned the panel it had, correctly reported
`selectivity below minimum`, and recommended `not_found`. The new one finishes
the beam inside its evaluation budget, reaches the floor, and recommends the
panel. A user asking for a density of 40 on this pair previously got a truthful
refusal; they now get a design.

The delivered panel differs in that case, as it must: five of twelve primers
change and the Jaccard against the old answer is 0.412. Coverage falls from
0.7240 to 0.7123, which is the trade the constraint asks for and the reason the
old panel did not qualify.

## What changed in the code

`neoswga/core/occupancy_coverage.py` holds the computation.
`BaseOptimizer._compute_effective_coverage` delegates to it and supplies the
reach and geometry from its config, so every caller is unchanged. The
extraction was forced by the module size ratchet, which the plan says to answer
by splitting rather than by raising a budget.

`coverage.merged_window_intervals` expresses a primer's window union as
intervals, beside the `_mark_window` that marks the same union into an array.
`tests/test_effective_coverage_is_computed_over_intervals.py` compares the two
base by base across five reach and length combinations, both geometries, and
forty random layouts each; a one-base error in either fails 31 of them.

## What this does not establish

- One pair, one panel size for the binding case, one machine, two runs each.
  The 1.87x will move with the target length, the panel size and how much of a
  run is objective evaluation.
- The density floor of 40.0 was chosen because it binds on this pair, not
  because it is a sensible design constraint. That it now succeeds says the
  search got further in the same time, not that 40 is achievable in general.
- Neither implementation confines a coverage window to the record holding its
  site: `_compute_effective_coverage` never passed `record_starts` to
  `_mark_window`, `_union_coverage` does not either, and the rewrite preserves
  that deliberately so the change could be verified by exact agreement. On a
  multi-record reference both therefore let a window cross a contig boundary.
  That is Phase 6's subject and is open, not fixed here.
- The saturation branch, for a primer bound all of the time, cannot be reached
  by any valid reaction: it needs a melting temperature near 200 C and a 12-mer
  tops out around 70. It is handled because log(1 - 1) is negative infinity and
  the alternative failure is a silent NaN, and it is exercised through a
  conditions stub rather than by a real design.
- Making a modelled quantity faster does not make it more accurate. Coverage
  here is still site geometry weighted by a two-state occupancy approximation,
  not measured amplification. See
  [evidence_matrix_2026-09-15.md](evidence_matrix_2026-09-15.md).
