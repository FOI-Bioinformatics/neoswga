# Benchmarking scripts

## Current (added by the 2026-08 audit)

These drive the **CLI**, because that is what users run: a library-level harness
misses the dispatch, the position-cache build and the metric computation that
dominate a real invocation. See
[docs/reports/AUDIT_2026-08_alternatives_and_scaling.md](../../docs/reports/AUDIT_2026-08_alternatives_and_scaling.md)
for the results they produced.

### `sweep_optimize.py` — how optimize scales in set size

```bash
python sweep_optimize.py <workdir> <sizes> <methods>
python sweep_optimize.py ./tierM 6,12,32,64,128 dominating-set,hybrid
```

`<workdir>` needs a `params.json`, the genomes it names, and a `step3_df.csv`.
Records wall time, peak child RSS and the reported metrics to
`<workdir>/sweep_results.jsonl`, one line per run, and keeps each run's summary
under `<workdir>/run_<method>_S<size>/`.

The set size goes on the CLI rather than into `params.json` deliberately:
`param_validator.py` caps `num_primers` at 50 while `params.schema.json` allows
200, so a params.json asking for 64 is rejected while `--num-primers 64` runs.

### `max_coverage_bound.py` — exact and LP bounds on coverage

The optimizers answer "given a budget of S primers, which S maximise coverage?".
This module builds that formulation over the same `BipartiteGraph` bins and the
same `extension_reach`, and reports both the integer optimum and the LP
relaxation.

It was written because the ILP that then shipped in
`dominating_set_optimizer.optimize_ilp` answered a different question (fewest
primers covering every reachable bin), so it was infeasible below the minimum
cover and could not bound the greedy result. That has since been replaced with
this formulation, so `optimize_ilp` and `coverage_upper_bound` now give the same
answers as this harness. The harness remains the sweep driver, and is a useful
independent check on the library implementation.

Needs a solver: `pip install mip`.

### `optimality_gap.py` — greedy vs exact vs random

```bash
cd <workdir> && python ../../scripts/benchmarking/optimality_gap.py . 6,12,32,64 100
```

Writes `gap_results.json`: greedy coverage, ILP optimum, LP bound, the gap, and
the distribution of N random sets of the same size as a floor. Run it from the
working directory, since `PositionCache` resolves HDF5 prefixes relative to cwd.

### `compare_refinement.py` -- saved-pool refinement comparison

Runs each network/swap comparison in a fresh sequential process, without
modifying input pools or pipeline outputs. Records complete primer panels,
independent coverage and background metrics, optimizer and process timing,
worker peak RSS, source/input hashes and per-run logs. This benchmarks the
library optimizer, not the complete CLI pipeline.

```bash
python scripts/benchmarking/compare_refinement.py \
  --params tests/validation/genomes/params.json \
  --pools tests/validation/genomes/out10/step3_df.csv \
  --sizes 6 12 --repeats 2 --background-aware --timeout 60 \
  --output /tmp/refinement-comparison-new
```

Use `--modes swap --evaluations 100000` for a larger-budget check. Output
directories must be new. `--timeout` covers a complete worker;
`--search-seconds` covers only the cooperative swap search.

See [the measured comparison](../../docs/validation/refinement_real_pools_2026-09.md)
for the available real-genome pools and the limits of those measurements.

## What one panel evaluation costs (added 2026-09-17)

Four scripts, all reading a design directory that has already been through
`count-kmers`, `filter` and `score`. They read its position indexes and
candidate inventory, which are hundreds of megabytes and are not committed, so
the directory has to be built before any of them will run.

### `measure_objective_cost.py` -- the price of one evaluation

Times `PoolObjective.coverage(panel + [candidate])` at panel size 24 and
projects a full-universe greedy scan from it. Produced the 14.9 ms figure in
[objective_evaluation_cost_2026-09-17.md](../../docs/validation/objective_evaluation_cost_2026-09-17.md).

### `measure_objective_attribution.py` -- where that price goes

Same evaluation, one variable at a time: occupancy weighting on and off,
background present and absent, and four panel sizes. Written because the cost
note blamed the 144 Mb background; the background is 16% of it and the
occupancy-weighted coverage term is 95%.

### `count_objective_calls.py` -- how often it is actually called

Counts objective cache misses during a real `plan_pool` size row. One when no
constraint binds, about 3,200 when the repair fires. Also shows that
`swap_max_seconds` rather than `--swap-max-evaluations` is what stops the loop.
See [what_actually_bounds_the_search_2026-09-17.md](../../docs/validation/what_actually_bounds_the_search_2026-09-17.md).

### `check_interval_sweep.py` -- a cheaper way to compute it

Compares `_compute_effective_coverage` against the same quantity accumulated
over window endpoints rather than over bases: 19x to 21x faster, agreeing to
6e-9, on both linear and circular geometry. A prototype for measurement, not a
replacement; it does not confine windows to records, and neither does the loop
it reproduces. See
[parallelism_opportunities_2026-09-17.md](../../docs/validation/parallelism_opportunities_2026-09-17.md).

### `scan_width_sweep.py` -- how narrow the objective-scored scan can be

Phase 4 increment 4's acceptance measurement. Runs `plan_pool` at widths 16, 64,
256 and unbounded over the real pool and reports Jaccard, coverage, density,
objective evaluations and pairs seen per size.

```bash
python scripts/benchmarking/scan_width_sweep.py 40.0 <design_dir>
```

The first argument is the selectivity density floor. A floor nothing violates
leaves the repair unattempted and every width identical, which measures nothing,
so it defaults to a value that binds on the bundled Wolbachia design. See
[scan_width_2026-09-17.md](../../docs/validation/scan_width_2026-09-17.md).

## Stale

`benchmark_suite.py`, `run_benchmarks.py` and `benchmark_improvements.py` predate
the optimizer consolidation. `benchmark_suite.py` imports
`neoswga.core.milp_optimizer` and `neoswga.core.moea_optimizer`, both deleted,
and its config axes (greedy / milp / moea) no longer match the shipped method
set. They are kept for reference; they do not run as written.
