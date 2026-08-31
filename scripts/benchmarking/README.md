# Benchmarking scripts

## Current (added by the 2026-08 audit)

These drive the **CLI**, because that is what users run: a library-level harness
misses the dispatch, the position-cache build and the metric computation that
dominate a real invocation. See
[docs/AUDIT_2026-08_alternatives_and_scaling.md](../../docs/AUDIT_2026-08_alternatives_and_scaling.md)
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
The ILP that ships in `dominating_set_optimizer.optimize_ilp` answers a
different question (fewest primers covering every reachable bin), so it is
infeasible below the minimum cover and cannot bound the greedy result. This
module builds the max-coverage formulation over the same `BipartiteGraph` bins
and the same `extension_reach`, and reports both the integer optimum and the
LP relaxation.

Needs a solver: `pip install mip`.

### `optimality_gap.py` — greedy vs exact vs random

```bash
cd <workdir> && python ../../scripts/benchmarking/optimality_gap.py . 6,12,32,64 100
```

Writes `gap_results.json`: greedy coverage, ILP optimum, LP bound, the gap, and
the distribution of N random sets of the same size as a floor. Run it from the
working directory, since `PositionCache` resolves HDF5 prefixes relative to cwd.

## Stale

`benchmark_suite.py`, `run_benchmarks.py` and `benchmark_improvements.py` predate
the optimizer consolidation. `benchmark_suite.py` imports
`neoswga.core.milp_optimizer` and `neoswga.core.moea_optimizer`, both deleted,
and its config axes (greedy / milp / moea) no longer match the shipped method
set. They are kept for reference; they do not run as written.
