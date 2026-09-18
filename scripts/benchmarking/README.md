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

### `frontier_refill_sweep.py` -- does searching more candidates change the panel

Phase 4 increment 5's acceptance measurement. Runs one size row with refilling
off and on and reports refills used, which exhaustion ended it, how much of the
universe was examined, and the delivered coverage and density.

```bash
python scripts/benchmarking/frontier_refill_sweep.py 100.0 <design_dir>
```

The floor matters: a row that qualifies never refills, so a floor the shortlist
already meets measures nothing. See
[frontier_refill_2026-09-17.md](../../docs/validation/frontier_refill_2026-09-17.md).

### `floor_boundary.py` -- what does retaining candidates buy

Phase 5's measurement. Raises a selectivity floor until each candidate universe
fails, and compares where the shortlist, the post-Gini inventory and all hard-QC
candidates stop. Budgets are held equal so only the universe varies.

```bash
python scripts/benchmarking/floor_boundary.py <design_dir> post_gini 40,60,80,100,140
python scripts/benchmarking/floor_boundary.py <all_qc_design_dir> all_qc 60,80
```

The answer on the Wolbachia pair at panel size 12 is that retention buys
nothing: all three universes stop between 60 and 80, and the two larger ones
deliver identical panels. See
[retention_changes_no_delivered_panel_2026-09-17.md](../../docs/validation/retention_changes_no_delivered_panel_2026-09-17.md).

### `achievable_density.py` -- a density-only ceiling, and a caution

Exact arithmetic for the highest selectivity density an N-primer panel from a
pool can reach, using the fact that the occupancy loads are additive over
primers. Useful for one narrow question: whether retaining more candidates
raises that ceiling at all. On the Wolbachia pool it does not, which confirms
Phase 5 by an exact method.

**Do not read its number as headroom.** It ignores the coverage target and the
dimer screen. On the measured pool the panel achieving 79.807 has coverage
0.4042 against a 0.5 target and 30 dimerising pairs out of 66, and chasing that
apparent gap cost three Stage 1 search rules that could not help. With every
constraint in force nothing beats the search's 60.112 at coverage 0.6535. See
[stage_one_constraint_awareness_2026-09-18.md](../../docs/validation/stage_one_constraint_awareness_2026-09-18.md).

`selectivity_budget.py` holds the accounting both this and that note use.

### `selectivity_compatibility_tension.py` -- is selectivity at odds with compatibility

Compares the most selective candidates against random controls on pairwise dimer
compatibility, largest mutually compatible subset, and GC content. Answers
whether the dimer screen caps achievable specificity.

```bash
python scripts/benchmarking/selectivity_compatibility_tension.py <design_dir> 70
```

On the Wolbachia pool: pairwise compatibility is 8 to 10 points lower among
selective candidates and they are GC-richer, but the largest compatible subset
is not smaller, so compatibility is not the barrier. See
[no_search_headroom_on_this_pool_2026-09-18.md](../../docs/validation/no_search_headroom_on_this_pool_2026-09-18.md),
which also records three claims of mine that controls and constraints refuted.

### `candidate_gate_audit.py`

What the candidate gate admits, and how informative it is. Needs no design
directory and no external tool: it samples random k-mers at a fixed GC fraction
and runs them through the real `ReactionConditions` and occupancy code.

```bash
python scripts/benchmarking/candidate_gate_audit.py
python scripts/benchmarking/candidate_gate_audit.py --gc 0.41 --n 40000 --k 12
```

Produces every figure in
[pool_selection_audit_2026-09-18.md](../../docs/validation/pool_selection_audit_2026-09-18.md):
occupancy spread within the admitted pool, how far each polymerase's Tm floor
sits below its reaction temperature, which bound of the window binds as a
function of primer length, what an additive gains and loses per GC class, and
the shipped Tm window against a band on occupancy.

Sampling is over sequence space, not over a genome's distinct k-mers. It
measures a property of the gate and is not a pool composition.

### `published_gap_thresholds.py`

Tests whether a hard coverage-hole constraint derived from the polymerase reach
separates the wet-lab winners among the 18 published sets already in
`tests/validation/data/`. Needs no design directory and no external tool.

```bash
python scripts/benchmarking/published_gap_thresholds.py
```

It does not. Every reach-derived threshold rejects all 18, winners included, and
the only separating threshold implies a 20 kb reach, outside the calibrated 3.0
to 6.2 kb band and fitted to these points. Read it as a description instead:
panels enriching 96 to 120 fold carry 31 to 33 kb holes. See
[getting_ahead_on_spacing_2026-09-18.md](../../docs/validation/getting_ahead_on_spacing_2026-09-18.md).

## Stale

`benchmark_suite.py`, `run_benchmarks.py` and `benchmark_improvements.py` predate
the optimizer consolidation. `benchmark_suite.py` imports
`neoswga.core.milp_optimizer` and `neoswga.core.moea_optimizer`, both deleted,
and its config axes (greedy / milp / moea) no longer match the shipped method
set. They are kept for reference; they do not run as written.
