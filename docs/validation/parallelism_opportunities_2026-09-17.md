# Where parallel execution would buy wall clock, and where it would not

Reproduce the attribution measurements with
`scripts/benchmarking/measure_objective_attribution.py`.

Read of the source on 2026-09-17, branch `audit-phase-4-increment-3`, with short
supporting measurements on one machine (11 cores, Python 3.11.14, numpy 2.x,
macOS, `multiprocessing` start method `spawn`).

**The headline is that the most expensive thing this tool does is not waiting on
parallelism.** One objective evaluation costs 14.9 ms
(`docs/validation/objective_evaluation_cost_2026-09-17.md`) and the dominant term
is an occupancy-weighted coverage loop that touches the whole target genome once
per primer. Replacing that loop with an interval sweep is 19x to 21x faster on
one core, measured against the real implementation on the real index, and agrees
with it to 6e-9. Eleven cores cannot buy 20x on this loop; eight threads buy
3.7x. Every parallel option below competes against that number.

A second result: the attribution in the cost note is wrong about which function
dominates, and the reason is recoverable from the two benchmark scripts. See
**What contradicts the premises**.

## 1. What is already parallel

`neoswga/core/utility.py:37` (`create_pool`) and `:58`
(`create_pool_with_progress`) are the shared process-pool primitive. Every other
site builds its own executor.

| Site | Primitive | Unit of work | Reached on a default run |
|---|---|---|---|
| `core/kmer_counter.py:513` | `ThreadPoolExecutor` | one jellyfish subprocess per k | yes |
| `core/filter.py:682` | `ThreadPoolExecutor` | one count table per (prefix, k) | yes |
| `core/string_search.py:800` | `multiprocessing.Pool` | one (prefix, k) position scan | no, fallback only |
| `core/primer_attributes.py:332` | `multiprocessing.Pool` | one (prefix, k) Gini batch | no, see below |
| `core/dimer.py:307` | `ThreadPoolExecutor` | one primer pair | only above 200 primers, via `dimer_validator.py:88` |
| `core/thermodynamic_filter.py:498` | `ProcessPoolExecutor` | one primer pair, bounded window | no, see below |
| `core/rf_preprocessing.py:505`, `:555` | `create_pool` | one primer | only under `--amp-model` |
| `core/dimer_matrix.py:195` | numpy boolean matmul | whole pool at once | yes |

Four of those eight deserve a note.

**The Gini pool is dead on the `filter` path, and deliberately so.**
`primer_attributes.get_gini_from_txt` takes the multiprocessing branch only when
`position_cache is None`. `pipeline.step2` always supplies one:
`_scan_foreground_positions` returns `string_search.get_positions(...)`
unconditionally and the result is passed to `filter.get_gini` as
`position_cache=`. So on every run of `neoswga filter` the Gini stage is
single-threaded, and the pool exists for standalone callers only.

Re-checked 2026-09-17 against the source. Unlike the other two dead sites, this
one is a recorded choice rather than an oversight: the comment above the branch
says "no multiprocessing needed since the data is already loaded -- avoids
pickle overhead", and the data in question is the foreground position arrays. So
item 4 in the ranking below is not "reach the pool" but "measure whether the
serial in-cache path is actually the cheaper one at 491,836 candidates", which
is a different and unanswered question. Reaching the pool by withholding the
cache would reintroduce the pickle cost the comment is about.

**The per-k position pool is dead whenever pyahocorasick is installed.**
`string_search.get_positions` takes the single-pass Aho-Corasick branch
(`core/string_search.py:704`) and falls through to the pool at `:800` only when
the package is absent or `overwrite=True`. pyahocorasick is installed in this
environment. The Aho-Corasick branch iterates prefixes and k values serially.

**The heterodimer process pool has no production caller.**
`ThermodynamicFilter.filter_candidates` defaults to `check_heterodimers=True`
(`core/thermodynamic_filter.py:277`) but both callers pass False:
`core/hybrid_thermo_screen.py:107` and `core/multi_genome_pipeline.py:373`, each
with a comment saying why. The bounded sliding window at `:498` is careful work
that nothing reaches.

**The pairwise dimer thread pool buys close to nothing.** `_check_dimer_pair`
runs `dimer.is_dimer_fast`, a pure-Python substring test, so the threads at
`core/dimer.py:307` serialise on the GIL. It is also superseded:
`dimer_matrix.build` computes the same relation as one boolean matmul, measured
in that module's own docstring at 0.018 s for 449 primers.

One arithmetic defect while in the area. `pipeline.py:840` calls
`run_jellyfish(genome, prefix, min_k, max_k)` without `cpus`, so the default of 4
stands, and `kmer_counter.py:511` computes
`max_workers = min(num_k, cpu_count // cpus)`. On this 11-core machine that is 2
concurrent k values at 4 threads each, 8 threads of 11. On a 64-core machine it
is 7 k values at 4 threads, 28 of 64. The formula never sees `parameter.cpus`.

## 2. Classification

(a) embarrassingly parallel and worth it, (b) parallel but dominated by transfer
or setup, (c) inherently sequential, (d) already fast enough.

| Work | Class | Unit of work | Reason |
|---|---|---|---|
| Per-candidate objective scan in a greedy step | (d), after a fix | `dominating_set_optimizer.py:556` | No production caller passes `objective=`; and the per-call cost is about 20x reducible on one core |
| Objective evaluation itself | (a) as threads, but change the algorithm first | `base_optimizer.py:1419` | 3.7x measured with 8 threads; 19x to 21x measured by changing the algorithm |
| Occupancy-weighted swap and beam evaluation | (b) | `panel_beam.py:128`, `swap_refinement.py:97` | Independent per pair, but the budget check sits inside the loop, so parallelising changes the delivered panel |
| Background position scan (`filter`) | (a) | `string_search.py:236` | About 70 s of a 153 s filter at `all_qc`; chunking machinery already present |
| Foreground position scan (`filter`) | (b) | `string_search.py:236` | Same loop, but the foreground is 1.27 Mb and the pattern set the same, so the constant dominates |
| HDF5 index build (`PositionCache`) | (c) in practice | `position_cache.py:324` | h5py serialises on a global lock; processes would have to ship the arrays back |
| k-mer counting | (d) | `kmer_counter.py:513` | Already threaded over k, and jellyfish is a subprocess, so the GIL is irrelevant. Only the worker arithmetic is wrong |
| Filter per-k count-file scan | (b) | `filter.py:682` | Inner work is Python `line.split()`, GIL-bound; processes would work but the stage is not the cost |
| Gini computation | (a) | `primer_attributes.py:201` | Pure per-primer arithmetic over cached positions; the pool exists and is simply not reached |
| Tm and sequence-quality gate | (a), modest | `pipeline.py:1268` | 25.2 us per primer measured, so 12.4 s for 491,836 candidates on one core |
| Pairwise dimer screen | (d) | `dimer_matrix.py:195` | Already one matmul. At inventory scale it is a representation problem, not a parallelism one: 491,836 squared booleans is 242 TB |
| Ensemble over methods | (a) | `unified_optimizer.py:355` | Independent by construction, and `_reseed` was written to make them order-independent |
| `plan-pool` size sweep | (a) | `pool_planner.py:302` | Each size is an independent `optimizer.optimize` |
| `design_sweep` condition grid | (a), but unreachable | `pool_design_sweep.py:83` | Coarsest grain in the repo, and no command calls it |
| `set_size_optimizer` frontier | (a) | `set_size_optimizer.py:358` | At most `num_refine` (default 5) independent optimizations; phase 1 is closed form |

### The objective evaluation, measured

`compute_pool_metrics` (`core/pool_metrics.py:45`) computes five fields. Timing
its parts with a stub position cache, a 1.27 Mb target, a 144 Mb background,
20 foreground and 9 background sites per primer per strand, 24 primers, reach
3000:

| Term | Per evaluation |
|---|---|
| `_compute_effective_coverage` (`base_optimizer.py:1419`) | 12 to 15 ms |
| `_effective_site_load` -> `weighted_site_load` | 2.7 ms |
| `_union_coverage` for raw `fg_coverage` | 1.4 ms |
| foreground position set union | 0.08 ms |
| background position set union | 0.03 ms |

Reproduced end to end through the real objective, same cache and panel for both
rows:

| Configuration | Per call |
|---|---|
| `conditions=ReactionConditions(temp=30)`, `coverage_metric="effective"` | 16.83 ms |
| `conditions=None`, `coverage_metric="raw"` | 0.95 ms |

`_compute_effective_coverage` allocates nothing per call but does, per primer,
one `window[:] = False` over the target length and one
`not_covered[window] *= 1 - theta` over the same length. The masked multiply
alone measured 0.482 ms at 1.27 Mb, so 24 primers is 11.6 ms before any site is
marked. The cost is therefore linear in panel size and linear in target length,
and nearly independent of the background.

The same quantity computed as an interval sweep, sorting all window endpoints and
accumulating `log(1 - theta)` across segments, costs 0.277 ms against 13.12 ms
for the existing loop, a factor of 47, and agrees to 5.5e-9 absolute on the same
input. An intermediate form that keeps a full-length difference array and one
`cumsum` gives 6.02 ms, a factor of 2.3, and is the easier of the two to write.

`weighted_site_load`'s 2.7 ms is a sum of per-primer terms that do not depend on
the panel at all: `mismatch_class_counts(primer, prefixes, max_mismatches)` is a
pure function measured at 50 us per 12-mer, of which 25 us is
`_variants_at_distance`, and `calculate_effective_tm` is 6 us. Table size does
not matter: `weighted_site_load` over 24 primers measured 1.36, 1.33 and 1.38 ms
against count tables of 200,000, 2,000,000 and 8,000,000 entries. So this term is
removable by memoising per primer, not by parallelising.

### The position scans, measured

Aho-Corasick throughput over a synthetic 20 Mb sequence of uniform base
composition:

| Patterns | Build | Scan | Throughput |
|---|---|---|---|
| 20,000 | 0.01 s | 0.64 s | 31.4 Mb/s |
| 100,000 | 0.10 s | 1.99 s | 10.0 Mb/s |

Throughput falls as the automaton grows, which is the expected cache behaviour
and is the reason splitting matters. A `filter` run at
`candidate_retention="all_qc"` on the Wolbachia pair indexes 491,836 candidates
against a 144 Mb host, so the automaton carries about 983,672 patterns across
both strands. The measured `filter` runtimes of 153.4 s for `all_qc` against
86.4 s for `post_gini`
(`docs/validation/wolbachia_retention_benchmark_2026-09-16.md`) put roughly 70 s
of that difference on this one loop, and it runs on one core.

The loop at `string_search.py:236` already walks the sequence in overlapping
windows of `MAX_SCAN_CHUNK`, with each match attributed to the chunk its start
falls in. That is exactly the decomposition a process pool needs, so the
correctness argument for splitting is already written and already tested.

## 3. Which primitive fits

| Opportunity | Primitive | Inner work | Data a worker needs |
|---|---|---|---|
| Objective evaluation | vectorise first, then threads | numpy masked ops on float32 and bool arrays | nothing new; shares the cache in-process |
| Background position scan | processes | pyahocorasick `iter()` yielding into a Python loop with `bisect` | the genome string (144 Mb) and the pattern list, or a saved automaton |
| Gini computation | processes | Python arithmetic over `int64` lists | position arrays for one (prefix, k) slice |
| Tm and quality gate | processes | Python string work | the primer strings only, a few MB |
| Ensemble over methods | processes | mixed, mostly Python | the whole `PositionCache` |
| `plan-pool` size sweep | processes | mixed, mostly Python | the whole `PositionCache` |
| `set_size_optimizer` frontier | processes | mixed, mostly Python | the whole `PositionCache` |

Threads are viable for the objective and useless everywhere else in that table.
Measured on the occupancy coverage loop:

| Workers | Per evaluation | Speedup |
|---|---|---|
| 1 | 13.97 ms | 1.00x |
| 2 | 7.05 ms | 1.96x |
| 4 | 3.97 ms | 3.47x |
| 8 | 3.76 ms | 3.67x |

numpy releases the GIL for the large masked operations, so the scaling is real up
to about four workers and then flattens, which is consistent with the loop being
memory-bandwidth bound rather than compute bound. That is also why the interval
sweep wins: it stops touching the array.

**Picklability and size decide the process cases.** `PositionCache` holds a dict
of `np.int64` arrays keyed on `(prefix, primer, strand)`, and
`position_cache.POSITION_DTYPE` must stay `int64` (Known Issue 7). The measured
index sizes on the Wolbachia design are 443.1 MB for the `all_qc` background and
18.9 MB for `post_gini`, plus 381.3 MB foreground. The start method here is
`spawn`, so a process pool over sizes, methods or frontier points pickles the
whole cache once per worker. At `all_qc` that is roughly 800 MB per worker
against a `plan-pool` run that takes 20 s: the transfer would dominate. At
`post_gini` it is about 400 MB, still poor value. Either case argues for a shared
memory-mapped index or for `StreamingPositionCache`
(`position_cache.py:745`), not for pickling.

One caveat on those figures. At the time of writing, uncommitted work in this
tree adds `PositionCache.load` (`position_cache.py:359`) and
`PositionCache.release` (`:406`), so the index a worker would need becomes the
frontier rather than the whole inventory. That changes the arithmetic above in
the direction that favours a process pool, and the estimate should be retaken
against the incremental cache once it lands rather than read off the 443 MB
figure.

The scans and the per-primer gates do not have this problem. `filter_extra` takes
a string. `mismatch_class_counts` takes a string and a path. The Gini worker
takes one k's slice of positions, and `get_gini_from_txt_for_one_k` already
threads `min_sites` through as an argument for precisely this reason.

## 4. What parallelism would break

**The search budget is part of the answer, not an accuracy knob.**
`docs/validation/wolbachia_search_budget_2026-09-16.md` measures panels at
budgets 1,000 and 100,000 sharing as little as 0.38 of their primers. In
`refine_by_swaps` (`swap_refinement.py:104`) and `beam_search`
(`panel_beam.py:138`) the budget is checked inside the candidate loop and
`evaluations` increments in scan order. Any parallel scan changes which
candidates fit inside the budget, so it changes the delivered panel. A
result-preserving parallel version has to fix the work set per layer first, score
it all, and only then apply the budget, which is a change to the search
definition rather than to its implementation.

**Order, not seeding, is what makes a run reproducible.** Greedy selection
involves no RNG at all: `dominating_set_optimizer.py` returns `order` rather than
`list(selected)` because per-process string hashing once produced five different
answers at seed 42. `step3_ordering.order_step3_rows` establishes the pool order
the optimizer reads, and reversing it moved a delivered set to Jaccard 0.600. Any
parallel scan must therefore recombine results in the original scan order, not in
completion order, and must break ties the same way `_is_better` does.

**Threads would break the ensemble's reproducibility.** `_reseed`
(`unified_optimizer.py:194`) sets the global `random` and `numpy.random` state
before each method. Two methods in two threads share that state, so the
re-seeding that was written to make methods order-independent would instead make
them interfere. Processes give each method its own interpreter state and preserve
the property.

**The memory ceiling binds before the core count does.** A single `filter`
against hg38 peaks at about 8.5 GB, so parallel filters over genomes are not
available on a normal machine; the retention benchmark also reports 2650 to
3349 MB peak resident for the Wolbachia pair. Splitting one scan across
processes multiplies the resident genome string unless the workers map it from
disk. `mismatch_counts.load_kmer_counts` holds one dict per (prefix, k) and
measured 3.6 s to build for 8,000,000 entries; four workers each holding their
own copy of a host-sized table is the same trap in a new place.

**Spawn against fork.** `rf_preprocessing.py:233` says worker data is "shared via
fork", which is false here: the start method is `spawn` on macOS and the code
works only because the pool passes an explicit `initializer`.
`primer_attributes.py:133` records the same hazard correctly and is the pattern to
copy: resolve configuration in the parent and pass it as data. Any new worker
must be checked for module-global reads, because a spawned worker sees
`parameter`'s import-time defaults, and the Gini default of 3 differing from a
configured 2 is not an error, only a different answer.

**A wrong incremental union is invisible.**
`docs/validation/pool_plan_profile_2026-09-16.md` declined an incremental
`_union_coverage` on exactly this ground, and the caution applies to the interval
sweep recommended here: circular wrap and the `record_starts` confinement that
`_mark_window` implements are where an off-by-one would produce a coverage figure
that is wrong without looking wrong. The remedy is the one already in the repo,
`tests/test_pool_metrics_agree_with_the_full_evaluation.py`, extended to
randomised circular and multi-record cases.

## 5. Ranking

Saving is per affected run, not per pipeline. Measured means measured here or in
a cited note; estimated means arithmetic over a measured rate.

| # | Change | Saving | Risk | Basis |
|---|---|---|---|---|
| 1 | Interval sweep for `_compute_effective_coverage` | 19x to 21x on the dominant term, measured on the real pair | Medium: record-boundary confinement, which neither version does | Measured against the real implementation; wrapping verified on five fixtures |
| 2 | Memoise the per-primer terms of `_effective_site_load` | 2.7 ms of 14.9 ms | Low: pure function of (primer, prefixes, max_mismatches) | Measured, panel-independence measured across three table sizes |
| 3 | Processes over genome chunks in the background position scan | About 70 s of a 153 s `all_qc` filter, times core count less overhead | Medium: worker must receive the genome or map it | Measured throughput, estimated saving |
| 4 | Parallelise the cached Gini path | Unmeasured share of an 86 to 153 s filter | Low: the spawn-safe argument threading exists | Read of the source only; the serial branch is a recorded choice, see above |
| 5 | Processes over sizes in `plan_pool` and over frontier points | Up to the number of sizes, 9 on the published sweep | Medium: 400 to 800 MB per worker under spawn | 20 s run measured; saving estimated |
| 6 | Processes over ensemble methods | Up to 4x on an ensemble run | Low for correctness, high for memory | Estimated |
| 7 | Processes over `filter_extra` | 12.4 s of a 153 s filter | Low | Measured rate, estimated total |
| 8 | Threads inside one objective evaluation | 3.7x, and only if item 1 is not done | High: contends with item 1 and with the budget semantics | Measured |

**Do item 1 first, and it is not parallelism.** It is the largest measured factor
available, it applies to every path that evaluates a panel, and it removes the
reason to parallelise the same code. Item 2 is a natural companion and is close
to free. The first genuinely parallel move is item 3, the background position
scan, because it is the largest single-threaded block in the pipeline, the chunk
decomposition and its correctness argument already exist in
`string_search.get_all_positions_multi_k`, and it is the one opportunity where the
worker does not need the position index.

Before item 6, note that the ensemble's redundancy is a larger effect than its
parallelism: `hybrid` returns a set identical to `dominating-set` at Jaccard 1.000
(Known Issue 8) and `background-aware` wraps the same `HybridOptimizer`.
Not running two of the four beats running four in parallel.

## What contradicts the premises

**The cost note attributes the 14.9 ms to the wrong functions.**
`docs/validation/objective_evaluation_cost_2026-09-17.md` states that the cost
"is dominated by the background, specifically `_effective_site_load` and a
background set union". Measured at 12-mer site counts, those two are about 2.7 ms
and 0.03 ms. A random 12-mer has an expected 17.2 sites in a 144 Mb host and 0.15
in a 1.27 Mb target, so the background union is tens of microseconds, not
milliseconds. The dominant term is `_compute_effective_coverage`, and it scales
with the target length and the panel size.

**The order of magnitude between the two profiles is the occupancy weighting, not
the background.** `scripts/benchmarking/profile_pool_plan.py:81` passes
`conditions=None` and `:93` sets `coverage_metric="raw"`. With `conditions=None`,
`_effective_site_load` returns immediately and `_compute_effective_coverage` is
never called, so the 1.4 ms synthetic figure never entered the occupancy path at
all. `scripts/benchmarking/measure_objective_cost.py:71` passes
`ReactionConditions(temp=30.0)` and the default `PoolConstraints()`, whose
`coverage_metric` is `"effective"`. Holding the cache, panel and both genomes
fixed and changing only that, the same objective measured 16.83 ms against
0.95 ms. The two published figures are measurements of two different quantities.

This does not weaken the note's conclusion. A bounded per-step scan is still
necessary; the reason is that the objective costs milliseconds, and that stands
whichever function spends them. What changes is where to spend effort: on the
target-length loop, not on the background.

**The 48.9 h greedy scan has no production caller.**
`optimize_greedy(objective=...)` at `dominating_set_optimizer.py:860` is passed an
objective by nothing in `neoswga/`; the only reachable objective-scored searches
are `refine_by_swaps` (`pool_planner.py:126` and `swap_refinement.py:195`) and
`beam_search` (`pool_planner.py:179`). The arithmetic in the note is correct and
the architectural conclusion follows from it, but the scan it prices is not a
path a user can currently run.

**A naive reading of the swap loop and the 14.9 ms figure does not reconcile with
the published sweep.** `refine_by_swaps` increments `evaluations` before the dimer
guard (`swap_refinement.py:108`) and calls `_score`, and therefore
`objective.metrics`, on every pair that survives the guard. At 14.9 ms a budget of
100,000 would be about 25 minutes per size, while the search budget note measures
23.0 s for nine sizes. So either the dimer guard rejects almost every pair before
it is scored, or the loop reaches a local optimum far inside its budget. Both are
plausible and neither is measured. It matters for item 8: if the surviving pair
count is small, parallelising that loop has little to work with.

**`design_sweep` is still unreachable.**
`tests/test_no_capability_is_unreachable.py:59` lists it with the reason "no
production caller; Phase 4 wires it", and `CandidateProvider`, `load_design_grid`
and `ensure_positions` sit beside it. The condition grid is the coarsest parallel
grain in the repository and there is no command that runs it.

## The interval sweep, verified against the real implementation

**Added 2026-09-17, and it corrects the figure above.** The 47x in the ranking
was a synthetic comparison. Re-measured against the real
`_compute_effective_coverage` on the real Wolbachia index, the speedup is 19x to
21x, not 47x. Reproduce with `scripts/benchmarking/check_interval_sweep.py`.

| Geometry | Panel | Loop | Sweep | Speedup | Absolute difference |
|---|---|---|---|---|---|
| linear | 12 | 8.69 ms | 0.46 ms | 18.8x | 4.0e-9 |
| linear | 24 | 14.63 ms | 0.68 ms | 21.5x | 6.2e-9 |
| circular | 12 | 8.19 ms | 0.44 ms | 18.6x | 4.0e-9 |
| circular | 24 | 14.44 ms | 0.78 ms | 18.5x | 6.2e-9 |

The agreement holds, and it holds on the geometry that carries the risk. A first
attempt to check wrapping used sites near both ends of a circular target and
found circular and linear identical, which looked like the flag not reaching the
function. It was not: with a site within one reach of each end, each site's
window already covers the other end, so wrapping changes nothing and the test
proved nothing. On five fixtures where wrapping does change the union, circular
and linear differ and the sweep matches the loop in every case:

| Case | Circular | Linear | Agrees |
|---|---|---|---|
| one site at 10 | 0.0297410 | 0.0153662 | yes |
| one site at 19,990 | 0.0297410 | 0.0153662 | yes |
| site at 10 plus one interior | 0.0556687 | 0.0412939 | yes |
| site exactly at 0 | 0.0297410 | 0.0148705 | yes |
| two primers both near the origin | 0.0329636 | 0.0205044 | yes |

**The remaining risk is record confinement, not wrapping.**
`_compute_effective_coverage` calls `_mark_window` without `record_starts`, so
it does not confine a window to the record holding its site. The prototype
matches that, which is right for a comparison and wrong for a replacement: a
real one has to decide the question rather than inherit it. That is the same
concern Phase 6 of the plan raises for the coverage path generally.

So the recommendation stands with a smaller number. 19x to 21x on one core still
exceeds the 3.7x that eight threads buy on the same loop, and it still applies
to every path that evaluates a panel.

## The attribution finding, re-measured on the real pair

Everything above was measured on synthetic arrays. The one finding that
contradicts a published note was re-taken on the real Wolbachia index built for
the increment 3 acceptance check: wMel `NC_002978.6` with the Drosophila
`GCF_000001215.4` background, the 2,000-candidate shortlist, panel size 24,
reach 3000, one variable at a time.

| Configuration | Per call |
|---|---|
| Occupancy-weighted coverage, background present | 18.91 ms |
| Raw coverage, background present | 0.92 ms |
| Occupancy-weighted coverage, no background at all | 15.81 ms |

| Panel size | Per call |
|---|---|
| 6 | 6.38 ms |
| 12 | 10.96 ms |
| 24 | 18.72 ms |
| 48 | 35.94 ms |

Removing the 144 Mb host costs 16% of the time; removing the occupancy
weighting costs 95%. The panel scaling is close to linear, about 0.62 ms per
primer over a constant, which is what a per-primer pass over a 1.27 Mb target
looks like and is not what a background term looks like. The synthetic
breakdown and the real pair agree, so the attribution correction stands on real
data rather than on a model of it.
`docs/validation/objective_evaluation_cost_2026-09-17.md` has been corrected.

The interval sweep HAS now been re-measured against the real implementation on
the real index, and the synthetic 47x does not survive it: see the section above.
The real figure is 19x to 21x.

## What this analysis does not establish

- No measurement here ran against a real genome pair. Every number attributed to
  this note was taken on synthetic arrays, synthetic count tables or a synthetic
  sequence of uniform base composition, sized to match the Wolbachia design. Real
  site distributions are clustered and real primers are selected for abundance in
  the target, so the per-primer site counts assumed in the objective breakdown
  (20 foreground, 9 background per strand) are plausible rather than measured.
- Single runs on one machine, one core count, one start method. No spread. The
  thread-scaling table in particular was taken while nothing else ran; under load
  it will be worse.
- The interval sweep has since been verified against the real implementation on
  the real index and on five wrapping fixtures. It still does not handle
  `record_starts` confinement, and neither does the loop it reproduces, so that
  question is open rather than answered.
- Nothing here measures the end-to-end effect of any proposed change. The 19x to
  21x is on one function; what a `plan-pool` run would actually save depends on
  how much of it is objective evaluation, which was not profiled after the
  focused evaluator landed. The call counts in
  `what_actually_bounds_the_search_2026-09-17.md` suggest that share is small
  unless a constraint binds.
- The background scan saving is an estimate from a throughput curve on random
  sequence, extrapolated across an order of magnitude in pattern count. A real
  host genome has repeat structure that changes both the automaton's cache
  behaviour and the match count.
- No claim is made that any parallel change here would preserve a delivered
  panel. Where a change would alter one, this note says so; where it says nothing,
  it has not been checked.
- The working tree carried uncommitted changes to `candidate_provider.py`,
  `candidate_source.py`, `pool_planner.py` and `position_cache.py` while this was
  written, and `pool_planner.py` moved under the read. Line references are as of
  the tree state on 2026-09-17 and will drift when that work lands.
- Coverage and selectivity remain modelled site geometry, not measured
  amplification. Making them faster does not make them more accurate. See
  `docs/validation/evidence_matrix_2026-09-15.md`.
