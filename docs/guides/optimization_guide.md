# Optimization Method Selection Guide

This guide helps you choose the right optimization method for your SWGA primer design.

## Plan pool size against coverage and specificity

Configured chemistry must be valid before a new optimization starts. Invalid
reaction settings and failures in the configured network Tm calculation now
stop that design instead of substituting a different buffer or ignoring
additives. Condition-free library evaluators remain available explicitly.

Standard optimization now attaches a shared panel objective before generating
proposals and applies bounded swap refinement under that objective. With
reaction conditions, the objective uses effective coverage and the configured
panel limits. `--minimize-primers` preserves the same metric and limits when
removing oligos. This can change previously delivered panels and adds bounded
refinement work; it does not establish improved experimental recovery.

`plan-pool` carries supported search settings from params.json, including
`stage1_objective_width`, refill controls, and swap budgets. An explicit CLI
swap-evaluation budget overrides the parameter-file value. Stage-one scoring
remains geometric unless its objective width is configured; refinement and
final acceptance use the selected panel objective. Chemistry-independent
sequence dimer limits remain mandatory.

After `count-kmers`, `filter`, and `score`, use `plan-pool` to compare oligo
counts and find small panels meeting explicit targets:

```bash
neoswga plan-pool -j params.json --primer-length 12 --max-size 64 \
  --coverage-targets 0.8 0.9 0.95 --min-selectivity-density 10 \
  -o pool_plan
```

Candidate generation must include the requested oligo length. Each size is
optimized separately, then coverage, background binding and dimer compatibility
are recomputed for the delivered panel. The command recommends the smallest
qualifying panel found; this is not a proof of the global minimum. A partial
panel is assessed using its actual count. Larger requests need not improve
coverage because the search is heuristic and panels are not necessarily nested.

With background indexes, specify `--min-selectivity-density`,
`--max-background-sites`, or both. The density ratio uses the application's
occupancy-weighted binding loads per reference base; it is not predicted fold
enrichment. The exact background-site limit counts matches separately from that
weighted metric. A ratio of 10 above is an illustrative comparison threshold,
not a validated experimental cutoff. `--no-background` explicitly permits
planning without assessing specificity.

A panel that misses a constraint gets one bounded second attempt before the row
is reported as ineligible. First a swap refinement, which keeps the panel size
and accepts only swaps that reduce the number of violated constraints or, once
none are violated, raise coverage. If that reaches a local optimum still in
violation, a bounded beam rebuilds at the same size, keeping four partial panels
alive rather than one. A beam is worth the cost here because the selectivity
density floor is not monotonic under additions: a primer that binds the target
well and the background little raises the ratio, so the locally best first
choice can be the one that puts the qualifying panel out of reach.

Both attempts run inside the swap evaluation and time budgets, and the beam is
skipped when its bound does not fit what is left of them. The repaired panel is
re-evaluated through the same objective the row is judged on, so the reported
coverage describes the panel actually delivered. A dimer violation is not
repaired: the objective does not see dimers, so nothing steers the swap score
towards removing one, and a panel arriving with a dimerising pair means an
upstream relaxation fired. The report names the rows that needed a repair, in
the `repaired_by` column and beside the reason in the HTML table. Pass
`--no-repair` for the earlier behaviour, where the optimizer's first answer is
the one reported.

The default coverage metric is occupancy-weighted (`effective`); `--coverage-metric
raw` uses the union of windows around binding sites. `--coverage-reach` controls
their radius (3,000 bp by default for phi29). These estimates do not establish
sequencing recovery or coverage uniformity. Dimer limits remain strict; a target
that is not reached is reported as such.

Outputs include an HTML report, coverage-versus-count figure, all evaluated
panels in CSV/JSON, and FASTA files for qualifying recommendations. Use a new
output directory for each run. See the [Wolbachia wMel example](../../examples/wolbachia_pool_design/README.md)
for a comparison against the full Drosophila reference.

### Use reports in design scripts

The report is organism-independent and is generated automatically by every
`plan-pool` run. It includes the target and background filenames, design settings,
coverage and specificity curves, and oligo exports. Backgrounds can be host,
environmental, or other reference sequences supplied in the design parameters.

Regenerate a report from saved results without loading genomes or rerunning
optimization:

```bash
neoswga report-pool --input pool_plan/pool_plan.json \
  --title 'Target against sample background' -o design_report
```

`--input` also accepts the directory containing `pool_plan.json`. This command
preserves the saved coverage estimates, constraints, recommendations and
provenance; changing the title does not rerun or change the design. Use
`plan-pool` for a new design or different coverage/specificity requirements.
Both commands require a new or empty report directory, preventing stale FASTA
exports from an earlier design. Reports made from older saved results display
missing metadata as `not recorded`.

Python scripts can use the same renderer:

```python
import json
from pathlib import Path
from neoswga.core.pool_plan_report import write_pool_plan

plan = json.loads(Path("pool_plan/pool_plan.json").read_text())
report = write_pool_plan(plan, "design_report")
```

Keep the generated HTML, PNG, JSON, CSV and FASTA files together when sharing
the report. Report generation happens after the pool-size search completes;
it does not display live optimization progress.

## Quick Decision Tree

```
                    Start Here
                        |
                        v
              Is this for clinical/
              diagnostic use?
                   /        \
                 Yes         No
                  |           |
                  v           v
           background-    How many
              aware       candidates?
                             |
                       /          \
                    >500         <=500
                     |             |
                     v             v
              dominating-       hybrid
                  set          (default)
```

**Summary:**
1. **Clinical/diagnostic?** Use `background-aware`
2. **Large candidate pool (>500)?** Use `dominating-set`
3. **Otherwise:** Use `hybrid` (default)

## Pairwise dimer policy

Greedy selection in `dominating-set`, `hybrid`, `background-aware`, and
`network` now stops when no compatible candidate can be added under
`max_dimer_bp`. It may return fewer primers than requested. A stalled greedy
search does not establish that a larger compatible set is impossible.

To allow relaxation when the greedy search stalls, pass
`--allow-dimer-relaxation` or set `"allow_dimer_relaxation": true` in
params.json. Each admission made without the pairwise screen is logged. The
screen is restored for the next pick. `clique` always keeps its compatibility
constraint. This option does not change upstream thermodynamic filtering.

The validation report records the policy, requested and delivered counts, and
existing dimer checks on the delivered panel. A panel below the requested size
is reported as partial. Unsupported dimer thresholds raise an error instead of
disabling the screen.

Effective coverage is reported as `null` when reaction conditions are absent;
a computed zero remains zero in scoring. Coverage is calculated separately for
each target and weighted by sequence length. Gap statistics use a wrap-around
gap for circular targets and separate terminal gaps for linear targets.

## Bounded swap refinement

For `hybrid` and `background-aware`, opt into a bounded one-for-one search:

```bash
neoswga optimize -j params.json -m hybrid --refinement-method swap \
    --swap-max-evaluations 10000 --swap-max-seconds 10
```

The equivalent params.json keys are `refinement_method`,
`swap_max_evaluations`, and `swap_max_seconds`. Network refinement remains the
default. Swap mode starts with the greedy panel at the requested size and
considers the full thermodynamically filtered candidate pool. It preserves
fixed primers, does not increase panel size, and screens incoming primers
against the retained panel. It does not repair all conflicts in a baseline
created with explicit dimer relaxation, or fill a smaller baseline to budget.

The objective first increases covered bases on the optimizer's bins. In
`background-aware` mode, equal-coverage swaps prefer fewer background binding
sites. This mode replaces background pruning and network removal; it does not
optimize network connectivity or occupancy-weighted coverage. Final network
metrics and any requested simulation validation are still computed.

The search stops at a local optimum, the evaluation limit, or the cooperative
time limit. Time is checked between candidate swaps. Pool preprocessing and
final reporting are outside that limit, so it is not a whole-run timeout.
The log records accepted swaps, evaluations, and the stopping reason. A zero
budget retains the initial panel.

A [36-run comparison on saved Prevotella pools](../validation/refinement_real_pools_2026-09.md)
found that the default 10,000 evaluations did not complete a full pass and made
no swaps. An explicit `--swap-max-evaluations 100000` improved coverage in five
of six pool/size combinations, but increased exact background sites in two.
Network remains the default; inspect background binding as well as coverage
when using the larger swap budget.

## Constrained coverage benchmarks

`DominatingSetOptimizer.optimize_ilp` and `coverage_upper_bound` now enforce
`max_dimer_bp` by default. `fixed_primers` are mandatory and count toward the
**total** `max_primers` budget. Pass `enforce_dimers=False` explicitly to obtain
a coverage-only comparison. The objective remains base-weighted binned
coverage, not predicted experimental amplification.

An exact-solver result distinguishes `coverage` (the incumbent),
`coverage_upper_bound`, `proven_optimal`, `feasible`, and `mip_gap`. An
infeasible model has no coverage measurement. The LP helper returns the
solver's upper bound, not a possibly suboptimal incumbent. Solver time limits
exclude model construction.

`scripts/benchmarking/optimality_gap.py` now compares strict greedy selection,
swaps, and constrained ILP/LP bounds. Its random panels are explicitly labelled
unconstrained and are not compatible-panel baselines. A percentage optimality
gap is reported only when the integer optimum is proven.

## Method Comparison

| Method | Speed | Coverage | Specificity | Best For |
|--------|-------|----------|-------------|----------|
| `hybrid` | Slow, and superlinear in set size | Excellent | Good | General use (default) |
| `dominating-set` | Fast, and flat in set size | Excellent | Fair | Large pools, quick results |
| `background-aware` | Slowest | Good | Excellent | Clinical, low background |
| `network` | Slow | Fair | Good | Tm-balanced sets |
| `clique` | Medium | Fair | Best measured | Dimer-free primer sets |
| `ensemble` | Slowest | Best of the above | Varies | When unsure which to use |

> **Read this table with the measurements below.** On the one target measured
> so far, `hybrid`, `dominating-set` and `background-aware` returned the
> *identical* primer set at every size tested, so the speed column is the only
> column separating them. `genetic`, `moea`, `milp` and `greedy` were removed
> from this guide on 2026-08-31: those optimizers no longer exist and the CLI
> rejects the names.

## Detailed Method Descriptions

### hybrid (Default)

**Two-stage approach combining coverage and connectivity.**

```bash
neoswga optimize -j params.json --optimization-method=hybrid
```

**How it works:**
1. Stage 1: Uses dominating-set to select primers with good coverage
2. Stage 2: Refines selection using network connectivity analysis

**Strengths:**
- Good balance of speed and quality
- Both coverage and amplification optimized
- Suitable for most applications

**Weaknesses:**
- Not the fastest option for very large pools
- Background minimization not explicit

**When to use:**
- General primer design
- When you want reliable results without tuning


### dominating-set

**Graph-based coverage optimization (8x faster).**

```bash
neoswga optimize -j params.json --optimization-method=dominating-set
```

**How it works:**
- Models primer selection as a set cover problem
- Greedy selection of primers covering the most uncovered regions
- Provable ln(n) approximation to optimal

**Strengths:**
- Very fast (8x faster than hybrid)
- Guaranteed coverage bounds
- Works well with large candidate pools

**Weaknesses:**
- Ignores amplification network structure
- Does not consider background binding

**When to use:**
- Large candidate pools (>500 primers)
- Quick screening runs
- When coverage is primary concern


### background-aware

**Three-stage optimizer with explicit background minimization.**

```bash
neoswga optimize -j params.json --optimization-method=background-aware
```

**How it works:**
1. Stage 1: Filter by background binding frequency
2. Stage 2: Score by selectivity ratio (fg/bg)
3. Stage 3: Optimize coverage subject to background constraint

**Strengths:**
- Lower background binding than `hybrid`: host sites in the delivered panel fall 7-35% against `hybrid` at n=24 and n=36, measured against hg38 on the three GC-tier designs, for 0.1-3.1 points of coverage. The figure of
  10-20x that this guide used to quote was never reproduced
- Designed for clinical samples with high background
- Explicit selectivity optimization

**Weaknesses:**
- Slower than other methods
- May sacrifice some coverage for specificity

**When to use:**
- Clinical/diagnostic applications
- Samples with high host DNA
- When background contamination is a concern


### network

**Tm-weighted amplification network optimization.**

```bash
neoswga optimize -j params.json --optimization-method=network
```

**How it works:**
- Builds amplification network from primer binding sites
- Weights edges by Tm compatibility
- Penalizes primer pairs with dimer potential
- Selects primers maximizing network connectivity

**Strengths:**
- Considers thermodynamic compatibility
- Dimer-aware selection
- Good for uniform amplification

**Weaknesses:**
- Slower for large networks
- May not maximize coverage

**When to use:**
- When Tm uniformity is important
- Primer sets with dimer concerns
- Uniform amplification desired


### clique

**Dimer-free primer sets via clique finding.**

```bash
neoswga optimize -j params.json --optimization-method=clique
```

**How it works:**
- Builds a compatibility graph where edges connect dimer-free primer pairs
- Finds maximum clique (largest set with no primer dimers)

**Strengths:**
- Guarantees dimer-free primer sets
- Suitable when dimer avoidance is critical

**Weaknesses:**
- May sacrifice coverage for dimer-freedom
- Computationally intensive for large candidate pools

**When to use:**
- When primer-dimer formation is a primary concern
- Multiplex applications requiring dimer-free sets


### ensemble

**Runs several methods on one shared position cache and keeps the best.**

```bash
neoswga optimize -j params.json --optimization-method=ensemble
neoswga optimize -j params.json --optimization-method=ensemble --ensemble-combine=union
```

With a shared panel evaluator, selection prioritizes constraints, then objective
coverage and background load. The reported `normalized_score` is a [0,1] descriptive value across optimizers,
weighted by `--application`; the raw `score` is not comparable between methods.
The runner-up table is written to `step4_improved_df_summary.json` as
`ensemble_comparison`. `--ensemble-combine union` additionally re-optimizes over
the pooled primers from all methods, guarded so it never returns a worse set.

Note that on a target where the methods converge — as they did on the one
measured here — the comparison table will show several indistinguishable
runners-up. That is information about the target, not a failure of the ensemble.

### Host-Free Optimization

When no background genome is available, use the `--no-background` flag:

```bash
neoswga optimize -j params.json --no-background
```

This skips all background-related scoring and focuses on target genome coverage and primer compatibility.


## Command Line Examples

### Basic Usage

```bash
# Default (hybrid)
neoswga optimize -j params.json

# Explicit method selection
neoswga optimize -j params.json --optimization-method=dominating-set

# View detailed method comparison
neoswga optimize --method-guide
```

### Clinical Workflow

```bash
# For clinical samples with high host DNA
neoswga optimize -j params.json \
  --optimization-method=background-aware \
  --num-primers 10
```

### Fast Screening

```bash
# Quick screening with large candidate pool
neoswga optimize -j params.json \
  --optimization-method=dominating-set \
  --num-primers 15
```

### High-Quality Design

```bash
# Thorough optimization for important designs
neoswga optimize -j params.json --optimization-method=hybrid
```

`iterations` is a params.json key, not a flag. It bounds the search for
ALTERNATIVE sets only, so raising it offers more alternatives; it does not make
the primary selection more thorough.

## Iterative Design Workflow

For iterative wet-lab optimization:

```bash
# 1. Design initial set
neoswga optimize -j params.json --optimization-method=hybrid

# 2. Predict efficiency before synthesis
neoswga predict-efficiency -j params.json \
  --primers SEQ1 SEQ2 SEQ3 SEQ4 SEQ5 SEQ6

# 3. After testing, expand with additional primers
neoswga expand-primers -j params.json \
  --fixed-primers SEQ1 SEQ2 SEQ3 \
  --failed-primers SEQ4 \
  --num-new 4 \
  --output expanded/
```

## Performance Benchmarks

Measured 2026-08-31 on a 3.17 Mb bacterial target against human chr21, 2000
candidates, seed 42, phi29 at 30 C, one run each (repeat noise about 9%). Driven
through the CLI, so the figures include dispatch, position-cache build and
metric computation. Reproduce with `scripts/benchmarking/sweep_optimize.py`.

| Method | 12 primers | 32 primers | 64 primers | 128 primers | fg_coverage at 32 |
|--------|---:|---:|---:|---:|---:|
| dominating-set | 8.3 s | 8.2 s | 9.0 s | 8.6 s | 0.5376 |
| hybrid (default) | 13.7 s | 63.3 s | 347 s | 2239 s | 0.5376 |
| background-aware | 20.1 s | 133 s | — | — | 0.5376 |
| clique | 13.1 s | 34.7 s | — | — | 0.4499 |
| network | 64.9 s | 377 s | — | — | 0.4283 |

Three things this table is for. **`dominating-set` is flat in set size** while
coverage climbs; every other method's cost curve is superlinear. **The first
three rows returned the identical primer set** at every size tested (Jaccard
1.000), so the extra runtime bought nothing on this target. **`network` is worse
on coverage and costs 46x more** at 32 primers.

`clique` is the one method that is not dominated: it gives up coverage but was
the only one to return `dimer_risk` 0.0000, and it had the best selectivity
density of any method here (42.4 against 31.7).

Caveats that belong with these numbers: one target, one k, one run each. Enough
to show that a method returning an identical set at 38x the cost has a defect;
not enough to change a shipped default on. See
[`AUDIT_2026-08_alternatives_and_scaling.md`](../reports/AUDIT_2026-08_alternatives_and_scaling.md)
(F5, F5b) for the full measurement and its limits.

## Troubleshooting

### Poor coverage

Try:
1. Use `dominating-set` for coverage-focused optimization
2. Increase `--num-primers`
3. Check if candidates have sufficient binding sites

### High background

Try:
1. Use `background-aware` method
2. Pre-filter with stricter `max_bg_freq`
3. Check background genome is correctly specified

### Slow optimization

Try:
1. Use `dominating-set` for speed
2. Reduce candidate pool size
3. Decrease `iterations` in params.json

### Poor Tm uniformity

Try:
1. Use `network` method
2. Adjust `min_tm` and `max_tm` parameters
3. Consider polymerase-specific presets

### Shared panel stages

`OptimizationRequest` and `run_panel_search` in `core/optimization_service.py`
provide the panel-stage entry point used by standard optimization, pool planning
(including condition-grid rows), contraction, and expansion with a resolved design context.
Reference loading and candidate-frontier orchestration remain with the front ends.
The optimizer carries the resolved chemistry, indexes and evaluator; a request
adds the candidate frontier, fixed/excluded oligos, size bounds and stage choices.

The sequence is proposal selection, bounded repair of unmet limits or an explicit
coverage target, shared-objective refinement, and optional pool reduction.
After a reduction, repair/refinement can run again before further deletion.
Deletion preserves the requested coverage, panel limits and fixed oligos.
These are bounded searches for smaller qualifying pools, not proofs of minimum
pool size. Swap and deletion passes each use the configured evaluation/time
budget; multiple stages may consume more total computation than one pass.

Ensemble proposals use constraint-first objective ranking when a common evaluator
is available. Result-only library optimizers retain their normalized-score
comparison. The normalized score remains a reported metric, but does not override
an explicit panel constraint. Sequence-based dimer checks and the optional
stability floor are checked on changed panels at the configured temperature;
this does not establish a full additive-dependent dimer equilibrium model.

Results include stage history with pool sizes, objective coverage, failed limits,
time, and available search counts/stopping reasons. Pool reports show this history
under each row. Expansion labels deficit recovery separately from genome coverage
and retains the panel constraints when using a deficit objective. It refuses a
configured expansion that cannot produce a qualifying pool.

Pool planning refills a candidate frontier when constraints or the requested
coverage remain unmet, within `max_frontier_refills`. It retains the best feasible
panel found and records all frontier attempts. Standard optimization and expansion
still need a common frontier/refill policy; this stage runner does not by itself
make every retained candidate reachable from those front ends.

`contract-set` uses the same bounded deletion stage and effective-coverage metric.
It scans the supplied primers against the reference FASTAs without changing
pipeline indexes, keeps target/background circularity separate, and reports raw
coverage alongside effective coverage. Its existing `--min-coverage 0` option is
supported. When no qualifying reduction is found, the original panel is retained
and `meets_target` states whether it actually satisfies the request. Removal
traces retain per-primer coverage; stage history also includes timing and budgets.
