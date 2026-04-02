# NeoSWGA Optimizer Selection & Multi-Optimizer Workflow Audit

**Date:** 2026-03-31
**Team:** 6 agents (selection guidance, serial workflows, parallel workflows, architecture, workflow builder, gap analysis)

---

## Executive Summary

NeoSWGA has well-engineered multi-optimizer infrastructure (factory pattern,
cascade chains, parallel orchestrator, normalized scoring) that is largely
**hidden from CLI users**. The tool has 18 registered optimizers but the
wizard only exposes 3. Multi-agent parallel execution and 4 aggregation
strategies exist but have no CLI surface. The architecture scores 7.5/10
but the user experience for optimizer selection and comparison scores 5/10.

**Core finding: Ferrari engine behind a bicycle interface.**

---

## 1. Optimizer Selection: Is It Easy?

### Discovery Paths (rated 1-10)

| Path | Rating | Notes |
|------|:------:|-------|
| `--method-guide` flag | 9 | Excellent decision tree, comparison tables. Poorly discoverable. |
| `docs/optimization_guide.md` | 9 | Comprehensive per-method analysis. Decision tree only covers 3/18. |
| `--help` text | 8 | All 17 CLI methods listed with embedded decision tree. |
| Strategy presets (--strategy) | 7 | clinical/discovery/fast/balanced/enrichment. Only for normalized method. |
| Workflow selector | 5 | Only shows 3 methods, no descriptions. |
| Wizard | 5 | Only offers hybrid/dominating-set/background-aware. |
| Condition suggester | 4 | Barely mentions optimizer selection. |

**Average: 6.7/10**

### Key Gaps in Selection
1. Wizard shows 3 of 18 methods — users miss 15 options
2. No genome-aware auto-recommendation (GC-rich -> equiphi29, clinical -> background-aware)
3. `multi-agent` registered but missing from CLI choices
4. `weighted-set-cover` in CLI but absent from docs
5. Strategy presets only work with `normalized` — not cross-cutting

---

## 2. Serial (Chained) Optimizer Workflows

### What Exists
- **SerialCascadeOptimizer**: Clean architecture with 3 pre-built chains:
  - `coverage-then-dimerfree` (DS -> Clique)
  - `dimerfree-scored` (Clique -> Network)
  - `bg-prefilter-hybrid` (BG prune -> Hybrid)
- Supports arbitrary chains programmatically
- Per-stage candidate limits protect against expensive algorithms

### What's Missing
- **No --candidates-file on optimize**: Users can't feed one optimizer's output to another
- **No --cascade CLI flag**: Can't define custom chains from command line
- **No --fixed-primers on optimize**: Only `expand-primers` supports this
- **Scores not normalized across methods**: Manual comparison is unreliable
- **step4_improved_df.csv overwritten each run**: No --output-prefix

### Serial Workflow Grade: Architecture B+, User Experience C+

---

## 3. Parallel (Simultaneous) Optimizer Workflows

### What Exists
- **MultiAgentOrchestrator**: Well-engineered parallel runner
  - ThreadPoolExecutor with up to 8 workers
  - 4 aggregation strategies: BEST_SCORE, CONSENSUS, WEIGHTED_VOTE, ENSEMBLE
  - Per-agent timeouts (300s), fault tolerance, performance tracking
  - Registered as `multi-agent` (aliases: `ensemble`, `parallel`)
  - All 4 strategies tested and functional

- **CompositeOptimizer**: Exists in base_optimizer.py but runs sequentially. Superseded by MultiAgent.

### What's Missing
- **No CLI surface**: `multi-agent` not in --optimization-method choices
- **No --compare flag**: Can't compare methods from CLI
- **No --aggregation flag**: Can't choose CONSENSUS vs BEST_SCORE from CLI
- **File locking**: Two concurrent CLI runs overwrite same output file
- **ENSEMBLE strategy is a placeholder**: Just returns BEST_SCORE

### Parallel Workflow Grade: Architecture A-, User Experience D

---

## 4. Architecture Quality

| Aspect | Rating | Notes |
|--------|:------:|-------|
| Extensibility | 8/10 | 2 files, ~30 lines to add new optimizer |
| Clarity | 8/10 | Clean BaseOptimizer contract, decorator registration |
| Result Comparability | 7/10 | normalized_score exists but not used everywhere |
| Composability | 7/10 | Cascade + MultiAgent work but CompositeOptimizer uses raw score |

### BUG: Raw score vs normalized_score
CompositeOptimizer and WEIGHTED_VOTE aggregation compare raw `score` values
across optimizers. These are NOT comparable (enrichment-based 5.0 vs
coverage-based 0.999). Should use `normalized_score` instead.

### Architecture Strengths
- Thread-safe factory with lazy registration
- Frozen OptimizationResult dataclass prevents accidental modification
- normalized_score: weighted [0,1] composite of 5 clamped metrics
- fixed_primers parameter in BaseOptimizer for warm-starting

---

## 5. Gap Analysis: What's Missing

### 5 Recommended CLI Features

1. **`--compare hybrid,network,dominating-set`** (HIGH priority)
   - Runs specified methods, prints comparison table, saves all results
   - Uses MultiAgentOrchestrator internally
   - Example: `neoswga optimize -j params.json --compare hybrid,network,clique`

2. **`--output-prefix`** (HIGH priority)
   - Prevents step4_improved_df.csv overwrite
   - Example: `neoswga optimize -j params.json -m network --output-prefix step4_network`

3. **`--fixed-primers` / `--from-step4`** on optimize command (MEDIUM)
   - Enables iterative refinement: take hybrid result, refine with network
   - Example: `neoswga optimize -j params.json -m network --from-step4 step4_hybrid.csv`

4. **`--optimization-method=auto`** (MEDIUM)
   - Auto-selects based on genome size, GC, candidate count, background
   - Decision tree: clinical->background-aware, GC-rich->equiphi29,
     >1000 candidates->dominating-set, <200->coverage-then-dimerfree, else->hybrid

5. **`--optimizer-args key=value`** pass-through (LOW)
   - Exposes method-specific tuning (GA population, MILP timeout, DS bin_size)
   - Example: `neoswga optimize -j params.json -m genetic --optimizer-args population=200,seed=42`

### Auto-Select Decision Tree

```
Is this a clinical/diagnostic application?
  YES -> background-aware
  NO -> continue

Is genome GC > 65% and using equiphi29?
  YES -> equiphi29
  NO -> continue

Are there > 1000 candidates?
  YES -> dominating-set (fast) or hybrid (balanced)
  NO -> continue

Are there < 200 candidates?
  YES -> coverage-then-dimerfree or clique
  NO -> hybrid (default)

Want Pareto front / multi-objective?
  YES -> multi-agent with CONSENSUS aggregation
```

---

## 6. Verdict

### Can users select the correct optimizer?
**Partially.** The `--method-guide` and docs are excellent but poorly
discoverable. The wizard/workflow selector only show 3 options. A new user
will use the default (hybrid) and never discover that network, normalized,
or multi-agent might be better for their case.

### Can users run multiple optimizers?
**Programmatically yes, via CLI no.** The MultiAgentOrchestrator is
production-quality but has zero CLI exposure. Manual comparison requires
copying files between runs. There's no comparison command.

### What would make this excellent?
1. Add `--compare` flag (biggest impact, exposes existing infrastructure)
2. Add `multi-agent` to CLI method choices (one-line fix)
3. Add `--output-prefix` (prevents file overwrite)
4. Expand wizard to show all methods with guidance
5. Fix normalized_score usage in cross-optimizer comparison paths
