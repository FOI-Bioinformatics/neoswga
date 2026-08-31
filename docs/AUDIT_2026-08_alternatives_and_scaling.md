# Audit: step alternatives, scaling, and oligo-set quality

**Date:** 2026-08-31
**Baseline:** commit `9e7026e`, Python 3.11.14, macOS 26.2 arm64, 11 cores,
19.3 GB RAM, jellyfish 2.3.1, `spawn` start method.
**Status:** interim. Measured: the optimize step (all five methods), the option
inventory, the oligo-set quality question, the effect of the score step, and a
regression pass over the March 2026 audits. Not yet done: the count-kmers and
filter timing sweeps and the 3.1 Gb tier — listed under
[Not yet measured](#not-yet-measured) so the scope of these conclusions is not
overstated.

**Scope of the evidence.** Every optimizer conclusion rests on one target
(Prevotella vs human chr21) at one k. That is enough to show a defect exists —
a method that returns the identical set at 38x the cost does so regardless of
dataset — but not enough to fix a default on. Confirm on a second target before
changing shipped behaviour.

The harnesses are in
[`scripts/benchmarking/`](../scripts/benchmarking/README.md) —
`sweep_optimize.py` (set-size scaling), `max_coverage_bound.py` (exact and LP
bounds) and `optimality_gap.py` (greedy vs exact vs random). Every figure below
is reproducible from them; see [Method](#method) for the datasets and caveats.

---

## Summary

The optimizer is sound and its output is good: on a real bacterial target the
greedy set cover lands within 3.4-7.2% of a proven exact optimum, and **no
random set out of 100 beat it at any budget**. The March 2026 STEP4 finding
that the optimizer performed below random does not reproduce on real genomes.

The problems are elsewhere, and they are mostly about work that is done but not
used, and options that are offered but not wired.

Four results carry the audit:

1. **The `score` step does not change the oligo set.** Its `min_amp_pred=10.0`
   gate passed 2000 of 2000 candidates on tier M, removed a median of zero
   primers across 20 archived runs (maximum 4), and step 4 then reads only the
   primer column and discards the scores. The threshold equals the baseline of
   the hand-written rule the model was fit to, so a filter-surviving primer
   clears it by construction.
2. **The default optimizer does 38x the work for the same answer.** `hybrid`
   returned a set identical to `dominating-set` at every budget tested — the
   same primers, not merely the same coverage — while taking 347 s against
   8.95 s at 64 primers. Its Stage 2 optimises a network-connectivity term that
   measures 0.00 in almost every run on record, produces a worse set, and is
   then discarded by a guard.
3. **There is 3-7% more coverage available, and it costs seconds.** An exact
   max-coverage ILP solves this instance in 0.8-17 s. At 32 primers the default
   method spends 63 s to return a set the ILP beats by 5.45% in 3.3 s.
4. **Several documented configuration routes are silently inert.** Setting
   `coverage_reach` in params.json changes nothing while `--coverage-reach`
   changes the headline coverage figure by 2.5x; `num_primers` above 50 is
   rejected from params.json but accepted on the CLI; six flags advertise
   capability and are read by no module.

## Findings

Severity: **HIGH** = wrong result or hard failure on a documented path;
**MEDIUM** = misleading output or wasted work; **LOW** = cosmetic.

### F0 (HIGH) The `score` step does not change the oligo set

Three independent facts compose:

**Its gate passes everything.** `min_amp_pred` defaults to 10.0 and is applied
at `pipeline.py:1108` (`results["on.target.pred"] >= threshold`). On tier M the
step counted 2000 primers in and wrote 2000 out — its own log line reads
"Scoring 2000 primers... Filtered 0 primers based on efficacy". Across the 20
archived score runs on disk the gate removed 0 primers in twelve of them, 1 in
five and 4 in three. Its median effect is zero, and its maximum observed effect
is 4 primers out of 2000.

(The post-gate `step3_df.csv` cannot be used to measure this — it contains only
the primers that passed, so its minimum is above the threshold by construction.
The counts above come from the pre-gate input size and the step's own log.)

**Its scores are then discarded.** `unified_optimizer.py:598` reads
`candidates = step3_df["primer"].tolist()` — only the primer column. `amp_pred`
is never a ranking signal inside optimization.

**The threshold is the model's own baseline.** The model was fit to
`scripts/retrain_rf_model.py:208 compute_target_score`, a hand-written rule that
starts at `score = 10.0` and mostly adds for the properties the `filter` step has
already selected for (Tm 30-42 C, GC 0.4-0.6, GC clamp). A filter-surviving
primer scores above 10 essentially by construction, so the gate cannot bite.

So a whole pipeline stage — with a shipped `.skops` artifact, a SHA-256
allowlist, retraining scripts and its own documentation — has no effect on the
delivered set on any dataset measured here.

**The provenance is also documented two ways, and they contradict.**
`docs/training_amp_pred_on_real_data.md:3-4` says the shipped model "was trained
on synthetic data via `scripts/retrain_rf_model.py`". `rf_preprocessing.py:31-32`
says "Model trained on experimental SWGA primer sets with known amplification
outcomes." The training script settles it: `compute_target_score` is a rule, not
a measurement. The docstring should be corrected regardless of what else is done.

A related oddity: the rule the model was fit to is close kin to
`_heuristic_primer_score` (`rf_preprocessing.py:763`), the fallback used when the
model fails to load — same 10.0 baseline, same Tm and GC bands, different
coefficients. On the tier M pool the two rank primers only loosely alike
(Spearman 0.32, mean absolute difference 4.77 on a 0-20 scale), so they are not
interchangeable. Treat that figure as indicative only: it was computed on the
post-gate frame, whose restricted range depresses correlation. The point that
survives is structural rather than numerical — the shipped model and its
"fallback" are fits to variants of one hand-written rule, and neither reaches
the output.

**What to do.** This is a decision, not a bug fix. Either give `amp_pred` real
influence — feed it into candidate ranking in step 4, and set a threshold that
discriminates — or retire the stage. What should not persist is a stage that
looks load-bearing, costs a pipeline step, and changes nothing. Settling the
provenance question is the prerequisite for either path.

### F1 (HIGH) `coverage_reach` in params.json is ignored

Tier M, `dominating-set`, S=8, seed 42, everything else identical:

| route | reach used | fg_coverage |
|---|---:|---:|
| `"coverage_reach": 20000` in params.json | 3000 | 0.3816 |
| key absent entirely | 3000 | 0.3816 |
| `--coverage-reach 20000` on the CLI | 20000 | 0.9513 |

The params.json route gives output identical to omitting the key. The headline
coverage differs 2.5-fold between two routes the documentation presents as
equivalent: `coverage.py:236` and `reaction_conditions.py:95,172` both describe
the params.json route, and `coverage_reach` is a declared property in
`params.schema.json`. `unified_optimizer.py:669` reads it with
`getattr(parameter, "coverage_reach", None)`, and `parameter.py` never assigns
it. `additionalProperties: true` means nothing warns.

Three more keys have the same shape — read via `getattr` on a global
`get_params` never populates:

- `occupancy_ranking` (`pipeline.py:729`) and `occupancy_shortlist`
  (`pipeline.py:740`), both declared in the schema and reachable only by
  monkeypatch; `max_mismatches` (`pipeline.py:804`) likewise.
- `sampled_index_path` (`filter.py:276`) — and `filter.py:301`'s own error
  message instructs the user to set it in params.json, a route that does not
  exist.

### F2 (HIGH) `num_primers` above 50 is rejected from params.json, accepted on the CLI

`param_validator.py:95-96` caps `num_primers` and `target_set_size` at `(1, 50)`;
`params.schema.json` allows 200. The validator runs first:

```
params.json num_primers=64  ->  "Value 64 outside valid range [1, 50]"  (exit 1)
--num-primers 64            ->  runs, returns 64 primers
```

The 2026-08-19 audit raised the schema ceiling to 200 precisely because a >95%
design on a 3.2 Mb target needs 31-96 primers. That change did not reach the
validator, so the capability is unreachable through the documented configuration
file. The blocked range is the useful one: on tier M, 64 primers reach 0.661
coverage against 0.538 at 32, and 128 reach 0.810.

The error message also names neither the key nor the file — it prints the same
line twice for two different parameters.

### F3 (HIGH) `--enable-qa` crashes on all four steps

```
neoswga filter -j params.json --enable-qa
ERROR: Step 2 failed: module 'neoswga.core.pipeline_qa_integration'
       has no attribute 'run_step2_with_qa'          (cli/pipeline.py:322)
```

The module's real API is `apply_post_step2_qa_filter`, `combine_rf_qa_scores`
and `QAAwareOptimizer`. A live crash on a documented flag, not a no-op.

### F4 (HIGH) `--application` silently loses a third of its effect

`optimize_step4` always forwards `uniformity_weight` (argparse default `0.0`,
`unified_optimizer.py:1200,1246`), so the `kwargs.setdefault` at `:726` can never
install the profile's own value. From the CLI, `metagenomics` (0.25) and
`discovery` (0.20) both run at 0.0 while the log at `:729-732` reports the
profile as applied. `tm_weight` and `dimer_penalty` do take effect. Silent, and
it changes the design.

### F5 (HIGH) The default optimizer does 38x the work for an identical set

Tier M, identical params, same 2000-candidate pool, seed 42. Set equality is
exact (Jaccard 1.000), not merely equal coverage.

| S | dominating-set | hybrid | slowdown | same set | fg_coverage |
|---:|---:|---:|---:|:--:|---:|
| 6 | 7.99 s | 10.64 s | 1.3x | yes | 0.3688 |
| 12 | 8.28 s | 13.65 s | 1.7x | yes | 0.4116 |
| 32 | 8.17 s | 63.30 s | 7.8x | yes | 0.5376 |
| 64 | 8.95 s | 347.02 s | **38.8x** | yes | 0.6607 |
| 128 | 8.61 s | >25 min, did not finish | >175x | — | 0.8097 |

The S=128 hybrid run was still going after 25 minutes when this was written,
against `dominating-set`'s 8.61 s for the same budget — already past 175x. It is
reported as a lower bound rather than a measurement.

`dominating-set` is flat in set size while coverage climbs from 0.369 to 0.810.
The entire cost curve belongs to hybrid's Stage 2, and its exponent is rising:
1.28x over 6->12 (exponent 0.36), 4.64x over 12->32 (1.56), 5.48x over 32->64
(**2.45**). Extrapolating `S^2.45` to the schema ceiling of 200 puts a
default-method run near **95 minutes** for a set `dominating-set` returns in
about 9 seconds — and nothing interrupts it, because no optimizer has a
wall-clock budget (`docs/development/optimizers.md:343` documents
`OptimizerConfig.timeout_seconds`; the field does not exist).

**Why Stage 2 cannot help.** Captured at S=12:

```
Stage 1 target: 20 primers (coverage)    Stage 1 runtime: 0.17s
Stage 2 target: 12 primers (network)     Stage 2 runtime: 5.29s
Connectivity improved: 0.00 -> 0.00
Refinement covered 42.8% against 45.0% for the stage-1 prefix; keeping the prefix.
```

Stage 2 is not a silent no-op: it produces a worse set and a guard discards it.
The identical output is the guard working; the time is the cost of a refinement
that has not yet been better. Two structural reasons:

- **Its objective is degenerate.** Algebraic connectivity reads 0.00 in every
  run measured here and in **14 of the 17 archived optimize logs** on disk; the
  other three read 0.01. Across k=10/11/12, set sizes 8/16/32, two polymerases
  and both the chr21 and hg38 backgrounds, Stage 2 is maximising a quantity that
  is indistinguishable from zero.
- **Greedy set cover is nested.** Stage 1 takes `max(final+8, final*1.67)`
  primers and Stage 2 prunes back to `final`, so removing by coverage loss
  reverses the greedy order and lands on the prefix again.

This corroborates, on a second dataset, the note already in
`hybrid_optimizer.py:49-56` recording an MTB run where Stage 2 spent 150 s and
*reduced* coverage from 54.8% to 53.2%.

### F5b Three of the five methods return the same set; the one that differs is worse

Same pool, same seed, tier M. Jaccard is against `dominating-set`'s selection at
the same budget, so 1.000 means the identical primers, not merely equal
coverage.

All five methods at S=32, same pool, same seed:

| method | time | cost | fg_coverage | dimer_risk | selectivity_density | gap_gini | Jaccard vs ds |
|---|---:|---:|---:|---:|---:|---:|---:|
| dominating-set | 8.2 s | 1x | 0.5376 | 0.3468 | 31.7 | 0.726 | — |
| hybrid (default) | 63.3 s | 7.8x | 0.5376 | 0.3468 | 31.7 | 0.726 | 1.000 |
| background-aware | 133.3 s | 16.3x | 0.5376 | 0.3468 | 31.7 | 0.726 | 1.000 |
| network | 377.2 s | 46.2x | 0.4283 | 0.3327 | 36.3 | 0.886 | 0.208 |
| clique | 34.7 s | 4.2x | 0.4499 | **0.0000** | **42.4** | 0.827 | (different by design) |

The first three agree on **every metric**, not just coverage. At S=12 the same
holds: dominating-set 8.3 s, hybrid 13.7 s, background-aware 20.1 s, all at
0.4116 with Jaccard 1.000; network 64.9 s at 0.3644 (Jaccard 0.333); clique
13.1 s at 0.3037.

`background-aware` delegates to `HybridOptimizer(background_pruning=True)`
(`background_aware_optimizer.py:664-700`), so its agreement with hybrid is
structural rather than coincidental — but its agreement with plain
`dominating-set` means the background-pruning stage also changed nothing here,
while adding a further 2.1x on top of hybrid's cost at S=32. It is documented
for clinical use on the strength of a "10-20x background reduction"; on this
target it removed no primer.

`network` is the only method that genuinely searches elsewhere, and on this
dataset it is worse on both axes at once: 20% less coverage at 46x the cost at
S=32. The published guidance in `docs/optimization_guide.md` puts it at
"20-60s"; measured, it is 377 s.

**`clique` is the one method that is not dominated**, and the result suggests a
recommendation change. It is documented as the choice when a set must be
dimer-free, and it delivers that (`dimer_risk 0.0000` against 0.3468) — but it
*also* returns the best selectivity density of any method, 42.4 against 31.7,
in a quarter of `background-aware`'s runtime. `background-aware` is the method
currently recommended for clinical work on the strength of a "10-20x background
reduction", and here it returned a set indistinguishable from plain
`dominating-set`. **On this target, for specificity-led design, `clique`
outperforms the method documented for it on both specificity and dimer safety,
at 26% of the runtime**, giving up 16% coverage and some binding evenness
(gap_gini 0.827 against 0.726).

Two caveats belong with any use of `clique`: it saw 200 of the 2000 candidates
(`max_pool_size`), and it stopped at the `max_cliques` ceiling — its own message
reads "Evaluated 10000 dimer-free sets from 200 candidates" — so it is a
truncated search over a truncated pool, not an optimum. That its truncated
search still wins on specificity is itself worth following up.

The practical reading: on this target the five documented methods collapse to
**two distinct answers plus one strictly worse one**. Choosing among hybrid,
dominating-set and background-aware buys nothing but runtime; `network` is worse
on coverage and evenness for 46x the cost; and `clique` is the genuine
alternative. That may not generalise — a dataset where the amplification network
is not degenerate could separate the first three — but no such dataset has been
demonstrated, and the connectivity evidence below suggests they are rare.

### F6 How good is the oligo set? Within 3-7% of optimal, and far above random

Tier M, 2000 candidates, reach 3000, bin 750. All quantities measured on the
same bins in covered bases, so they are comparable to each other.

| S | greedy | ILP optimum | gap | LP bound | ILP solve | random median | random beating greedy |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 6 | 0.4048 | 0.4192 | 3.44% | 0.4192 | 0.8 s | 0.0298 | 0 / 100 |
| 12 | 0.4500 | 0.4817 | 6.59% | 0.4821 | 2.0 s | 0.0634 | 0 / 100 |
| 32 | 0.5833 | 0.6169 | 5.45% | 0.6170 | 3.3 s | 0.1612 | 0 / 100 |
| 64 | 0.7109 | 0.7664 | 7.24% | 0.7670 | 17.1 s | 0.2961 | 0 / 100 |

Every ILP was proven optimal, and the LP relaxation is tight to within 0.06% of
the integer optimum, so the LP is a usable bound wherever the ILP will not
finish. The gap widens with budget, which is expected of a greedy submodular
cover and means the headroom is largest exactly where the tool was recently
extended to operate.

These figures are in the bin space the set cover selects on. They are not
directly comparable to the base-level `fg_coverage` the pipeline reports, which
is a different measure — see F8.

### F7 (MEDIUM) The shipped ILP cannot bound anything, and is miswired

`dominating_set_optimizer.py:500 optimize_ilp` has no callers and three faults:

1. **Wrong question.** It minimises primer *count* subject to covering every
   reachable bin (`:546-555`), so below the minimum cover it is infeasible and
   returns `None` with only a warning. On tier S, budgets of 4 and 6 return
   `None`; 8 and above solve. Every greedy optimizer answers the opposite
   question — best coverage within a fixed budget.
2. **Wrong argument type.** `Model(sense=minimize)` passes the `minimize` helper
   where mip expects the string `'MIN'`. It survives only by falling through to
   the default; the same mistake with `maximize` yields a model CBC reports as
   "Empty problem - 0 rows, 0 columns" with objective 0 and no error raised.
3. **Its solver is not installed.** `mip` sits in the `improved` extra, so the
   path could never have run as shipped.

### F8 (MEDIUM) Greedy coverage is a fraction of bins, not of bases

`optimize_greedy` divides covered bins by bin *count*, weighting every bin
equally regardless of the span it represents. On the 6157 bp plasmid at
`bin_size=500` that reads 0.1538 where the base-weighted value is 0.1624 — 5.3%
apart with no primer changed. The same split shows on tier M at S=64, which
reports `score: 0.5938` (bin fraction) beside `fg_coverage: 0.6607` (measured
bases). This belongs with the wider coverage-reconciliation work: the tree holds
at least 14 implementations of "coverage".

A related trap sits in the same file. `_total_bins` (`:281-297`) carries a
careful docstring explaining why coverage must be divided by the bins in the
*genome* rather than the bins some candidate happens to reach — "a pool touching
a tenth of the target and covering all of it scored 1.000". But
`BipartiteGraph.get_coverage_score` (`:234-240`) still divides by
`len(self.regions)`, the old reachable-pool denominator. It has no production
caller; the only callers are four tests in
`tests/test_dominating_set_optimizer.py:174-200`, which pin the superseded
semantics. Dead code retaining a fixed defect, with tests asserting the defect,
is a trap for whoever wires it up next. Delete it or align it.

### F9 (MEDIUM) The headline `score` is saturated, but reaches reporting only

`hybrid_optimizer.py:1348` sets `score = final_predicted_amplification`, capped
at `2**20` (`network_optimizer.py:82`). Every hybrid run reports exactly
`1048576.0` regardless of set size or quality. It does **not** reach selection:
the ensemble picks by `metrics.normalized_score(application=...)`
(`unified_optimizer.py:288-291`) and nothing compares two results by `.score`.
So the number in `step4_improved_df_summary.json` is uninformative, but no
design is chosen because of it.

### F10 (MEDIUM) Six flags advertise capability and change nothing

None crash; all log encouraging messages. `--use-gpu` logs "GPU acceleration
enabled" and `--use-enhanced-features` logs "Using enhanced 120+ feature model",
while no module reads either attribute.

| Flag | Set at | Readers |
|---|---|---|
| `--use-gpu` / `--no-gpu` / `--gpu-device` | `cli/_common.py:418-451` | none |
| `--use-enhanced-features` / `--enhanced-model-path` | `cli/pipeline.py:396-409` | none |
| `--use-background-filter` (optimize) | `cli/pipeline.py:1549` | none |
| `--background-bloom-path` / `--background-sampled-path` | `cli/pipeline.py:1562,1565` | none |
| `--max-optimization-time` | `cli/pipeline.py:1606` | none |
| `--fast-score` | `cli/pipeline.py:1450` | self-documented no-op |

Also `--max-extension` is honoured only inside the `--validate-simulation` block
(`cli/pipeline.py:984`) and never sizes the optimizer's amplification network;
and `--show-frontier`'s MOEA table cannot fire, because
`OptimizationResult.pareto_front` has no assignment anywhere and no MOEA
optimizer is registered. This is the class
`tests/test_design_options_have_effect.py` exists to prevent, and it covers none
of them.

The GPU module is worth a separate note: only `batch_binding_probability` is
genuinely vectorised, while `batch_calculate_tm` loops in Python writing
element-by-element into a CuPy array, which is slower than NumPy. CLAUDE.md
advertises a "10-100x speedup". That line should be corrected.

### F11 (MEDIUM) The optimizer documentation describes a different tool

`docs/optimization_guide.md:385-397` gives a "Performance Benchmarks" table
listing `greedy`, `genetic`, `moea` and `milp`. None of the four can be
selected: the CLI accepts only
`{hybrid, dominating-set, network, background-aware, clique, ensemble}`, and
`unified_optimizer._ensure_optimizers_registered` (`:103`, imports at
`:122-126`) pulls in five modules — background-aware, clique, dominating-set,
hybrid, network — none of them a genetic or multi-objective optimizer. Half the
rows in the table describe methods that were retired.

The guidance is also inverted by measurement. `neoswga optimize --help`
positions `dominating-set` as the choice for "speed-critical, large pools" and
`hybrid` as the default for "general use". On tier M, `dominating-set` returned
the same set as `hybrid` everywhere tested, faster at every size, and was the
*only* method that stayed flat as the set grew. On this evidence the two
descriptions should be swapped, with `hybrid` documented as the special case
rather than the default.

### F12 (LOW) Dead params keys ship in every example

`max_sets` and `iterations` appear in all nine shipped example and integration
`params.json` files and are read by nothing on any user path. `max_sets` is
range-checked by `param_validator.py:98` and written by `wizard.py:628`, which
makes both look load-bearing. Also dead: `selection_metric`, `retries`,
`drop_iterations`, `top_set_count`.

### Checked and cleared

Recorded so they are not re-investigated:

- **`bg_coverage: 0.0` on every tier M run** is consistent with the
  `max_bg_freq` gate, which selects primers with near-zero chr21 sites; the hg38
  run does report a nonzero value. Not a repeat of the silent-empty-scan bug.
- **The GC-adaptive polymerase override works.** `pipeline.py:495-508` preserves
  a user-set polymerase; a current run logs "Keeping user-specified polymerase
  'phi29'". The August validation logs that appear to show an override predate
  the 2026-08-19 fix (`011df74`) and must not be used as correctness baselines.
- **`long_primer_mode` is live** (`parameter.py:1039-1052`), gating the RF k-mer
  sampling rate.
- **`efficiency_predictor.py` is reachable** via `neoswga predict-efficiency`;
  `amplicon_network.py` is reachable but only under `--validate-simulation`.
- **The missing-position warning is well designed.** When position files are
  absent the tool says the coverage number "is meaningless, not low" rather than
  reporting zero. That is the right behaviour and caught an error in this
  audit's own harness.

## Regression pass over the March 2026 audits

Every checkable finding in `STEP1`-`STEP4_AUDIT_REPORT.md`,
`AUDIT_REPORT_v3.0.md` and `OPTIMIZER_UX_AUDIT.md` was verified against HEAD.
Those documents now carry a supersession header pointing here.

**The March audits have largely been acted on.** Fixed since: the dead
`filter_extra()`, the sort that kept the worst primers, the `min_tm/max_tm`
crash, the RF 0-20 vs 0-1 scale error in all three report consumers, the
`--seed` irreproducibility (now pinned by a golden snapshot test), the silently
ignored `--use-mechanistic-model`, `normalized_score=0.0` from the
dominating-set adapter, the Jarzynski-Sarkar loop penalty missing RT, the
missing dependency upper bounds, all four `pickle.load` calls, the absent audit
trail, and the CI gaps. The broken GA crossover, the unregistered MILP and the
whole 18-optimizer UX surface are **superseded** — those modules were deleted.

**Still open, ranked by how directly they change which oligos are selected:**

1. **The `score` step's provenance and inertness** — see F0. The March STEP3
   finding about delta-G features was addressed by *skipping* them, not
   retiring them: `fast_score=True` feeds 27 zeros into a 53-feature model
   (`rf_preprocessing.py:726-730`). Their measured combined importance on the
   shipped artifact is 2.13%, so the risk is bounded, but the honest fix is a
   26-feature retrain rather than zeros.
2. **`--auto-size` reads the un-adapted polymerase.** `cli/pipeline.py:625-627`
   uses `json_data.get("polymerase", "phi29")`, bypassing the GC-adaptive value
   installed at `pipeline.py:501`. On a GC-rich target with no explicit
   `polymerase` in params.json, the set size is recommended for phi29 at 30 C
   while the run is designed for equiphi29. This changes the number of primers
   delivered. It is the one surviving instance of a defect fixed everywhere else.
3. **`neoswga validate model` cannot fail.**
   `model_validation.py:186` sets `passed = True` and never clears it; its one
   check logs at INFO (`:217`) and `equiphi_at_42` (`:221`) is computed and
   discarded; `validate_binding_kinetics` (`:245`) has no conditional at all.
   The command reports a green mechanistic model unconditionally. This is the
   open finding most likely to mask a future regression in the chemistry that
   drives primer weighting.
4. **`StreamingPositionCache` hardcodes k=6-12** (`position_cache.py:612`). Any
   equiphi29 (12-18 bp) or bst (15-25 bp) design gets no positions from that
   path — the same silent-zero class as the two `2**31` bugs already in
   CLAUDE.md's Known Issues.
5. **Step-3 scores discarded by step 4** (`unified_optimizer.py:598`) — the
   mechanism behind F0.
6. **`MultiGenomeKmerCounter` has no subprocess timeout**
   (`kmer_counter.py:205,210`) while the main `run_jellyfish` path does. A hung
   jellyfish in a multi-genome run blocks indefinitely.
7. **The vendored `melting_temp.temp()` pre-filter ignores user salt**
   (`kmer_counter.py:536,565`, Na=10/Mg=20 fixed, with a deliberate GC bug).
   Much less severe than in March, because `filter_extra` now applies the
   correct Owczarzy/NN Tm as a real gate (`filter.py:523`) — but it is pinned to
   the shipped RF model (`melting_temp.py:8-15`), so it and F0 have to be
   retired together.

Also still open and worth recording: `except Exception` has grown from 21
occurrences to 134; there is no `--output-prefix` or file locking, so two
concurrent optimize runs clobber `step4_improved_df.csv`
(`unified_optimizer.py:1290`); and CLAUDE.md still advertises the "10-100x" GPU
speedup that `gpu_acceleration.py:110-112` now disclaims in its own docstring.

## SWOT

Each entry cites a measurement or a file.

### Strengths

- **The selection is genuinely good.** Within 3.4-7.2% of a proven optimum, and
  0 of 100 random sets beat it at any budget (F6).
- **`dominating-set` scales flat in set size** — 8-9 s from 6 to 128 primers
  while coverage rises from 0.369 to 0.810 (F5).
- **A calibrated reach.** The 3 kb selection default is fitted to a measured
  outcome at 3.0-6.2 kb (`docs/validation/reach_calibration.md`). No other SWGA
  designer fits this from data; swga 2.0 uses 70 kb, at which almost any set
  scores ~1.0.
- **A real wet-lab benchmark suite** — seven published datasets in
  `tests/validation/data/`, including 310 qPCR measurements over 284 primers.
- **An unusually honest validation record.** `docs/validation/` documents
  rejected ideas, contradictory evidence and uncalibrated thresholds rather than
  only successes, and the codebase carries guards (the Stage-2 prefix guard, the
  meaningless-coverage warning) that fail safe.
- **Capabilities no competitor has**: multi-genome pan-target design, BAM-driven
  gap analysis and iterative set expansion.

### Weaknesses

- **A whole pipeline stage changes nothing.** `score` gates on a threshold that
  removes a median of zero primers and produces values step 4 discards (F0). It
  is the most elaborate inert surface in the tool: a shipped model artifact, a
  checksum allowlist, retraining scripts and a documentation page, none of which
  reach the delivered oligos.
- **The one ML component is a fit to a hand-written rule, described as a fit to
  lab data.** `rf_preprocessing.py:31-32` claims experimental training;
  `docs/training_amp_pred_on_real_data.md:3-4` says synthetic; the training
  script's label function settles it (F0). Of everything here this is the
  finding most likely to matter to a reviewer or a collaborator.
- **The default method is the worst-value one** — 38.8x the cost of
  `dominating-set` for an identical set at 64 primers (F5).
- **A large inert surface**: six no-op flags, one crashing flag, four params
  keys read from a global nothing assigns, dead keys in every shipped example
  (F1, F3, F10, F11).
- **Two configuration routes that disagree.** `coverage_reach` and `num_primers`
  behave differently depending on whether they are set in params.json or on the
  CLI, in one case changing the headline number 2.5-fold (F1, F2).
- **At least 14 implementations of "coverage"**, including a fallback in
  `report/metrics.py:503` that estimates it as `n_primers x 30000 / fg_size`.
- **No runtime guardrail at all**: no wall-clock budget in any optimizer,
  `make benchmark` matches no files, `pytest-benchmark` is not a dependency, and
  the only automated check is a nightly 2 GB RSS assertion.
- **The dominant cost is single-threaded.** Position scanning was 838 s of a
  924 s filter step on the 3.1 Gb background, on one of 11 cores.

### Opportunities

- **A real primer-quality model is an open lane, not a repair.** Because
  `amp_pred` currently reaches nothing (F0), wiring it into step-4 ranking and
  fitting it to the seven wet-lab datasets already in `tests/validation/data/`
  would be a genuine capability gain rather than a fix — and the plumbing,
  artifact handling and threshold machinery already exist.
- **Take the 3-7%.** The exact solve costs 0.8-17 s at this scale and the LP
  bound is tight enough to report routinely in place of the saturated `2**20`.
- **Delete Stage 2, or gate it on a non-degenerate connectivity.** The cheapest
  large win in the tool: it removes the entire set-size cost curve at no
  measured cost in quality.
- **Parallelise the genome scan.** 91% of the dominant step on one core, over
  windows the chunking already proves independent.
- **Extend `test_design_options_have_effect.py` to every option.** The
  vary-the-option test already exists; the inert flags are the ones it does not
  cover.

### Threats

- **Silent divergence between documentation and behaviour.** The reach and set
  size cases both arise from a fix reaching one route and not the other, and
  `additionalProperties: true` guarantees no warning. This class recurs: it is
  the same shape as the two independent `2**31` bugs in Known Issues #5 and #7.
- **Breadth outrunning evidence.** Options and modules accumulate faster than
  tests and calibration. In `quality_thresholds.py`, `TM_RANGE` and
  `DIMER_RISK` are documented as operating points "not calibrated against an
  outcome" (`:99`), with recorded counter-evidence that the Tm term runs
  backwards on the Leishmania sets; `GINI` is calibrated to the published
  distribution but labelled descriptive rather than predictive. The project is
  candid about this, which is a strength — but it means several grades a user
  reads as verdicts are preferences.
- **External dependencies with sharp edges** — jellyfish as a hard requirement,
  and the pyahocorasick 2 Gb limit that the chunking works around but does not
  remove.

## Recommendations, in order of measured payoff

1. **Decide what the `score` step is for.** It currently changes nothing:
   gate passes everything, scores discarded downstream (F0). Either wire
   `amp_pred` into step-4 candidate ranking and calibrate a threshold that
   discriminates, or retire the stage. First settle the provenance
   contradiction between `docs/training_amp_pred_on_real_data.md:3-4` and
   `rf_preprocessing.py:31-32` — that is a docstring fix owed either way.
2. **Stop running Stage 2 by default.** Either make `dominating-set` the default
   or gate Stage 2 on a connectivity that is not zero. Removes a 38x cost for no
   measured loss (F5), and the same change retires `background-aware`'s extra
   16x, since it returned the identical set (F5b).
3. **Re-check the clinical recommendation.** `clique` beat `background-aware`
   on both selectivity density (42.4 vs 31.7) and dimer risk (0.0 vs 0.347) at
   26% of the runtime on this target (F5b). Confirm on a second dataset before
   changing the guidance, but the current recommendation is not supported by
   this measurement.
4. **Fix the two split configuration routes** — `coverage_reach` and
   `num_primers`. Both are small changes with large, silent consequences (F1,
   F2). Assign the params keys in `get_params` and align
   `param_validator.py:95-96` with the schema.
5. **Fix or remove `--enable-qa`.** It crashes (F3).
6. **Fix `--application`'s `uniformity_weight` shadowing** (F4).
7. **Replace `optimize_ilp` with the max-coverage formulation** and report the
   LP bound as a routine output, in place of the saturated `2**20` (F7, F9).
   Consider using the exact solve outright at bacterial scale, where it costs
   seconds.
8. **Delete the six inert flags**, or wire them, and extend
   `test_design_options_have_effect.py` to cover every option (F10).
9. **Add a wall-clock budget** — implement `OptimizerConfig.timeout_seconds` or
   remove the documentation that promises it.
10. **Parallelise the position scan**, guided by a filter-step measurement
    rather than by inspection.

## Method

Tiers: **S** = pcDNA plasmid (6157 bp) vs pLTR; **M** = Prevotella (3,168,282 bp)
vs human chr21 (46,709,983 bp), k=12, 2000 candidates from the pre-built
`step3_df.csv`, phi29 at 30 C, `--seed 42`.

Timings drive the CLI rather than the library, because that is what users run:
a library harness would miss the dispatch, the position-cache build and the
metric computation that dominate a real invocation. The existing
`scripts/benchmarking/benchmark_suite.py` was not reused — it imports
`neoswga.core.milp_optimizer` and `neoswga.core.moea_optimizer`, both deleted,
and its axes no longer match the shipped method set. It and its two siblings are
left in place but marked stale in the directory's README; removing them is a
separate call.

The optimality bound is a max-coverage ILP built over the same `BipartiteGraph`
bins the greedy optimizer sees, with the same `extension_reach`, so it bounds
the quantity actually being optimised:

```
maximise    sum_b w_b * y_b
subject to  y_b <= sum_{i covers b} x_i     for every bin b
            sum_i x_i <= S
```

with `w_b` the bases the bin spans, so the objective is covered bases rather
than covered bins.

**Caveats on the numbers.** Timings are from a single run each. One
configuration was run twice by accident (`dominating-set` at S=6, 8.0 s and
8.7 s), which puts the repeat noise around 9%, so treat any difference under
about 10% as nothing. The effects reported here — 7.8x, 16.3x, 38.8x, 46.2x —
are orders of magnitude outside that. The S=128 hybrid run overlapped a second
sweep for part of its life and had not finished after 25 minutes, so it is
reported as a lower bound on the slowdown rather than as a timing. Peak RSS is measured on the child
process and includes the position-cache build, which is why every run reads
~1.35-1.5 GB regardless of method; it is not a discriminating figure here.

## Not yet measured

- `count-kmers`, `filter` and `score` sweeps across k and genome size, including
  the exact / Bloom / sampled-index background paths.
- Whether the score step's fast path changes primer ranking, which decides
  whether the delta-G features should be retired rather than merely skipped.
- Tier L (3.1 Gb hg38 background) confirmation of the optimize results.
- The `ensemble` method, which was not run. Its selection criterion
  (`normalized_score`) is understood, but since three of the four methods it
  runs by default return the identical set here, it is worth confirming what it
  reports in a comparison table where the runners-up are indistinguishable.
- The published-set comparison against the seven wet-lab fixtures.
- Two findings the regression pass could not settle by inspection: whether a
  genome named in both `params.json` and `--blacklist` is counted twice, and
  whether the replication simulator's inter-binding intervals are Poisson. Both
  need a run rather than a read.
