# Plan: completing the alternatives, scaling and oligo-set-quality audit

**Companion to** [`AUDIT_2026-08_alternatives_and_scaling.md`](AUDIT_2026-08_alternatives_and_scaling.md),
which is the interim report. That document answers the question for the
`optimize` step and for oligo-set quality against a proven optimum, on one
target. This plan covers what it could not settle: the three earlier pipeline
steps, a second dataset, the scaling tier that matters for host-background work,
and the improvement lane.

**Ground rule carried over from the audit.** Every optimizer conclusion on
record rests on one target (Prevotella vs human chr21) at one k. That is enough
to prove a defect exists but not enough to change a shipped default. Any phase
below that ends in a default change requires a second target first, and that
requirement is written into the phase gate rather than left to judgement.

---

## Phase A - close the measurement gaps

The audit measured `optimize` in full and left the other three steps as
inspection only. Nothing in Phase A changes behaviour; it produces the evidence
the later phases are gated on.

### A1 Scaling sweep of count-kmers, filter and score

**Question.** How does each step scale in genome size and in k, and which of the
three background paths (exact jellyfish, Bloom filter, sampled index) is
actually the right default at each size?

**Method.** Extend `scripts/benchmarking/sweep_optimize.py` into a
`sweep_pipeline.py` that drives the CLI, not the library, for the same reason
the optimize sweep did: dispatch, cache build and metric computation dominate a
real invocation. Axes:

- k in {8, 10, 12, 14, 16} - spans phi29 (6-12) and equiphi29 (12-18) regimes.
- Target: tier S (6157 bp plasmid), tier M (3.17 Mb Prevotella).
- Background: none, chr21 (46.7 Mb), hg38 (3.1 Gb).
- Background path: exact, `build-filter` Bloom, sampled index.

Record wall time per step, peak RSS on the child process, and the per-stage
funnel counts already written to `filter_stats.json`.

**Two things this is expected to expose.**

1. The single-threaded position scan. It was 838 s of a 924 s filter step on the
   3.1 Gb background, on one of 11 cores. The sweep should establish whether
   that fraction holds across k, which decides whether parallelising it is the
   single highest-value performance change or only a hg38-scale one.
2. Whether the Bloom and sampled-index paths earn their complexity. CLAUDE.md
   already records that exact counting of hg38 at k=12 costs about 7 minutes and
   138 MB, which is well inside reach, and that the sampled path has a
   resolution trap (`_warn_if_sample_too_sparse`). If exact counting wins up to
   k=14, the recommendation in the docs should narrow to where it does not.

**Output.** A scaling table per step, and a statement of the cost model for each
(linear in genome size, in candidates, in k).

### A2 Tier L confirmation of the optimize results

**Question.** Do F5 (hybrid returns the identical set at 38x-260x the cost) and
F5b (three of five methods agree exactly) hold against a 3.1 Gb background?

Run the five methods at S in {12, 32} against hg38. This is the tier where the
two `2**31` bugs in Known Issues #5 and #7 lived, so it is also a standing
regression check on the chunked scan and the int64 positions.

**Cost.** Hybrid at S=32 took 63 s on chr21; the position cache build dominates
at tier L, so budget a working day of machine time and run it detached.

### A3 A second target

The gate for every default change below. Pick a target that differs from
Prevotella on the axis that matters: an organism whose amplification network is
not degenerate, since Stage 2's objective read 0.00 in 14 of 17 archived logs
and the whole F5 argument is that it is maximising nothing. A GC-extreme target
(the Wolbachia fixtures, 35.2% GC) is the natural second case because it also
exercises the GC-adaptive strategy.

If algebraic connectivity is non-zero on some target, Stage 2 becomes a
conditional rather than dead weight, and recommendation 2 changes shape. If it
is zero there too, the case for deleting it is closed.

### A4 The `ensemble` method, and the published-set comparison

Two smaller items the audit lists as not run.

- `ensemble` was never executed. Its selection is by `normalized_score`, but
  three of the four methods it runs by default return the identical set on tier
  M, so its comparison table should be inspected for what it reports when the
  runners-up are indistinguishable.
- The seven wet-lab datasets in `tests/validation/data/` (including 310 qPCR
  measurements over 284 primers) have never been used to compare a neoswga set
  against a published one. This is the only measurement in the plan that speaks
  to whether the tool designs good primers rather than good coverage numbers,
  and it is the prerequisite for C4.

---

## Phase B - the decisions that change which oligos are delivered

Each of these is a decision, not a bug fix, and each is already supported by
measurement. They are ordered by how directly they change the output set.

### B1 Decide what the `score` step is for (F0)

The step gates on a threshold that removed a median of zero primers across 20
archived runs, and step 4 then reads only the primer column
(`unified_optimizer.py:598`). Two paths:

- **Wire it.** Feed `amp_pred` into step-4 candidate ranking and calibrate a
  threshold that discriminates. This is the honest version only if the model is
  refit to something other than the hand-written rule it currently reproduces -
  see C4, which is the same work.
- **Retire it.** Remove the stage, the `.skops` artifact, the checksum allowlist
  and the retraining scripts. Note the coupling recorded in the audit's
  regression pass: the vendored `melting_temp.temp()` pre-filter is pinned to
  the shipped RF model, so it and the score step have to be retired together.

**Owed either way, and immediately:** the provenance contradiction.
`rf_preprocessing.py:31-32` claimed training on experimental outcomes;
`docs/training_amp_pred_on_real_data.md:3-4` said synthetic. The training
script's `compute_target_score` settles it - it is a rule, over features drawn
from `np.random.multinomial`.

> **Done.** Both places now say synthetic, name the labelling rule, and record
> why the default gate cannot bite. The decision above is untouched by this: the
> model is described correctly now, and still reaches no delivered oligo.

### B2 Stop running Stage 2 by default (F5)

**Gated on A3.** If connectivity is degenerate on the second target too, either
make `dominating-set` the default or gate Stage 2 on a non-zero connectivity.
This removes the entire set-size cost curve - 38.8x at 64 primers, 260x at 128 -
at no measured cost in quality, and retires `background-aware`'s further 16.3x
along with it, since it returned the identical set.

Prefer the gate over deletion if A3 finds any non-degenerate case: a conditional
that fires rarely is cheaper to justify than removing a documented capability.

### B3 Re-check the clinical recommendation (F5b)

**Gated on A3.** On tier M, `clique` beat `background-aware` - the method
currently documented for clinical work - on selectivity density (42.4 vs 31.7)
and on dimer risk (0.0000 vs 0.3468) at 26% of the runtime, giving up 16%
coverage. Two caveats belong with any use of that result: `clique` saw 200 of
2000 candidates and stopped at its `max_cliques` ceiling, so it is a truncated
search over a truncated pool. Confirm on the second target, then either change
the guidance or record why not.

---

## Phase C - improved methods for the oligo set

This is the "could we do better" lane. The audit establishes the headroom
exactly: greedy lands 3.4-7.2% below a proven exact optimum, and the gap widens
with budget, so the headroom is largest where the tool was recently extended to
operate.

### C1 Replace `optimize_ilp` with the max-coverage formulation (F7) - DONE, except the wiring

> **Done.** `optimize_ilp` now solves the max-coverage program; a fourth fault
> found while fixing it (the reach was never passed into the graph, so a
> non-zero reach *lowered* reported coverage) is fixed with it, and
> `coverage_upper_bound` adds the LP relaxation. Pinned by
> `tests/test_dominating_set_max_coverage_ilp.py`. **Still open:** it has no
> production caller, so the boundary measurement below is what remains.

The shipped ILP had no callers and three faults: it minimised primer count
subject to full coverage (the opposite question), passed the `minimize` helper
where mip expects `'MIN'`, and its solver ships only in the `improved` extra.

The correct formulation was already written and measured in
`scripts/benchmarking/max_coverage_bound.py`:

```
maximise    sum_b w_b * y_b
subject to  y_b <= sum_{i covers b} x_i     for every bin b
            sum_i x_i <= S
```

It solved tier M in 0.8-17 s and every instance was proven optimal. At bacterial
scale this is not a bound to report, it is a method to use - taking the full
3-7% for seconds of runtime. Report the LP relaxation (tight to within 0.06%)
wherever the ILP will not finish, and put it where the saturated `2**20` score
currently sits (F9).

Decide the boundary by measurement: find the genome size and candidate count at
which the ILP stops finishing in reasonable time, and dispatch to greedy above
it.

### C2 Cheap improvements to greedy, for the range the ILP cannot reach

If C1 establishes a ceiling above which exact solving is out, two standard
options are worth measuring against the same bound:

- **Lazy-greedy (accelerated greedy).** Exploits submodularity to skip marginal
  gain recomputation. Identical output to greedy, materially faster - relevant
  because it makes the fixed budget go further, not because it improves quality.
- **Local search / GRASP over the greedy solution.** Swap-based improvement on
  the greedy set. Known to close a meaningful part of the greedy gap on
  max-coverage. Measure the closed fraction against the ILP optimum directly -
  `optimality_gap.py` already computes exactly this comparison.

Neither is worth implementing before C1 shows where the exact method runs out.

### C3 Lift `clique`'s two truncations

`clique` produced the best selectivity density in the audit while seeing 10% of
the candidate pool and stopping at its clique ceiling. That a truncated search
wins on specificity is the strongest hint in the audit that the method is
under-exploited. Measure quality against `max_pool_size` and `max_cliques` to
see whether the win grows or was luck.

### C4 A primer-quality model fit to the wet-lab data

Depends on A4. Because `amp_pred` currently reaches nothing, fitting a model to
the seven datasets in `tests/validation/data/` and wiring it into step-4 ranking
is a capability gain rather than a repair, and the plumbing, artifact handling
and threshold machinery already exist. Its success criterion is the A4
comparison: a set ranked with the new model should beat the current set on the
qPCR measurements, or the model does not ship.

---

## Phase D - the inert surface

Not scientific findings, but each is a documented capability that does not work,
and together they are the largest source of the "documentation describes a
different tool" class.

| Item | Finding | Status |
|---|---|---|
| ~~`--enable-qa` crashes on all four steps~~ | F3 | **Done.** Wired to the API that exists, through three wrappers owning the CSV I/O. Two further defects fixed with it: `combine_rf_qa_scores` read a column step3_df.csv does not have, and `parameter.enable_qa` was never cleared, so the flag leaked across steps in one process. |
| ~~`--application` loses `uniformity_weight`~~ | F4 | **Done.** `None` sentinel from argparse to `_resolve_selection_weights`; an explicit `0.0` stays a real choice. Changes which primers five application profiles select. |
| ~~`occupancy_ranking`, `occupancy_shortlist`, `max_mismatches`, `sampled_index_path`~~ | F1 | **Done.** Fixed as `coverage_reach` was, in one extracted `_apply_params_only_keys`; the two keys missing from the schema are now declared. |
| ~~Six no-op flags~~ | F10 | **Done, without deleting any.** `--max-extension` wired; the rest report honestly instead of logging success. Which five to delete outright is left to the owner. |
| ~~CLAUDE.md's "10-100x" GPU claim~~ | F10 | **Done.** Removed from CLAUDE.md and from the help text. |
| `optimization_method` in params.json is ignored | **F1b** | **Open, and new.** Found while assessing the workflow, after the sweep above. Same defect on the key that decides which optimizer runs, so it carries F5's cost: up to 260x for an identical set. Needs the F4 sentinel, since `'hybrid'` is both the default and a valid explicit choice. |
| No wall-clock budget in any optimizer | F5 | Documentation half done (`optimizers.md:343` corrected, `--max-optimization-time` no longer advertises an unenforced default). No optimizer checks a clock. |
| Two vacuous tests in `test_design_options_have_effect.py` | F10 | `test_removed_flags_are_not_offered` and `test_scoring_weights_is_not_offered` assert against `create_parser().format_help()`, which lists subcommand names only, so neither inspects the flags it guards. |
| `get_coverage_score` dead code with tests pinning the superseded denominator | F8 | Delete or align; tests currently assert the defect. |
| Retired methods in the benchmark table | F11 | `docs/optimization_guide.md` lists `greedy`, `genetic`, `moea`, `milp`; none are selectable. |

`sampled_index_path` was worth doing early regardless of its severity, because
`filter.py:301`'s own error message instructs the user to set it in params.json
- a route that did not exist until now, and one A1 exercises.

**What Phase D taught, beyond the individual fixes.** The sweep above closed
every known instance of the inert-option class, and F1b was found immediately
afterwards by asking a different question — how the pipeline behaves as a
workflow rather than what each step's options do. That suggests the class is not
exhausted by enumerating flags, and that the vary-the-option test is the durable
fix rather than any particular repair. Pointing the two vacuous tests at the
real help text is the cheapest step toward that.

---

## SWOT

The full SWOT is in the audit, with each entry cited to a measurement or a file.
It is not restated here. Two notes on how this plan interacts with it:

- **The SWOT is provisional in one direction.** Every weakness about the
  optimizer rests on one target. A3 is the phase that either confirms them or
  narrows them to a dataset property, and the SWOT should be revised after it -
  in particular "the default method is the worst-value one", which is the
  strongest claim in the document and the one most dependent on connectivity
  being degenerate.
- **Phase A adds no SWOT entries by design.** It measures three steps that have
  so far been assessed only by inspection, so entries about the filter and
  count-kmers steps should be added only after A1, not before.

---

## Sequencing and effort

```
A1 sweep ---------+
A2 tier L --------+---> revise SWOT ---> B2, B3 (default changes)
A3 second target -+
A4 wet-lab + ensemble ---> C4

B1 provenance docstring   (do now, independent)
D items                   (independent, any time)
C1 ILP                    (independent of A; gated only on its own boundary measurement)
C2, C3                    (after C1)
```

Rough machine cost: A1 is the largest, dominated by hg38 filter runs at several
k; A2 is a working day detached; A3 is a tier-M-scale rerun. C1 is
implementation rather than measurement, and the formulation is already written
and validated in the harness.

The independent items - B1's docstring correction, the Phase D fixes, and C1 -
need no further measurement and can proceed while Phase A runs.

## Reproducibility

All measurement uses `scripts/benchmarking/` and drives the CLI with `--seed 42`.
Repeat noise on this machine is about 9%, established from one accidental
duplicate run, so treat any difference under 10% as nothing. Record the
environment line the audit uses (commit, Python, OS, cores, RAM, jellyfish
version, start method) with every new table.
