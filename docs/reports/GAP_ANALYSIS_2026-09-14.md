# Gap analysis, 14 September 2026

Snapshot taken at 07:20:51 local time against `main` at `298534a`, with the
working tree as it stood then. Every figure below was measured on that tree;
the command or the file and line is given wherever a claim is made.

**The tree was being edited by another process while this analysis ran.**
Source files were written at 07:15:31, 07:17:19, 07:19:48 and 07:19:50, and a
`pytest tests/ -n 8` that was not part of this analysis started at 07:20:04.
The working-tree fingerprint at snapshot time was `291123c`
(`git diff | shasum`). Findings are a snapshot, not a description of a settled
tree.

Scope: plan execution against code, test coverage, packaging and release
readiness, configuration and flag routing, scientific validation, dead code, and
the evidence base.

---

## Summary

The largest gap is not in the code. A substantial body of optimizer work exists
in exactly one working tree, has never reached CI, fails nine of its own tests,
and cannot be committed as it stands without breaking the command-line entry
point.

| Section | Gap | Severity |
|---|---|---|
| 1 | Committing the tracked edits without the new files breaks every command | Critical |
| 2 | `--num-primers` has become a request rather than a guarantee | High, decision needed |
| 3 | No test relates a delivered panel to the configured dimer limit | High |
| 4 | The shared dimer policy module was never built; nine files keep the looser default | High |
| 5 | The quantitative evidence base is absent from this machine and unversioned | High |
| 6 | Three untested paths can each return a confident wrong answer | High |
| 8 | Six config keys and four flags are accepted and never reach behaviour | High |
| 9 | The headline score of the default method is never compared to the measured data in the repo | High |
| 12 | The documented output columns and JSON keys do not match what is written | Medium |
| 7 | Plans and audits record execution state that is wrong or absent | Medium |
| 10 | Roughly 1,800 lines of unreachable code, some documented as live | Medium |
| 11 | Quality gates advisory only; 2,793 lint findings, 182 type errors | Low |

---

## 1. The entry point does not import without the untracked files

`neoswga/cli_unified.py` is tracked and modified. It registers a `plan-pool`
command and imports `neoswga.cli.plan_pool` at module scope
(`neoswga/cli_unified.py:269`). That module, together with
`neoswga/core/pool_planner.py` and `neoswga/core/pool_plan_report.py`, is
untracked.

Reproduced by copying only the files git tracks into a scratch directory and
importing the entry point:

```
IMPORT FAILED: ModuleNotFoundError No module named 'neoswga.cli.plan_pool'
```

Committing the tracked modifications without `git add` on the new modules
therefore breaks `neoswga --help`, not only the new command. This is an
all-or-nothing commit.

A second instance is milder. `neoswga/core/hybrid_optimizer.py:927` imports
`.swap_refinement`, also untracked, but the import sits inside `_swap_refine`
and that method is reached only when `refinement_method` is `swap`, which is not
the default. That one breaks a flag rather than the program.

Untracked source and test files at snapshot time:

```
neoswga/cli/plan_pool.py
neoswga/core/pool_planner.py
neoswga/core/pool_plan_report.py
neoswga/core/swap_refinement.py
tests/test_dimer_policy.py
tests/test_pool_planner.py
tests/test_swap_refinement.py
```

`main` is level with `origin/main`, zero commits either way, so none of the 29
modified tracked files has been through CI.

## 2. The set-size contract has been reversed

`python -m pytest tests/ -n 8`, run three times across two machines' worth of
sessions, agrees: 9 failed, 4360 passed, 75 skipped, out of 4442 collected.

| Failing tests | Cause |
|---|---|
| 4 set-size guards | Dimer relaxation now defaults off, so the greedy stops short of the requested count |
| 3 module-size ratchets | `hybrid_optimizer.py` 1641 lines against a 1600 budget, `unified_optimizer.py` 1606 |
| 1 function-length ratchet | `hybrid_optimizer.__init__` 205 lines against 200, `optimize` 356 against 350 |
| 1 shipped-config check | `examples/wolbachia_pool_design/params.json` carries an unknown key, `verbose` |

The cause of the four is a single default flip.
`dominating_set_optimizer.py` previously hardcoded
`self.relax_dimer_constraint_when_stuck = True`; it now reads
`= allow_dimer_relaxation`, a new configuration field defaulting to false.
Causality was confirmed directly by forcing that field true at run time, after
which all four pass.

Both failing test files are unmodified in the working tree. They encode the
previous contract, that a request for N primers delivers N.
`tests/test_design_options_have_effect.py` describes itself as fixing the
headline defect that `--num-primers` was not honoured, under the commit message
"fix: stop the inert CLI flags reporting success". The new default reverses
that: a request for 8 delivers 4 and a request for 10 delivers 6. The panel is
marked partial and a warning is logged, so it is not silent, but the flag has
gone from a guarantee to a request.

The two contracts are mutually exclusive on the shipped pools. At threshold 3
those pools support 29, 31 and 26 primers for S. aureus, E. coli and
M. tuberculosis, against shipped panels of 96 to 160.

**The safety argument for strictness was not built.** The plan justified
strict-by-default with a backtrack-before-relax mechanism, which would try
alternative earlier choices before either relaxing or stopping. No such
mechanism exists anywhere in the tree. Strictness currently means stopping, not
searching harder.

This is a product decision and it blocks the commit.

## 3. Nothing checks a delivered panel against its dimer limit

This is the missing test that turned the section above into a conflict rather
than a fix. Both halves are tested apart and never together: the validator
against a hand-built primer list
(`tests/test_delivered_pool_is_screened_for_dimers.py:120` writes the validation
JSON by hand), and the optimizers without the validator. The production call
site is `neoswga/core/unified_optimizer.py:1249`. The documented 11 bp delivered
violation against a configured 3 appears only in docstrings.

Had that test existed it would have been failing all along, and the default flip
would read as a fix. As it stands the suite can tell you that the behaviour
changed, not which behaviour is wanted. Write this test before resolving the
four failures, because it is the artifact that records the decision.

## 4. The shared dimer policy module was never built

Of the nine tasks in the 2026-09-12 optimizer plan, four are done, two partial
and three not started. Every checkbox in every plan file is unchecked, including
those whose work shipped, so all states below were verified against code.

| # | Task | State | Evidence |
|---|---|---|---|
| 1 | Two metric edge cases | Done, by a different mechanism | `base_optimizer.py:226`, `:386`, `:1407` |
| 2 | Shared dimer-constraint policy module | Not done | `core/dimer_constraint_policy.py` absent |
| 3 | Strict by default, backtrack before relax | Partial | strict at `dominating_set_optimizer.py:331`; no backtrack |
| 4 | Route three optimizers through the resolver | Partial | relaxation routed, threshold untouched |
| 5 | Structured relaxation events | Not done | no `relaxation_events` symbol anywhere |
| 6 | Covered-base-weighted Stage-1 selection | Not started | `dominating_set_optimizer.py:446`, `:458` still count bins |
| 7 | MILP dimer and fixed-primer constraints | Done | `dominating_set_optimizer.py:1001`, `:1007` |
| 8 | Bounded swap refinement behind a flag | Done | `swap_refinement.py:19`, `hybrid_optimizer.py:816`, `:925` |
| 9 | Benchmark swap refinement, decide the default | Done, different artifacts | `docs/validation/refinement_real_pools_2026-09.md` |

Task 2 is the largest omission. Its stated purpose was to make 3 the single
canonical `max_dimer_bp` default. Every scattered 4 it named survives verbatim
across nine files, including `base_optimizer.py:682`,
`unified_optimizer.py:567`, `network_optimizer.py:671`,
`hybrid_optimizer.py:1573`, `background_aware_optimizer.py:158`,
`optimizer_factory.py:375`, `clique_optimizer.py:85`, `search_context.py:135`
and `improved_pipeline.py:58`. A caller who omits the key still gets a looser
threshold than the documentation states. What shipped instead is a narrower
relaxation boolean covering half the policy.

Task 1 carries a fix the plan did not ask for: foreground positions are now
keyed per prefix, so effective coverage no longer pools multi-target
coordinates.

Task 9's measurement supports keeping `network` as the default. At a
10,000-evaluation budget swaps made zero replacements in all six cases; at
100,000 coverage rose in five of six while host binding rose in both 10-mer
cases.

## 5. The evidence base is absent

`runs/` does not exist on this machine and is gitignored. The three GC-tier
reference designs that most quantitative claims in CLAUDE.md rest on, including
the dimer sweep, the background-aware comparison at n=12, 24 and 36, and the
redundancy-threshold measurements, were never versioned and are gone.

The uncommitted swap-refinement evaluation records the same problem in its own
method section: "The earlier GC-tier run directories were absent." It fell back
to three candidate pools on a single target, Prevotella against human
chromosome 21, which it states plainly is "three pools on one target, not three
independent species."

Re-deriving is feasible. `hg38.analysisSet.fa.gz` and a human 12-mer count table
are present under `/Volumes/sekvens2/neoswga/`.

## 6. Three untested paths that can return a confident wrong answer

Each repeats a failure mode the project has already been bitten by, where a
measurement silently reads as zero or as truncated rather than raising.

- **The exclusion-genome filter stage has no test at all.** It removes
  candidates from the pool the optimizer sees and writes the
  `after_exclusion_blacklist` funnel row. Defined at
  `neoswga/core/pipeline.py:25`, applied at `:1138-1143`. No test references
  `filter_by_exclusion_genome`, `excl_genomes`, `excl_prefixes` or
  `excl_threshold`. The single funnel mention,
  `tests/test_filter_funnel_stages.py:126`, is a literal in a report-rendering
  assertion.
- **The Known Issue 13 fix went one class deep.**
  `StreamingPositionCache.get_positions` still returns an empty array for a
  prefix it holds no file for, silently
  (`neoswga/core/position_cache.py:713-746`, return at `:730`), while its
  sibling raises. `unified_optimizer.py:492-495` chooses between them on a flag.
  The AST guard at `tests/test_expansion_counts_background.py:210` matches only
  the literal name `PositionCache` and so skips the streaming call site.
- **The HDF5 equal-length write branch inherits an existing dataset's dtype.**
  An older int32 dataset plus a rescan yielding the same site count bakes
  coordinate truncation into the file (`neoswga/core/string_search.py:307`).
  The far-coordinate round-trip test takes the create branch instead
  (`tests/test_position_cache.py:508` against `string_search.py:304`).

## 7. Plans and audits misreport their own state

`docs/reports/AUDIT_pipeline_four_steps_2026-09-06.md:1204` closes with
"Nothing in this list has been implemented. The audit changed no code." That is
now false: 17 of its 32 findings are done. It is the most misleading sentence in
the documentation set.

`.gitignore:174` ignores `docs/superpowers/`. Six plan documents live only on
this machine and hold 349 unchecked boxes between them:

| Plan | Open boxes |
|---|---:|
| 2026-09-06 observability and startup | 98 |
| 2026-09-06 candidate pool quality | 71 |
| 2026-09-12 optimizer review fixes | 67 |
| 2026-09-06 optimizer cost and dimer criterion | 44 |
| 2026-09-06 correctness, silent wrong answers | 39 |
| 2026-09-06 filter throughput and memory | 30 |

Not one box is ticked in any of them although much of the work landed, so the
plans cannot be used to tell what remains.

Sixteen items were verified still open against code, including a circular genome
held twice during the scan (`string_search.py:202`), a memory guard whose
exception is raised nowhere (`exceptions.py:334`), two prerequisite validators
with no call site or no checks (`pipeline.py:128`, `:260`), three filter keys
read in `filter.py` but never assigned a global in `parameter.py`, an
efficiency predictor that counts a 10 kb bin as covered when a single primer
binds anywhere in it, with no extension reach applied, so its figure is not
comparable to the base-by-base coverage the summary reports
(`efficiency_predictor.py:488`), and a max-coverage ILP with no production
caller (`dominating_set_optimizer.py:1070`).

`tests/test_12bp_pipeline.py` collects 16 tests and skips all 16, because the
fixture directory `tests/integration/wolbachia_e2e/12bp/` neither exists nor is
tracked.

## 8. The inert-option class is still live in ten places

CLAUDE.md states that Known Issue 8 was "the last known instance of one class:
a config key or flag that is documented, accepted, and read by nothing." Ten
instances remain, in three shapes.

**Five keys never bind a module global**, so every reader takes its fallback.
Confirmed by execution: `hasattr` on the parameter module returns false for
`mismatch_penalty`, `retries`, `drop_iterations`, `top_set_count` and
`selection_metric`. The last four appear only as entries in a defaults
dictionary at `neoswga/core/pipeline.py:454-461` that nothing consumes.

`mismatch_penalty` is the sharpest case. A consumer was written for it,
`occupancy.default_mismatch_penalty` at `neoswga/core/occupancy.py:60`, whose
docstring states "This is its first consumer, and its intended meaning." It
reads `getattr(parameter, "mismatch_penalty", None)`, which is always None, so
it always returns the hardcoded 4.0. The consumer exists and the configured
value still never arrives.

`bl_penalty` binds a global, is assignable from params.json and from
`--bl-penalty`, and is range-validated at `neoswga/core/param_validator.py:106`,
but no scoring code reads it.

**Three argparse defaults beat params.json**, the exact mechanism Known Issue 8
describes. The fix used elsewhere is a `None` sentinel.

- `filter --gc-tolerance` defaults to 0.15 at `neoswga/cli/pipeline.py:1322`,
  so the block at `:201-208` always fires and recomputes `gc_min` and `gc_max`
  from the genome GC, overwriting whatever params.json configured.
- `plan-pool --swap-max-evaluations` defaults to 100,000 and is passed straight
  through, while the schema and module default is 10,000, so the two disagree.
- `filter --excl-threshold` assigns its default 0 unconditionally whenever an
  exclusion genome is given (`neoswga/cli/pipeline.py:281-283`).

**`expand-primers --optimization-method` was never fixed.** Known Issue 8 says
the argparse default is now the `None` sentinel in three places. On this
subcommand it still defaults to `"hybrid"` (`neoswga/cli/iterate.py:887`) and is
forwarded directly at `:225`, bypassing `optimization_method_from_params`. A
configured method cannot reach `expand-primers`.

`analyze-set --simulate` is accepted, never read, and absent from the
`UNIMPLEMENTED_OPTIONS` registry, so unlike the six honest inert flags it warns
about nothing.

`long_primer_mode` is partially routed. Its comment claims it enables k-mer
sampling, GPU acceleration and memory optimisations; its only effect is
selecting a `sample_rate` default.

Two smaller inconsistencies. The user guide's example spells `top_set_count` as
`top_sets_count` (`docs/guides/user-guide.md:417`), which is silently accepted.
The root `VERSION` file reads 3.0.0 while `pyproject.toml` and
`neoswga/__init__.py` read 3.6.0, and nothing reads the file.

The four new optimizer flags reach `optimize` but not `design`, which is the
same asymmetry Known Issue 8 recorded for `--optimization-method`.

**The guard does not cover any of them.** CLAUDE.md names three tests as the
ones that hold this line: `tests/test_design_options_have_effect.py`,
`tests/test_params_json_routes_optional_keys.py` and
`tests/test_optimizer_config_reaches_optimizers.py`. Between them they exercise
roughly fifty keys. None of the thirteen above appears in any of the three:
not the six inert keys, not the three argparse overrides, not
`expand-primers --optimization-method`, and not the four optimizer fields added
in the working tree this session. That is the structural reason the class
survived being declared closed. Extending these three tests to every schema key
and every subcommand flag would convert the whole class from an audit finding
into a build failure.

## 9. The default method's headline score is never validated

The project has more external validation than its documentation suggests. A
published-dataset benchmark suite exists: seven datasets from five papers, 101
tests, all passing in 8 seconds. It records negative results honestly, including
an occupancy selectivity model that does not predict on-target yield, with a
test that fails if a future change makes it appear to work.

The gap is specific. For the default `hybrid` method the `score` written to
`step4_improved_df_summary.json` is a predicted fold-amplification
(`hybrid_optimizer.py:1625`). Its shape constants, five-fold per site and a
crossover at component size ten (`network_optimizer.py:92-93`), appear nowhere
in `docs/SCIENCE_CITATIONS.md`. Meanwhile 310 qPCR fold measurements sit in
`tests/validation/data/dwivedi_yu_2023_plasmid_amplification.json`. Three test
files read them, and all three test sequence-level properties or record what the
dataset cannot answer. None compares the shipped amplification model to the
measured amplification.

Other items on the same axis:

- **No external oracle for melting temperature.** The nearest-neighbour tables
  are cited and pinned against transcribed published values, but no test
  compares an assembled Tm to a reference implementation.
- **Additive behaviour at the operating temperatures is unmeasured.** The
  coefficients at 37 C are individually sourced. The activation energies that
  extrapolate them to 30 C and 42 C are all annotated as estimated, two of them
  noting that no multi-temperature data was available
  (`mechanistic_params.py:52-138`).
- **The 3 kb realistic reach was fitted once**, to one breadth proxy on one
  organism, with a band of 3.0 to 6.2 kb. Its test skips on a clean clone
  because `tests/validation/genomes/` is untracked.
- **The equiphi29 and bst reach figures are reasoned in comments but not
  sourced.** The comments explain the inheritance and the processivity bound;
  neither points to a measurement.

## 10. Unreachable code

| Module | Lines | Package callers | Test files |
|---|---:|---:|---:|
| `core/swga_simulator.py` | 723 | 0 | 7 |
| `core/gpu_acceleration.py` | 585 | 0 dispatching | - |
| `core/minimal_primer_selector.py` | 492 | 0 | 1 |

`swga_simulator` is kept alive only by its own tests; the `simulate` command
uses `replication_simulator` instead. A fourth item is a single function rather
than a module: the module-level `optimize()` in
`neoswga/core/dominating_set_optimizer.py:1153` writes its own
`step4_improved_df.csv` at `:1241-1252` with six columns, `set_index` hardcoded
to zero and none of the metrics the real writer emits. Nothing calls it; every
importer takes the class instead. It is dead rather than divergent, but it would
produce a second incompatible schema if it were ever wired up. `minimal_primer_selector` is listed in
CLAUDE.md as an optimizer post-process, but `unified_optimizer.py:769` records
that it has no caller and that its coverage semantics are wrong: it counts
binding coordinates, so thirty sites on a thirty-kilobase genome read as 0.1
percent coverage and no target is reachable.

## 11. Quality gates and test quality

| Measure | Value |
|---|---:|
| ruff findings across `neoswga` and `tests` | 2,793 |
| of those, auto-fixable | 2,175 |
| mypy errors | 182 in 50 files |
| in-body `pytest.skip` calls | 92 |
| tests skipped in a local run | 75 |
| functions with "coverage" in the name | 56 |

ruff, mypy and pip-audit all run with `continue-on-error: true`, the coverage
gate is 50 percent, and `black --check` covers `neoswga/` only. The nine
competing coverage implementations reported in TECH_DEBT_2026-09 are unresolved.

Five tests assert against source text and would pass against a
constant-returning function, at `tests/test_hybrid_optimizer_run.py:218`,
`tests/test_byo_oligos.py:252`, `tests/test_design_options_have_effect.py:265`,
`tests/test_host_profile.py:321` and
`tests/test_dimer_screen_runs_once_per_pool.py:242`. Three assertions cannot
fail at all, at `tests/test_delivered_pool_is_screened_for_dimers.py:57`,
`tests/test_position_cache.py:370` and `tests/test_hybrid_optimizer.py:312`.

One suite race is real but undiagnosed: `tests/test_schema_renderer.py` passed
on one full run and failed on another under eight workers.

## 12. Smaller observations

- **The reference documentation is honestly labelled but thin.** Checking the
  module names in `docs/reference/MODULE_REFERENCE.md`,
  `docs/development/architecture.md` and `docs/development/optimizers.md`
  against the package turns up 13, 7 and 10 names that do not exist. Every one
  of them appears inside a dated accuracy warning that names it as removed, so
  this is disclosure rather than a false claim. The residual gap is coverage,
  which those warnings state themselves: 27 of 92 core modules have no entry in
  the module reference and 86 of 92 have none in the architecture document.
  `docs/development/optimizers.md` is marked superseded in full and still
  ships, describing a 17-optimizer architecture against the six the CLI accepts.
- **The documented output schemas do not match what is written.** Every file
  name in CLAUDE.md's data-flow section is correct; the columns and JSON keys
  are not. `step2_df.csv` is documented at CLAUDE.md:66 as holding "fg_freq,
  bg_freq, gini, Tm"; a file written by the current code holds `primer,
  fg_count, bg_count, gini, gini_bool, ratio, occupancy_ratio`, so only `gini`
  matches and no melting temperature is written at all. In
  `step4_improved_df_summary.json` the key documented as `coverage` is actually
  `fg_coverage`, and the ten metrics listed at CLAUDE.md:87-88 as though they
  were top level are all nested under `metrics`
  (`neoswga/core/base_optimizer.py:260`, `:536-555`). A script reading either
  file by the documented names fails.
- **`filter_stats.json` is written on a best-effort basis.** The write is
  wrapped and a failure is logged at debug level
  (`neoswga/core/pipeline.py:1333-1336`), so at default verbosity a missing
  funnel file looks like a file that was never meant to be there.
- `plan-pool` is a substantial new command that sweeps panel size against
  coverage and selectivity targets with no size cap, which is the gap CLAUDE.md
  describes for `--auto-size` and `--show-frontier`, both limited to 20 primers.
  It is documented only in the uncommitted `docs/guides/optimization_guide.md`,
  not in CLAUDE.md, README or the `neoswga-cli` skill, and appears in no plan or
  audit. Its planner logic is well covered but against a fake optimizer with
  hardcoded metrics; `run_plan_pool` itself is tested only for parser
  registration.
- On the Wolbachia run in the working tree the foreground frequency gate removed
  nothing, `after_fg_frequency` equalling `total_kmers` at 874,596, because
  `min_fg_freq` of 1e-6 against a 1.27 Mb target admits any k-mer present at
  all. The Gini stage was the largest single cut, 491,836 to 20,670, where
  CLAUDE.md states `after_max_primer_cut` usually is.
- That same configuration requests six primers for a 1.27 Mb target. At the
  three-kilobase selection reach six primers cannot approach the coverage such a
  design usually aims for.
- `examples/wolbachia_pool_design/` is safe to commit: its own `.gitignore`
  excludes `input/` and `work/`, which hold all 684 MB.
  `docs/validation/refinement_comparison_2026-09/` is 280 KB of measurements and
  should be committed. The untracked provenance sidecars in
  `examples/plasmid_example/` need a keep-or-ignore decision.
- `validate --quick` and `validate --smoke` both pass, the latter in 4.8 s. The
  smoke target is small enough that one primer covers it completely, so it
  exercises the optimizer only shallowly.

## Suggested order

1. Decide section 2, whether `--num-primers` is a guarantee or a request. It
   governs whether the four failing tests are stale or the default is wrong, and
   nothing else can land until it is settled.
2. Write the delivered-panel dimer test from section 3 first, so the decision is
   recorded as an assertion rather than as a default.
3. `git add` the four new modules and three new tests, or move them aside. Do
   not commit the tracked edits alone.
4. Bring the ratchets back under budget, or record a reviewed budget increase.
5. Remove the unknown `verbose` key from the Wolbachia example config.
6. Push the branch so CI sees it.
7. Correct the false closing line in the pipeline four-steps audit, and the
   claim in CLAUDE.md that the inert-option class is closed.
8. Fix the three argparse defaults that beat params.json, and
   `expand-primers --optimization-method`, using the `None` sentinel already
   established elsewhere. These are small and independent of the decision above.
9. Decide, for each of the six inert keys, whether to wire or delete. Deleting
   is the cheaper honest answer for the four that appear only in an unread
   defaults dictionary.
10. Re-derive the GC-tier designs and record where they live, or move the claims
    that depend on them into documents carrying their own measurements.
11. Correlate the shipped amplification model against the 310 qPCR measurements
    already in the repository, or stop reporting its output as the headline
    score.
