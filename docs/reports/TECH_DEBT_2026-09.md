# Technical debt assessment, September 2026

Branch `audit-alternatives-and-scaling`, commit `86f8b8f`. Every number below was
measured on 2026-09-02 against that tree; the commands are listed in
[Method](#method) so each one can be re-run.

This report covers *structural* debt: what makes the code slow to change and
prone to a recurring class of defect. It is deliberately complementary to
[AUDIT_2026-08_alternatives_and_scaling.md](AUDIT_2026-08_alternatives_and_scaling.md),
which covers behavioural findings (F0-F12). Where a structural cause explains one
of those findings, the link is drawn rather than the finding restated.

---

## Summary

The codebase is in better shape than its size suggests. Duplication is 0.5%,
comment and docstring density is 36%, there are no circular imports, the test
suite is green in 3.5 minutes, and two size ratchets already prevent the
god-file from regrowing. Most of the classic debt categories are simply absent
here, and [What is not debt](#what-is-not-debt-here) says so explicitly so that
effort is not spent there.

The debt that does exist is concentrated, and almost all of it traces to one
design decision.

**`neoswga/core/parameter.py` is a mutable module-global singleton.** A single
539-line function, `get_params`, assigns 62 module globals with a cyclomatic
complexity of 103, and downstream code reads those globals through
`getattr(parameter, name, default)`. When a key is added to `params.schema.json`
and documented but no matching global is assigned, the read site silently
returns its default. Nothing warns, because the schema sets
`additionalProperties: true`.

That is not a hypothetical. Nine instances have shipped. Five have since been
fixed, one was already documented as open, and this audit found three more that
were not.

| Debt category | Verdict | Headline measure |
|---|---|---|
| Configuration architecture | **Severe** | 1 function, 62 globals, complexity 103; 9 inert-key bugs shipped, 4 open |
| Coverage rule consistency | **Severe** | 9 implementations of one quantity; 1 stale copy ships a wrong number |
| Complexity concentration | **High** | 7 F-grade blocks; the top 4 hold 21% of all package complexity |
| Test verification gaps | **High** | 10 tests skip themselves when the feature under test produces nothing |
| Dead code | Moderate | 5 modules unreachable (2,407 LOC); 26 of 29 exception classes unused |
| Startup performance | Moderate | 0.93 s of a 1.37 s `--help` is an unnecessary sklearn import |
| Quality gate enforcement | Moderate | ruff, mypy and pip-audit all run `continue-on-error: true` |
| Duplicated business rules | Moderate | 3 GC classifiers, 2 threshold sets, disagreeing today |
| Code duplication | **Low** | 349 LOC exact structural duplication (0.5%) |
| Documentation | **Low** | 36.1% comment and docstring density; 58 documents |
| Dependency health | **Low** | 1 advisory in the declared dependency closure |

---

## What is not debt here

Stating the negatives matters as much as the findings, because each one is a
place where remediation effort would be wasted.

**Copy-paste duplication is low.** An AST-normalised scan of every function over
12 lines found 13 clusters totalling 349 redundant lines, 0.5% of the package. A
separate 25-line token-window scan found only two cross-file near-duplicate
regions, one of which is legitimate optimizer-factory registration boilerplate.
There is no copy-paste problem to solve, and a de-duplication pass aimed at
repeated text would find almost nothing.

That is a narrower statement than "duplication is low", and the distinction is
the point. The expensive duplication here is one rule reimplemented several times
with drifting semantics, which no text-similarity tool detects because the copies
do not resemble each other. Findings 2 and 6 are both of that kind, and Finding 2
is the most consequential item in this report.

**There are no circular imports.** Zero cycles in the module-level import graph.

**Documentation is not missing.** 36.1% of all lines are comments or docstrings,
across 58 documents in `docs/`. The docstrings explain *why*, not just *what*:
`_configured_int` in `neoswga/core/filter.py:395` documents the MagicMock
interaction that motivated its type check rather than a truthiness check. The
documentation problem in this codebase is accuracy in a few places, not volume.

**The test suite is fast and honest.** 3,789 tests pass in 208 seconds with no
failures, no unconditional skips, no `time.sleep`, and no hardcoded absolute
paths. Every one of the 85 core modules is referenced somewhere in `tests/`.

**Self-labelled debt is zero.** Not one `TODO`, `FIXME`, `HACK` or `XXX` marker
in the package. The debt here is structural, not deferred-and-annotated.

**Dependency health is good.** Of the ten declared runtime dependencies, one
carries an advisory: biopython 1.86, PYSEC-2026-1221, fixed in 1.87. Every other
vulnerability reported by `pip-audit` belongs to the surrounding conda
environment, not to the package's dependency closure. Two optional extras
(`deep-learning-torch`, `deep-learning-tf`) declare packages that are imported
nowhere, and `openpyxl` and `seaborn` in the `viz` extra are likewise unimported.

---

## Finding 1: the parameter module is a mutable global singleton

**Severity: severe. This is the root cause of a recurring defect class.**

`neoswga/core/parameter.py:857` defines `get_params`, which spans 539 lines and
has a cyclomatic complexity of 103. It declares 62 `global` statements and
assigns 62 module-level names. Downstream modules read configuration by
attribute lookup on the module object, with a fallback:

```python
# neoswga/core/filter.py:408
value = getattr(parameter, name, None)
```

A key that reaches `get_params`'s returned dictionary but is never assigned to a
module global therefore resolves to its default, silently. `params.schema.json`
sets `additionalProperties: true`, so a misspelled or unrouted key produces no
warning either.

### Three new instances, demonstrated

A dynamic probe set every one of the 79 schema keys to a distinctive, type-valid
value, ran `get_params`, and read back what the code actually sees. Seventy-two
keys route correctly. Three that do not were previously unrecorded:

| params.json key | value set | value the filter used | read site |
|---|---|---|---|
| `max_homopolymer_run` | 9 | 5 (default) | `core/filter.py:422` |
| `gc_clamp_window` | 8 | 5 (default) | `core/filter.py:438` |
| `max_gc_in_clamp` | 4 | 3 (default) | `core/filter.py:439` |

None of the three has a module global on `parameter` at all. All three are
declared in `params.schema.json`, accepted without complaint, and inert. They
control primer rejection rules, so a user tuning them gets an unchanged candidate
pool and no indication why.

The fourth unrouted key, `optimization_method`, is already recorded as Known
Issue #8 and audit finding F1b. `schema_version`, `adaptive_gc` and
`bg_seq_lengths` also lack globals but are consumed correctly from the returned
dictionary, so they are not defects.

### The class, not the instances

Nine instances of this defect have shipped. The distinction between fixed and
open matters, so both lists are given in full:

| Key | Status |
|---|---|
| `coverage_reach` | Fixed, assigned in `_apply_params_only_keys` |
| `occupancy_ranking` | Fixed, same helper |
| `occupancy_shortlist` | Fixed, same helper |
| `max_mismatches` | Fixed, same helper |
| `sampled_index_path` | Fixed, same helper |
| `optimization_method` | **Open** (Known Issue #8, finding F1b) |
| `max_homopolymer_run` | **Open**, found here |
| `gc_clamp_window` | **Open**, found here |
| `max_gc_in_clamp` | **Open**, found here |

The three found here are worse than `optimization_method` in one respect: they
have a live read site that runs on every filter invocation, and the docstrings
around it assert the opposite of the truth. `core/filter.py:415` says the value
"Was a module constant... Read per call now" and `core/filter.py:434` says
"Configurable so that conclusion can be re-tested." Neither is reachable from
params.json. The only thing that has ever set them is
`tests/test_sequence_rule_configuration.py`, which monkeypatches the module
attribute directly and therefore never exercises the params.json route at all.

Git history shows at least eight separate commits fixing members of this family,
among them `8982254` ("stop the inert CLI flags reporting success"), `bd67d37`
("four bugs found sweeping the pipeline for defaults that never fire") and
`946250d` ("warn when `--use-mechanistic-model` flag is silently ignored").

The same shape reaches the CLI. `neoswga analyze-set --simulate` is declared at
`cli/analysis.py:508` with the help text "Run replication simulation". The string
`simulate` appears nowhere else in that file: the handler never reads it, and no
simulation runs. Unlike `--fg` and `--fg-kmers` on the same subcommand, which are
honestly labelled "(unused)" and warned about at runtime, this one makes a
promise it does not keep.

Four regression tests guard the class today:
`test_design_options_have_effect.py`, `test_params_json_routes_optional_keys.py`,
`test_params_json_routes_reach_the_pipeline.py` and
`test_unimplemented_cli_options_are_honest.py`. They are hand-maintained
allowlists covering roughly 30 named keys out of a surface of 79 schema keys and
372 CLI flags. They catch the instances someone thought to list. They cannot
catch the tenth.

### Why it keeps recurring

Adding a parameter correctly requires four coordinated edits in three files: a
schema property, a `global` declaration, an assignment inside a 539-line
function, and a read site. Three of the four are silent when omitted. The
function's complexity of 103 means no reviewer holds all its branches in mind,
and its length is already pinned by a ratchet whose own comment concedes the
point:

> "Adding a parameter to the function whose job is assigning parameters is the
> fix; the real remedy for the size is splitting this function, not this line."
> (`tests/test_function_length_ratchet.py:33`)

---

## Finding 2: nine implementations of coverage, and one of them ships a wrong number

**Severity: severe. This one changes a number a user acts on.**

"Fraction of the genome within reach of a binding site" is implemented nine times
with four incompatible definitions:

| Implementation | Method |
|---|---|
| `core/base_optimizer.py:121` `_union_coverage` | Bitmask, symmetric window, circular-aware. Canonical. |
| `core/coverage.py:19` `compute_per_prefix_coverage` | Shares `_mark_window` with the above |
| `core/hybrid_optimizer.py:973` `_calculate_coverage` | Binned, via `BipartiteGraph` |
| `core/background_aware_optimizer.py:375` `_calculate_coverage` | Interval merge, strand-directional |
| `core/set_size_optimizer.py:648` `_compute_coverage` | Sorted-run merge, downstream-only processivity |
| `core/dominating_set_optimizer.py:296` `_genome_fraction` | Binned |
| `core/primer_expansion.py:332` `_calculate_coverage` | Binned |
| `core/efficiency_predictor.py:467` `_calculate_coverage` | Binned |
| `core/results_interpreter.py:331` `_calculate_coverage` | Binned |

The codebase already knows this is a hazard. The canonical implementation says so
in its own docstring:

> Three implementations of one quantity is how this codebase has produced
> disagreeing coverage numbers before; one audit found three different semantics
> for "coverage" at once.
> (`core/base_optimizer.py:126`)

There are nine now.

### The stale copy

`BipartiteGraph.add_primer_coverage` marks a bin covered when any part of it is
within reach of a binding site. That approximation only holds while a bin is no
larger than the reach, which is why `coverage_bin_size` exists at
`core/dominating_set_optimizer.py:62`. Its docstring records what happened at the
shipped defaults of 10 kb bins against a realistic 3 kb reach: greedy set cover
saw the genome covered after two primers and stopped, and "the coverage reported
for that set was 1.000 where the honest figure was 0.433."

Seven of the eight `add_primer_coverage` call sites now pass `extension_reach`
and route their bin size through `coverage_bin_size`. One does not:

```python
# neoswga/core/efficiency_predictor.py:477
bin_size = 10000
graph = BipartiteGraph(bin_size=bin_size)
...
graph.add_primer_coverage(primer, positions, prefix, length)   # no extension_reach
```

`extension_reach` defaults to 0, so every binding site claims a full 10 kb bin.
This is exactly the configuration the fix was written to eliminate, in the one
place the fix never reached. It is live code: `predict-efficiency` at
`cli/commands.py:457` builds an `EfficiencyPredictor`, whose `predict` calls
`_calculate_coverage` at line 249.

### What it costs

Measured on the *Wolbachia* wMel genome (1,267,782 bp) with 40 binding sites per
primer and a 3 kb reach, comparing the three figures:

| Primers | `predict-efficiency` reports | Every other call site | Base-by-base truth |
|---|---|---|---|
| 2 | 0.449 | 0.350 | 0.176 |
| 4 | 0.685 | 0.544 | 0.309 |
| 8 | 0.929 | 0.816 | 0.529 |
| 12 | **0.992** | 0.922 | **0.673** |

The error grows with set size and converges on 1.000, which is the direction a
reader will not question. `predict-efficiency` is the pre-synthesis go/no-go
check, described in its own help text as "Predict efficiency of primer set before
synthesis". It currently reports near-total coverage for a set that covers two
thirds of the target.

The fix is one line: pass `extension_reach` and route `bin_size` through
`coverage_bin_size`, as the other seven call sites already do.

---

## Finding 3: complexity is concentrated in four files

Average cyclomatic complexity across 1,678 blocks is 4.80, a grade A. The
distribution is what matters: 1,237 blocks are A, and 12 are E or F.

| Block | File | Complexity |
|---|---|---|
| `get_params` | `core/parameter.py:857` | 103 |
| `run_step4` | `cli/pipeline.py:550` | 68 |
| `run_optimization` | `core/unified_optimizer.py:535` | 62 |
| `run_step2` | `cli/pipeline.py:150` | 49 |
| `render_technical_report` | `core/report/technical_report.py:901` | 48 |
| `HybridOptimizer.optimize` | `core/hybrid_optimizer.py:376` | 44 |
| `ThermodynamicFilter.filter_candidates` | `core/thermodynamic_filter.py:259` | 42 |

Weighting each file by how often it changed in the last twelve months gives the
risk ranking. High churn against high complexity is where defects concentrate.

| File | Commits (12 mo) | Max block complexity | Risk score | LOC |
|---|---|---|---|---|
| `core/parameter.py` | 26 | 103 | 2678 | 1411 |
| `core/unified_optimizer.py` | 30 | 62 | 1860 | 1374 |
| `cli/pipeline.py` | 18 | 68 | 1224 | 1618 |
| `core/hybrid_optimizer.py` | 22 | 44 | 968 | 1361 |
| `cli_unified.py` | 57 | 13 | 741 | 431 |
| `core/pipeline.py` | 22 | 29 | 638 | 1152 |

`cli/pipeline.py` is the only file in the package with a maintainability index
below grade A, at B (15.23). It also holds 11 of the 134 blind `except Exception`
handlers and 41 of the 500 function-local imports.

Three of the top four risk files have no dedicated test module:
`unified_optimizer.py`, `base_optimizer.py` and `background_aware_optimizer.py`
are exercised only incidentally through other test files.

---

## Finding 4: ten tests skip themselves when the feature fails

A test that skips when its subject produces no output reports green on the exact
failure it was written to detect. There are ten such sites:

| Site | Skips when |
|---|---|
| `tests/cli/conftest.py:225` | `step3_df.csv` was not produced |
| `tests/cli/test_diagnostic_commands.py:108` | `validate-model --output-json` printed nothing |
| `tests/cli/test_diagnostic_commands.py:138` | `simulate` produced no output directory |
| `tests/cli/test_diagnostic_commands.py:169, 204` | `simulate` wrote no JSON |
| `tests/test_12bp_pipeline.py:98` | the `score` column is absent |
| `tests/test_12bp_pipeline.py:105` | the `selectivity` column is absent |
| `tests/test_string_search_positions.py:55, 140, 152` | the fallback path returns `None` |

The failure mode is concrete. `tests/cli/conftest.py:225` is a session-scoped
fixture: if the `score` step stops writing `step3_df.csv`, every test depending
on `scored_primers` skips rather than fails. Audit finding F0 records that the
`score` step currently changes nothing downstream. This is the shape of test that
would let that go unnoticed.

The remaining 88 skip markers are legitimate: absent optional dependencies
(Jellyfish, pyahocorasick, plotly, jsonschema, PyYAML) or absent example data.
There are no unconditional skips and no expected-failure markers hiding
regressions.

---

## Finding 5: dead code, and an abstraction that was never adopted

Transitive reachability from the entry point `neoswga.cli_unified`, walking every
`import` node including the 500 function-local ones, reaches 109 of 114 modules.

| Module | LOC | Test references | Note |
|---|---|---|---|
| `core/swga_simulator.py` | 723 | 13 | Largest orphan |
| `core/simulation_analysis.py` | 505 | 6 | Only comment mentions in production |
| `core/minimal_primer_selector.py` | 492 | 4 | **CLAUDE.md:26 says this runs as a post-process** |
| `core/search_context.py` | 399 | 1 | |
| `logging_config.py` | 288 | 2 | See below |

Two of these deserve attention beyond their line count.

`CLAUDE.md:26` states that `minimal_primer_selector` runs "as a post-process" in
the optimizer dispatch. It has zero production references. The documentation
asserts a wiring that does not exist.

`logging_config.py` defines seven public entities including `setup_logging`,
`get_logger`, `ProgressLogger` and a colour formatter. `get_logger` has zero
callers in the package. Meanwhile the package calls `logging.getLogger` directly
at 93 sites and `print()` at 892 sites, 122 of them in
`core/workflow_selector.py` and 111 in `core/wizard.py`. A logging abstraction
was written and never adopted; library modules print to stdout instead, which
makes their output impossible for a caller to redirect or suppress.

### An exception hierarchy nothing raises

`neoswga/core/exceptions.py` is 413 lines defining 29 exception classes. Three
are referenced anywhere outside the file: `NoCandidatesError`,
`OptimizerNotFoundError` and `InvalidParameterError`. The other 26, including
`OptimizerConvergenceError`, `InsufficientCoverageError`, `PipelineStateError`
and `MemoryLimitError`, are referenced by no production code and by no test.

This sits alongside the 134 `except Exception` handlers, 29 of which log only at
DEBUG level and 8 of which pass silently. A precise exception hierarchy was
designed and then not used, while the code that would have used it catches
everything instead.

### Totals

| Category | Count | LOC |
|---|---|---|
| Unreachable modules | 5 | 2,407 |
| Public definitions referenced nowhere | 41 | 1,002 |
| Referenced only by documentation | 12 | 451 |
| Referenced only by tests | 36 | 877 |

Roughly 4,700 lines, about 7% of the package. Around 1,700 lines are safe to
delete outright: `search_context.py` and `logging_config.py` (687), the 26 unused
exception classes (roughly 234), `expand_primers` and `ExpansionInput` (150), and
the three benchmark helpers `benchmark_network_vs_ratio`, `benchmark_gpu_vs_cpu`
and `log_gpu_status` (90).

Five names look dead to a text search and are not: `CliqueOptimizer`,
`HybridBaseOptimizer`, `NetworkBaseOptimizer`, `BackgroundAwareBaseOptimizer` and
`DominatingSetAdapter` all reach production through `@OptimizerFactory.register`
decorators. Building the registry at runtime lists all five. Any deletion pass
must resolve decorator registration before trusting a reference count.

The rest is not safe to delete unexamined. Each is a feature that was built and
tested but never connected, so the choice for each is to wire it up or retire it,
and that is a product decision rather than a cleanup. `swga_simulator.py` and
`simulation_analysis.py` are 1,228 lines with four test files and a user-facing
guide at `docs/simulator-guide.md`; that is an unfinished feature, not debris.

---

## Finding 6: one duplicated algorithm and one duplicated scientific rule

Duplication by volume is negligible: 0.5% by exact structural match, 1.5% by a
looser reading that counts divergent reimplementations of one rule. The costly
duplication in this codebase is not copied text. The largest instance is the
coverage rule in Finding 2. Two more are real.

**`_prune_background`, 74 lines, duplicated verbatim.** Present at
`core/hybrid_optimizer.py:1124` and `core/background_aware_optimizer.py:289`. The
two bodies differ in 57 diff lines, almost all comments and docstring wording,
plus one substantive divergence: the hybrid copy reads
`self.min_coverage_threshold` and the background-aware copy reads
`self.min_coverage`. The same concept already carries two names. A correction to
the pruning rule must be made twice, and the copies have begun to drift.

**GC classification is implemented three times, with two threshold sets, and they
disagree today.**

| GC | `gc_adaptive_strategy._classify_genome` | `condition_suggester.classify_gc` | `genome_analysis.get_gc_class` |
|---|---|---|---|
| 0.30 | AT_RICH | at_rich | at_rich |
| **0.3523** | **BALANCED** | **at_rich** | **at_rich** |
| 0.38 | **BALANCED** | **at_rich** | **at_rich** |
| 0.50 | BALANCED | balanced | balanced |
| 0.62 | **BALANCED** | **gc_rich** | **gc_rich** |
| 0.66 | GC_RICH | gc_rich | gc_rich |

`gc_adaptive_strategy.py:146` uses boundaries of 0.35 and 0.65; the other two use
0.40 and 0.60. The classifiers disagree across two bands, 0.35 to 0.40 and 0.60 to
0.65.

This is live, not theoretical. *Wolbachia* wMel, the project's current worked
example, is 35.23% GC. It falls in the first disagreement band: the GC-adaptive
strategy classifies it as balanced while the other two classify it as AT-rich.
Whether that changes the chosen conditions was not traced here, but a genome
being AT-rich and not AT-rich at the same time within one run is a defect
regardless of its downstream effect.

Related: 451 magic-value comparisons remain across the package, concentrated in
`reaction_conditions.py` (44), `genome_analysis.py` (32) and
`condition_suggester.py` (24). For a tool whose thresholds are scientific
calibrations, an unnamed literal is a calibration nobody can find or cite.

---

## Finding 7: every command pays 0.93 s for an import it does not need

`neoswga --help` takes 1.37 s. Import profiling attributes 0.93 s of that, 68%,
to a single chain:

```
neoswga/cli_unified.py:60   from neoswga.core.pipeline import StepPrerequisiteError
neoswga/core/pipeline.py:16 from neoswga.core import parameter, rf_preprocessing, ...
neoswga/core/rf_preprocessing.py:61-62  import sklearn / import sklearn.ensemble
                                        -> sklearn.utils -> scipy.stats
```

`cli_unified.py` imports the whole pipeline module at module scope to obtain
`StepPrerequisiteError`, a four-line exception class. That drags in
`rf_preprocessing`, which imports scikit-learn eagerly, which imports
`scipy.stats`. All 38 subcommands pay it, including `--help`, `show-presets` and
`schema`, none of which score anything.

The fix is consistent with what the codebase already does elsewhere: 500
function-local imports exist across the package, 245 of them internal. The one
place where deferring the import would pay for itself is the one place it was not
done.

---

## Finding 8: the quality gates do not gate

`.github/workflows/ci.yml` runs six test matrix cells across three Python
versions and two operating systems, installs Jellyfish, and smoke-imports the
built wheel in a clean virtualenv. That build job is genuinely well constructed
and guards a real packaging-omission class.

The static-analysis half does not gate anything:

| Check | CI status | Current count |
|---|---|---|
| `black --check` | blocking | clean |
| `isort --check` | blocking | clean |
| `ruff check` | `continue-on-error: true` | 11,334 findings |
| `mypy` | `continue-on-error: true` | 172 errors in 50 files |
| `pip-audit` | `continue-on-error: true` | 1 in-scope advisory |
| coverage | blocking at 50% | actual figure not tracked anywhere |

Three checks run, report, and are ignored. Nothing prevents the counts from
growing, and no trend is recorded, so nobody can tell whether they are improving.
The coverage floor of 50% is far below the level that would have caught the
findings above, and the actual percentage is not written to any file in the repo,
so there is no baseline to ratchet against.

The two ratchets that *do* exist, `test_module_size_ratchet.py` and
`test_function_length_ratchet.py`, are the right pattern. Both pin current values
as ceilings, both refuse stale entries, and both carry comments explaining why
each pin exists. They are the model the other gates should follow.

Selected ruff counts, for context on where the bulk sits:

| Rule | Count | Reading |
|---|---|---|
| `UP006` non-PEP-585 annotation | 1136 | Mechanical, autofixable |
| `G004` f-string in logging call | 1065 | Mechanical |
| `T201` `print()` | 892 | Real: library code writing to stdout |
| `PLC0415` import outside top level | 500 | Real: hidden coupling |
| `PLR2004` magic value comparison | 451 | Real: uncited calibrations |
| `F401` unused import | 174 | Mechanical, autofixable |
| `BLE001` blind except | 132 | Real: see below |
| `PLW0603` global statement | 116 | Real: 75 of them in `parameter.py` |

The 134 `except Exception` handlers include 29 that log only at DEBUG level and 8
that pass silently. There are no bare `except:` clauses, and the package defines
34 custom exception types, so the machinery for precise handling exists and is
partly unused.

---

## Remediation roadmap

Effort figures are estimates for a developer familiar with this codebase. The
measured facts are the counts and timings above; the hours are judgement.

### Quick wins: high value, under a day each

**1. Pass `extension_reach` in `efficiency_predictor._calculate_coverage`.**
One line, plus routing `bin_size` through `coverage_bin_size`, exactly as the
other seven `add_primer_coverage` call sites already do. Stops
`predict-efficiency` reporting 0.992 coverage where the measured figure is 0.673.
This is a wrong number a user acts on before spending money on synthesis, so it
outranks everything else on speed of payback. Estimated 1 hour, plus a test that
pins the three coverage figures against each other.

**2. Move the sklearn import behind a function boundary.** Import
`StepPrerequisiteError` from a leaf module, or defer `rf_preprocessing`'s
scikit-learn import into the functions that use it. Removes 0.93 s from every
invocation of all 38 subcommands. Estimated 1 to 2 hours. Verify by re-running
the import profile in [Method](#method).

**3. Turn on parallel test execution.** `pytest -n 8` runs the same 3,789 tests
in 103 s against 208 s serial, with identical results, so the suite is already
parallel-safe. `pytest-xdist` 3.8.0 is installed locally but is absent from the
`dev` extra and from CI. Add it to both. Halves the feedback loop on six matrix
cells. Estimated 1 hour.

**4. Fix the three inert filter keys.** Assign `max_homopolymer_run`,
`gc_clamp_window` and `max_gc_in_clamp` in `_apply_params_only_keys`, which
already exists at `core/parameter.py:830` for exactly this purpose and already
handles five keys the same way. Extend
`tests/test_params_json_routes_optional_keys.py` to cover them. Estimated 2 hours.

**5. Upgrade biopython to 1.87.** Clears PYSEC-2026-1221, the only advisory in
the declared dependency closure. Estimated 1 hour including a test run.

**6. Delete the four unimported optional extras.** `deep-learning-torch`,
`deep-learning-tf`, and `openpyxl` and `seaborn` from `viz`, none of which is
imported anywhere. Estimated 30 minutes.

**7. Autofix the mechanical lint.** `ruff check --fix` clears `UP006`, `F401`,
`F541`, `COM812` and `RET505`, roughly 2,400 findings, with no behaviour change.
Do it as one commit so it never obscures a real diff. Estimated 2 hours including
review of the diff.

### Structural work: the ordered sequence

**8. Replace the hand-maintained routing allowlists with a generated test.**
Do this *before* the parameter refactor, because it is the safety net for it. The
dynamic probe used in this audit is roughly 60 lines: build a params.json setting
every schema key to a distinctive value, run `get_params`, assert each key either
reaches a module global or appears on a small explicit exemption list with a
stated reason. That test would have caught all nine instances of the bug class
including the three found here, and it catches the tenth without anyone
remembering to add it. Estimated 4 to 6 hours. **This is the single highest-value
item in the report.**

**9. Set `additionalProperties: false` in `params.schema.json`.** A misspelled key
becomes a validation error rather than silence. This is a breaking change for any
params.json carrying extra keys, so it wants a deprecation window: warn on unknown
keys for one release, then reject. Estimated 3 hours plus the release cycle.

**10. Split `get_params`.** With the generated test in place this becomes safe.
The natural seams are already visible in the function: filtering parameters,
thermodynamic parameters, reaction conditions, optimizer parameters, and
GC-adaptive resolution. Extract each into a function taking the parsed dictionary
and returning a typed dataclass, keeping the module globals as a thin
compatibility shim assigned from those dataclasses so no read site has to change
in the same commit. Then lower the ratchet entry from 539 as each piece lands.
Estimated 3 to 5 days. Retires the defect class at its root rather than
instance by instance.

**11. Unify the three GC classifiers.** One function, one threshold set, one
citation for where the boundaries come from. Decide deliberately whether the
boundary is 0.35 or 0.40 rather than inheriting both. Check what changes for
*Wolbachia* at 35.23% GC before and after. Estimated 4 hours plus the scientific
decision.

**12. Extract `_prune_background` to `BaseOptimizer`.** Reconcile
`min_coverage_threshold` and `min_coverage` to one name. 74 lines removed, one
algorithm with one home. Estimated 3 hours.

**13. Convert the ten self-skipping tests into assertions.** Where the artefact is
genuinely optional, assert on the condition explicitly. Where it is not, as with
`step3_df.csv`, fail. Estimated 4 hours.

**14. Decide the fate of the five orphan modules.** For each of the 2,407 lines,
wire it up or retire it. Start with `minimal_primer_selector`, because CLAUDE.md
currently claims it is wired up, so either the code or the documentation is wrong
today. Estimated 1 day for the decisions, more if any is to be connected.

**15. Consolidate the nine coverage implementations.** With the stale copy fixed
by item 1, the remaining work is to reduce nine implementations of one quantity
to one canonical function plus deliberate, documented approximations of it.
`_union_coverage` at `core/base_optimizer.py:121` is already the canonical one and
its docstring already argues the case. Every binned variant should call
`coverage_bin_size` and state in one place that it is an approximation used for
progress reporting, not the authoritative figure. Estimated 3 to 4 days.

**16. Fix the inert `--simulate` flag on `analyze-set`.** Either wire it to the
simulator, or add it to `UNIMPLEMENTED_OPTIONS` in `cli/_common.py:428` so it
warns like the other six. Estimated 1 hour.

### Longer term

**17. Split `cli/pipeline.py`.** 1,618 lines, the only file below maintainability
grade A, holding two of the seven F-grade blocks. `run_step4` at complexity 68 and
`run_step2` at 49 are the targets. The `add_parsers` half is declarative argparse
and can move wholesale. Estimated 1 week.

**18. Give the top-risk optimizers dedicated test modules.**
`unified_optimizer.py`, `base_optimizer.py` and `background_aware_optimizer.py`
are the three highest-churn modules without one. Estimated 1 week.

**19. Adopt `logging_config` or delete it.** 892 `print()` calls in library code
make output uncontrollable by callers. If the module is to be adopted, the
`workflow_selector` and `wizard` modules hold 233 of those calls between them and
are the place to start. If not, delete 288 lines. Estimated 3 days to adopt.

---

## Prevention

The project already has the right instinct: two ratchets, four routing guard
tests, and CI comments explaining *why* each pin exists. The gap is that the
guards are hand-maintained lists, so they cover the instances someone
remembered, and that three of the six CI checks do not fail the build.

**Generate the guards instead of listing them.** Item 8 above is the pattern:
enumerate the schema, assert the property, keep an explicit exemption list with a
written reason per entry. Apply the same shape to the 372 CLI flags, asserting
that each flag's `dest` is read somewhere outside the parser definition.

**Ratchet the lint counts.** Record today's ruff and mypy totals in a file, add a
test that fails if either grows, and lower the numbers as they improve. This is
exactly what the two size ratchets already do, applied to two more metrics. It
converts 11,334 ruff findings from an ignored wall into a number that can only go
down.

**Record the coverage figure.** CI enforces a 50% floor but writes the actual
percentage nowhere, so there is no trend and no basis for raising the floor. Emit
it to a tracked file, then raise the floor to just below the current value and
repeat.

**Make `mypy` blocking for new files.** 172 errors across 50 files is too many to
fix at once, but a per-file allowlist, same shape as the ratchets, stops file 51
from joining them.

**Give `pip-audit` a scope.** It currently reports the whole conda environment,
which is why one genuine biopython advisory sits among two dozen irrelevant ones.
Auditing the declared dependency closure would make the check actionable, and then
it can be made blocking.

---

## Method

Every figure is reproducible from the repository root on commit `86f8b8f`.

```bash
# Size and complexity
find neoswga -name '*.py' | xargs wc -l | sort -rn | head -25
radon cc neoswga/ -s --total-average
radon cc neoswga/ -n E -s          # E and F grade blocks
radon mi neoswga/ -s | sort -t'(' -k2 -n | head

# Lint and types
ruff check neoswga/ --select ALL --statistics
mypy neoswga --ignore-missing-imports --no-strict-optional

# Risk ranking: churn x complexity
git log --since='12 months ago' --name-only --pretty=format: -- 'neoswga/*.py' \
  | grep -v '^$' | sort | uniq -c | sort -rn | head -20

# Startup cost
python -X importtime -c 'import neoswga.cli_unified' 2>&1 \
  | awk -F'|' 'NR>1 {gsub(/ /,"",$2); print $2, $3}' | sort -rn | head

# Test suite
python -m pytest tests/ -q -o addopts='' --ignore=tests/integration          # 208 s
python -m pytest tests/ -q -o addopts='' --ignore=tests/integration -n 8     # 103 s

# Dependency audit
pip-audit --progress-spinner off
```

The coverage divergence in Finding 2 was measured by building a `BipartiteGraph`
twice over the same synthetic binding positions on the *Wolbachia* genome length,
once with the arguments `efficiency_predictor` passes and once with the arguments
the other seven call sites pass, and comparing both against a base-by-base
NumPy mask. The three columns of that table come from one script run.

Duplication was measured with two purpose-written scripts: an AST-normalised
function hash over every function of 12 lines or more, and a 25-line token-window
shingle scan. Every cluster reported above was read at both sites before being
recorded; clusters that the normalisation matched spuriously were discarded.

Module reachability was computed by walking every `ast.Import` and
`ast.ImportFrom` node anywhere in each tree, not only at module level, since the
package contains 500 function-local imports that a top-level-only walk would miss.

The parameter routing table was produced dynamically rather than statically. A
static check of `get_params` reports 27 unrouted keys and is wrong, because five
keys are assigned in the `_apply_params_only_keys` helper and several others are
consumed from the returned dictionary. The dynamic probe sets each of the 79
schema keys to a distinctive type-valid value, runs `get_params`, and compares
what the code actually reads. That is the table in Finding 1, and it is the test
recommended as item 8.

### Not measured

- Actual coverage percentage. A `--cov` run reached 3% in nine minutes and was
  abandoned as too slow to be worth the wall-clock here; the CI floor of 50% is
  the only enforced figure and the real number is recorded nowhere in the repo.
- Whether the GC classifier disagreement changes the conditions chosen for
  *Wolbachia*. The disagreement itself is demonstrated; its downstream effect is
  not traced.
- Runtime profiling of the pipeline steps. Optimizer timings are already covered
  in the August audit, findings F5 and F5b.
