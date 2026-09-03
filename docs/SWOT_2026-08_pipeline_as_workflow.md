# SWOT: the neoswga pipeline as a workflow

**Date:** 2026-08-31
**Subject:** the four-step pipeline as something a scientist operates end to end
— data flow between steps, provenance, reproducibility, resumability, failure
modes and the surface around them.

**Relationship to the other audit.**
[`AUDIT_2026-08_alternatives_and_scaling.md`](AUDIT_2026-08_alternatives_and_scaling.md)
carries a SWOT of the *methods*: which optimizer to use, how they scale, how far
the returned oligo set is from optimal. This document deliberately does not
repeat it. The question here is different: given that the methods are sound,
what happens to a result as it moves through four commands, three intermediate
files and a config file that can change underneath them?

**Evidence.** Every entry cites `file:line` or a measurement. Claims about
behaviour were checked by running the code, not by reading it, wherever a run
was cheap. Two entries below correct or sharpen findings from the investigation
pass that produced them; those are marked.

---

## The workflow as it actually runs

```
count-kmers          filter              score               optimize
     |                  |                  |                    |
     v                  v                  v                    v
 *_Xmer_all.txt --> step2_df.csv  --> step3_df.csv  --> step4_improved_df.csv
 (cached per k)      positions.h5                          + _summary.json
                     filter_stats.json                     + _validation.json

 every step appends to run_manifest.json  ...  on success only
 all of it lands in one flat data_dir with fixed filenames
```

`design` runs all four in one process, with `--start-from` / `--stop-at`
(`cli/commands.py:567-575`).

---

## Strengths

- **Provenance records the reaction that ran, not the one that was configured.**
  `run_manifest.py:80-173` writes a per-step entry with neoswga version, git
  SHA, jellyfish version, Python version, platform, the raw `sys.argv`, the
  seed, and input checksums. The load-bearing field is `effective_conditions`
  (`cli/_common.py:699-716`), built from the live parameter module *after* CLI
  merge and after the GC-adaptive strategy and `retune_for_polymerase` have had
  their say. Those can change the enzyme and temperature at run time and never
  write back to params.json, so a manifest that recorded only the file would
  misreport the chemistry. Few pipelines make this distinction; making it is
  correct.

- **The reporting layer reads authoritative sources rather than re-deriving.**
  `report/metrics.py:404` prefers `step4_improved_df_summary.json` over
  recomputing from CSV, and `:676-679` reads `effective_conditions` back from
  the manifest, which is what makes a report's Tm agree with the optimizer's.
  Values are badged MEASURED or ESTIMATED rather than presented uniformly.

- **The expensive step is genuinely cached.** `kmer_counter.py:336` guards the
  jellyfish subprocess with `if not os.path.exists(txt_file)`, per k. Re-running
  `count-kmers` reuses existing counts instead of repeating a job that costs
  about seven minutes and 138 MB for hg38 at k=12. This is real resumability,
  not a completion marker.

- **Missing inputs stop the run instead of producing an empty answer.**
  `StepPrerequisiteError` (`pipeline.py:116-144`) is raised by steps 2 and 3
  when the predecessor's outputs are absent (`:845-852`, `:1046-1050`), with a
  remediation message. In a codebase whose recorded failure mode is the
  plausible-looking zero (Known Issues #5 and #7), a step that refuses to run is
  the right instinct.

- **Determinism is tested where it is achievable, and honestly scoped.**
  `tests/integration/test_plasmid_golden_snapshot.py:59` pins metric bounds
  rather than exact sequences — deliberately, because tied optimizers may
  permute — and `:178` asserts two seed=42 runs return identical primer sets.
  Step 2 contains no RNG at all, so it has nothing to seed rather than a seed
  that is ignored.

- **A wide operational surface.** 38 commands, including a setup wizard
  (`init`), pre-flight validation (`validate params`), result interpretation
  (`interpret`), quality reporting (`report`), and an iterative design loop
  (`analyze-coverage`, `expand-primers`, `swap-primer`, `contract-set`) that
  accepts real BAM depth. Most SWGA tooling stops at "here are some primers".

---

## Weaknesses

- **One namespace, no prefix, no lock.** Every step writes fixed filenames into
  `data_dir` (`pipeline.py:995`, `:1138`; `unified_optimizer.py:1041`), there is
  no `--output-prefix` anywhere in the CLI, and no file locking on any pipeline
  output. Two parameter settings cannot coexist in one directory, and two
  concurrent `optimize` runs silently clobber each other — last writer wins,
  with no error on either side.

- **Nothing detects a stale intermediate.** A `params_hash` is computed at
  `unified_optimizer.py:1107` and written into the audit record. It is never
  read: `grep` finds no comparison anywhere in `neoswga/` or `tests/`. There is
  no mtime check either. Edit params.json, re-run only `optimize`, and the run
  silently mixes a step-2 pool built under the old settings with a step-4
  selection made under the new ones. The ingredient for detecting this exists
  and is unused.

- **A config key that pins every params.json user to the slowest method.**
  `optimization_method` is schema-declared and documented, and inert:
  `get_params` assigns no global, and `run_step4` passes
  `args.optimization_method` through (`cli/pipeline.py:766`) with argparse
  default `'hybrid'`. Verified end to end — params.json asking for
  `dominating-set` yields `parameter.optimization_method` unset and hybrid run.
  F5 measured hybrid returning an *identical* set to dominating-set at up to
  260x the cost. Recorded as F1b in the companion audit.

- **`design` cannot reach the parameters that decide the design.** Its
  subparser offers 14 options against `optimize`'s 40, and omits `--seed`,
  `--optimization-method`, `--application`, `--coverage-reach`,
  `--num-primers` and the rest;
  `run_design` then hardcodes `optimization_method="hybrid"` and
  `application="enrichment"` for the missing attributes
  (`cli/commands.py:577-595`). Since `seed` is also absent from
  `params.schema.json`, **the one-shot command cannot be seeded by any route** —
  the path most likely to be used for a whole design is the one that cannot be
  made reproducible.

- **The manifest is a success log, not a run log.** In steps 2 and 3 the write
  sits after the `try/except`, every branch of which ends in `sys.exit(1)`
  (`cli/pipeline.py:348-356`, `:439-444`); in step 4 it is near the end of the
  try block, after the no-results exit. So a `StepPrerequisiteError` — the most
  common resumability failure — leaves no trace at all. The record covers
  exactly the runs least in need of forensics.

- **A corrupt manifest is discarded in silence.** `run_manifest.py:154-161`
  reads the existing file and, on `OSError` or `JSONDecodeError`, falls through
  to `existing = {"steps": []}`. The next write then persists a manifest
  containing one entry, and every prior step's provenance is gone with no
  warning. Combined with the unlocked read-modify-write at `:154-167`, two
  concurrent steps can lose entries. *(Sharpened here: the investigation noted
  the missing lock but not the silent history loss, which is the worse half.)*

- **`resolved_params` is not resolved.** It is a verbatim copy of params.json
  (`run_manifest.py:123-129`); CLI overrides survive only inside the raw argv
  string. The repo knows — `tests/test_run_manifest.py:142-148` documents it —
  but the field name invites exactly the wrong reading, and it sits beside
  `effective_conditions`, which *is* resolved. Two adjacent fields with opposite
  semantics and no naming cue.

- **Step 1 has no prerequisite check.** `validate_step1_prerequisites` is
  defined at `pipeline.py:147-190` and has **zero call sites**. A missing or
  unreadable genome surfaces later as a raw `FileNotFoundError` from
  `kmer_counter.py:418`. The validator that would give the good message exists
  and is dead.

- **134 `except Exception` blocks**, three of them in the step path where a
  swallowed exception lets a step exit 0 with a correctness guard silently
  disabled: exclusion/blacklist k-mer file reads (`pipeline.py:55`, `:103`) —
  an unreadable file scores as "no hits", so a primer passes a screen it should
  fail; blacklist re-injection collection (`unified_optimizer.py:450-466`) —
  returns `[]` on any failure, disabling the check; and the post-optimization
  validator itself (`:895-906`) — if it crashes, duplicates, size drift, zero
  coverage and blacklist re-injection all go unchecked and the run still reports
  success.

- **The manifest omits the dependencies that move the numbers.** No numpy,
  pandas or scikit-learn version; no hash or version of the random forest model
  used in scoring; no checksum of params.json itself, only its path — so if the
  file is edited afterwards, the manifest cannot recover what was used. For a
  tool whose own CLAUDE.md warns that a major sklearn upgrade may warrant
  re-validating the model, the model's identity is the obvious thing to record.

---

## Opportunities

- **`--output-prefix`, or a run-id subdirectory.** The single highest
  ratio of risk removed to work required. It ends the clobbering, lets a
  parameter sweep run in one `data_dir`, and makes concurrent runs safe. The
  benchmarking harness in `scripts/benchmarking/` already works around its
  absence externally, which is evidence both that it is needed and that the
  shape is known.

- **Make `params_hash` load-bearing.** It is already computed. Write it beside
  each intermediate, compare on read, and warn when a step consumes an input
  produced under different parameters. That converts a decorative field into the
  staleness guard the workflow currently lacks entirely.

- **Move the manifest write into a `finally`, and record the exception.** The
  manifest becomes a run log rather than a success log, at the cost of a few
  lines. Failed runs are precisely the ones someone will later need to
  reconstruct.

- **Give `design` a `--seed` and the four flags that change the design**, or
  document it explicitly as a fixed-recipe convenience command. Either is
  defensible; the current state — a full-pipeline command that silently pins
  method and application and cannot be seeded — is not.

- **Narrow the three step-path `except Exception` blocks** to the exceptions
  they actually mean. Each currently converts an I/O problem into a silently
  weakened correctness guarantee.

---

## Threats

- **The house failure mode is a plausible number, not a crash.** Two
  independent `2**31` truncations (Known Issues #5 and #7) each reported a
  confident, wrong, low figure; the inert-config class (F1, F2, F1b, and the
  four keys fixed on 2026-08-31) silently ignores what the user asked for; the
  swallowed guards above turn an I/O error into a passed screen. None of these
  announce themselves. A workflow whose errors look like results is the hardest
  kind to catch in review and the most damaging to publish from.

- **A result can look reproducible without being reproducible.** A shared
  `data_dir`, no staleness detection, no failure record, and a one-shot command
  that cannot be seeded compose into a specific risk: re-running the pipeline
  and getting the same numbers proves less than it appears to, because the
  inputs may not be the ones the manifest describes. This is the threat that
  matters most for a paper, and it arises from the interaction of the weaknesses
  above rather than from any single one.

- **Provenance that is trusted further than it is true.** `effective_conditions`
  is genuinely good, which lends credibility to `resolved_params` beside it,
  which is not what its name says — and to a manifest that reads as a complete
  history while recording only successes. Partial provenance presented without
  qualification is worse than none, because it is relied upon.

- **Surface growing faster than the evidence that it works.** 38 commands, and
  this session alone found one crashing flag, six inert ones, five inert config
  keys and a dead validator. `tests/test_design_options_have_effect.py` exists
  precisely to catch this class and, until today, covered almost none of it.

---

## What to do first

Ranked by risk removed per unit of work, not by severity:

1. **F1b `optimization_method`** — a documented config key that pins params.json
   users to a method measured at up to 260x the cost for an identical set. Same
   fix as the four keys already done, plus a "was the flag given" sentinel,
   because `'hybrid'` is both the default and a legitimate explicit choice.
2. **`--output-prefix` plus a `params_hash` check on read** — together these
   close the clobbering and the staleness holes, which are the two that make a
   result unreproducible without saying so.
3. **Manifest on failure, and a params.json checksum** — turns the provenance
   record into something that can be relied on, including when it matters most.
4. **Narrow the three swallowed guards** in the step path.
5. **`design`: add the flags or document the recipe.** Cheap either way, and it
   removes a trap from the command most likely to be a user's first.

Items 1, 3 and 4 are small and independent. Item 2 is the one worth planning.

---

## Method

Behaviour was checked by running it wherever a run was cheap: the
`optimization_method` finding comes from calling `get_params` on a params.json
that sets it and reading back the module global, not from tracing the code. Flag
surfaces were enumerated from the live argparse tree rather than from the
source, which is how the `design` gaps surfaced. Counts (`134`, zero call sites)
are `grep` totals over `neoswga/` excluding caches, re-run at the time of
writing. Claims that could not be settled cheaply are not in this document.
