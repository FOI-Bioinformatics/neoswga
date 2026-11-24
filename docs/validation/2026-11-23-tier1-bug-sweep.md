# Tier-1 defects from the 2026-11-23 sweep, and their fixes

Date: 2026-11-23. Base: `5343bc8` (main after #151). Branch:
`tier1-bug-sweep`.

## Source

Six read-only reviews of the CLI arguments, the gaps between pipeline steps,
silent substitution of unknowns, params.json resolution, the optimize step,
and the commands that consume a delivered set. About 60 findings were
reported. The twelve below were rated tier 1: a stale or failed run presented
as current, two commands resolving one params.json differently, or an
unmeasured quantity reaching selection as a zero. Each was confirmed by
reading the code, and most were reproduced on the plasmid example before
being fixed. Tiers 2 and 3 are listed at the end and are not addressed here.

## The twelve, with the test that pins each

| # | Defect | Fix | Test |
|---|---|---|---|
| 1 | `design` caught every error and exited, so a `DesignError` wrote no `design_failure.json` and `export` shipped the previous panel | `run_design` re-raises; the record names the step `design` was running | `tests/test_failure_records_name_the_step.py`, `tests/integration/test_strict_design_pipeline.py::test_a_failed_design_rerun_blocks_the_previous_panel` |
| 2 | `design` back-filled `mechanistic_weight=0.3`, which `run_step4` reads as a request for the mechanistic model; `optimize` on the same file left it off | back-fill is None | `tests/test_design_optimize_parity.py` |
| 3 | `interpret` and `report` did not read the failure record or the set index; both called a set ready that `export` refused | both read `export_is_blocked`; the report keeps its grade and changes only the recommendation | `tests/test_interpret_and_report_agree_with_export.py` |
| 4a | Only a `DesignError` wrote a failure record; `exit_on_step_failure`, refused prerequisites, `ValueError`s and schema refusals wrote none | every nonzero exit of a design step writes a record naming the step; `prepare-candidates` and `optimize` refuse after an earlier step's record; a step's success clears only its own | `tests/test_failure_records_name_the_step.py` |
| 4b | The manifest recorded step-2 and step-3 checksums and nothing compared them; `optimize` ran on a pool prepared from an earlier `filter` | `optimize` refuses a demonstrated mismatch and warns when no entry describes the pool | `tests/test_step3_provenance.py` |
| 5 | `get_primer_positions` returned an empty array on any exception, so an unindexed or unreadable index counted as zero sites | raises `ReferenceDataError` | `tests/test_position_lookup_is_not_zero.py` |
| 6 | `get_params` defaulted the k range from the polymerase, `resolve_design_request` to 6-12; an equiphi29, bst or klenow file without k passed `filter` and was refused by `optimize` | the request defaults from the same registry | `tests/test_k_range_defaults_agree.py` |
| 7 | Absent `min_fg_freq`, `max_bg_freq`, `max_gini`, `max_primer` resolved to None (`filter` crashed); absent `max_sets` delivered 1 set against a hashed 5, `iterations` 100 against 8 | the seven keys take `PipelineParameters`' default; the schema states the four it did not | `tests/test_absent_keys_take_documented_defaults.py` |
| 8 | `count-kmers` and `design` required `jellyfish` even with KMC3 installed | the pre-flight resolves the counter as counting does | `tests/test_kmer_counter_check.py` |
| 9 | `network` used a 70 kb connectivity reach for every polymerase | resolved from the reaction's polymerase, as on `hybrid` | `tests/test_optimizer_config_reaches_optimizers.py::test_network_reach_follows_the_polymerase_when_unset` |
| 10 | Vendor CSV length, Tm and GC were computed on the sequence with modification syntax (`TGTACTGCC*A*A`: 13 nt) | computed on the bare primer | `tests/test_export_vendor_csv_figures.py` |
| 11 | `export` and `report` took conditions from the latest manifest entry of any step | the `optimize` entry | `tests/test_export_conditions_from_optimize_step.py` |
| 12 | `hybrid`, `background-aware`, `network`, the composite member loop and the ensemble union-combine turned a `DesignError` into a failure result | re-raised | `tests/test_optimizers_propagate_design_errors.py` |

## Final review

A fresh review of the branch found two defects in the new record rules, fixed
before merge with tests in `tests/test_failure_records_name_the_step.py`:

- A record naming `design` (written before records named steps, or by a
  `design` refused before any step ran) or an unreadable record blocked
  `prepare-candidates` and `optimize` and nothing could clear it. A `design`
  that ran no step now writes none, and `count-kmers` success clears both.
- The checks resolved `data_dir` against the params file's directory and the
  run resolves it against the working directory, so with a params file in
  another directory the check read the wrong one; and `--data-dir`, which the
  steps ignore, placed the record. Both now use the run's resolution.

## Behaviour that changes for a user

- A params file without `max_sets` now receives up to five alternative sets,
  as documented. Set 0 is unaffected.
- A params file with `polymerase` equiphi29, bst or klenow and no k range now
  reaches `optimize` instead of being refused there.
- A `network` design with a non-phi29 polymerase and no `max_extension` now
  scores connectivity at that polymerase's processivity (bst: 2 kb). phi29 is
  unchanged. The effect on a delivered panel was not measured.
- After any failed design step, `export`, `interpret` and `report` refuse or
  say "DO NOT ORDER" until the step is re-run successfully, and
  `prepare-candidates`/`optimize` refuse after a failed earlier step.
- `optimize` refuses a `step3_df.csv` that its manifest shows was prepared
  from a different `step2_df.csv`.
- A schema refusal of `optimize` now leaves a failure record. It left none
  before, which let an earlier panel in the directory remain exportable.

## Suite

`pytest -n 8 tests/`: 7733 passed, 7 failed. Six are `test_docs_links_resolve`
on gitignored local files. One
(`test_the_host_term_has_one_calibration::test_worst_case_with_background_aware_and_two_hosts_is_refused`)
fails on `main` too, when it runs after `tests/test_polymerase_defaults.py`,
which leaves `parameter.polymerase = "custom_enzyme"`.

Serial, in six interleaved chunks: the only failures besides the link checks
were `test_validate_smoke.py::test_smoke_runs_all_four_steps_and_passes` and
seven tests in `test_application_selection_weights.py`. The same chunks fail
identically on `main` (a leaked `parameter.data_dir` pointing at a deleted
temporary directory, and a leaked reaction fingerprint). Every one passes on
its own.

## Not addressed (tiers 2 and 3)

Tier 2 was addressed on 2026-11-24; see
`docs/validation/2026-11-24-tier2-bug-sweep.md`. Tier 3 remains.

Tier 2: a corrupt `run_manifest.json` is erased on the next append;
`step2_df.csv`/`step3_df.csv` are written non-atomically; `simulate
--duration` is inert and `--polymerase` defaults to phi29; `evaluate-set`
ignores `fg_circular`; `rescore-set` ignores `coverage_reach` and reports a
`fg_cov*1000` selectivity; unindexed primers become whole-genome gaps in
`analyze-coverage` and `expand-primers`; `plan-pool` defaults beat params.json;
`seed` in params.json seeds nothing and is missing from the schema with
`fixed_oligos`/`excluded_oligos`; the report renders a NaN Gini as 0 and the
1e6 selectivity sentinel as an enrichment in the technical report; `export`
drops `na_conc`/`k_conc`/`nh4_conc`/`dntp_conc`/`dtt_mm`; alternative sets are
not checked against panel limits; `analyze-set` accepts non-ACGT input;
`-n 0` and negative `-n`; several search flags are outside the request hash;
`resolve_design_request` defaults `min_tm`/`max_tm` to 20/50 where `get_params`
derives them from the polymerase; `design -z DIR` writes a failure record into
`DIR`, which the run itself ignores.

Tier 3: documented `min_tm`/`max_tm`/`dtt_mm` defaults; `--max-gini` help
text; inert `--bl-penalty`; `show-presets` optimizer hints; stale entries in
the option ratchet's allowlist; `design` lacks `--seed`, `--preset`, `-n`;
`expand-primers` next-step advice; BED export ignores `-j`;
`--no-modifications` drops `--pto-bonds`; `validate` reports jellyfish only.
