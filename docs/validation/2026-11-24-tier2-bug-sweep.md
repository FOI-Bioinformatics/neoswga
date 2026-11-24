# Tier-2 defects from the 2026-11-23 sweep, and their fixes

Date: 2026-11-24. Base: `d633e3b` (main after #152). Branch:
`tier2-bug-sweep`.

## Source

The tier-2 list in `docs/validation/2026-11-23-tier1-bug-sweep.md` ("Not
addressed"), and two minors the tier-1 review deferred (M3, M5). The
`design -z DIR` item there had already been fixed by the tier-1 review (8c463ad).
Each fix has a test that failed before the change.

| # | Defect | Fix | Test |
|---|---|---|---|
| 1 | An unparsable `run_manifest.json` was replaced by one entry, erasing every earlier step and its `effective_conditions` | kept as `run_manifest.json.corrupt-<time>`, with a warning | `tests/test_run_manifest_corruption.py` |
| 2 | `step2_df.csv` and `step3_df.csv` were written with a plain `to_csv`; an interrupted write left a short file step 4 accepted | written through `atomic_output.atomic_write_dataframe` | `tests/test_step_tables_are_atomic.py` |
| 3 | `simulate`: the `--duration` config was never passed (every run 8 h); `--polymerase` defaulted to phi29; lower-case primers found no sites | config passed; polymerase from the flag, else the reaction the optimize step recorded, else phi29; primers upper-cased and refused unless DNA | `tests/test_simulate_cli_inputs.py` |
| 4 | `evaluate-set` ignored `fg_circular` (every site read `not args.linear`) | `--linear`, else params.json, else circular | `tests/cli/test_evaluate_set_geometry.py` |
| 5 | `rescore-set` ignored `coverage_reach` and reported `fg_cov * 1000` as the selectivity when host coverage was 0, including with no host | reach through `resolve_coverage_reach`; selectivity None without a host, `MAX_SELECTIVITY` at zero host coverage | `tests/test_rescore_set_reach_and_selectivity.py` |
| 6 | `analyze-coverage` and `expand-primers` read a primer missing from the index as binding nowhere (one whole-genome gap for a primer with two sites) | the cache scans the run's genomes for missing primers and refuses without them | `tests/test_unindexed_primers_are_scanned.py` |
| 7 | `plan-pool`: `--primer-length` 12, `--coverage-metric` effective and the method fall-through beat params.json | flag, else params.json, else the old default; a multi-length k range needs `--primer-length` | `tests/test_plan_pool_defaults.py` |
| 8 | `seed` in params.json was hashed into the request and seeded nothing; `seed`, `fixed_oligos`, `excluded_oligos` were reported as "will be ignored" | step 4 seeds from `--seed`, else the file; request-only keys are known keys | `tests/test_params_seed_and_request_keys.py` |
| 9 | The report rendered a NaN Gini as 0.0 (uniformity 1.0), the 1e6 selectivity sentinel as "1000000x", and a genome-length placeholder as a measured strand gap | Gini None and dropped from the contribution score; "no background binding detected"; no gap note without a measured strand ratio | `tests/report/test_unknowns_render_as_unknown.py` |
| 10 | `export`, `report` and `interpret` dropped `na_conc`, `k_conc`, `nh4_conc`, `dntp_conc`, `dtt_mm`; `interpret` substituted 50 mM Na+ and 2.5 mM Mg2+ | the recorded buffer reaches `ReactionConditions`; `interpret` uses the report's reconstruction | `tests/test_export_carries_the_buffer.py` |
| 11 | An alternative set that missed a configured panel limit was written to `step4_improved_df.csv` | not delivered, with a log line naming the limit; dimer checks left to the search | `tests/test_alternatives_respect_panel_limits.py` |
| 12 | `analyze-set` printed thermodynamics for `ACGTXXQQ`; the set-editing loader kept duplicates and accepted anything | `validate_primer_sequence`, duplicates dropped | `tests/test_primer_lists_are_validated.py` |
| 13 | `-n 0` read as no flag; `-n -3` reached the parameter module | a `positive_int` argparse type | `tests/test_num_primers_flag_is_positive.py` |
| 14 | Search flags with no params.json key (`--minimize-primers`, `--application`, `--mechanistic-weight`, `--no-background`, ...) never reached the request hash | recorded in `DesignRequest.search_flags` when set to a non-default value, hashed only when non-empty; every optimize flag classified, checked against the live parser | `tests/test_a_flag_that_changes_a_design_moves_its_identity.py::TestSearchFlagsMoveTheHash` |
| 15 | The request's Tm window was fixed at 20-50 C | `parameter.default_tm_range(polymerase)`; phi29 unchanged | `tests/test_k_range_defaults_agree.py` |
| 16 | M5: `_dispatch_optimizer` turned a `DesignError` raised while building an optimizer into a failure result; M3: with no counter installed the refusal named only KMC | re-raised; the message names both | `tests/test_optimizers_propagate_design_errors.py`, `tests/test_kmer_counter_check.py` |

## Behaviour that changes for a user

- `simulate --from-results DIR` simulates the recorded reaction rather than
  phi29 at 30 C, for the `--duration` asked for.
- `evaluate-set -j params.json` follows the file's `fg_circular`.
- `rescore-set` reports `selectivity_ratio: null` with no host configured.
- `analyze-coverage` and `expand-primers` read the FASTA for primers the index
  does not hold, and refuse when no genome is available.
- `plan-pool` refuses a params file whose k range spans several lengths unless
  `--primer-length` is given.
- A params file's `seed` now seeds `optimize`.
- Alternatives that miss a configured panel limit are no longer written. With
  no limits configured nothing changes.
- A run that sets any of the search flags records a different `request_hash`
  from one that does not; a run that sets none keeps its hash.
- The request for bst, equiphi29 or klenow records the polymerase's Tm window,
  so the hash of such a run with no `min_tm`/`max_tm` changes; phi29's does not.
  `plan-pool` and the library configs now optimize under that window too.
- `plan-pool` without `--method` follows params.json's `optimization_method`;
  the Wolbachia example's README now passes `--method background-aware`, which
  is what its saved result used.
- `simulate` with no `--duration` now simulates its documented 60 minutes; it
  ran 8 hours before. Default coverage figures change.
- A params.json `seed` that is not an integer is refused rather than hashed and
  ignored; the manifest records a seed taken from the file.

## Final review

A fresh review found four defects in this branch, fixed before merge
(`tests/test_tier2_review_fixes.py`): the interactive report crashed on an
unknown Gini; `analyze-coverage` without primers scanned every genome for a
placeholder or refused naming it; `plan-pool` and the library configs kept the
20-50 C Tm window; and the Wolbachia example no longer reproduced its result.
Two further points were fixed with them: the manifest seed, and the set-editing
commands now refuse ambiguity codes (A, C, G, T only), which would score as zero
sites.

## Suite

`pytest -n 8 tests/` with a private `--basetemp`: 7813 passed. The failures were
`test_docs_links_resolve` on gitignored local files, and integration tests that
failed with `ENOSPC` while the data volume was full and pass when run again
serially.

Serially, in six interleaved chunks: the only other failures are tests that
fail after an earlier test leaves `parameter.data_dir` pointing at a deleted
temporary directory (the 2026-10-02 defect list, item 4). For the two chunks
that showed them, the same files run on `main` failed on exactly the same
tests.
