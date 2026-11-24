# Tier-3 defects from the 2026-11-23 sweep, and their fixes

Date: 2026-11-24. Base: `bd14b77` (main after #153). Branch: `tier3-bug-sweep`.

## Source

The tier-3 list in `docs/validation/2026-11-23-tier1-bug-sweep.md` ("Not
addressed"), and five of the six minors the final review of tier 2 deferred
(`docs/validation/2026-11-24-tier2-bug-sweep.md`). Each fix has a test that
failed before the change.

| # | Defect | Fix | Test |
|---|---|---|---|
| 1 | The schema documented `min_tm`/`max_tm` as 15/45 and `dtt_mm` as 0.0; all three, and `reaction_temp`, resolve from the polymerase. A copied `dtt_mm: 0.0` reads as DTT omitted | no fixed default; "Polymerase-aware default is used if absent" | `tests/test_polymerase_derived_keys_document_no_fixed_default.py` |
| 2 | `filter --max-gini` help said 0.6 (the default is 0.7); `--min-tm`/`--max-tm` said 15/45 | help states the resolved default | `tests/test_help_states_the_resolved_default.py` |
| 3 | `--bl-penalty` on `count-kmers` and `filter` wrote a global nothing reads (the key was retired 2026-09-14) | flag removed; argparse refuses it | `tests/test_bl_penalty_is_retired.py` |
| 4 | `show-presets` printed an optimization method and set size for four presets that `--preset` never applies, and listed the other seven by a second path | one `PRESET_NAMES` tuple feeds `filter --preset` and the listing, which reads each preset through `load_preset_conditions` | `tests/test_show_presets_matches_what_preset_applies.py` |
| 5 | The option ratchet excused 40 options as unread. Fourteen were read: `-j` on three commands through `@params_command(seed=True)`'s default merge; nine `rescore-set` chemistry flags plus `evaluate-set --reaction-temp` through `build_reaction_conditions`, which reads by computed name; and `optimize --polymerase` through a module-level alias (`_polymerase_from_params = polymerase_from_params`). Of the rest, most read as working: `count-kmers --polymerase bst` counted at the params.json enzyme's k range | the walker sees all three routes and those entries are gone; the 17 options below refuse when given and name the params.json key that does what they say; the nine whose help already says they do nothing stay accepted | `tests/test_every_cli_option_has_an_effect.py` (`test_an_excused_option_refuses_or_says_it_does_nothing`, `test_every_refused_option_is_one_no_path_reads`) |
| 6 | `design` had no `-n`, `-m` or `--seed`; step 3 seeded from `--seed` alone | the three flags, `-m` with the optimize parser's choices; step 3 uses `resolve_seed` | `tests/test_design_takes_the_optimize_size_seed_and_method.py` |
| 7 | `expand-primers` suggested `interpret`/`report`/`export -d OUT`, which read an optimize result it does not write | it also writes `expanded_primers.txt` and suggests `evaluate-set --primers-file` on it | `tests/test_expand_primers_advice_runs.py` |
| 8 | BED/BedGraph export globbed every `*_positions.h5` in the directory (host indexes included, under the target's name), ignored `-j`, wrote no line for a primer the index lacked, and used length 0 when `fg_seq_lengths` was absent | target from `-j`, else `DIR/params.json`, else `./params.json` when its `data_dir` is `DIR` (the `init` layout); paths read relative to the working directory, as the pipeline reads them; unindexed primers scanned from the target FASTA or refused; a second target refused; length recorded, else measured, else refused | `tests/test_bed_export_reads_the_run.py` |
| 9 | `export --no-modifications --pto-bonds 3` dropped the bonds silently | refused | `tests/test_export_modification_flags_conflict.py` |
| 10 | `validate` reported jellyfish alone, so a KMC3-only machine read as unable to count | both counters, and the one counting uses | `tests/test_validate_reports_both_counters.py` |
| 11 | Tier-2 minor: `simulate` took only polymerase and temperature from the record, and logged the enzyme optimum | the whole recorded reaction, through the report's reconstruction; refuses when it cannot be rebuilt | `tests/test_simulate_cli_inputs.py` |
| 12 | Tier-2 minor: with no panel limit set, an alternative whose coverage could not be measured was dropped | checked only when a limit is configured | `tests/test_alternatives_respect_panel_limits.py` |
| 13 | Tier-2 minor: the coverage-reach ratchet in `test_byo_oligos.py` passed on a docstring | reads the AST; carries a prose-only module it must fail | `tests/test_byo_oligos.py` |
| 14 | Tier-2 minor: `plan-pool --design-grid` was refused for a multi-length k range and needed a single-length pool and index it never used; a configured method plan-pool cannot run fell back silently | the grid runs before those checks; the fallback is logged | `tests/test_plan_pool_defaults.py` |
| 15 | Tier-2 minor: with no measured Gini, the set's mean and max Gini were 0.0 and the grade called uniformity Excellent | None; the grade leaves the component out and renormalises; both reports say "not measured" | `tests/report/test_unknowns_render_as_unknown.py` |

## Options that now refuse

`count-kmers`: `-z/--data-dir`, `--polymerase`, `-x/--fasta-fore`,
`-y/--fasta-back`, `-k/--kmer-fore`, `-l/--kmer-back`. `filter`: `-z`,
`--min-gini-sites`. `prepare-candidates`: `-z`, `--polymerase`. `optimize`: `-z`.
`design`: `-z`. `evaluate-set`: `-z`, `--enable-qa`. `analyze-dimers`:
`--visualize`. `analyze-set`: `--simulate`. `analyze-stability`: `--temp`.
Plus `--bl-penalty` (item 3). Each was accepted and changed nothing; a script
passing one now stops with exit 2 and the route to use. The Task-5 commit
message gives 22 for this list; that count predates finding five of them read.
`-z` on `evaluate-set` had one effect, the directory a failure record is
written to; it is refused with the rest.

## Other behaviour that changes for a user

- `design -n/-m/--seed` move the request hash exactly as on `optimize`.
- `prepare-candidates` seeds from params.json's `seed` when `--seed` is absent.
- `simulate --from-results` simulates the recorded additives and buffer.
- The grade of a set with no measured Gini is computed over the other four
  components.
- BED/BedGraph export refuses a run with more than one target.

## Final review

A fresh review found that `optimize --polymerase` is read (through a
module-level alias the option walker did not follow) and had been refused; the
walker now follows aliases and the flag is accepted again. It also found BED
export resolving paths against the params file's directory, which the pipeline
does not, and the JSON report writing `NaN` for an unmeasured uniformity; both
fixed. `optimize --polymerase` sets the coverage reach and `--auto-size`
chemistry while Tm and dimer conditions still come from params.json; that
half-application predates this branch and is not addressed here.

## Not done

- `design --preset`: a preset reaches step 2 only (steps 3 and 4 warn on drift
  from it), so on a one-shot command it would read as configuring the whole
  design. The route is the preset's values in params.json.
- Circular wrap-around in `_common.design_position_cache`: `PositionCache` has
  one circular flag for all genomes, and `fg_circular` and `bg_circular` can
  differ. Needs a per-genome flag.
- `analyze-stability --temp` was refused rather than wired: the analyzer takes a
  `ReactionConditions`, which needs a polymerase the command does not ask for.
