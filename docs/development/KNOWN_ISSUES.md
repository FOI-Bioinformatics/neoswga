# Known issues

> Moved verbatim from `CLAUDE.md` on 2026-10-01, when that file was reduced to a
> working summary. The numbering is stable: code comments and tests
> cite "Known Issue N", so do not renumber. Entries are dated records: a figure is what was
> measured on the date given, on the references named, and a later entry may
> correct an earlier one. `CLAUDE.md` keeps the short form and links here.

Most entries are fixed defects kept for the mechanism and the test that now
guards it. The unnumbered entries between 7 and 8 were recorded in that
position and are left there.

1. **Large background genomes**: Use `neoswga build-filter` to pre-build a Bloom filter for human genome. Note this is not always needed: at k=12 exact jellyfish counting of the whole human genome costs about 7 minutes and a 138 MB table (8,368,418 canonical 12-mers), which is well within reach. The Bloom path matters at longer k, where the count table stops being small. Measure before reaching for it -- the sampled-index path has its own resolution trap (see `_warn_if_sample_too_sparse`).

2. **sklearn compatibility**: The RF model ships in skops format (version-tolerant, no arbitrary-code deserialization), so minor sklearn upgrades no longer require retraining. A major sklearn upgrade may still warrant re-validating the model.

   **The format is version-tolerant; its default trust list is not**, and the
   two are easy to conflate. skops 0.15.0 stopped implicitly trusting
   `sklearn.tree._tree.Tree`, and every model-loading test went red on CI while
   passing locally on 0.14.0 -- `pip install -e ".[dev]"` resolves
   `skops>=0.11,<1` to whatever is newest, and `requirements-dev.lock` is not
   used by the test job.

   `rf_preprocessing._TRUSTED_MODEL_TYPES` now names the types a forest
   legitimately needs and `unexpected_model_types` refuses anything else,
   naming it. That is narrower than trusting the archive wholesale and does not
   depend on a future release keeping today's defaults; pinning skops would
   have worked until the next release did the same thing. The digest check
   against `models/checksums.json` still runs first, so this narrows what an
   already-vouched-for file may reconstruct rather than replacing provenance.

   The refusal rule is a pure function because which types skops reports as
   untrusted depends on the installed version: a test driving the loader would
   exercise it on 0.15 and skip straight past it on 0.14.

3. **Memory usage**: The filter command loads all background k-mers into memory. Use Bloom filter for large backgrounds.

4. **PositionCache strand parameter**: Uses 'forward', 'reverse', 'both' (not '+' or '-').

5. **pyahocorasick silently finds nothing past 2 Gb** (upstream, unfixed): `Automaton.iter()`
   indexes with a 32-bit int, so for a string longer than `2**31 - 1` its scan loop never runs.
   It yields nothing, raises nothing and warns nothing. Verified on **2.3.1, the latest release
   as of 2026-08**, with the needle planted at offset 1000: length `2**31 - 10` finds it,
   `2**31 + 10` does not. Not found in the upstream tracker.

   The re-check is automated rather than left to memory. `tests/test_pyahocorasick_limit_canary.py`
   fails as soon as the installed version leaves the set that has actually been measured, and
   names the command that re-measures it:

   ```bash
   NEOSWGA_VERIFY_AHOCORASICK_LIMIT=1 pytest tests/test_pyahocorasick_limit_canary.py -k still_present
   ```

   That run allocates over 2 GB, so it is opt-in. If the limit is gone, the chunking is still
   correct but no longer load-bearing, and keeping or dropping it is a deliberate choice.

   `string_search.MAX_SCAN_CHUNK` (2**30) works around this by scanning in overlapping windows,
   so this is **load-bearing, not defensive** -- removing the chunking silently breaks any
   background above 2.147 Gb, which includes human (3.1 Gb) and mouse (2.7 Gb).

   Symptom before the fix: `total_bg_sites` and `bg_coverage` read 0 for a whole-genome host,
   which is indistinguishable downstream from a perfectly specific primer set. A 27-primer panel
   scored 0 against hg38 where the jellyfish counts put the true figure at 860.

6. **A partial background is not a specific design**: `selectivity_ratio` is a ratio of counts
   with no genome length in it, so it moves with how much background sequence you supply --
   about 66x between human chr21 and whole hg38, with nothing about the primers changed. Read
   `selectivity_density` (added beside it) when comparing designs scored against different
   backgrounds. At k=12 note that 99.7% of all canonical 12-mers occur in the human genome, so
   a gate demanding zero host sites returns a single-site candidate pool rather than failing.
   See [docs/validation/additive_specificity.md](../validation/additive_specificity.md#the-background-was-one-chromosome-and-that-mattered).

7. **Genome coordinates are int64** (`position_cache.POSITION_DTYPE`), and must stay so.
   Positions were cast to `np.int32` on the way out of HDF5. int32 tops out at 2,147,483,647,
   so for human (3.1 Gb) and mouse (2.7 Gb) every site past that offset **saturated at the
   ceiling**; the `strand="both"` path then calls `np.unique`, which collapsed all of them into
   a single site.

   This is a second, independent instance of the 2**31 failure in Known Issue #5 -- different
   place, same symptom, and neither one is visible on a single chromosome. Measured on an
   *M. tuberculosis*-vs-hg38 run: `CACCGACGACGA` occurs 48 times in hg38 (jellyfish and a
   direct string count agree), the scan stored all 48 correctly, and the cache returned 5.
   Across the twelve-primer set `total_bg_sites` read 48 against a true 114.

   It changed the design, not just the report: with the true background visible the optimizer
   **drops** `CACCGACGACGA`, the primer whose host load int32 had been hiding. After the fix
   `total_bg_sites` matches the jellyfish count exactly (77 = 77).

   The lesson both issues share: **test against a whole genome, not a chromosome.** chr21 is
   46 Mb and cannot reach either limit, so both bugs sat behind a passing test suite.
   `tests/test_position_cache.py::TestPositionsPastTheInt32Ceiling` pins this one.

**`--design-grid` on `plan-pool`** designs once per condition and length in a
JSON grid and writes `design_sweep.json` beside the usual report, rather than
one `pool_plan`. The grid names `lengths` and `conditions`, where each condition
names only the fields it changes: the baseline is the reaction this run
resolved, so a grid varying DMSO alone keeps the buffer, salts and oligo
concentration, and the comparison is between chemistries rather than against
library defaults. A cache and optimizer are rebuilt per condition and length,
since the index is per length and the chemistry is what varies.

It needs the candidate inventory, and it looks each condition up by reaction
fingerprint, so a condition the filter never recorded is reported as having no
eligible candidate rather than designed with an empty pool. Wired on 2026-09-17
in Phase 4 increment 6; it was audit finding F4, parsed and documented and read
by nothing, and its entries are now gone from both the inert-option and
unreachable-capability allowlists.

**`--min-fg-bg-ratio` was read and then overruled** -- FIXED 2026-09-17.
`optimize`'s background prefilter kept every candidate at or above the ratio,
then, if that removed more than `max_removal_fraction` of them, discarded the
threshold and kept the top 80% by ratio instead. On the 2,000-candidate
Wolbachia shortlist the threshold removes 64.8% at its default of 1.0, so the
clause fired at 1.0, 2.0, 5.0 and 20.0 and removed exactly 400 every time. The
flag changed nothing above about 1.0 and the rule in force was "drop the worst
20%".

That is the Known Issue 8 class in a shape none of its ratchets look for: not a
flag nobody reads, but a flag that is read and then overruled by a second rule
on the same decision. `max_removal_fraction` was also a bound on the fraction
of a BATCH, so which candidates survived depended on how many others were below
the threshold alongside them.

`order_candidates_by_background` replaces it. Candidates at or above the ratio
are searched first and the rest are searched last; nothing is deleted, so the
400 the old path made unreachable at every setting are reachable again, which
matters because increment 5's refill can now reach them. The partition is
stable, preserving the inventory's `search_rank` traversal. Delivered panel on
the measured design: 11 of 12 primers shared, Jaccard 0.846
([measurement](../validation/background_ordering_2026-09-17.md)).
`bg_max_removal` is retired with the clause.

**The objective never reached the stage that refines** -- FIXED 2026-09-18.
`plan_pool` attached `pool_objective` to the optimizer it was handed, which on
every command-line path is a wrapper (`HybridBaseOptimizer` or
`BackgroundAwareBaseOptimizer`) that delegates the search to an inner
`HybridOptimizer`. `_swap_refine` is a method of the INNER one, so
`refine_hybrid_stage2` read the attribute off an object nobody had set. Measured
through a real design, the refinement ran once and received None: Stage 2
refined on raw covered bases while the row was accepted on occupancy-weighted
coverage under a specificity floor. Delivered density on a failing row went from
28.78 to 42.62 once connected.

Two tests covered it and neither could see it -- one asserted by AST that
`plan_pool` assigns an attribute of that name, the other by source text that the
refinement reads one. Both ends existed and the path did not. Use
`swap_refinement.attach_search_config` for anything a delegate must read, and
assert the PATH: `tests/test_the_objective_reaches_the_stage_that_refines.py`
drives a real factory-built optimizer under both methods.

**Stage 1 is deliberately NOT constraint-aware**, and there is no demonstrated
reason for it to be. The specificity density floor is exactly additive over
primers, so a density-only ceiling is computable and comes to 79.807 for a
12-primer panel on the shipped pool. **That figure bounds nothing deliverable**:
the panel achieving it has coverage 0.4042 against a 0.5 target, and the
constructions that appear to beat the search carry 30 dimerising pairs out of 66
at the configured `max_dimer_bp` of 3. Only 6 of the 16 most selective
candidates are mutually compatible, so a 12-primer panel cannot be built from
them at all.

With dimers and coverage both in force, no deterministic construction beats the
search: the best reach density 50.6 at coverage 0.527, or 64.0 at coverage 0.385
which fails the target, against the search's **60.112 at 0.6535**. So its
failure at a floor of 65 is probably correct. Three Stage 1 rules built on the
accounting each failed to improve a delivered panel, which is best explained by
there being nothing to find.

Settled with controls on 2026-09-18: **no construction respecting the dimer
screen beats the search.** Nine of them, over three reference densities and
three pool sizes, all land between 0.375 and 0.387 coverage against a 0.5 target
and none exceeds 64.0 density. Selective candidates are GC-richer (0.44-0.48
against 0.39) and pairwise less compatible (58-63% against 69-71%), but the
largest mutually compatible subset is NOT smaller -- 18 among the top 64 by
slack, more than a 12-primer panel needs -- so compatibility is not the barrier.
The specificity against coverage trade-off is
([measurement](../validation/no_search_headroom_on_this_pool_2026-09-18.md)).

The lesson is the reusable part: **an achievability figure that omits a
constraint bounds nothing, and a striking ratio without a control is not a
finding.** Three claims of mine died in sequence here -- a 79.807 ceiling that
ignored coverage and dimers, an existence proof carrying 30 dimerising pairs out
of 66, and a "6 of 16" compatibility barrier that sits inside the random range
of 6 to 8. Three Stage 1 search rules were built on the first two before the
third was tested. The check that would have killed all three at the outset is
the same one: evaluate a candidate panel through the acceptance path a delivered
panel takes, and compare it against a control. Quote 60.112 at coverage 0.6535
as the reference for this pool. The accounting lives in
`scripts/benchmarking/selectivity_budget.py` as a diagnostic, not in the
package, because nothing in the search uses it.

8. **`optimization_method` in params.json did nothing** — FIXED 2026-09-05
   (audit finding F1b). The key was declared in `params.schema.json`,
   documented above, accepted by the validator, and read by nothing:
   `get_params` assigned no module global, and `run_step4` passed
   `args.optimization_method` straight through with an argparse default of
   `'hybrid'`, so the flag's default beat the config every time.

   ```
   params.json "optimization_method": "dominating-set"
     -> parameter module global: <UNSET>,  optimizer actually run: hybrid
   ```

   It cost more than provenance: `hybrid` returns a set **identical** to
   `dominating-set` (Jaccard 1.000) at 7.8x the cost at 32 primers and 260x at
   128 (2239 s against 8.6 s), so every params.json user ran the slowest method
   for the same answer. `design` pinned it a second way — that subparser has no
   `--optimization-method` and `run_design` hardcoded `"hybrid"`.

   Fixed in three places, because one alone was not enough. The global is
   assigned in `_apply_params_only_keys`. The flag's argparse default is now
   `None`, the sentinel that distinguishes an explicit `--optimization-method
   hybrid` — which must beat a configured `dominating-set` — from an absent
   flag, which must not; do not give it a real default again.
   And the lookup goes through `optimization_method_from_params`, beside the
   other pre-read resolvers, because `run_step4` builds its argument list
   before `optimize_step4` triggers `get_params`, so reading the global at that
   point sees nothing. Tests:
   `tests/test_optimization_method_routes_from_params.py`.

   It was NOT the last instance of its class -- a config key or flag that is
   documented, accepted, and read by nothing. An audit on 2026-09-14 found ten
   more. The class is now closed, and held closed by a ratchet.

   Closed on 2026-09-14 in the two ways available. Wired, because each already
   had a reader taking its fallback: `mismatch_penalty` (whose consumer
   `occupancy.default_mismatch_penalty` was written for it and received None on
   every call), `max_homopolymer_run`, `gc_clamp_window` and `max_gc_in_clamp`.
   Retired from the schema, because nothing implemented what they named:
   `retries`, `drop_iterations`, `top_set_count`, `selection_metric` and
   `bl_penalty`. The first four appeared only in a module-level `defaults` dict
   in `core/pipeline.py` that itself had no reader; the dict is gone. Setting
   any of the five now produces the unknown-key warning rather than silence.

   `tests/test_no_schema_key_is_inert.py` is the ratchet: every schema key must
   bind a `parameter` global or appear on a short list of keys consumed during
   loading, each with its reason.

   Fixed on 2026-09-14: `filter --gc-tolerance`, `filter --excl-threshold`,
   `expand-primers --optimization-method` and
   `plan-pool --swap-max-evaluations` all carried a real argparse default and so
   beat params.json on every run. They now use the `None` sentinel, and
   `expand-primers` routes through `resolve_optimization_method` rather than
   reading the attribute. `--gc-tolerance` was the costly one: the block it fed
   also computed its own GC window, clamping the lower bound at 0.20 where
   `adaptive_gc_window` releases it to zero below the extreme-AT threshold, so
   on a 19% GC target it excluded exactly the zero-GC primers published AT-rich
   designs are built from. It now routes through `adaptive_gc_window`.

   `tests/test_design_options_have_effect.py`,
   `tests/test_params_json_routes_optional_keys.py` and
   `tests/test_optimizer_config_reaches_optimizers.py` are named as the tests
   that hold the line, but between them they cover about fifty keys and none of
   the ten above -- which is why the class survived being declared closed.
   `tests/test_cli_defaults_do_not_beat_params_json.py` now pins the four flag
   defaults. Extend all four when adding an option.

   **Declared closed twice, and closed neither time.** The audit of 2026-09-16
   found `--design-grid` on `plan-pool` parsed, documented in the help text, and
   never read -- added *after* the second closure. Checking for more found
   `--data-dir` inert on `count-kmers`, `filter`, `prepare-candidates`,
   `optimize`, `design`
   and `evaluate-set`, and `--min-gini-sites` inert on `filter`, which the Key
   Parameters section above documented as working. Both were verified by
   resolving a config whose flag value differed from the file value: the file
   won each time.

   The reason none of the four ratchets caught them is structural.
   `test_no_schema_key_is_inert.py` and `test_params_json_routes_optional_keys.py`
   iterate params.json **schema keys**, so a CLI flag is invisible to them.
   `test_cli_defaults_do_not_beat_params_json.py` checks four argparse
   **defaults** and asserts nothing about whether a flag is read.
   `test_design_options_have_effect.py` calls `run_optimization` directly, so it
   covers the **optimize path only**. None of them asks "does this flag do
   anything".

   `tests/test_every_cli_option_has_an_effect.py` now does. It reads the dispatch
   table out of `main()`, walks from each handler through the functions it calls,
   collects every attribute read off the argparse namespace -- including the
   `merge_args_to_parameter` and `@params_command(merge=...)` routes, which are
   reads by another name -- and fails on any declared option it cannot account
   for. Currently-inert options are listed in `KNOWN_INERT` with a reason, and a
   second test fails on an entry that has since been wired, so the list can only
   shrink.

   The same defect exists one layer down, where a capability is built and tested
   and no command can reach it. The audit found six at once, all from the
   condition-aware pool design work: `CandidateProvider`, `design_sweep`,
   `load_design_grid`, `load_grid_file`, `ensure_positions` and, in practice,
   `beam_search`. Unit tests cannot see this, because a test that constructs the
   thing directly and asserts it behaves passes whether or not anything calls it.
   `tests/test_no_capability_is_unreachable.py` walks transitive reach from the
   dispatch table -- not bare references, since `design_sweep` calling
   `provider.expand` must not make `expand` count -- and holds the same kind of
   shrinking allowlist.

   Do not record this class as closed again. Record what the ratchets cover.

   A third variant surfaced on 2026-09-17, and none of the five ratchets looks
   for it: not a CLI flag and not a schema key nobody reads, but a schema key
   read in some places and not in the one that would spend it. `cpus` reached
   `create_pool` and did not reach step 1. `kmer_counter.run_jellyfish` declares
   `cpus: int = 4` and computes
   `max_workers = min(num_k, cpu_count // max(cpus, 1))`, so that default set
   both the threads per jellyfish process AND how many k values ran at once,
   while the configured value set neither. All four of step 1's call sites
   omitted it. On a 64-core machine with `cpus: 16` the run used 4 threads per
   process and 16 concurrent k values, the transpose of what was asked for.
   Fixed by passing `parameter.cpus`; pinned by
   `tests/test_the_configured_cpu_count_reaches_jellyfish.py`, which walks the
   AST of step 1 rather than asserting a thread count.

9. **Evenness is not measurable from one or two sites** -- FIXED 2026-09-10
   (audit finding B4). `filter.get_gini` keeps a primer when
   `gini.notna() & (gini < max_gini)`. The `.notna()` half was written for
   exactly the case where evenness cannot be measured, but a single-site primer
   produced 0.0, the best value available, so the guard never fired and the
   gate ranked that primer first. 86% of the shipped Prevotella-against-chr21
   pool and 96.2% of the plasmid example sat at 0.0, all with two or fewer
   foreground sites; the three whole-genome GC-tier pools have none at 0.0,
   which is why this was invisible in the runs this project usually inspects.
   `min_gini_sites` is the threshold, default 3, settable in params.json and as
   `--min-gini-sites` on `filter`; `primer_attributes.DEFAULT_MIN_GINI_SITES`
   holds the default. It is threaded into `get_gini_from_txt_for_one_k` as an
   argument rather than read from a module global, because that function runs in
   a spawned multiprocessing worker which would otherwise see the default
   instead of the configured value. `pipeline.check_gini_stage_kept_something`
   refuses to write an empty pool and names the threshold in force.
   Tests: `tests/test_gini_needs_enough_sites.py`,
   `tests/test_min_gini_sites_is_configurable.py`.

10. **The candidate loader filtered on a Tm known to be wrong** -- FIXED
    2026-09-10 (audit finding B5). `melting_temp.py:153` keeps the original melt
    package's GC-fraction bug for compatibility with the random forest retired
    on 2026-09-05. It reads 10.12 C high at k=12 (measured over 20,000 random
    12-mers, sd 0.30 C), so with the old symmetric 15 C margin the loader's
    window was `[min_tm - 25, max_tm + 5]` in true-Tm terms. Nothing was lost on
    plain phi29 and 9.6% of k=12 candidates were lost under DMSO 10% plus
    betaine 1.5 M, silently. `kmer_counter.get_primer_list_from_kmers` now calls
    the same `ReactionConditions.calculate_effective_tm` the gate calls, on the
    window `filter._resolve_tm_window` resolves, with 2 C of stated headroom.
    The shim itself is still used by `rf_preprocessing` and is correct to leave
    there: it is what the bundled model was fitted against.

11. **Near-duplicate primers reached the delivered panel** -- measured, and the
    remedy is OFF by default (audit finding A6). Delivered E. coli set 0 held 9
    pairs at Hamming distance 1 or less and 5 primers sharing the 3' hexamer
    GCGAAA. No optimizer had a similarity rejection test: the greedy picked the
    largest ABSOLUTE new coverage, so a primer with fifty sites and forty-eight
    already covered beat one with five sites all new.
    `dominating_set_optimizer.DEFAULT_REDUNDANCY_THRESHOLD` can skip a candidate
    whose covered bins are already covered above the threshold, on site sets
    rather than on sequences.

    It defaults to **1.0, which disables it**, because measurement did not
    support switching it on. Across three real pools at ten combinations of tier
    and panel size, a 0.9 threshold fired 328,846 times and changed nothing:
    coverage identical in five of six cases and 0.03 points lower in the sixth,
    with the Hamming-1 pair count and the duplicate 3' hexamer count -- the two
    things it was built to reduce -- identical in all six. The reason is
    structural: a candidate more than 90% already covered has a small marginal
    gain by construction, so the greedy's argmax was never going to pick it. The
    criterion and the objective are nearly the same signal. No threshold beats
    disabled on average, and gains sit beside large losses in the same tier:
    M. tuberculosis at n=36 gains 5.53 coverage points at 0.05 and loses 10.72
    at 0.00. The mechanism and its tests are kept; pass an explicit threshold to
    use it. It affects `hybrid` and `background-aware` too, which both call
    `optimize_greedy` for their Stage-1 set cover.

    Taken together, entries 9 to 11 moved the three shipped whole-genome designs
    by almost nothing. Coverage changed by at most 0.17 percentage points, the
    delivered panels have Jaccard 0.993, 0.976 and 1.000 against their
    baselines, and the candidate pool size is identical on all three. The
    evenness rule bites on small targets, which is where the defect was
    measurable in the first place.

12. **CLI startup used to import scikit-learn** -- FIXED 2026-09-10 (audit
    finding E6). `cli_unified.py` imported `core.pipeline` at module scope only
    to make `StepPrerequisiteError` catchable, and `core/pipeline.py` imports
    `rf_preprocessing`, which imported sklearn at module scope. Every invocation
    paid it, `--help` and `show-presets` included, for a model retired from the
    default path on 2026-09-05.

    `StepPrerequisiteError` and `StepValidationResult` now live in
    `core/exceptions.py` (re-exported from `core/pipeline.py`, the same objects,
    so `except` clauses and the two importing tests are unaffected), and the
    sklearn alias fix runs inside `load_model_safely` instead of at import.

    `neoswga --help` now costs **about a third of what it did**. Quote the ratio
    rather than an absolute pair: measured twice hours apart the saving was
    about 3.3x both times, while the before figure itself moved from 1.21 s to
    1.40 s between sessions with no code change, purely with machine load.
    `python -X importtime -c "import neoswga.cli_unified"` reports no sklearn
    entry at all. `tests/test_cli_import_is_light.py` fails if either import
    comes back. What remains is not sklearn: it is pandas, reached through
    `cli/_common.py` -> `reaction_conditions` -> `thermodynamics` -> `utility`.
    That chain is pre-existing and is the obvious next target if CLI startup is
    worth more work.

13. **Five commands measured the host genome they were told to ignore** --
    FIXED 2026-09-10. Each read `bg_prefixes` from params, built a
    `PositionCache` over `fg_prefixes` alone, then handed the background
    prefixes to something that queries that cache by prefix.
    `PositionCache.get_positions` answered an unindexed prefix with an empty
    array, silently, so every background lookup read zero.

    The manifestation is `NetworkOptimizer._evaluate_primer_addition`, whose
    score is `fg_improvement / (1.0 + bg_added)`. With an fg-only cache
    `bg_added` is always 0.0, so every candidate scored as perfectly selective.
    This is the same symptom as Known Issues 5 and 6, reached by a third route:
    not a scan that found nothing and not an integer that saturated, but a cache
    asked for something it does not hold.

    It mattered most in `expand-primers`, which exists to add primers to an
    existing panel, so specificity is the property the user is asking it to
    preserve.

    `get_positions` now raises `MissingPositionsError` for a prefix the cache
    was not built over. That uses a separate `on_unindexed_prefix` knob, not the
    existing `on_missing`: a primer with no hits on an INDEXED prefix is a
    plausible measurement of zero and warns, while a prefix nobody indexed is a
    caller error and raises.

    `tests/test_expansion_counts_background.py` walks the AST for any function
    that forwards `bg_prefixes` while building a cache without them. That check
    found the fifth site after a manual review had settled on four.

14. **Reading the background is not acting on it.** A background-aware stage
    that does not choose the panel changes nothing useful. `expand-primers` was
    fixed on 2026-09-10 to build its `PositionCache` over
    `fg_prefixes + bg_prefixes`, and a real run afterwards still queried the
    host prefix zero times. Three further seams had to be closed before the data
    was read at all: `background_pruning` defaulted to False on the expansion
    path, `PrimerExpander.expand` silently substituted `hybrid` for every method
    it did not recognise including `background-aware`, and `_prune_background`
    would have removed primers from the very panel the user asked to extend.

    Even then it read the host without acting on it. Stage 1.5 pruning is not
    the stage that picks the panel; Stage 2 `_network_refine` is, and it ranked
    on amplification connectivity and unique coverage bins alone. Enabling
    pruning therefore only shrank the pool Stage 2 drew from, and on a
    40-candidate expansion over a 300 kb synthetic pair it moved delivered host
    binding the wrong way, 32 sites to 45. `_STAGE2_BACKGROUND_WEIGHT` adds the
    host as a third normalised axis in that stage, gated on `background_pruning`
    so `hybrid` panels are unchanged (verified identical on all three GC tiers
    at n=12/24/36).

    The general lesson: check which stage produces the delivered result before
    concluding that a measurement reaching the code means it reached the user. A
    query count answers "was it read", not "did it matter".

    **The two Stage 2s carry different things and neither carries both** --
    found 2026-09-19 while wiring Phase 6's deficit objective. The host term
    above lives in `_network_refine`. The objective a search can be steered by,
    `pool_objective`, is read only by `_swap_refine`. So a host-aware expansion
    cannot rank by recovered deficit, and a deficit-targeted one is not
    host-aware. `PrimerExpander._expand_hybrid` chooses between them on
    `background_pruning` and WARNS when target gaps are present but cannot
    steer selection, rather than narrowing the pool to the gaps and then
    ranking by something else. Switching expansion to `swap` wholesale was the
    first attempt and `tests/test_expansion_uses_the_background.py` caught it
    immediately: background-aware and hybrid returned the same panel, because
    the host term had been left behind. Combining them means putting the host
    term into the swap score as a weighted axis rather than its current
    lexicographic tie-break, which is unmeasured.

    `examples/plasmid_example` cannot demonstrate any of this. Six primers
    already cover its 5.4 kb target completely at 3 kb reach, so expansion adds
    nothing and `optimize` early-returns before Stage 1.5. Pin this behaviour on
    a target large enough that Stage 1 over-selects;
    `tests/test_expansion_uses_the_background.py` builds one at 300 kb with no
    external tool.

15. **The guard against the silent zero was itself silent** -- FIXED 2026-09-17
    (Phase 4 increment 3 of the 2026-09-16 pipeline audit).
    `CandidateProvider.ensure_positions` exists to refuse a candidate whose
    binding data is absent, so an unmeasured primer cannot be scored as though
    it bound nothing. It had two defects and each one alone made it useless.

    It returned quietly when no position cache was attached, and nothing in
    production attached one. So the single configuration it was written to
    catch was the configuration in which it did not run.

    And its predicate asked whether the candidate had a hit on ANY prefix:

    ```
    not any(len(cache.get_positions(prefix, sequence, "both"))
            for prefix in cache.fname_prefixes)
    ```

    A candidate with fifty foreground sites and no background entry at all
    therefore passed, which is unknown specificity reported as perfect
    specificity. A candidate indexed against a host it binds nowhere failed,
    though that zero is a measurement and a good one. The question is whether
    there is an ENTRY on EVERY prefix the design scores against, which is the
    distinction `_resolve_missing` already drew for the constructor's primer
    list and `PositionCache.has_entry` now exposes.

    `PositionCache.require_entries` holds the rule once, for both the inventory
    provider and a `--candidates` list, and `pool_planner._prepare_candidate_pool`
    is where `plan-pool` attaches the cache and runs the check, before any panel
    is evaluated. Tests:
    `tests/test_positions_arrive_on_demand.py` for the behaviour and
    `tests/test_the_frontier_is_vouched_for_before_it_is_scored.py` for the
    wiring, the second because the first would have passed throughout the years
    the check was inert.

    `load` and `release` arrive with it. The cache took a fixed primer list at
    construction, which was sufficient only while a design never looked past
    the `max_primer` shortlist. `load` also drops the memoized `both` key for
    the primers it admits: `get_positions` writes one for any primer it is
    asked about, including the empty one it returns for a primer the cache does
    not hold, so without that invalidation a candidate would keep answering
    with the zero it gave before its positions arrived.

16. **The stage that picks the panel is the least informed one** -- WIRED and
    measured 2026-09-19, and deliberately OFF by default (audit
    [pool_selection_audit_2026-09-18.md](../validation/pool_selection_audit_2026-09-18.md)).
    `optimize_greedy` takes an `objective`, and supplying it makes Stage 1
    select on occupancy-weighted coverage with a background tie-break instead
    of on unweighted coverage bins. Commit `59a4ee3` added it under the heading
    "The greedy now chooses on the quantity the design is judged on". **No
    production caller passes it.** `hybrid_optimizer.py:711`,
    `dominating_set_adapter.py:164` and `primer_expansion.py:652` all omit it;
    the only caller that supplies it is
    `tests/test_partial_panel_pruning.py:214`.

    It matters exactly where additives matter. Occupancy depends only on the
    primer, so an unweighted bin count misranks two candidates by the ratio of
    their occupancies, and across the pool the Tm gate admits that ratio is 1.8
    on phi29 at 30 C, 7.8 on equiphi29 at 42 C and 8.3 under DMSO 5% plus
    betaine 1 M. On phi29 occupancy is saturated and the unweighted count is
    nearly right, which is the same reason phi29 offers no discrimination.

    The ratchets cannot see this class.
    `tests/test_no_capability_is_unreachable.py` walks reach to FUNCTIONS and
    `optimize_greedy` is reachable; a PARAMETER no caller supplies is invisible
    to it. This is a fifth route into Known Issue 8's class, and the list there
    should be read as covering options and capabilities but not arguments.

    **Now reachable, as `stage1_objective_width` in params.json, defaulting to
    None.** An integer turns the objective on and bounds its cost: the cheap
    bin gain ranks every candidate and only that many leaders are scored. The
    bound is not optional -- one `compute_metrics` call costs 36 ms on the
    Wolbachia design, so a full scan is 14.4 minutes against 39 s for the whole
    run, which is the likeliest reason this stayed unwired.

    **It is off by default because measurement does not support switching it
    on.** At n=6/12/24 it improves the metric it now selects on (effective
    coverage +0.0073, +0.0139, +0.0583) and costs specificity every time
    (density -1.64, -6.50, -7.60; host sites 261 to 456 at n=24) for 3.5x to
    8.4x the runtime. With no width set the delivered panel is identical to
    before, verified at n=12 to the last digit. Same resolution as Known Issue
    11, for the same reason
    ([measurement](../validation/stage_one_objective_2026-09-19.md)).

    Two things the wiring clarified. `PoolObjective.coverage()` is coverage and
    nothing else, and `_objective_gain` is a coverage delta, so constraints
    reach Stage 1 only as a STOP rule and never as a selection criterion --
    "select on the quantity the design is judged on" changes what COVERAGE
    means, not whether specificity is weighed. And occupancy weighting favours
    primers whose Tm sits near the reaction temperature, which is the same
    property that makes them bind the host, so the unweighted bin count was
    accidentally the more specific rule. That is Known Issue 17's axis.

    Stage 1 uses `stage1_pool_objective`, NOT `pool_objective`. Writing it into
    the latter overwrote the objective `plan_pool` attaches for Stage 2 -- with
    None on every default run -- silently undoing the Phase 6 fix recorded in
    `attach_search_config`. Four tests caught it; keep the two names apart.

17. **The pool cannot discriminate, and the candidate filter is not the fix**
    -- measured 2026-09-19. Half of this entry's original diagnosis does not
    survive measurement, and the remedy it implied makes panels worse.

    **The floor is not padding the pool.** phi29's default floor does sit 10 C
    below its reaction temperature, but at k = 12 there is nothing down there:
    of 40,000 random 12-mers, 8 fall in the Tm 20-25 band (0.02%). Moving the
    floor changes essentially nothing.

    **Saturation is real and severe.** 86% of random 12-mers sit at or above
    0.998 occupancy at phi29 30 C, where a 4 C mismatch penalty leaves
    discrimination -- matched over single-mismatch occupancy -- at 1.01 or
    less. On the real Wolbachia shortlist, 65% of 2,000 candidates are above
    0.99 occupancy and mean discrimination is 1.11. Specificity in such a pool
    is a property of where sites fall, not of binding.

    **But gating on occupancy delivers a worse panel.** Measured at n=12 with
    the candidate list authoritative: capping occupancy at 0.95 raises the
    delivered panel's discrimination 1.098 to 1.400 and costs coverage 0.7334
    to 0.4397, selectivity density 25.62 to 6.60, with host sites RISING 149 to
    237. Discrimination lives in a tail too small to build a panel from --
    candidates above 2 are 1.6% of the space at k = 12 and 30 C.

    **The lever that works is the reaction.** Same 20,000 12-mers: phi29 30 C
    gives mean discrimination 1.065 with 1.6% above 2; DMSO 10% plus betaine
    1.5 M gives 1.374 and 11.3%; equiphi29 at 42 C gives 2.190 and 36.0%.
    Nothing about the candidates changes in any row. Saturation is a
    phi29-at-30-C problem, not a filtering problem.

    So what ships is a measurement, not a gate. `occupancy.discrimination_profile`
    computes the regime and `log_discrimination_profile` reports it at the end
    of `filter`, warning below `DISCRIMINATION_FLOOR` (1.5, between the two
    measured regimes) and naming the lever that works and the one that does
    not. No candidate is filtered and no delivered panel moves
    ([measurement](../validation/occupancy_and_discrimination_2026-09-19.md)).

    There is still no occupancy gate anywhere, and that is now a decision
    rather than an omission; `occupancy_ranking` remains the nearest thing and
    ranks on background load rather than on whether the candidate binds the
    target. Untested: whether a discrimination TERM in selection, as opposed to
    a gate on the pool, would help.

    A methodological note worth keeping. The first run of the gate experiment
    appeared to show density IMPROVING to 44.28. It did not:
    `open_source_or_list` prefers the inventory over the supplied list, so
    swapping `step3_df.csv` only set the frontier SIZE and the run searched the
    inventory as usual. The apparent gain was the smaller frontier. It was
    caught by checking that the delivered primers were actually in the capped
    pool -- none of them were.

    The consequence for the additive lever is measured in the audit. Occupancy
    and mismatch discrimination move in opposite directions along the Tm axis,
    so an additive improves every GC class at or above 6 of 12 and degrades
    every class below it, moving the best class up one step. The two routes to
    specificity conflict: compositional rarity favours GC-rich against an AT-rich
    host, thermodynamic discrimination favours AT-rich, and their correlation at
    k = 12 is about -0.89. An additive is the only lever that moves a candidate
    along the thermodynamic axis without changing its composition, which is why
    the best design measured in `docs/validation/additive_specificity.md` is an
    additive design at k = 12 rather than a longer-primer one.

    Do not conclude from the pool-size table that longer primers help. Above
    k = 15 an additive admits more candidates rather than fewer, and that regime
    is saturated: mean discrimination is 1.09 at k = 18 against 2.99 at k = 12,
    and occupancy spread across the admitted pool collapses to 1.0. A draft of
    the audit recommended k >= 15 before the discrimination column was measured.

18. **`max_gap` and `bg_coverage` are computed and read by nothing that
    selects** -- found 2026-09-18, partly acted on. Both reach
    `step4_improved_df_summary.json` and the reports. Neither appears in
    `normalized_score` or in any optimizer's scoring, and deliberately still
    does not.

    What changed the same day: both are now **constrainable** via
    `max_worst_hole` and `max_host_coverage` (see **Panel limits** above) and
    both are **reported** by the "What limits this panel" table, which also
    names them as having no reference. So a user can hold a panel to either,
    and neither has acquired a default -- no threshold derived from the reach
    separates the published wet-lab winners, so picking one would be the
    scoring change that evidence refuses.

    Also fixed on 2026-09-18: three of the five strand quantities
    `PositionCache.compute_strand_alternation_stats` returns were computed and
    discarded at the call site, and the loop stopped after the first foreground
    prefix so the host was never measured at all.
    `core/strand_metrics.py` collects all five for every foreground genome AND
    the background onto `PrimerSetMetrics.strand_stats`, keyed by prefix.
    `strand_alternation_gap_max` is the one that mattered: exponential
    amplification needs two sites in convergent orientation within the
    polymerase's reach, so the widest gap between opposite-strand sites is the
    closest quantity here to the mechanism, and on the host it is what swga 2.0
    approximates with `within_mean_gap_ratio` and fits against measured
    sequencing breadth. It reaches the report as `convergent_gap` and
    `host_convergent_gap`.

    A prefix the cache cannot answer for is now ABSENT from `strand_stats`
    rather than zero, and the two headline scalars are `None` rather than 0.0.
    They were initialised to 0.0 and left there, so a zero meant either
    "measured zero" or "never asked".

    **The source was fixed on 2026-09-19 and the consumers audited.** A
    one-site panel does NOT genuinely score 0.0 for alternation, which an
    earlier note here got wrong: alternation is the fraction of ADJACENT site
    pairs on opposite strands, so below two sites there is no pair and the
    fraction is 0/0. `compute_strand_alternation_stats` now returns None there,
    while two same-strand sites still return a measured 0.0.
    `strand_coverage_ratio` is min/max over the two strand counts and needs
    only one site, so it is None only with no sites at all; a lone forward site
    really is maximally unbalanced. The gap figures keep the genome length,
    which encodes "no convergent pair anywhere" rather than a missing
    measurement, and `worst_convergent_gap` reads it that way.

    Every consumer was checked and none needed changing: `panel_regime._as_float`
    and `report/metrics._safe_float` preserve None, `headline_strand_scalars`
    already returned `(None, None)`, and the technical report skips a row whose
    value is None. Pinned by
    `tests/test_strand_scores_say_when_they_are_unmeasurable.py`.

    They are still NOT constrainable, but the reason is now different: no
    threshold has a reference, which is why `max_worst_hole` and
    `max_host_coverage` ship unset rather than defaulted. Adding a strand limit
    is a decision about evidence, not a blocked repair.

    `bg_coverage` is the only computed quantity that sees background site
    POSITION. `selectivity_density` and `total_bg_sites` are additive in
    per-primer counts -- `occupancy.weighted_site_load` sums `count * theta` per
    mismatch class and no position enters -- so two backgrounds with identical
    per-primer counts score identically whether their sites are clustered or
    dispersed. That distinction is most of off-target amplification, since SWGA
    needs two convergent sites within the polymerase's reach.

    swga 1.0 made both criteria hard in 2017, and they are the only two things
    that constrain its selection: the clique search runs `--unweighted --all`
    and stores every clique passing the `max_fg_bind_dist` gap cut, so its score
    expression ranks the output and never steers the search. The background side
    is a pruning budget inside the recursion -- each vertex weight is the
    primer's raw background site count (`weight = primer.bg_freq` in
    `graph.py`), and a partial clique is pruned once the summed weight exceeds
    `bg_length / min_bg_bind_dist`. It is NOT a ranking by mean background
    binding distance; that quantity is what the search emits, as
    `bg_len / graph_subgraph_weight`, and an earlier draft of this entry
    conflated the two.

    **The field does not agree that gap statistics belong in the objective.**
    swga 2.0 fits both as `on_gap_gini` and `off_gap_gini`, where `off_gap_gini`
    carries the second largest recorded weight in the only set-level model
    fitted against measured sequencing breadth -- a value nobody has been able
    to verify from a source that opens, since it sits in a CAPTCHA-gated table.
    COATswga (2025) computes no Gini and no gap statistic at all, on the stated
    ground that a per-primer Gini cannot speak for a whole set, which is an
    argument against `max_gini` as much as for the interval-union objective this
    project already uses. So treat background evenness as a measurement to make
    before it is a term to add. No background amplification network is built
    here, in contrast to the foreground network the hybrid and network methods
    build at about 70 kb.

    Worth knowing about the ancestry: swga 2.0 is this project's direct
    ancestor, and the Known Issue 8 class is partly inherited. In its shipped
    master `filter.filter_extra` implements the GC, homopolymer, GC-clamp and
    self-dimer rules the paper describes and nothing calls it, and it would
    raise if called, reading a `default_max_self_dimer_bp` that `parameter.py`
    never assigns. Its step 2 also computes `ratio = bg_count / fg_count`, where
    lower is more specific, then keeps `sort_values(by=["ratio"],
    ascending=False)[:max_primer]`, retaining the LEAST specific survivors.
    NeoSWGA sorts that ascending. All three were reported from that repository's
    source and were not re-verified here.

19. **Four commands meant "every alternative set" where they should have meant
    one** -- FIXED 2026-09-21. `step4_improved_df.csv` holds up to `max_sets`
    (default 5) primer sets, one per `set_index`, and they are ALTERNATIVES:
    each is found by excluding the primers already chosen and selecting again.
    Set 0 is the one the summary describes. `export`, `interpret`, `report` and
    `simulate` all read every row and treated the union as one panel.

    The costly one is `export`, whose output the tool calls "Primers ready for
    ordering!". Measured on the bundled plasmid example at 300 bp reach, where
    the pool is large enough for alternatives to be found: set 0 is 8 oligos
    with no pair above the configured `max_dimer_bp` of 3, while the exported
    FASTA was 18 oligos with five pairs above it, the worst a 9 bp
    complementary run. All five join oligos from DIFFERENT sets, so no screen
    had ever compared them and by construction none could. Nothing in the file
    marked a set boundary: records run SWGA_001 upward straight through.
    `interpret` reported 18 primers, a count matching no orderable set.

    This is the 11 bp delivered heterodimer of the `max_dimer_bp` entry reached
    by a second route, with selection behaving correctly throughout. No saved
    run in this repository exhibits it -- every one holds set 0 alone -- which
    is why it survived. It needs only a pool big enough for a second set.

    `core/delivered_set.py` holds the rule once. Default set 0; `--set N` on
    `export` and `interpret`; a requested set the file does not hold raises
    `ReferenceDataError` rather than returning an empty panel; a file with no
    `set_index` column is older output holding one set and every row is
    returned. Tests: `tests/test_commands_read_one_primer_set.py`.

20. **Mixed oligo lengths work, and on the one pair measured they bought
    nothing.** A design may mix lengths: the schema admits k of 4 to 30, the
    scan writes one index per length, `PositionCache.load` groups by length and
    the dimer screen codes t-mers. No delivered panel here had ever mixed,
    so this was an argument until
    `tests/integration/test_variable_oligo_length.py` took a k 7-11 design
    through all four steps and out to an ordering file.

    Measured on Prevotella against human chr21 at equiphi29 42 C, a mixed
    k 10-12 design lands BETWEEN the single-length designs on both axes at
    panel sizes 8 and 20. It is beaten on coverage by the all-11-mer panel and
    on specificity by the all-12-mer panel, which binds the host once against
    the mixed panel's seventeen at n=8. **Choosing k matters far more than
    choosing whether to mix.**

    The occupancy spread across lengths is governed by the Tm WINDOW, not by
    the length range, which a control established after a first draft concluded
    otherwise. A 30 C window gives a 23.5 C median Tm spread across lengths and
    occupancy from 0.042 to 0.994; a 12 C window gives 2.0 C and 0.580 to
    0.794. `core/length_occupancy.py` reports this per length from `filter` and
    `optimize` and never enforces it, the resolution Known Issue 17 reached for
    the same quantity. Silent on a single-length design
    ([measurement](../validation/variable_oligo_length_2026-09-21.md)).

21. **The HDF5 lock failure is two runs sharing a data directory, not mixed
    length** -- established 2026-09-21. `tests/validation/genomes/f_mixed.log`
    records `BlockingIOError: [Errno 35] unable to lock file` from a
    mixed-length run, and mixed length was blamed because that is what the
    config changed.

    Four measurements settle it. **The errno is cross-process by
    construction**: a foreign process holding the file, even read-only,
    produces exactly the recorded message, while a second handle inside ONE
    process produces `OSError: ... file is already open for read-only` with no
    errno 35, so a leaked handle cannot be the cause of this message. **A
    reader is enough to stop a writer**, so a command that only reads the
    index can stop a `filter`. How wide that window is depends on the cache:
    the default `PositionCache` opens each file in a `with` block and closes
    it promptly, making the collision a race, while `StreamingPositionCache`
    holds handles until `close()` and is selected only when the in-memory
    cache is DISABLED. An earlier version of this entry said `optimize` holds
    them for a whole run, which is true only of that non-default path. **Concurrent pairs failed 10 of 11
    across two trials and single-process runs 0 of 11**, the latter including
    multi-k runs at realistic scale; a mixed k 10-12 `filter` in a clean
    directory takes 146 s and exits 0. **The recorded directory shows the
    collision directly**: its `run_manifest.json` has `score` on that config
    completing 2.2 s before the failing filter's last log write and `optimize`
    8.7 s after, with no `filter` entry at all because the manifest is written
    on completion, and a `prevotella_13mer_positions.h5` of 800 bytes holding
    zero datasets sits there with the same mtime for a k that run never
    requested.

    The recorded traceback differs from the reproduction only in `h5f.open`
    against `h5f.create`, which is whether the target file already existed.
    The frames match the sequential Aho-Corasick branch, so the per-k
    multiprocessing fallback -- the first guess -- did not run.

    **Which process held the lock is NOT determined, and probably cannot be.**
    At least three runs were live in that window. Because the default cache
    holds each file only briefly the collision is a race, so the blocker is
    whichever process had that one file open at that one instant, and no
    durable fact about "the holder" exists for the artifacts to have recorded.
    Read it as "not determined" rather than "not determined yet": better
    evidence would not settle it. The cause is settled and the mechanism is
    not fully reconstructed, which are different claims. Nothing depends on
    the answer, since the remedy is the same either way.

    `core/concurrent_runs.py` translates it at the step boundary: steps 2, 3
    and 4 all write HDF5 and each consults it. No lock is taken and no retry is
    attempted -- a retry loop would hide a genuine second process, and two
    designs writing one directory have a provenance problem that outlasts the
    lock. The message steers AWAY from `HDF5_USE_FILE_LOCKING=FALSE`, which is
    the usual first hit for this error and risks a corrupt index. The predicate
    matches on errno AND message, because EAGAIN alone is raised by unrelated
    things. Tests: `tests/test_two_runs_sharing_a_directory_say_so.py`.

    **This entry took three corrections and they were all the same shape**,
    which is worth keeping because it is the silent-zero family in a register
    this file does not otherwise cover: not a value, but a SENTENCE that reads
    as more definite than its evidence. "Another writer" when a reader
    suffices. "`optimize` holds these for its whole run" when only the
    non-default cache does. "Nothing in the artifacts separates the three
    processes" when a race means there is nothing to separate. None was wrong
    about the cause; each was wrong about the size of the claim, and the first
    two shipped to users in an error message. A diagnostic written from a
    finding is a claim about someone's machine, so name the mechanism that
    generalises rather than the command that happened to be involved, and say
    "not determined" only when better evidence would in fact settle it.

22. **The Bloom path could not be built at the scale it exists for** -- FIXED
    2026-09-24. Seven defects, and the shape they share is the one this file
    records throughout: every one of them passed every small-genome test.

    **Capacity bounded the wrong quantity.** pybloom allocates its bit array
    upfront from `capacity` and its `add` counts only items the filter did not
    already hold, so capacity bounds DISTINCT k-mers. `genome_size * 10`
    bounded INSERTIONS -- the multiplier was chosen for the seven k-mer lengths
    each position contributes. Measured at 9.59 bits per item:

    | genome | asked for | allocation | distinct k-mers |
    |---|---|---|---|
    | plasmid 5.4 kb | 53,860 | 0.06 MB | 36,361 |
    | E. coli 4.64 Mb | 46,416,520 | 55.6 MB | 10,232,681 |
    | Drosophila 144 Mb | 1,440,000,000 | 1.73 GB | 22,368,256 |
    | hg38 3.3 Gb | 33,000,000,000 | **39.56 GB** | 22,368,256 |

    The last two agree because each term saturates at `4**k`. hg38 is the
    documented reason the module exists and it is the row that cannot be
    allocated. `distinct_kmer_capacity` holds the bound and all four call sites
    ask it.

    **A saved filter reloaded only at one geometry.** `make_hashfuncs` selects
    the hash from `num_slices` and `num_bits`, so the constructor in the pickle
    moves with capacity and error rate. At error rate 0.01 the bands are
    capacity below about 3,400 (xxh3_128), up to about 224 million (sha256),
    above that (sha512); other error rates reach sha384 and sha1. The
    safe-pickle allowlist named sha256 alone, which is what one observed filter
    carried. `save()` succeeded and `load()` raised for a small background, and
    for any long-oligo host design. The ratchet asserts the RULE -- every
    constructor `make_hashfuncs` can select must be listed -- because reaching
    the sha512 arm behaviourally costs a 268 MB allocation.

    **A length the filter never indexed read as absent.** `contains` answers
    False for a k-mer of a length nobody inserted, the count is then zero, and
    zero clears any frequency gate, so a design at k 13-18 screened against a
    phi29-range filter passed its whole pool. Both artifacts now record
    `min_k`/`max_k` and `get_bg_rates_via_bloom` refuses outside them. A filter
    with no recorded range predates the field and warns rather than refusing,
    the rule `digest_algorithm` established. `BackgroundFilter.build_from_genome`
    could not have built a correct filter anyway: it took `add_genome`'s 6-12
    defaults regardless of configuration.

    **`use_bloom_filter` without a path screened nothing.** It took neither
    branch and fell through to exact counting over `bg_prefixes`, which that
    same flag leaves empty, so every background count was absent and an absent
    count passes. Now `InvalidDesignRequest`, raised before any counting.

    **One filename, two quantities.** `bg_sampled.pkl` holds sampled positions
    at rate 100 from the FASTA route and exact jellyfish counts at rate 1 from
    `--from-kmers`. `_warn_if_sample_too_sparse` reasons about sampling and was
    skipped on the second only because that route left `genome_size` at 0.
    `source` now records which quantity an index holds.

    **The library's auto-built filter was wrong in both directions.**
    `genome_library.add_genome` took the 3e9 default capacity -- 3.6 GB, far
    more than a k 6-12 filter needs and far less than the 14.6 billion distinct
    k-mers a k 6-18 filter holds, which is the range that path computes. The
    `except Exception` turned pybloom's IndexError into "Bloom filter build
    failed" with no filter registered. Capacity is now required, which is what
    stops a third such call site appearing.

    **The companion index undoes the memory argument.** Measured at 125 bytes
    per entry (`scripts/benchmarking/sampled_index_rss.py`, ru_maxrss, one size
    per process): hg38 at rate 100 over k 6-12 projects to 22.4 million entries
    and about 2.8 GB, against 26.8 MB for the filter beside it. So the
    structure the filter was chosen to avoid reappears at about a hundred times
    its size. `warn_if_sampled_index_is_large` says so and names `--from-kmers`.

    **What the FASTA route costs, and why `--from-kmers` is the answer.** The
    scan now slides inside maximal ACGT runs instead of revalidating every
    position; measured 1.21x to 1.33x on E. coli at k 6, 12, 18 and on
    Drosophila at k 12, because the dominant cost is pybloom's insert, which
    neither version changes. At 1.65 us per position, hg38 over k 6-12 is 23.1
    billion inserts, about 10.6 hours EXTRAPOLATED, against about 37 seconds
    for the 22.4 million unique k-mers `--from-kmers` reads. That route is
    therefore the one to use for a host background, and it is now validated:
    it previously indexed the first whitespace-delimited field of every line
    with no base or length check, while `add_genome` skipped ambiguous k-mers,
    so the two artifacts one command writes already disagreed.

    **Still not measured.** No Bloom filter has ever been built against a host
    genome in this repository, so every hg38 figure above is arithmetic from a
    measured per-item or per-entry constant, not a build. The `--from-kmers`
    route is also the only one whose artifacts `filter` can use, since
    `genome_library` writes a filter with no sampled index beside it and
    `get_bg_rates_via_bloom` requires one.

23. **Python 3.13 only, and the code modernised to match** -- 2026-09-25.
    `requires-python` is `>=3.13`; black, ruff and mypy all target it and both
    workflow matrices are one interpreter on two operating systems. There were
    no `sys.version_info` guards anywhere in the package, so nothing had to be
    unwound.

    2,294 modernisation sites were rewritten across the package and tests,
    almost all `typing.List` to `list` and `Optional[X]` to `X | None`. Ruff
    fixed 1,993; the orphaned `typing` imports took two passes and the last
    five were done by hand. **Ruff reports zero UP findings now**, so those
    rules could be made blocking if anyone wants them to be.

    The import removal is the part worth knowing about. A blanket unused-import
    fix would have deleted deliberate re-exports, so it was scoped to `typing`
    names only. The first pass then kept `List` imported wherever the word
    appeared anywhere in the file, including docstring prose like
    "candidates: List of candidate primers"; the second asks whether the name
    is used as CODE via the syntax tree, while still keeping forward references
    inside string annotations.

    **One test broke, and it was asserting a spelling rather than a property.**
    `test_multi_genome_result_allows_none_metrics` required the literal string
    "Optional" in a dataclass field's annotation, so it failed on `float | None`
    while its subject had not changed. It now checks that NoneType is in
    `typing.get_args`, which holds under either spelling and still fails if a
    field is made non-optional.

    **`mip`'s default solver kills Python 3.13.** Constructing a CBC model
    terminates the interpreter with SIGKILL -- no exception, no traceback, no
    stderr -- measured on macOS arm64 with mip 2.0.0 and cbcbox 2.935. The same
    versions on 3.11 solve the same models in 0.34 s. `core/ilp_solver.py`
    prefers HiGHS, which works on both, and **refuses rather than falling back**
    as of 2026-09-27: a warning that precedes a SIGKILL is never read in
    context, because the process dies with no traceback and the user has no
    reason to connect a killed command to a log line. Refusing costs one
    `pip install highsbox` and names it. The fatality measurement is ONE
    platform, so `NEOSWGA_ALLOW_CBC=1` asks for CBC explicitly, which is a
    request rather than a substitution nobody made.
    The fallback is deliberately unprobed: a probe would construct a CBC model,
    which is the operation that kills the process, so nothing in process can
    tell a working CBC from a fatal one.

    Two construction sites needed it, and the second is easy to miss:
    `dominating_set_optimizer.py` and
    `scripts/benchmarking/max_coverage_bound.py`, which a test loads and
    executes. Wiring only the library left 16 of 17 tests passing and the
    seventeenth still killing the run.

    **Those 17 tests did not run in CI until 2026-09-30.** `mip` lives only in
    the `improved` and `all` extras and CI installed `.[dev]`, so they skipped
    there, and nightly installs `improved` but runs only the scale-marked
    tests. It was not only the ILP path: 19 modules `importorskip` an optional
    dependency at module level, so 185 tests were never collected in CI, every
    Bloom and BAM test among them, while the job stayed green. The test job
    now installs `.[dev,improved,bam,viz,interactive]` and passes `-rs`, so a
    skip is in the log with its reason. Reading those reasons found 17 more
    tests skipped as "needs jellyfish" on a runner that had it: a
    `skipif(not plasmid_example_ready())` is evaluated at collection, before
    the session fixture primes the example. Request the
    `primed_plasmid_example` fixture instead; a ratchet in
    `tests/test_plasmid_example_dependency_is_visible.py` refuses an
    import-time call. CI went from 6,416 passed and 111 skipped to 6,643
    passed and 33 skipped, and every remaining skip is opt-in, by design, or
    needs a reference genome the repository does not hold.

    **bioconda has no Python 3.13 build of `kmer-jellyfish`.** 2.3.1 exists for
    3.9 to 3.12 only, so `conda create ... python=3.13 kmer-jellyfish` silently
    resolves to 1.1.12, whose CLI this project refuses. CI is unaffected
    because it installs Jellyfish through apt and brew rather than conda.
    Recorded in docs/guides/TROUBLESHOOTING.md, since it costs an hour to
    diagnose from the symptom.

    **PyYAML was reaching CI only as a transitive dependency of pre-commit**,
    while `tests/test_workflows_invoke_real_commands.py` `importorskip`s it.
    A ratchet that can silently vanish is the defect class that ratchet exists
    to catch, so it is now declared in the `dev` extra.

    Measured after all of it: 6,424 passed and 27 skipped on BOTH 3.13 and
    3.11, with identical collection of 6,444 tests on each, so nothing is
    quietly missing from either.

24. **KMC3 is preferred, and tables are read as databases** -- 2026-09-25.
    `kmer_counter` in params.json takes "kmc" or "jellyfish" and, when set,
    REQUIRES that counter. Unset, KMC3 is used when installed and jellyfish
    otherwise. `core/kmer_backend.py` owns the invocation;
    `core/kmer_tables.py` is the one place anything asks about a table, and
    callers name a prefix and a k rather than a file.

    **Unset falls back rather than failing**, which departs from a strict
    "KMC is the default". A hard requirement would have broken CI, which
    installs only jellyfish, and every jellyfish-only installation. The
    fallback is acceptable to do unasked only because the choice does not
    change any result -- and that was FALSE until a fix described below.

    **The first version of this change was inert.** `kmer_counter` was
    declared, validated, defaulted and documented, and `count-kmers` still ran
    jellyfish at every call site: Known Issue 8's class, on the branch that
    introduced the key. Tests covered the backend by calling it directly, so
    none could see it. `tests/test_count_kmers_uses_the_configured_counter.py`
    now walks from the two counting entry points, `run_jellyfish` and
    `MultiGenomeKmerCounter`, and five of its tests fail against the old
    behaviour.

    **The counter choice used to change results.** Step 2 sorted by `ratio`
    then `fg_count`, and ties kept their INPUT order -- the order the counter
    emitted k-mers in, hash order for jellyfish and sorted order for KMC. That
    order leads step 3 into an order-sensitive optimizer, and through
    `[:max_primer]` it decided which tied candidates survived the shortlist.
    The primer sequence is now the final sort key, and step 2's written index
    is reset so it records rank rather than emission order. Verified by running
    the plasmid example end to end under each counter from clean copies:
    `step2_df.csv`, `step3_df.csv` and `step4_improved_df.csv` are
    byte-identical, and the KMC run writes no text table. Found by comparing
    files after the claim "changes speed, never a result" had already been
    written into the code, which is why the comparison is worth keeping as a
    habit. Tie order now differs from the old jellyfish-only order, so an
    existing design with ties at a boundary can see a different shortlist; no
    test-pinned output moved.

    **What this buys is the QUERY path, not the counter.** `filter` asks one
    question of a table: the counts of a known candidate list. It answered by
    streaming every line into Python and testing set membership, which is a
    scan answering a set question. Measured on *Drosophila* at k=18 with 2,000
    candidates:

    | approach | time | intermediate |
    |---|---|---|
    | dump to text, then Python scan | 17.9 s | 2.3 GB |
    | stream the dump, then Python scan | 15.6 s | none |
    | `kmc_tools simple ... intersect` | **2.5 s** | 2,000 lines |

    **`-ocleft` is load-bearing and silent when wrong.** It keeps the counters
    of the FIRST database. Without it the output carries the candidate
    database's counters, which are all 1, so every background count reads 1 --
    a wrong answer that looks entirely plausible. Verified by removing it,
    which fails the test asserting the database and text paths agree.

    **What it does NOT buy is memory, and on hg38 the gap is 18.6x.** Measured
    on the 3.1 GB human genome at k=12 through the shipping backend: KMC
    counts in 12.2 s against jellyfish's 89.7 s, and peaks at 1,950 MB against
    105 MB. So KMC is 7.2x faster and uses 18.6x the memory, which is the
    trade in one line. `kmc -m1` also refuses outright; the floor is 2 GB,
    where jellyfish counted wMel in 18 MB.

    The lookup on that hg38 database beats scanning its 8.4 million line text
    table by 4.8x at 2,000 candidates and 2.1x at 500,000, the advantage
    narrowing because building the candidate database is itself work. Every
    absolute figure there is under two seconds, so at k=12 the win is real and
    the stakes are modest; k=18 on *Drosophila* is where the 7x lives
    ([measurement](../validation/kmer_counter_comparison_2026-09-25.md)).

    **A host with no table is now scanned rather than refused.**
    `core/query_scan.py` counts a KNOWN k-mer set directly in a reference,
    which is the question `filter` asks of a background, and it reaches that
    module only where a prefix has no table. A counted prefix takes the table
    path unchanged. Counting is otherwise FASTER and stays the default: on
    Drosophila at k=12, jellyfish counts in 0.9 s and answers a 2,000-candidate
    batch in 0.2 s, against 4.0 s to scan. What the scan avoids is the table --
    33.7 MB at k=12 and 818 MB at k=18 on that same 144 Mb reference -- so it
    earns its place only where the table is the problem.

    Agreement is checked three ways, because a fallback that disagrees with
    the path it replaces is worse than none: against a brute-force oracle
    written in the test file, against both counters, and against the
    production `counts_for` on real references. One difference is known and
    pinned. A canonical table stores one spelling of each reverse-complement
    pair, so `counts_for` answers 0 for the other spelling while the scan
    answers the pair's count -- 13 of 1,272 queries on Drosophila at k=18, all
    non-canonical, with all 608 canonical queries agreeing exactly. Every
    caller here reads its k-mers from a canonical table, so the pipeline never
    asks the other spelling; the 0 is still not a measurement, and a test
    holds both behaviours so a deliberate fix would show as a change
    ([measurement](../validation/query_scan_2026-09-25.md)).

    Memory is bounded by the largest RECORD and the chunk, NOT by the query
    set: 332 MB on Drosophila, 46.6 MB on wMel. The module's first draft
    claimed otherwise and the first measurement refuted it. The chunk size was
    measured rather than chosen -- the prototype's 8,000,000 positions was the
    worst value on both axes, 6.6 s and 1,194 MB against 3.9 s and 351 MB at
    the shipped 125,000.

    **CI installs both counters as of 2026-09-25.** It had only jellyfish, so
    the preferred path never ran there, and four tests had quietly encoded
    that: two asserted jellyfish's `*mer_all.txt`, which KMC does not write,
    and two asserted an absent jellyfish is fatal, which stopped being true
    when KMC became preferred. All four passed only because no runner had KMC.
    Five more skipped on jellyfish alone, so a KMC-only machine skipped most
    of the counting tests. **Installing it found three more fixtures asking
    for a table by jellyfish's filename**, which is the same defect one layer
    down: with KMC installed they found nothing, so `plasmid_example_ready()`
    reported the example unprepared and 48 tests skipped as unavailable while
    the directory was ready. Ask `kmer_tables.table_exists` or
    `discover_prefixes`, never a glob -- and note a KMC database is
    `{prefix}_{k}mer.kmc_pre`, so cutting at the last underscore cuts inside
    `kmc_pre`, which is how the first attempt at that guard silently answered
    "not prepared" for a prepared directory.

    Unskipping those 48 exposed a real defect they had been hiding.
    `_filter_blacklist_penalty` called `counts_for` once per PRIMER, and a
    lookup against a KMC database builds a database of the query set,
    intersects and dumps it -- three processes per call. A four-prefix design
    spent minutes in that gate; batched by k it is seconds. A count lookup has
    a fixed cost per CALL, so ask once per group, which is why
    `get_rates_for_one_species` groups by k before asking. It now lives in
    `core/blacklist_penalty.py`, extracted because `pipeline.py` had reached
    its size budget; `pipeline` re-exports the name so its five importers are
    unaffected.

    KMC comes from the upstream release tarball pinned
    to 3.2.4, not from a package manager: it is in neither apt nor the default
    brew taps, and the conda route would put a second Python on PATH.

    **Above k=12 a host genome cannot have a text table at all.** The distinct
    count stops being bounded by the k-mer space and becomes bounded by the
    genome, so hg38 at k=16 or k=18 would dump about 78 to 84 GB of text. That
    is the sharpest argument for reading databases rather than dumps, and it
    is also why those k could not be benchmarked here.

    **Two figures this file recorded did not reproduce.** It says hg38 at k=12
    costs "about 7 minutes and a 138 MB table (8,368,418 canonical 12-mers)".
    The table size matches at 138.4 MB, which is what confirms it is the same
    quantity; the time was 91 s here, which one machine against another
    explains; and the distinct count was 8,368,476, which does not have an
    explanation. Both counters agree with each other on that number, so it is
    not a tool artifact. A different hg38 assembly is the likeliest cause and
    is unestablished.

    **KMC's defaults compute a different quantity and both differences are
    failures this file already carries.** `-ci2` excludes k-mers occurring
    once, which at k=18 is 96.9% of wMel's and 95.6% of *Drosophila*'s: a
    background counted that way reports a host as almost k-mer-free and every
    candidate as specific, which is Known Issues 5, 6, 13 and 15's shape.
    `-cs255` saturates the counter at 255, which is Known Issue 7 exactly. The
    backend passes `-ci1` and `-cs1000000`, and a test asserts neither default
    returns.

    **`py_kmc_api` is not usable here.** It exists, and the bioconda package
    even ships `py_kmc_api.so`, but that build is compiled for Python 3.10: it
    fails on 3.11 with an explicit version mismatch and on 3.13+ because
    `__PyThreadState_UncheckedGet` was removed from CPython. Upstream's README
    also warns the wrapper is "much slower than native C++ API". So Python
    does not read `.kmc_suf` directly; everything goes through the C++ tools.

    **A table is any of three forms**, and `table_exists` is the one predicate
    that says so. Half a KMC database is not a table: it writes two files, and
    one alone is an interrupted run. Two things a find-and-replace would have
    broken: the genome library symlinked the text table by name, so a database
    source linked nothing and the pipeline recounted a genome it already had;
    and auto-discovery globbed `*_6mer_all.txt`, so a directory of databases
    looked empty.

    **Every existing data directory still works.** They hold text tables and
    no database, and that path is preserved. It is also the only path CI
    exercises, since CI installs no KMC -- which is why the text-fallback
    tests are written to run without a counter rather than skipping with the
    rest.
