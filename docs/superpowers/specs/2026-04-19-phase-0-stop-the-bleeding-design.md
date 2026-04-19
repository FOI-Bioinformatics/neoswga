# Phase 0: Stop the Bleeding - Design

**Status:** Draft for review
**Date:** 2026-04-19
**Parent roadmap:** `2026-04-19-production-readiness-roadmap.md`
**Target version:** neoswga 3.7.0

## Purpose

Close the five correctness bugs identified in the production-readiness audit
that are capable of producing incorrect results today. No new features. No
reproducibility work. No model work. Everything downstream in the roadmap rests
on this phase.

## Out of scope

- RF model polymerase-awareness (Phase 3).
- Run manifest and seed control (Phase 1).
- Integration tests in CI, Jellyfish install in CI (Phase 2).
- BAM ingestion and anything feedback-loop related (Phase 4).
- Windows CI, conda-forge, deprecation decorator (Phase 5).

## Work streams

Four independent work streams. They can be sequenced or parallelized. Each
ends with tests committed.

### 0.1 Wire AdaptiveGCFilter into the pipeline

**Problem.** `neoswga/core/adaptive_filters.py` defines `AdaptiveGCFilter` with
a tolerance band around the genome's own GC content. The class is never
instantiated in `filter.py` or `pipeline.py`. The actual GC filter uses fixed
37.5-62.5 percent thresholds, so any genome with mean GC outside that band
(Francisella 33 percent, Plasmodium 20 percent, Streptomyces 70 percent)
loses all candidate primers at the filter stage.

**Changes.**

- `neoswga/core/parameter.py`: add to `PipelineParameters`:
  - `use_adaptive_gc: bool = True`
  - `gc_tolerance: float = 0.15`
- `neoswga/core/filter.py`: in the GC filter pass, when `use_adaptive_gc` is
  true, instantiate `AdaptiveGCFilter(genome_gc, tolerance=gc_tolerance)` and
  use its min/max thresholds. Fall back to fixed thresholds when false.
- `neoswga/core/pipeline.py`: compute foreground genome GC in step 2 from
  the target FASTA(s) via `genome_io.GenomeLoader`. For multi-target, use
  the length-weighted mean of target GCs. Cache the value for reuse across
  stages. No new `fg_gc` parameter field - GC is always derived.
- `neoswga/core/param_validator.py`: after step 2, if `use_adaptive_gc=False`
  and the derived foreground GC is outside 0.35-0.65, emit a warning
  pointing at `use_adaptive_gc` as the fix.

**Tests.**

- Unit: `tests/test_adaptive_filter_wiring.py` - given `fg_gc=0.33, tolerance=0.15`,
  filter accepts a primer at 40 percent GC and rejects one at 60 percent.
- Integration: add a low-GC target (e.g. truncated synthetic 33 percent GC
  FASTA) under `tests/fixtures/low_gc/`. Run the filter stage and assert
  `len(step2_df) > 0`. Running with `use_adaptive_gc=False` should return
  approximately zero - this is the bug being fixed; the contrast is the
  regression guard.
- Regression: existing `tests/test_plasmid_e2e.py` still passes.

**Acceptance.** Low-GC integration test passes. Plasmid E2E unchanged.

### 0.2 Complete blacklist integration

**Problem.** Three untracked files exist:
`neoswga/core/genome_library.py` (new, well-formed),
`tests/test_blacklist_filter.py`, `tests/test_blacklist_cli.py`,
`tests/test_genome_library.py`, `tests/test_gc_adaptive_pipeline.py`.
The tests expect fields that are partially present, partially missing:

- `PipelineParameters.bl_genomes` and `.bl_prefixes` exist.
- `PipelineParameters.bl_seq_lengths`, `.bl_penalty`, `.max_bl_freq` are
  referenced by pipeline code and tests but not declared in the dataclass.
- `pipeline.py:659-660` reads `bl_seq_lengths` via `getattr(..., [])`,
  silently passing an empty list downstream, which risks division-by-zero
  in frequency calculation and masks the misconfiguration.
- `MultiGenomeFilter` does not receive blacklist k-mer counts - there are
  two parallel code paths, pipeline-legacy and multi-genome, and only the
  legacy one sees the blacklist.

**Changes.**

- `neoswga/core/parameter.py`: add to `PipelineParameters`:
  - `bl_seq_lengths: list[int] = field(default_factory=list)`
  - `bl_penalty: float = 5.0`
  - `max_bl_freq: float = 0.0`
- `neoswga/core/pipeline.py`: in the step that reads blacklist prefixes, if
  `bl_seq_lengths` is empty and `bl_genomes` is non-empty, compute lengths
  from the FASTAs via `genome_io.GenomeLoader`. Raise
  `ConfigurationError` if `bl_genomes` is non-empty but prefixes cannot be
  resolved.
- `neoswga/core/pipeline.py:_filter_blacklist_penalty`: guard against
  `sum(bl_seq_lengths) == 0`, raise a clear error rather than returning
  `nan`.
- `neoswga/core/multi_genome_filter.py`: accept blacklist counts via
  `load_genome_counts(..., blacklist_counts=...)`. Integrate blacklist
  scoring into `score_primer` with `max_bl_freq` enforcement.
- `neoswga/core/wizard.py:598-605`: verify the generated config round-trips
  through `PipelineParameters(**config)` without AttributeError. Write a
  round-trip test.
- Commit the four untracked test files after they pass.

**Tests.**

- Unit: the four untracked tests (`test_blacklist_filter.py`,
  `test_blacklist_cli.py`, `test_genome_library.py`,
  `test_gc_adaptive_pipeline.py`) pass as-is or with minimal adjustments.
- Unit: `tests/test_parameter.py` - new test asserts all five `bl_*` fields
  exist with documented defaults.
- Unit: `PipelineParameters(**wizard_config).bl_penalty == 5.0` round-trip.
- Integration: add a blacklist genome scenario under
  `tests/fixtures/with_blacklist/`. Compare `step4_improved_df` with and
  without the blacklist and assert fewer blacklist-hitting primers in the
  final set with it.

**Acceptance.** Four committed tests pass. Blacklist integration test
demonstrates fewer blacklist hits when blacklist is supplied.

### 0.3 Collapse the Tm paths to one canonical implementation

**Problem.** Two Tm-with-additives code paths exist:

1. `neoswga/core/thermodynamics.py:910` `calculate_tm_with_additives()` uses
   fixed coefficients, for example -2.3 C/M for betaine regardless of GC,
   which is only correct at ~50 percent GC.
2. `neoswga/core/reaction_conditions.py` +
   `neoswga/core/additives.py:554-637` implement a GC-aware sigmoid model
   with length scaling, which is correct.

Both paths are callable. `reaction_conditions.py:331` also carries a
linear saturation shortcut (`min(1.0, betaine_m / 5.2)`) that conflicts
with the sigmoid model.

**Changes.**

- `neoswga/core/thermodynamics.py`: mark `calculate_tm_with_additives`
  deprecated with a `DeprecationWarning` pointing to
  `ReactionConditions.calculate_tm_correction`. Keep the function body so
  downstream code does not hard-break.
- `neoswga/core/reaction_conditions.py:331` (`_calculate_gc_normalization`):
  remove the linear saturation shortcut. Route through the sigmoid in
  `additives.py` for all additives.
- Audit every caller: `grep -rn calculate_tm_with_additives neoswga/`.
  Route each to the canonical path. If a caller wants raw Tm, it should go
  through `ReactionConditions.calculate_effective_tm`.
- `neoswga/core/adaptive_filters.py` uses `calculate_tm_with_additives` per
  commit 1491807 replacement - verify the canonical path is wired and the
  legacy call is gone.

**Tests.**

- Property-based (hypothesis): for AT-rich primers (gc <= 0.35), betaine Tm
  correction is stronger (more negative) than for GC-rich (gc >= 0.65).
  Run over all four supported additives.
- Property-based: Tm correction is monotone in betaine concentration for a
  fixed primer.
- Unit: `DeprecationWarning` emitted when `calculate_tm_with_additives` is
  called directly.
- Regression: `tests/test_thermodynamics_properties.py` and
  `tests/test_reaction_conditions.py` continue to pass.

**Acceptance.** Single canonical Tm path used by all pipeline stages.
Property tests green. Deprecation warning emitted for legacy call.

### 0.4 Multi-target aggregation and max_gini enforcement

**Problem.** `neoswga/core/multi_genome_filter.py:345` aggregates multi-target
foreground frequencies with `np.mean`. A primer present in one target at
1e-3 and absent in another at 0 has mean 5e-4 and still passes the
`min_fg_freq` filter, defeating the intent of pan-primer design. `max_gini`
is a parameter but `score_primer` does not enforce it.

**Changes.**

- `neoswga/core/parameter.py`: add to `PipelineParameters`:
  - `multi_target_aggregation: Literal["min", "mean", "geomean"] = "min"`
- `neoswga/core/multi_genome_filter.py`: dispatch `target_freq` computation
  on the new field. `min` is the default (every target must meet the
  threshold). `geomean` is softer but still requires all targets to be
  non-zero. `mean` restores the prior behavior.
- Actually enforce `max_gini` in `score_primer`: reject primers with
  per-target Gini exceeding the threshold.
- `neoswga/core/param_validator.py`: warn when `multi_target_aggregation=
  "mean"` is set with two or more targets (describe the footgun).

**Tests.**

- Unit: `tests/test_multi_genome_filter.py` - given two targets, primer
  hits only target 1:
  - `aggregation="min"`: rejected.
  - `aggregation="geomean"`: rejected (contains zero).
  - `aggregation="mean"`: accepted.
- Unit: `max_gini=0.3`, primer with Gini 0.5 rejected; 0.2 accepted.
- Integration: existing multi-genome tests - set `multi_target_aggregation=
  "mean"` to preserve previous behavior. Add CHANGELOG note that the
  default has changed.

**Acceptance.** New aggregation test passes. `max_gini` enforcement test
passes. Existing multi-genome tests pass (with the explicit `mean`
opt-in where they relied on the old behavior).

## Cross-cutting acceptance

- All previously-passing tests still pass.
- Four new integration scenarios exist: low-GC target, blacklist,
  multi-target-min, and a combined smoke scenario exercising all four fixes
  together.
- `CHANGELOG.md` entry for v3.7.0 with:
  - Bug fixes for adaptive GC wiring and blacklist integration.
  - Breaking behavior change: `multi_target_aggregation` default is now
    `min`. Users with two or more targets should review results.
  - Deprecation notice for `calculate_tm_with_additives`.
- Version bump in `neoswga/__init__.py`, `pyproject.toml`, `VERSION`,
  `docs/CHANGELOG.md` (the single-source-of-truth pattern from commit
  5d3f2aa).

## Testing strategy

Use test-driven development per the superpowers workflow.

- Start each work stream by writing the failing test, then fix.
- For 0.1, 0.2, 0.4, the failing integration test is the most important
  deliverable - it is the regression guard for the rest of the roadmap.
- For 0.3, property-based tests (hypothesis) are the most valuable because
  the bug is in numeric values, not control flow.

Fixture files:

- `tests/fixtures/low_gc/target.fna` - synthetic ~33 percent GC 50 kbp
  sequence.
- `tests/fixtures/with_blacklist/blacklist.fna` - a short sequence sharing
  some k-mers with the plasmid example target.
- `tests/fixtures/multi_target/target_a.fna`,
  `tests/fixtures/multi_target/target_b.fna` - two targets with partial
  k-mer overlap.

Pre-compute k-mer counts for these fixtures so CI does not need Jellyfish
(Phase 2 fixes the CI Jellyfish story; Phase 0 tests stay Jellyfish-free).

## Backwards compatibility

- All new `PipelineParameters` fields are additive with defaults, so existing
  `params.json` files load.
- `multi_target_aggregation` default change is a behavior change. Document in
  CHANGELOG. Single-target users see no change. Multi-target users should
  re-run; results may be stricter.
- `calculate_tm_with_additives` is deprecated, not removed. Scheduled for
  removal in v4.0.0 (documented in Phase 5 deprecation enforcement).

## Risks and mitigations

- **Risk.** Collapsing Tm paths reveals callers in experimental modules.
  *Mitigation.* Keep the deprecated function callable with a warning;
  concrete removal is Phase 5.
- **Risk.** Min-aggregation rejects primers that users previously relied on.
  *Mitigation.* `mean` is an opt-in escape hatch via `params.json`.
  Document in CHANGELOG.
- **Risk.** Adaptive GC filter admits too many candidates for balanced-GC
  genomes.
  *Mitigation.* `gc_tolerance=0.15` matches the implemented default;
  `use_adaptive_gc=False` restores previous behavior exactly.
- **Risk.** New fixture files bloat the repo.
  *Mitigation.* Target fixtures under 100 kbp each; pre-computed k-mer
  files are small text.

## Exit criteria (restated, for the plan)

1. `tests/fixtures/low_gc/` scenario passes and returns non-empty
   `step2_df`.
2. `tests/fixtures/with_blacklist/` scenario produces fewer
   blacklist-hitting primers with the blacklist than without.
3. Four committed blacklist/library tests green.
4. Property-based Tm tests assert GC-asymmetric betaine effect and
   monotonicity in concentration.
5. Multi-target `min` aggregation rejects primers absent in one target.
6. `max_gini` enforcement rejects high-Gini primers.
7. All pre-existing tests pass.
8. `CHANGELOG.md` and version bumped to 3.7.0.
9. No new CLI commands. No new external dependencies.

## Next step after approval

Invoke the writing-plans skill to convert this design into a stepwise
implementation plan with explicit file-by-file TDD iterations.
