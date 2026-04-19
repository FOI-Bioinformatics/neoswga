# neoswga Production-Readiness Roadmap

**Status:** Draft for review
**Date:** 2026-04-19
**Current version:** 3.6.0 (Development Status: 4 - Beta)
**Target audience:** Tier 1 — single research lab / solo bioinformatician via PyPI (now).
Tier 2 — public community tool with stable API and reproducibility (later).

## Scope framing

"Production-ready" is decomposed into six ordered phases. Each phase has a bounded scope
and a testable exit criterion. Each becomes its own design spec and implementation plan.

The phases address, in order: correctness bugs that invalidate results today; reproducibility
so future iteration is possible; a validation harness that locks in correctness across the
parameter matrix the tool claims to support; a polymerase-aware scoring model; the
wet-lab feedback loop (the distinctive capability the user asked for); and finally the
community-tool hardening needed to graduate from Beta.

The parameter matrix the tool must eventually handle correctly:

- Polymerases: phi29 (30 C, 6-12 bp primers), equiphi29 (42-45 C, 12-18 bp primers).
- Additives: DMSO, betaine, formamide, urea, TMAC, trehalose, ethanol, Mg2+, Na+.
- GC content: 32-68 percent (AT-rich, balanced, GC-rich bacteria).
- Genome sizes: ~5 kbp plasmids through ~10 Mbp bacteria, with ~3 Gbp human background.
- Multiple foreground targets (pan-primer design).
- Multiple background genomes.
- Blacklist genomes (must not be amplified).

## Audit summary

Evidence was collected by three parallel audits against the v3.6.0 codebase.
Key findings, grouped by severity.

### Critical (produces incorrect results)

1. `random_forest_filter.p` is polymerase-agnostic and its training provenance is
   unclear. equiphi29 predictions at 42 C may be systematically biased.
   `neoswga/core/rf_preprocessing.py:18-34`.
2. `AdaptiveGCFilter` exists but is never instantiated by the pipeline. Extreme-GC
   genomes (Francisella, Plasmodium) are rejected wholesale against the fixed
   37.5-62.5 percent thresholds. `neoswga/core/adaptive_filters.py` vs `neoswga/core/filter.py`.
3. Blacklist half-integrated: `bl_seq_lengths`, `bl_penalty`, `max_bl_freq` missing
   from `PipelineParameters` defaults; `pipeline.py:659-660` reads a field never set;
   `wizard.py` writes fields the param system does not accept.
4. Two conflicting Tm paths: `thermodynamics.py:910` uses GC-blind betaine coefficient;
   `additives.py:554-637` uses the correct GC-aware sigmoid. Both callable.
5. Multi-target aggregation uses `np.mean` for foreground
   (`multi_genome_filter.py:345`). A primer absent in one target still passes.
   `max_gini` advertised but not enforced in `score_primer`.

### Engineering

6. No run manifest. No persisted version, git SHA, seed, input checksums, or
   Jellyfish version beside outputs.
7. `seed` not a `params.json` field. Only `random`/`numpy` seeded in GA;
   `networkx`, `scipy`, `scikit-learn` unseeded. Determinism is aspirational.
8. CI does not install Jellyfish and does not run integration tests; the E2E
   plasmid test is not exercised in CI. `.github/workflows/ci.yml:48`.
9. No Windows CI. No conda-forge recipe. Deprecation policy is prose in README,
   not enforced by a decorator.

### Missing capability

10. No BAM/FASTQ ingestion. No pysam. `active_learning.py` exists but is
    experimental, in-silico, expects hand-entered `ExperimentalResult`, and is
    not CLI-wired.

### Solid foundations (do not disturb)

Exception hierarchy (`neoswga/core/exceptions.py`), polymerase database
(`reaction_conditions.py:52-97`), logging infrastructure, unsafe-pickle removal,
property-based thermodynamics tests, E2E plasmid test fixture, pyproject
versioning, new `genome_library.py`.

## Phased plan

```
Phase 0  -> Phase 1  -> Phase 2  -> Phase 3 (RF)
                                 \
                                  -> Phase 4 (feedback loop)
                                 -> Phase 5 -> Phase 6
```

### Phase 0 - Stop the bleeding

**Scope.** Fix the five [BUG]/[CRITICAL] correctness issues above, minus the RF
model (that is Phase 3 because it requires Phase 2's harness to validate
changes). No new features.

**Work streams.**

- 0.1 Wire `AdaptiveGCFilter` into the filter stage; add
  `use_adaptive_gc` and `gc_tolerance` to `PipelineParameters`.
- 0.2 Complete blacklist integration: add missing fields to `PipelineParameters`,
  fix `bl_seq_lengths` computation, wire counts into `MultiGenomeFilter`,
  ensure wizard-produced configs load, make the untracked
  `test_blacklist_*.py` tests pass and commit them.
- 0.3 Collapse thermodynamics to a single canonical Tm path.
  Deprecate `calculate_tm_with_additives` with a `DeprecationWarning`;
  remove the linear saturation shortcut in `reaction_conditions.py:331`;
  ensure every caller routes through `ReactionConditions.calculate_tm_correction`.
- 0.4 Multi-target aggregation: add
  `multi_target_aggregation: Literal["min", "mean", "geomean"] = "min"` to
  `PipelineParameters`. Default is `min` (a pan-primer must hit every target);
  `mean` is an opt-in escape hatch for backwards compatibility.
  Actually enforce `max_gini` in `score_primer`.

**Exit criteria.**

- All five issues closed with new regression tests.
- A Francisella-like (low-GC) target yields a non-empty `step2_df` (it would
  return zero primers today).
- A run with a blacklist genome produces strictly fewer blacklist-hitting
  primers in `step4_improved_df` than the same run without.
- Every existing test still passes.
- CHANGELOG entry for v3.7.0 with the behavior change for
  `multi_target_aggregation` clearly documented.

**Out of scope.** RF model changes, reproducibility work, new CLI commands.

### Phase 1 - Reproducibility baseline

**Scope.** Make runs identifiable and repeatable before any iterative science
work rests on them.

- Write `run_manifest.json` next to every pipeline output. Contents: neoswga
  version, git SHA if dev install, Jellyfish version, resolved RNG seed,
  sha256 of every input FASTA and k-mer file, resolved `params.json`, CLI
  invocation, UTC timestamp.
- Add `seed: int | None = None` to `PipelineParameters`.
  When None, generate one and record it. Thread the seed through
  `genetic_algorithm.py`, `moea_optimizer.py`, `network_optimizer.py`,
  `dominating_set_optimizer.py`, and into numpy, random, networkx, and scipy.
- CI check: same seed + same inputs -> byte-identical `step4_improved_df.csv`
  for the plasmid example. Run twice in CI, diff.

**Exit criteria.** Determinism CI check green. Manifest schema documented
in `docs/reproducibility.md`.

### Phase 2 - Scientific validation harness

**Scope.** A test matrix that locks in correctness across the parameter space
the tool claims to support. Without this, every subsequent change is guesswork.

- Matrix fixtures under `tests/matrix/`:
  `{phi29, equiphi29} x {low_gc, balanced, high_gc} x {small, medium, large} x
  {single, multi_target} x {no_bg, small_bg, large_bg} x {no_bl, with_bl}`.
  Not the full Cartesian product initially - pick ~20 cells that cover the
  edges and each dimension independently.
- Each cell: runs end-to-end (count-kmers -> optimize), compares against a
  committed golden `step4_improved_df.csv`, asserts scientific invariants.
- Invariants (partial list): Tm monotone in Na+ and Mg2+; additive corrections
  opposite-signed for AT-rich vs GC-rich; enrichment finite in host-free mode;
  non-empty filter output on legitimate inputs; multi-target coverage sums
  match per-target contributions.
- CI installs Jellyfish (Linux via apt, macOS via brew). Integration tier runs
  on a separate CI job, not gated on the lint+unit job.

**Exit criteria.** Matrix green in CI. A sample regression PR (e.g. "change
Tm formula by 0.5 C") turns it red.

### Phase 3 - RF model polymerase-awareness

**Scope.** Produce trustworthy scoring across the polymerase matrix.

- Audit `models/random_forest_filter.p` training data provenance. Document
  the validity envelope (which polymerase, temp, GC range, primer lengths).
- Per the user's direction: **ship separate models per polymerase**, not a
  unified model with a polymerase indicator.
- Produce `random_forest_phi29.p` and `random_forest_equiphi29.p`.
  Select at runtime based on `PipelineParameters.polymerase`.
- Warn when users run outside the validated envelope.
- Use the Phase 2 harness to verify no regression on phi29 cells and
  measurable improvement on equiphi29 cells.
- Document retraining procedure in `scripts/retrain_rf_model.py` and
  `docs/rf_model.md`.

**Blocker.** Requires training data. If retraining data is not available,
fall back to a validity-envelope warning + document as a known limitation
and escalate Phase 3 as a separate wet-lab project.

**Exit criteria.** Matrix harness green with per-polymerase models.
Validity warning fires on out-of-envelope runs.

### Phase 4 - Wet-lab feedback loop

**Scope.** Close the design -> wet-lab -> sequence -> re-design loop. New
subsystem. This is the flagship capability.

- `neoswga/core/coverage_ingest.py`: BAM/FASTQ -> per-position depth and
  per-region coverage. pysam dependency under a new `[feedback]` optional
  extra.
- `neoswga/core/coverage_attribution.py`: attribute observed coverage to
  primers via `PositionCache`.
- `neoswga/core/observed_analysis.py`: predicted vs observed report;
  attribute discrepancies to GC bias, repeats, secondary structure, etc.
- `neoswga/core/empirical_rescorer.py`: boost and penalize primers in
  `step3_df` based on observed outcomes.
- `PipelineParameters` extensions: `observed_coverage_bam`, `feedback_weight`,
  `rejected_primers_file`, `boosted_primers_file`.
- CLI: `neoswga ingest-coverage`, `neoswga analyze-feedback`,
  `neoswga redesign --feedback`.
- Wire `active_learning.py` to the CLI; ingest observed metrics
  automatically instead of expecting hand-entered results.
- Worked example under `examples/feedback_loop/` with a toy BAM,
  feedback report, and redesign params.
- Depends on Phase 1 (runs must be identifiable across iterations) and
  Phase 2 (regression safety net).

**Exit criteria.** End-to-end feedback roundtrip demonstrated on the
plasmid example. Redesigned primer set measurably different from the
initial set under a synthetic observed-underperformance scenario.

### Phase 5 - Tier 2 community readiness

**Parked per user direction.** Deferred until after Phase 4 ships. When
revived, scope is: Windows CI, coverage floor 50 -> 80 percent, conda-forge
recipe, `@deprecated` decorator and enforcement, public-API Sphinx docs
with stability markers, benchmark budgets in CI.

### Phase 6 - UX and docs

**Continuous.** Not sequenced. As each of Phase 0-4 lands, a small UX/docs
pass follows: error-message audit for the affected surface, docs
reconciled against code, wizard updates, tutorial updates.

## Dependency graph

```
Phase 0 -> Phase 1 -> Phase 2 -> Phase 3
                              \-> Phase 4
```

Phase 3 and Phase 4 can proceed in parallel after Phase 2 exits, but Phase 4
alone depends on Phase 1.

## Rough effort (solo)

- Phase 0: 2-3 weeks
- Phase 1: 1-2 weeks
- Phase 2: 3-4 weeks
- Phase 3: 3-4 weeks (assuming retraining data available; otherwise escalates)
- Phase 4: 6-8 weeks
- Phase 5: 3-4 weeks (parked)

Serial total through Phase 4: ~16-21 weeks. Parallelizing Phase 3 and Phase 4
after Phase 2 brings it to ~13-17 weeks.

## Decisions locked at brainstorming

- Phase boundaries as above.
- Multi-target aggregation default: `min`, with `mean`/`geomean` selectable.
- Phase 3: separate models per polymerase, not a unified indicator model.
- Phase 5 parked until after Phase 4.

## Open items deferred to phase-level design

- Exact schema of `run_manifest.json` (Phase 1).
- Matrix cell selection and golden-file strategy (Phase 2).
- Retraining data source for Phase 3.
- Feedback report JSON schema (Phase 4).

## Next step

Phase 0 design spec is written as a companion doc. After user review, it
proceeds to an implementation plan via the writing-plans skill.
