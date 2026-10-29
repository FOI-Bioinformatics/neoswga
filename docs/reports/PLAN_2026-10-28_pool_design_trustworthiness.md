# Plan: acting on the pool-design trustworthiness SWOT

Against [the SWOT](SWOT_2026-10-28_pool_design_trustworthiness.md), written at
`dec9047`. Its weaknesses, opportunities and threats become the seven phases
below, resequenced by dependency rather than by the order they appear there.

**What was reproduced before planning, and what was not.** W3's table, W5's
line counts and T1's 299/500 were measured for the SWOT and the probes are in
the tree (`docs/validation/hybrid_reaction_2026-10-26/probe.py`,
`docs/validation/clique_dimer_dg_2026-10-26/probe.py`). **W4 is the exception
and it governs Phase 5**: the coverage-inflation argument is read off
`coverage.py:28`, `:117` and `:234` and has never been measured, so Phase 5
begins with a measurement and not a fix.

**Status, 28 October 2026.** Phases 2 and 1 are implemented, in that order, for
the reason given under *Sequencing* below. Everything else is unstarted; the
table is the contract, and the `State` column is the only part of it that is a
progress report.

## Startability in this container

Three phases need something this container does not have. Stating it up front so
no phase is begun and abandoned:

| Needs | Phases blocked |
|---|---|
| A k-mer counter (KMC3 or jellyfish) to generate a candidate pool | 5 (magnitude), 6 |
| A wet-lab outcome | none — deliberately. No phase below promises prospective validation (W6), which stays open |

Phases 1, 2, 3, 4 and 7 are fully startable here. Phase 5's measurement script
can be written and unit-tested here and must be RUN where a counter exists.

## The contract

| Task | Phase | State | Pinned by |
|---|---|---|---|
| 1.1 Restore headroom in `parameter.py` | 1 | **Done** `e8962e1` — 1596→1453, headroom 148 | size + headroom ratchets |
| 1.2 Restore headroom in `unified_optimizer.py`, `hybrid_optimizer.py` | 1 | **Done** `e8962e1` — headroom 48 and 59 | as above |
| 2.1 Disclose the pre-2026-10-26 reaction defect | 2 | **Done** `3f648f3` | `test_reaction_provenance_is_reported.py`, CHANGELOG, LIMITATIONS |
| 3.1 Ratchet the adapter/inner parameter diff | 3 | **Done** `925916f` — covers 8 seams, not 4 | `test_no_adapter_drops_an_inner_parameter.py` |
| 3.2 Audit the two unaudited adapters against it | 3 | **Done** `925916f` — found 4 dropped values in `background-aware`; allowlist empty for all method adapters | same test + behavioural pin |
| 3.3 Stop the run log asserting weights the engine lacks | 3 | **Done** `f4fd1af` — also found all three inert on `dominating-set` and `clique` | `test_the_weights_log_does_not_overclaim.py` |
| 4.1 Ratchet the doc link/orphan sweep | 4 | **Done** `b0dc6ed` — found 12 more dead references; 7 repointed, 5 reworded | `test_docs_links_resolve.py` |
| 5.1 Measure the directional-coverage inflation | 5 | **Done** `58cdc00` — 1.27x-1.85x at 3 kb; synthetic genome, real panels | `docs/validation/2026-10-29-directional-coverage-measured.md` |
| 5.2 Decide and, if warranted, fix + re-baseline | 5 | **Done** `fbf2190` — default stays `symmetric`, both figures now reported and labelled, directional without an explicit reach refused. No re-baselining: nothing moved | `test_both_coverage_geometries_are_reported.py`, Known Issue 25 |
| 6.1 Measure the two Stage-2 terms in a discriminating regime | 6 | **Done** `0fd2c2d` — Tm term: measured, benefit not demonstrated. Host term: could not run, now refuses | two dated records |
| 7.1 Strand-balance report | 7 | **Done** `162a699` — all ten figures per run, labelled target/host; also found a `null` alternation score rendering as 0.00 | `test_the_strand_report_reads_every_genome.py` |
| 7.2 Compute and report `off_gap_gini` | 7 | **Done** `366df9b` — computed from positions the strand pass already fetched; 0.497 dispersed against 0.966 clustered at equal site count; in no score | `test_the_background_evenness_is_computed.py` |
| 7.3 Record what `max_dimer_dg` does at the default | 7 | **Done** `58cac68` | schema description + regenerated params-reference |

---

## Phase 1 — Restore headroom, because the next option does not fit (W5)

First because it blocks later phases. `parameter.py` is 1596/1601 and the last
new params.json key cost 6 lines there (`40af6e0`). Phase 7.2 adds a reported
metric and Phase 5.2 may add a coverage-geometry key; neither fits today.

**1.1** Extract from `parameter.py`. The module is a flat list of module-level
defaults plus `_apply_params_only_keys`. The search-control block already has a
natural home: `search_control.SEARCH_SETTING_DEFAULTS` holds the same keys'
defaults. Move the search-control globals' *definitions* behind that dict and
re-export, so `parameter.stage2_tm_weight` still resolves — CLAUDE.md requires
the module-level default to exist, not to be literal.

**1.2** `unified_optimizer.py` and `hybrid_optimizer.py` are both 1646/1650.
`OptimizationConfig` is a large dataclass in the former; the Stage-2 ranking
helpers already moved out of the latter into `stage2_ranking.py` and more can
follow the same way.

**Acceptance.** Every touched module at headroom ≥ 20, the existing ratchets
green, and no re-export removed: `test_public_api_exports.py` and the
`parameter` globals a schema key binds must still resolve.

**Risk.** This is the phase most likely to break something invisibly, because
`parameter` globals are read by name across the tree. Do it with no behaviour
change in the same commit, and diff the full suite against baseline.

## Phase 2 — Disclose the reaction defect (T1)

Independent of every other phase and the only one with reach outside the
repository. Not a code fix: a disclosure.

Every `equiphi29`, `bst` or `klenow` design produced before 2026-10-26 ran as
phi29 at 30 C, with the Stage-0 thermodynamic pre-screen skipped for the two
presets that set `thermo_filter=True`. Under the shipped `equiphi29_scenario`
reaction that screen keeps 299 of 500 candidates.

**2.1** Three places, because three audiences:
- `docs/CHANGELOG.md` — a fixed-defect entry naming the affected polymerases
  and the date boundary.
- `docs/LIMITATIONS.md` — what a panel designed before that date does not
  establish. This is the page the README tells users to read first.
- `neoswga doctor` (or `validate`) — when a `run_manifest.json` in `data_dir`
  records a non-phi29 polymerase and a step that predates the fix, say the
  design should be re-run. The manifest already stores `effective_conditions`
  per step, which is exactly the evidence needed.

**Acceptance.** A test that the doctor check fires on a manifest fixture
carrying a pre-fix hot-polymerase step and stays silent on a phi29 one.

**Scope.** Re-running affected designs is the user's call, not this plan's. The
deliverable is that they can tell whether they have one.

## Phase 3 — Make the defect class structural (W1, T2, O4, T5)

Four reachability defects in one pass, all through seams the six existing
ratchets do not inspect. The check that found two of them is cheap and
mechanical, which is the argument for encoding it.

**3.1** A ratchet that, for every adapter that constructs an inner engine, diffs
the inner's `__init__` parameters against what the adapter passes or `setattr`s.
Scope is four modules: `hybrid_optimizer`, `background_aware_optimizer`,
`dominating_set_adapter`, `network_optimizer`. Any parameter neither passed nor
set is a finding; the allowlist starts empty and may only shrink, like the other
reachability ratchets.

**3.2** Run it against `network_optimizer` and `background_aware_optimizer`,
which this session did not audit by hand. Fix what it finds, or allowlist with a
stated reason.

**3.3** The `--application` run log asserted "selection weights applied:
tm=0.30, uniformity=0.05, dimer_penalty=0.45" while the engine used 0.0 — a log
line that actively misleads is worse than silence. Either read the value off the
engine that will use it, or say the weights are inert on this method. A test
should pin that the figure logged is the figure the engine holds.

**Acceptance.** 3.1's test fails when the `resolve_reaction` call or either
`setattr` in `HybridBaseOptimizer` is removed — the same negative check the four
fixes used.

## Phase 4 — Guard the documentation (T4)

Three stale entries and a five-document dead link were found by hand this week,
all incidentally. The sweep that found them is twenty lines.

**4.1** A test asserting zero dangling relative `.md` links across the tree.
Orphans are NOT an error — `docs/README.md` states that `validation/` holds
records the curated table does not list — so the test checks links only. Inline
code-span paths (`` `docs/...md` ``) must be checked too: four of the five dead
references were code spans, which a link-only regex misses.

**Acceptance.** The test fails if any of the five dead references from
`92849cb` is reintroduced.

## Phase 5 — Measure the directional-coverage inflation, then decide (W4, O2)

The headline metric. `fg_coverage` is authoritative in every
`step4_improved_df_summary.json` and across `validation/`, and the argument that
it is inflated is structural and unmeasured. **Measurement precedes any fix, and
the fix is not pre-approved by this plan.**

**5.1** A script that, for a delivered panel, computes coverage both ways: the
shipped symmetric `[pos-r, pos+r]` window, and the directional model —
a literal occurrence covering `[pos, pos+r]` and a reverse-complement occurrence
`[pos-r+len, pos+len]`, one direction per occurrence. Report both per prefix, on
a dense panel and a sparse one, since the expectation is that the error grows as
panels sparsify. The multi-record references are the ones that also exercise
record-join geometry.

**5.2** Then decide. If the inflation is material, fixing it moves every
recorded coverage figure in this repository, so the change needs: a geometry key
so the old behaviour stays reachable for comparison, a re-baselining of the
summaries and validation records that quote the figure, and a note in
`reach_calibration.md`, whose reach was fitted *against* the symmetric model and
may need refitting. If it is immaterial on both panels, record that and close
W4 with a measurement instead of a fix.

**Blocked here.** 5.1 needs a candidate pool and position indexes.

## Phase 6 — Measure the two Stage-2 terms where they can work (O1, W2)

`stage2_tm_weight` and `stage2_background_weight` both ship off, each with
"unmeasured against a reaction" as the stated reason. Until 2026-10-26 no
discriminating reaction was reachable through the default method, so the
measurement was not merely undone but impossible.

**6.1** With `polymerase: equiphi29` or `bst` now actually reaching the
optimizer, measure both terms on a hot-polymerase pool: delivered coverage,
selectivity density, mean Tm deviation and host coverage, at several weights,
against the weight-0 control. Known Issue 17 predicts the phi29 pool cannot
separate; the question is whether a hot pool can.

**Outcome is open.** A result showing no benefit is a result: it keeps both keys
off and converts "unmeasured" into "measured, no benefit demonstrated", which is
the honest end state the project's own convention asks for.

**Blocked here.** Needs a counter.

## Phase 7 — The cheap gaps and the two decisions (O3, W3, T3)

**7.1 Strand-balance report.** No new data: `strand_metrics.collect_strand_stats`
already produces all five figures per genome onto `PrimerSetMetrics.strand_stats`.
Report a balance figure and flag an imbalanced panel. COATswga ranks on the
*mean* of its two strand coverages, so 0.99/0.30 beats 0.60/0.60 — that is the
failure this report would surface. Known Issue 18 keeps the strand scalars
unconstrainable because an unmeasured zero cannot be told from a measured one;
that reason does not extend to reporting.

**7.2 `off_gap_gini`.** The second-largest weight in swga 2.0's fitted model and
computed nowhere here. Compute and report it. Do **not** fold it into
`normalized_score`: the sample size does not support a weight, and
`tool_comparison.md` already states that position.

**7.3 Record what `max_dimer_dg` actually does.** At the shipped
`max_dimer_bp=3` a −6.0 and a −3.0 floor reject nothing the base-pair screen has
not already rejected (SWOT W3). The key is honest and nearly inert. Say so where
it is documented, so a user does not set it believing it tightened anything, and
note the regime where it bites.

**7.4 The grade caveat (T3).** Deferred, not planned. The grades are calibrated
retrospectively and documented as descriptive, and NeoSWGA is the only one of
the four tools that grades at all. Softening them is a judgement about
presentation that needs the maintainer's call, not a task.

---

## Sequencing, and what to do first

1. **Phase 2** before anything else if any hot-polymerase design has left this
   repository. It is a disclosure, it is cheap, and it is the only item with a
   consequence already in the world.
2. **Phase 1** next, because Phases 7.2 and 5.2 cannot add what they need until
   it lands.
3. **Phase 3** and **Phase 4** are independent, cheap, and each prevents a
   recurrence rather than fixing an instance. Phase 3 is the higher value: it
   encodes the judgement that found four defects in a day.
4. **Phase 5** is the largest and the one with the most downstream churn. It is
   also the one where the project's own comparative review already concedes the
   defect, so it is the strongest claim available — once measured.
5. **Phase 6** is the measurement that was impossible last week.

## What this plan does not do

- It does not promise prospective validation (W6). Nothing here produces a
  wet-lab outcome, and T3's gap between presentation and evidence is narrowed
  only by Phase 7.4, which is deferred to a decision.
- It does not pre-approve the directional-coverage fix. Phase 5 can end in a
  recorded measurement and no code change.
- It does not address W7 (the suite is green only in CI). A thin container is a
  contributor-experience problem, and the baseline-diff habit this session used
  is a workaround that works; encoding it is lower value than any phase above.
- It assigns no dates. Every phase's acceptance criterion is a test or a dated
  record, which is the project's unit of done.
