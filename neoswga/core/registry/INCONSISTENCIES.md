# Known inconsistencies in the seeded polymerase data

The registry was seeded **verbatim** from the tables it replaces, so that landing it
changed no behaviour. That faithfulness means it also inherited the contradictions
those tables had accumulated.

This file is the ledger. Each entry is tracked by an assertion in
`tests/test_registry_invariants.py`.

**Status: entries 1-6 are fixed; entry 7 is open by deliberate choice.** They
were not fixed during the registry refactor, because each changes numeric output
and several change which primers get selected; they were carried as
`xfail(strict=True)` until they could be addressed with the whole suite
watching. No strict-xfail entries remain — #7 is documented rather than
pinned, because it is a gap in the model rather than a demonstrated defect
(see below).

Kept as a record of what was wrong and why, because several of these were
long-standing and the reasoning matters more than the diffs.

---

## 1. Network reach disagreed with processivity — FIXED

`hybrid_optimizer.POLYMERASE_PRESETS[...].max_extension` was a hand-maintained
number that disagreed with `processivity_bp`: bst 10000 against a processivity of
2000 (a reach the enzyme cannot achieve), klenow 5000 against a corrected 40 nt.

The original ledger entry described this as the *set-cover* reach. That was
wrong. `max_extension` feeds the Stage-2 amplification-**network** connectivity
question — "could two primers connect via one uninterrupted extension?" — for
which single-molecule processivity is exactly the right quantity. The Stage-1
set-cover **coverage** objective uses `coverage_reach`, threaded separately by
`unified_optimizer` and already set to the realistic per-primer reach.

**Fixed by** deriving `max_extension` from `processivity_bp` and deleting the
separate `legacy_hybrid_max_extension` field, so the two cannot diverge again.
Affects bst, bst3.0, bsu and klenow; phi29 and equiphi29 already agreed.

## 2. klenow warn band was wider than the hard-reject band — FIXED

`temp_warn_range` was (20.0, 42.0) while `temp_hard_range` is (25.0, 40.0), so a
klenow run at 41 C passed `neoswga validate-params` cleanly and then raised
`ValueError` inside `ReactionConditions.__init__`. Catching exactly that before
the run starts is the validator's job.

**Fixed by** narrowing `temp_warn_range` to (25.0, 40.0). An invariant test now
asserts the warn band sits inside the hard band for every polymerase.

## 3. bst operating temperature disagreed with its optimum — FIXED

`preset_reaction_temp` was 60.0 while `optimal_temp` is 63.0, so
`hybrid_optimizer` and `additive_optimizer` ran bst 3 C below the optimum the
rest of the codebase advertised. No comment explained the difference and no other
polymerase diverged this way, so it read as drift rather than a deliberate
conservative operating point.

**Fixed by** collapsing `preset_reaction_temp` to `optimal_temp` (63.0).

## 4. klenow processivity contradicted the literature by ~250x — FIXED

`processivity_bp` was 10000 bp citing Bambara et al. (1978) JBC 253:413, which
does not support it. Klenow is **distributive**: it dissociates after roughly
5-40 nucleotides, with 42 nt measured by single-molecule nanocircuit recording
(*J. Am. Chem. Soc.* 2013, PMC3738269).

**Corrected to** `processivity_bp = 40`, `typical_amplicon_bp = 500`, plus a new
`distributive = True` flag.

The two fields now differ in kind, deliberately. A distributive enzyme re-binds
and continues, so its effective reach over a long isothermal incubation
legitimately exceeds single-event processivity — `typical_amplicon_bp >
processivity_bp` is correct here, not a contradiction. The invariant test skips
distributive enzymes for that comparison. The 500 bp figure is a conservative
estimate for a multi-primer reaction, not a measurement.

## 5. phi29 extension rate was ~3x too fast — FIXED

`extension_rate_nt_s` was 150, annotated "median of 100-170 range". Blanco et al.
(1989) — the same paper cited for the 70 kb processivity — reports ~53 nt/s, and
single-molecule work agrees at ~50 nt/s. **Corrected to 53** in both the registry
and `replication_simulator.ReplicationConfig`.

Affects simulation timings only, not primer selection.

## 6. The magnesium model was calibrated for PCR, not MDA — FIXED

`mg_optimal` was 2.5 mM with `mg_high_threshold` 6.0 mM — PCR figures. Since
`mg_default_mm` correctly puts phi29 at the vendor 10 mM, the model penalised the
pipeline's own default configuration as having too much magnesium.

Worse, the test guarding the old value cited **"Rahman et al. (2014) PLoS One
9:e112515"**. That attribution is wrong twice over: the DOI resolves to an
unrelated rice chromatin-remodelling paper, and the real Rahman et al. (2014)
(*Anwer Khan Modern Medical College Journal* 4(1):30-36) is a general **PCR**
review — not an authority for isothermal strand-displacement chemistry.

**Corrected to** `mg_optimal = 10.0`, `mg_low_threshold = 4.0`,
`mg_high_threshold = 15.0`, matching the standard phi29 buffer (Thermo
MAN0030290, NEB phi29 buffer). The replacement test asserts that no polymerase's
own default falls outside the model's usable band.

## 7. `normalized_score` computes `max_gap` and never uses it — OPEN (by choice)

`PrimerSetMetrics.normalized_score` weights coverage, selectivity, dimer risk,
evenness and Tm. It computes `max_gap` and then ignores it — even though on the
only wet-lab-validated benchmark available (six *Prevotella* sets, Dwivedi-Yu
et al. 2023) `max_gap` is the strongest single predictor of measured enrichment
(Spearman rho -0.83, and ranking by it alone reproduces the measured top three in
order), while mean binding distance has essentially none (rho -0.09).

**This is a gap in the model, not a demonstrated defect in the ranking.** The
shipped scoring does rank the two winners first on that benchmark. An earlier
version of this entry claimed otherwise; that came from a bug in the validation
harness (background sites counted over the 3.2 Mb target rather than the 3.1 Gb
human background) and has been corrected.

A coverage-hole term based on `max_gap` was implemented and measured against the
benchmark. It changed nothing in the ranking, so it was **not** shipped:
re-ranking all optimizer output requires demonstrating a benefit, and on the
available data there is none to show.

**Recommended resolution:** extend the benchmark first (Clarke et al. 2017
*M. tuberculosis*, or in-house sets with outcomes), then re-test a hole term.
Full write-up in `docs/validation/published_primer_sets.md`.
