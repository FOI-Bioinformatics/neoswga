# Known inconsistencies in the seeded polymerase data

The registry was seeded **verbatim** from the tables it replaces, so that landing it
changed no behaviour. That faithfulness means it also inherited the contradictions
those tables had accumulated.

This file is the ledger. Each entry is tracked by an assertion in
`tests/test_registry_invariants.py`; the ones that fail today are marked
`xfail(strict=True)`, so fixing one turns the test green and forces the ledger to be
updated rather than letting the fix pass unnoticed.

**None of these should be "cleaned up" as part of a refactor.** Each changes numeric
output, and several change which primers get selected. They need a dedicated change
with integration baselines regenerated deliberately.

---

## 1. bst set-cover reach exceeds bst processivity — physically impossible

| Field | bst value | Source |
|---|---|---|
| `legacy_hybrid_max_extension` | 10000 | `hybrid_optimizer.POLYMERASE_PRESETS` |
| `processivity_bp` | 2000 | `reaction_conditions.POLYMERASE_CHARACTERISTICS` |
| `typical_amplicon_bp` | 1000 | same |

A polymerase cannot extend further than it stays bound. The set-cover stage
therefore credits each bst primer with 5x more reach than the enzyme has, and 10x
the realistic per-primer reach the result is later *scored* on.

klenow is inverted the other way: `legacy_hybrid_max_extension` 5000 vs
`processivity_bp` 10000.

**Why not fixed here:** `max_extension` (set-cover reach) and `processivity`
(single-molecule reach) are genuinely different quantities, and the project's own
coverage-reach doctrine says selection should use the *realistic* reach
(`coverage.polymerase_extension_reach('realistic')`, i.e. 1000 for bst, 1500 for
klenow) — which matches neither column. Changing bst from 10000 re-ranks every
hybrid-optimizer result for that polymerase.

**Recommended resolution:** point Stage-1 set-cover at `typical_amplicon_bp` for all
polymerases, so selection and scoring agree, and regenerate the integration
baselines.

---

## 2. klenow warn band is wider than the hard-reject band

| Band | klenow value | Source |
|---|---|---|
| `temp_warn_range` | (20.0, 42.0) | `param_validator.POLYMERASE_TEMP_RANGES` |
| `temp_hard_range` | (25.0, 40.0) | `ReactionConditions._validate` |

A klenow run at 41 °C passes `neoswga validate-params` cleanly and then raises
`ValueError` inside `ReactionConditions.__init__`. Same at 21–24 °C. The validator's
job is to catch exactly this before the run starts.

**Recommended resolution:** narrow `temp_warn_range` to sit inside
`temp_hard_range`. Low risk — it only affects warning text — but it does change
validator output, so it is grouped here.

---

## 3. bst operating temperature disagrees with its optimum

| Field | bst value | Source |
|---|---|---|
| `optimal_temp` | 63.0 | `reaction_conditions`, `mechanistic_params` |
| `preset_reaction_temp` | 60.0 | `hybrid_optimizer`, `additive_optimizer` |

Two subsystems run bst 3 °C below the optimum the rest of the codebase advertises.
Possibly deliberate (a conservative operating point), possibly drift; there is no
comment either way. phi29, equiphi29 and klenow agree across both fields.

**Recommended resolution:** decide which is intended and collapse to one field.

---

## 4. klenow processivity contradicts the literature by ~250x

`processivity_bp = 10000` cites Bambara et al. (1978) JBC 253:413. That paper is
about the processive mechanism of DNA polymerase I and does not support a 10 kb
value. Klenow is a distributive enzyme: it dissociates after roughly 5–40
nucleotides, with 42 nt measured directly by single-molecule nanocircuit recording
(*J. Am. Chem. Soc.* 2013, PMC3738269).

This is not an internal inconsistency but an external one, and it is the largest
single numeric error found in the audit. It feeds `coverage.polymerase_extension_reach`
and hence klenow coverage estimates.

**Recommended resolution:** correct to a literature value and re-derive
`typical_amplicon_bp`, as part of the scientific-corrections phase with a
`schema_version` bump.

---

## 5. phi29 extension rate is ~3x too fast

`extension_rate_nt_s = 150`, annotated "median of 100-170 range". Blanco et al.
(1989) — the same paper cited for the 70 kb processivity — reports ~53 nt/s. This
affects simulation timings in `replication_simulator` and `stochastic_simulator`, not
primer selection.

---

## 6. The magnesium model is calibrated for PCR, not MDA

`mechanistic_params.MECHANISTIC_MODEL_PARAMS["enzyme"]` sets `mg_optimal = 2.5` mM
and `mg_high_threshold = 6.0` mM, which are PCR figures. The standard phi29 buffer is
10 mM MgCl₂ (with 10 mM (NH₄)₂SO₄ and 4 mM DTT), and `mg_default_mm` correctly
defaults phi29 to 10.0.

The consequence is self-contradictory: the pipeline's own default configuration is
penalised by its own enzyme-activity model as having too much magnesium.

**Recommended resolution:** recalibrate the Mg response for strand-displacement
isothermal amplification.
