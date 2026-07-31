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

**Klenow now violates this too.** Correcting its processivity from 10000 bp to
40 nt (entry #4) exposed a violation the wrong figure had been hiding: the
set-cover stage credits Klenow with 5000 bp against ~500 bp of effective reach,
a 10x over-credit.

A polymerase cannot extend further than it stays bound. The set-cover stage
therefore credits each bst primer with 5x more reach than the enzyme has, and 10x
the realistic per-primer reach the result is later *scored* on.

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

Knock-on effect: this exposed the Klenow half of entry #1 above, which the wrong
10000 bp figure had been masking.

---

## 5. phi29 extension rate was ~3x too fast — FIXED

`extension_rate_nt_s` was 150, annotated "median of 100-170 range". Blanco et al.
(1989) — the same paper cited for the 70 kb processivity — reports ~53 nt/s, and
single-molecule work agrees at ~50 nt/s. **Corrected to 53** in both the registry
and `replication_simulator.ReplicationConfig`.

Affects simulation timings only, not primer selection.

---

## 6. The magnesium model was calibrated for PCR, not MDA — FIXED

`mg_optimal` was 2.5 mM with `mg_high_threshold` 6.0 mM — PCR figures. Since
`mg_default_mm` correctly puts phi29 at the vendor 10 mM, the model penalised the
pipeline's own default configuration as having too much magnesium.

Worse, the test guarding the old value cited **"Rahman et al. (2014) PLoS One
9:e112515"**. That attribution is wrong twice over: the DOI resolves to an
unrelated rice chromatin-remodelling paper, and the real Rahman et al. (2014)
(*Anwer Khan Modern Medical College Journal* 4(1):30-36) is a general **PCR**
review — not an authority for isothermal strand-displacement chemistry. A
mis-attributed, domain-inappropriate citation was lending false authority to a
wrong number.

**Corrected to** `mg_optimal = 10.0`, `mg_low_threshold = 4.0`,
`mg_high_threshold = 15.0`, matching the standard phi29 buffer (Thermo
MAN0030290, NEB phi29 buffer). The band is wide because dNTPs chelate Mg2+
roughly 1:1, so free Mg2+ sits several mM below nominal total.

The replacement test asserts that no polymerase's own default falls outside the
model's usable band — the invariant the old value violated.
