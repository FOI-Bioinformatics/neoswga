# Audit, gap and SWOT: chemistry, additives and temperature in pool design

**Date:** 2026-10-26
**Subject:** how reaction chemistry — the enzyme, the temperature, the buffer and
the additives — reaches the **selection** of an oligo pool, and what stands
between the code as it is and "select and optimise oligo pools with chemistry,
additives and temperature taken into account during design".

**Tree:** `d27e137`, working tree clean at snapshot apart from this report.

**How this was measured.** The container holds the source only. No reference
genome is present (`tests/validation/genomes/` does not exist) and neither k-mer
counter is on PATH, so no step of the four-step pipeline can run here and no
design against a real target could be produced. The Python dependencies were
absent when the audit began and were installed partway through, which is why
gap 2 carries a measurement and the surrounding sections do not.

Everything below is either (a) read from code with a `file:line`, (b) quoted
from a dated record already in the repository, (c) measured by
`chemistry_audit_2026-10-26/probes.py`, which is genome-free by construction, or
(d) measured by running a test, which is stated where it happens. The probes:

```
PYTHONPATH=. python docs/reports/chemistry_audit_2026-10-26/probes.py
```

Where a question needs a genome or a test run, it is marked open rather than
answered. Section 5 lists those explicitly.

**Relationship to the existing record.** This is not a re-audit of the additive
coefficients.
[`additives_coverage_audit_2026-09-14`](../validation/additives_coverage_audit_2026-09-14/README.md)
did that and [`PLAN_2026-09-14_additives_and_coverage.md`](PLAN_2026-09-14_additives_and_coverage.md)
closed its phases 1 to 3. The question here is the one that plan left open: given
a chemistry model of known and limited strength, how much of it reaches the
**choice of oligos**, and what would it take for chemistry to be a design
variable rather than a design input.
[`additive_specificity.md`](../validation/additive_specificity.md) is the
load-bearing prior measurement and is quoted rather than repeated.

---

## 1. How chemistry reaches a design today

Five doors, and it matters which is which, because the promise of "design under
my chemistry" is carried by only some of them.

```
params.json  (chemistry flags on `filter` only -- see gap 1)
      |
      v
resolve_design_request --> model_evidence.require_model_support   (door 1: refusal)
      |                    unknown enzyme, out-of-range oligo length,
      |                    additive with no model -> UnsupportedModelError
      v
ReactionConditions  (the frozen reaction, hashed into the request)
      |
      +--> filter            Tm window gate on each candidate            (door 2)
      |
      +--> optimize Stage 0  thermodynamic screen: Tm bounds + structure (door 3)
      |                      hybrid_thermo_screen.py:26-40
      |
      +--> compute_metrics   occupancy-weighted selectivity and          (door 4)
      |                      effective_fg_coverage -- the figure the
      |                      objective actually judges on
      |                      base_optimizer.py:313-322
      |
      +--> run_manifest      effective_conditions: the reaction that ran (door 5)
```

Doors 2, 3 and 4 are real and load-bearing. Door 4 is the strongest part of the
design: `PoolObjective.coverage` returns `effective_fg_coverage` by default, and
returns `None` rather than zero when no conditions were attached, because "no
reaction conditions means no temperature at which to evaluate occupancy, and a
fabricated number there would be indistinguishable from a measured one"
(`pool_objective.py:206-213`).

What no door does is **choose** a chemistry. In `design` and `optimize` the
reaction is an input, fixed before the first index opens. Chemistry is varied in
only two places, both outside that path:

| Path | Varies | Compares | Reachable by |
|---|---|---|---|
| `condition_sweep.sweep_conditions` | temperature x additive grid | **one fixed panel**, rescored | `optimize-conditions` with `--primers` and `--bg` |
| `pool_design_sweep.design_sweep` | conditions x oligo length | **a separate design per cell** | `plan-pool --design-grid grid.json` |

`design_sweep` is the right construction and its docstring says why: a panel
chosen under one chemistry "carries that chemistry's decisions", so rescoring it
elsewhere "measures the transplant rather than the alternative"
(`pool_design_sweep.py:5-16`). It redesigns per condition against the shared
inventory, with per-condition eligibility, and returns a non-dominated set over
size, coverage and background load rather than one winner.

---

## 2. Gap analysis

| # | Gap | Severity |
|---|---|---|
| 1 | Chemistry flags exist on `filter` and not on `optimize`/`design`: the pool is selected under one Tm model and scored under another | **Critical** |
| 2 | No test anywhere asserts that a chemistry change alters the delivered panel | **High** |
| 3 | The default chemistry cannot use the chemistry lever at all | High, decision |
| 4 | Chemistry is never a decision variable; the one redesign-per-chemistry path needs a hand-written grid on another command | High |
| 5 | The melting-temperature term is inert on the default method | High, recorded |
| 6 | 7 of 27 chemistry constants are `measured`; the temperature-domain check guards exactly one additive | High |
| 7 | Nothing is calibrated against a reaction, and no recovery target is defined | High, blocked |
| 8 | `rescore-set` accepts nine chemistry flags and reads none; two further commands carry inert temperature flags | Medium-high |
| 9 | No panel limit constrains Tm spread, although Tm tightness is the precondition for chemistry having leverage | Medium-high |
| 10 | The primer-primer duplex model is additive-blind while the primer-template model is not | Medium-high, new |
| 11 | Dimer free energy is condition-blind, so `max_dimer_dg` cannot track chemistry; the reported dimer Tm is not a temperature | Medium, new |
| 12 | `dtt_mm` and `ssb` have no Tm term and no evidence record, so they sit outside the mechanism that refuses glycerol | Medium, new |
| 13 | The Arrhenius temperature extrapolation — eight bare activation energies — is entirely outside the evidence registry | Medium |
| 14 | The capability ratchet watches no chemistry module; two further ratchets cover no chemistry key | Medium |
| 15 | On the only wet-lab-validated benchmark, the strongest single predictor enters no score | Medium, open by choice |

### Gap 1 — the chemistry a pool is selected under is not the chemistry it is scored under

This is the gap that most directly blocks the stated goal, and the project
already knows it. Additive, salt and temperature flags are declared on
`filter_parser` (`cli/pipeline.py:1415-1475`) and on `rescore_parser`
(`cli/iterate.py:1340-1349`). `optimize` and `design` carry no chemistry flag
other than `--polymerase`. The warning the code emits says the consequence
plainly (`cli/_common.py:770-776`):

> "Reaction conditions differ from the recorded '%s' step... Chemistry flags
> (`--preset`, `--betaine-m`, `--dmso-percent` and the rest) exist on `filter`
> and not on this step, so a preset applied there does not carry over. **The
> candidate pool was selected under one Tm model and is being scored under
> another.**"

`tests/test_one_reaction_one_fingerprint.py:1-11` states the asymmetry in the
same terms: "`filter` has sixteen chemistry flags and a `--preset`; `optimize`
has none."

So the mitigation is a log line. A user who wants a pool designed for 5% DMSO at
44 C must put it in params.json and must not use the flags, and nothing enforces
that: the drift is detected after the fact and the run continues. For a tool
whose purpose is chemistry-aware design, the chemistry is not expressible on the
command that performs the selection.

### Gap 2 — nothing tests that chemistry changes the panel

There are roughly 16,000 lines of chemistry-subject tests across ~68 files
(additives 14 files/2,816 lines; thermodynamics and conditions 22/5,166;
polymerase and mechanistic 15/4,008; dimers 17/3,649). Volume is not the
problem. What they assert is.

`tests/test_design_options_have_effect.py` is the only ratchet that runs the
optimizer twice and compares `result_fingerprint`, defined as the sorted primer
tuple plus the rounded `normalized_score`. Its `BEHAVIOURAL_OPTIONS` list, at
`:255-266`, holds exactly two entries: `minimize_primers` and `target_coverage`.
Neither is chemistry. Its own comment names the principle it then does not apply
to the reaction: "Stored is not consumed. Only running the thing twice can tell
the difference."

Three tests look like they close this and do not:

- `tests/integration/test_additive_ranking_e2e.py` has a docstring claiming
  "additives change the chosen primer set... must now surface a measurable
  downstream difference (different set membership OR different Tm distribution
  among selected primers)". The file asserts that `conditions is not None`, that
  `betaine_m` arrived as 1.5, a Tm delta on one primer, and that the inner
  optimizer received conditions. There is no second run at `betaine_m = 0`
  anywhere in it, and no set comparison.
- `tests/integration/test_polymerase_scenarios_e2e.py:128`,
  `test_phi29_vs_equiphi29_produce_different_primer_counts`, asserts only
  `len(phi_df) > 0 and len(equi_df) > 0` (`:171`), at the filter stage, with an
  inline note that identical counts are tolerated.
- `tests/test_network_optimizer_mechanistic.py:339`,
  `test_mechanistic_weight_affects_selection`, asserts only that both results
  have length 2, closing with "They may or may not be different".

The consequence is structural, not cosmetic: gaps 5 and 10, where a chemistry
term is wired to something that does not select, are exactly the failures this
class of test exists to catch, and the inert Tm weight (gap 5) was in fact found
by hand in September rather than by a ratchet.

**Closed, with a result that matters (2026-10-26).** A `betaine_m` case was added
to `BEHAVIOURAL_OPTIONS`, varying 0 against 1.5 M at equiphi29 42 C, and it
passes. What it found is the more important half:

| Reaction | Delivered | `normalized_score` |
|---|---|---|
| equiphi29 42 C, no betaine | 4 primers | 0.755929 |
| equiphi29 42 C, 1.5 M betaine | **the same 4 primers** | 0.695307 |

Deterministic over three runs. The reaction reaches the delivered result, and on
this pool it does not choose different oligos: the score moves by 0.061 and
membership is identical. At the phi29 30 C default the same comparison moves the
score by 0.0004, which is gap 3's regime argument appearing as a number.

This is direct confirmation of gap 5 from the selection side rather than from
reading the wiring. The scope of the claim is one synthetic 30 kb target with no
background and a four-primer panel that stops early on the dimer limit, so it
does not establish that chemistry never moves a panel; occupancy-weighted
selectivity needs a background to bite, and this harness has none. Membership is
deliberately **not** asserted, because wiring the Tm term into Stage 2 is the
change that would start moving it, and a test demanding sameness would fail on
the fix.

### Gap 3 — the default chemistry cannot use the chemistry lever

The additive lever is not a property of the additive. It is a property of
operating near Tm, and the shipped defaults do not operate there.
`parameter.py` gives `polymerase = phi29`, `min_k = 6`, `max_k = 12`,
`reaction_temp = None` (so the phi29 preset, 30 C). Against that,
[`additive_specificity.md`](../validation/additive_specificity.md) measures:

| Platform | Selectivity | With 5% DMSO + 1 M betaine | Fold |
|---|---|---|---|
| phi29, 30 C | 1.82 | 1.84 | **1.01x** |
| equiphi29, 42 C | 2.35 | 3.54 | **1.51x** |

By oligo length at equiphi29 across 38-50 C: 7-8mers flat to two decimal places
(**1.00x**), 10mers 1.36x, 11mers 1.74x, 12mers 3.03x.
[`occupancy_and_discrimination_2026-09-19.md`](../validation/occupancy_and_discrimination_2026-09-19.md)
gives the mechanism: at phi29 30 C, 86% of 40,000 random 12-mers sit at Tm >= 40 C
with occupancy 0.998-1.000 and discrimination 1.00-1.01. Far below Tm, perfect
and mismatched sites are both saturated and the occupancy terms cancel in the
ratio.

A user who sets `dmso_percent` on a default run is therefore adjusting a knob
the repository has already measured at about 1%, and nothing in the tool says
so. `filter` reports the regime and does not enforce it (Known Issue 17). The
lever lives at equiphi29 and k >= 10, which is a configuration choice rather
than a code change — and that is what makes a preset the right remedy.

### Gap 4 — chemistry is not a decision variable

`--design-grid` exists only on `plan-pool` (`cli/plan_pool.py:426`), requires a
hand-authored JSON grid naming both `lengths` and `conditions`
(`pool_design_sweep.load_design_grid:158-166`), and requires the candidate
inventory (`cli/plan_pool.py:138`). There is no chemistry axis on `optimize`, no
generated grid, and no iteration between the pool and the reaction.

The machinery for the user's goal already exists as `design_sweep` plus
`nondominated` plus `frontier`. It is gated behind a different command and a
file the user must write by hand, and nothing derives a grid from the enzyme's
own supported envelope, which the registry already knows
(`registry/views.hard_temp_ranges`).

### Gap 5 — the Tm term is inert on the default method

Measured in
[`selection_weights_are_inert_2026-09-21.md`](../validation/selection_weights_are_inert_2026-09-21.md):
all four `--application` profiles deliver an identical panel from `hybrid`, and
`tm_weight` at 0.25 and at 1.0 both change nothing. `HybridOptimizer.__init__`
constructs a `NetworkOptimizer` and passes it both weights; nothing ever reads
`self.network_optimizer` again. Stage 2 is a method on a different object and
carries no Tm term.

A Tm-tightness term does exist — `tm_tight = max(0, 1 - tm_spread / 20)` at
`base_optimizer.py:355-356` — inside `PrimerSetMetrics.normalized_score`, which
scores an **already-built** panel. Its consumers are the ensemble winner choice,
panel contraction and the report (`optimization_service.py:236,296`,
`panel_contraction.py:52`). So thermodynamic tightness influences which method's
panel wins under `ensemble`, and never which oligos a search picks. CLAUDE.md
states the position: wiring these weights "moves every delivered panel and is a
decision, not a repair."

### Gap 6 — the evidence base, and the one guarded additive

`registry/model_evidence.json` is `schema_version: 1`, compiled 2026-09-21, 27
records (25 at audit time; `tm_dtt` and `tm_ssb` added 2026-10-26, both `assumed`), with a five-level vocabulary mirrored at `model_evidence.py:53` and
ordered from strongest to weakest claim:

| Status | Meaning (verbatim, abridged) | Count |
|---|---|---|
| `measured` | "the cited work reports this value for a case the model applies it to" | 7 |
| `estimated` | "extrapolated from cited data taken at another temperature, on longer DNA, or in another buffer" | 7 |
| `empirical` | "chosen so the model behaves plausibly; no source reports it" | 1 |
| `assumed` | "a modelling decision with a stated reason and no supporting measurement" | 11 |
| `absent` | "no model computes this effect; the code must not report zero for it" | 1 |

Restricted to the twelve additive Tm coefficients: 3 measured, 5 estimated, 1
empirical, 2 assumed-absence (BSA, PEG), 1 absent (glycerol). The file's own
`provenance_note` is the most important line in it: "The primary literature was
NOT re-read while compiling this file, so a `source` records who the value is
attributed to and not that the attribution was checked."

Enforcement is correctly placed and deliberately narrow. `require_model_support`
runs at `design_request.py:663`, before any index opens, and every failure is a
hard `UnsupportedModelError` with no warning tier. It does **not** refuse
estimates, for a stated reason: "a model that declined to run on an estimate
would decline to run" (`model_evidence.py:204-209`).

The gap is the temperature domain. Only records carrying `temperature_range_c`
are checked, and **exactly one does** — `tm_dmso`, 30-45 C. The defect it was
written for is quoted in the module: a Bst design at 63 C with DMSO returned an
effective Tm computed from a 37 C coefficient and nothing said so. That defect
remains open for betaine, urea, trehalose, TMAC, formamide and ethanol: a Bst
design at 63 C with 1 M betaine is unrefused and unflagged today.

Two further specifics:

- **A second, contradicting vocabulary.** `mechanistic_params.ADDITIVE_TM_PARAMS`
  carries its own `"evidence"` key whose values (`extrapolated`, `measured`,
  `unsupported`) are not `EVIDENCE_STATUSES`, and it disagrees with the registry
  on four constants. Sharpest: urea is `measured` in the registry — one of only
  seven — on a source (`Lesnick & Bhalla 1995`) that `mechanistic_params.py`
  says "could not be located from the details given".
- **A live value drift.** The registry records propanediol at -5.3 C/M; all three
  code paths use -5.4 (`additives.py:576`, `mechanistic_params.py:143,176`).
  `tests/test_model_evidence_contract.py:103-126` pins only three quantities,
  all processivity or reach, so no additive coefficient is pinned and nothing
  catches it.

### Gap 7 — nothing is calibrated against a reaction

[`additive_specificity.md`](../validation/additive_specificity.md) states it
plainest: "none of this is validated against a measured wet-lab outcome."
[`evidence_matrix_2026-09-15.md`](../validation/evidence_matrix_2026-09-15.md)
reports seven datasets, all phi29, with "conditions available: partial" in every
row and no observation that could confirm or refute a reaction-condition
prediction. Phase 5 of the 2026-09-14 plan is blocked on data not in this
repository; Phase 4 — defining a coverage target as a fraction of reference bases
at a stated depth and sequencing budget — is an open decision.

Two consequences belong in any chemistry-facing output. Without Phase 4, a
coverage figure under a chemistry names a model output and not an experimental
outcome. Without Phase 5, a condition ranking is internally consistent and
mechanistically grounded, and is not evidence about reactions.

One positive precedent: [`reach_calibration.md`](../validation/reach_calibration.md)
calibrated the per-primer reach against Clarke 2017 Wolbachia breadth outcomes
at 3.0-6.2 kb, corroborating the 3 kb selection reach and **not** the 10 kb
reporting reach. The calibration pattern works; it has been run once, on
geometry rather than chemistry.

### Gap 8 — chemistry flags that are accepted and read by nothing

`tests/test_every_cli_option_has_an_effect.py` carries a `KNOWN_INERT`
allowlist, and **14 of its roughly 40 entries are chemistry**. The concentrated
case is `rescore-set`, which declares nine additive and salt flags and reads none
of them (`:90-99`), each excused as "additive flags are not merged on this
path":

```
betaine_m, dmso_percent, ethanol_percent, formamide_percent,
mg_conc, na_conc, tmac_m, trehalose_m, urea_m
```

`rescore-set` is the command whose job is to re-score a panel, and it is the one
place a user would naturally expect to ask "what does this panel look like under
a different buffer". Two further entries: `evaluate-set` has `reaction_temp`
with "no reader on the evaluate-set path" (`:86`), and `analyze-stability` has
`temp` with "no reader on the analyze-stability path" (`:105`).

A related finding is that one chemistry excuse is probably false and the ratchet
cannot retire it. `optimize: polymerase` is excused as "EmptyOptions;
params.json wins" (`:78`), but `tests/test_optimize_polymerase_flag_and_reach.py:52`
asserts the flag beats the params file and `cli/pipeline.py:737` reads it. The
stale-entry guard misses it because the read goes through a module-level alias
(`cli/pipeline.py:565`) and the walker indexes functions by their own name. So
the one ratchet built to find false inertness carries a false excuse
indefinitely, on a chemistry flag.

### Gap 9 — nothing constrains Tm spread

`pool_objective._LIMITS` (`:67-80`) holds six panel limits:
`min_selectivity_density`, `max_background_sites`, `max_worst_hole`,
`max_mean_gap`, `max_evenness`, `max_host_coverage`. None is thermodynamic.

This is the gap most at odds with the goal, because `condition_sweep`'s own
rationale names Tm tightness as the precondition for chemistry being able to
help: "A set whose primers melt far apart has no single temperature where all of
them sit in that window, which is why Tm tightness is a precondition for
conditions being able to help at all" (`condition_sweep.py:27-30`). The tool can
deliver a panel whose Tm spread puts the lever permanently out of reach, satisfy
every configured limit, and then report a condition sweep over it. The quantity
is already computed — `report/metrics.py:771` records `tm_range = max(tms) -
min(tms)` and `report/quality.py:324-333` rates it — so this is wiring and a
threshold decision, not new science.

### Gap 10 — the dimer model is additive-blind, the template model is not

New here, measured by probe 4. Homodimer prediction for `GCGCGCGCGCGC` at
equiphi29 42 C:

| Reaction | Dimer energy (kcal/mol) | Dimer Tm (C) |
|---|---|---|
| baseline, Na 50 / Mg 10 | -41.5093 | 225.6163 |
| Mg 0 mM | -41.5093 | 214.6062 |
| Mg 20 mM | -41.5093 | 227.2844 |
| DMSO 10% | -41.5093 | **225.6163** |
| betaine 2 M | -41.5093 | **225.6163** |
| DMSO 10% + betaine 2 M | -41.5093 | **225.6163** |

Salt moves the dimer Tm by 12.7 C across 0 to 20 mM magnesium. DMSO at 10% and
betaine at 2 M move it by nothing, to the last digit. The same two additives
move a primer-template effective Tm by 2.79 C and 1.98 C at a *fifth* of those
doses (probe 1).

So under an additive-bearing reaction the model destabilises primer-template
binding and leaves primer-primer binding at its no-additive stability. The
asymmetry biases every screen that weighs the two against each other, in the
direction of making dimers look relatively more stable than they should.
Nothing in the repository records this as a choice.

### Gap 11 — dimer free energy is condition-blind entirely

`thermodynamics.calculate_free_energy(seq, temperature=37.0)`
(`thermodynamics.py:555`) takes a sequence and a temperature: no salt, no
additive, no concentration. That is the quantity the optional `max_dimer_dg`
floor gates on, and probe 4 confirms the energy column is identical across all
six reactions including both magnesium extremes.

Severity is limited by context.
[`dimer_stability_floor_2026-09-18.md`](../validation/dimer_stability_floor_2026-09-18.md)
measured that at the shipped `max_dimer_bp = 3` a -6.0 kcal/mol floor rejects
nothing the run-length screen already rejects, and the floor ships unset. The
default gate, `lazy_dimer.dimer_screen(pool, max_dimer_bp)`, is pure run length
and carries no chemistry by design. So this bites a user who sets
`max_dimer_dg` expecting it to track their buffer.

A secondary observation from the same probe: the dimer Tm values are not
physical temperatures — 225 C for a 12-mer, -29 C for a six-base overlap. The
`forms_dimer` decision is sane in both cases (the second has a positive energy,
+0.2036, and is rejected), so the decision is not wrong; the reported number is
not a temperature and should not be printed as one.

### Gap 12 — two knobs outside the refusal mechanism

Probe 1 finds five schema knobs that move an effective Tm by exactly zero:
`glycerol_percent`, `bsa_ug_ml`, `peg_percent`, `dtt_mm`, `ssb`.

Three are handled. `tm_glycerol` is `absent`, the one status that refuses, and
the registry records the same measurement this probe reproduces. `tm_bsa` and
`tm_peg` are `assumed`, a stated decision.

`dtt_mm` and `ssb` appear nowhere in `model_evidence.json` — zero matches in the
file. Both are accepted by the schema, range-validated, carried into
`ReactionConditions`, and move no Tm. A design setting either gets silence where
a design setting glycerol gets a refusal. BSA and SSB have no physical reason to
shift a duplex Tm, so the right answer is probably a record saying so rather
than a coefficient: the gap is the missing record, not the missing number.

**Closed 2026-10-26.** `tm_dtt` and `tm_ssb` records were added (status
`assumed`, no duplex model, like BSA and PEG) and wired into
`_ADDITIVE_EVIDENCE`, so both are now examined by `require_model_support` like
every other additive. A completeness ratchet
(`test_every_reaction_additive_has_an_evidence_record`) pins that every future
`ReactionConditions` additive field carries a record or is named a buffer term,
so this cannot recur. `docs/validation/2026-10-26-suggested-additives-reach-design.md`.

### Gap 13 — the temperature extrapolation is outside the registry

Every additive coefficient in `additives.py` carries an inline citation; there
are no bare Tm coefficients in that file. The bare numbers are one module over,
in `mechanistic_params.py`, which is where `ArrheniusTmCorrector` reads its
parameters (`additives.py:148`).

Eight activation energies carry no source and no registry record: DMSO 2500
J/mol (`:52`), betaine 1800 (`:70`), formamide 3000 (`:81`), trehalose 1500
(`:97`), urea 2000 (`:114`), TMAC 1000 (`:125`), propanediol 2000 (`:145`),
ethanol 2200 (`:155`). Several are self-labelled "estimated; no multi-temp data
available". Alongside them sit bare shape constants: `gc_preference 1.3`
(`:117`), four betaine and TMAC midpoint/steepness values (`:183-187`),
`dmso_melt_rate` and `dmso_max_effect` (`:193-194`).

This is the machinery that makes every additive coefficient
temperature-dependent — precisely the coupling the request is about — and it is
the least evidenced part of the chemistry model. The 2026-09-14 audit's
objection is structural rather than numerical: "A kinetic Arrhenius law applied
to an equilibrium Tm shift needs an empirical or thermodynamic justification."

### Gap 14 — the ratchets do not watch the chemistry

`tests/test_no_capability_is_unreachable.py` has an empty `KNOWN_UNREACHABLE`
(`:87`), which looks like full coverage and is not: its `WATCHED` scope
(`:48-77`) lists sixteen modules and **none of the twelve chemistry modules**.
The scope limit is stated deliberately in the file, and the practical effect is
that a chemistry capability going orphaned is invisible to the ratchet. Nothing
is orphaned today; it is unguarded.

Two further ratchets cover no chemistry key. `test_cli_defaults_do_not_beat_params_json.py`
is an inclusion list of four pairs, none chemistry, so the `None`-sentinel
property is untested for `reaction_temp`, `polymerase` and every additive and
salt flag (the property does hold today by inspection of the declarations).
`test_params_json_routes_optional_keys.py` covers seven keys, none chemistry;
`polymerase` and `reaction_temp` appear in it only as fixture scaffolding.
`test_no_schema_key_is_inert.py` excuses no chemistry key, but asserts only that
a key binds a `parameter` global, so a chemistry key binding a global nothing
reads would pass.

### Gap 15 — the strongest known predictor enters no score

`neoswga/core/registry/INCONSISTENCIES.md` entry 7, open by deliberate choice. On the only
wet-lab-validated benchmark in the repository (six Prevotella sets, Dwivedi-Yu
2023), `max_gap` is the strongest single predictor of measured enrichment
(Spearman rho **-0.83**, and ranking by it alone reproduces the measured top
three in order), while mean binding distance has essentially none (rho -0.09).
`normalized_score` computes `max_gap` and never uses it. A coverage-hole term was
implemented, measured, changed nothing, and was not shipped.

This belongs in a chemistry audit because it bounds how much chemistry work is
worth doing now: the one quantity with a measured link to a real outcome is
geometric, not thermodynamic, and it is unused.

---

## 3. SWOT

### Strengths

- **"Unknown is not zero" is enforced, not aspirational.** `PoolObjective.coverage`
  returns `None` rather than 0.0 when no conditions were attached
  (`pool_objective.py:206-213`), `_unmet` treats an unmeasurable value as a miss
  (`:83-94`), and a limit on a background-measured quantity is refused outright
  when no background exists (`:150-158`). This is the discipline that makes a
  chemistry claim auditable at all.
- **The refusal is placed before any index opens**, and it raises rather than
  substituting: "the alternative is a number produced with no evidence behind
  it, which looks exactly like one produced with evidence"
  (`design_request.py:646-648`).
- **Both halves of the objective sit on one footing.** Occupancy-weighted
  selectivity would have rewarded a panel for reaching sites it does not occupy
  if coverage had stayed a site count; `effective_fg_coverage` unions the same
  windows weighted by occupancy, so the coverage cost of an additive is visible
  rather than hidden (-9% at equiphi29 42 C with DMSO and betaine).
- **The condition sweep refuses to search specificity alone**, on a correct
  argument: stringency buys discrimination by the same motion that pushes
  duplexes toward not being bound, so "a search on selectivity alone walks to
  the hottest, most-additive corner of any grid and recommends a reaction that
  amplifies nothing" (`condition_sweep.py:18-25`). The optimum is interior.
- **`test_condition_sweep.py` is the one chemistry test that proves a chemistry
  effect end to end.** It builds a real 60 kb foreground/background world and
  asserts that selectivity is not flat across temperature and is monotone,
  that `by_dmso[5.0] != by_dmso[0.0]`, that amplification falls as stringency
  rises, that the recommendation is not simply the hottest point, and that the
  sweep's selectivity equals the optimizer's `compute_metrics` figure to
  `rel=1e-6`. That last assertion is the model that gap 2 needs copying.
- **`design_sweep` redesigns rather than rescores**, avoiding the transplant
  fallacy, and returns a non-dominated set rather than a winner chosen by an
  arbitrary weighting (`pool_design_sweep.py:5-23,35-45`).
- **Provenance records the reaction that ran.** `effective_conditions` is built
  after CLI merge and after the GC-adaptive strategy and `retune_for_polymerase`
  have had their say, which can change the enzyme and temperature at run time
  without writing back to params.json.
- **Negative results ship as defaults.** The position-dependent mismatch model
  scored *below* a panel-size control (within-study ROC AUC 0.39 against 0.69)
  and ships off; gating candidates on occupancy made panels worse and was not
  adopted. A project that ships its own negative findings is one whose positive
  findings are worth more.

### Weaknesses

Gaps 1 to 15. Ranked by what they cost the stated goal:

1. The chemistry cannot be expressed on the command that selects, and the
   mitigation is a warning (gap 1).
2. No test relates any chemistry change to a delivered panel, which is why gaps
   5 and 10 could exist unnoticed (gap 2).
3. The lever is out of reach in the default configuration and the tool does not
   say so (gap 3); chemistry can be evaluated but not optimised from the design
   commands (gap 4).
4. Thermodynamics does not weight oligo choice on the default method (gap 5) and
   nothing keeps a panel inside a single usable temperature window (gap 9).
5. Template and dimer duplexes are treated inconsistently under additives
   (gaps 10, 11).
6. The evidence base is thin where the coupling is strongest, its two
   vocabularies disagree, and the ratchets do not watch it (gaps 6, 13, 14).

### Opportunities

Ordered by benefit per unit of risk. The first five need no new science.

- **Add one chemistry case to `BEHAVIOURAL_OPTIONS`** (gap 2). Two runs differing
  only in `betaine_m`, compared on `result_fingerprint`. It is a handful of
  lines in a file that already has the fixtures, and it converts the central
  question of this audit from an opinion into a pinned fact — in whichever
  direction it comes out. If the panel does not move, that is the most important
  single finding available about chemistry-aware design in this tool, and it is
  currently unknown.
- **Put the chemistry flags on `optimize` and `design`** (gap 1), or refuse the
  run instead of warning when the recorded and current reactions differ. Either
  removes a silent two-chemistry design; the second is smaller and is the safer
  first move.
- **A `max_tm_spread` panel limit** (gap 9), as a seventh entry in `_LIMITS`, fed
  by the `tm_range` already computed at `report/metrics.py:771`, unset by default
  like every other panel limit, repaired through the existing `repair_panel`. It
  encodes the precondition `condition_sweep` names.
- **Apply the additive correction to the primer-primer duplex** (gap 10). The
  correction already exists and is already called for the template duplex; this
  removes an asymmetry rather than adding a model.
- **Pin the additive coefficients in the contract test** (gap 6). The test
  already pins three processivity and reach constants; extending it to the twelve
  `tm_*` values would have caught the propanediol drift.
- **Extend `temperature_range_c` to every additive record** (gap 6). One record
  carries a domain and the refusal machinery is written; the work is reading each
  source for the range it covers. This is the highest-value evidence task,
  because it converts a silent extrapolation into a refusal exactly where users
  will cross it when they move from phi29 to Bst.
- **Generate the design grid from the enzyme envelope** (gap 4). `plan-pool
  --design-grid` already does the hard part; a grid derived from the polymerase's
  `hard_temp_ranges` plus the additive doses whose evidence covers that range
  turns a hand-written file into a default.
- **Ship an equiphi29 discrimination preset** (gap 3): equiphi29, 42-45 C,
  k >= 10, with the measured 1.51x and its -9% coverage cost stated in the
  preset's own description. It converts a buried measurement into the default
  answer for users who want the lever, without changing anyone's default.
- **Reconcile the two evidence vocabularies** and bring the eight Arrhenius
  activation energies into the registry even if every one enters as `assumed`
  (gaps 6, 13). A recorded assumption is auditable; a bare module constant is
  not.
- **Calibrate what the data can carry.** The 310 qPCR fold measurements in
  `tests/validation/data/dwivedi_yu_2023_plasmid_amplification.json` cannot
  calibrate coverage or selectivity, since site geometry is constant by
  construction, and can calibrate a sequence-level amplification term.

### Threats

- **Overclaiming.** A command reading "optimise conditions" invites the inference
  that the recipe is validated. With 7 of 25 constants `measured`, the primary
  literature unverified by the registry's own note, and nothing calibrated
  against a reaction, a chemistry recommendation presented as an optimum is the
  project's largest reputational exposure. The existing documents are scrupulous
  about this; a new feature is where the discipline would slip.
- **The temperature dependence may be wrong in form, not just in value.** Applying
  a kinetic Arrhenius law to an equilibrium Tm shift is unjustified (gap 13). If
  the form is wrong, every temperature-dependent additive prediction is wrong in
  a way coefficient tuning cannot fix, and the Bst and klenow ends of the range
  are where it would show.
- **The magnesium model runs outside its characterised range.** The registry
  records that the sodium-equivalence was characterised over PCR magnesium
  concentrations while SWGA runs 10 mM, and that "the code applies total, not
  free" magnesium with no dNTP chelation term. Probe 2 shows magnesium is the
  dominant salt term, moving mean effective Tm 9.13 C from 0 to 10 mM, so a bias
  here propagates to every Tm gate and every occupancy weight.
- **Misallocated effort.** The one quantity with a measured link to a real outcome
  is `max_gap` at rho -0.83, and it is unused (gap 15). Chemistry work that
  displaces it is work spent on the less evidenced of two available levers.
- **Silent divergence between the two vocabularies.** Propanediol has already
  drifted -5.3 against -5.4 with no test to catch it. Two registries describing
  the same constants will keep drifting.
- **Tests that read as coverage they do not provide.** Three named in gap 2. A
  future reader who greps for "additives change the chosen primer set" finds a
  docstring that says the work was done.

---

## 4. What to do, sequenced

**First, and before anything else is decided.**

1. ~~Add a chemistry case to `BEHAVIOURAL_OPTIONS` and run it (gap 2).~~ **Done
   2026-10-26.** The panel's score moves and its membership does not, on the
   synthetic harness. That answer sizes everything below: chemistry reaches the
   result but does not steer selection, so item 11 stops being optional. Still
   worth doing in the same area: correct the three misleading test names, and
   repeat the comparison on a real pool with a background, where
   occupancy-weighted selectivity can actually bite.

**Then, no new science, each pinned by a test.**

2. Refuse rather than warn on condition drift between `filter` and `optimize`,
   or carry the chemistry flags onto the design commands (gap 1).
3. Record `dtt_mm` and `ssb` in `model_evidence.json` with the status their
   evidence supports (gap 12). Smallest change here, and it closes the only case
   where a no-Tm-term knob escapes the refusal mechanism.
4. Pin the twelve additive Tm coefficients in the contract test and resolve
   propanediol to one value (gap 6).
5. Apply the additive correction to the primer-primer duplex, and stop printing a
   dimer Tm that is not a temperature (gaps 10, 11).
6. Add `max_tm_spread` to the panel limits, unset by default (gap 9).
7. Merge the chemistry flags on `rescore-set` or remove them (gap 8), and fix the
   walker's alias blindness so the `optimize: polymerase` excuse can retire.
8. Extend the capability ratchet's `WATCHED` scope to the chemistry modules
   (gap 14).

**Next, a measurement each.**

9. Extend `temperature_range_c` across the additive records, then measure how
   many shipped example configurations the new refusals would reject (gap 6). The
   2026-09-14 precedent applies: label and measure the blast radius before
   changing a delivered design.
10. Derive a default design grid from the enzyme envelope and measure, on a
    prepared example, whether the non-dominated set differs from the single
    baseline design (gap 4). If it does not, chemistry is not the binding
    constraint on that pool and the output should say so.

**Decisions, not branches.**

11. Whether to wire the Tm term into Stage 2 (gap 5). It moves every delivered
    panel, so it is a decision with a measurement attached, not a repair.
12. Whether the shipped default should remain phi29 at 30 C given the lever is
    measured at 1.01x there (gap 3). A preset is the low-risk form; changing the
    default is not.
13. Phase 4 of the 2026-09-14 plan: define the coverage target as a fraction of
    reference bases at a stated depth and sequencing budget (gap 7).

**Blocked on data.**

14. Phase 5. Independent SWGA reactions spanning real pools and relevant
    conditions, held out as whole reactions. Until then the honest position is
    the one the repository already takes: the tool can express many chemistries
    and rank them, and its numbers do not establish an optimal recipe or a
    fraction of genome recovered in a laboratory.

---

## 5. What this audit did not establish

- **Whether a chemistry change moves the delivered panel on a real pool.**
  Partly settled after this report was first written. On the synthetic harness in
  `test_design_options_have_effect.py` the score moves and membership does not
  (gap 2). That harness has no background genome, so occupancy-weighted
  selectivity cannot contribute, and the panel is four primers stopped early by
  the dimer limit. The same comparison on a prepared example with a background
  is the measurement that would settle it, and it needs `step3_df.csv`, the
  position indexes and a reference genome. None is in a fresh clone. The
  adjacent question — whether chemistry changes candidate *eligibility* — is
  settled affirmatively by probes 1 and 2, since the Tm gate moves with every
  additive and salt term measured.
- **A clean full-suite run.** The dependencies were installed partway through
  this work, after which `pytest tests/ -n 8` reports 6,454 passed, 282 skipped,
  55 failed and 4 errors. Every failure is
  environmental and pre-existing: no k-mer counter is on PATH, so
  `examples/plasmid_example/step3_df.csv` is never built and the fixtures that
  need it error rather than skip, and no reference genome is present. Verified by
  running the suite with this report's changes stashed, which gives 6,451 passed
  and 56 failed, so the chemistry case added for gap 2 contributes two passing
  tests and no failures. At least one of the environmental failures is flaky
  across runs. `test_design_options_have_effect.py` and the five ratchet files
  pass in full.
- **Whether any additive coefficient is numerically right.** Unchanged from
  2026-09-14 and out of scope: this audit reports evidence status, not new
  values. That plan's rule still holds, that replacing one unsourced number with
  another repeats the error being corrected.
