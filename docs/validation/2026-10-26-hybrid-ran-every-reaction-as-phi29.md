# `hybrid` ran every reaction as phi29 at 30 C

Found and fixed 26 October 2026, on `main` at `7fa718f`, by following the
dominating-set dimer-temperature fix
([`2026-10-26-dominating-set-dimer-temperature.md`](2026-10-26-dominating-set-dimer-temperature.md))
back to where the reaction is supposed to reach an optimizer. Reproduce with:

```
PYTHONPATH=. python docs/validation/hybrid_reaction_2026-10-26/probe.py
```

The regressions are `tests/test_optimizer_config_reaches_optimizers.py`
(`test_the_configured_reaction_reaches_the_hybrid_inner`,
`test_an_explicit_polymerase_still_beats_the_attached_conditions`).

## The defect

`unified_optimizer` hands every optimizer its `ReactionConditions` and hands it
**neither `polymerase` nor `reaction_temp`** -- the module contains no mention
of either name. Both hybrid wrappers then built the inner with:

```python
polymerase=kwargs.get("polymerase", "phi29"),   # a REAL default
reaction_temp=kwargs.get("reaction_temp"),       # always None
```

So the inner always built a **phi29** `PolymeraseConfig`, and
`resolve_tm_inputs` fell back to that preset's 30 C. Measured through the real
factory dispatch, before the fix:

| configured | inner.polymerase | inner.reaction_temp | thermo_filter |
|---|---|---|---|
| phi29 | phi29 | 30.0 | False |
| equiphi29 (45 C) | **phi29** | **30.0** | **False** |
| bst (65 C) | **phi29** | **30.0** | **False** |

Identical on `hybrid` and on `background-aware`, which build the same inner from
two separate argument lists.

This is Known Issue 8's most-repeated shape -- a real default beating the
configured value -- and the fourth variant recorded under it, the wrapper/inner
seam. It is also the *second* defect found at this exact seam: the comment
beside the fix records an earlier one where "six values arrived at `hybrid` and
were dropped here", costing two thirds of a user-approved equiphi29 pool.

## What rode on it

**1. The Stage-0 thermodynamic pre-screen never ran for the polymerases that ask
for it.** `optimize` gates that whole pre-stage on `poly_config.thermo_filter`:

```python
if self.poly_config.thermo_filter and candidates_filtered:
    candidates_filtered = self._thermo_filter_with_cache(...)
```

and the presets disagree on that flag:

| polymerase | preset reaction_temp | thermo_filter |
|---|---|---|
| phi29 | 30.0 | False |
| equiphi29 | 42.0 | **True** |
| bst | 63.0 | **True** |
| klenow | 37.0 | False |

Resolved as phi29, the flag is False, so for an equiphi29 or bst design the
screen that rejects candidates on hairpin, homodimer and heterodimer energies
and on the Tm window was skipped entirely.

Measured under the configuration of the **shipped**
`examples/equiphi29_scenario/params.json` (equiphi29, `reaction_temp` 43 C, Tm
window 40-65, Na 50 mM, Mg 10 mM, DMSO 5%, betaine 1.0 M), over 500 synthetic
candidates in that example's `long_primer_mode` k range:

- pre-screen allowed to run: keeps **299/500**
- resolved as phi29, pre-screen skipped: **500/500** reach selection

So 201 of 500 candidates -- 40% -- entered selection although the configured
reaction rejects them. A shipped example was in the affected configuration.

**2. Both Stage-2 temperatures were wrong.** Stage 2's dimer screen reads
`optimizer.reaction_temp` (`refine_hybrid_stage2`), and the Stage-2 Tm term
ranks on effective Tm against `reaction_temp + 5`. Both ran at 30 C against a
63 C or 42 C reaction. The dimer half is the same permissive-direction error
quantified in the dominating-set record.

## The fix

One door, in the inner that both wrappers construct, mirroring the fallback
`network_optimizer` already carried for `reaction_temp` (added there with the
note that its Tm score had "previously been a silent no-op"):

```python
def resolve_reaction(polymerase, reaction_temp, conditions):
    polymerase = polymerase or getattr(conditions, "polymerase", None) or "phi29"
    if reaction_temp is None:
        reaction_temp = getattr(conditions, "temp", None)
    return polymerase, reaction_temp
```

`stage2_ranking.resolve_reaction`, called once in `HybridOptimizer.__init__`.
The inner's `polymerase` default became a `None` sentinel and both wrappers now
pass `kwargs.get("polymerase")` without a literal default, so an absent value no
longer beats the file -- the remedy this project's conventions already name for
CLI flags, applied to a constructor argument.

Precedence is **explicit > conditions > preset**. A caller that states a
polymerase still gets it even when the attached conditions disagree, pinned by
`test_an_explicit_polymerase_still_beats_the_attached_conditions`; without that
test, fixing the default would merely have inverted which value was ignored. A
`reaction_temp` left None still falls through to the preset in
`resolve_tm_inputs` -- now the preset of the polymerase actually configured.

## Scope, and what is unchanged

- **phi29 is unchanged, and it is the default.** Its conditions resolve to the
  same preset and the same 30 C the bug produced, so no default-path run moves.
  phi29 is carried as the control row in the parameterised test for that reason.
- **equiphi29, bst and klenow designs change**, which is the point: they begin
  running the screen and the temperatures they configured. Their delivered panels
  will differ from previous releases. This is a correction, not a new capability,
  so it is on rather than behind a key -- "defaults follow measurement" governs
  capabilities with no demonstrated benefit, not a wrong polymerase.
- Full suite before and after: the same 60 failures, all environmental in this
  container (no k-mer counter, no `pysam`). Zero added, zero removed.

## Scope of the claim

The 299/500 figure is **synthetic candidates**, not a delivered panel: no
prepared shortlist exists in this container, since generating one needs a k-mer
counter. It bounds what the pre-screen rejects under a real shipped
configuration; it does not say how many primers a finished equiphi29 design
gains or loses, which needs a run against a genome and is not determined here.
