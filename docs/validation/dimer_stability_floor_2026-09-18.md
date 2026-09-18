# A stability floor for the dimer screen, and what measurement decided

Item 4 of [getting_ahead_on_spacing_2026-09-18.md](getting_ahead_on_spacing_2026-09-18.md).
Measured 2026-09-18 before anything was built, which changed the design twice
and corrected two claims in the audit that proposed it.

Reproduce with `scripts/benchmarking/dimer_policy_comparison.py`, which needs no
design directory and no external tool.

## What the audit claimed, and what is actually true

The audit said `dimer.is_dimer_thermodynamic` was condition-aware and that
COATswga was "genuinely ahead here". Both are overstated.

**It is temperature-aware, not condition-aware.** It accepts a
`ReactionConditions` and reads only its `temp`, because
`thermodynamics.calculate_free_energy(seq, temperature)` takes no salt and no
additive term. So wiring it buys a screen that moves with the reaction
TEMPERATURE. That is still more than any published SWGA tool has, since all
three call the same melting routine with no arguments at all, but it is not the
full chemistry the Tm path applies.

**NeoSWGA's default screen is already stricter than COATswga's.** COATswga
rejects a pair scoring below -2.79 by PrimerROC and applies no length cap.
Measured on the same pools, a -2.79 free-energy floor with no length cap admits
complementary runs of 5 to 6 bp, while NeoSWGA's default `max_dimer_bp` of 3
admits 3. The models are not identical, so read the direction rather than the
digits; the direction is not in COATswga's favour.

## Cost is not the obstacle

| Screen | Per pair |
|---|---|
| longest complementary run | 10.8 us |
| free energy of that run | 11.5 us |

A ratio of 1.1. The audit said this "is on the O(n^2) hot path, so it needs the
lazy screen and measurement before it becomes a default". The measurement says
cost is irrelevant, and the reason is that the run screen is not cheap either:
both walk the same dynamic-programming table.

## At the shipped default, a floor decides nothing

Every pair a -6 kcal/mol floor rejects, `max_dimer_bp` 3 already rejects, at
35%, 50% and 65% GC alike. So the floor cannot be sold as catching what the run
screen misses. It is not an accuracy improvement at the defaults.

## What it is for: the run screen is what bounds panel size

200 primers, 19,900 pairs, greedy largest compatible subset.

| Policy | Pairs rejected, 35% GC | Panel | 50% GC | Panel | 65% GC | Panel |
|---|---|---|---|---|---|---|
| run <= 3 (default) | 29.1% | 20 | 22.9% | 17 | 28.9% | 21 |
| run <= 4 | 7.2% | 42 | 4.9% | 51 | 7.2% | 45 |
| run <= 5 | 1.6% | 93 | 0.9% | 110 | 1.7% | 95 |
| run <= 5 and dG > -6 | 1.7% | 91 | 1.1% | 104 | 3.5% | 83 |
| run <= 5 and dG > -4 | 2.5% | 77 | 4.6% | 70 | 13.7% | 50 |
| dG > -6 alone | 0.2% | 172 | 0.7% | 139 | 3.4% | 91 |

The default supports a panel of 17 to 21 out of 200. That is the constraint
CLAUDE.md records as usually binding: the shipped pools support 29, 31 and 26
primers at `max_dimer_bp` 3 against panels of 200, 160 and 36, and
`--allow-dimer-relaxation` is the only escape hatch. It is blunt, admitting
violating pairs with a warning.

`run <= 5 and dG > -4` supports 50 to 77 primers with a hard 5 bp length cap
still in force. That is the case for the floor: **a principled way to raise
`max_dimer_bp` rather than a way to replace it.**

One row is worth noticing on its own. At 65% GC a -6 floor rejects twice as many
pairs as `run <= 5` (3.4% against 1.7%) and supports a panel of 91 against 95,
so the pairs it rejects are spread across more primers rather than concentrated
in a few hubs. The two screens are not simply stricter and looser versions of
each other.

## Why a floor is never offered alone

A floor bounds STABILITY, not LENGTH, so it can admit a long AT-rich run.

| Policy | Worst run admitted, 35% GC | 50% GC | 65% GC | Pairs at 8 bp or more |
|---|---|---|---|---|
| run <= 3 (default) | 3 | 3 | 3 | 0 |
| run <= 5 | 5 | 5 | 5 | 0 |
| dG > -6 alone | **8** | 7 | 7 | 3 at 35% GC |
| dG > -4 alone | 7 | 6 | 5 | 0 |
| dG > -2.79 alone | 6 | 6 | 5 | 0 |

A -6 floor with no length cap admits an 8 bp complementary run. This project
delivered an 11 bp heterodimer against a configured 3 once already, and
`max_dimer_bp` is what stops that.

So `max_dimer_dg` is applied only to pairs the length screen has already
passed. It can make the screen stricter and never looser, and
`TestTheFloorOnlyEverAdds` pins that a lenient floor cannot readmit a long run.

## What shipped

- `max_dimer_dg` in params.json, unset by default, refused above 0.
- The screen consults it only after `max_dimer_bp` passes, at the reaction
  temperature, and forwards it from all seven construction sites. A test walks
  the source for a site that builds a screen without it.
- A configured floor forces the pairwise screen, because the dense matrix codes
  t-mers and cannot express free energy, so using it would silently ignore the
  floor.
- `is_dimer_thermodynamic` gained a plain `temperature` argument.
  `ReactionConditions` validates its temperature against the polymerase and
  refuses 63 C under the phi29 default, so a screen that built one to pass a
  temperature through would raise on a bst design.

## What this does not establish

- -6.0 kcal/mol follows Rychlik (1995) Mol Biotechnol 3:129-134. Nothing here
  validates it, or any other threshold, against a reaction. No dataset in this
  repository carries per-panel dimer outcomes.
- Random GC-matched pools, not a genome's k-mers, and one length (k = 12). The
  panel figures are a greedy largest-compatible-subset, so they are lower
  bounds.
- The comparison with COATswga's PrimerROC threshold sets two different models
  against each other on one scale. The direction is clear; the digits are not
  transferable.
- No delivered panel was re-derived under a floor. The claim is that the
  capability exists and is safe by construction, not that it improves a design.
- The floor is temperature-aware only. Salt and additives move Tm in this
  codebase and do not move this delta-G.
