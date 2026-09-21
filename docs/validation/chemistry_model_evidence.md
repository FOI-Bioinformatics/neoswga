# What the chemistry model supports, and what it only appears to

Compiled 21 September 2026, for Task 4 of the valid-design plan. The
machine-readable form is `neoswga/core/registry/model_evidence.json`, which
`neoswga/core/model_evidence.py` loads and
`tests/test_model_evidence_contract.py` checks against the shipped constants.

## What this document is not

The plan asks for an audit of the original papers. **That did not happen here.**
Every `source` in the registry is the attribution this repository already
carried, in `docs/SCIENCE_CITATIONS.md`, `registry/INCONSISTENCIES.md` or the
module comments. The primary literature was not re-read while compiling it.

Saying so is the point rather than a caveat on it. The rule the registry exists
to enforce is that an assumption is not promoted to a measurement because it has
a citation, and a registry that quietly implied verification would break that
rule in the act of recording it. What the registry does add is the second
judgement, which the prose ledger did not make: whether the cited work covers
**this** case, meaning short oligos, this buffer and this temperature.

A real audit would re-read each source and record the exact table or equation
and its experimental domain. Until then `status` is a statement about fit, not
about having checked the arithmetic.

## The five statuses

| Status | Meaning |
|---|---|
| `measured` | the cited work reports this value for a case the model applies it to |
| `estimated` | extrapolated from cited data at another temperature, on longer DNA, or in another buffer |
| `empirical` | chosen so the model behaves plausibly; no source reports it |
| `assumed` | a modelling decision with a stated reason and no supporting measurement |
| `absent` | no model computes this effect, and the code must not report zero for it |

The line that carries the weight is between the first two. Most additive
coefficients are 37 °C figures for PCR-length duplexes, applied to 12-mers in a
30 °C isothermal reaction. That extrapolation may well be fine. It is not a
measurement of it, and only `tm_urea` was selected specifically because its
source concerns short oligos.

## What is refused

`require_model_support(request)` runs inside `resolve_design_request`, so it
fires before any index is opened. It refuses three things and deliberately not
a fourth.

**An unknown polymerase.** No parameter set, so nothing downstream has evidence.

**An oligo length outside the enzyme's modelled range.** Bst is modelled for
15–25 nt. A 6-mer under Bst at 63 °C is not wrong in a way anyone can see: it is
a number with no evidence behind it. This is the recorded defect in which a Bst
design was filtered through phi29's 6–12 bp window.

**An additive whose duplex effect nothing computes, where the literature expects
one.** Currently glycerol alone.

**It does not refuse estimates.** A model that declined to run on an
extrapolated coefficient would decline to run at all, since that is most of the
additive chemistry in this field. The registry records the strength; it does not
gate on it.

## Two findings from compiling it

**Glycerol is accepted, validated and has no melting-temperature term.**
Measured here on 21 September 2026: a 12-mer at 10 % glycerol returns the same
effective Tm as at 0 %, to the last digit, while the value is range-validated to
0–15 %, carried through `to_dict`, and printed in the condition summary. The
mechanical reason it is accepted at all is that the mechanistic model does carry
glycerol enzyme-stability and speed terms. The literature expects a real
destabilisation of duplex DNA, so zero here is an absent model rather than a
result, and the shipped `q_solution` preset sets 10 %.

No coefficient is invented to close this. Inventing one is exactly the promotion
of an assumption the registry exists to prevent. A design that sets glycerol is
refused, naming the missing model.

Measured for comparison, on the same 12-mer, as ΔTm against no additive:

| Additive | ΔTm (°C) |
|---|---|
| DMSO 5 % | −2.689 |
| betaine 1.0 M | −1.565 |
| trehalose 0.5 M | −1.480 |
| formamide 5 % | −3.164 |
| ethanol 2 % | −0.784 |
| urea 1.0 M | −2.581 |
| TMAC 0.1 M | −0.116 |
| propanediol 0.5 M | −2.652 |
| glycerol 10 % | 0.000 |
| BSA 200 µg/mL | 0.000 |
| PEG 4 % | 0.000 |

BSA and PEG are recorded as `assumed` rather than `absent`: they act on the
enzyme and on crowding, not on duplex stability, so no Tm term is a decision
with a reason. Glycerol is the one where the absence is a gap.

**The registry package was not installed.** `neoswga.core.registry` was missing
from `pyproject.toml`'s explicit `packages` list, so it was absent from a built
wheel entirely. Verified by building one: the only entries matching "registry"
were two unrelated modules. `core/parameter.py` imports
`neoswga.core.registry.views` at module scope, so an installed neoswga could not
import its own configuration module.

Package data could not have rescued this, since `include-package-data` applies
to packages that are being installed and this one was not.
`tests/test_model_evidence_contract.py` now compares the declared package list
against the packages on disk, so a new subpackage cannot be added without being
installable.

## The drift the registry replaces

`docs/SCIENCE_CITATIONS.md` states Klenow processivity as 10,000 bp, citing
Bambara et al. (1978). The shipped registry says 40 bp, Klenow being
distributive. The prose was correct when written and the code moved.

That is the argument for a machine-readable registry over a document: this one
is checked against `registry/views.as_characteristics()` by a test, so the two
cannot disagree silently. `SCIENCE_CITATIONS.md` remains useful as the long-form
discussion and should be read as commentary rather than as the record.

## Known gaps not closed here

- The coefficients exist in five places (`mechanistic_params`, the legacy
  methods in `additives.py`, the non-Arrhenius branch of `reaction_conditions`,
  `mechanistic_model`, and a knowingly wrong set in `thermodynamics.py`).
  Nothing asserts they agree, and this registry records one of them.
- The DMSO coefficient of −0.55 °C/% sits below its own cited range of
  0.6–0.75, so DMSO looks gentler on Tm than the source says.
- The magnesium salt equivalence is applied to **total** magnesium. Free
  magnesium is lower, since dNTPs chelate it roughly 1:1, so 10 mM total with
  1.6 mM dNTPs leaves nearer 8.4 mM free.
- Neither coverage reach has been measured. `reach_phi29` is a design-density
  convention taken from sets with wet-lab success, recorded as `assumed`.
  `calibrate-reach --bam` exists and there is no BAM in this repository.
- The mismatch penalty is uniform in identity and position. Any specificity
  claim resting on it is an exact-match claim with a uniform correction.

## Concentration policy: measured, and not yet propagated

`DesignRequest.concentrations_molar` implements both declared allocation modes.
It is **not** propagated into per-candidate evaluation, and that is a decision
rather than an omission.

Under a fixed total of 4 µM, the per-oligo concentration falls with panel size
and the melting temperature falls with it. Measured on a 12-mer, 21 September
2026:

| Panel | Per-oligo (M) | Tm (°C) |
|---|---|---|
| 1 | 4.000e-06 | 59.362 |
| 6 | 6.667e-07 | 55.191 |
| 24 | 1.667e-07 | 52.035 |
| 96 | 4.167e-08 | 48.940 |

Ten degrees across the range, which is large. The quantity selection actually
uses is occupancy, and occupancy barely moves, for the reason recorded in Known
Issue 17: at phi29 30 °C almost everything is saturated.

| Reaction | Occupancy at n=6 | at n=96 | Ratio |
|---|---|---|---|
| phi29 30 °C | 0.999993 | 0.999888 | 1.000 |
| equiphi29 42 °C | 0.997476 | 0.961089 | 1.038 |
| phi29 30 °C, DMSO 10 % + betaine 1.5 M | 0.999781 | 0.996007 | 1.004 |

So a fixed-total policy changes the occupancy-weighted objective by under 4 % at
the warmest supported reaction and by nothing measurable at the default one.
Threading concentration through every evaluation and cache key is a large change
for that, and this project's standard is that a change without a demonstrated
benefit does not ship on by default. The policy is therefore recorded on the
request, reported, and available to a caller, while selection continues to use
the configured `primer_conc`.

What this does **not** license is quoting a Tm from a fixed-total design without
its panel size. That figure moves by ten degrees and the report should say which
size it was computed at.
