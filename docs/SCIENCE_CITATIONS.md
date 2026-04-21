# Scientific constants and their citations

Every numeric scientific constant used by neoswga, its primary literature
source (where one exists), and any deviation from the published value.
Values flagged `EMPIRICAL` are model tuning parameters rather than
measured constants; they are calibrated to match qualitative behaviour
observed in practitioner protocols and are documented explicitly so
users know which numbers to treat as ground truth vs as heuristics.

## Nearest-neighbor thermodynamics (SantaLucia 1998)

`neoswga/core/thermodynamics.py:34-80` — all 10 canonical
Watson-Crick NN stacks match SantaLucia (1998) PNAS 95:1460-1465
Table 1 exactly:

| Pair | dH (kcal/mol) | dS (cal/mol/K) |
|------|---------------|----------------|
| AA/TT | -7.9 | -22.2 |
| AT/TA | -7.2 | -20.4 |
| TA/AT | -7.2 | -21.3 |
| CA/GT | -8.5 | -22.7 |
| GT/CA | -8.4 | -22.4 |
| CT/GA | -7.8 | -21.0 |
| GA/CT | -8.2 | -22.2 |
| CG/GC | -10.6 | -27.2 |
| GC/CG | -9.8 | -24.4 |
| GG/CC | -8.0 | -19.9 |

Initiation, symmetry, and terminal-AT corrections also match
Table 1 exactly (`thermodynamics.py:88-96`).

## Salt correction

- SantaLucia 1998 additive form: `Tm(salt) = Tm(1M) + 12.5 * log10([Na+])` (`thermodynamics.py:~305`).
- Owczarzy 2004 entropy-based: `dS_corrected = dS + 0.368 * N * ln([Na+])` — the 0.368 coefficient matches Owczarzy et al. (2004) Biochemistry 43:3537-3554 exactly (`thermodynamics.py:~318`).
- Owczarzy 2008 Mg2+ equivalence: `[Na+]_equiv = [Na+] + 3.3 * sqrt([Mg2+])` — matches Owczarzy et al. (2008) Biochemistry 47:5336-5353 (`thermodynamics.py:~328`).

## Additive Tm corrections (at 37 C reference)

Values below are implemented both in `neoswga/core/additives.py` and
mirrored in `mechanistic_params.ADDITIVE_TM_PARAMS`.

| Additive | Coefficient | Unit | Citation | Status |
|---|---|---|---|---|
| DMSO | -0.55 | C per % | Chester & Marshak (1993); Varadaraj & Skinner (1994) | literature midpoint |
| Betaine (uniform component) | -1.2 | C per M | Rees et al. (1993), Henke et al. (1997) | consensus of two studies |
| Betaine (full GC equalization) | 5.2 | M | Rees et al. (1993) | primary literature |
| Trehalose | -3.0 | C per M | Spiess et al. (2004) | midpoint of 2-4 C/M range they reported |
| Formamide | -0.65 | C per % | Blake & Delcourt (1996); McConaughy (1969) | primary literature |
| Ethanol | -0.4 | C per % | Cheng et al. (1994) | primary literature |
| Urea | -2.5 | C per M | Lesnick & Bhalla (1995) for short oligos | selected over Hutton (1977) -5 C/M which measured long DNA |
| TMAC (uniform) | -0.5 | C per M | Melchior & von Hippel (1973) | minimal uniform component; GC-dependent effect is primary |

## Arrhenius activation energies for temperature-dependent corrections

`neoswga/core/mechanistic_params.py:ADDITIVE_TM_PARAMS` — activation
energies (J/mol) used to extrapolate the 37 C coefficients above to
other reaction temperatures (30 C phi29, 42-45 C equiphi29).

| Additive | Ea (J/mol) | Status |
|---|---|---|
| DMSO | 2500 | ESTIMATED from Chester & Marshak (1993) multi-temp data |
| Betaine | 1800 | ESTIMATED from Rees et al. (1993) multi-temp data — verify against lab data before absolute predictions |
| Formamide | 3000 | ESTIMATED from McConaughy (1969) |
| Trehalose | 1500 | EMPIRICAL — no multi-temperature data in the literature for this compound on short oligos |
| Urea | 2000 | ESTIMATED from Hutton (1977) multi-temp data |
| TMAC | 1000 | EMPIRICAL — no direct literature source; conservative low value |
| Ethanol | 2200 | ESTIMATED |

Consequence: Tm corrections AT 37 C are well-validated; corrections at
the SWGA operating temperatures (30 C phi29, 42-45 C equiphi29) carry
additional uncertainty from the Arrhenius extrapolation. This uncertainty
is small enough that relative ranking between primer candidates stays
meaningful, but absolute Tm predictions can drift 1-2 C from measured
values at the temperature extremes.

## Polymerase characteristics

`neoswga/core/reaction_conditions.py:POLYMERASE_CHARACTERISTICS`
(lines 52-97).

| Polymerase | Optimal T (C) | Processivity (bp) | Citation |
|---|---|---|---|
| phi29 | 30 | 70,000 | Blanco et al. (1989) JBC 264:8935 |
| equiphi29 | 42 | 80,000 | NEB technical data |
| bst | 63 | 2,000 | Notomi et al. (2000) NAR 28:e63 |
| klenow | 37 | 10,000 | Bambara et al. (1978) JBC 253:413 |

Mg2+ recommendations by GC content
(`reaction_conditions.py:~740`): Rahman et al. (2014) PLoS One
9:e112515.

### Typical amplicon length vs. processivity

The processivity values above are single-molecule theoretical maxima
measured under optimal lab conditions (infinite time, no competition).
In a real MDA / SWGA reaction, the mean fragment length is substantially
shorter because of reaction-time limits, primer competition, strand-
displacement kinetics, and secondary-structure stalling. Practitioner
reports (Leichty & Brisson 2014; Clarke et al. 2017; Cowell et al. 2017)
consistently cluster amplicon size distributions in the 1-10 kb range,
with medians around 2-5 kb for phi29.

NeoSWGA exposes both:

| Polymerase | `processivity` (bp) | `typical_amplicon_length` (bp) | Used for |
|---|---|---|---|
| phi29 | 70,000 | 3,000 (±1 kb) | Realistic coverage metric |
| equiphi29 | 80,000 | 4,000 (±1 kb) | Realistic coverage metric |
| bst | 2,000 | 1,000 (±0.5 kb) | Realistic coverage metric |
| klenow | 10,000 | 1,500 (±0.5 kb) | Realistic coverage metric |

`reaction_conditions.POLYMERASE_CHARACTERISTICS` carries both fields.
Helpers:

- `get_polymerase_processivity(name)` → theoretical maximum (legacy).
- `get_typical_amplicon_length(name)` → practitioner-observed mean
  fragment length (new in Phase 16).

`coverage.compute_per_prefix_coverage` defaults to the realistic
amplicon length (3 kb for phi29) so user-facing `fg_coverage` /
`per_target_coverage` numbers reflect what a lab will actually see.
The amplicon-network reachability path still uses processivity because
the question there is "can two primers connect in principle?" not
"how much DNA will I see?".

Pass `coverage_metric='processivity'` to `polymerase_extension_reach`
if you specifically want the theoretical upper bound.

## Mechanistic-model empirical tuning

These constants drive the four-pathway mechanistic model (Tm, accessibility,
enzyme activity, kinetics) and are **not** primary-literature values.
They are calibrated so model predictions reproduce the qualitative
behaviour reported in SWGA protocols. For absolute-yield prediction,
calibrate against your own wet-lab measurements via
`ExperimentalTracker` and `predict-efficiency --track` (see
[`active-learning-guide.md`](active-learning-guide.md)).

| Parameter | File:Line | Status |
|---|---|---|
| phi29 DMSO threshold 5% | `mechanistic_params.py:182` | EMPIRICAL — operating-point chosen from vendor protocol guidance |
| phi29 dmso_mild_coef 0.02 | `mechanistic_params.py:184` | EMPIRICAL — tuned to model |
| phi29 dmso_steep_coef 0.12 | `mechanistic_params.py:185` | EMPIRICAL — tuned to model |
| equiphi29 DMSO threshold 4% | `mechanistic_params.py:192` | EMPIRICAL — slightly tighter than phi29 per EquiPhi29 vendor notes |
| betaine_peak 1.0 M | `mechanistic_params.py:218` | EMPIRICAL — peak enhancement, not full equalization (5.2 M per Rees 1993) |
| betaine_enhancement 0.12 | `mechanistic_params.py:219` | EMPIRICAL |
| betaine_inhibition_start 1.5 M | `mechanistic_params.py:220` | EMPIRICAL — practitioner consensus |
| mg_optimal 2.5 mM | `mechanistic_params.py:227` | EMPIRICAL — matches Rahman (2014) PLoS One 9:e112515 balanced-GC recommendation |
| Gaussian delta_t width 5 C | `mechanistic_params.py:~243` | EMPIRICAL — Gaussian for effective-Tm scoring |

## Validation ranges (ParamValidator)

The bounds in `neoswga/core/param_validator.py:PARAM_RANGES` mirror the
authoritative physical bounds in `reaction_conditions._validate()`:

- DMSO 0-10% (10% is PCR-safe upper; higher inhibits polymerases)
- Betaine 0-2.5 M (near saturation)
- Trehalose 0-1.0 M (practical SWGA range)
- Formamide 0-10%
- Ethanol 0-5% (>5% inhibits polymerases)
- Urea 0-2.0 M (practical range for short oligos)
- TMAC 0-0.1 M (higher is excessive)
- Mg2+ 0-20 mM (10-15 mM typical max; 0.5 mM practical min)
- Reaction temp 20-70 C (covers all four polymerases)

## How to audit a constant

Before changing any value in the tables above:

1. Find it in `tests/test_scientific_constants.py` and update the
   asserted expected value along with its citation link.
2. Update this document with the new value and source.
3. If the value is empirical, add a comment in-code explaining the
   calibration basis (dataset name, issue ID, or "in-house titration
   YYYY-MM-DD").
4. Run the nightly CI workflow which exercises the full regression
   suite under the new value to confirm no downstream test breaks.
