# Additive and coverage audit

Audited 2026-09-14 at commit `54165d5`. This is a code and literature audit,
not an experimental validation. Production calculations were not changed.

## Assessment

The pipeline can express many experimental conditions, and the common condition
builder forwards most buffer and additive settings. However, its numerical
predictions combine empirical melting-temperature corrections, estimated
coefficients, software heuristics and geometric coverage scores. They should
not currently be used to claim an optimal reaction recipe or a specified
fraction of genome recovered in the laboratory.

Some literature anchors are valid, but several exact numerical attributions are
incorrect or could not be verified. Successful unit tests establish software
behavior; they do not establish the biological accuracy of these models.

## Pipeline trace

| Stage | Current behavior | Consequence when changing conditions |
|---|---|---|
| Configuration | `parameter.py` loads additive/buffer settings; `build_reaction_conditions` forwards constructor fields, with explicit overrides before configuration and global defaults. | Most settings are connected, but `primer_conc` is not a constructor field. |
| Candidate counting | Exact sequence counts and positions are sequence-dependent. | Additives should not alter these counts. Topology and references do affect indexes. |
| Filtering | `filter_extra` calls additive- and salt-corrected Tm, then sequence-quality and complementarity checks. | Tm changes can remove candidates; redesign requires rerunning filtering and downstream stages. A saved candidate pool cannot recover candidates excluded under previous conditions. |
| Candidate scoring | Current default stage 3 preserves filtered candidates; legacy amplification scoring is opt-in with `--amp-model`. Optional quality scorers use conditions. | An old empirical amplification model is not validation for changed chemistry. |
| Set evaluation | `BaseOptimizer` computes raw windows, occupancy-weighted windows, and occupancy-weighted target/background loads. | Tm-active additives alter effective coverage and selectivity; raw geometric coverage is unchanged for fixed sites/reach. |
| Additive optimization | `AdditiveOptimizer` combines mechanistic heuristics, a synthetic representative primer, optional specificity callback, and bounded recipe searches. | A high score is a heuristic ranking, not demonstrated experimental performance. |
| Report | `plan-pool` recommends panels against the chosen metric and constraints; `report-pool` renders saved results. | Regenerating a report does not recompute conditions. Output must distinguish raw coverage, effective coverage and any separate simulation. |

## Findings, ordered by priority

### F1 — P1: configured primer concentration is ignored in the main Tm/occupancy path

`ReactionConditions.calculate_effective_tm` defaults to 0.5 uM
(`reaction_conditions.py:514`). The constructor does not retain `primer_conc`,
and filtering, occupancy and effective coverage call this method without passing
the configured concentration (`filter.py:522`, `occupancy.py:159`,
`base_optimizer.py:1380`). The common builder therefore does not fix this gap.

Probe: configuring 0.05 versus 5 uM returns the same 47.717 C for the test oligo.
Passing those concentrations explicitly returns 42.678 versus 52.921 C. This
can affect candidate eligibility, binding estimates and the pool-size curve.
The model also needs to state whether concentration is per oligo or total pool;
increasing pool size at fixed per-oligo concentration is a different experiment
from increasing it at fixed total concentration.

**Action:** store and propagate per-oligo concentration in the condition object;
support or explicitly exclude fixed-total-pool designs. Test through filtering
and final metrics, not only the low-level Tm function.

### F2 — P1: coverage percentage is an uncalibrated geometric/occupancy proxy

`base_optimizer.py:1334` marks symmetric windows around exact sites and weights
them by a sigmoid derived from Tm and duplex enthalpy. It does not resolve
strand-directed extension, competing template molecules, repeated priming,
enzyme activity, reaction duration, depletion, sequencing depth or breadth at a
specified read-depth threshold. An equilibrium binding fraction is not by
itself a probability of eventual amplification or read recovery.

The default 3 kb phi29 reach derives from design spacing choices, not a measured
extension distribution. `coverage.py:224` still claims neighboring products
truncate extension, whereas the explanatory text in `reaction_conditions.py`
explicitly rejects that explanation. Also, swga2's published `coverage_ratio`
counts opposite-strand neighboring sites within 70 kb; it is not the same
quantity as NeoSWGA's fraction of bases inside symmetric windows. See
[Dwivedi-Yu et al. 2023](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010137).

For the same saved 26-oligo wMel panel and identical conditions:

| Window radius | Effective coverage estimate |
|---|---:|
| 1,000 bp | 41.3% |
| 3,000 bp | 80.1% |
| 5,000 bp | 93.5% |
| 10,000 bp | 99.6% |

These are sensitivity calculations, not new evidence of recovery. The previous
80% recommendation is valid only as a result of the stated 3 kb model.

**Action:** label it a coverage proxy, always show reach sensitivity, and fit
observed sequencing breadth/depth across held-out reactions before presenting
recovery probabilities. `calibrate-reach` is a useful existing diagnostic but
fits a symmetric triangular depth profile on a selected contig; it does not
validate the full additive/occupancy model or identify physical processivity.

### F3 — P1: FASTA record boundaries are lost before position scanning

`string_search.py:111` concatenates records without separators and the scanning
functions use the resulting string. A synthetic two-record FASTA produces the
12-mer `CCCCCCGGGGGG` across a boundary despite neither record containing it.
Coverage windows also have no record boundaries within a prefix, so they can
reach from one contig into another. Per-prefix aggregation fixes separate
files, but does not fix multiple contigs in one file.

**Impact:** target coverage and exact background matches can be wrong for
fragmented/multichromosome references. The wMel target is a single record, but
its Drosophila background has many records; background comparisons remain
exposed to this issue. The magnitude on that reference was not quantified here.

**Action:** index record IDs and local coordinates, preserve contig lengths,
and clip/wrap within each record according to its topology.

### F4 — P1: a small circular target can receive 100% coverage without a site

`coverage.py:77` marks the entire target when `2 * extension >= length` before
checking for any primer positions. An empty-position cache returns coverage
1.0 for a 1 kb circular target and 1 kb reach. The marginal-curve helper has the
same shortcut. This is a concrete bug in these helpers, not a claim that every
coverage implementation has it; the main optimizer has a separate path.

**Action:** require at least one valid position before the circular shortcut;
add empty-site, absent-contig and short-circle regression cases.

### F5 — P2: the documented occupancy formula differs from the implementation

The docstring in `base_optimizer.py:1334` describes a product over sites. The
implementation at `:1378` forms a union for each primer and applies its weight
once, then combines distinct primers. At T = Tm with 100 bp windows on a 1 kb
target, sites at 500 and 510 give 10.5% for one primer sequence but 15.25% for
two distinct sequences with the same imposed occupancy. A site-independent
model would treat these cases differently from the present grouping.

Neither interpretation has been validated here. The comments' claim of an
upper bound on single-molecule recovery is not established by this arithmetic.

**Action:** specify and test the intended approximation. Do not silently replace
it with per-site independence and present the larger number as validated.

### F6 — P2: additive conversion and interaction paths are inconsistent

- `ReactionConditions.from_additives` omits `propanediol_m`
  (`reaction_conditions.py:1024`). A 1 M input becomes 0 M.
- The legacy `ReactionConditions.calculate_tm_correction(use_arrhenius=False)`
  also omits propanediol, while the additive dataclass's legacy method includes it.
- `MechanisticModel.calculate_effects` changes `tm_correction` for interactions,
  but `_calculate_effective_tm` normally delegates to the condition object and
  ignores that argument (`mechanistic_model.py:313`). Supplying corrections of
  -1 and -10 C gives the same effective Tm. The fallback uses a separate
  correction model, including a different GC-equalization curve.

**Action:** use one canonical Tm calculation, explicitly apply supported
interaction deltas, and test condition conversion round trips for every field.

### F7 — P2: non-Tm additive effects do not enter pool-plan effective coverage

Glycerol, BSA, PEG, SSB and DTT have effects in the optional mechanistic model
but no direct term in the main occupancy-window coverage calculation. Enzyme
inhibition by a Tm-active additive is likewise not incorporated there. Merely
forwarding a setting into `ReactionConditions` does not mean all predictions
use its biological effects. This is a model limitation, not evidence that an
additive has no experimental effect.

**Action:** report which effects each metric includes. A condition optimizer
must not optimize only binding discrimination while implying that enzyme
activity and genome recovery were assessed.

### F8 — P2: confidence and primer-length recommendations overstate evidence

`additive_optimizer.py:496` maps score >0.5 to `high` confidence without an
uncertainty model. Its representative sequence is assembled from blocks of
G/C/A/T based on length and GC (`:432`), rather than the actual designed pool.
`reaction_conditions.py:703` uses additive thresholds to permit longer primers
and calls 18 bp a literature-validated limit. The cited PCR papers do not
validate these SWGA length/dose thresholds. Musso 2006 concerns GC-rich PCR,
not the SWGA optimization attributed to it in `additive_optimizer.py`.

**Action:** replace confidence with a clearly labeled heuristic score; separate
experimental constraints from software search bounds. Evaluate actual pools
where sequence-specific claims are made.

## Literature values

| Model component | Assessment of current evidence |
|---|---|
| DMSO, -0.55 C/% | Plausible empirical approximation, not a universal constant. [Chester & Marshak 1993](https://pubmed.ncbi.nlm.nih.gov/8470801/) concerns primer Tm/PCR. [von Ahsen et al. 2001](https://pubmed.ncbi.nlm.nih.gov/11673362/) fits -0.75 C/% under its conditions. The code's exact 37 C anchor and activation energy need an identifiable derivation. |
| Betaine, -1.2 C/M plus GC equalization | [Rees et al. 1993](https://pubmed.ncbi.nlm.nih.gov/8418834/) supports approximately 5.2 M isostabilization. It does not by itself validate a linear dose response, Wallace-rule correction, or the uniform term for short SWGA oligos. [Henke et al. 1997](https://academic.oup.com/nar/article/25/19/3957/2549226) is a PCR enhancement study and tests practical betaine concentrations, contradicting the code comment that neither cited study did so. The claimed -1.3 C/M measurement from Henke was not substantiated. |
| Trehalose, -3 C/M | General Tm/protein effects have literature support, but the cited Spiess paper is [Clinical Chemistry 50:1256-1259](https://pubmed.ncbi.nlm.nih.gov/15229160/), not BioTechniques 36:732-736. The exact claimed 2-4 C/M range was not verified. Horakova reports 0.7-1.5 C reduction at 0.2 M, showing dependence on assay conditions rather than establishing the current coefficient. |
| 1,2-propanediol, -5.4 C/M | The strongest direct numerical match: [Horakova et al. 2011](https://link.springer.com/article/10.1186/1472-6750-11-41) reports a 4.9-5.9 C shift at 1 M on a short duplex. Treat -5.4 as a local empirical midpoint, not a validated linear model across doses, temperatures, buffers and polymerases. |
| Urea, -2.5 C/M | The attribution needs repair. [Hutton 1977](https://academic.oup.com/nar/article/4/10/3537/2380527) reports -2.25 C/M over 0-8 M, not the -5 C/M claimed in `mechanistic_params.py:95`. The cited Lesnick & Bhalla 1995 source was not located from the supplied details. The number is near Hutton's result, but its claimed provenance and the additional GC preference are not validated. |
| Formamide, -0.65 C/% | Broad scale is plausible: Hutton reports -0.60 C/% over specified salt/concentration ranges. Extrapolation into a phi29 reaction and the added temperature dependence require separate evidence. |
| TMAC, -0.5 C/M plus equalization | High-concentration isostabilization is a literature anchor; a uniform term and linear interpolation down to the configured 0-0.1 M range are extrapolations, not validated short-oligo SWGA corrections. |
| Ethanol, -0.4 C/% | Not substantiated from the cited [Cheng et al. 1994](https://doi.org/10.1073/pnas.91.12.5695) long-PCR study, which describes glycerol and DMSO. Retain as unsupported until a relevant measurement is provided. |
| Arrhenius scaling for all these Tm terms | Activation energies in `mechanistic_params.py:48-140` are marked estimated, with no reproducible fits or uncertainties. A kinetic Arrhenius law applied to an equilibrium Tm shift needs an empirical or thermodynamic justification. The default path should not be described as more accurate merely because it is temperature-dependent. |
| Additive combinations / enzyme-response multipliers | Effects can depend on polymerase, buffer, template and inhibitors. [Musso et al. 2006](https://pmc.ncbi.nlm.nih.gov/articles/PMC1876170/) supports a specific PCR combination including 7-deaza-dGTP; it does not validate generic phi29 multipliers or recipe confidence. |

Not all full texts were accessible. An unverified coefficient is not asserted
to be false; it is insufficiently supported for its current level of certainty.
No quantitative literature evidence found here calibrates the combined model
to recovery from the Wolbachia/Drosophila example.

## Recommended sequence of work

1. Fix F1, F3, F4 and F6 with regression cases before further condition searches.
2. Correct citations and attach source, experimental domain, evidence category
   and uncertainty to each coefficient. Mark unsupported corrections as
   experimental or disable them by default; do not substitute new point values
   merely because they look plausible.
3. Use explicit names for raw window coverage, occupancy-weighted coverage,
   simulated amplification and measured sequencing breadth. Show sensitivity
   to reach and concentration alongside pool-size recommendations.
4. Define the intended coverage target: fraction of reference bases at a
   specified read depth and sequencing budget, plus uniformity and specificity.
5. Calibrate and validate against independent SWGA reactions spanning actual
   primer pools and relevant conditions. Hold out complete reactions/pools,
   not random nearby bases from the same depth profile. Report prediction
   errors and uncertainty before using model thresholds as recovery targets.

## Verification

147 existing tests passed across additive wiring, condition completeness,
corrections, occupancy, effective coverage, explicit additive overrides and
scoring propagation. These mostly test internal consistency/directionality.

The read-only [probes](probes.py) reproduce concentration loss, propanediol
loss, ignored correction arguments, empty circular coverage, primer grouping,
record concatenation and (when local indexes exist) wMel reach sensitivity.
Observed values are saved in [probe_results.json](probe_results.json).

```bash
PYTHONPATH=. python docs/validation/additives_coverage_audit_2026-09-14/probes.py
```

No coefficients, production behavior or prior design results were modified by
this audit. The probe script writes only its numerical evidence file.
