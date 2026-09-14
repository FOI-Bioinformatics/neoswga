# Plan: additive and coverage audit remediation

Against [the audit](../validation/additives_coverage_audit_2026-09-14/README.md),
written at `54165d5`. It follows the audit's own **Recommended sequence of
work**; its five steps become the five phases below.

Every finding was reproduced before planning. `PYTHONPATH=. python
docs/validation/additives_coverage_audit_2026-09-14/probes.py` reruns the
evidence, and F1 was additionally confirmed by hand: `ReactionConditions` has no
`primer_conc` attribute, and `calculate_effective_tm` takes it as an argument
defaulting to 0.5 uM.

**Line numbers in the audit have drifted.** It was written before the formatting
pass in `1bd28e1`, so anchor on symbol names rather than the numbers it quotes.

## Status, 14 September 2026

Phases 1 to 3 are implemented. Phase 4 is a decision and phase 5 needs data that
does not exist in this repository; neither is startable here.

| Task | State | Pinned by |
|---|---|---|
| 1.1 F4 circular shortcut | Done | `test_circular_coverage_needs_a_site.py` |
| 1.2 F6 one canonical Tm path | Done | `test_one_canonical_tm_path.py` |
| 1.3 F1 primer concentration | Done | `test_primer_concentration_reaches_the_tm.py` |
| 1.4 F3 record boundaries | Done, both halves | `test_kmers_do_not_span_record_boundaries.py`, `test_coverage_windows_stop_at_record_edges.py` |
| 2.1 coefficient provenance | Done | `test_additive_coefficients_carry_their_evidence.py` |
| 2.2 citation corrections | Done | prose, in `mechanistic_params.py` and `additives.py` |
| 2.3 unsupported corrections | Labelled, not disabled | as above |
| 2.4 F8 confidence and length claims | Done | `test_recommendation_confidence_is_labelled_heuristic.py` |
| 3.1 distinct names | Partial: reports only, JSON keys unchanged | - |
| 3.2 F2 proxy label and sensitivity | Done | `test_pool_plan_shows_reach_sensitivity.py` |
| 3.3 F7 metric scope | Done | rendered in the pool-plan report |
| 3.4 F5 occupancy specification | Done | `test_occupancy_grouping_is_specified.py` |

Suite: 4458 passing, 75 skipped, no failures. 65 tests added.

**Measurements the plan asked for.**

F3's impact on the shipped pair was measured and is nil. The Drosophila
background holds 1870 records across 143.7 Mb, and of the 148 background sites
the delivered 26-oligo wMel panel has, none spanned a record join. The defect is
real and now guarded; it did not move this design.

The window half of F3 was measured the same way and is also small on that
pair: 3 of the delivered panel's 148 background sites sit within one reach of a
record edge, and the overcount is 11,289 bases, 0.008% of the assembly. The
exposure is in draft references rather than this one. 92.3% of that background's
1870 records are shorter than twice the 3 kb reach, but they hold only 2.28% of
the sequence, and six chromosome-scale records hold 94.8%. On an assembly made
only of contigs that size every window would cross an edge.

Record starts are now stored in the position index under a reserved key that
cannot collide with a primer. An index written before this has none and the
helpers behave exactly as they did, which was verified against the shipped
Drosophila index.

F1's blast radius is the opposite. No shipped configuration sets `ethanol_percent`,
`urea_m`, `formamide_percent` or `propanediol_m` at all, and one sets `tmac_m`,
so labelling the unsupported coefficients rather than disabling them changes no
delivered design today.

**Two decisions taken, both recorded in code.**

`primer_conc` is per oligo. Nearest-neighbour Tm depends on the concentration of
the annealing strand, so a 96-oligo panel at 0.5 uM each is what is modelled. A
fixed total concentration shared across a growing panel is a different
experiment and is explicitly not supported.

Unsupported coefficients are labelled, not disabled. Ethanol carries
`evidence: "unsupported"` and keeps its value. Disabling it would change
delivered designs for anyone who sets it; labelling states the position without
that cost, and the measurement above says nobody currently does.

**Three things deliberately not done.**

The audit's F6 probe still shows `tm_correction` of -1 C and -10 C producing the
same effective Tm, and that is now by design. Both the conditions object and the
mechanistic model compute a full additive correction, so adding the model's on
top of the canonical one would count every additive twice. What was genuinely
missing -- the delta an additive INTERACTION contributes -- is now applied
explicitly and verified end to end: a 1.5x Tm interaction modifier moves the
effective Tm by 3.65 C where it previously moved nothing.

The four quantities are named distinctly in the reports but the summary JSON
keys are unchanged. `fg_coverage` and `effective_fg_coverage` were already
distinct; renaming them would break every existing report consumer for a
cosmetic gain.

---

## What kind of work this is

The eight findings are not one kind of problem, and treating them as one is the
main way this plan could go wrong.

| Kind | Findings | What "done" means |
|---|---|---|
| Defect with a right answer | F1, F3, F4, F6 | The code does what it already claims. Pinned by a regression test |
| Claim not supported by its evidence | F8, and the literature table | The text matches what the source says |
| Modelling choice presented as settled | F2, F5, F7 | The choice is stated, its scope is labelled, and a decision is recorded |

Only the first kind can be finished without a decision. The second needs sources
re-read. The third needs someone to say what the model is for, and phases 4 and
5 need laboratory data that does not exist in this repository.

---

## Phase 1 — the four defects, with regression cases

The audit puts F1, F3, F4 and F6 first and says to land them "before further
condition searches". Ordered here by dependency rather than by number, because
two of them touch the same call sites.

### 1.1 F4: a circular target cannot be covered by nothing

`coverage.py` marks the whole target when `circular and 2 * extension >= length`,
and its own comment says "covered by any single site" while the code never
checks that a site exists. An empty position cache returns 1.0 for a 1 kb
circular target at 1 kb reach.

- Require at least one position from at least one primer before the shortcut.
- Apply the same guard to the marginal-coverage helper, which shares it.
- Tests: empty cache, absent contig, and a circle shorter than twice the reach
  with exactly one site. The first must return 0.0, not 1.0.

Smallest and wholly isolated. Do it first so the phase starts green.

### 1.2 F6: one canonical melting-temperature path

Three separate inconsistencies, all in the same subsystem, so fix them together
or fix the same call sites twice.

- `ReactionConditions.from_additives` does not mention `propanediol_m` at all,
  so 1 M becomes 0 M. Probe confirms.
- The legacy `calculate_tm_correction(use_arrhenius=False)` omits propanediol
  while the additive dataclass's legacy method includes it.
- `MechanisticModel.calculate_effects` computes a `tm_correction` for additive
  interactions, and `_calculate_effective_tm` delegates to the condition object
  and ignores it. Corrections of -1 C and -10 C both return 47.717 C.

Work: make one function the single place an effective Tm is computed; have the
mechanistic path apply its interaction delta explicitly through that function
rather than alongside it; add a round-trip test asserting **every** field of the
additive dataclass survives conversion, so the next omission fails rather than
silently zeroing.

### 1.3 F1: the configured primer concentration must reach the calculation

Depends on 1.2. Propagating a concentration through three call sites is wasted
work if there are still three Tm paths to propagate it through.

`calculate_effective_tm` defaults to 0.5 uM, the constructor discards
`primer_conc`, and filtering, occupancy and effective coverage all call it
without the configured value. Configuring 0.05 uM against 5 uM returns an
identical 47.717 C; passing them explicitly returns 42.678 C and 52.921 C, a
10 C spread that decides which candidates survive the Tm gate.

- Retain `primer_conc` on the condition object and pass it at every call site.
- Decide and document whether the figure is **per oligo or total pool**. These
  are different experiments: raising panel size at fixed per-oligo
  concentration is not the same reaction as raising it at fixed total. A
  96-oligo panel at 5 uM each is not a 96-oligo panel sharing 5 uM.
- Test through `filter` and the delivered metrics, not only the Tm function, so
  the test fails if the value stops arriving anywhere in between.

This is a seventh instance of the inert-key class recorded in CLAUDE.md Known
Issue 8, and the most consequential one found so far: the other six change
nothing when ignored, this one moves Tm.

### 1.4 F3: record boundaries inside one FASTA

`string_search.py` concatenates records with no separator, so a synthetic
two-record file yields the 12-mer `CCCCCCGGGGGG` spanning the join, present in
neither record. Coverage windows have no record boundaries either, so a window
can reach out of one contig and into the next.

Largest change of the four, because it alters what the position index stores.
Per-prefix aggregation already separates distinct files; this is about multiple
contigs within one.

- Index record id and local coordinate; keep per-record lengths.
- Clip or wrap within each record according to that record's topology.
- Tests: the planted cross-boundary k-mer must not be found; a window at a
  record end must not extend into the next record.

**Not yet quantified:** the shipped wMel target is a single record, but its
Drosophila background has many, so background counts are exposed. Measure the
error on that pair before and after, and record it. If it is negligible the fix
is still right, but the size of the claim should match the measurement.

---

## Phase 2 — citations, provenance and overstated confidence

The audit's step 2. No production number should change merely because a better
one looks plausible; the output of this phase is provenance, not new constants.

### 2.1 Attach provenance to every coefficient

For each additive correction record source, experimental domain, evidence
category and uncertainty, beside the value. The audit's literature table is the
starting inventory.

### 2.2 Correct the attributions that are wrong

These are errors of citation, not of value, and each is cheap:

- Trehalose cites BioTechniques 36:732-736; the Spiess paper is Clinical
  Chemistry 50:1256-1259.
- Urea cites -5 C/M to Hutton 1977, which reports -2.25 C/M over 0-8 M. The
  shipped -2.5 C/M is near Hutton; its stated provenance is not.
- Ethanol's -0.4 C/% is attributed to Cheng 1994, a long-PCR study describing
  glycerol and DMSO. Unsupported as cited.
- A code comment says neither cited betaine study tested practical
  concentrations; Henke 1997 did. The -1.3 C/M attributed to Henke was not
  substantiated.
- Musso 2006 concerns GC-rich PCR including 7-deaza-dGTP, not the SWGA
  optimisation `additive_optimizer` attributes to it.

### 2.3 Decide what happens to unsupported corrections

**Needs a decision.** The audit's instruction is to mark them experimental or
disable them by default. Ethanol is the clearest case: no located source.
Disabling changes delivered designs for anyone setting `ethanol_percent`;
labelling does not. Propose labelling first, with a measured statement of how
many shipped configurations set each affected field.

### 2.4 F8: confidence, and the primer-length claim

- `additive_optimizer` maps a score above 0.5 to `high` confidence with no
  uncertainty model. Replace with a plainly labelled heuristic score.
- Its representative sequence is assembled from G/C/A/T blocks by length and GC
  rather than taken from the designed pool. Where a claim is sequence-specific,
  evaluate the actual pool.
- `reaction_conditions` calls 18 bp a literature-validated limit on the strength
  of PCR papers that do not address SWGA oligo length. Separate the experimental
  constraint from the software search bound and cite each for what it is.

---

## Phase 3 — say which quantity is being reported

The audit's step 3. Mostly naming and output, but it is what makes F2, F5 and F7
honest, and none of it needs new science.

### 3.1 Four different quantities, four names

Raw window coverage, occupancy-weighted coverage, simulated amplification and
measured sequencing breadth are currently discussed as though interchangeable.
Name them distinctly in code, in the summary JSON and in the reports.

### 3.2 F2: label the proxy and fix the contradiction

Effective coverage marks symmetric windows around exact sites and weights them
by a sigmoid of Tm and duplex enthalpy. It resolves neither strand direction nor
competing templates, repeated priming, enzyme activity, reaction duration,
depletion, nor depth at any threshold. An equilibrium binding fraction is not a
probability of recovery, and it should not be printed as though it were.

`coverage.py` also contradicts itself: one docstring says extension "is
truncated by neighbouring primers' strand-displacement", another says "strand
displacement means a downstream primer does not truncate" it. One of them is
wrong and the default 3 kb reach rests on the answer.

Always show reach sensitivity beside any pool-size recommendation. On the saved
26-oligo wMel panel under identical conditions the same design reads 41.3% at
1 kb, 80.1% at 3 kb, 93.5% at 5 kb and 99.6% at 10 kb. A single number from that
family, printed without its reach, is the whole problem in one line.

Note that swga2's published `coverage_ratio` counts opposite-strand neighbouring
sites within 70 kb. It is not this quantity and the two should not be compared.

### 3.3 F7: state which effects each metric includes

Glycerol, BSA, PEG, SSB and DTT have terms in the optional mechanistic model and
none in the occupancy-window coverage that `plan-pool` reports. Enzyme
inhibition by a Tm-active additive is likewise absent there. Forwarding a
setting into `ReactionConditions` is not the same as every prediction using it.
A condition optimiser must not optimise binding discrimination while implying
that enzyme activity and genome recovery were assessed.

### 3.4 F5: specify the occupancy approximation, then test it

The docstring describes a product over sites; the implementation forms a union
per primer, applies that primer's weight once, and then combines distinct
primers. At T = Tm with 100 bp windows on a 1 kb target, sites at 500 and 510
give 10.5% for one primer and 15.25% for two distinct primers at the same
imposed occupancy.

Neither reading is validated. Write down which is intended and test that; do not
swap in per-site independence and present the larger number as validated. The
comment claiming an upper bound on single-molecule recovery is not established
by this arithmetic and should go or be proved.

---

## Phase 4 — define the target (decision, not code)

The audit's step 4. A coverage target has to name a fraction of reference bases
**at a stated read depth and sequencing budget**, together with uniformity and
specificity. Until that exists, "80% coverage" names a model output and not an
experimental outcome, and no threshold in the tool can be called a recovery
target.

This is a decision for the maintainer. Nothing in phases 1 to 3 depends on it,
and everything in phase 5 does.

## Phase 5 — calibrate against real reactions (blocked)

The audit's step 5, and **not implementable from this repository.** It requires
independent SWGA reactions spanning actual primer pools and relevant conditions,
held out as complete reactions rather than as random nearby bases from one depth
profile, with prediction errors and uncertainty reported.

What exists today: 310 qPCR fold measurements in
`tests/validation/data/dwivedi_yu_2023_plasmid_amplification.json`, on a plasmid
with site geometry constant by construction. They can calibrate a
sequence-level amplification term and cannot calibrate coverage or selectivity.
`calibrate-reach` fits a symmetric triangular depth profile on a selected
contig; it is a useful diagnostic and does not identify physical processivity or
validate the additive and occupancy model.

Until phase 5 has data, the honest position is the audit's own: the tool can
express many conditions and rank them, and its numbers should not be used to
claim an optimal recipe or a fraction of genome recovered in the laboratory.

---

## Sequencing and checkpoints

1. **Phase 1** in the order 1.1, 1.2, 1.3, 1.4. Each task lands with its
   regression test and leaves `pytest tests/ -n 8` green. 1.4 additionally
   reports the measured error on the wMel/Drosophila pair.
2. **Phase 2** after phase 1, because 2.3 changes delivered designs and should
   not be entangled with defect fixes.
3. **Phase 3** can run in parallel with phase 2; it touches output and naming,
   not coefficients.
4. **Phase 4** is a conversation, not a branch.
5. **Phase 5** starts when there is data.

Re-run the audit probes after each phase. They are read-only apart from their
own results file, and they are the cheapest check that a fix did what it claimed.

## What this plan does not do

It does not change any additive coefficient. The audit is explicit that an
unverified coefficient is not thereby false, only insufficiently supported for
the certainty currently attached to it, and replacing one unsourced number with
another would repeat the error being corrected.
