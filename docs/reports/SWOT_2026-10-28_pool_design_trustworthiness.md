# SWOT: can this tool be trusted to deliver an oligo pool?

**Date:** 2026-10-28
**Subject:** the tool's capacity to deliver a *trustworthy* pool, read through
the four defects found in the 26 October optimizer and reachability audit. Not
"are the methods sound" and not "does the workflow hold together" — those are
the two SWOTs below. The question here is narrower and harder: when a run
prints a panel and a coverage figure, what is actually known about them?

**Relationship to the other reports.**
[`AUDIT_2026-08_alternatives_and_scaling.md`](AUDIT_2026-08_alternatives_and_scaling.md)
SWOTs the *methods*;
[`SWOT_2026-08_pipeline_as_workflow.md`](SWOT_2026-08_pipeline_as_workflow.md)
SWOTs the *workflow*;
[`AUDIT_2026-10-26_chemistry_in_pool_design.md`](AUDIT_2026-10-26_chemistry_in_pool_design.md)
asks how chemistry reaches selection. This document deliberately repeats none
of them. It is about the gap between what a run asserts and what it establishes.

**Evidence.** Every entry cites `file:line`, a dated validation record, or a
measurement made for this document. Where a figure is synthetic or a single
sample, it says so. Nothing here is a projection of benefit.

---

## Strengths

**S1. The project grades its own evidence, and refuses to model what it cannot.**
`registry/model_evidence.json` carries a status per chemistry constant; counted
on this tree, 27 entries: **7 measured, 7 estimated, 11 assumed, 1 empirical, 1
absent.** `absent` is load-bearing rather than decorative —
`model_evidence.py:92` makes it the one status that is not a model, and a design
setting glycerol is refused outright because nothing computes its Tm effect.
Declining to produce a number is rarer in this field than producing one.

**S2. "Unknown is not zero" is enforced in code, not aspired to in prose.**
The newest example is the one written this week: the swap Stage-2 score applies
its weighted host term only when both coverages are measurable and otherwise
falls back to the lexicographic form, rather than substituting a zero host
coverage (`swap_refinement.py:383-384`). A fabricated zero there would have been
invisible and would have read as the best possible specificity.

**S3. A measurement log that can be re-run, not a changelog of assertions.**
124 markdown files and ~44,500 lines under `docs/`, the bulk of it one record
per claim under `validation/`, most carrying a probe script that reproduces the
figure. The four defects of 26 October each shipped with one.

**S4. One input gives one panel.** `hybrid` Stage 2 seeds its
`algebraic_connectivity` calls (`_CONNECTIVITY_SEED`), and a source-level test
asserts every call carries the seed. Before 2026-10-26 one input delivered two
or three different panels
([`2026-10-26-stage2-determinism.md`](../validation/2026-10-26-stage2-determinism.md)).
Reproducibility is a precondition for every other claim here.

**S5. The documentation argues against the tool where the evidence does.**
[`2026-09-22-swga-comparative-review.md`](../validation/2026-09-22-swga-comparative-review.md)
records NeoSWGA's own directional-coverage defect in the sharpest available
terms — "`polymerase_extension_reach` returns an EXTENSION distance and the code
spends it as a RADIUS" — and credits COATswga for getting it right. A
comparison document that concedes the central metric is a strength of process.

---

## Weaknesses

**W1. Options failing to reach the deciding code is endemic, not residual.**
Known Issue 8 is explicitly never recorded as closed, and the 26 October audit
found **four fresh instances in one pass**, each shipped and each silent:

| Instance | Consequence |
|---|---|
| `stage2_tm_weight` dropped by the `HybridBaseOptimizer` wrapper | the Tm term shipped days earlier was inert on every CLI `hybrid` run |
| `dimer_temp` never forwarded by `DominatingSetAdapter` | free-energy dimer screen ran at 37 C under every reaction |
| `polymerase` / `reaction_temp` never passed by `unified_optimizer` | **every** `hybrid` and `background-aware` run built a phi29 preset at 30 C |
| `max_dimer_dg` read nowhere in `clique_optimizer` | the method sold on a dimer-free guarantee applied no free-energy floor |

All six reachability ratchets passed through all four, because they assert the
two ends and not the seam between an adapter and the engine that decides. The
class is now four variants deep (flags, schema keys, capabilities, arguments,
and the wrapper/inner seam).

**W2. The flagship default cannot discriminate, and the escape hatch was broken.**
Known Issue 17: at phi29 30 C the host and target k-mer distributions overlap so
heavily that no selection term separates them; the lever is the reaction, not
the score. That lever — a hotter polymerase — did not work until 2026-10-26,
because the configured polymerase never reached the optimizer (W1, row 3). The
two facts compound: the regime that cannot discriminate was also the only regime
anyone was actually running.

**W3. The dimer-tightening lever decides nothing at the shipped default.**
Measured for this document over 20,000 random 12-mer pairs, counting pairs that
pass the base-pair screen but fail a free-energy floor:

| `max_dimer_bp` | floor −6.0 | floor −3.0 | floor −1.0 |
|---|---|---|---|
| **3 (shipped default)** | **0** | **0** | 4,879 |
| 4 | 0 | 1,237 | 8,317 |
| 5 | 53 | 1,894 | 9,102 |
| 7 | 149 | 2,088 | 9,297 |

So `max_dimer_dg` is a real option that changes nothing unless the user also
loosens `max_dimer_bp` past 4 or picks an unusually shallow floor. This bounds
the practical reach of two of the four fixes above: both are correct, and at the
shipped default neither moves a panel. Consistent with
[`dimer_stability_floor_2026-09-18.md`](../validation/dimer_stability_floor_2026-09-18.md),
which reached the same conclusion from a different direction.

**W4. The headline number is inflated by 1.3x to 1.9x** (measured 2026-10-29;
this entry said "probably inflated, and by how much is unknown" when written).
At the shipped 3 kb reach a 12-primer panel reported at 0.72 coverage covers
0.46 under the directional model. See
[`2026-10-29-directional-coverage-measured.md`](../validation/2026-10-29-directional-coverage-measured.md);
the magnitude is for synthetic site geometry, and the mechanism below is what
generalises.
Both coverage entry points default to `strand="both"` (`coverage.py:28`, `:117`)
and `_mark_window` marks `occupied[pos-extension:pos+extension]`
(`coverage.py:234`) — a symmetric window centred on each site, so every
occurrence earns coverage in the direction it cannot extend in. `fg_coverage` is
the authoritative figure in every `step4_improved_df_summary.json` and across
`validation/`. The argument is structural and predicts inflation without
quantifying it; the expectation is small on dense panels and larger on sparse
ones, and sparse is where SWGA designs live.

**W5. The extension points are at their size ceilings simultaneously.**
`hybrid_optimizer.py` 1646/1650, `unified_optimizer.py` 1646/1650,
`parameter.py` 1596/1601 — headroom 4, 4 and 5, against a ratchet requiring 4.
The last new params.json key (`stage2_background_weight`, `40af6e0`) cost 6
lines in `parameter.py` and 3 in `unified_optimizer.py`. **The next one does not
fit**: it needs 6 lines where 5 remain. The architecture's documented extension
recipe is blocked on an extraction nobody has done.

**W6. Nothing is prospectively validated.** 11 of 27 chemistry constants are
`assumed`, both coverage reaches are `assumed`, and the betaine, trehalose and
urea coefficients each carry an in-code retraction. `docs/LIMITATIONS.md` and
`EVIDENCE.md` say so plainly, which belongs under Strengths — but the
underlying fact belongs here.

**W7. The suite is green only in CI.** 60 tests fail in a container lacking a
k-mer counter and `pysam` — 50 of them because a `ReferenceDataError` about
pysam arrives where a specific refusal was asserted. A contributor cannot
distinguish "my change broke something" from "my container is thin" without
running a baseline diff first, which is how this session verified every commit.

---

## Opportunities

**O1. The polymerase fix unlocks the measurement programme Known Issue 17
blocked.** Now that `bst` and `equiphi29` designs actually run their configured
reaction and their Stage-0 screen, the two new selection terms
(`stage2_tm_weight`, `stage2_background_weight`) can finally be measured in a
regime that *can* discriminate. Both currently ship off with "unmeasured against
a reaction" as the stated reason — and until 2026-10-26 no such reaction was
reachable through the default method. This is the highest-value measurement
available, and it was impossible last week.

**O2. Directional coverage is a bounded fix with a worked reference.**
COATswga's model is correct and short: a forward site covers `[pos, pos+frag]`,
a reverse-complement site `[pos+len-frag, pos+len]`, one direction per
occurrence. Fixing W4 corrects the tool's headline metric. It will move every
recorded `fg_coverage`, so it needs the measurement first and a re-baselining
plan after — but the direction of the error is already known.

**O3. Two cheap gaps with external justification.** `off_gap_gini` carries the
second-largest weight in swga 2.0's fitted model and is computed nowhere here. A
strand-balance *report* needs no new data — `strand_stats` already holds all five
figures per genome — and COATswga's mean-of-two-strands ranking shows the
failure it would catch (0.99/0.30 beating 0.60/0.60).

**O4. One resolution path instead of a seam per adapter.**
`stage2_ranking.resolve_reaction` is the template: explicit > conditions >
preset, in one place both wrappers pass through. Generalising it, and adding the
cheap structural check the audit found (diff an inner's `__init__` parameters
against what its adapter passes), would convert W1 from a recurring discovery
into a ratchet.

---

## Threats

**T1. Provenance: every hot-polymerase design produced before 2026-10-26 ran as
phi29 at 30 C.** This is the sharpest consequence of the audit and the only one
that reaches outside the repository. For `bst` or `equiphi29` the Stage-0
thermodynamic pre-screen was skipped entirely — phi29's preset sets
`thermo_filter=False` where theirs set it True — and both Stage-2 temperatures
ran at 30 C against a 63 C or 42 C reaction. Measured under the reaction of the
*shipped* `examples/equiphi29_scenario`, the pre-screen keeps 299 of 500
candidates, so 500/500 reached selection unscreened. Any equiphi29 or bst panel
in a notebook, an order, or a freezer predates the fix and was not designed
under the reaction its params.json declared.

**T2. The defect class produces plausible output, so users cannot detect it.**
An inert option changes nothing visible: the panel looks reasonable, the report
is fully populated, and in the `--application` case the run log actively
asserted the opposite ("selection weights applied: tm=0.30, uniformity=0.05,
dimer_penalty=0.45") while the object doing the selecting used 0.0. Silent
wrong answers are the hardest class to retire because no user report will ever
surface them.

**T3. Presentation outruns evidence.** The tool grades designs, badges values
MEASURED or ESTIMATED, and prints quality verdicts, on a model with 7 of 27
constants measured and nothing validated prospectively. The grade tables are
calibrated against the seven-dataset retrospective suite and documented as
descriptive rather than predictive — but a grade is read as a prediction
whatever the documentation says, and NeoSWGA is the only one of the four tools
that grades at all, so it is the surface with the least external calibration.

**T4. Drift is found incidentally, and it runs in both directions.** Three
entries in `tool_comparison.md` were stale this week, **all understating the
tool** (the k-mer counter, the retired RF, `bg_coverage`'s reader), and one was
falsified by a commit made two days earlier in this same session. A dead link to
a plan file that never existed in git survived in five documents. Across 124
documents and 149 `core/` modules, correctness of prose depends on someone
happening to read the right line.

**T5. The discipline is the asset, and it is concentrated.** Everything under
Strengths — the evidence registry, the ratchets, the probe-per-claim log, the
willingness to record a defect against the tool's own headline metric — is a
practice rather than a mechanism. The ratchets encode part of it; the judgement
that found four reachability defects in a day is not encoded anywhere.

---

## What this document does not establish

- The coverage inflation in W4 is argued from source and unmeasured; no figure
  here bounds it.
- W3's table is random 12-mers, not a candidate pool from a real design, and it
  speaks to what a screen rejects rather than to a delivered panel.
- T1's 299/500 is synthetic candidates under a real shipped configuration. No
  equiphi29 or bst design was re-run end to end, because that needs a k-mer
  counter this container does not have.
- No wet-lab outcome is cited anywhere above. The seven-dataset suite is
  retrospective, and T3 is about what that cannot support.
