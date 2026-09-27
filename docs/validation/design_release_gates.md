# Two gates, and what each one licenses

2026-09-27. Task 10 of the valid-design plan. This is a study design and a
release policy. It is not a reaction recipe, and it does not report a
validation outcome, because none has been attempted.

## Why two gates and not one

A design tool can be correct and still be useless, and it can be useful while
computing something other than what it claims. Those are different failures and
they are established by different evidence, so they get separate gates.

- **The software gate** asks whether the tool computes what it says it
  computes. It is settled by tests, ratchets and recorded measurements inside
  this repository, and it is the only gate this project has ever passed.
- **The empirical gate** asks whether a panel this tool designs recovers target
  sequence in a reaction. It is settled by sequencing, against endpoints
  declared before the outcome is seen. **It has never been attempted here.**

Passing the first licenses proxy claims. It licenses nothing about recovery.

## Claim tiers

Every figure a report shows belongs to exactly one tier, and the tier bounds
what may be said about it.

| Tier | What it is | What licenses it | Present today |
|---|---|---|---|
| **proxy** | a geometric or thermodynamic quantity computed from binding positions under a declared reach, such as `fg_coverage` | the software gate | yes, and it is all there is |
| **calibrated** | a prediction of sequencing breadth from a model fitted to measured breadth, with a recorded domain | the empirical gate, for that model, in its domain | **no** |
| **observed** | breadth or depth measured in a named experiment | the experiment, reported with its own uncertainty | **no** |

`PanelAssessment` already states the tier of its headline quantity in prose:
coverage carries the note "a geometric proxy at the declared reach. It is not a
predicted sequencing breadth and has not been calibrated against one." There is
no `claim_tier` field in the code. Adding one is worth doing and is not done.

## What the software gate covers

Concretely, and only this:

- The four pipeline steps run end to end on a packaged reference, and each
  required refusal survives argument parsing, parameter resolution, the step's
  own `except` clause and the command boundary
  (`tests/integration/test_strict_design_pipeline.py`).
- A required calculation that fails fails the run rather than returning a
  substitute (the design-failure contract).
- Coverage agrees with an independent base-by-base oracle written without
  calling the production interval code.
- Positions agree with a brute-force scan, and the two production scanners
  agree with each other.
- Every chemistry constant carries an evidence status, and a computation
  outside a recorded domain is refused. Of the 23 statuses recorded, 7 are
  `measured`, 6 `estimated`, 8 `assumed`, 1 `empirical` and 1 `absent`.
- The counted search allowance sees every panel evaluation made in a loop.

What it does not cover: whether either coverage reach is right. Neither has
been measured against a reaction, and `calibrate-reach` has never been run
against measured depth because there is no BAM or CRAM in this repository.

## The empirical study design

Predeclared, because an endpoint chosen after seeing the data is not an
endpoint. None of this has been run.

### Endpoints

Primary:

1. **Callable breadth at a stated depth.** Fraction of target bases at or above
   a depth named before sequencing, reported per replicate.

Secondary, all reported whether or not they favour the panel:

2. **Target read fraction.** Reads mapping to the target over reads passing
   quality filters.
3. **Breadth uniformity.** Gini of covered-interval gaps, on the same
   definition the design uses, so the proxy and the observation are comparable.
4. **Oligo count and cost.**

### Design

- **Independent amplification replicates**, so within-reaction variance is
  estimable rather than assumed.
- **Held-out experiments.** `sequencing_feedback.require_disjoint_experiments`
  already refuses a validation set sharing an experiment with the training set,
  and refuses an empty set on either side.
- **Comparators:** the previously published pool for the same target, and a
  baseline appropriate to the question.
- **One factor at a time.** A design differing in both chemistry and oligo
  length cannot attribute a combined change to either. Additive and length
  changes are separate arms.

### Acceptance margins

**No margin is stated here, and that is the point.** A margin must be:

1. derived from the coverage the intended application actually requires, not
   from what a tool happens to achieve;
2. sized by a pilot-based precision and power assessment, so the study can
   distinguish the margin from replicate noise;
3. versioned and frozen with the study, before any outcome is observed.

Software cannot invent a universal acceptable recovery threshold, and this
document will not supply one in its place. Until a pilot exists there is no
defensible margin, so there is no empirical gate to pass.

## Release policy

1. Ship proxy calculations once the software gate passes. Say proxy.
2. Promote a model to **calibrated** only after its own predeclared empirical
   gate passes, in the domain the pilot covered, citing the frozen record.
3. Report negative findings. Revise or retire a model that does not transfer.
4. Never substitute another model to make a validation pass. A model that fails
   its gate has told you something.

## Historical results are proxy results

Two figures in this repository are sometimes quoted as evidence that this work
improves recovery. They are not.

- The **9-oligo panel at 70.07%** from
  [sequential_panel_service_2026-09-19.md](sequential_panel_service_2026-09-19.md)
  is *predicted effective coverage* under the design's own occupancy model, at
  the realistic reach. It is a proxy figure. No sequencing was done, and the
  same document says so.
- Every coverage figure anywhere in `docs/validation/` is proxy, for the same
  reason.

Read them as what they are: evidence about the search, measured against the
tool's own objective. A search that finds a smaller panel at the same proxy
coverage has demonstrated something real about the search, and nothing about a
reaction.

## What would move this document

One pilot, with a depth threshold and a margin fixed in advance, on a target
this repository already carries a prepared design for. That is the smallest
thing that would turn any figure here from proxy into calibrated, and it needs
sequencing rather than code.
