# What one panel evaluation costs, and what that rules out

Measured 2026-09-17, on the real Wolbachia design rather than a synthetic one.
Phase 4 of the plan for `docs/validation/pipeline_audit_2026-09-16/` names this
as the measurement to take before committing to a search architecture.

Reproduce with `scripts/benchmarking/measure_objective_cost.py`.

**A full-universe objective scan is not available at any scale this tool works
at.** Bounding the per-step scan is the only viable design, not an optimisation
to add later.

## The measurement

`PoolObjective.coverage(panel + [candidate])` at panel size 24, against wMel
`NC_002978.6` with the Drosophila `GCF_000001215.4` background, k = 12, phi29 at
30 C.

| Calls | Wall clock | Per call |
|---|---|---|
| 200 | 3.7 s | 18.5 ms |
| 1,000 | 14.9 s | 14.9 ms |

The 1,000-call figure is the one to use; the shorter run carries warm-up. Single
measurements on one machine, so treat these as an order of magnitude rather than
a constant.

The plan's threshold was 200 microseconds. The measurement is 75 times that.

## What it rules out

Scoring every candidate with the full objective at every greedy step, at 14.9 ms
per candidate:

| Candidates | One greedy step | A 24-primer panel |
|---|---|---|
| 2,000 (today's shortlist) | 30 s | 0.2 h |
| 20,670 (post-Gini) | 308 s | 2.1 h |
| 491,836 (all hard QC) | 7,329 s | 48.9 h |

Forty-nine hours for one panel, and a size sweep evaluates twenty of them. The
decision to search the full retained inventory stands; what cannot stand is
doing it by scoring everything at every step.

## Why the earlier profile was optimistic

`docs/validation/pool_plan_profile_2026-09-16.md` recorded about 1.4 ms per
evaluation. That was a 2 Mb synthetic target with a synthetic background of the
same size. This is a 1.27 Mb target against a 144 Mb host, and both
`_effective_site_load` and the background set union scale with the background.
An order of magnitude separates them, and the difference is the background, not
the target.

The general point is the one the retention benchmark already made: measure on
the pair you design against. A synthetic profile is the right tool for finding
which function dominates and the wrong one for deciding what is affordable.

## What this does not say

- It does not say the greedy is slow today. The greedy scores on coverage bins,
  not on the objective; the objective enters through the repair swaps and the
  beam, both of which carry budgets. A `plan-pool` run over nine sizes takes
  about 20 s.
- It does not measure the bin-based scan, which is the cheap prescreen any
  bounded design would put in front of the objective.
- One machine, one pair, one panel size. The panel size matters: the cost grows
  with the panel, since every evaluation re-unions the whole panel's sites.

## Consequence for the search design

A bounded frontier is necessary but not sufficient. At 14.9 ms even a 2,000
candidate frontier costs 30 s per greedy step if every candidate is scored, so
the objective cannot be the per-candidate criterion at all. The design the plan
sketches is the one the measurement supports: a cheap bin-gain prescreen over
the frontier, then the full objective over a small number of leaders.

At a scan width of 64 the objective costs about 1 s per step and 24 s per
24-primer panel, which is affordable. The width that is actually needed is a
separate question, settled by comparing delivered panels against an unbounded
scan on a pool small enough to afford one.


## A defect this measurement uncovered

The first attempt to read the inventory back from `plan-pool` found nothing and
fell back to the CSV, reporting "no eligible candidates under this reaction" --
which was true, and was not the reason.

The inventory files every verdict under a reaction fingerprint. The filter
writes it through `filter._get_reaction_conditions`; `plan-pool` computed it
through `build_reaction_conditions(SimpleNamespace(**params))`. From the same
params.json those disagreed:

    filter     tm-2026-09-14:2b3536d4305777e8
    plan-pool  tm-2026-09-14:6dc7004c72a51941

One field, `dtt_mm`, 4.0 against 0.0.

`ReactionConditions` resolves `mg_conc=None` and `dtt_mm=None` to the
polymerase's own buffer, and both carry comments explaining why a literal zero
is wrong: no polymerase runs without magnesium, and the mechanistic model scores
absent DTT as a deficiency worth 20% of stability.
`build_reaction_conditions` passes a value through whenever it is not None, and
`parameter.mg_conc` and `parameter.dtt_mm` were 2.0 and 0.0 at module scope. A
path that had run `get_params` got the resolved values; one that had not got the
stale module defaults. `parameter.reaction_temp` was already `None` for exactly
this reason; these two now match it.

2.0 mM is the old PCR magnesium figure that schema version 2 corrected to 10 mM.
So any path building conditions without `get_params`, and any params file not
setting `mg_conc` explicitly, was designing at the pre-correction value.

The point for this note is narrower: with nothing in production reading the
inventory, a key that could never match was invisible. Wiring the reader is what
found it.
