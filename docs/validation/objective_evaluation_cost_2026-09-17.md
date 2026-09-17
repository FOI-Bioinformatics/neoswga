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

**Corrected 2026-09-17.** This section first attributed the gap to the
background, naming `_effective_site_load` and the background set union. That was
wrong. The two profiles differ because one computes occupancy-weighted coverage
and the other does not, and that term scales with the TARGET length and the
panel size, not with the background.

`docs/validation/pool_plan_profile_2026-09-16.md` recorded about 1.4 ms per
evaluation on a 2 Mb synthetic pair. Its driver passes `conditions=None` and
`coverage_metric="raw"`, so `_effective_site_load` returns immediately and
`_compute_effective_coverage` is never called. The two published figures are
measurements of two different quantities.

Measured on the real Wolbachia pair at panel size 24, changing one thing at a
time and holding the cache, the panel and both genomes fixed:

| Configuration | Per call |
|---|---|
| Occupancy-weighted coverage, background present | 18.9 ms |
| Raw coverage, background present | 0.9 ms |
| Occupancy-weighted coverage, no background at all | 15.8 ms |

Removing the 144 Mb host costs 16% of the time. Turning off the occupancy
weighting costs 95% of it. The cost also grows about linearly in the panel,
which is the signature of a per-primer pass over the target rather than anything
to do with the host:

| Panel size | Per call |
|---|---|
| 6 | 6.4 ms |
| 12 | 11.0 ms |
| 24 | 18.7 ms |
| 48 | 35.9 ms |

The conclusion above is unaffected: the objective costs milliseconds, so a
bounded per-step scan is still mandatory. What changes is where the effort
belongs. `_compute_effective_coverage` does, per primer, one full-length window
reset and one masked multiply over the target, and computing the same quantity
as an interval sweep over window endpoints measured 19x to 21x faster against
the real implementation on this pair, agreeing to 6e-9, with wrapping verified
separately. At under 1 ms for that term a 2,000-candidate frontier costs a few
seconds per scan rather than 30 s, which changes what scan width increment 4 can
afford. See
[parallelism_opportunities_2026-09-17.md](parallelism_opportunities_2026-09-17.md),
which found this while looking for something else.

The general point is the one the retention benchmark already made: measure on
the pair you design against. A synthetic profile is the right tool for finding
which function dominates and the wrong one for deciding what is affordable --
and it is the wrong tool for the first job too when it silently runs a different
code path, which is what happened here.

## The scan this prices has no production caller

**Added 2026-09-17.** Nothing in `neoswga/` passes `objective=` to
`optimize_greedy`. The greedy scores on coverage bins, and the two reachable
objective-scored searches are the swap repair and the beam, both already
budgeted. The arithmetic above is correct and the architectural conclusion
follows from it, but it prices a scan a user cannot currently run.

Measured on the same pair, one `plan-pool` size row makes ONE objective
evaluation when no constraint binds, and about 3,200 when the repair fires.
Neither is 491,836. See
[what_actually_bounds_the_search_2026-09-17.md](what_actually_bounds_the_search_2026-09-17.md),
which also finds that the exposed evaluation budget is not what stops the loop.

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
