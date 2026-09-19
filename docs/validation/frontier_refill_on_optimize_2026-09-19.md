# Widening `optimize`'s frontier changes the panel, and mostly for the worse

Measured 2026-09-19 on the Wolbachia/Drosophila pair, before doing the work
rather than after.

Audit finding F1 is only partly closed. All three commands open the shared
inventory-backed candidate source, but the frontier opens at the size of the
list the command would have read anyway and `pool_planner.py` is the only
caller of `source.advance()`. So `optimize` and `expand-primers` still search a
`max_primer`-sized universe, which is the effect F1 named.

The obvious remedy is to give `optimize` a refill too. This measures whether
that would buy anything.

**It would not, and the prior stated when this was proposed was wrong.** That
prior -- that a refill changes nothing, from `plan-pool`'s own refill
measurement and the retention benchmark -- does not transfer, and the reason it
does not is the point of this document.

## What was run

`count-kmers`, `filter` and `score` from
`examples/wolbachia_pool_design/params.json` with `candidate_retention` at
`post_gini`: wMel `NC_002978.6` as the target, Drosophila `GCF_000001215.4` as
the background, k = 12, phi29 at 30 C, `max_primer` 2000, `max_gini` 0.7. The
funnel retains 491,836 hard-QC survivors, 20,670 clear the evenness gate and
are indexed, and 2,000 are shortlisted.

Then `neoswga optimize` twice per panel size, differing in ONE file:

- **narrow** -- `step3_df.csv` as written, the 2,000-primer shortlist. This is
  what `optimize` searches today.
- **wide** -- `step3_df.csv` replaced by all 20,670 eligible candidates in the
  inventory's `search_rank` order. This is what a frontier refill would reach.

Same params, same seed, same position files, same method (`hybrid`). No panel
limits are configured, which is the default and is load-bearing here.

Raw numbers are in
[frontier_refill_on_optimize_2026-09-19/results.json](frontier_refill_on_optimize_2026-09-19/results.json),
including every delivered panel. Reproduce with `build_wide.py` (writes the
20,670-row `step3_df.csv` from the inventory) and `run_pair.sh` (runs the two
`optimize` invocations and copies the outputs out), both beside it.

## The panels differ, and they differ more the larger the panel

| n | seed | Jaccard | primers outside the shortlist | coverage | selectivity density | host sites | worst hole (bp) |
|---|---|---|---|---|---|---|---|
| 6 | 20260919 | 0.500 | 2 of 6 | 0.6503 -> 0.6614 | 33.594 -> 29.415 (-12.4%) | 92 -> 93 | 65,860 -> 40,962 |
| 12 | 20260919 | 0.263 | 7 of 12 | 0.7334 -> 0.7408 | 25.616 -> 14.760 (-42.4%) | 149 -> 209 | 43,615 -> 28,070 |
| 12 | 777 | 0.263 | 7 of 12 | 0.7334 -> 0.7408 | 25.616 -> 14.760 (-42.4%) | 149 -> 209 | 43,615 -> 28,070 |
| 24 | 20260919 | 0.171 | 15 of 24 | 0.7772 -> 0.8378 | 19.463 -> 11.735 (-39.7%) | 261 -> 375 | 21,935 -> 21,935 |
| 24 | 777 | 0.171 | 15 of 24 | 0.7772 -> 0.8378 | 19.463 -> 11.735 (-39.7%) | 261 -> 375 | 21,935 -> 21,935 |

The two seeds return identical panels at both sizes. `hybrid`'s deciding stages
are deterministic here; the seed reaches network refinement, which is not what
selects on this path. So the differences below are the frontier's doing and not
run-to-run noise.

Cost is not the objection. Narrow against wide: 22.1 s -> 39.8 s at n=6,
39 s -> 58 s at n=12, and 213 s -> 201 s at n=24, where the wider pool was
slightly FASTER.

## The trade is coverage and gaps for specificity

Consistent in direction at every size: coverage rises by 1.1, 0.7 and 6.1
points, the worst hole shrinks or holds, and selectivity density falls by 12 to
42 per cent. Host sites rise by 40% at n=12 and 44% at n=24.

A user who runs `optimize` with no panel limits has asked for neither trade.
Wiring a refill in unconditionally would make that trade for them.

## Why, and the control that settles it

Two causes, and they compound.

**The shortlist is the most specific 2,000, so widening necessarily admits
worse candidates.** `max_primer` cuts on the step-2 ranking, which sorts
`ratio = bg_count / fg_count` ascending. Measured over the two groups:

| | median `bg_count` | median `fg_count` | implied ratio |
|---|---|---|---|
| the 2,000 shortlisted | 10 | 6 | ~1.7 |
| the 18,670 beyond it | 20 | 3 | ~6.7 |

The candidates a refill reaches bind the host twice as often and the target
half as often. This is not a defect in the shortlist; it is what the shortlist
is for.

**And nothing in Stage 1 resists them.** The set-cover greedy maximises
coverage bins and carries no specificity term: `optimize_greedy` accepts an
`objective` and no production caller supplies one, which is Known Issue 16.
Given ten times the candidates it finds more coverage, and pays for it in
specificity because specificity is not in the quantity it is maximising.

That second cause is also why the earlier `plan-pool` measurements do not
transfer. There the refill fired only when a size row could not satisfy its
constraints, and a row that already qualified never refilled -- so the panel
did not move wherever the shortlist was enough. `optimize` with no panel limits
has no notion of a row that fails, so an unconditional refill has nothing
holding it to the specificity the shortlist was selected for.

## What to do instead

**Do not give `optimize` an unconditional refill.** Either of these first:

1. Make Stage 1 specificity-aware, by passing the objective `optimize_greedy`
   already accepts (Known Issue 16). Then a wider universe is searched on the
   quantity the design is judged on rather than on coverage alone.
2. Trigger a refill only on an unmet CONFIGURED panel limit, which is what
   `plan-pool` does. `optimize` can now carry those limits, so this is
   well-defined: a user who sets `min_selectivity_density` and misses it has
   asked for a wider search, and a user who sets nothing has not.

Option 2 is the smaller change and matches the existing semantics. Option 1 is
the one that would make a wider universe worth searching at all.

## What this does not establish

One organism pair, one retention mode, one method, three panel sizes. The
`all_qc` universe of 491,836 was not measured here; the direction should be the
same and stronger, since its extra candidates are the ones the evenness gate
also rejected, but that is an expectation and not a measurement.

The reusable part is the same lesson this repository keeps relearning: a result
measured under constraints does not transfer to a path that has none. The prior
quoted against this work came from `plan-pool` runs with a specificity floor in
force, and the floor was doing the work being credited to the refill.
