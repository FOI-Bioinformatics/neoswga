# Looking past the first qualifying frontier costs specificity

Measured 22 September 2026 on a Prevotella melaninogenica design against human
chr21, phi29 at 30 C, 12-mers, effective coverage, selectivity density floor 10.
Inventory of 365,073 hard-QC candidates behind a 2,000-candidate opening
frontier.

## The question

`search_frontiers` broke on the first qualifying attempt. Refills therefore
existed only to rescue a search that had not qualified, so a run qualifying on
its opening frontier examined about 2,000 of 365,073 candidates and never
looked wider. Is that costing anything?

## The answer: it is buying specificity, not costing coverage

`improve_until_budget` keeps widening while the refill cap allows, retaining
the best-ranked incumbent.

| panel size | coverage | density | runtime | Jaccard |
|---|---|---|---|---|
| 8, first_feasible | 0.4041 | **29.18** | 9.2 s | |
| 8, improve_until_budget | 0.4071 | 25.64 | 36.4 s | 0.778 |
| 12, first_feasible | 0.4339 | **28.36** | 7.4 s | |
| 12, improve_until_budget | 0.4403 | 25.68 | 29.0 s | 0.500 |
| 16, first_feasible | 0.4635 | **24.72** | 10.7 s | |
| 16, improve_until_budget | 0.4647 | 18.23 | 54.5 s | 0.778 |

Every size tells the same story. Coverage rises by 0.3% to 1.5% relative.
Density falls by 9% to 26%. Runtime rises about fourfold. Half to a quarter of
the panel changes.

**So the default stays `first_feasible`.** Not because the difference is noise,
but because the trade is the wrong way round for a tool whose purpose is
selective amplification.

## Why, and it is already documented

This is not a surprise once the mechanism is read. `max_primer` cuts the
shortlist on `bg_count / fg_count` ascending, most specific first. So the
candidates a refill reaches are by construction the ones that bind the host
more and the target less, and CLAUDE.md already records the consequence:

> the candidates a refill reaches bind the host twice as often and the target
> half as often, and Stage 1's greedy has no specificity term to resist them

The measurement confirms that mechanism rather than discovering a new one. It
is the same shape as Known Issue 16, where supplying Stage 1 an objective
improves the metric it selects on and costs specificity every time.

## What would change the answer

A Stage 1 that weighs specificity. Until then, widening the frontier reaches
worse candidates and the greedy has nothing to stop it choosing them, so a
larger search is a worse search on the axis that matters.

## Two ways to get this measurement wrong

Both produce a confident number about a search that never widened, and both
were hit while taking it.

**The reaction fingerprint must be the one `filter` wrote under.** Reading
`params.json` is not enough and neither is calling `get_params`:

    params.json, fresh process     tm-2026-09-14:2b3536d4305777e8
    after get_params()             tm-2026-09-14:2b3536d4305777e8
    what filter actually wrote      tm-2026-09-14:2fb49b12c595c147

The difference is betaine 1.0 M that the params file never mentions. The
GC-adaptive strategy added it at run time and never writes back, which is
exactly what `effective_conditions` in the run manifest exists to record, and
what `export` and `report` already read in preference to the file. A mismatch
makes the inventory lookup refuse and the run fall back to the supplied list.

**`plan_pool` takes a source, not a list.** Handing it a list wraps it in a
`ListCandidateSource`, whose `advance()` returns False unconditionally. Every
refill reads zero and the run stops with `inventory_exhausted`, which looks
like a finding and is an artefact. The CLI passes `InventoryCandidateSource`.

`scripts/benchmarking/compare_delivered_panels.py` does both correctly and
prints the reachable universe and opening frontier, so a run that is about to
measure nothing says so first.

## Caveats

One organism pair, one reaction, three panel sizes, one density floor, one
seed. The direction is consistent across all three sizes and has a mechanism
behind it, which is why it is enough to leave a default alone; it would not be
enough to change one.
