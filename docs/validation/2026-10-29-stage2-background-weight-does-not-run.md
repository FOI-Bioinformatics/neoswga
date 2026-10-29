# `stage2_background_weight` never ran on a real design, and now refuses

Found 29 October 2026, on `main` at `e47b90b`, by the first attempt to MEASURE
the key (Phase 6.1 of
[the plan](../reports/PLAN_2026-10-28_pool_design_trustworthiness.md)). It
corrects [`2026-10-26-stage2-background-weight.md`](2026-10-26-stage2-background-weight.md),
which is this repository's own record of shipping it three days earlier.

## What happened

Setting `stage2_background_weight` on a real design crashed inside Stage 2:

```
Hybrid optimization failed: 'PoolMetrics' object has no attribute 'bg_coverage'
```

The weighted axis reads `objective.metrics(panel).bg_coverage`. No production
objective supplies that. `PoolMetrics` carries exactly five fields —
`fg_coverage`, `effective_fg_coverage`, `selectivity_density`,
`total_bg_sites`, `max_gap` — and `pool_metrics.py`'s own module docstring says
why:

> The returned object carries only those five fields on purpose. A caller that
> reaches for a sixth gets an `AttributeError` rather than a default that would
> be indistinguishable from a measurement.

Background coverage is absent **deliberately**: that module exists because
`compute_metrics` was 27.6 s of the 29.7 s a repair-heavy `plan_pool` spent, and
computing a host-coverage fraction per candidate panel is part of that cost. So
the quantity the weighted axis needs is the quantity the search objective was
built to stop computing.

## Why the tests did not catch it

The three behavioural tests shipped with the key used a stub objective built on
`SimpleNamespace`, and that stub supplied `bg_coverage` because **the author
wrote the interface the implementation wanted rather than reading the one the
objective provides**. `SimpleNamespace` accepts any attribute, so the stub could
not have failed. The feature was tested against a fiction.

The regression test added here uses a metrics object with `__slots__` naming the
five real fields, so an absent sixth is expressible and the test would fail if
the refusal were removed.

## What the key actually did between 2026-10-26 and 2026-10-29

Nothing, and then it crashed. On `hybrid`'s default `network` refinement the
swap score is not used at all, so the key was inert; on `plan-pool`,
`expand-primers`, or `hybrid` with `refinement_method="swap"`, any value above
zero raised `AttributeError` mid-search. The default — unset, the lexicographic
host tie-break — was never affected, which is why nothing failed until someone
switched the key on.

The earlier record's claim that a weight above zero "trades foreground coverage
against host coverage" described the stub's behaviour, not the product's. Its
discriminator table is a measurement of the stub.

## What changed here

The key now **refuses** before any search, rather than crashing inside one:

```
Invalid design request field 'stage2_background_weight': it requires a
host-coverage fraction and PoolMetrics does not carry one: the search
objective reads five quantities and background coverage is deliberately not
among them, because computing it dominated the search loop. Unset the key to
keep the lexicographic host tie-break, which is the shipped default.
```

This is the rule the project already states for this situation — a required
calculation that fails raises a `DesignError` subclass and is not caught and
substituted. The alternatives were both worse:

- **Add `bg_coverage` to `PoolMetrics`.** That pays, inside the search loop, the
  cost the module was created to remove. For an off-by-default term with no
  demonstrated benefit, that is the wrong trade.
- **Re-base the axis on `total_bg_sites`**, which the objective does carry. But
  that is a count against a 0..1 coverage fraction, and the original record
  claimed in as many words that the weight was "a genuine trade and not a
  cross-scale fudge". Making it exactly that fudge to keep the feature alive
  would be worse than admitting it does not work.

An objective that genuinely carries `bg_coverage` still gets the weighted axis;
the code checks rather than assumes.

## Scope, and what is still open

Known Issue 14 is **not** closed by this. The gap it names — neither Stage 2
carrying both the shared objective and the host term — is where it was before
2026-10-26, except that the attempt is now documented and the broken path
refuses instead of crashing. Closing it properly needs a host quantity the
search can afford, which is a measurement about cost, not a weight.

Phase 6.1 could therefore only measure the other term; see
[`2026-10-29-stage2-tm-weight-in-a-discriminating-reaction.md`](2026-10-29-stage2-tm-weight-in-a-discriminating-reaction.md).
