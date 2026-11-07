# The weighted host axis runs, and does not reliably reduce host coverage

Measured 2 November 2026, on `main` at `c3c1bc7` plus the change this document
records. Phase 2.2 of [the plan](../reports/PLAN_2026-11-02_known_issues.md),
against Known Issue 14. Reproduce with:

```
PYTHONPATH=. python docs/validation/host_weight_2026-11-02/measure.py \
    --params <work>/params.json --out results.json --size 12
```

Results in [`host_weight_2026-11-02/`](host_weight_2026-11-02/) with the pool's
configuration. The cost measurement that unblocked this is
[`2026-11-02-what-host-coverage-would-cost-the-search.md`](2026-11-02-what-host-coverage-would-cost-the-search.md).

## What changed to make the measurement possible

`stage2_background_weight` was REFUSED above zero from 2026-10-29, because the
weighted axis needs a host-coverage fraction and `PoolMetrics` carried five
fields without one. Phase 2.1 measured what the field costs: **+48% per
evaluation**, and **2.05x cheaper** than the full `compute_metrics` that
`panel_refinement` already falls back to whenever a host-coverage limit is set.
The positions were already being gathered for `total_bg_sites` and thrown away,
so the cost objection applied to the position lookups it was written about and
not to this field.

So `PoolMetrics` carries `bg_coverage`, the blanket refusal is gone, and the
key runs. `tests/test_pool_metrics_agree_with_the_full_evaluation.py` pins the
new field against `compute_metrics` like the other five, which is what stops
the focused and full evaluators drifting.

**One refusal remains**, and it is the part the old one was right about: a
weight above zero with **no background genome**. `bg_coverage` is `None` there,
not 0.0, and trading foreground coverage against a fabricated zero would make
every panel look equally specific. That replaces a silent fallback — the
previous behaviour returned a panel with the weight ignored and nothing said,
which is a key accepted and read by nothing.

## The measurement

equiphi29 at 43 C, k 10-12, 1.27 Mb target at 35% GC, 4 Mb host at 42% GC,
300 candidates, `refinement_method="swap"` in every arm. The control leaves the
key UNSET rather than setting 0.0, so it is the shipped path and not a weighted
axis that happens to weigh nothing. Stage 2 has been deterministic since
2026-10-26, so one run per weight is the whole measurement.

**12 primers**

| weight | fg_cov | **bg_cov** | bg sites | sel_dens | panel moved |
|---|---|---|---|---|---|
| unset | 0.4010 | 0.0220 | 11 | 4.795 | — |
| 0.25 | 0.3906 | **0.0160** | 8 | 4.809 | yes |
| 0.5 | 0.3906 | **0.0160** | 8 | 4.809 | yes |
| 1.0 | 0.4010 | 0.0220 | 11 | 4.795 | **no** |

**24 primers**

| weight | fg_cov | **bg_cov** | bg sites | sel_dens | panel moved |
|---|---|---|---|---|---|
| unset | 0.6114 | 0.0360 | 18 | 4.882 | — |
| 0.25 | 0.6100 | **0.0380** | 19 | 5.051 | yes |
| 0.5 | 0.6100 | **0.0380** | 19 | 5.051 | yes |
| 1.0 | 0.6018 | **0.0300** | 15 | 5.083 | yes |

## The finding: the direction does not replicate

Read the n=12 table alone and the term works: host coverage falls 0.0220 to
0.0160, a 27% reduction, and host sites 11 to 8, for 1.0 pp of foreground
coverage. That is the trade the key exists to offer.

**The second panel size contradicts it.** At 24 primers the same weights of
0.25 and 0.5 move host coverage the WRONG way — 0.0360 up to 0.0380, and host
sites 18 up to 19 — so the axis made its own objective worse. Only weight 1.0
reduces it there, and weight 1.0 is the setting that did nothing at all at
twelve primers, returning the control panel unchanged.

So across two sizes there is no weight that reliably reduces host coverage:

| weight | n=12 | n=24 |
|---|---|---|
| 0.25, 0.5 | host coverage down | host coverage **up** |
| 1.0 | **no change at all** | host coverage down |

This is a stronger negative than the Tm term's
([`2026-10-29-stage2-tm-weight-in-a-discriminating-reaction.md`](2026-10-29-stage2-tm-weight-in-a-discriminating-reaction.md)),
which at least moved its own axis in the claimed direction at full weight on
both sizes. Here the sign itself depends on the panel size.

**Why is not established.** A weighted axis changes which swaps the search
accepts, so it changes which local optimum it reaches; a term that reorders
candidates can land on a worse point for its own objective without any of the
arithmetic being wrong. Attributing it would need the search trace, which this
measurement does not collect.

**Selectivity density tells no consistent story either.** Flat at twelve
primers (4.795 to 4.809) and up at twenty-four (4.882 to 5.051/5.083). Since
that is the figure a design is accepted on, and the pool it moves on is the one
where host coverage went the wrong way, nothing here supports switching the
term on for specificity.

## Consequences

- **`stage2_background_weight` stays off by default**, with its reason upgraded
  from "cannot run" to **"runs, measured, and the direction on its own axis
  does not replicate across panel size"**. That is a better place for the key
  to sit than a blanket refusal: a user can now set it, and the schema says
  what the measurement found.
- **Known Issue 14 is still not closed, and the reason has changed.** It is no
  longer "the host term reaches no selection" — it reaches it now. It is that
  reaching it does not demonstrably help.
- **The host-coverage field is worth keeping regardless of the weight.** Any
  design setting `max_host_coverage`, `max_mean_gap` or `max_evenness` was on
  the full `compute_metrics` path; those runs can now use the focused evaluator
  at about half the per-evaluation cost. That is a speedup independent of this
  result.

## Scope of the claim

One synthetic pool, one reaction, two panel sizes, four weights, one refinement
method, one run per arm (which is the whole sample, since Stage 2 is
deterministic). The genome is synthetic because no real reference is reachable
from this container (`tests/validation/genomes/` is gitignored and absent, NCBI
is refused by the egress proxy), so site geometry is Poisson where a real
genome clusters — and the non-replication found here is precisely a
geometry-sensitive effect, so a real genome could behave differently in either
direction. Nothing here is a wet-lab outcome: `bg_cov` is a modelled union of
binding windows, not measured off-target amplification.
