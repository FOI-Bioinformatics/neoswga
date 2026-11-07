# The Stage-2 Tm term, measured in a reaction that can discriminate

Measured 29 October 2026, on `main` at `e47b90b`. Phase 6.1 of
[the plan](../reports/PLAN_2026-10-28_pool_design_trustworthiness.md).
Reproduce with:

```
PYTHONPATH=. python docs/validation/stage2_terms_2026-10-29/measure.py \
    --params <work>/params.json --out results.json
```

## Why this was impossible until three days ago

`stage2_tm_weight` shipped on 2026-10-26 off by default, with "unmeasured
against a reaction" as the stated reason. That was not merely undone but
unreachable: `unified_optimizer` passed no polymerase to any optimizer, so every
CLI run was phi29 at 30 C
([`2026-10-26-hybrid-ran-every-reaction-as-phi29.md`](2026-10-26-hybrid-ran-every-reaction-as-phi29.md)),
and Known Issue 17 records that the phi29 pool cannot discriminate. The term's
own earlier measurement
([`2026-10-26-wiring-the-tm-term.md`](2026-10-26-wiring-the-tm-term.md)) was
made by constructing a `HybridOptimizer` directly, which is the only path on
which it worked at all.

**The regime here does discriminate**, which `filter` reports directly:

```
Pool discrimination: mean matched/mismatched occupancy 4.92,
0% of candidates saturated.
```

Against Known Issue 17's phi29 30 C case, where 86% of random 12-mers are bound
essentially all the time. equiphi29 at 43 C, k 10-12, 1.27 Mb target at 35% GC,
4 Mb host at 42% GC, 500 candidates.

Stage 2 has been deterministic since 2026-10-26, so one configuration gives one
panel and a single run per weight is the whole measurement rather than a sample
from a distribution.

## The measurement

`tm_dev` is mean |effective Tm − (reaction_temp + 5)| over the delivered panel:
the quantity the term's removal ranking is built to reduce.

**12 primers**

| weight | fg_cov | eff_cov | bg_cov | bg sites | sel_dens | tm_dev |
|---|---|---|---|---|---|---|
| off | 0.2785 | 0.0888 | 0.0060 | 3 | 4.25 | 7.01 |
| 0.25 | 0.2831 | 0.0913 | 0.0020 | 1 | 4.04 | 6.57 |
| 0.5 | 0.2781 | 0.0955 | 0.0020 | 1 | 4.10 | 6.33 |
| 1.0 | 0.2462 | 0.1012 | 0.0040 | 2 | 4.12 | **5.89** |

**24 primers**

| weight | fg_cov | eff_cov | bg_cov | bg sites | sel_dens | tm_dev |
|---|---|---|---|---|---|---|
| off | 0.4560 | 0.1293 | 0.0100 | 5 | 4.16 | 7.85 |
| 0.25 | 0.4656 | 0.1255 | 0.0080 | 4 | 4.08 | 7.88 |
| 0.5 | 0.4421 | 0.1279 | 0.0100 | 5 | 4.22 | 7.73 |
| 1.0 | 0.4111 | 0.1331 | 0.0080 | 4 | 4.18 | **7.40** |

## What holds on both panel sizes

**At full weight the term does what it claims, and charges for it.** Mean Tm
deviation falls on both sizes — 7.01 to 5.89 C (−16%) at twelve primers, 7.85 to
7.40 C (−5.7%) at twenty-four — and raw foreground coverage falls on both,
−3.2 pp and −4.5 pp. That is the trade the term exists to offer, and it is
real.

**Nothing holds at intermediate weights.** At 0.25 and 0.5 the figures are
non-monotone on the larger panel: Tm deviation is slightly WORSE at 0.25 than
off (7.88 against 7.85), and raw coverage is slightly better. Whatever happens
between zero and full weight is below the resolution of one pool.

## What does NOT replicate, and would have been the headline

At twelve primers, **effective** (occupancy-weighted) coverage rose
monotonically with the weight, 0.0888 to 0.1012, +14%. Read alone that is a
striking result: the term would be trading raw coverage for effective coverage,
buying the quantity that accounts for whether a primer is actually bound at the
reaction temperature, and so arguably not a trade at all but a win.

**It does not survive the second panel size.** At twenty-four primers the same
column is 0.1293 → 0.1255 → 0.1279 → 0.1331: non-monotone, dipping below the
off case at 0.25, and +2.9% net. One panel size would have produced a claim the
other refutes, which is why this project's own convention treats a single panel
size as "no benefit demonstrated".

So: **measured, benefit not demonstrated.** The term's effect on its own
objective is established at full weight; no coverage or specificity benefit is.

## Consequences

- `stage2_tm_weight` **stays off by default**, with its reason upgraded from
  "unmeasured against a reaction" to "measured in a discriminating reaction; it
  reduces Tm deviation at full weight and costs raw coverage, and no further
  benefit is demonstrated."
- Selectivity density is flat across every weight on both sizes (4.04–4.25), so
  the term neither helps nor harms the quantity a design is accepted on. The
  background-site counts are too small here (1–5) to read as anything.
- The 2026-10-26 record's figures stand for the pool they were taken on. They
  were measured under phi29 because that was the only reachable reaction, and
  the direction of its finding — lower Tm deviation, lower coverage — is
  reproduced here under a reaction that discriminates.

## Scope

One synthetic pool, one reaction, two panel sizes, four weights, one
refinement method. The genome is synthetic because no real reference is
reachable from this container (`tests/validation/genomes/` is absent and NCBI
is refused by the egress proxy), so site geometry is Poisson where a real
genome clusters. Nothing here is a wet-lab outcome: `eff_cov` is a modelled
occupancy weighting, not measured amplification.

The companion term could not be measured at all: see
[`2026-10-29-stage2-background-weight-does-not-run.md`](2026-10-29-stage2-background-weight-does-not-run.md).
