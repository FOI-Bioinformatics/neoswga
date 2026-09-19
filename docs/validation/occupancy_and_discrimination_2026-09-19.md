# The pool cannot discriminate, and the candidate filter is not the fix

Measured 2026-09-19 on the Wolbachia/Drosophila pair and on random k-mers,
revisiting Known Issue 17. **One half of that entry's diagnosis does not
survive measurement, and the remedy it implies makes delivered panels worse.**

Occupancy is the fraction of the time a site is bound at the reaction
temperature. Discrimination is the ratio of matched to single-mismatch
occupancy: how much better a primer binds its true site than a near-miss,
which is what specificity ultimately rests on. The mismatch penalty is 4 C
(`occupancy.default_mismatch_penalty`).

## What does not survive: the Tm floor is not padding the pool

Known Issue 17 says the default floors sit 5-17 C below the reaction
temperature, "so the pool is padded with primers whose occupancy at the
reaction temperature is 0.002 to 0.13".

phi29's `primer_tm_range` floor is indeed 10 C below its 30 C reaction
temperature. But at k = 12 there is almost nothing down there to admit. Over
40,000 random 12-mers at phi29 30 C:

| Tm band | n | share | mean occupancy | mean discrimination |
|---|---|---|---|---|
| 20-25 | 8 | 0.02% | 0.081 | 5.599 |
| 25-30 | 130 | 0.3% | 0.365 | 4.181 |
| 30-35 | 1,148 | 2.9% | 0.785 | 2.070 |
| 35-40 | 4,165 | 10.4% | 0.967 | 1.163 |
| 40-50 | 20,304 | 50.8% | 0.998 | 1.009 |
| 50-70 | 14,235 | 35.6% | 1.000 | 1.000 |

The band between the default floor and the reaction temperature holds 8 of
40,000 candidates. Moving that floor changes essentially nothing at k = 12.

## What does survive, and it is severe: nothing excludes the saturated

The other half of the entry is right. 86% of random 12-mers sit at or above
0.998 occupancy, where a 4 C mismatch penalty moves occupancy so little that
discrimination is 1.01 or less. On the real Wolbachia shortlist: **65% of the
2,000 candidates are above 0.99 occupancy and mean discrimination is 1.11.**

A primer that binds its true site 1.11 times better than a single-mismatch
site is not selecting for anything. Specificity in such a pool is a property
of where sites happen to fall, not of binding.

## But gating on occupancy delivers a worse panel

The obvious remedy is to exclude the saturated candidates. Measured, with the
candidate list authoritative (the inventory moved aside, so the pool is the
pool), `optimize` at n=12 on the same data and seed:

| pool | size | coverage | selectivity density | host sites | panel discrimination |
|---|---|---|---|---|---|
| no cap | 2,000 | 0.7334 | **25.62** | 149 | 1.098 |
| occupancy < 0.99 | 700 | 0.4530 | 6.59 | 174 | 1.189 |
| occupancy < 0.95 | 260 | 0.4397 | 6.60 | **237** | 1.400 |

The gate does what it says -- the delivered panel's discrimination rises -- and
it is a clear loss on everything that matters. Coverage falls by 29 points,
selectivity density falls by three quarters, and host binding RISES by 59%.

Discrimination lives in a tail too small to build a panel from: at k = 12 and
30 C, candidates with discrimination above 2 are 1.6% of the space. Forcing
selection into that tail costs far more than it returns.

**A methodological note, because it nearly produced the opposite conclusion.**
The first run of this experiment swapped `step3_df.csv` and appeared to show
the cap IMPROVING density to 44.28. It did not: `open_source_or_list` prefers
the candidate inventory over the supplied list, so the swapped file only set
the frontier SIZE and the run searched the inventory as usual. The apparent
gain was the smaller frontier, not the occupancy cap. None of the delivered
primers were in the capped pool, which is how it was caught.

## What does move the pool is the reaction

Same 20,000 random 12-mers, four chemistries:

| reaction | occupancy > 0.99 | mean discrimination | discrimination > 2 |
|---|---|---|---|
| phi29 30 C | 85.6% | 1.065 | 1.6% |
| phi29 30 C + DMSO 5% | 75.6% | 1.153 | 4.5% |
| phi29 30 C + DMSO 10% + betaine 1.5 M | 53.5% | 1.374 | 11.3% |
| equiphi29 42 C | 24.2% | **2.190** | **36.0%** |

Saturation is a phi29-at-30-C problem, not a filtering problem. A warmer
polymerase moves 36% of the same sequence space into real discrimination; an
additive moves phi29 part of the way without changing the enzyme. Nothing
about the candidates changed in any row of that table.

That is the point Known Issue 17's second paragraph already made -- "an
additive is the only lever that moves a candidate along the thermodynamic axis
without changing its composition" -- and this measurement gives it numbers and
rules out the alternative.

Across reaction temperature the relationship is monotonic over phi29's
permitted range, on 300 random 12-mers:

| reaction temp | mean discrimination | fraction saturated |
|---|---|---|
| 25 C | 1.012 | 0.97 |
| 30 C | 1.078 | 0.83 |
| 35 C | 1.320 | 0.57 |
| 40 C | 1.939 | 0.32 |

## What ships

`occupancy.discrimination_profile` measures the regime and
`log_discrimination_profile` reports it at the end of `filter`. On the
Wolbachia design it prints the profile and then warns, naming the lever that
was measured to work and the one that was measured not to:

```
Pool discrimination: mean matched/mismatched occupancy 1.11, 65% of candidates saturated.
WARNING: This pool cannot discriminate: ... The lever that moves this is the
reaction, not the candidate filter: a warmer polymerase (equiphi29 at 42 C
measures 2.19 against phi29's 1.07 at 30 C) or an additive (DMSO 10% with
betaine 1.5 M reaches 1.37). Narrowing the Tm window instead was measured and
delivers a worse panel.
```

`DISCRIMINATION_FLOOR` is 1.5, between the two measured regimes: phi29 at 30 C
is 1.065 and equiphi29 at 42 C is 2.190. No candidate is filtered and no panel
moves.

## What this does not establish

One organism pair for the panel measurement, k = 12 throughout, and a 4 C
mismatch penalty that is a modelling choice rather than a measurement. The
discrimination figures scale with that penalty, so they are a comparison
between chemistries rather than an absolute prediction of off-target rate.

It does not test whether a discrimination TERM in selection -- as opposed to a
gate on the candidate pool -- would help. That is a different experiment and
the one worth doing next.
