# The host term as a weighted axis in swap Stage 2 (Known Issue 14)

Changed 26 October 2026, on `main` at `5dbe5dc` plus the change this document
records, during the optimizer audit that also produced
[`2026-10-26-clique-fixed-primer-dimers.md`](2026-10-26-clique-fixed-primer-dimers.md).
Reproduce with:

```
PYTHONPATH=. python docs/validation/stage2_background_weight_2026-10-26/probe.py
```

The regressions are `tests/test_swap_refinement.py`
(`test_off_by_default_keeps_the_higher_coverage_primer`,
`test_a_weighted_host_term_trades_coverage_for_a_lower_host_panel`,
`test_unmeasurable_host_coverage_falls_back_without_inventing_a_zero`) and
`tests/test_optimizer_config_reaches_optimizers.py`
(`test_stage2_weights_reach_the_inner_hybrid`).

## What Known Issue 14 said, and what changed

Known Issue 14: "Reading the background is not acting on it. `_network_refine`
carries the host term and `_swap_refine` reads the objective; neither carries
both." The swap Stage-2 score was a lexicographic tuple:

```
(-shortfall, coverage, -total_bg_sites)
```

so the host load (`total_bg_sites`) only ever decided between panels of **equal**
foreground coverage. A swap that shed a great deal of host binding at the cost of
a little foreground coverage was never taken, because coverage outranks the host
term absolutely.

This change adds a **new params.json key**, `stage2_background_weight`, that
turns the host term into a weighted axis instead:

| | |
|---|---|
| Swap score carries a weighted host axis | yes, in `refine_by_swaps._score` |
| Driven by | `stage2_background_weight` (params.json), None by default |
| Which methods | the swap refinement: `plan-pool`, `expand-primers`, and `hybrid` when `refinement_method` is `swap` |
| Off by default | yes -- weight 0/None is the old lexicographic tuple, bit for bit |

It is deliberately off by default, for the same reason the Tm term is
([`2026-10-26-wiring-the-tm-term.md`](2026-10-26-wiring-the-tm-term.md)): it moves
the delivered panel, and the trade it makes -- foreground coverage against host
coverage -- is **unmeasured against a reaction**. On the phi29 30 C regime the
host and target k-mer distributions overlap so heavily that no selection term
discriminates (Known Issue 17); the lever there is the reaction, not the score.
So this is a knob a user opts into with a background genome and a hotter
polymerase, not a default.

## How the weighted score works

```python
shortfall = -objective.shortfall(panel)
coverage = objective.coverage(panel)
metrics = objective.metrics(panel)
if background_weight > 0 and coverage is not None and metrics.bg_coverage is not None:
    return (shortfall, coverage - background_weight * metrics.bg_coverage, 0.0)
return (shortfall, coverage, -metrics.total_bg_sites)
```

Constraints still come first (`shortfall`), so a feasible panel is never traded
for an infeasible one. On the coverage axis, `coverage` and `bg_coverage` are
both 0..1 fractions, so `background_weight` is a genuine trade and not a
cross-scale fudge: a weight of 0.8 says "give up 0.8 points of host coverage to
gain 1 point of foreground coverage, or vice versa."

**Unknown is not zero.** The weighted form runs only when the objective can
measure both coverages. With no background genome `bg_coverage` is `None`, and
the score falls back to the lexicographic tuple rather than substituting a zero
host coverage -- which would both be a fabricated measurement and, done
naively, a `TypeError`. `test_unmeasurable_host_coverage_falls_back_without_inventing_a_zero`
pins that it neither raises nor invents.

## The discriminator (from the probe)

A panel starts as `('A',)` with candidate `B` available:

| primer | fg coverage | bg coverage |
|---|---|---|
| A | 0.60 | 0.90 |
| B | 0.50 | 0.10 |

- **weight 0.0 (default):** delivered `('A',)`. On pure coverage A (0.60) beats
  B (0.50); the host difference is invisible to the lexicographic score.
- **weight 0.8:** delivered `('B',)`. A scores `0.60 - 0.8*0.90 = -0.12`, B scores
  `0.50 - 0.8*0.10 = 0.42`; the weighted axis prefers B even though it covers
  less foreground.

## A latent Known Issue 8 closed on the way

Wiring `stage2_background_weight` required the `HybridBaseOptimizer` wrapper to
forward it to the inner `HybridOptimizer`. The wrapper builds the inner and sets
the two Stage-2 weights on it afterward (the inner's `__init__` is at its
line-length budget, so they cannot be added as constructor parameters):

```python
self._hybrid.stage2_tm_weight = getattr(self.config, "stage2_tm_weight", None)
self._hybrid.stage2_background_weight = getattr(self.config, "stage2_background_weight", None)
```

Adding the forwarding exposed that `stage2_tm_weight` had the same gap: the
`OptimizerConfig` carried it, but the wrapper never set it on the inner, so the
Tm term shipped in [`2026-10-26-wiring-the-tm-term.md`](2026-10-26-wiring-the-tm-term.md)
was **inert on every factory/CLI hybrid run** -- it only ever fired when a
`HybridOptimizer` was constructed directly (as that record's own measurement
script does). This is the textbook Known Issue 8 shape: both ends of the option
existed and the path between them did not.

`test_stage2_weights_reach_the_inner_hybrid` now builds a hybrid through the
real factory dispatch and asserts both weights reach `optimizer._hybrid`. With
the forwarding removed, both parametrisations fail
(`stage2_background_weight=0.5 ... is using None`), so the test guards the path,
not the two ends.

## Scope of the claim

No benefit is demonstrated. This records that the key exists, is off by default,
changes the swap pick in the direction it claims when switched on, and reaches
the code that decides. Whether a weighted host axis produces a better panel than
the lexicographic tie-break on a real design, under a reaction that can actually
discriminate, is not determined here and needs a measurement against a control.
