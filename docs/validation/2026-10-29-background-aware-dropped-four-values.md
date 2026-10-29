# `background-aware` dropped four configured values, and a ratchet now checks the seam

Found and fixed 29 October 2026, on `main` at `334443c`, by the ratchet Phase
3.1 of [the plan](../reports/PLAN_2026-10-28_pool_design_trustworthiness.md)
asked for. The ratchet is
`tests/test_no_adapter_drops_an_inner_parameter.py`; the behavioural pin is
`test_stage2_weights_reach_the_inner_hybrid`, now parameterised over both
wrappers.

## What the ratchet does

For every class that constructs an inner optimizer, it diffs the inner's
`__init__` parameters against what the constructing class passes as a keyword,
sets with `setattr`, or assigns on the instance. A parameter supplied by none of
those routes silently takes the inner's default, and no configuration can reach
it.

That check is mechanical and takes no judgement, which is the argument for
having a machine do it: of the five instances of this defect found by hand in
October, it would have caught two directly.

## What it found immediately

`BackgroundAwareBaseOptimizer` builds the same `HybridOptimizer` class as the
`hybrid` wrapper does, from a second, separately maintained argument list. Four
configured values reached its `OptimizerConfig` and never reached the engine:

| Value | What was inert |
|---|---|
| `max_dimer_dg` | the free-energy dimer floor |
| `stage1_objective_width` | Stage 1's objective-scored set cover |
| `stage2_tm_weight` | the Stage-2 Tm term |
| `stage2_background_weight` | the Stage-2 weighted host axis |

Measured through the real factory dispatch, before the fix, with all four set on
the config:

| method | `max_dimer_dg` | Stage-1 width | `stage2_tm_weight` | `stage2_background_weight` |
|---|---|---|---|---|
| `hybrid` | -6.0 | 32 | 1.0 | 0.5 |
| `background-aware` | **None** | **None** | **None** | **None** |

Two of those four are keys this repository shipped three days earlier, with
their own validation records
([`2026-10-26-stage2-background-weight.md`](2026-10-26-stage2-background-weight.md)).
Both were fixed for the `hybrid` wrapper and left broken here, so "off by
default" was true of `background-aware` in a way nothing recorded: the keys could
not be switched on at all.

**This is the second time these two lists have diverged.** The comment above the
fix records the first: "six values arrived at `hybrid` and were dropped here",
costing two thirds of a user-approved equiphi29 pool. A seam that has drifted
twice is one a test should hold, which is what this phase was for.

## The fix

The three that are inner constructor parameters are passed; `stage2_background_weight`
is not one, so it is set on the instance the way `HybridBaseOptimizer` sets it.
Both wrappers now agree, and a default run is unchanged — all four values are
`None` unless params.json supplies them, verified on a plain config.

## What the ratchet does NOT catch

Stated here and in the test, because Known Issue 8 is never recorded as closed:

- **A real default passed instead of the configured value.** This rule sees
  `polymerase=kwargs.get("polymerase", "phi29")` as "passed", and that was the
  widest of the five October instances
  ([`2026-10-26-hybrid-ran-every-reaction-as-phi29.md`](2026-10-26-hybrid-ran-every-reaction-as-phi29.md)).
  Pinned separately by `test_the_configured_reaction_reaches_the_hybrid_inner`.
- **Resolution inside the inner.** `resolve_reaction` is not part of an
  adapter's argument list, so removing it leaves this ratchet green. The plan's
  acceptance criterion claimed otherwise and was wrong about that; the reaction
  test above is what holds it.
- **Whether the inner then uses the value.** Reaching the constructor is not the
  same as entering a decision, which is the distinction Known Issue 14 is about.

## A third seam, recorded rather than fixed

The ratchet's completeness test — which fails if a class constructs an inner
optimizer and is not listed — found three further seams beyond the four method
adapters. `HybridOptimizer`'s own Stage 1 wraps the dominating-set engine and is
clean. The other two carry documented `ALLOWED_UNREACHED` entries, and a
staleness test fails on an entry that has since been wired, so the list can only
shrink:

- **`PrimerExpander` drops 20 of `HybridOptimizer`'s 32 parameters**, including
  `polymerase`, `reaction_temp`, `min_tm`, `max_tm`, `max_dimer_dg`,
  `stage1_objective_width` and `stage2_tm_weight`, and 4 of 10 on its
  dominating-set engine including `dimer_temp`. So `expand-primers` builds a
  lean optimizer that does not know the reaction, the Tm window or the dimer
  floor. **This is a live command and a real finding.** It is not fixed here
  because what a full forwarding changes about a delivered expansion is not
  measured, and Known Issue 14 records this path's own history of reading the
  background without acting on it. It wants its own phase.
- `ImprovedPipeline` drops 10 of 16 on a `NetworkOptimizer`, but it is reached
  only from `core/validation.py` and no design command, so a dropped parameter
  there moves no delivered panel.

## Verification

Three negative checks, each confirmed to fire:

1. Removing `HybridBaseOptimizer`'s `setattr` loop → fails naming
   `stage2_tm_weight`.
2. Reverting this fix → fails naming `max_dimer_dg`, `stage1_objective_width`
   and `stage2_tm_weight`.
3. An `ALLOWED_UNREACHED` entry for a parameter that is forwarded → the
   staleness test and the no-holes-on-method-adapters test both fire.

Full suite identical to baseline: the same 60 failures, all environmental in
this container. Zero added.
