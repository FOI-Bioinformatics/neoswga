# Wiring the Tm term into hybrid's Stage 2

Measured 26 October 2026, on `main` at `44ed69b` plus the change this document
records. Reproduce with:

```
PYTHONPATH=. python docs/validation/tm_term_2026-10-26/measure.py
```

## What was asked and what shipped

[`selection_weights_are_inert_2026-09-21.md`](selection_weights_are_inert_2026-09-21.md)
measured that `tm_weight` and `uniformity_weight` reach nothing on the default
method: `HybridOptimizer` hands both to a `NetworkOptimizer` that nothing reads
back, and Stage 2 is `_network_refine`, a method on the class itself with no Tm
term. CLAUDE.md records the consequence: wiring them "moves every delivered
panel and is a decision, not a repair."

The decision was taken. What shipped is the Tm term, live, driven by a **new
key** and **off by default**:

| | |
|---|---|
| Stage 2 carries a Tm axis | yes, in `_rank_removal_candidates` |
| Driven by | `stage2_tm_weight` (params.json), None by default |
| `tm_weight` from `--application` | still inert on `hybrid`, still warned about |
| `uniformity_weight` | still inert on `hybrid`, still warned about |

The two keys are deliberately separate, and that separation is the whole of the
decision. All five `--application` profiles set `tm_weight` between 0.10 and
0.30, and `--application` defaults to `enrichment` on the optimize parser, so
feeding the profile into Stage 2 would have switched the term on for every CLI
run. The measurement below does not support that.

## How it works

Stage 2 is a greedy removal loop: from the Stage-1 set it drops one primer at a
time until the target size, ranking candidates per step on rank-normalised
axes. It had three: network score, unique-bin coverage cost, and host load. It
now has four.

The Tm axis scores each primer by `selection_weights.tm_fitness`, a Gaussian on
the distance between the primer's effective Tm and `reaction_temp + 5 C`, which
halves 3.7 C away. The axis is inverted relative to the fitness, because every
axis in that loop reads "higher is better to drop": the primer furthest from the
reaction's optimum scores 1.

Three properties were made to hold deliberately:

- **A zero weight is bit-identical.** The weight budget is arranged so that at
  `tm_share == 0` the other weights are exactly what they were. Checked against
  a reimplementation of the previous formula over 4,000 randomised candidate
  sets, both background-aware and not: zero disagreements. Pinned by
  `test_a_zero_tm_weight_cannot_rank`.
- **One definition of the Tm preference.** `NetworkOptimizer._calculate_tm_score`
  now delegates to the same `tm_fitness`, and the effective-Tm lookup moved to
  `selection_weights.effective_tm`, which `BaseOptimizer._estimate_tm` also
  delegates to. `HybridOptimizer` does not inherit `BaseOptimizer` -- its MRO is
  `(HybridOptimizer, ThermoScreenMixin, object)` -- so without that extraction
  the two stages would have held two copies of the same quantity, which is the
  failure this codebase has produced most often.
- **No temperature means no ranking.** `tm_fitness` returns 1.0 uniformly when
  there is no reaction temperature, so the term cannot rank, and
  `_network_refine` logs that the configured weight is not being applied. An
  absent temperature is not evidence that a primer melts badly. An explicit
  `reaction_temp` wins and the polymerase preset supplies the default, matching
  how `_thermo_criteria` resolves the Tm window.

Two inputs had to be retained on the instance to make any of this possible.
Neither was: `reaction_temp` was a constructor parameter forwarded to the unread
`NetworkOptimizer` and never stored, which is why its docstring has read "for
Stage-2 Tm scoring" since long before any Stage-2 code could consult it. Known
Issue 8 in its usual shape, and the reason the first version of this measurement
reported noise as signal -- see below.

## Setup

| | |
|---|---|
| target | `examples/plasmid_example/pcDNA.fasta`, 6,157 bp, one record |
| pool | the real step-2 shortlist, `step2_df.csv.original`, sampled to 150 per k at seed 20261026 |
| reaction | phi29, 30 C, so the optimum is 35 C |
| panel | 6 primers, `coverage_reach` 500 |
| repeats | 7 per weight |

No k-mer counter is needed: position indexes are written straight from the
sequence, as `tests/test_hybrid_optimizer_run.py` does.

## Stage 2 is not deterministic, and that had to be measured first

The first version of this measurement ran each weight once and reported
differences of 0.1 to 0.6 C as the term's effect. They were noise. Two
independent sources, both pre-existing and neither introduced here:

- **Set and dict iteration order.** The removal loop iterates sets of nodes and
  primers, and iteration order varies with `PYTHONHASHSEED`, which breaks ties
  in the ranking differently per process. Three runs at a fixed seed agree;
  three runs at different seeds do not.
- **`nx.algebraic_connectivity(subgraph, method="tracemin_lu")` is called with
  no `seed`.** Its starting basis is random, so the connectivity axis varies
  per call within one process.

Over 7 repeats at the shipped default the baseline returns **2 to 3 distinct
panels** on each pool, varying by run. Any claim about a
selection change smaller than that spread is unsupportable. This is a
reproducibility defect in its own right, on a tool that records a seed in its
run manifest, and it is **not fixed here**: fixing it changes delivered panels
for everyone and is its own decision with its own measurement.

## Result

Medians over 7 repeats. **The medians reproduce exactly across independent
runs of the script; the ranges and the jaccard column do not**, because the
baseline they are measured against is itself one of 2 or 3 panels (see above).
Two full runs of the script gave identical medians in every cell and different
spreads, so read the medians as the result and the spread as an indication of
the noise. `jaccard` is the best overlap with any panel that run's baseline
produced, which is why it moves when the baseline does.

### Single-length pool, k=10, 150 candidates

| `stage2_tm_weight` | distinct panels | mean &#124;dTm&#124; | range | coverage | range | jaccard |
|---|---|---|---|---|---|---|
| 0.00 (default) | 2-3 | 5.229 | 5.23-5.86 | 0.7400 | 0.7400-0.7600 | 1.000 |
| 0.10 | 1 | 5.229 | 5.23-5.23 | 0.7400 | 0.7400-0.7400 | 1.000 |
| **0.25** | 1 | **4.857** | 4.86-4.86 | 0.7400 | 0.7400-0.7400 | 0.714 |
| **0.30** | 1 | **4.857** | 4.86-4.86 | 0.7400 | 0.7400-0.7400 | 0.714 |
| 1.00 | 1 | 5.615 | 5.62-5.62 | **0.6800** | 0.6800-0.6800 | 0.714 |

**The term does its job here.** At 0.25 and 0.30 it improves mean Tm deviation
from 5.23 to 4.86 C, a 7% reduction, and the coverage it delivers sits inside
the baseline's own range. One primer in six changes. At 0.10 the term is present
but too weak to outrank the other axes, and the panel is the baseline's.

At 1.00 it is harmful: with the network and coverage weights driven to zero the
removal is Tm-only, coverage falls from 0.74 to 0.68, and mean deviation ends up
*worse* than the baseline's because the panel it reaches is then floored
differently. A pure-Tm Stage 2 is not a useful setting.

### Mixed-length pool, k=9 to 12, 600 candidates

| `stage2_tm_weight` | distinct panels | mean &#124;dTm&#124; | range | coverage | range | jaccard |
|---|---|---|---|---|---|---|
| 0.00 (default) | 2-3 | 6.481 | 5.67-7.33 | 0.8600 | 0.8200-0.8600 | 1.000 |
| 0.10 | 1 | 7.334 | 7.33-7.33 | 0.8200 | 0.8200-0.8200 | 0.09-1.00 |
| 0.25 | 1 | 7.334 | 7.33-7.33 | 0.8200 | 0.8200-0.8200 | 0.09-1.00 |
| 0.30 | 1 | 7.334 | 7.33-7.33 | 0.8200 | 0.8200-0.8200 | 0.09-1.00 |
| 1.00 | 1 | 7.334 | 7.33-7.33 | 0.8200 | 0.8200-0.8200 | 0.09-1.00 |

**No benefit, at any weight.** Every nonzero weight converges on the same panel,
which is the Stage-1 prefix: mean deviation 7.33 at or beyond the worst end of
the baseline's own range, and coverage at the bottom of it. The jaccard column
spans 0.09 to 1.00 across runs precisely because the baseline sometimes lands on
that prefix by itself and sometimes does not.

The mechanism, traced by instrumenting the refinement:

| `stage2_tm_weight` 0.25, mixed pool | mean &#124;dTm&#124; | coverage |
|---|---|---|
| what the removal loop produced | 7.33 | 0.8200 |
| the Stage-1 prefix it was floored against | 7.33 | 0.8200 |
| what a Tm-optimal panel of this size would be | **2.81** | **0.7400** |

Stage 2 floors its output at the Stage-1 greedy prefix of the same size and
**compares on coverage alone**. The Tm term trades coverage for Tm; that trips
the floor; and the floor discards the refinement including its Tm gain. The
Tm-optimal panel exists in this pool and is 3.7 C better, and the floor would
reject it at every weight because it covers 12 points less.

So on a mixed-length pool the term cannot deliver as currently wired. It is not
that the Tm preference is wrong; it is that a Tm term inside a stage gated on
coverage can only help where it is nearly coverage-neutral, which is the
single-length case.

A side effect worth recording: the term removes the nondeterminism it runs
under, collapsing 2 or 3 distinct panels to 1, because the Tm axis breaks the
ties that iteration order was breaking. That is not an argument for enabling it.

## Why the default is off

CLAUDE.md: "A capability with no demonstrated benefit ships off by default. Do
not switch one on without a new measurement."

The benefit is demonstrated on one pool shape and absent on the other, and the
absent one is the shape the shipped `min_k=6`/`max_k=12` range produces. With
`--application` defaulting to `enrichment` at `tm_weight 0.25`, routing the
profile into this term would have changed every CLI optimize run on the
strength of a single-length result. So:

- `stage2_tm_weight` defaults to None, and None is bit-identical to the
  pre-change Stage 2.
- The profile's `tm_weight` keeps feeding `NetworkOptimizer`, where it was
  always live, and the inert-weight warning now names `stage2_tm_weight` as the
  supported route rather than only saying the weight does nothing.
- `tests/test_design_options_have_effect.py::test_stage2_tm_weight_defaults_to_off`
  pins the default on the delivered panel, not on the attribute.

## Size of the claim

One target, one 6,157 bp plasmid, one panel size, one reaction, no background
genome, two pool shapes, 7 repeats each. That is enough to say the term is live
and to show the floor interaction, and it is **not** enough to recommend a
weight. Per this project's own standard, one pool and one panel size is "no
benefit demonstrated", not "no benefit".

What would settle it: the same grid on a real target with a background, at
several panel sizes, where occupancy-weighted selectivity contributes and the
Tm term's effect on specificity rather than only on Tm deviation can be seen.
That needs a reference genome and a k-mer counter, neither of which is in this
container.

## What this suggests next, and does not do

- **The floor is the binding constraint, not the term.** A Stage-2 floor that
  compared on the accepted objective rather than on raw coverage would let a
  Tm-improving, slightly-coverage-costing refinement survive. That is a change
  to the floor's contract and wants its own measurement.
- **Stage 2's nondeterminism should be fixed on its own.** Passing a seed to
  `algebraic_connectivity` and ordering the loop's set iterations would make
  `optimize` reproducible. It changes delivered panels, so it is a decision.
- **`uniformity_weight` remains inert** and is a separate problem: it is a
  property of a whole panel's gap distribution, not of one primer, so it does
  not fit the per-primer removal axis the Tm term uses.
