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
feeding the profile into Stage 2 would switch the term on for every CLI run. The
measurement below shows the term trades coverage for Tm by the amount of the
weight, and leaves specificity -- the quantity a design is accepted on --
unmeasured, so it is a trade a user should opt into, not a default.

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

> **Fixed later the same day**, after this and the floor record were written.
> `algebraic_connectivity` is now seeded, so Stage 2 is reproducible and the
> "measure around the noise" machinery below was no longer strictly necessary
> by the end of the day. The measurements here stand as taken; the fix and its
> own measurement are in
> [2026-10-26-stage2-determinism.md](2026-10-26-stage2-determinism.md). Of the
> two sources named just below, only the eigensolver turned out to matter;
> set-iteration order did not.

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

End-to-end delivered panels at HEAD (the wiring AND the floor change below),
medians over 7 repeats, `PYTHONHASHSEED=0`. **Medians reproduce exactly across
independent runs; the per-cell ranges do not**, because Stage 2 is
nondeterministic (above), so read the median as the result and the spread as
the noise. Reproduce with `measure.py`.

### Single-length pool, k=10, 150 candidates

| `stage2_tm_weight` | mean &#124;dTm&#124; (med) | coverage (med) |
|---|---|---|
| 0.00 (default) | 5.229 | 0.7400 |
| 0.10 | 5.229 | 0.7400 |
| **0.25** | **4.857** | 0.7400 |
| **0.30** | **4.857** | 0.7400 |
| 1.00 | **3.687** | **0.5000** |

### Mixed-length pool, k=9 to 12, 600 candidates

| `stage2_tm_weight` | mean &#124;dTm&#124; (med) | coverage (med) |
|---|---|---|
| 0.00 (default) | 6.481 | 0.8600 |
| 0.10 | **5.030** | 0.7600 |
| 0.25 | **5.030** | 0.7600 |
| 0.30 | **5.030** | 0.7600 |
| 1.00 | **2.813** | 0.7400 |

**The term lowers Tm deviation on both pools, and pays for it in coverage.**
The weight sets the exchange rate. On the mixed pool a light 0.10 already moves
mean deviation from 6.48 to 5.03 C for 10 points of coverage (0.86 to 0.76); a
full 1.0 reaches 2.81 C at 0.74. On the single pool the trade is milder until
1.0, which melts to 3.69 C at the cost of half the coverage (0.74 to 0.50). A
pure-Tm Stage 2 (weight 1.0) is a corner, not a recommendation: it is what
"optimise Tm and nothing else" delivers.

**This is HEAD, and it differs from the wiring commit alone.** At commit
`91d7c9a` -- the Tm term wired but the floor still on raw coverage -- every
nonzero weight on the mixed pool delivered the Stage-1 prefix unchanged (6.48
to 7.33 C, no benefit), because the floor discarded the trade. The next section
is the floor change that produced the mixed-pool numbers above, measured in
isolation so the nondeterminism cannot confound it.

A side effect worth recording: the term removes the nondeterminism it runs
under, collapsing 2 or 3 distinct panels to 1, because the Tm axis breaks the
ties that iteration order was breaking. That is not an argument for enabling it.

## Why the default is off

CLAUDE.md: "A capability with no demonstrated benefit ships off by default. Do
not switch one on without a new measurement."

After the floor change below, the term lowers mean Tm deviation on BOTH pool
shapes when it is weighted (mixed pool: 6.48 C at weight 0 to 2.81 at weight
1.0), and it does so by trading coverage, by the amount of the weight. That is a
trade a user should choose, not one a default should make for them: with
`--application` defaulting to `enrichment` at `tm_weight 0.25`, reading the
profile into this term would move every CLI optimize run off coverage toward Tm
without anyone asking. The effect on specificity -- the quantity a design is
actually accepted on -- is still unmeasured here (no background genome). So the
rule stands: demonstrated to move a panel, not demonstrated to improve a design,
therefore off by default. Concretely:

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

## The floor was the binding constraint, and was changed

**Done, same day.** The section above traced the mixed-length "no benefit" to
the Stage-1 floor: the removal loop produced a Tm-improved panel, the floor
compared it to the Stage-1 prefix on raw coverage alone, saw only the coverage
drop, and returned the prefix. So the floor, not the term, was discarding the
result.

The floor now keeps the prefix only when it PARETO-DOMINATES the refinement on
the axes Stage 2 optimized -- the accepted coverage (raw, or an attached
`pool_objective`'s), that objective's shortfall when one is attached, and mean
Tm fitness when the Tm term is active. A refinement that trades one axis for
another is incomparable and survives; one that is worse on every axis is still
floored. With a single axis (a run that declares neither an objective nor the
Tm term) dominance is the strict `>` it replaced, so that path is unchanged.
`neoswga/core/stage2_ranking.prefix_floors_refinement`.

### The floor decision, in isolation

A/B-ing "floor on" against "floor off" end to end compares two different loop
outputs, because Stage 2 is nondeterministic (below), so it measures the noise.
`docs/validation/tm_term_2026-10-26/measure_floor.py` instead takes the ONE
refinement the loop produced and asks what each floor rule does with that same
pair. Reproduce with `PYTHONHASHSEED=0`, which the script pins along with an
`algebraic_connectivity` seed; production stays unseeded.

Mixed-length pool, the refinement the loop produced vs the Stage-1 prefix, and
what each floor rule delivers (mean |dTm| to the 35 C optimum, lower better):

| `stage2_tm_weight` | refinement cov / dev | prefix cov / dev | old floor delivers | new floor delivers |
|---|---|---|---|---|
| 0.00 | 0.86 / 6.48 | 0.82 / 7.33 | 6.48 (refined) | 6.48 (refined) |
| 0.10 | 0.76 / 5.03 | 0.82 / 7.33 | **7.33 (prefix)** | **5.03 (refined)** |
| 0.25 | 0.76 / 5.03 | 0.82 / 7.33 | **7.33 (prefix)** | **5.03 (refined)** |
| 0.50 | 0.70 / 4.49 | 0.82 / 7.33 | 7.33 (prefix) | 4.49 (refined) |
| 1.00 | 0.74 / 2.81 | 0.82 / 7.33 | 7.33 (prefix) | 2.81 (refined) |

At weight 0 the two rules agree and keep the higher-coverage refinement: the
single coverage axis makes the new rule identical to the old. At every nonzero
weight the old rule returned the prefix and discarded the Tm gain; the new rule
keeps the refinement. The single-length pool behaves the same way, with the old
rule flooring the refinement at weight 0.50 and 1.0 where it costs coverage, and
the new rule keeping it.

**The coverage cost is real and the weight governs it.** On the mixed pool the
delivered panel moves from 0.86 coverage / 6.48 C deviation at weight 0 to 0.74
/ 2.81 at weight 1.0. A `stage2_tm_weight` run is choosing Tm over coverage by
the amount of the weight; the floor change stops the floor from silently
overriding that choice, it does not decide the trade. The term is still off by
default, so a default run is unaffected.

### What the floor change does not reach

The floor only stops discarding a refinement; it cannot improve one the loop
never produced. On a pool where the removal loop converges to the prefix on its
own, the floor change is a no-op. And it is confined to the one place in
`_network_refine` that was still on raw coverage: the swap path
(`_swap_refine`) already compared on the attached objective.

## What this still suggests, and does not do

- **Stage 2's nondeterminism should be fixed on its own.** DONE, same day:
  `algebraic_connectivity` is seeded and `optimize` is now reproducible. The
  set-iteration ordering this bullet also proposed proved unnecessary -- node
  order was already deterministic. [2026-10-26-stage2-determinism.md](2026-10-26-stage2-determinism.md).
- **The dominance floor relaxes the coverage-monotonicity guarantee when
  another axis is active.** `test_more_primers_do_not_reduce_coverage` still
  holds because it sets no Tm term and attaches no objective, so it stays on
  the single coverage axis. A run that sets `stage2_tm_weight` has opted into
  trading coverage and is no longer promised that more primers never lower
  coverage.
- **`uniformity_weight` remains inert** and is a separate problem: it is a
  property of a whole panel's gap distribution, not of one primer, so it fits
  neither the per-primer removal axis nor this per-panel floor axis directly.
