# Getting ahead on spacing, and what the evidence rules out

Written 2026-09-18, following
[pool_selection_audit_2026-09-18.md](pool_selection_audit_2026-09-18.md), which
found NeoSWGA behind swga 1.0 and swga 2.0 on foreground and background site
spacing.

**The two obvious ways to close that gap are both already refuted, one by
earlier work here and one by a measurement in this note.** What is left is not a
modelling gap. NeoSWGA already computes more spacing information than any of the
three published tools and discards or ignores almost all of it.

Reproduce the new measurement with
`scripts/benchmarking/published_gap_thresholds.py`.

## Ruled out: a fitted spacing weight

[published_primer_sets.md](published_primer_sets.md) settled this with two
benchmarks that disagree.

| Metric | Clarke *M. tuberculosis*, n=12 | Dwivedi-Yu *Prevotella*, n=6 |
|---|---|---|
| mean binding distance | separates winners | does not |
| worst gap | does not | separates winners |
| Gini evenness | does not, runs slightly backwards | separates winners |

Clarke pre-filtered the candidate pool for even binding, so evenness had little
residual variance left to explain the outcome and density dominated. Dwivedi-Yu
did not, so coverage holes dominated. The conclusion recorded there is the
useful one: **what predicts is whichever property is currently limiting**, and a
weight tuned for one dataset is wrong for the other. A `max_gap` term was built
and measured against the Prevotella benchmark, changed nothing, and was not
shipped.

That also disposes of copying swga 2.0. Its five fitted coefficients assume one
regime for every design. Two benchmarks here say the regime moves.

## Ruled out: a hard gap constraint derived from the reach

This was the next idea and it is worse than the first. If a hole wider than
twice the coverage reach cannot be bridged, a panel carrying one is defective
whatever regime it is in, and a constraint says so without fitting anything.

The reach is calibrated to 3.0 to 6.2 kb ([reach_calibration.md](reach_calibration.md)).
Every published set in both fixtures carries `fg_max_distance`, so the claim is
directly testable on 18 sets with wet-lab outcomes, 7 of them flagged successful
by their own papers.

| Reach | Threshold at 2x | Winners above it | Others above it | Separates |
|---|---|---|---|---|
| 3,000 | 6,000 | 7 of 7 | 11 of 11 | no |
| 4,000 | 8,000 | 7 of 7 | 11 of 11 | no |
| 6,200 | 12,400 | 7 of 7 | 11 of 11 | no |
| 10,000 | 20,000 | 3 of 7 | 10 of 11 | no |
| 20,000 | 40,000 | 0 of 7 | 4 of 11 | yes |
| 70,000 | 140,000 | 0 of 7 | 0 of 11 | no |

**At every reach-derived threshold the constraint rejects all 18 published
panels, the wet-lab winners included.** The only threshold that separates is
40 kb, which implies a 20 kb reach, outside the calibrated band and fitted to
these 18 points, which is the thing this project has twice declined to do. Two
of the four sets it catches are `MtbUneven` and `MtbSparse`, chosen by Clarke
et al. as deliberate negative controls, so the genuine separation is 2 of 11.

## What the working panels actually look like

The measurement is more useful as a description than as a test.

| | Worst hole |
|---|---|
| Tightest of all 18 published sets | 13,876 bp, and it is not a winner |
| The four most effective *M. tuberculosis* sets | 14,949 to 18,721 bp |
| The three successful *Prevotella* sets | 31,300 to 33,300 bp |
| The two deliberate negative controls | 92,071 and 99,986 bp |

**Panels that enrich 96 to 120 fold carry 31 to 33 kb coverage holes.** At the
3 kb selection reach that is a hole ten times the per-primer reach, in a
published success. Enrichment is a ratio and does not require closing holes, so
hole-closing is not what makes an SWGA panel work. This is worth holding against
COATswga's `target_coverage` of 0.95 and against this project's own application
profiles, which set coverage targets of 70 to 95 percent: no published outcome
supports treating coverage completeness as the binding criterion.

It also explains, after the fact, why the `max_gap` term built here changed
nothing, and why an occupancy-weighted version of the same statistic should not
be expected to do better. Weighting a quantity that does not separate winners
gives a better-founded quantity that still does not separate them.

## What is actually missing

Not modelling. NeoSWGA computes, per panel: `mean_gap`, `max_gap`, `gap_gini`,
`gap_entropy`, `strand_alternation_score`, `strand_coverage_ratio` and
`bg_coverage`. `PositionCache.compute_strand_alternation_stats` additionally
computes `strand_alternation_gap_mean`, `strand_alternation_gap_max` and
`longest_same_strand_run`. That is more spacing information than swga 1.0,
swga 2.0 and COATswga hold between them.

Four things are wrong with what happens to it.

**Three of the five strand quantities are thrown away at the call site.**
`base_optimizer._compute_metrics` reads `strand_alternation_score` and
`strand_coverage_ratio` and drops the rest. `strand_alternation_gap_max` is the
mechanistically closest quantity to the thing that matters, since exponential
amplification needs convergent sites within reach rather than merely many sites,
and it is computed and discarded on every run.

**The strand statistics cover one foreground genome and never the background.**
The loop `break`s after the first `fg_prefix`, so a multi-target design reports
one target's strand structure, and the host's is never computed although the
same method would compute it.

**`bg_coverage` has no reader**, which Known Issue 18 records. It is the only
computed quantity that sees background site POSITION, and it is the direct form
of what swga 2.0's `off_gap_gini` proxies: clustered host sites leave most of
the host unreachable, which is why an uneven host scores well in their fitted
model. NeoSWGA can measure the thing rather than the proxy.

**The constraint framework reaches one command.** `PoolConstraints` carries
`min_selectivity_density` and `max_background_sites`, is enforced inside the
search through `shortfall`, and is repaired by a beam. Only `plan-pool`
constructs one. This is the one piece of architecture no published tool has, and
it is where being ahead is available.

## What would put NeoSWGA ahead, by axis

Ordered by evidence, not by appeal.

**1. A regime diagnostic, which nothing else offers.** Two benchmarks support
exactly one conclusion: the limiting property varies by design. No tool tells a
user which property limits theirs. NeoSWGA has every quantity and has
`shortfall`, which already measures per-constraint distance in comparable units.
Reporting which criterion binds a delivered panel, and how far the others are
from binding, requires no fitting, no threshold and no new science. It is the
direct implementation of the only thing the outcome data supports.

**2. Generalise `PoolConstraints` and make it reachable from `optimize`.** The
spacing quantities above become constrainable rather than scored, so a user in a
hole-limited regime constrains the hole and a user in a density-limited regime
does not. swga 1.0 has hard constraints with fixed magic numbers; swga 2.0 has
fitted weights that assume one regime; this would have neither limitation.

**3. Stop discarding the strand quantities, and compute them for the host.**
Pure wiring. `strand_alternation_gap_max` on the background is the convergent
pair term swga 2.0 approximates with `within_mean_gap_ratio`.

**4. Wire the condition-aware free-energy dimer model that already exists.**
`dimer.is_dimer_thermodynamic` takes a delta-G threshold and a
`ReactionConditions`, is tested, and has no production caller: the shipped
screen is longest complementary run. COATswga is genuinely ahead here, using
PrimerROC free energy, and its model is NOT condition-aware. Wiring this one
would put NeoSWGA ahead rather than level. It is on the O(n^2) hot path, so it
needs the lazy screen and measurement before it becomes a default.

**5. A worst-target term for multi-genome designs.** `per_target_coverage` is
already computed and used only as a post-hoc warning. No published tool has
multi-target designs at all, so this extends a lead rather than closing a gap.

**6. Occupancy-weighted spacing, with expectations set low.** A gap bounded by a
site occupied 10 percent of the time is not really bounded, and no other tool
can express that. It is the most novel item here and the least supported: the
unweighted statistic does not separate winners on 18 sets, so weighting it is
unlikely to. Worth computing as a diagnostic under item 1, not as a scoring
term.

**Not on this list: fitting a set-level score.** The position in
[tool_comparison.md](tool_comparison.md) stands and is now better supported.
Eighteen sets across two benchmarks that disagree is not a basis for weights.

## A correction to the audit

The audit's comparison table says dimer-freedom is hard "in `clique` only;
penalised in the others". That understates the default. With
`allow_dimer_relaxation` false, which is the default, selection STOPS rather
than admitting a violating pair, and
`tests/test_delivered_panel_honours_the_dimer_limit.py` pins that every
delivered pair is within the limit. `clique` differs by finding the MAXIMUM
dimer-free set rather than a greedy one. NeoSWGA is at parity with the two tools
that guarantee it, not behind them.

## What this does not establish

- The 18 sets carry PUBLISHED statistics, not statistics recomputed from the
  genomes, so `fg_max_distance` is each paper's own figure under its own
  definition. The two papers agree on the definition of `fg_bg_ratio` and that
  quantity separates winners in opposite directions in the two datasets, which
  is a warning about how far to trust cross-dataset comparison here.
- Winner status is the fixtures' own flag: `most_effective` from Clarke Table 2
  and `successful` from Dwivedi-Yu. The Prevotella flag includes a 13x set
  alongside the 96x and 120x sets.
- Seven winners against eleven others cannot support any threshold, including
  the 40 kb one that happens to separate them. The table is evidence AGAINST a
  reach-derived constraint, not evidence for any other.
- Ruling out a hard gap constraint at the panel level says nothing about the
  per-primer `max_gini` gate, which is a different rule applied to a different
  object and was re-derived against delivered coverage in 2026-09-10.
- Items 1 to 6 are unmeasured. None has been shown to change a delivered panel,
  and this project's recent record on plausible search improvements is four
  attempts that did not.
- Coverage, occupancy and selectivity remain modelled site geometry rather than
  measured amplification. See
  [evidence_matrix_2026-09-15.md](evidence_matrix_2026-09-15.md).
