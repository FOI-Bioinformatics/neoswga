# Making the search see the specificity floor

Measured 2026-09-17 and 2026-09-18 on the real Wolbachia pair. Follows the
defect [Phase 5](retention_changes_no_delivered_panel_2026-09-17.md) recorded:
a design given a selectivity floor spends its candidates on coverage and
delivers a panel that misses the floor.

Two things shipped, one did not, and the one that did not is the more useful
finding.

## What shipped: the objective never reached the stage that refines

`plan_pool` attaches its objective to the optimizer it is handed. Every
command-line path is handed a wrapper -- `OptimizerFactory` returns
`HybridBaseOptimizer` or `BackgroundAwareBaseOptimizer` -- and each delegates
the search to an inner `HybridOptimizer`. `_swap_refine` is a method of the
INNER one, so `refine_hybrid_stage2` read the attribute off an object nobody had
set it on.

Measured through a real design: the refinement ran once and received None. So
Stage 2 refined on raw covered bases while the row was accepted on
occupancy-weighted coverage under a specificity floor, which is exactly the
two-rule split Phase 2 set out to remove.

| Delivered density on a failing row | Before | After |
|---|---|---|
| Floors 65 to 79, panel size 12 | 28.78 | 42.62 |

Two tests covered this and neither could see it. One asserted by AST that
`plan_pool` assigns an attribute of that name; the other asserted by source text
that the refinement reads one. Both ends existed and the path did not.
`attach_search_config` now sets it on the wrapper and the delegate, and
`tests/test_the_objective_reaches_the_stage_that_refines.py` drives a real
factory-built optimizer under both methods and asserts the object the refinement
receives is the object the row is accepted on.

## What the exact accounting shows

`occupancy.weighted_site_load` is a SUM of per-primer terms that depend only on
the primer, so a panel's loads are additive and the density floor rearranges
into a linear condition:

    density(panel) >= D   <=>   SUM_i (f_i * bg_len - D * b_i * fg_len) >= 0

That makes the achievable ceiling computable rather than searchable. On the
2,000-candidate shortlist, the best 12-primer panel reaches a density of
**79.807** -- and so does the best from the 20,670 post-Gini universe, to three
decimals, which independently confirms Phase 5's finding that retention adds no
achievable specificity here.

**But that ceiling ignores coverage, and the panel achieving it has coverage
0.4042, below the 0.5 target.** So it is not attainable for a qualifying row,
and an earlier draft of this work quoted 19.7 points of headroom that partly
did not exist.

The honest headroom comes from a constructive existence proof. Restricting to
the candidates that meet the floor on their own and taking the best-coverage
twelve:

| Floor | Individually affordable | Density | Coverage | Both met |
|---|---|---|---|---|
| 60 | 19 | 75.593 | 0.6100 | yes |
| 65 | 16 | 78.178 | 0.5702 | yes |
| 70 | 12 | 79.477 | 0.5706 | yes |
| 75 | 8 | too few for a 12-panel | | |
| 79 | 6 | too few | | |

Floors of 65 and 70 are satisfiable with coverage above target, by a rule with
no parameters. **The shipped search reaches 60.11 and fails at 65.** That gap is
real, and it is about the search rather than the pool.

## What did not ship, and why

Three Stage 1 rules were built on the accounting above. None shipped, because
none reliably improved a delivered panel.

| Rule | Density at floors 65 / 70 / 75 / 79 |
|---|---|
| Stage 2 connected only | 42.62 / 42.62 / 42.62 / 42.62 |
| Plus feasible-at-every-step | 46.64 / 41.56 / 39.53 / 24.97 |
| Plus feasibility-still-reachable | 45.85 / 45.77 / 41.86 / 39.53 |
| Plus preferring affordable candidates | 34.12 / 39.53 / 37.54 / 42.20 |

The first rule was actively harmful at strict floors: requiring the panel to
satisfy the floor at every intermediate step forbids taking a candidate with
negative slack however much later headroom would restore it, so the greedy spent
its headroom early and could then afford nothing good. Lookahead fixed that and
is better than Stage 2 alone at two floors and worse at two. Preferring
affordable candidates, which is what the existence proof does, was worse than
either.

**The structural reason is that Stage 1 is not the stage that picks the panel.**
It selects `max(final_count + 8, final_count * 1.67)` primers -- 20 for a
requested 12 -- and Stage 2 narrows them. Making Stage 1 build a feasible
20-primer panel does not make the delivered 12 feasible, and the accounting's
`slots_left` refers to the wrong panel size throughout. That is why a rule which
works offline over a fixed set fails inside the greedy.

The other half is the selection criterion. The existence proof picks by TOTAL
coverage contribution within the affordable set; the greedy picks by MARGINAL
gain, and after a handful of affordable candidates the marginal gains collapse,
so it falls through to looser candidates and dilutes the pool it hands Stage 2.

So `SelectivityBudget` lives in `scripts/benchmarking/selectivity_budget.py`
rather than in the package. It is a diagnostic: it produced every number above
and none of it steers a search. Keeping it in `neoswga.core` unused would be the
defect this audit is named after, and shipping a rule with the results in that
table would be the other one.

## What the next attempt should do

The evidence points at Stage 2's choice of twelve from twenty, not at Stage 1's
choice of twenty. Specifically:

- Constrain the SIZE the accounting is done at. The feasibility test must be
  over the delivered panel, which means either Stage 2 carrying the budget or
  Stage 1 selecting at the final size.
- Select within the affordable set by total contribution, not marginal gain.
  The parameter-free rule above reaches 78.178 at floor 65 where the search
  reaches nothing.
- Measure against the existence proof, not against the density-only ceiling.
  79.807 is unattainable under the coverage target; 78.178 at floor 65 is not.

## What this does not establish

- One pair, one panel size, one floor sweep, single runs. The existence proof is
  a lower bound on what is achievable, not a claim about what is optimal.
- The four rules were each measured once per floor. The differences between the
  middle two are within a few density points and no spread was measured, so
  "better at two floors and worse at two" should not be read as a reliable
  ranking.
- The shipped improvement, 28.78 to 42.62, is one quantity on rows that are
  rejected either way. It does not make a design qualify that did not.
- Nothing here changes a delivered panel on a row that already qualifies: at
  floor 60 the panel and its density are unchanged.
- Coverage and selectivity remain modelled site geometry, not measured
  amplification. See [evidence_matrix_2026-09-15.md](evidence_matrix_2026-09-15.md).
