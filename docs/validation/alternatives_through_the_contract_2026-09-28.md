# Alternative primer sets through the shared search contract

2026-09-28. Task 6 of the valid-design plan: "route standard, planning, grid,
expansion, contraction and alternatives through the same contract." Measured on
the prepared Wolbachia design (wMel against *Drosophila*, 2,000 candidate
12-mers), requested size 12, a density floor of 20, one shared allowance of
4,000 evaluations.

## What changed

`collect_alternative_sets` called `optimizer.optimize` directly. Every other
path -- standard, planning, expansion, contraction -- goes through
`run_panel_search`. Alternatives now do too.

## The accounting was the point, and it was broken

| route | ledger evaluations spent by 5 alternatives | seconds |
|---|---|---|
| bare `optimize` | **0** | 1.8 |
| through the contract | **2,988** | 8.9 |

Zero is the number that matters. The ledger was bound around the alternative
search on 2026-09-27 and bounded nothing, because `optimize` consults no
objective: a run declaring `total_search_evaluations` could spend past it once
per alternative, and `collect_alternative_sets`'s own
`except SearchBudgetExhausted` clause was a handler for something that could not
happen. That correction is recorded in CLAUDE.md beside the claim it corrects.

Through the contract the work is counted, so a declared total now bounds the
whole run rather than the primary search alone.

## Panel quality is a wash

| set | bare: coverage / density | contract: coverage / density |
|---|---|---|
| 0 (primary) | 0.7200 / 28.78 | 0.7200 / 28.78 |
| 1 | 0.6710 / 20.80 | 0.6868 / 23.07 |
| 2 | 0.6632 / 21.48 | 0.6796 / 20.97 |
| 3 | 0.6997 / 24.61 | 0.6903 / 22.08 |
| 4 | 0.6518 / 17.05 | 0.6517 / 14.29 |

Two alternatives improve on both axes, one trades coverage for density, one is
slightly worse on density. **No claim is made that routing produces better
alternatives.** It produces honestly counted ones, assessed the way the primary
is assessed, at about 5x the wall clock -- a cost paid only when `max_sets` is
above 1, which is not the default.

## A separate defect this measurement exposed

**Set 4 violates the configured density floor on both routes**, at 17.05 and
14.29 against a floor of 20, and is offered to the user either way. Routing did
not fix that and was never going to: the repair could not resolve it, and
nothing filters or marks a violating alternative.

The primary is held to the limits -- `limit_violation_issue` turns an unmet
configured limit into a blocking finding that stops `export`. Alternatives are
not. `export_is_blocked(results_dir)` takes a directory and no set index, so the
findings it reads describe set 0 while `export --set 4` delivers a different
panel.

That is Known Issue 19's family reached from the other side: the commands were
taught to read ONE set, and the sets after the first are never assessed.

**The export half is fixed.** `export_is_blocked` now takes the set index and
refuses a non-zero set, naming why: nothing evaluated it. Verified end to end on
the plasmid example with five sets -- set 0 exports, set 1 is refused and writes
no files. `--allow-unqualified` overrides it, as it does every other refusal
there.

**The offering half is not**, and is a decision rather than a repair. Suppressing
a violating alternative, marking it, and assessing every set are three different
answers about what `max_sets` is for.

## What is still not routed

`primer_expansion._expand_hybrid` calls `optimizer.optimize` directly, with
expansion-specific arguments (`fixed_primers`, `apply_polymerase_multiplier`)
that the contract does not currently carry. Expansion's other path already goes
through `run_panel_search`.

## Reproducing

The measurement script is not kept; it builds one optimizer per route, runs the
primary through `run_panel_search`, then calls `collect_alternative_sets` with
`through_contract` False and True against the same seeded state, and reports the
ledger delta and per-set metrics. `through_contract=False` is still a parameter,
so the comparison can be re-run.
