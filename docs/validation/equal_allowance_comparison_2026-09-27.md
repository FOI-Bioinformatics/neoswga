# Panel searches under equal allowances

2026-09-27. Task 8 of the valid-design plan. Measured on the prepared
Wolbachia design (`examples/wolbachia_pool_design`, wMel against *Drosophila*,
2,000 candidate 12-mers), macOS arm64, Python 3.13.

Reproduce with `scripts/benchmarking/equal_allowance_comparison.py`, which
replaces `sequential_panel_search.py`.

## What was wrong with the comparison this replaces

The previous script conceded it in its own docstring: "Budgets are per stage,
not equal total compute". An arm running more stages was handed more compute,
so a difference in its panel could not be attributed to the search. It also ran
one unseeded configuration, and it built one proposal OUTSIDE the timed region
and shared it between arms.

What is equal now: one shared `SearchBudget` allowance per arm, the same integer
for every arm and seed, spent from the moment the arm begins, with proposal
generation inside it. References, inventory, constraints, reaction and requested
size are identical. What is not equal, and is reported rather than fixed: wall
clock, which is an outcome.

## The finding that only equal allowances could show

**At 400 evaluations the reduction arm looks useless. It is starved.**

| allowance | delivered | coverage | evaluations spent | stop reason |
|---|---|---|---|---|
| 400 | 12 | 0.7130 | 400 | budget exhausted |
| 2,000 | 3 | 0.5613 | 2,000 | budget exhausted |
| 8,000 | 3 | 0.5613 | 4,785 | finished inside |

Requested size 12, coverage target 0.50. The stage history says why: at 400,
`refinement` spent all 400 and `reduction` then ran for 0.19 ms with nothing
left to spend, changing nothing. Reporting that as "sequential reduction buys
nothing" would have been a statement about the allowance wearing the clothes of
a statement about the search.

A shared ledger starving a later stage is by design — the plan requires that a
later stage cannot get a fresh allowance — so the lesson is about reading the
result, not about the ledger. **An arm whose extra stage runs last is only
comparable once the allowance is large enough that the earlier stages do not
exhaust it.** `stop_reason` is what distinguishes the two cases, and it is why
the harness records it per stage.

This is also why the harness now records stage history. Without it a sweep
cannot tell a stage that ran and found nothing from a stage that never got to
run, and those support opposite conclusions.

## The comparison, at an allowance every arm finishes inside

Allowance 8,000, seeds 1-3, coverage target 0.50, density floor 20. Every arm
returned `stop_reason: null`, so none was truncated.

| requested | arm | delivered | effective coverage | density | host sites | seconds | evaluations |
|---|---|---|---|---|---|---|---|
| 6 | single_pass | 6 | 0.6495 | 36.82 | 90 | 4.2 | 2,551 |
| 6 | hybrid_single_pass | 6 | 0.6495 | 36.82 | 90 | 6.6 | 2,645 |
| 6 | sequential_reduction | **3** | 0.5613 | **61.49** | **45** | 7.6 | 6,270 |
| 12 | single_pass | 12 | 0.7200 | 28.78 | 148 | 3.1 | 1,012 |
| 12 | hybrid_single_pass | 12 | 0.7299 | 26.69 | 147 | 9.5 | 1,744 |
| 12 | sequential_reduction | **3** | 0.5613 | **61.49** | **45** | 6.5 | 4,785 |

**Reduction converges on the same three-oligo panel from either starting size.**
Identical coverage, density and host-site count from a request of 6 and of 12.
Against the twelve-primer single pass that is a quarter of the oligos, 2.1x the
selectivity density and 3.3x fewer host sites, for 0.159 fewer coverage points.

That is a trade, not a win. Which end of it is right depends on the application,
which is the decision `--application` exists for and which no measurement here
settles.

**`hybrid` against `dominating-set` is close to a wash.** At size 12 it buys
0.0099 coverage points for 2.09 density and costs 3.1x the runtime; at size 6
the two deliver the identical panel. That is consistent with what this
repository already records for the pair.

## What the numbers are, and are not

- **Coverage is a geometric proxy** at the realistic reach, occupancy-weighted.
  Nothing here has been calibrated against sequencing breadth. See
  [design_release_gates.md](design_release_gates.md).
- **The evaluation count** covers uncached evaluations of the shared objective
  plus the `clique` method's own scoring loop. It excludes one final assessment
  per stage, deliberately, so reporting cannot consume a search's allowance.
- **Every arm was deterministic across three seeds**: one distinct panel each.
  So the medians above are exact values and the seed range is zero. Three
  identical answers are evidence of determinism, not an uncertainty estimate,
  and the harness reports `distinct_panels` so a reader is not invited to read a
  range that does not exist.
- **The position cache is shared across arms.** It holds binding positions,
  which are reference data identical for every arm, and building it per arm
  would dominate every timing while changing no answer. Each arm builds a fresh
  optimizer, so per-optimizer memoisation starts cold.
- **Memory is one figure for the process**, 1,321 to 1,466 MB across runs, and
  deliberately not per arm. `ru_maxrss` is a high-water mark that never falls,
  so the first arm absorbs the shared cache build: an early version of this
  harness reported 1,112 MB for the first arm and 54 MB for the next doing the
  same work. A per-arm figure needs one process per arm.

## A defect the sweep then found in the search itself

The allowances above straddle the point where `reduction` begins. Running them
showed the stage delivering 12 primers at 1,020 and at 1,060 alike, which did
not fit: at 1,060 it had 48 evaluations to spend and the first removal needs
about 12.

It was discarding them. `reduce_result` removes one primer per iteration and
each smaller panel is a valid incumbent -- it met the coverage target and
violated nothing, which is why it was accepted -- but the shared allowance can
run out inside that loop, since `objective.coverage` is the budgeted call. The
exception propagated out of the function, and `run_panel_search`'s `execute`
returned the PRE-STAGE incumbent because it has no access to the stage's partial
state.

Measured on a controlled eight-primer panel with a twenty-evaluation allowance:
two removals accepted, both thrown away. On the real design, before and after:

| allowance | delivered before | delivered after |
|---|---|---|
| 1,020 | 12 | 12 (nothing found in 8 evaluations) |
| 1,040 | 12 | **10** |
| 1,060 | 12 | **8** |
| 1,200 | 12 | **3** |

So the answer that needed 2,000 evaluations now arrives at 1,200, because
progress is no longer thrown away and re-attempted. **No default run changes**:
`total_search_evaluations` is None by default, so there is no allowance to
exhaust.

The stop reason now distinguishes `shared_search_allowance` from
`evaluation_budget`. One says the run is out of compute, the other that this
loop reached its own per-stage bound while the run had more to give, and a
reader deciding whether to raise a limit needs to know which.

## Not covered

The task asks for more than this, and the rest is not done:

- **The full host reference.** This is wMel against *Drosophila* (144 Mb), not
  against hg38.
- **Synthetic adversarial geometry**, multiple oligo lengths, and the
  fixed-total against per-oligo concentration policies.
- **More than two requested sizes**, and more than two coverage targets. The
  0.70 target was also run and changed nothing, because it gates removals inside
  reduction only and at 400 evaluations reduction never ran.

## A defect this harness had, found by running it

`--target` was accepted, recorded in the output JSON, and read by nothing, so
the first sweep ran every arm at `OptimizationRequest`'s default of 0.7 while
the saved record said otherwise. That is the inert-option class this repository
keeps a ratchet for, in the script written to measure the search.
`tests/test_the_benchmark_compares_equal_allowances.py` now pins that the target
reaches the request, along with the rest of the method: one shared allowance, a
required rather than defaulted allowance, several seeds, proposal generation
inside the timed region, and the warm cache disclosed.
