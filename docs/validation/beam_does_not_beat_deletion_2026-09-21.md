# The beam finds no smaller panel on five instances, including the one that started this

Measured 21 September 2026. **This supersedes
[smallest_pool_on_wolbachia_2026-09-21.md](smallest_pool_on_wolbachia_2026-09-21.md),
written earlier the same day, whose central claim does not reproduce.**

## What was claimed

That deletion stops one oligo short of the smallest qualifying panel on the
shipped Wolbachia design: `--minimize-primers` returns 12 where a beam finds 11
at higher coverage. The conclusion drawn was that switching
`--minimize-primers` to the beam "would deliver 11 oligos where it currently
delivers 12, on every design".

## What five instances say

One method throughout: greedy at n=12 through the normal pipeline, deletion at
the greedy's own `fg_coverage`, then `panel_beam.beam_search` at width 4 over
the delivered panel plus the first 100 candidates in the pool's own order. All
against human chr21 except the last row.

| target | genome | greedy coverage | deletion | beam best at n=11 | saving |
|---|---|---|---|---|---|
| E. coli | 4.64 Mb | 0.2551 | 12 | 0.2497 | none |
| S. aureus | 2.82 Mb | 0.1985 | 12 | 0.1859 | none |
| M. tuberculosis | 4.41 Mb | 0.1661 | 12 | 0.1619 | none |
| Wolbachia wMel | 1.27 Mb | 0.6652 | 12 | 0.6553 | none |
| wMel against *Drosophila* | 1.27 Mb | 0.7334 | 12 | 0.7239 | none |

Deletion removes nothing on any of the five, which the earlier document also
found. What does not hold is the other half: the beam finds nothing smaller
either. In every case its best eleven-primer panel falls short of the greedy's
own coverage.

The last row is the original instance, the same pool and the same position
indexes, with the greedy panel and the beam both produced inside one script so
that neither can differ from the other by setup.

## Why the earlier measurement was wrong

Its baseline does not match this design. It quotes the greedy's coverage as
0.7529. Two independent records put that panel at **0.7334**:
[occupancy_and_discrimination_2026-09-19.md](occupancy_and_discrimination_2026-09-19.md)
and Known Issue 17 in CLAUDE.md, both from before the beam work. Re-running the
design here returns 0.7334 to the last digit.

So the panel the beam was asked to beat was not the panel this design
delivers. Whatever that starting set was, it was a weaker twelve than the
greedy produces, and beating a weaker twelve with an eleven is a much easier
problem. What exactly the earlier run did differently is not recovered: it was
measured inline with no script kept, which is why it cannot be re-examined.

That is the lesson worth keeping. **A measurement whose baseline disagrees with
the project's own recorded figure for the same quantity is reporting on a
different object**, and the check that would have caught it costs nothing:
compare the baseline against what the repository already says before drawing a
conclusion from the delta.

## A plausible explanation, tested and refuted

The three bacteria sit at 17 to 26 percent coverage while Wolbachia was
reported at 75, which suggested the beam needs the redundancy that only appears
near saturation. Holding genome, pool and panel size fixed and widening only
the reach put S. aureus at 0.8017 coverage, comparable to Wolbachia:

| S. aureus, n=12 | reach | greedy coverage | beam best at n=11 |
|---|---|---|---|
| shipped | 3 kb | 0.1985 | 0.1859 |
| widened | 12 kb | 0.6067 | not run |
| widened | 20 kb | 0.8017 | 0.7773 |

Still no saving at 0.80. The coverage regime does not explain the result,
because there was nothing to explain.

## What still stands

Nothing here rehabilitates deletion as an optimal procedure, and two findings
from the earlier work are untouched:

- **Deletion is a local optimum, not the smallest pool.**
  `tests/test_smallest_pool_search.py` proves it by enumerating every subset of
  a four-candidate fixture where one candidate covers what two others cover
  between them. Deletion stops at 3; the beam finds 2. That is a true statement
  about the algorithm.
- **`beam_search` is not reachable from `optimize`.** It is called only from
  `pool_planner`, so `plan-pool` can escape the local optimum and
  `optimize --minimize-primers` cannot.

What changes is the reason for not wiring it. It was recorded as a pending
decision backed by a measured saving. There is no measured saving, so it
returns to what it was before: a capability with no demonstrated benefit, the
same resolution Known Issues 11 and 16 reached.

## Caveats

The beam pool is bounded at roughly 110 oligos for cost, so each row is "no
saving found within this bound" rather than "no saving exists". A wider pool
can only help the beam. That asymmetry is what makes the earlier positive
suspect rather than these negatives weak: it claimed a saving within the same
bound that five runs cannot now find.

Five instances, one panel size, one beam width. Whether some target has the
dominance structure the fixture demonstrates is still unknown; what is now
measured is that four real bacterial and intracellular genomes do not.
