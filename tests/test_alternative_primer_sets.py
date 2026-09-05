"""`max_sets` and `iterations` must change what a run produces.

Both are declared in `params.schema.json` and documented in CLAUDE.md under
Optimization -- "Parallel primer sets to build (default: 5)" and "Search
iterations (default: 8)" -- and neither reached anything. `max_sets` was read
only by `search_context.BFSConfig`, which has no callers anywhere in the
package, and no optimizer reads `config.max_iterations` at all; it is validated
by `OptimizerConfig.validate` and then never consulted.

They belong to one feature rather than two. The output format already
anticipates it: `step4_improved_df.csv` carries a `set_index` column that was
hardcoded to 0 on every row, and `optimize_step4` returns a `primer_sets` list
that always held exactly one entry.

So `max_sets` is implemented as "offer up to N distinct alternative sets", and
`iterations` bounds the search for them. Alternatives are produced by excluding
the primers already chosen and re-running selection, which is why they are
genuinely different sets rather than reorderings, and why the supply runs out on
a small candidate pool -- fewer than `max_sets` is a legitimate outcome, and the
first set is always the best one.

`iterations` deliberately does NOT bound the primary selection. Bounding that
would cap how many primers a run can choose, so `iterations: 8` would silently
truncate a 96-oligo panel.
"""

from dataclasses import replace

import pytest

from neoswga.core.base_optimizer import OptimizationResult, OptimizationStatus, PrimerSetMetrics

POOL = [
    "AAACCCGGGT",
    "ACCCGGGTTT",
    "AAGGCCTTAC",
    "GTAAGGCCTT",
    "TTGACCAGTC",
    "GACTGGTCAA",
    "CACACAACAT",
    "ATGTTGTGTG",
]


class _Optimizer:
    """Picks the first `target_size` candidates it is offered.

    Deterministic and order-sensitive, so excluding a set's primers genuinely
    yields a different answer -- which is what makes the alternatives real.
    """

    name = "stub"

    def __init__(self):
        self.calls = 0

    def optimize(self, candidates, target_size=None, **kwargs):
        self.calls += 1
        chosen = list(candidates)[: (target_size or 2)]
        return OptimizationResult(
            primers=tuple(chosen),
            score=1.0,
            metrics=replace(PrimerSetMetrics.empty(), fg_coverage=0.9),
            status=OptimizationStatus.SUCCESS if chosen else OptimizationStatus.NO_CONVERGENCE,
            optimizer_name="stub",
            iterations=1,
        )


def _alternatives(**kwargs):
    from neoswga.core.unified_optimizer import collect_alternative_sets

    optimizer = _Optimizer()
    primary = optimizer.optimize(POOL, 2)
    sets = collect_alternative_sets(
        primary=primary, optimizer=optimizer, candidates=POOL, target_size=2, **kwargs
    )
    return sets, optimizer


def test_one_set_is_the_default_behaviour():
    """`max_sets=1` must leave the run exactly as it was."""
    sets, optimizer = _alternatives(max_sets=1, max_iterations=8)
    assert len(sets) == 1
    assert optimizer.calls == 1, "no extra search should run for a single set"


@pytest.mark.parametrize("max_sets", [2, 3, 4])
def test_max_sets_controls_how_many_sets_come_back(max_sets):
    sets, _ = _alternatives(max_sets=max_sets, max_iterations=20)
    assert len(sets) == max_sets, f"asked for {max_sets} sets, got {len(sets)}"


def test_the_sets_are_distinct():
    sets, _ = _alternatives(max_sets=4, max_iterations=20)
    as_frozen = {frozenset(s) for s in sets}
    assert len(as_frozen) == len(sets), f"duplicate sets returned: {sets}"


def test_the_first_set_is_the_primary_result():
    """Alternatives are alternatives; the best set stays first."""
    sets, _ = _alternatives(max_sets=3, max_iterations=20)
    assert list(sets[0]) == POOL[:2]


def test_running_out_of_candidates_is_not_an_error():
    """Eight candidates and sets of two cannot yield ten distinct sets."""
    sets, _ = _alternatives(max_sets=10, max_iterations=50)
    assert 1 <= len(sets) <= 4
    assert all(sets)


def test_iterations_bounds_the_search_for_alternatives():
    """The parameter has to do something observable."""
    few, few_opt = _alternatives(max_sets=5, max_iterations=1)
    many, many_opt = _alternatives(max_sets=5, max_iterations=20)

    assert few_opt.calls < many_opt.calls, (
        "a tighter iteration budget did not reduce the number of searches run"
    )
    assert len(few) <= len(many)
