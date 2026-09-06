"""When ensemble methods tie, the winner must not be decided by flag order.

`_run_ensemble` picked the winner with `max()` over a dict keyed by method name.
`max()` returns the FIRST key attaining the maximum and the dict is built by
`results[m] = res` inside `for m in methods`, so the tie-break was
`--ensemble-methods` order.

Ties are the common case rather than a corner. Measured on the plasmid example,
hybrid, dominating-set and background-aware return the identical set -- the
latter two because background-aware wraps the same HybridOptimizer -- so their
metrics are identical and their scores equal to the last digit in every
application profile. Reordering the flag moved the winner every time:

    --ensemble-methods hybrid dominating-set background-aware   -> hybrid
    --ensemble-methods dominating-set background-aware hybrid   -> dominating-set
    --ensemble-methods background-aware hybrid dominating-set   -> background-aware

all three rows scoring 0.9003. CLAUDE.md claims the ensemble is
order-independent; `_reseed` makes the RUNS order-independent and does nothing
about the CHOICE.

On that example only provenance moved, because the tied methods returned the
same primers. The tie is on the composite score, though, not on the set, so two
methods returning DIFFERENT primers can tie and then flag order picks the pool
that gets delivered.
"""

from dataclasses import replace

import pytest

from neoswga.core.base_optimizer import OptimizationResult, PrimerSetMetrics


def _result(name, primers, coverage=0.9):
    return OptimizationResult(
        primers=tuple(primers),
        score=1.0,
        metrics=replace(PrimerSetMetrics.empty(), fg_coverage=coverage),
        status=__import__(
            "neoswga.core.base_optimizer", fromlist=["OptimizationStatus"]
        ).OptimizationStatus.SUCCESS,
        optimizer_name=name,
        iterations=1,
    )


TIED = {
    "hybrid": _result("hybrid", ["AAACCCGGGT", "ACCCGGGTTT", "AAGGCCTTAC"]),
    "dominating-set": _result("dominating-set", ["AAACCCGGGT", "ACCCGGGTTT"]),
    "background-aware": _result("background-aware", ["AAACCCGGGT", "ACCCGGGTTT", "AAGGCCTTAC"]),
}


def _pick(order):
    from neoswga.core.unified_optimizer import _select_ensemble_winner

    results = {name: TIED[name] for name in order}
    return _select_ensemble_winner(results, application="balanced")


@pytest.mark.parametrize(
    "order",
    [
        ["hybrid", "dominating-set", "background-aware"],
        ["dominating-set", "background-aware", "hybrid"],
        ["background-aware", "hybrid", "dominating-set"],
    ],
)
def test_the_winner_does_not_move_with_the_flag_order(order):
    """Same three tied methods, three orders, one answer."""
    assert _pick(order) == "dominating-set", (
        f"order {order} chose {_pick(order)}; the winner is following "
        "--ensemble-methods order rather than a stated rule"
    )


def test_a_tie_is_broken_towards_the_smaller_set():
    """Fewer oligos at equal quality is cheaper to synthesise and simpler to
    run, so it is the defensible tie-break."""
    assert _pick(["hybrid", "dominating-set"]) == "dominating-set"


def test_a_genuine_winner_still_wins():
    """The tie-break must only apply to ties."""
    from neoswga.core.unified_optimizer import _select_ensemble_winner

    results = {
        "hybrid": _result("hybrid", ["AAACCCGGGT"], coverage=0.99),
        "dominating-set": _result("dominating-set", ["AAACCCGGGT", "ACCCGGGTTT"], coverage=0.10),
    }
    assert _select_ensemble_winner(results, application="balanced") == "hybrid"
