"""The smallest qualifying pool, checked against every subset there is.

Task 7 of the 2026-09-21 valid-design plan. A design that stops at the first
feasible panel reports a size, not the smallest size, and the two are only the
same when the first thing found happens to be optimal.

The oracle here enumerates all 2**n subsets of a tiny candidate set inside this
file and finds the smallest that qualifies. It does not call the production
solver to compute the expected answer, which is what makes the comparison
evidence rather than a restatement. Ten candidates is 1,024 subsets, which is
free; the point of keeping it small is that exhaustive is then possible at all.

The fixture is a set-cover instance with an exact structure, so the answer can
be reasoned about rather than trusted:

    P1 covers bins {0,1,2}
    P2 covers bins {3,4,5}
    P3 covers bins {6,7,8,9}
    P4 covers bins {0,1,2,3,4,5}

Full coverage of the ten bins needs P3, plus either {P1,P2} or {P4}. So
{P1,P2,P3} qualifies at size 3 and {P4,P3} qualifies at size 2, and **no
single deletion from {P1,P2,P3} qualifies**: dropping any one of them loses
bins nothing else in that panel covers. A search that only deletes is stuck at
3 while the answer is 2, and reaching 2 needs a swap before a deletion.

**The gap is constructed, and did not reproduce on either real instance
available here** (measured 2026-09-21). On `examples/plasmid_example` one
primer covers the target completely at 3 kb reach, so there is nothing to
reduce. On a 300 kb random sequence with 60 candidates, deletion stopped at
the requested 20 and a beam over the same pool found nothing smaller at the
same coverage. Random sequence rarely produces the structure this fixture
has, where one candidate dominates the union of two others.

So the deletion-only local optimum is real and the beam does escape it, both
shown below. What is NOT established is that it costs anything on a candidate
pool anyone would design from, and that is why `optimize --minimize-primers`
is left on the deletion path rather than switched to the beam. Same resolution
as Known Issues 11 and 16: a mechanism that fires without a demonstrated
benefit does not ship on by default.
"""

import itertools
from dataclasses import replace
from types import SimpleNamespace

import pytest

from neoswga.core.base_optimizer import OptimizationResult, OptimizationStatus, OptimizerConfig
from neoswga.core.base_optimizer import PrimerSetMetrics

# Distinct 12-mers, so nothing is rejected by a self-dimer or length rule.
P1, P2, P3, P4 = "ACGGACGGACGG", "AGGAGGAGGAGG", "ACAACAACAACA", "AGCAGCAGCAGC"

BINS = {
    P1: frozenset({0, 1, 2}),
    P2: frozenset({3, 4, 5}),
    P3: frozenset({6, 7, 8, 9}),
    P4: frozenset({0, 1, 2, 3, 4, 5}),
}
TOTAL_BINS = 10
TARGET = 1.0


# ---------------------------------------------------------------------------
# The oracle
# ---------------------------------------------------------------------------


def oracle_coverage(panel):
    """Fraction of bins the panel covers. A union, computed directly."""
    covered = set()
    for primer in panel:
        covered |= BINS[primer]
    return len(covered) / TOTAL_BINS


def qualifying_subsets(candidates, target=TARGET):
    """Every subset that reaches the target. All of them, by enumeration."""
    found = []
    for size in range(len(candidates) + 1):
        for subset in itertools.combinations(sorted(candidates), size):
            if oracle_coverage(subset) >= target:
                found.append(subset)
    return found


def assert_matches_enumerated_optimum(result, qualifying):
    best_size = min(map(len, qualifying))
    assert len(result.primers) == best_size
    assert frozenset(result.primers) in {frozenset(p) for p in qualifying}


# ---------------------------------------------------------------------------
# The fixture is what the docstring says it is
# ---------------------------------------------------------------------------


def test_the_enumerated_optimum_is_two():
    qualifying = qualifying_subsets(BINS)

    assert min(map(len, qualifying)) == 2
    assert frozenset({P4, P3}) in {frozenset(p) for p in qualifying}


def test_no_single_deletion_from_the_three_primer_panel_qualifies():
    """This is what makes deletion-only insufficient rather than merely slow."""
    incumbent = (P1, P2, P3)

    assert oracle_coverage(incumbent) >= TARGET
    for primer in incumbent:
        smaller = tuple(p for p in incumbent if p != primer)
        assert oracle_coverage(smaller) < TARGET, primer


def test_a_swap_then_a_deletion_reaches_the_optimum():
    """The two-step move the search has to be able to make.

    Swap P1 for P4, which does not shrink the panel and does not lose
    coverage, then delete P2, which is now redundant.
    """
    after_swap = (P4, P2, P3)
    assert oracle_coverage(after_swap) >= TARGET

    after_delete = (P4, P3)
    assert oracle_coverage(after_delete) >= TARGET
    assert len(after_delete) == 2


# ---------------------------------------------------------------------------
# What the production reduction does
# ---------------------------------------------------------------------------


def make_optimizer(candidates):
    """An optimizer whose metrics come from the oracle's coverage.

    A fake, deliberately: this file is about which subsets the SEARCH reaches,
    and a real evaluator would make the expected answer depend on binding
    geometry rather than on a table a reader can check.
    """

    def compute_metrics(panel):
        coverage = oracle_coverage(panel)
        return replace(
            PrimerSetMetrics.empty(),
            fg_coverage=coverage,
            effective_fg_coverage=coverage,
            selectivity_density=50.0,
            total_bg_sites=1,
        )

    return SimpleNamespace(
        compute_metrics=compute_metrics,
        config=OptimizerConfig(target_set_size=len(candidates)),
        conditions=SimpleNamespace(fingerprint=lambda: "test", temp=30.0),
        bg_prefixes=["bg"],
        bg_seq_lengths=[1000],
        fg_prefixes=["fg"],
        fg_seq_lengths=[1000],
    )


def result_for(panel):
    metrics = make_optimizer(BINS).compute_metrics(panel)
    return OptimizationResult(
        tuple(panel),
        metrics.normalized_score(),
        OptimizationStatus.SUCCESS,
        metrics,
        0,
        "fixture",
    )


def test_deletion_alone_stops_above_the_enumerated_optimum():
    """Measured, and the gap this task exists to close.

    `reduce_result` removes one primer at a time and keeps the panel when no
    single removal qualifies. On this fixture that is a local optimum at 3
    while the enumerated answer is 2. Nothing is wrong with the deletion loop;
    it is doing what it says. What it is not is a search for the smallest pool.
    """
    from neoswga.core.optimization_service import reduce_result

    optimizer = make_optimizer(BINS)
    diagnostics = {}
    reduced = reduce_result(
        result_for((P1, P2, P3)),
        optimizer,
        TARGET,
        diagnostics=diagnostics,
    )

    qualifying = qualifying_subsets(BINS)
    assert min(map(len, qualifying)) == 2
    assert len(reduced.primers) == 3, "deletion alone cannot reach the two-primer answer"
    assert diagnostics["stop_reason"] == "no_qualifying_deletion"


def test_deletion_does_reach_the_optimum_when_one_exists_below_it():
    """The deletion loop is not broken, which is why the gap is a design one.

    Given a panel with a genuinely redundant member, it removes it.
    """
    from neoswga.core.optimization_service import reduce_result

    optimizer = make_optimizer(BINS)
    reduced = reduce_result(result_for((P4, P2, P3)), optimizer, TARGET)

    assert_matches_enumerated_optimum(reduced, qualifying_subsets(BINS))


def test_a_wider_frontier_holds_a_smaller_qualifying_pool():
    """First feasibility is not smallest.

    The initial frontier {P1,P2,P3} already meets the target, so a search that
    stops there reports 3. The candidate that makes 2 possible is P4, which is
    only in the wider frontier. Continuing past feasibility is the only way to
    find it.
    """
    narrow = qualifying_subsets({P1, P2, P3})
    wide = qualifying_subsets(BINS)

    assert min(map(len, narrow)) == 3
    assert min(map(len, wide)) == 2


# ---------------------------------------------------------------------------
# The beam can make the two-step move
# ---------------------------------------------------------------------------


def test_the_beam_reaches_the_enumerated_optimum_at_the_right_size():
    """`panel_beam` builds panels of a requested size from the whole pool.

    It is the piece that can reach {P4, P3}, because it does not start from an
    incumbent and delete: it grows a panel of the size it is asked for. Asked
    for 2, it must find the one qualifying pair.
    """
    from neoswga.core.panel_beam import beam_search
    from neoswga.core.panel_acceptance import PoolConstraints
    from neoswga.core.panel_refinement import objective_for_optimizer

    optimizer = make_optimizer(BINS)
    objective = objective_for_optimizer(optimizer, PoolConstraints())

    outcome = beam_search(
        objective=objective,
        candidates=sorted(BINS),
        size=2,
        beam_width=8,
    )

    assert oracle_coverage(outcome.primers) >= TARGET, outcome.primers
    assert frozenset(outcome.primers) == frozenset({P4, P3})


def test_the_beam_at_the_incumbent_size_does_not_do_worse():
    """A larger request must not return less coverage than a smaller one."""
    from neoswga.core.panel_beam import beam_search
    from neoswga.core.panel_acceptance import PoolConstraints
    from neoswga.core.panel_refinement import objective_for_optimizer

    optimizer = make_optimizer(BINS)
    objective = objective_for_optimizer(optimizer, PoolConstraints())

    at_two = beam_search(objective=objective, candidates=sorted(BINS), size=2, beam_width=8)
    at_three = beam_search(objective=objective, candidates=sorted(BINS), size=3, beam_width=8)

    assert oracle_coverage(at_three.primers) >= oracle_coverage(at_two.primers)


# ---------------------------------------------------------------------------
# Claims about minimality
# ---------------------------------------------------------------------------


def test_exhausting_candidates_is_not_exhausting_subsets():
    """The claim the plan forbids, stated as arithmetic.

    A search that has examined every candidate has examined `n` things. The
    subsets number `2**n`. "No pool found" after the first is not a statement
    about the second, and only an exhaustive run over a toy case like this one
    can produce a minimum certificate.
    """
    candidates = sorted(BINS)
    subsets = 2 ** len(candidates)

    assert len(candidates) == 4
    assert subsets == 16
    assert len(qualifying_subsets(BINS)) < subsets
