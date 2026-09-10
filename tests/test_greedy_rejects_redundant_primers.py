"""Near-duplicate primers reached the delivered panel.

Audit finding A6. Delivered E. coli set 0 -- 160 primers drawn from a
449-candidate pool -- holds 9 pairs at Hamming distance 1 or less, and 5 primers
share the 3' hexamer GCGAAA, 4 share GCCAGA, 4 share GGATAA. The pool itself
holds 65 pairs at Hamming distance 1 or less. Primers differing by one base bind
largely the same sites, so the second buys little coverage while adding
synthesis cost and dimer surface, and shared 3' ends correlate mispriming
behaviour across the panel.

The greedy picked the candidate with the largest number of NEW covered bins,
which is not the same as picking a candidate that is not redundant: a primer
with fifty sites, forty-eight of them already covered, can still beat a primer
with five sites all of which are new.

The criterion is bin containment, not sequence similarity. The covered bin sets
are already in the graph the optimizer holds, two primers can differ by one base
and bind different sites, and two unrelated sequences can bind the same places.
Hamming distance is the symptom the audit counted, not the thing being tested.

EVERY SEQUENCE IN THIS FIXTURE CONTAINS ONLY A AND G. That is deliberate and
load-bearing. `optimize_greedy` also carries a dimer rejection guard, added by
`docs/superpowers/plans/2026-09-06-optimizer-cost-and-dimer-criterion.md`, which
builds a DimerMatrix over the candidates and skips any primer that pairs with
the selected set. A sequence with no C and no T cannot form a Watson-Crick pair
with another such sequence at any `max_dimer_bp`, so this fixture isolates the
redundancy rule from that guard rather than colliding with it.

Audit finding D2 tested and refuted the alternative explanation, that the
candidate ORDER lets the greedy take runs of near-duplicates. Near-duplicates
sit 12% to 35% closer together than chance on the median gap, but the
near-adjacency excess is tiny and on M. tuberculosis it reverses. The absence
of a rejection test in the selection rule is the cause.
"""

import numpy as np
import pytest

GENOME = 1_000_000
REACH = 3000
BIN = 750

# A and B differ by one base, which is the pattern the audit counted, but the
# rule under test does not look at that.
#
# A covers 25 well-separated sites. B repeats twenty of them, offset by 10 bp so
# they land in the same bins, and adds two sites of its own. C has one site,
# entirely its own. Measured against the real BipartiteGraph at these settings:
#
#   A: 125 regions
#   B: 110 regions, 10 new after A, redundancy 0.909
#   C:   5 regions,  5 new after A, redundancy 0.000
#
# So the plain greedy takes B second (10 new beats 5) and a redundancy rule at
# 0.9 takes C.
A = "AAGAGGAAGAGG"
B = "AAGAGGAAGAGA"  # Hamming distance 1 from A
C = "GGAAGAGGAAGA"
A_SITES = [20_000 + i * 30_000 for i in range(25)]
SITES = {
    A: A_SITES,
    B: [p + 10 for p in A_SITES[:20]] + [850_000, 880_000],
    C: [950_000],
}


class _Cache:
    def get_positions(self, prefix, primer, strand="both"):
        if strand == "reverse":
            return np.array([], dtype=np.int64)
        return np.asarray(SITES.get(primer, []), dtype=np.int64)


@pytest.fixture
def optimizer():
    from neoswga.core.dominating_set_optimizer import DominatingSetOptimizer

    return DominatingSetOptimizer(
        cache=_Cache(),
        fg_prefixes=["fg"],
        fg_seq_lengths=[GENOME],
        bin_size=BIN,
        extension_reach=REACH,
    )


def test_the_fixture_cannot_trip_the_dimer_guard():
    """Guard the guard, part one.

    `optimize_greedy` also skips a candidate that dimerises with the selected
    set. If this fixture tripped that guard, every assertion below would pass
    for the wrong reason. A sequence with no C and no T cannot form a
    Watson-Crick pair with another such sequence at any max_dimer_bp.
    """
    for sequence in (A, B, C):
        assert set(sequence) <= {"A", "G"}, sequence


def test_a_and_b_are_the_near_duplicate_pair_the_audit_counted():
    """The rule does not measure this, but the fixture should still be the
    shape A6 describes."""
    assert sum(x != y for x, y in zip(A, B)) == 1


def test_the_default_threshold_is_stated():
    from neoswga.core.dominating_set_optimizer import DEFAULT_REDUNDANCY_THRESHOLD

    # Disabled by default on measurement, 2026-09-10: at 0.9 the rule fired
    # 328,846 times across a ten-case sweep of the three GC-tier pools and
    # changed neither coverage nor either redundancy measure it exists to
    # reduce. The mechanism is kept and every test below drives it with an
    # explicit threshold.
    assert DEFAULT_REDUNDANCY_THRESHOLD == 1.0


def test_the_fixture_isolates_the_fault(optimizer):
    """Guard the guard. B must genuinely beat C on absolute new coverage, or
    this file proves nothing about the redundancy rule."""
    result = optimizer.optimize_greedy(
        candidates=list(SITES), max_primers=1, verbose=False, redundancy_threshold=1.0
    )
    graph = result["graph"]
    covered = set(graph.primer_to_regions[A])
    new_b = len(graph.primer_to_regions[B] - covered)
    new_c = len(graph.primer_to_regions[C] - covered)
    assert result["ordered_primers"] == [A]
    assert new_b > new_c
    assert 1 - new_b / len(graph.primer_to_regions[B]) > 0.9


def test_a_redundant_candidate_is_skipped_for_an_independent_one(optimizer):
    """The regression."""
    result = optimizer.optimize_greedy(
        candidates=list(SITES), max_primers=2, verbose=False, redundancy_threshold=0.9
    )
    assert result["ordered_primers"] == [A, C]


def test_a_threshold_of_one_reproduces_the_old_behaviour(optimizer):
    """An escape hatch that is also the before-picture."""
    result = optimizer.optimize_greedy(
        candidates=list(SITES), max_primers=2, verbose=False, redundancy_threshold=1.0
    )
    assert result["ordered_primers"] == [A, B]


def test_a_redundant_candidate_is_still_taken_when_nothing_else_is_left(optimizer):
    """The fallback. Without it a long panel truncates and loses coverage that
    the plain greedy would have delivered."""
    result = optimizer.optimize_greedy(
        candidates=[A, B], max_primers=2, verbose=False, redundancy_threshold=0.9
    )
    assert result["ordered_primers"] == [A, B]


def test_the_rule_does_not_lower_coverage_at_the_full_set_size(optimizer):
    """Taking all three candidates must cover the same bins either way."""
    strict = optimizer.optimize_greedy(
        candidates=list(SITES), max_primers=3, verbose=False, redundancy_threshold=0.9
    )
    loose = optimizer.optimize_greedy(
        candidates=list(SITES), max_primers=3, verbose=False, redundancy_threshold=1.0
    )
    assert set(strict["ordered_primers"]) == set(loose["ordered_primers"])
    assert strict["coverage"] == pytest.approx(loose["coverage"])


def test_the_first_pick_is_unaffected(optimizer):
    """Nothing is covered yet, so no candidate can be redundant."""
    result = optimizer.optimize_greedy(
        candidates=list(SITES), max_primers=1, verbose=False, redundancy_threshold=0.9
    )
    assert result["ordered_primers"] == [A]


def test_fixed_primers_count_toward_the_covered_union(optimizer):
    """A fixed primer's coverage is pre-filled, so a candidate redundant
    against it must be skipped."""
    result = optimizer.optimize_greedy(
        candidates=[B, C],
        max_primers=1,
        fixed_primers=[A],
        verbose=False,
        redundancy_threshold=0.9,
    )
    assert result["new_primers"] == [C]


# ---------------------------------------------------------------------------
# The two hazards in sharing the scan with the dimer guard
# ---------------------------------------------------------------------------


def test_the_redundancy_fallback_cannot_readmit_a_dimerising_primer(optimizer, monkeypatch):
    """The hazard. If the redundancy fallback were tracked BEFORE the dimer
    guard, a primer the dimer screen rejected could come back as the fallback
    and silently undo that constraint.

    B is the only candidate the fallback could reach here, so forcing the dimer
    guard to reject it must leave the fallback empty rather than selecting it.
    """
    monkeypatch.setattr(
        type(optimizer),
        "_would_dimerise",
        lambda self, candidate, selected, matrix: candidate == B and bool(selected),
    )
    result = optimizer.optimize_greedy(
        candidates=[A, B], max_primers=2, verbose=False, redundancy_threshold=0.9
    )
    assert B not in result["ordered_primers"]


def test_a_redundancy_skip_does_not_look_like_a_dimer_skip(optimizer):
    """`skipped_for_dimer` is the optimizer plan's signal that the DIMER
    constraint stalled the loop, and it drives that plan's relaxation. A
    redundancy skip must not set it, or the dimer constraint would be dropped
    for a reason that has nothing to do with dimers.

    Observed through the outcome: B is skipped for redundancy here, and if that
    had set the dimer flag the relaxation would have fired and logged a warning
    naming max_dimer_bp.
    """
    import logging

    with pytest.MonkeyPatch.context():
        records = []
        handler = logging.Handler()
        handler.emit = records.append
        logger = logging.getLogger("neoswga.core.dominating_set_optimizer")
        logger.addHandler(handler)
        try:
            optimizer.optimize_greedy(
                candidates=list(SITES), max_primers=2, verbose=True, redundancy_threshold=0.9
            )
        finally:
            logger.removeHandler(handler)

    assert not any("max_dimer_bp" in record.getMessage() for record in records)


def test_the_two_guards_cannot_deadlock(optimizer, monkeypatch):
    """The roadmap asks for this explicitly.

    A candidate rejected by one guard and admitted by the other must still be
    reachable, and the two guards together must never deliver a smaller set
    than either alone. Here the dimer guard rejects C and the redundancy guard
    rejects B, which between them cover every candidate but A. The run must
    still return as many primers as the dimer guard alone would.
    """
    monkeypatch.setattr(
        type(optimizer),
        "_would_dimerise",
        lambda self, candidate, selected, matrix: candidate == C and bool(selected),
    )
    both = optimizer.optimize_greedy(
        candidates=list(SITES), max_primers=3, verbose=False, redundancy_threshold=0.9
    )
    dimer_only = optimizer.optimize_greedy(
        candidates=list(SITES), max_primers=3, verbose=False, redundancy_threshold=1.0
    )
    assert len(both["ordered_primers"]) == len(dimer_only["ordered_primers"])
    assert both["coverage"] == pytest.approx(dimer_only["coverage"])
