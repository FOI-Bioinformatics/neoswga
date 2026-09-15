"""The search budget bounds what is examined, not what exists.

Task 4 of the condition-aware pool design plan. Task 3 made every hard-QC
candidate durable and indexed; this makes them reachable.

`max_primer` seeds a working shortlist. When that shortlist cannot meet the
constraints, the search asks for more rather than concluding the design is
infeasible. The distinction matters because the two look identical from inside a
greedy that has run out of candidates: "nothing left that helps" and "nothing
left in the batch I was given" produce the same stall.

Expansion is deterministic. The optimizers this feeds are order-sensitive, so an
unordered traversal would make an otherwise reproducible run depend on which
candidates happened to come back first.
"""

import pytest

from neoswga.core.candidate_inventory import CandidateInventory
from neoswga.core.candidate_provider import CandidateProvider

COND = "tm-2026-09-14:abc123"


def _sequence(index: int) -> str:
    """Distinct 12-mers, ordered so rank and sequence order differ."""
    letters = "ACGT"
    out = []
    n = index
    for _ in range(6):
        out.append(letters[n % 4])
        n //= 4
    return "".join(out) + "AAAAAA"


@pytest.fixture
def inventory(tmp_path):
    """Twenty eligible candidates, ranked so the useful one ranks last."""
    with CandidateInventory(tmp_path / "inv.sqlite") as inv:
        for i in range(20):
            inv.record_candidate(_sequence(i), {"step2_rank": 19 - i})
            inv.record_assessment(_sequence(i), COND, passed=True, reasons=[], metrics={})
        inv.commit()
    with CandidateInventory(tmp_path / "inv.sqlite") as inv:
        yield inv


def test_the_initial_shortlist_respects_its_limit(inventory):
    provider = CandidateProvider(inventory, COND, [12])

    assert len(provider.initial(limit=5)) == 5


def test_the_shortlist_is_seeded_by_the_existing_ranking(inventory):
    """Rank still orders the search; it just no longer bounds it."""
    provider = CandidateProvider(inventory, COND, [12])

    first = provider.initial(limit=3)
    assert first == [_sequence(19), _sequence(18), _sequence(17)]


def test_expansion_reaches_a_candidate_the_shortlist_left_out(inventory):
    """The fixture's point: the useful candidate ranks below the cut."""
    provider = CandidateProvider(inventory, COND, [12])
    examined = set(provider.initial(limit=5))
    wanted = _sequence(0)  # ranked last
    assert wanted not in examined

    while True:
        batch = provider.expand(excluded=examined, limit=5)
        if not batch:
            break
        examined.update(batch)

    assert wanted in examined, "a candidate below the cut was unreachable"


def test_expansion_eventually_exhausts_the_eligible_inventory(inventory):
    provider = CandidateProvider(inventory, COND, [12])
    examined = set(provider.initial(limit=5))
    for _ in range(10):
        batch = provider.expand(excluded=examined, limit=5)
        if not batch:
            break
        examined.update(batch)

    assert len(examined) == 20
    assert provider.expand(excluded=examined, limit=5) == []


def test_expansion_never_returns_what_was_already_examined(inventory):
    provider = CandidateProvider(inventory, COND, [12])
    examined = set(provider.initial(limit=5))

    batch = provider.expand(excluded=examined, limit=5)
    assert not (set(batch) & examined)


def test_expansion_is_deterministic(inventory):
    """An order-sensitive optimizer must not depend on traversal order."""

    def run():
        provider = CandidateProvider(inventory, COND, [12])
        examined = list(provider.initial(limit=4))
        for _ in range(4):
            batch = provider.expand(excluded=set(examined), limit=4)
            if not batch:
                break
            examined.extend(batch)
        return examined

    assert run() == run()


def test_a_budget_stop_reports_what_was_not_examined(inventory):
    """Stopping early is not the same as proving nothing else would help."""
    provider = CandidateProvider(inventory, COND, [12])
    examined = set(provider.initial(limit=5))

    remaining = provider.unexamined(excluded=examined)
    assert remaining == 15
    assert provider.exhausted(excluded=examined) is False

    everything = set(provider.initial(limit=100))
    assert provider.unexamined(excluded=everything) == 0
    assert provider.exhausted(excluded=everything) is True


def test_a_primer_is_not_removed_for_conflicting_with_candidates_nobody_uses():
    """The pool-wide heterodimer count is not a property of a panel.

    The Stage-0 screen removed a primer when it conflicted with more than a
    fraction of the WHOLE pool. What matters is whether it conflicts with the
    primers actually selected, which the greedy already enforces pairwise
    against the panel it is building.

    Counting across the pool also made the verdict depend on pool size. Under
    `candidate_retention="all_qc"` the pool is an order of magnitude larger, so
    a primer that survived in a 2,000-candidate pool can be removed from a
    20,000-candidate one without anything about the primer changing.
    """
    from types import SimpleNamespace

    from neoswga.core.hybrid_thermo_screen import ThermoScreenMixin
    from neoswga.core.reaction_conditions import ReactionConditions
    from neoswga.core.thermodynamic_filter import ThermodynamicCriteria

    seen = {}

    class _Screen(ThermoScreenMixin):
        polymerase = "phi29"
        poly_config = SimpleNamespace(reaction_temp=30.0)
        conditions = ReactionConditions(temp=30.0, polymerase="phi29", mg_conc=10.0)
        max_dimer_bp = 3
        _thermo_filter_cache = None

        def _thermo_criteria(self):
            return ThermodynamicCriteria(na_conc=50.0, mg_conc=10.0)

    import neoswga.core.thermodynamic_filter as tf

    real = tf.ThermodynamicFilter.filter_candidates

    def _spy(self, candidates, **kwargs):
        seen.update(kwargs)
        return real(self, candidates, **kwargs)

    tf.ThermodynamicFilter.filter_candidates = _spy
    try:
        _Screen()._thermo_filter_candidates(["ACGTACGTACGT", "CCCCGGGGAAAA"], verbose=False)
    finally:
        tf.ThermodynamicFilter.filter_candidates = real

    assert seen.get("check_heterodimers") is False, (
        "the pool-wide heterodimer hub count is still acting as a hard gate; "
        "pairwise compatibility against the selected panel is the real constraint"
    )
