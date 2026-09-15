"""Candidates survive ranking; only declared hard gates remove them.

Task 3 of the condition-aware pool design plan. The measured Wolbachia funnel
shows the problem: 20,670 candidates survive the Gini filter and the ranking cap
keeps 2,000, so 90.3% are deleted by a QUALITY ORDER rather than by a stated
requirement. An optimizer cannot later select them, and the saved CSV cannot be
described as the candidate universe because it is a truncation of one.

The distinction this module draws:

hard gate
    A declared requirement. Sequence composition, self-dimer limits, exclusion
    genomes and the configured frequency and background limits. Failing one is a
    reason to be ineligible, and the reason is recorded.
measurement
    Gini, abundance, rank, occupancy. Useful for ordering a search; not grounds
    for permanent deletion. A primer that only covers a difficult region ranks
    poorly and may still be the one a pool needs.

Assessments are keyed by condition and QC-policy version, because eligibility is
a property of a candidate UNDER A REACTION. The same oligo can fail a Tm window
under one chemistry and pass under another, and both verdicts are worth keeping.
"""

import pytest

from neoswga.core.candidate_inventory import CandidateInventory

A, C, G = "AAAAAAAAAAAA", "CCCCCCCCCCCC", "GGGGGGGGGGGG"
COND = "tm-2026-09-14:abc123"
OTHER = "tm-2026-09-14:def456"


@pytest.fixture
def inventory(tmp_path):
    with CandidateInventory(tmp_path / "inv.sqlite") as inv:
        yield inv


def test_a_poorly_ranked_candidate_stays_eligible(inventory):
    """Rank is a measurement. It orders a search; it does not delete."""
    inventory.record_candidate(A, {"gini": 0.98, "rank": 19_999, "fg_count": 2})
    inventory.record_assessment(A, COND, passed=True, reasons=[], metrics={})

    assert list(inventory.iter_eligible(COND, [12])) == [A]


def test_a_hard_gate_failure_is_not_eligible(inventory):
    inventory.record_candidate(C, {"gini": 0.1})
    inventory.record_assessment(C, COND, passed=False, reasons=["self_dimer"], metrics={})

    assert list(inventory.iter_eligible(COND, [12])) == []


def test_the_recorded_reason_survives(inventory):
    inventory.record_candidate(C, {})
    inventory.record_assessment(
        C, COND, passed=False, reasons=["self_dimer", "exclusion_genome"], metrics={}
    )

    assert inventory.reasons(C, COND) == ["self_dimer", "exclusion_genome"]


def test_eligibility_is_a_property_of_the_candidate_under_a_reaction(inventory):
    """The same oligo, two chemistries, two verdicts, both kept."""
    inventory.record_candidate(A, {})
    inventory.record_assessment(A, COND, passed=False, reasons=["tm_window"], metrics={})
    inventory.record_assessment(A, OTHER, passed=True, reasons=[], metrics={})

    assert list(inventory.iter_eligible(COND, [12])) == []
    assert list(inventory.iter_eligible(OTHER, [12])) == [A]


def test_a_different_policy_version_is_a_separate_verdict(inventory):
    inventory.record_candidate(A, {})
    inventory.record_assessment(
        A, COND, passed=False, reasons=["tm"], metrics={}, policy_version="v1"
    )
    inventory.record_assessment(A, COND, passed=True, reasons=[], metrics={}, policy_version="v2")

    assert list(inventory.iter_eligible(COND, [12], policy_version="v1")) == []
    assert list(inventory.iter_eligible(COND, [12], policy_version="v2")) == [A]


def test_importing_the_same_candidate_twice_is_harmless(inventory):
    inventory.record_candidate(A, {"fg_count": 2})
    inventory.record_candidate(A, {"fg_count": 5})
    inventory.record_assessment(A, COND, passed=True, reasons=[], metrics={})
    inventory.record_assessment(A, COND, passed=True, reasons=[], metrics={})

    assert list(inventory.iter_eligible(COND, [12])) == [A]
    assert inventory.metrics(A)["fg_count"] == 5, "the later measurement should win"


def test_iteration_is_deterministic(inventory):
    for seq in (G, A, C):
        inventory.record_candidate(seq, {})
        inventory.record_assessment(seq, COND, passed=True, reasons=[], metrics={})

    first = list(inventory.iter_eligible(COND, [12]))
    assert first == sorted(first), "an unordered scan makes a run irreproducible"
    assert first == list(inventory.iter_eligible(COND, [12]))


def test_only_the_requested_lengths_come_back(inventory):
    inventory.record_candidate(A, {})
    inventory.record_candidate("ACGTACGTAC", {})
    for seq in (A, "ACGTACGTAC"):
        inventory.record_assessment(seq, COND, passed=True, reasons=[], metrics={})

    assert list(inventory.iter_eligible(COND, [10])) == ["ACGTACGTAC"]
    assert sorted(inventory.iter_eligible(COND, [10, 12])) == sorted([A, "ACGTACGTAC"])


def test_a_committed_batch_survives_reopening(tmp_path):
    """Resume: an interrupted design should not rescan what it already counted."""
    path = tmp_path / "inv.sqlite"
    with CandidateInventory(path) as inv:
        inv.record_candidate(A, {"fg_count": 3})
        inv.record_assessment(A, COND, passed=True, reasons=[], metrics={})

    with CandidateInventory(path) as reopened:
        assert list(reopened.iter_eligible(COND, [12])) == [A]
        assert reopened.metrics(A)["fg_count"] == 3


def test_an_unassessed_candidate_is_counted_but_not_eligible(inventory):
    """Counted and assessed are different numbers, and both are reported."""
    inventory.record_candidate(A, {})

    assert inventory.counts()["candidates"] == 1
    assert inventory.counts()["assessed"] == 0
    assert list(inventory.iter_eligible(COND, [12])) == []


def test_the_filter_step_records_what_it_enumerated(tmp_path, monkeypatch):
    """The inventory has to be written, not merely available.

    `record_stage2_inventory` is what stage 2 calls. Candidates cut by the
    ranking cap stay ELIGIBLE: `max_primer` orders a search, it is not a
    declared requirement, and the measured Wolbachia design lost 90.3% of its
    Gini survivors to it.
    """
    import pandas as pd

    from neoswga.core.candidate_inventory import record_stage2_inventory

    cleared_hard_gates = pd.DataFrame(
        {"primer": [A, C, G], "fg_count": [9, 5, 1], "bg_count": [0, 1, 2]}
    )
    after_gini = cleared_hard_gates[cleared_hard_gates["primer"] != G]
    shortlisted = after_gini[after_gini["primer"] == A]

    path = record_stage2_inventory(
        tmp_path,
        condition_id=COND,
        cleared_hard_gates=cleared_hard_gates,
        after_gini=after_gini,
        shortlisted=shortlisted,
    )

    with CandidateInventory(path) as inv:
        counts = inv.counts()
        assert counts["candidates"] == 3
        assert counts["hard_qc_passed"] == 3, "the ranking cap must not make a candidate ineligible"
        eligible = list(inv.iter_eligible(COND, [12]))
        assert eligible == sorted([A, C, G])
        # The measurements that ordered the search are kept, and say so.
        assert inv.metrics(A)["shortlisted"] is True
        assert inv.metrics(G)["shortlisted"] is False
        assert inv.metrics(G)["passed_gini"] is False
        assert inv.metrics(A)["passed_gini"] is True
