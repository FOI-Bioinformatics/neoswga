"""Every candidate that cleared hard QC stays reachable by a search.

Task 3 of the 2026-09-21 valid-design plan, the inventory half. The plan's
global constraint is that ranking and a bounded frontier may postpone when a
candidate is examined but must never delete it from the universe.

That constraint has been broken twice in ways that were invisible from inside.
`InventoryCandidateSource.advance` returned False from the day it was written,
so the 18,670 candidates beyond the Wolbachia shortlist could not affect any
panel. And the retired `legacy` retention mode gave a background index to the
shortlist alone, so 489,836 candidates that cleared hard QC scored against an
empty background and read as perfectly specific.

The expected sequences here are enumerated from the input FASTA by this file,
not read back from the inventory under test. An expectation derived from the
thing it checks agrees with it by construction, which is how a source that
silently dropped a third of its candidates would still have passed.
"""

import itertools

import pytest

from neoswga.core.candidate_inventory import CandidateInventory
from neoswga.core.candidate_provider import CandidateProvider
from neoswga.core.candidate_source import InventoryCandidateSource

CONDITION = "tm-test:completeness"
OTHER_CONDITION = "tm-test:a-warmer-reaction"

# A short sequence with no repeats, so every k-mer of a given length is
# distinct and the expected set is exactly the sliding window.
GENOME = "ACGTTGCAAGGCTTACCGATGCATGGCTAACGTCAGTCCAAGTTGCACTGA"


def enumerated_kmers(sequence, length):
    """Every distinct k-mer of one length, counted here rather than read back.

    Deliberately a sliding window over the string rather than a call into the
    package: an expectation computed by the code under test is not a check on
    it. This is the same reason the plan asks for an independent oracle for
    coverage.
    """
    windows = [sequence[i : i + length] for i in range(len(sequence) - length + 1)]
    return sorted({window for window in windows if len(window) == length})


def assert_inventory_preserved(qc_eligible, visited, rejection_reasons):
    assert set(visited) == set(qc_eligible)
    assert not set(qc_eligible).intersection(rejection_reasons)


def stock(path, sequences, condition=CONDITION, rejected=()):
    """An inventory holding `sequences` as eligible and `rejected` as not."""
    inventory = CandidateInventory(path)
    for rank, sequence in enumerate(sorted(sequences) + sorted(rejected)):
        inventory.record_candidate(sequence, {"fg_count": 1, "bg_count": 0}, search_rank=rank)
    for sequence in sequences:
        inventory.record_assessment(sequence, condition, True, [], {})
    for sequence in rejected:
        inventory.record_assessment(sequence, condition, False, ["gini_above_max"], {})
    inventory.commit()
    return inventory


def exhaust(source):
    """Every candidate a search would ever see, frontier by frontier."""
    visited = list(source.initial())
    while source.advance(keep=visited):
        for sequence in source.frontier():
            if sequence not in visited:
                visited.append(sequence)
    return visited


@pytest.fixture
def small_inventory(tmp_path):
    def build(lengths=(10,), condition=CONDITION, rejected=()):
        expected = sorted(
            itertools.chain.from_iterable(enumerated_kmers(GENOME, k) for k in lengths)
        )
        path = tmp_path / f"inventory_{'_'.join(map(str, lengths))}_{condition[-6:]}.sqlite"
        inventory = stock(path, expected, condition=condition, rejected=rejected)
        return expected, inventory

    return build


# ---------------------------------------------------------------------------
# Nothing eligible is unreachable
# ---------------------------------------------------------------------------


def test_a_bounded_frontier_postpones_rather_than_deletes(small_inventory):
    """More eligible candidates than the frontier holds. All must be reachable.

    This is the shape `advance` returned False for. A frontier smaller than the
    universe is normal and wanted; a frontier that is also a ceiling is not.
    """
    expected, inventory = small_inventory()
    provider = CandidateProvider(inventory, CONDITION, [10])
    source = InventoryCandidateSource(provider, frontier=5)

    assert len(source.initial()) == 5, "the frontier must actually be bounded"
    assert source.universe_size() == len(expected)

    visited = exhaust(source)
    assert_inventory_preserved(expected, visited, rejection_reasons=[])
    inventory.close()


def test_the_universe_is_reported_whole_even_while_the_frontier_is_not(small_inventory):
    expected, inventory = small_inventory()
    source = InventoryCandidateSource(CandidateProvider(inventory, CONDITION, [10]), frontier=3)

    assert sorted(source.universe()) == expected
    assert len(source.frontier()) == 0, "nothing is examined before initial()"
    inventory.close()


@pytest.mark.parametrize("lengths", [(10,), (10, 11), (9, 10, 11, 12)])
def test_every_configured_length_survives(small_inventory, lengths):
    """A multi-length design must not lose a length to the traversal."""
    expected, inventory = small_inventory(lengths=lengths)
    source = InventoryCandidateSource(
        CandidateProvider(inventory, CONDITION, list(lengths)), frontier=4
    )

    visited = exhaust(source)
    assert_inventory_preserved(expected, visited, rejection_reasons=[])
    assert {len(sequence) for sequence in visited} == set(lengths)
    inventory.close()


def test_a_changed_chemistry_reaches_its_own_eligible_set(small_inventory):
    """Eligibility is per reaction; completeness must hold under each one."""
    expected, inventory = small_inventory(condition=OTHER_CONDITION)
    source = InventoryCandidateSource(
        CandidateProvider(inventory, OTHER_CONDITION, [10]), frontier=6
    )

    visited = exhaust(source)
    assert_inventory_preserved(expected, visited, rejection_reasons=[])
    inventory.close()


# ---------------------------------------------------------------------------
# What is excluded is excluded for a recorded reason
# ---------------------------------------------------------------------------


def test_an_ineligible_candidate_is_absent_and_carries_its_reason(small_inventory):
    """Removal from the universe has to be traceable to a named hard gate.

    A candidate that is simply missing is indistinguishable from one the
    traversal dropped, which is the failure the other tests in this file are
    about. The reason code is what tells the two apart.
    """
    rejected = ["TTTTTTTTTT", "AAAAAAAAAA"]
    expected, inventory = small_inventory(rejected=rejected)
    source = InventoryCandidateSource(CandidateProvider(inventory, CONDITION, [10]), frontier=5)

    visited = exhaust(source)

    assert_inventory_preserved(expected, visited, rejection_reasons=rejected)
    for sequence in rejected:
        assert sequence not in visited
        assert inventory.reasons(sequence, CONDITION) == ["gini_above_max"]
    inventory.close()


def test_a_candidate_nobody_assessed_is_not_silently_eligible(small_inventory):
    """Never assessed and assessed-and-passed must not read the same.

    `record_candidate` runs before any condition-specific admission, so the
    inventory knows about candidates no chemistry has judged. Treating those as
    eligible would admit primers to a design under a reaction they were never
    checked against.
    """
    expected, inventory = small_inventory()
    inventory.record_candidate("GGGGGGGGGG", {"fg_count": 1, "bg_count": 0}, search_rank=999)
    inventory.commit()

    source = InventoryCandidateSource(CandidateProvider(inventory, CONDITION, [10]), frontier=5)
    visited = exhaust(source)

    assert "GGGGGGGGGG" not in visited
    assert_inventory_preserved(expected, visited, rejection_reasons=[])
    assert inventory.reasons("GGGGGGGGGG", CONDITION) == []
    inventory.close()


# ---------------------------------------------------------------------------
# The traversal order is stable
# ---------------------------------------------------------------------------


def test_a_refill_extends_the_window_rather_than_reshuffling_it(small_inventory):
    """The optimizers are order-sensitive, so the order must not move.

    A traversal that reshuffled on refill would make an otherwise reproducible
    run depend on how many refills it happened to take.
    """
    _expected, inventory = small_inventory()
    source = InventoryCandidateSource(CandidateProvider(inventory, CONDITION, [10]), frontier=4)

    first = source.initial()
    source.advance(keep=first)
    widened = source.frontier()

    assert widened[: len(first)] == first
    assert len(widened) > len(first)
    inventory.close()


def test_advance_reports_false_only_when_the_universe_is_covered(small_inventory):
    expected, inventory = small_inventory()
    source = InventoryCandidateSource(CandidateProvider(inventory, CONDITION, [10]), frontier=5)
    source.initial()

    while source.advance(keep=()):
        pass

    assert len(source.frontier()) == len(expected)
    assert source.advance(keep=()) is False
    inventory.close()
