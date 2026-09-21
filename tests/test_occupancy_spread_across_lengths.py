"""A mixed-length design is told what it is asking of one reaction.

Occupancy rises steeply with length at a fixed temperature. Measured on 3,000
random k-mers per length at phi29 30 C, median occupancy runs 0.017 at k=7 to
0.9997 at k=12. So a panel spanning those lengths holds oligos bound under 2%
of the time beside oligos bound essentially always, and a panel size counts
them equally.

The pool-wide mean `filter` already reports cannot show this. On the bundled
plasmid pool it is 2.405 while the per-length values run 3.010 down to 1.015 --
a single number describing no length in the pool.

Reported, never enforced, which is the resolution Known Issue 17 reached after
measuring that an occupancy GATE made the delivered panel worse on both axes.
These tests therefore assert that nothing is removed, as much as that something
is said.
"""

import logging

import pytest

from neoswga.core.length_occupancy import (
    WEAK_OCCUPANCY,
    LengthOccupancy,
    log_occupancy_spread,
    occupancy_by_length,
)
from neoswga.core.reaction_conditions import ReactionConditions

LOGGER = "neoswga.core.length_occupancy"

SHORT = ["ACCGCAT", "GGCAACA", "TCAGCGA"]  # 7-mers
LONG = ["CCCATTGACG", "ACGTCAATGG"]  # 10-mers


@pytest.fixture
def phi29():
    return ReactionConditions(temp=30.0, polymerase="phi29")


# ---------------------------------------------------------------------------
# The measurement
# ---------------------------------------------------------------------------


def test_lengths_are_reported_separately_shortest_first(phi29):
    records = occupancy_by_length(SHORT + LONG, phi29)

    assert [record.length for record in records] == [7, 10]
    assert [record.n for record in records] == [3, 2]


def test_occupancy_rises_with_length(phi29):
    """The finding the module exists for, on real sequences."""
    short, long = occupancy_by_length(SHORT + LONG, phi29)

    assert short.median < long.median, (short.median, long.median)
    assert long.median > 0.9 and short.median < 0.5


def test_a_single_length_pool_still_measures(phi29):
    (record,) = occupancy_by_length(SHORT, phi29)

    assert record.length == 7
    assert record.lowest <= record.median <= record.highest


def test_no_temperature_means_no_measurement_rather_than_zero(phi29):
    """Absence and "never bound" must not share a representation."""
    assert occupancy_by_length(SHORT, None) == ()
    assert occupancy_by_length([], phi29) == ()


def test_a_primer_that_cannot_be_evaluated_is_counted_not_dropped():
    """Skipping it silently leaves a median over an unknown subset.

    That is the failure `occupancy.discrimination_profile` records at its own
    except clause, and the count is what keeps a reader able to see it.

    The unevaluable primer is made so by the conditions rather than by its
    bases, because an ambiguous sequence does NOT raise here:
    `ReactionConditions.calculate_effective_tm` returns a number for
    "NNNNNNN" computed from default penalty stacks, while
    `thermodynamics.calculate_tm_batch` refuses the same input with
    `InvalidSequenceError`. That disagreement is real and is not this file's
    subject.
    """

    class OneBadPrimer:
        temp = 30.0

        def calculate_effective_tm(self, sequence):
            if sequence == SHORT[0]:
                raise RuntimeError("no Tm for this one")
            return 25.0

    records = occupancy_by_length(SHORT, OneBadPrimer())

    assert len(records) == 1
    record = records[0]
    assert record.n == len(SHORT), "the primer must still be counted"
    assert record.unmeasured == 1
    assert record.median is not None, "the two that measured still have a median"


# ---------------------------------------------------------------------------
# What gets said
# ---------------------------------------------------------------------------


def test_a_single_length_design_gains_no_output(phi29, caplog):
    """Every design in this repository before 2026-09-21 is single-length."""
    with caplog.at_level(logging.INFO, logger=LOGGER):
        log_occupancy_spread(SHORT, phi29)

    assert caplog.text == ""


def test_a_mixed_design_is_told_the_spread(phi29, caplog):
    with caplog.at_level(logging.INFO, logger=LOGGER):
        log_occupancy_spread(SHORT + LONG, phi29, label="panel")

    assert "k=7" in caplog.text and "k=10" in caplog.text
    assert "panel" in caplog.text


def test_a_weak_length_beside_a_saturated_one_warns(phi29, caplog):
    """The actionable case, and the message must name the reaction.

    An occupancy gate was measured and makes panels worse, so the advice has
    to point at the temperature rather than at the candidate filter.
    """
    with caplog.at_level(logging.INFO, logger=LOGGER):
        log_occupancy_spread(["ATATATA", "ATATATAT"] + LONG, phi29, label="panel")

    warnings = [r for r in caplog.records if r.levelno >= logging.WARNING]
    assert warnings, caplog.text
    message = warnings[0].getMessage()
    assert "temperature" in message
    assert "Nothing is removed" in message


def test_the_report_removes_nothing(phi29):
    """Stated as a test because the whole design decision rests on it."""
    pool = SHORT + LONG
    before = list(pool)

    log_occupancy_spread(pool, phi29)

    assert pool == before


def test_a_failure_in_the_diagnostic_does_not_fail_the_step(phi29, caplog):
    """It describes a pool; it must never be the reason that pool is lost."""

    class Hostile:
        temp = 30.0

        def calculate_effective_tm(self, sequence):
            raise RuntimeError("boom")

    with caplog.at_level(logging.INFO, logger=LOGGER):
        log_occupancy_spread(SHORT + LONG, Hostile())  # must not raise


def test_the_weak_threshold_is_a_description_not_a_gate():
    """No caller may use it to reject a candidate."""
    assert 0.0 < WEAK_OCCUPANCY < 1.0

    import ast
    import pathlib

    root = pathlib.Path(__file__).resolve().parent.parent / "neoswga"
    users = [
        f"{path.relative_to(root)}:{node.lineno}"
        for path in root.rglob("*.py")
        if path.name != "length_occupancy.py"
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8")))
        if isinstance(node, ast.Name) and node.id == "WEAK_OCCUPANCY"
    ]
    assert not users, (
        "WEAK_OCCUPANCY is read outside the module that reports it. It is a "
        "description of a spread, not a threshold anything may be rejected "
        f"for -- see Known Issue 17. Read at: {users}"
    )


def test_records_are_immutable():
    record = LengthOccupancy(7, 3, 0.1, 0.05, 0.2, 0)

    with pytest.raises(Exception):
        record.median = 0.9
