"""`--minimize-primers` must not leave a panel unchanged without saying why.

Two routes to a silent no-op were found on 2026-09-30.

A candidate list handed to `run_optimization` directly was not screened for
self-dimers, because the screen lives with the candidate source and a plain
list has none. Selection could then deliver a self-dimerising primer, the
reduction stage counts a self-dimer as a violation, every deletion inherited
it, and the stage stopped having removed nothing.

And on any path, the coverage target is compared against the objective's
coverage, which is occupancy weighted when conditions are attached, while the
log prints the unweighted figure. On the plasmid example a 12-primer panel
printed 46.7% against a 30% target and was not reduced: the figure the target
was compared against was 19.7%, already below it. Nothing in the log said so.
"""

import logging
from dataclasses import dataclass

import pytest

from neoswga.core import dimer
from neoswga.core.base_optimizer import OptimizationStatus
from neoswga.core.exceptions import NoCandidatesError
from neoswga.core.optimization_service import (
    describe_reduction,
    screen_supplied_candidates,
)
from neoswga.core.result_validation import validate_result

SELF_DIMER_AT_3 = "GGTACGTGTC"
CLEAN = "AACAGGAACA"


@dataclass
class _Config:
    max_dimer_bp: int = 3
    max_self_dimer_bp: int = 3


def test_the_fixture_primers_are_what_they_are_named():
    assert dimer.is_dimer_fast(SELF_DIMER_AT_3, SELF_DIMER_AT_3, 3)
    assert not dimer.is_dimer_fast(SELF_DIMER_AT_3, SELF_DIMER_AT_3, 4)
    assert not dimer.is_dimer_fast(CLEAN, CLEAN, 3)


def test_a_supplied_list_is_screened_at_the_configured_limit(caplog):
    with caplog.at_level(logging.WARNING):
        kept = screen_supplied_candidates([CLEAN, SELF_DIMER_AT_3], _Config())

    assert kept == [CLEAN]
    assert "1 of 2" in caplog.text and "max_self_dimer_bp=3" in caplog.text


def test_a_looser_limit_keeps_the_same_primer_and_says_nothing(caplog):
    with caplog.at_level(logging.WARNING):
        kept = screen_supplied_candidates([CLEAN, SELF_DIMER_AT_3], _Config(max_self_dimer_bp=4))

    assert kept == [CLEAN, SELF_DIMER_AT_3]
    assert caplog.text == ""


def test_a_fixed_primer_is_the_users_instruction_and_is_kept():
    kept = screen_supplied_candidates(
        [CLEAN, SELF_DIMER_AT_3], _Config(), fixed_primers=[SELF_DIMER_AT_3]
    )

    assert kept == [CLEAN, SELF_DIMER_AT_3]


def test_a_screen_that_empties_the_list_names_itself():
    with pytest.raises(NoCandidatesError, match="max_self_dimer_bp=3"):
        screen_supplied_candidates([SELF_DIMER_AT_3], _Config())


def test_an_unchanged_panel_below_the_target_says_which_coverage_was_compared():
    message = describe_reduction(
        size_before=12,
        size_after=12,
        target_coverage=0.30,
        coverage=0.197,
        coverage_metric="effective",
        stop_reason="no_qualifying_deletion",
    )

    assert "removed no primer" in message
    assert "0.197" in message and "effective" in message and "0.300" in message
    assert "already below" in message


def test_an_unchanged_panel_above_the_target_names_the_stop_reason():
    message = describe_reduction(
        size_before=4,
        size_after=4,
        target_coverage=0.20,
        coverage=0.50,
        coverage_metric="effective",
        stop_reason="no_qualifying_deletion",
    )

    assert "removed no primer" in message and "no_qualifying_deletion" in message
    assert "already below" not in message


def test_a_reduced_panel_reports_how_many_were_removed():
    message = describe_reduction(
        size_before=12,
        size_after=8,
        target_coverage=0.30,
        coverage=0.31,
        coverage_metric="effective",
        stop_reason="no_qualifying_deletion",
    )

    assert "removed 4 of 12" in message and "0.310" in message


@dataclass
class _Result:
    primers: tuple
    status: OptimizationStatus = OptimizationStatus.SUCCESS
    stage_history: tuple = ()
    metrics: object = None
    optimizer_name: str = "stub"


def _size_issue(result, target_size):
    report = validate_result(result, target_size=target_size)
    return [i for i in report["issues"] if i["code"] == "set_size_mismatch"]


def test_a_panel_made_smaller_on_request_is_not_an_error():
    reduced = _Result(
        primers=("AAAACCCC", "CCCCGGGG"),
        stage_history=({"stage": "reduction", "before_size": 4, "after_size": 2},),
    )

    (issue,) = _size_issue(reduced, 4)

    assert issue["level"] == "warning"
    assert "minimi" in issue["detail"]


def test_a_short_panel_nobody_asked_to_shorten_is_still_an_error():
    (issue,) = _size_issue(_Result(primers=("AAAACCCC", "CCCCGGGG")), 4)

    assert issue["level"] == "error"


def test_a_reduction_stage_that_removed_nothing_excuses_nothing():
    unchanged = _Result(
        primers=("AAAACCCC", "CCCCGGGG"),
        stage_history=({"stage": "reduction", "before_size": 2, "after_size": 2},),
    )

    (issue,) = _size_issue(unchanged, 4)

    assert issue["level"] == "error"
