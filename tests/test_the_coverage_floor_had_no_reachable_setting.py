"""`validate_result`'s coverage floor could not fire, and is retired.

The check was `fg_cov < min_coverage`. `min_coverage` defaulted to 0.0, the one
production call site -- `unified_optimizer`, through `OptimizationResult.validate`
-- passed 0.0 with the comment "soft by default; caller can tighten", and no
caller tightened it. A coverage is never below zero, so `coverage_below_threshold`
could not be emitted on any run.

It was exercised, which is what made it look alive:
`tests/test_post_optimization_validator.py` passed real floors and the branch
behaved correctly. Known Issue 16 names this shape -- an argument no caller
supplies is invisible to every ratchet here, because they walk CLI options and
public capabilities, not parameters.

Retired rather than wired, for a reason that outlasts the reachability. A
coverage floor is an ACCEPTANCE criterion, and `validate_result` records
optimizer MISBEHAVIOUR: duplicates, drift from the requested size, blacklist
re-injection. A panel with low coverage is a design outcome, not a bug in the
search. Acceptance criteria have a home in this codebase --
`panel_acceptance.LIMIT_KEYS` and the objective `evaluate_panel` consults -- and
this was a second, quieter home that nobody could reach.

Wiring it to `optimize_step4`'s `target_coverage` was the alternative and is
not taken here: at its 0.70 default that newly fails validation on most measured
designs, which is a verdict change needing its own measurement.

**The gap this leaves is real and is recorded rather than filled.** None of the
six configured panel limits is a minimum foreground coverage: they are
`min_selectivity_density`, `max_background_sites`, `max_worst_hole`,
`max_mean_gap`, `max_evenness` and `max_host_coverage`. So after this there is
no way to say "refuse a panel below X coverage". Adding one is a separate
decision, and this project does not pick a threshold without evidence -- no
spacing threshold derived from the polymerase reach separates the 18 published
sets with wet-lab outcomes, which is why the other six ship unset.
"""

import inspect

from neoswga.core.base_optimizer import OptimizationResult, OptimizationStatus
from neoswga.core.panel_acceptance import LIMIT_KEYS
from neoswga.core.result_validation import validate_result


class _Metrics:
    fg_coverage = 0.0
    per_target_coverage: dict = {}


def _result(primers, status=OptimizationStatus.SUCCESS):
    return OptimizationResult(
        primers=list(primers),
        score=1.0,
        status=status,
        metrics=_Metrics(),
        iterations=1,
        optimizer_name="test",
    )


def test_the_floor_is_gone_from_both_signatures():
    """Both, because the delegate is where the default lived."""
    assert "min_coverage" not in inspect.signature(validate_result).parameters
    assert "min_coverage" not in inspect.signature(OptimizationResult.validate).parameters


def test_a_zero_coverage_panel_raises_no_coverage_issue():
    report = validate_result(_result(["ATCG", "GCTA", "TACG"]), target_size=3)

    codes = [issue["code"] for issue in report["issues"]]
    assert "coverage_below_threshold" not in codes


def test_the_checks_that_remain_still_fire():
    """Guard the guard. A validator that reports nothing would pass the test
    above, and the three checks below are the ones this record is for."""
    report = validate_result(_result(["ATCG", "ATCG"]), target_size=3, forbidden_primers=["ATCG"])

    codes = {issue["code"] for issue in report["issues"]}
    assert "duplicate_primers" in codes
    assert "set_size_mismatch" in codes
    assert "blacklist_primer_in_set" in codes
    assert report["ok"] is False


def test_no_configured_limit_is_a_foreground_coverage_floor():
    """The gap, asserted so that filling it is a deliberate act.

    If someone adds a minimum-coverage limit, this test fails and points them
    at this file's docstring, which says what the evidence bar is.
    """
    assert not [key for key in LIMIT_KEYS if "coverage" in key and key.startswith("min_")]
