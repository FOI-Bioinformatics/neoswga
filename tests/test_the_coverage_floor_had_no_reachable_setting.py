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

**The gap this leaves is narrower than it was, and is recorded rather than
filled.** No configured panel limit is read against the AGGREGATE foreground
coverage (`fg_coverage`, `effective_fg_coverage`, `headline_coverage`), so
there is no key that says "refuse a panel below X coverage" over several
targets pooled. Adding one is a separate decision, and this project does not
pick a threshold without evidence -- no spacing threshold derived from the
polymerase reach separates the 18 published sets with wet-lab outcomes, which
is why every panel limit ships unset.

Narrowed on 2026-10-02, not loosened. `min_per_target_coverage` became an
enforced limit in Phase 6c of the diversity plan, and it is a minimum coverage
-- per TARGET. It is read against `worst_target_coverage`, the record-confined
minimum over the foreground targets, and it has a measured loss behind it
(`docs/validation/2026-10-02-diversity-baseline.md` Q2 and Q3) rather than a
threshold someone picked. It also ships unset, so it refuses nothing by
default.

**On a SINGLE-TARGET design the two quantities coincide**, and the prose here
used to deny it. With one foreground prefix `worst_target_coverage` equals that
prefix's own coverage, which (absent record starts, so on every single-record
reference) equals `fg_coverage`. So `min_per_target_coverage: 0.7` on a
single-genome run IS the retired `min_coverage` floor at its old default, the
one whose verdict change this file says needs its own measurement. The
difference is that the user set it rather than the tool defaulting it, and the
schema and the CLI help now say so. What remains genuinely absent is a floor on
the POOLED figure of a multi-target design, which is what the assertion below
tests: no limit reads an aggregate coverage field. That is asserted on the
metrics field rather than the key's spelling, because the retired floor was
spelled `min_coverage` and would have slipped past a name check.
"""

import inspect
import pathlib

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


def test_no_configured_limit_is_an_aggregate_foreground_coverage_floor():
    """The gap, asserted so that filling it is a deliberate act.

    Asserted on the METRICS FIELD each limit is read against, not on the key's
    spelling: the retired floor was `min_coverage` and a name check would have
    let it back in. If someone adds a limit against the aggregate foreground
    coverage, this fails and points at this file's docstring, which says what
    the evidence bar is.
    """
    from neoswga.core.pool_objective import _LIMITS

    aggregate_fields = {"fg_coverage", "effective_fg_coverage", "headline_coverage"}
    offenders = [field for field, metric, *_rest in _LIMITS if metric in aggregate_fields]

    assert not offenders, offenders


def test_the_only_minimum_coverage_limit_is_the_per_target_one():
    """Guard the guard above, which would also pass with no limits at all.

    It also pins what Phase 6c added, so removing `min_per_target_coverage`
    from the enforced limits is as deliberate an act as adding one was.
    """
    minimums = {key for key in LIMIT_KEYS if "coverage" in key and key.startswith("min_")}

    assert minimums == {"min_per_target_coverage"}


def test_the_single_target_equivalence_is_documented_where_a_user_reads_it():
    """The prose above used to deny this. On one target the per-target floor IS
    a foreground coverage floor, and silence about that is what let the retired
    verdict change back in without saying so."""
    import json

    schema = json.loads(
        (pathlib.Path("neoswga") / "core" / "schema" / "params.schema.json").read_text()
    )
    description = schema["properties"]["min_per_target_coverage"]["description"].lower()

    assert "single-target" in description
    assert "foreground coverage" in description
