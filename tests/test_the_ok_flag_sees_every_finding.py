"""`validation["ok"]` must agree with the findings recorded beside it.

The saved validation file is read by three commands that tell a user whether a
pool may be ordered, and until now it could say `ok: true` while carrying a
finding `export` refuses on. Two independent mechanisms produced that, and
either one alone was enough.

**`ok` was computed before the issues were complete.**
`OptimizationResult.validate` folded `level == "error"` into `ok` and returned
the dict. `unified_optimizer` then appended the saturation warnings and the
delivered-panel dimer finding to that same dict. A finding appended after the
fold could not move the flag whatever level it carried, so the dimer breach was
invisible to `ok` by construction rather than by policy.

**A blocking code could be emitted at warning level.** It was:
`dimer_validation_issue` recorded `level="warning"` while
`BLOCKING_VALIDATOR_CODES` held its code. So `export` refused the pool,
`interpret` called it unfit, and the report rendered the same file as having no
errors.

This is the silent-zero family in its report-facing shape: not an absent
measurement standing in as zero, but a summary flag computed from a subset of
its own evidence and defaulting to the favourable answer.

`panel_validation_is_ok` checks both conditions, so the level and the code set
cannot drift apart again, and it is applied to the FINAL issue list.
"""

import ast
import pathlib

import pytest

from neoswga.core.design_result import BLOCKING_VALIDATOR_CODES, panel_validation_is_ok

ROOT = pathlib.Path(__file__).resolve().parent.parent


def issue(code, level="warning"):
    return {"level": level, "code": code, "detail": f"{code} happened"}


# ---------------------------------------------------------------------------
# The rule
# ---------------------------------------------------------------------------


def test_a_clean_pool_is_ok():
    assert panel_validation_is_ok([]) is True
    assert panel_validation_is_ok([issue("coverage_saturated_on_small_genome")]) is True


def test_an_error_blocks_whatever_its_code():
    assert panel_validation_is_ok([issue("something_new", level="error")]) is False


@pytest.mark.parametrize("code", sorted(BLOCKING_VALIDATOR_CODES))
def test_a_blocking_code_blocks_whatever_its_level(code):
    """The half that was wrong: a code `export` refuses on left `ok` true.

    Parametrised over the set rather than the one known code, so a future
    blocking code emitted at warning level is caught the day it is added.
    """
    assert panel_validation_is_ok([issue(code)]) is False


def test_saturation_still_leaves_a_pool_ok():
    """Deliberate, and the reason `ok` is not simply "no issues".

    It says a metric cannot be trusted on a small target, which is inherent to
    designing against a plasmid. Blocking on it would refuse every such design.
    """
    assert "coverage_saturated_on_small_genome" not in BLOCKING_VALIDATOR_CODES
    assert panel_validation_is_ok([issue("coverage_saturated_on_small_genome")]) is True


def test_a_malformed_entry_does_not_crash_the_verdict():
    assert panel_validation_is_ok([None, "text", issue("x")]) is True


# ---------------------------------------------------------------------------
# Both ends apply it
# ---------------------------------------------------------------------------


def test_the_validator_reports_a_dimer_breach_as_an_error():
    """`export` and `interpret` treat this code as blocking. A report rendering
    the same file as error-free is two commands disagreeing about one panel."""
    from neoswga.core.dimer import dimer_validation_issue

    # A pair with a long complementary run; any threshold below it will do.
    found = dimer_validation_issue(["AAAAAAAAAAAA", "TTTTTTTTTTTT"], max_dimer_bp=3)

    assert found is not None, "the fixture no longer breaks the threshold"
    assert found["code"] in BLOCKING_VALIDATOR_CODES
    assert found["level"] == "error"


def test_validate_builds_ok_through_the_shared_rule():
    """Driving a real result rather than reading source text: a test that
    matched a call by name would pass on a call that computed the wrong list."""
    from neoswga.core.base_optimizer import OptimizationResult, OptimizationStatus, PrimerSetMetrics

    result = OptimizationResult(
        primers=("ACGTACGTACGT", "ACGTACGTACGT"),
        score=1.0,
        status=OptimizationStatus.SUCCESS,
        metrics=PrimerSetMetrics.empty(),
        iterations=1,
        optimizer_name="test",
    )
    report = result.validate(target_size=2)

    assert report["ok"] is panel_validation_is_ok(report["issues"])
    assert report["ok"] is False, "a duplicated primer is an error-level finding"


def test_the_flag_is_recomputed_after_the_appends():
    """The defect that made the level irrelevant.

    `run_optimization` appends to `validation["issues"]` after `validate`
    returned, so the write site must recompute. Asserted by AST because the
    property is an ORDERING -- that the recompute follows the appends -- and a
    behavioural test would need a full optimizer run to reach it.
    """
    source = (ROOT / "neoswga" / "core" / "unified_optimizer.py").read_text()
    tree = ast.parse(source)

    appends, recompute = [], []
    for node in ast.walk(tree):
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute):
            if node.func.attr in {"append", "extend"} and "validation" in ast.dump(node.func):
                appends.append(node.lineno)
        if isinstance(node, ast.Call) and getattr(node.func, "id", None) == "panel_validation_is_ok":
            recompute.append(node.lineno)

    assert appends, "no issue is appended after validate(); this test is stale"
    assert recompute, "run_optimization no longer recomputes ok from the final list"
    assert max(recompute) > max(appends), (
        "ok is recomputed before the last append, so a finding recorded after "
        "it cannot move the flag -- which is the defect this closes"
    )
