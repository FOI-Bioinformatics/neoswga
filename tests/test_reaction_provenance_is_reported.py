"""A panel selected before the reaction reached the optimizer is reported.

Until 2026-10-26 `unified_optimizer` passed neither `polymerase` nor
`reaction_temp` to any optimizer, and both hybrid wrappers defaulted
`polymerase` to the literal "phi29", so every `hybrid` and `background-aware`
run built a phi29 preset at 30 C whatever params.json declared
(`docs/validation/2026-10-26-hybrid-ran-every-reaction-as-phi29.md`).

The code is fixed. These tests are about the designs that were already
produced: the manifest records the polymerase and temperature each step ran
under, so whether a given `data_dir` holds an affected panel is answerable, and
`doctor` answers it. A fix that leaves no way to tell which outputs predate it
leaves the user holding a panel they cannot evaluate.
"""

import io
import json
from contextlib import redirect_stdout

import pytest

from neoswga.core.run_manifest import steps_predating_the_reaction_fix

PRE_FIX = "2026-10-20T12:00:00+00:00"
POST_FIX = "2026-10-28T12:00:00+00:00"


def _manifest(tmp_path, entries):
    (tmp_path / "run_manifest.json").write_text(json.dumps({"steps": entries}))
    return str(tmp_path)


def _step(step="optimize", polymerase="bst", temp=65.0, stamp=PRE_FIX):
    conditions = {}
    if polymerase is not None:
        conditions["polymerase"] = polymerase
    if temp is not None:
        conditions["reaction_temp"] = temp
    return {"step": step, "timestamp_utc": stamp, "effective_conditions": conditions}


# ----------------------------------------------------------------------
# What must be reported
# ----------------------------------------------------------------------


@pytest.mark.parametrize(
    "polymerase,temp,screen_skipped",
    [
        # equiphi29 and bst set thermo_filter=True, so `optimize` skipped its
        # whole Stage-0 pre-screen when it resolved them as phi29.
        ("bst", 65.0, True),
        ("equiphi29", 45.0, True),
        # klenow's preset sets thermo_filter=False like phi29's, so no screen
        # was skipped -- but its 37 C reaction still ran at 30 C.
        ("klenow", 37.0, False),
    ],
)
def test_a_pre_fix_hot_polymerase_step_is_reported(tmp_path, polymerase, temp, screen_skipped):
    found = steps_predating_the_reaction_fix(
        _manifest(tmp_path, [_step(polymerase=polymerase, temp=temp)])
    )
    assert len(found) == 1, f"{polymerase} at {temp}C should be reported: {found}"
    assert found[0]["polymerase"] == polymerase
    assert found[0]["temp_intended"] == temp
    assert found[0]["temp_used"] == 30.0
    assert found[0]["screen_skipped"] is screen_skipped, (
        f"{polymerase}'s preset decides whether the Stage-0 screen was skipped; "
        f"expected {screen_skipped}"
    )


def test_a_pre_fix_phi29_run_at_a_configured_temperature_is_reported(tmp_path):
    """The non-obvious case, and the one a polymerase-only check would miss.

    An explicitly configured `reaction_temp` did not reach the optimizer either,
    so a phi29 design at any temperature other than 30 C was also selected under
    the wrong reaction -- no screen skipped, but Stage 2's dimer screen and Tm
    term both ran 5 C cold.
    """
    found = steps_predating_the_reaction_fix(
        _manifest(tmp_path, [_step(polymerase="phi29", temp=35.0)])
    )
    assert len(found) == 1, "a configured 35 C did not reach the optimizer either"
    assert found[0]["screen_skipped"] is False
    assert (found[0]["temp_intended"], found[0]["temp_used"]) == (35.0, 30.0)


def test_expand_primers_is_reported_too(tmp_path):
    """`expand-primers` also selects through the optimizer."""
    found = steps_predating_the_reaction_fix(_manifest(tmp_path, [_step(step="expand-primers")]))
    assert len(found) == 1 and found[0]["step"] == "expand-primers"


# ----------------------------------------------------------------------
# What must stay silent -- a provenance warning that cries wolf is ignored
# ----------------------------------------------------------------------


def test_a_pre_fix_phi29_design_at_30c_is_silent(tmp_path):
    """The default configuration is the one case the defect left correct.

    phi29 at 30 C is exactly what every pre-fix run resolved to, so those
    designs got the reaction they asked for and must not be flagged.
    """
    assert (
        steps_predating_the_reaction_fix(
            _manifest(tmp_path, [_step(polymerase="phi29", temp=30.0)])
        )
        == []
    )


def test_an_absent_polymerase_reads_as_the_default_and_is_silent(tmp_path):
    """Absence means the default, which was correct -- not unknown."""
    assert (
        steps_predating_the_reaction_fix(_manifest(tmp_path, [_step(polymerase=None, temp=None)]))
        == []
    )


def test_a_post_fix_step_is_silent(tmp_path):
    assert steps_predating_the_reaction_fix(_manifest(tmp_path, [_step(stamp=POST_FIX)])) == []


def test_steps_the_optimizer_does_not_decide_are_silent(tmp_path):
    """`filter` resolved its own Tm window from conditions and was unaffected."""
    entries = [_step(step=s) for s in ("count-kmers", "filter", "prepare-candidates")]
    assert steps_predating_the_reaction_fix(_manifest(tmp_path, entries)) == []


def test_an_unreadable_timestamp_is_left_out_rather_than_assumed_old(tmp_path):
    """Unknown is not old. A fabricated date would be a fabricated finding."""
    assert (
        steps_predating_the_reaction_fix(_manifest(tmp_path, [_step(stamp="not-a-timestamp")]))
        == []
    )


def test_a_missing_or_unreadable_manifest_reports_nothing(tmp_path):
    assert steps_predating_the_reaction_fix(str(tmp_path)) == []
    (tmp_path / "run_manifest.json").write_text("{not json")
    assert steps_predating_the_reaction_fix(str(tmp_path)) == []


# ----------------------------------------------------------------------
# The user-facing surface
# ----------------------------------------------------------------------


def _doctor_json(params_path):
    from neoswga.cli.setup import run_doctor

    class _Args:
        json_file = params_path
        output_json = True

    buf = io.StringIO()
    try:
        with redirect_stdout(buf):
            run_doctor(_Args())
    except SystemExit:
        pass
    return json.loads(buf.getvalue())


def test_doctor_reports_an_affected_design(tmp_path):
    data_dir = _manifest(tmp_path, [_step()])
    params = tmp_path / "params.json"
    params.write_text(json.dumps({"data_dir": data_dir, "fg_prefixes": []}))

    report = _doctor_json(str(params))

    assert report.get("reaction_provenance"), "doctor did not surface the affected step"
    joined = " ".join(report.get("warnings", []))
    assert "reached the optimizer" in joined, joined
    assert "phi29 at 30C" in joined, "the warning must say what the run actually used"


def test_doctor_is_silent_on_an_unaffected_design(tmp_path):
    data_dir = _manifest(tmp_path, [_step(polymerase="phi29", temp=30.0)])
    params = tmp_path / "params.json"
    params.write_text(json.dumps({"data_dir": data_dir, "fg_prefixes": []}))

    report = _doctor_json(str(params))

    assert "reaction_provenance" not in report
    assert not [w for w in report.get("warnings", []) if "reached the optimizer" in w]
