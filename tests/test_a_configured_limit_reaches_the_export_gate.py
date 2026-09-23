"""A configured limit the delivered panel misses must stop the export.

`optimize` already knows. On the packaged plasmid example with
`min_selectivity_density` set to 1e9, the run prints its "Configured limits"
table with `NO` beside the limit, prints "Repaired: no, the panel returned is
the one selected", and then recommends `neoswga export`. Export exits 0 and
prints "Primers ready for ordering!".

Two defects in `panel_acceptance.apply_configured_limits` put it there.

**The report is built only under `verbose`.** The `AcceptanceReport` carrying
`violations`, `shortfall` and the per-limit `values` is constructed as an
argument to `report_acceptance` inside `if verbose:`. On a non-verbose run the
violations are never computed at all, so there is nothing to record even in
principle.

**`return None` conflates two outcomes.** It returns None when
`updated is result`, and that module's own docstring states the rule that makes
the two identical: "A repair that does not resolve the violation returns its
input." So the caller cannot tell a panel that met every limit from one that
failed a limit the repair could not fix.

Nothing was serialized, and `export_is_blocked` selects issues whose `code` is
in `BLOCKING_VALIDATOR_CODES`, which held one member.

**Why not simply let the assessment decide.** `evaluate_panel` already records
`qualified: False` with `violations: ['selectivity below minimum', 'panel size 1
is below the requested 6']` for this very run. Blocking on that whole list would
also refuse every panel shorter than requested, and a short panel is a normal,
documented outcome -- `num_primers` is a request, not a guarantee. The
configured limits are narrower and different in kind: the user set them. So this
records the acceptance path's own list, scoped by `constraints`, and leaves the
assessment's broader list to be measured before anything gates on it.
"""

import json
import shutil
import subprocess
import sys

import pytest

from neoswga.core.design_result import VALIDATION_FILENAME
from tests.conftest import _EXAMPLE_DIR, plasmid_example_ready

pytestmark = pytest.mark.skipif(
    not plasmid_example_ready(), reason="needs jellyfish and the plasmid example"
)

#: Far outside anything this 6 kb pair can deliver, so the repair cannot
#: succeed and the panel is returned unchanged. That is the case the export
#: gate could not see.
UNREACHABLE_DENSITY = 1e9


def _run(tmp_path, **overrides):
    """Run `optimize` over a private copy of the primed example."""
    workdir = tmp_path / "plasmid"
    shutil.copytree(_EXAMPLE_DIR, workdir)

    params_path = workdir / "params.json"
    params = json.loads(params_path.read_text())
    params.update(overrides)
    params_path.write_text(json.dumps(params, indent=2))

    completed = subprocess.run(
        [sys.executable, "-m", "neoswga.cli_unified", "optimize", "-j", "params.json"],
        cwd=workdir,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        pytest.fail(f"optimize failed:\n{completed.stderr[-2000:]}")
    return workdir


def _issue_codes(workdir):
    payload = json.loads((workdir / VALIDATION_FILENAME).read_text())
    return [issue.get("code") for issue in payload.get("issues") or []]


def test_a_panel_that_misses_a_configured_limit_is_not_exportable(tmp_path):
    """The headline: this exported under "Primers ready for ordering!"."""
    from neoswga.core.export import export_is_blocked

    workdir = _run(tmp_path, min_selectivity_density=UNREACHABLE_DENSITY)

    blocked = export_is_blocked(str(workdir))
    assert blocked, "a panel that missed its configured limit was exportable"
    assert "selectivity" in blocked.lower(), blocked


def test_the_recorded_issue_names_the_limit(tmp_path):
    """A refusal a user cannot act on is barely better than none."""
    workdir = _run(tmp_path, min_selectivity_density=UNREACHABLE_DENSITY)

    payload = json.loads((workdir / VALIDATION_FILENAME).read_text())
    issues = [i for i in payload["issues"] if i["code"] == "panel_limit_not_met"]

    assert issues, _issue_codes(workdir)
    assert issues[0]["level"] == "error"
    assert "selectivity" in issues[0]["detail"].lower(), issues[0]["detail"]


def test_setting_no_limit_leaves_the_export_alone(tmp_path):
    """Guard the guard, and the documented contract.

    Set no panel limit and `constraints_from_parameter` returns None, no
    objective is built and nothing here runs. A change that blocked this run
    would be refusing designs on a limit nobody set.
    """
    from neoswga.core.export import export_is_blocked

    workdir = _run(tmp_path)

    assert "panel_limit_not_met" not in _issue_codes(workdir)
    assert export_is_blocked(str(workdir)) is None


def test_a_limit_the_panel_meets_does_not_block(tmp_path):
    """The other half of the guard: a configured limit is not a refusal by
    itself. Only missing one is."""
    from neoswga.core.export import export_is_blocked

    workdir = _run(tmp_path, min_selectivity_density=0.0001)

    assert "panel_limit_not_met" not in _issue_codes(workdir)
    assert export_is_blocked(str(workdir)) is None
