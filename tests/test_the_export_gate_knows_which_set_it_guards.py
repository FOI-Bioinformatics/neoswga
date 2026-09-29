"""The export gate must judge the panel about to be written.

`step4_improved_df.csv` holds up to `max_sets` primer sets and they are
ALTERNATIVES: separate answers to the same design, never pooled. Everything that
assesses a run assesses set 0 -- the validator findings, the panel assessment,
the summary the report reads. Nothing evaluates the sets after it.

`export_is_blocked` took a directory and no set index, so it passed judgement on
set 0 while `export --set 4` wrote a different panel. A user could therefore
order an alternative that violates their own configured limits under "Primers
ready for ordering!". Measured on the Wolbachia design with a density floor of
20: set 4 came back at 17.05 and was offered
(docs/validation/alternatives_through_the_contract_2026-09-28.md).

Refusing an unassessed set rather than passing it is the rule this function
already applies to a failure record it cannot parse. Unknown is not success, and
here the unknown is total: nothing evaluated that panel at all.
"""

import json

import pytest

from neoswga.core.design_result import VALIDATION_FILENAME
from neoswga.core.export import export_is_blocked


@pytest.fixture
def clean_run(tmp_path):
    """A directory whose last run finished with a clean set 0."""
    (tmp_path / VALIDATION_FILENAME).write_text(json.dumps({"ok": True, "issues": []}))
    return tmp_path


def test_the_primary_set_is_judged_as_before(clean_run):
    assert export_is_blocked(str(clean_run)) is None
    assert export_is_blocked(str(clean_run), 0) is None


@pytest.mark.parametrize("set_index", [1, 2, 4])
def test_an_alternative_set_is_refused_because_nothing_assessed_it(clean_run, set_index):
    blocked = export_is_blocked(str(clean_run), set_index)

    assert blocked, "an unassessed alternative left the gate as though it were verified"
    assert f"Set {set_index}" in blocked
    assert "not been assessed" in blocked
    assert "--allow-unqualified" in blocked, "the refusal must name its override"


def test_the_refusal_says_which_panel_the_findings_describe(clean_run):
    """A user has to know the verdict they are missing, not just that one is."""
    blocked = export_is_blocked(str(clean_run), 3)
    assert "describe set 0" in blocked


def test_a_defect_in_set_zero_still_blocks_set_zero(tmp_path):
    """The existing gate is unchanged for the set it was written for."""
    (tmp_path / VALIDATION_FILENAME).write_text(
        json.dumps(
            {
                "ok": False,
                "issues": [
                    {
                        "level": "error",
                        "code": "duplicate_primers",
                        "detail": "ACGT appears twice",
                    }
                ],
            }
        )
    )
    blocked = export_is_blocked(str(tmp_path), 0)
    # The gate quotes the finding's DETAIL, which is what a reader can act on;
    # the code is an internal label.
    assert blocked and "ACGT appears twice" in blocked


def test_the_default_is_the_primary_so_existing_callers_are_unaffected():
    import inspect

    default = inspect.signature(export_is_blocked).parameters["set_index"].default
    assert default == 0


def test_the_command_passes_the_set_it_is_about_to_write():
    """A source check: the parameter is useless if the CLI does not supply it."""
    import pathlib

    source = pathlib.Path("neoswga/cli/report.py").read_text()
    assert 'export_is_blocked(args.dir, getattr(args, "set_index", 0) or 0)' in source, (
        "export must hand the gate the set it is exporting, or the gate judges "
        "a different panel from the one written"
    )
