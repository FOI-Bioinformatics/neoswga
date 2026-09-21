"""The candidate-preparation stage is no longer called `score`.

Task 3 of the 2026-09-21 valid-design plan. The stage stopped scoring anything
on 2026-09-05, when the bundled amplification model was retired from the
default path: it was computing a prediction for every candidate and discarding
it, and every step-4 consumer reads only the primer column.

What it does is prepare `step3_df.csv`, the ordered candidate pool six modules
read. The name said otherwise, and a name that describes work the command does
not do is the same defect as a flag that is read and then overruled. Both make
someone reason about a pipeline they do not have.

No alias is kept. The plan's first global constraint is that backward
compatibility is not a requirement, and an alias would leave the misleading
name reachable and in every example someone copies.
"""

import subprocess
import sys

import pytest

from neoswga import cli_unified


def _run(*args):
    return subprocess.run(
        [sys.executable, "-m", "neoswga.cli_unified", *args],
        capture_output=True,
        text=True,
    )


def _subcommands():
    parser = cli_unified.create_parser()
    for action in parser._subparsers._group_actions:
        return set(action.choices)
    return set()


def test_the_command_is_named_for_what_it_does():
    assert "prepare-candidates" in _subcommands()


def test_the_old_name_is_gone_rather_than_aliased():
    """An alias would leave the misleading name in every copied example."""
    assert "score" not in _subcommands()


def test_the_old_name_fails_with_a_message_naming_the_new_one(tmp_path):
    result = _run("score", "-j", str(tmp_path / "params.json"))

    assert result.returncode != 0
    assert "prepare-candidates" in (result.stdout + result.stderr)


def test_it_is_dispatched(tmp_path):
    """A parser entry with no dispatch entry is a command that cannot run."""
    result = _run("prepare-candidates", "--help")

    assert result.returncode == 0
    assert "candidate pool" in result.stdout.lower()


def test_the_one_line_help_does_not_promise_scoring():
    """The line shown in the command list, not the full description.

    The description deliberately explains the old name and why it was wrong,
    which is the part a reader needs when an old example fails. The one-line
    help is what someone skims, and it must describe the work.
    """
    parser = cli_unified.create_parser()
    for action in parser._subparsers._group_actions:
        described = " ".join(
            choice.help or ""
            for choice in action._choices_actions
            if choice.dest == "prepare-candidates"
        )
        break
    assert described
    assert "score" not in described.lower(), described
    assert "candidate pool" in described.lower(), described


def test_the_retired_fast_score_flag_is_gone(tmp_path):
    """It was documented as accepted and doing nothing, which is Known Issue 8."""
    result = _run("prepare-candidates", "-j", str(tmp_path / "p.json"), "--fast-score")

    assert result.returncode != 0
    assert "--fast-score" in (result.stdout + result.stderr)


@pytest.mark.parametrize("flag", ["--amp-model", "--min-amp-pred"])
def test_the_flags_that_restore_the_old_behaviour_remain(flag):
    """`--amp-model` is not retired: it restores the score column and the gate."""
    parser = cli_unified.create_parser()
    for action in parser._subparsers._group_actions:
        options = {
            option
            for sub_action in action.choices["prepare-candidates"]._actions
            for option in sub_action.option_strings
        }
        break
    assert flag in options


def test_the_pipeline_documents_the_new_name():
    """No runnable example may still use the old name.

    Prose about the rename is fine and wanted: someone whose script broke needs
    to find out why. What must not survive is a line someone copies, so the
    check is on lines that read as commands rather than on the string anywhere.
    """
    import pathlib

    for relative in ("README.md", "CLAUDE.md"):
        text = (pathlib.Path(__file__).resolve().parent.parent / relative).read_text()
        assert "neoswga prepare-candidates" in text, relative
        runnable = [
            line
            for line in text.splitlines()
            if line.lstrip().lstrip("$ ").startswith("neoswga score")
        ]
        assert not runnable, (relative, runnable)
