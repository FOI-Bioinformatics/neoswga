"""A failure record is never written by creating a directory at a file path.

`_failure_artifact_path` used to read `args.data_dir or args.output` and then
`os.makedirs` the result. `--output` does not mean the same thing across
commands: `analyze-coverage` documents it as "Output directory for gap
BED/JSON", while `calibrate-reach` documents it as "Write the fit as JSON to
this path", and `predict` and `report` name a file too.

So a design failure under any of the file-valued ones left a DIRECTORY at the
path the user had asked for a file. That is destructive twice over:

- the next run, the successful one after the user fixes the problem, cannot
  write its own output there any more; and
- the record lands beside nothing, where `export.export_is_blocked` does not
  look, so the record whose whole purpose is to block a stale export blocks
  nothing.

`--output` is now consulted last and only when it is already a directory. The
params file's `data_dir` is a definite statement about where this run's
outputs live; `--output` is a flag whose meaning varies, so it is the weaker
evidence and goes last.

`calibrate-reach` is the one that made this reachable from the BAM work: it
takes `--bam`, `--output` names a JSON file, and a CRAM whose reference cannot
be resolved raises `ReferenceDataError` from inside it.
"""

import argparse
import json
import os

import pytest

from neoswga.cli._failure import _failure_artifact_path


@pytest.fixture(autouse=True)
def no_inherited_run_directory(monkeypatch):
    """`parameter.data_dir` is a module global any earlier test may have set,
    and the resolver reads it. Without this these tests would pass or fail on
    what else ran in the process."""
    from neoswga.core import parameter

    monkeypatch.setattr(parameter, "data_dir", None, raising=False)


def test_a_file_valued_output_is_not_turned_into_a_directory(tmp_path, monkeypatch):
    """The defect, stated as the thing a user would find on disk."""
    monkeypatch.chdir(tmp_path)

    _failure_artifact_path(argparse.Namespace(output="fit.json", json_file=None))

    assert not os.path.isdir("fit.json"), "a directory was created where a file was asked for"
    assert not os.path.exists("fit.json"), "nothing at all should be created there"


def test_the_params_file_data_dir_wins_over_output(tmp_path):
    """Where the record is USEFUL. `export_is_blocked` reads the run directory,
    so a record anywhere else is written and never read."""
    data_dir = tmp_path / "run"
    data_dir.mkdir()
    params = tmp_path / "params.json"
    params.write_text(json.dumps({"data_dir": str(data_dir)}))

    path = _failure_artifact_path(
        argparse.Namespace(output=str(tmp_path / "fit.json"), json_file=str(params))
    )

    assert path == str(data_dir / "design_failure.json")
    assert not os.path.exists(tmp_path / "fit.json")


def test_an_output_that_is_already_a_directory_is_still_used(tmp_path, monkeypatch):
    """The `analyze-coverage` shape, which was never broken and must not
    regress: `--output` there names a directory and is the right place when
    nothing else is known."""
    monkeypatch.chdir(tmp_path)
    (tmp_path / "gaps").mkdir()

    path = _failure_artifact_path(argparse.Namespace(output="gaps", json_file=None))

    assert path == os.path.join("gaps", "design_failure.json")


def test_no_run_directory_at_all_returns_none_rather_than_littering(tmp_path, monkeypatch):
    """A command with nowhere to write still fails loudly on stderr and exits
    nonzero. Writing into the working directory instead would litter."""
    monkeypatch.chdir(tmp_path)

    assert _failure_artifact_path(argparse.Namespace(output=None, json_file=None)) is None
    assert list(tmp_path.iterdir()) == []
