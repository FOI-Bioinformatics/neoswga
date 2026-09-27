"""A failed write must leave the previous result, not half of the new one.

An output directory is the only thing a later command sees. `export` and the
report read these files back, and neither can tell a truncated file from a
short one: a summary JSON cut off mid-object fails to parse, which the report
treats as "no optimizer summary" and answers from its own estimate instead, and
a results CSV cut off mid-row is a shorter panel that parses cleanly.

The same rule `position_index` follows, for the same reason.
"""

import json
import os

import pytest

from neoswga.core.atomic_output import (
    atomic_write_dataframe,
    atomic_write_json,
    atomic_write_text,
)

pd = pytest.importorskip("pandas")


def test_a_written_file_is_complete(tmp_path):
    path = str(tmp_path / "summary.json")
    atomic_write_json(path, {"metrics": {"fg_coverage": 0.5}}, indent=2)
    assert json.loads(open(path).read())["metrics"]["fg_coverage"] == 0.5


def test_a_failed_write_leaves_the_previous_file_intact(tmp_path):
    """The case that matters: the run that failed is the one whose directory
    must not be left describing nothing."""
    path = str(tmp_path / "summary.json")
    atomic_write_json(path, {"generation": 1})

    class Unserialisable:
        pass

    with pytest.raises(TypeError):
        atomic_write_json(path, {"generation": 2, "bad": Unserialisable()})

    assert json.loads(open(path).read()) == {"generation": 1}


def test_a_failed_write_leaves_no_temporary_file_behind(tmp_path):
    path = str(tmp_path / "summary.json")

    def explode(handle):
        handle.write('{"partial": ')
        raise OSError("disk full")

    with pytest.raises(OSError, match="disk full"):
        atomic_write_text(path, explode)

    assert os.listdir(tmp_path) == [], "a partial file must not survive under any name"


def test_nothing_partial_is_ever_visible_at_the_target_path(tmp_path):
    """The target must not exist at all until the content is complete."""
    path = str(tmp_path / "summary.json")
    seen = []

    def render(handle):
        handle.write('{"a": 1}')
        seen.append(os.path.exists(path))

    atomic_write_text(path, render)
    assert seen == [False], "the target existed while the write was in progress"
    assert os.path.exists(path)


def test_a_dataframe_is_published_atomically(tmp_path):
    path = str(tmp_path / "step4.csv")
    frame = pd.DataFrame({"primer": ["ACGT", "TTTT"], "set_index": [0, 0]})
    atomic_write_dataframe(path, frame, index=False)

    back = pd.read_csv(path)
    assert list(back["primer"]) == ["ACGT", "TTTT"]


def test_an_absent_directory_raises_rather_than_being_invented(tmp_path):
    """Writing must not create the location it was pointed at.

    An earlier draft of this helper called `makedirs`, and a shipped test
    caught it: `test_a_missing_verdict_is_not_a_passing_one` requires the
    validation record's write to fail loudly when its directory is missing.
    Inventing the location is how a failure record once left a DIRECTORY where
    the user had asked for a file.
    """
    path = str(tmp_path / "nested" / "deeper" / "summary.json")
    with pytest.raises(OSError):
        atomic_write_json(path, {"ok": True})
    assert not os.path.exists(path)


def test_step4_does_not_write_results_in_place():
    """A source check, because a new in-place write is invisible until a run
    fails at exactly the wrong moment."""
    import pathlib

    source = pathlib.Path("neoswga/core/step4_output.py").read_text()
    assert "df.to_csv(output_path" not in source
    assert "json.dump(summary" not in source
    assert "atomic_write_dataframe(" in source
    assert "atomic_write_json(" in source
