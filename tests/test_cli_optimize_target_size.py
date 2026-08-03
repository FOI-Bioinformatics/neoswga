"""params.json's requested set size has to survive to the optimizer.

`optimize` resolved its target primer count from `parameter.num_primers` before
anything had read params.json. That module global carries a default of 6 and is
only overwritten when `get_params()` runs, which happens lazily inside
`optimize_step4` -- after the target has already been resolved and passed down
as an explicit kwarg. So a params.json asking for 8 primers silently got 6.

The failure was quiet in the worst way. The completion check further down reads
`parameter.target_set_size`, and by then get_params() *has* run, so it compared
the 6 primers produced against the 8 requested and reported:

    WARNING: Found 6 primers but target was 8 (PARTIAL result: insufficient
    candidates)
    To improve, consider:
      - Relaxing filter thresholds (max_bg_freq, max_gini)
      - Increasing candidate pool (max_primer)

None of which was the problem. The pool held 2000 candidates and the optimizer
had simply been asked for 6. Acting on that advice -- widening the filters to
find primers that were never missing -- would degrade the design.

Found while running the pipeline at three primer lengths, where all three runs
reported the same false shortfall.
"""

import json

import pytest

from neoswga.cli.pipeline import _target_size_from_params


class FakeParameter:
    """Stands in for the parameter module, which is a module of globals."""

    def __init__(self, json_file=None, json_data=None):
        self.json_file = json_file
        self._json_data = json_data if json_data is not None else {}
        # The default that used to win over params.json.
        self.num_primers = 6
        self.target_set_size = 6


@pytest.fixture
def params_file(tmp_path):
    def write(payload):
        path = tmp_path / "params.json"
        path.write_text(json.dumps(payload))
        return str(path)

    return write


# ----------------------------------------------------------------------
# The regression
# ----------------------------------------------------------------------


def test_num_primers_from_params_json_is_used(params_file):
    """The bug, at its smallest: _json_data is empty at this point in the
    command, which is exactly when the old code consulted it."""
    parameter = FakeParameter(params_file({"num_primers": 8}))

    assert _target_size_from_params(parameter) == 8


def test_target_set_size_is_honoured_when_num_primers_is_absent(params_file):
    """Both spellings appear in params.json files in this repo and in the
    docs, and `get_params` accepts either, so this path has to as well."""
    parameter = FakeParameter(params_file({"target_set_size": 12}))

    assert _target_size_from_params(parameter) == 12


def test_num_primers_wins_over_target_set_size(params_file):
    """Matching parameter.get_params, which resolves them in this order. Two
    resolution orders for one value is how they come to disagree."""
    parameter = FakeParameter(params_file({"num_primers": 8, "target_set_size": 12}))

    assert _target_size_from_params(parameter) == 8


def test_an_already_loaded_json_is_preferred_over_re_reading(params_file):
    """When get_params() has run, `_json_data` is authoritative -- it is what
    the rest of the pipeline sees. Re-reading the file could differ if it were
    edited mid-run, and the two halves of one run disagreeing about the set
    size is worse than either answer."""
    parameter = FakeParameter(params_file({"num_primers": 8}), json_data={"num_primers": 20})

    assert _target_size_from_params(parameter) == 20


# ----------------------------------------------------------------------
# It falls back rather than failing
# ----------------------------------------------------------------------


def test_a_missing_params_file_falls_back_to_the_default():
    """`optimize` can run without -j. get_params() would raise here in a
    moment with a far better message than this helper could."""
    parameter = FakeParameter(json_file=None)

    assert _target_size_from_params(parameter) == 6


def test_an_unreadable_params_file_does_not_raise_here(tmp_path):
    """Pre-empting get_params()'s error with a partial one from the set-size
    lookup would report the wrong problem."""
    path = tmp_path / "params.json"
    path.write_text("{not json")
    parameter = FakeParameter(str(path))

    assert _target_size_from_params(parameter) == 6


@pytest.mark.parametrize("bad", ["eight", -1, 0, None, 3.5, True])
def test_a_nonsense_value_falls_back_and_says_so(params_file, caplog, bad):
    """A hand-edited params.json should not produce a quietly different set
    size. `True` is in the list because bool is an int in Python and would
    otherwise resolve to a one-primer set."""
    import logging

    parameter = FakeParameter(params_file({"num_primers": bad}))
    with caplog.at_level(logging.WARNING):
        assert _target_size_from_params(parameter) == 6

    assert "num_primers" in caplog.text
