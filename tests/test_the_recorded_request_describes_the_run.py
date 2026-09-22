"""The request recorded in the manifest must describe the run, not the file.

`request_hash` exists so "what was this designed under" has one answer to read
back. `design_request_for_run` reads the params FILE, and its docstring gives a
good reason: at that point in `run_step4` the `parameter` module has not been
resolved, because `get_params` runs inside `optimize_step4`.

But the panel size has already been overridden twice by then. `cli/pipeline.py`
sets `parameter.num_primers` from `args.num_primers` if given, and `--auto-size`
may replace it again with a recommendation, both BEFORE the request is built.
So `neoswga optimize -j params.json --num-primers 40` against a file saying 12
records the hash of a 12-oligo design and delivers a 40-oligo one.

The fix keeps both properties rather than trading one for the other:

- **Validation stays early and on the file.** A setting that cannot be applied
  is still named before the run has cost anything. Moving the only
  configuration gate past the position-cache build would mean a typo costs a
  full index.
- **Identity is recorded late, from the effective values.** A second resolve
  after the overrides, whose hash is what the manifest carries.

One resolver, two calls, two purposes.
"""

import json

import pytest

BASE = {
    "fg_genomes": ["a.fna"],
    "fg_prefixes": ["a"],
    "polymerase": "phi29",
    "reaction_temp": 30.0,
    "min_k": 12,
    "max_k": 12,
    "num_primers": 12,
}


@pytest.fixture
def params_file(tmp_path):
    path = tmp_path / "params.json"
    path.write_text(json.dumps(BASE))
    return path


class Args:
    """The attributes `run_step4` reads off its namespace, as a stand-in."""

    def __init__(self, json_file, **kwargs):
        self.json_file = str(json_file)
        for key, value in kwargs.items():
            setattr(self, key, value)


# ---------------------------------------------------------------------------
# What the file says
# ---------------------------------------------------------------------------


def test_the_file_request_is_unchanged_without_overrides(params_file):
    """Guard the guard: the two resolutions must agree when nothing overrides."""
    from neoswga.core.design_request import (
        design_request_for_run,
        effective_design_request,
    )

    args = Args(params_file)

    assert (
        effective_design_request(args, target_size=None).request_hash
        == design_request_for_run(args, None).request_hash
    )


# ---------------------------------------------------------------------------
# What the run does
# ---------------------------------------------------------------------------


def test_a_command_line_panel_size_reaches_the_recorded_request(params_file):
    """The defect. `--num-primers 40` against a file saying 12."""
    from neoswga.core.design_request import (
        design_request_for_run,
        effective_design_request,
    )

    args = Args(params_file)

    from_file = design_request_for_run(args, None)
    effective = effective_design_request(args, target_size=40)

    assert from_file.target_size == 12
    assert effective.target_size == 40, "the recorded request still says 12"
    assert effective.request_hash != from_file.request_hash


def test_an_auto_size_recommendation_reaches_the_recorded_request(params_file):
    """`--auto-size` replaces the size a second time, after the CLI value."""
    from neoswga.core.design_request import effective_design_request

    args = Args(params_file)

    assert effective_design_request(args, target_size=18).target_size == 18


def test_the_effective_request_still_refuses_what_the_file_request_refuses(tmp_path):
    """The late resolve must not become a way around the early gate."""
    from neoswga.core.design_request import effective_design_request
    from neoswga.core.exceptions import InvalidDesignRequest

    path = tmp_path / "params.json"
    path.write_text(json.dumps({**BASE, "concentration_mode": "fixed_total"}))

    with pytest.raises(InvalidDesignRequest):
        effective_design_request(Args(path), target_size=12)


def test_an_effective_size_of_zero_is_refused(params_file):
    """An auto-size recommending nothing is a failure, not a panel of none."""
    from neoswga.core.design_request import effective_design_request
    from neoswga.core.exceptions import InvalidDesignRequest

    with pytest.raises(InvalidDesignRequest):
        effective_design_request(Args(params_file), target_size=0)


def test_no_params_file_means_no_request(tmp_path):
    """The library and `--candidates` path. Absence, not a silent acceptance."""
    from neoswga.core.design_request import effective_design_request

    assert effective_design_request(Args(""), target_size=12) is None


# ---------------------------------------------------------------------------
# The gate stays where it is
# ---------------------------------------------------------------------------


def test_validation_still_happens_before_the_expensive_work():
    """`run_step4` must resolve the file BEFORE it builds the position cache.

    Asserted against the call order in the source rather than by timing. The
    point of validating early is that a setting which cannot be applied is
    named while the run has cost nothing; moving that gate after the cache
    build would mean a typo costs a full index.
    """
    import inspect

    from neoswga.cli import pipeline

    source = inspect.getsource(pipeline.run_step4)

    validate_at = source.index("design_request_for_run")
    cache_at = source.find("PositionCache(")

    if cache_at == -1:
        pytest.skip("run_step4 no longer builds the cache directly")
    assert validate_at < cache_at, (
        "the request is validated after the position cache is built, so a "
        "configuration mistake now costs a full index build"
    )
