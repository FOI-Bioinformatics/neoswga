"""`filter` has sixteen chemistry flags and a --preset; `optimize` has none.

So `neoswga filter -j params.json --preset high_gc_genome` followed by a plain
`neoswga optimize -j params.json` filters under one reaction and scores under
another. The optimizer builds its conditions from the parameter globals, which
come back from params.json, and nothing said the two disagreed.

This is a flag asymmetry, not an adaptation problem: the GC-adaptive strategy
runs in both steps, because unified_optimizer calls _initialize() before
building conditions.
"""

import json
import logging
import random
import shutil
import subprocess
import sys

import pytest

from neoswga.cli._common import warn_on_condition_drift
from neoswga.core import run_manifest as rm


def _drift_warnings(caplog):
    """Warning records the drift check itself emitted.

    Not ``caplog.text == ""``. caplog collects the whole call phase, so the
    INFO line write_manifest logs while the case is being set up lands in it
    too, and whether that line is captured at all depends on the root log
    level -- which other tests in the same xdist worker raise and lower. That
    made the silent-case assertions fail on worker assignment rather than on
    behaviour.
    """
    return [
        r
        for r in caplog.records
        if r.levelno >= logging.WARNING and "Reaction conditions differ" in r.getMessage()
    ]


class _Parameter:
    def __init__(self, data_dir, **conditions):
        self.data_dir = data_dir
        for name, value in conditions.items():
            setattr(self, name, value)


def test_drift_is_reported_field_by_field(tmp_path, caplog):
    rm.write_manifest(
        step="filter",
        data_dir=str(tmp_path),
        effective_conditions={
            "polymerase": "phi29",
            "reaction_temp": 30.0,
            "betaine_m": 2.0,
        },
    )

    parameter = _Parameter(str(tmp_path), polymerase="phi29", reaction_temp=30.0, betaine_m=0.0)

    with caplog.at_level(logging.WARNING):
        differing = warn_on_condition_drift(parameter)

    assert differing == ["betaine_m"]
    assert "betaine_m" in caplog.text
    assert "filter" in caplog.text


def test_the_warning_names_both_values_not_just_the_field(tmp_path, caplog):
    """A message naming only the field leaves the user to go and find both
    numbers themselves, and every test above still passes."""
    rm.write_manifest(
        step="filter",
        data_dir=str(tmp_path),
        effective_conditions={"betaine_m": 2.0},
    )
    parameter = _Parameter(str(tmp_path), betaine_m=0.0)

    with caplog.at_level(logging.WARNING):
        warn_on_condition_drift(parameter)

    assert "2.0" in caplog.text, "the filter step's value is missing"
    assert "0.0" in caplog.text, "the value this step would use is missing"


def test_agreeing_conditions_are_silent(tmp_path, caplog):
    rm.write_manifest(
        step="filter",
        data_dir=str(tmp_path),
        effective_conditions={"polymerase": "phi29", "reaction_temp": 30.0},
    )
    parameter = _Parameter(str(tmp_path), polymerase="phi29", reaction_temp=30.0)

    with caplog.at_level(logging.WARNING):
        differing = warn_on_condition_drift(parameter)

    assert differing == []
    assert _drift_warnings(caplog) == []


def test_no_filter_entry_is_silent(tmp_path, caplog):
    """A user optimizing a directory produced elsewhere has nothing to compare."""
    parameter = _Parameter(str(tmp_path), polymerase="phi29")

    with caplog.at_level(logging.WARNING):
        differing = warn_on_condition_drift(parameter)

    assert differing == []
    assert _drift_warnings(caplog) == []


def test_a_field_recorded_by_filter_and_absent_now_counts_as_drift(tmp_path, caplog):
    rm.write_manifest(
        step="filter",
        data_dir=str(tmp_path),
        effective_conditions={"polymerase": "bst", "dmso_percent": 5.0},
    )
    parameter = _Parameter(str(tmp_path), polymerase="bst")

    with caplog.at_level(logging.WARNING):
        differing = warn_on_condition_drift(parameter)

    assert differing == ["dmso_percent"]


def test_a_real_optimize_run_warns_after_a_preset_filter(tmp_path):
    """Drive the actual CLI, because the check's placement is the hard part.

    The first version of this test was a grep of cli/pipeline.py for the call.
    It passed while the call sat where it could never fire: `get_params` runs
    inside `optimize_step4`, so at every point before that call `data_dir` is
    None and every reaction global still holds its module default. On a real
    run with three differing fields the check returned an empty list and said
    nothing. Only running the steps catches that.

    Separate processes are required. Run in one process the parameter globals
    persist across steps, so a preset applied at `filter` is still in effect at
    `optimize` and there is genuinely no drift to report.
    """
    if not shutil.which("jellyfish"):
        pytest.skip("jellyfish not available (required for count-kmers)")

    rng = random.Random(20260910)
    seq = "".join(rng.choice("ACGT") for _ in range(20_000))
    fasta = tmp_path / "target.fasta"
    fasta.write_text(
        ">target\n" + "\n".join(seq[i : i + 70] for i in range(0, len(seq), 70)) + "\n"
    )

    params_file = tmp_path / "params.json"
    params_file.write_text(
        json.dumps(
            {
                "fg_genomes": [str(fasta)],
                "bg_genomes": [],
                "fg_prefixes": [str(tmp_path / "target")],
                "bg_prefixes": [],
                "data_dir": str(tmp_path / "results"),
                "min_k": 10,
                "max_k": 10,
                "polymerase": "phi29",
                "reaction_temp": 30.0,
                "min_fg_freq": 1e-6,
                "max_bg_freq": 1.0,
                "max_gini": 1.0,
                "max_primer": 60,
                "min_tm": 0,
                "max_tm": 100,
                "gc_min": 0.0,
                "gc_max": 1.0,
                "num_primers": 4,
                "target_set_size": 4,
                "max_sets": 2,
                "iterations": 2,
                "cpus": 1,
                "fg_circular": True,
                "schema_version": 2,
            }
        )
    )

    def _step(*extra, expect_failure=False):
        argv = [sys.executable, "-m", "neoswga.cli_unified", *extra, "-j", str(params_file)]
        proc = subprocess.run(argv, capture_output=True, text=True, cwd=str(tmp_path), timeout=900)
        if expect_failure:
            assert proc.returncode != 0, (
                f"{extra[0]} succeeded on a pool selected under another reaction:\n"
                f"{proc.stdout[-2000:]}"
            )
            return proc
        if proc.returncode != 0:
            # A step that RAN and failed is a failure, not a skip. Skipping
            # here turned a broken pipeline green: on 2026-09-10 two full-suite
            # runs reported one more skip and one fewer pass than the runs
            # either side, on a byte-identical tree, and these tests pass in
            # isolation. The genuinely absent prerequisite, jellyfish, is
            # checked before any of this runs.
            raise AssertionError(f"{extra[0]} exited {proc.returncode}:\n{proc.stderr[-2000:]}")
        return proc

    _step("count-kmers")
    _step("filter", "--preset", "enhanced_equiphi29")
    _step("prepare-candidates")
    optimize = _step("optimize", "--seed", "1", expect_failure=True)

    output = optimize.stdout + optimize.stderr
    # Strengthened on 2026-09-21. This used to assert a warning and a
    # successful run. A warning was as far as the old code could go, because
    # the inventory-open error quietly switched the search to `step3_df.csv`
    # -- so the run did proceed, over a shortlist selected at 42 C under
    # equiphi29 while scoring it at 30 C under phi29, and every candidate the
    # inventory held was unreachable.
    #
    # Under the valid-design contract that is a refusal with a named remedy.
    # The drift is the finding, not a footnote to a result.
    assert "candidate inventory" in output, output[-2000:]
    assert "neoswga filter" in output, output[-2000:]
    assert "polymerase" in output
