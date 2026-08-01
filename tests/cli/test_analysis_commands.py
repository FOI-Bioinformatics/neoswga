"""Execution tests for the analysis and reporting commands.

`analyze-coverage`, `predict-efficiency` and `report` are read-only commands
that summarise a design. Being read-only makes them easy to get wrong quietly:
nothing downstream fails when a number is meaningless, so a zero or an empty
section looks like a finding rather than a fault. That is exactly what happened
with coverage before the audit.
"""

import json
import os

import pytest

pytestmark = pytest.mark.usefixtures("pipeline_run")


@pytest.fixture
def primers(scored_primers):
    return scored_primers[:6]


# ----------------------------------------------------------------------
# analyze-coverage
# ----------------------------------------------------------------------


def test_analyze_coverage_writes_bed_and_json(run_cli, pipeline_run, primers, tmp_path):
    out = tmp_path / "cov"
    run_cli(
        [
            "analyze-coverage",
            "-j",
            pipeline_run["params_file"],
            "--primers",
            *primers,
            "-o",
            str(out),
        ]
    )
    written = set(os.listdir(out))
    assert any(f.endswith(".bed") for f in written), written
    assert any(f.endswith(".json") for f in written), written


def test_analyze_coverage_json_is_valid_and_populated(run_cli, pipeline_run, primers, tmp_path):
    out = tmp_path / "cov"
    run_cli(
        [
            "analyze-coverage",
            "-j",
            pipeline_run["params_file"],
            "--primers",
            *primers,
            "-o",
            str(out),
        ]
    )

    def reject(token):
        raise ValueError(f"non-standard JSON token {token!r}")

    files = [f for f in os.listdir(out) if f.endswith(".json")]
    payloads = [json.loads((out / f).read_text(), parse_constant=reject) for f in files]
    assert any(payload for payload in payloads), "every coverage JSON was empty"


def test_analyze_coverage_is_read_only(run_cli, pipeline_run, primers, tmp_path):
    """It must not mutate the workspace it is inspecting.

    An analysis command that rewrites the user's data_dir would silently
    invalidate a design in progress.
    """
    data_dir = pipeline_run["data_dir"]
    before = {f: os.path.getmtime(os.path.join(data_dir, f)) for f in os.listdir(data_dir)}

    run_cli(
        [
            "analyze-coverage",
            "-j",
            pipeline_run["params_file"],
            "--primers",
            *primers,
            "-o",
            str(tmp_path / "cov"),
        ]
    )

    after = {f: os.path.getmtime(os.path.join(data_dir, f)) for f in os.listdir(data_dir)}
    assert before == after, "analyze-coverage modified the pipeline data directory"


def test_min_gap_size_filters_reported_gaps(run_cli, pipeline_run, primers, tmp_path):
    """A larger minimum must not report more gaps than a smaller one."""

    def gap_count(min_gap, tag):
        out = tmp_path / tag
        run_cli(
            [
                "analyze-coverage",
                "-j",
                pipeline_run["params_file"],
                "--primers",
                *primers,
                "--min-gap-size",
                str(min_gap),
                "-o",
                str(out),
            ]
        )
        beds = [f for f in os.listdir(out) if f.endswith(".bed")]
        return sum(
            len([ln for ln in (out / b).read_text().splitlines() if ln.strip()]) for b in beds
        )

    assert gap_count(10_000, "big") <= gap_count(100, "small")


# ----------------------------------------------------------------------
# predict-efficiency
# ----------------------------------------------------------------------


def _predict(run_cli, pipeline_run, primers, caplog):
    """Run predict-efficiency and return everything it reported.

    It reports through `logger.info` rather than print, so stdout is empty --
    which is only observable now that the suite no longer disables logging
    globally at collection time.
    """
    import logging

    with caplog.at_level(logging.INFO):
        run_cli(["predict-efficiency", "-j", pipeline_run["params_file"], "--primers", *primers])
    return caplog.text


def test_predict_efficiency_runs_and_scores(run_cli, pipeline_run, primers, caplog):
    reported = _predict(run_cli, pipeline_run, primers, caplog)
    assert reported.strip(), "predict-efficiency reported nothing"
    assert "efficiency" in reported.lower()


def test_predict_efficiency_network_score_is_not_spuriously_zero(
    run_cli, pipeline_run, primers, caplog
):
    """The documented failure: a real ML score alongside a network score of zero.

    That zero came from primers missing in the position cache, not from a bad
    design, and nothing distinguished the two. With primers drawn from the
    pipeline's own scored table the positions exist, so a zero here means the
    lookup is broken again.
    """
    reported = _predict(run_cli, pipeline_run, primers, caplog).lower()

    # Any reported coverage/network figure should be non-zero for a set that
    # demonstrably binds the target.
    zero_markers = ("coverage: 0.0%", "coverage: 0%", "network score: 0.0")
    assert not any(marker in reported for marker in zero_markers), reported[-800:]


# ----------------------------------------------------------------------
# report
# ----------------------------------------------------------------------


def test_report_check_mode_validates_without_writing(run_cli, pipeline_run, tmp_path):
    """--check should report readiness, not produce a report."""
    before = set(os.listdir(pipeline_run["data_dir"]))
    try:
        run_cli(["report", "-d", pipeline_run["data_dir"], "--check"])
    except SystemExit as exc:
        # Exits non-zero when the run is incomplete, which it is here: the
        # fixture stops after `score`, so step4 output is absent. That is a
        # legitimate answer from a validation mode.
        assert exc.code in (0, 1)
    after = set(os.listdir(pipeline_run["data_dir"]))
    assert before == after, "--check wrote files into the data directory"


def test_report_on_an_incomplete_run_explains_itself(run_cli, pipeline_run, capsys, caplog):
    """A missing optimize step should produce a clear reason, not a traceback."""
    import logging

    with caplog.at_level(logging.ERROR):
        try:
            run_cli(["report", "-d", pipeline_run["data_dir"], "--check"])
        except SystemExit:
            pass

    combined = (capsys.readouterr().out + capsys.readouterr().err + caplog.text).lower()
    assert combined.strip(), "report --check said nothing about an incomplete run"
