"""The report shows the figures the acceptance record holds, not its own.

`PanelAssessment` carries every metric with the reach and denominator it was
computed at, and keeps a quantity that could not be measured apart from one
measured as zero. The report ignored it completely: it recomputed some figures
from the results CSV, read others from the summary JSON, and recomputed the
melting temperature under whatever conditions params.json happened to name.
`PanelAssessment.as_dict` meanwhile claimed "the report and the saved result
both read this", which was not true of the report.

The case that makes this more than tidiness is the selectivity sentinel.
`base_optimizer` reports MAX_SELECTIVITY (1e6) when the host carries no
exact-match site, because that is finite and JSON carries it, and its own
docstring concedes a reader cannot tell it from a measurement. The summary
therefore holds `selectivity_ratio: 1000000.0` beside `total_bg_sites: 0`, and
the report rendered the million. The assessment says undefined and why.
"""

import json
from pathlib import Path

import pytest

from neoswga.core.design_result import VALIDATION_FILENAME
from neoswga.core.report.metrics import _assessed, collect_pipeline_metrics

pytest.importorskip("pandas")


def _write_run(directory: Path, assessment: dict, summary_metrics: dict) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "step4_improved_df.csv").write_text(
        "primer,set_index,fg_count,bg_count\nACCACAGATAGC,0,3,0\n"
    )
    (directory / "step4_improved_df_summary.json").write_text(
        json.dumps({"metrics": summary_metrics})
    )
    (directory / VALIDATION_FILENAME).write_text(
        json.dumps({"ok": True, "issues": [], "assessment": assessment})
    )
    (directory / "params.json").write_text(
        json.dumps({"fg_seq_lengths": [5386], "bg_seq_lengths": [5386]})
    )
    return directory


def _measurement(value, unavailable=None):
    entry = {"value": value, "units": "dimensionless", "basis": "reach 3000 bp"}
    if unavailable:
        entry["unavailable"] = unavailable
    return entry


def _assessment(**metrics):
    return {
        "primers": ["ACCACAGATAGC"],
        "request_hash": "abc",
        "qualified": True,
        "acceptable": True,
        "violations": [],
        "blocking_violations": [],
        "advisory_violations": [],
        "metrics": metrics,
    }


def test_an_undefined_selectivity_is_not_rendered_as_a_million(tmp_path):
    """The sentinel must not reach a reader as though it were measured."""
    run = _write_run(
        tmp_path / "run",
        _assessment(
            fg_coverage=_measurement(1.0),
            total_bg_sites=_measurement(0.0),
            selectivity_ratio=_measurement(None, unavailable="zero denominator"),
        ),
        {"fg_coverage": 1.0, "total_bg_sites": 0, "selectivity_ratio": 1000000.0},
    )

    metrics = collect_pipeline_metrics(str(run))

    assert metrics.specificity.selectivity_ratio is None, (
        "the assessment says the ratio is undefined for want of a denominator; "
        "rendering the 1e6 sentinel states a specificity nobody measured"
    )
    assert metrics.specificity.background_sites == 0, "a measured zero is still a measurement"


def test_the_assessment_wins_over_the_summary(tmp_path):
    """Two saved figures for one panel: the acceptance record is the one shown."""
    run = _write_run(
        tmp_path / "run",
        _assessment(fg_coverage=_measurement(0.42), total_fg_sites=_measurement(7.0)),
        {"fg_coverage": 0.99, "total_fg_sites": 999},
    )

    metrics = collect_pipeline_metrics(str(run))

    assert metrics.coverage.overall_coverage == 0.42
    assert metrics.specificity.target_sites == 7
    assert metrics.coverage.from_optimizer is True


def test_a_directory_without_an_assessment_still_renders(tmp_path):
    """Directories written before the assessment existed are supported."""
    run = tmp_path / "old"
    run.mkdir()
    (run / "step4_improved_df.csv").write_text(
        "primer,set_index,fg_count,bg_count\nACCACAGATAGC,0,3,1\n"
    )
    (run / "step4_improved_df_summary.json").write_text(
        json.dumps({"metrics": {"fg_coverage": 0.77}})
    )
    (run / "params.json").write_text(json.dumps({"fg_seq_lengths": [5386]}))

    metrics = collect_pipeline_metrics(str(run))

    assert metrics.panel_assessment is None
    assert metrics.coverage.overall_coverage == 0.77, "the summary is still read"


def test_an_unmeasured_quantity_does_not_fall_back_to_an_estimate(tmp_path):
    """Unavailable must stay unavailable, not become the report's own guess.

    The report's own coverage estimate is `n_primers * 30000 / genome`, which
    for one primer on a 5.4 kb target saturates at 1.0. Substituting it for a
    quantity the assessment could not measure would render a perfect score for
    an absent measurement.
    """
    run = _write_run(
        tmp_path / "run",
        _assessment(fg_coverage=_measurement(None, unavailable="no position index")),
        {},
    )

    metrics = collect_pipeline_metrics(str(run))

    assert metrics.coverage.overall_coverage is None


def test_reading_one_measurement_separates_absent_from_unavailable():
    assessment = _assessment(
        measured=_measurement(0.5), missing=_measurement(None, unavailable="why")
    )
    assert _assessed(assessment, "measured") == (0.5, True)
    assert _assessed(assessment, "missing") == (None, True)
    assert _assessed(assessment, "never_recorded") == (None, False)
    assert _assessed(None, "measured") == (None, False)
