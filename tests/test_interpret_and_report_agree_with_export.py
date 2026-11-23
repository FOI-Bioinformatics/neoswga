"""`interpret` and `report` refuse what `export` refuses.

`export_is_blocked` reads `design_failure.json`, the validator findings and the
set index. `interpret` read only the validator findings and `report` read
nothing, so after a failed run `interpret` printed "Ready for synthesis" and
`report` printed "PROCEED WITH MONITORING" while `export` refused; and
`interpret --set 1` judged set 1 by set 0's findings
(2026-11-23 sweep, tier 1 item 3).
"""

import csv
import json

from neoswga.core.report.metrics import PipelineMetrics, _load_validation_report
from neoswga.core.report.quality import calculate_quality_grade
from neoswga.core.results_interpreter import ResultsInterpreter


def _finished_run(path, n=6):
    rows = [
        {
            "primer": f"ATCGATCG{i:02d}"[-8:],
            "fg_freq": "1e-3",
            "bg_freq": "1e-5",
            "gini": "0.35",
            "dimer_score": "0.1",
            "coverage": "0.92",
            "set_index": str(s),
        }
        for s in (0, 1)
        for i in range(n)
    ]
    with open(path / "step4_improved_df.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return path


def _fail(path):
    (path / "design_failure.json").write_text(
        json.dumps({"run_state": "failed", "stage": "optimize", "message": "refused"})
    )


def test_interpret_refuses_after_a_failed_run(tmp_path):
    _fail(_finished_run(tmp_path))
    report = ResultsInterpreter(str(tmp_path)).analyze()
    assert report.recommendation.startswith("Not ready for synthesis")
    assert "did not finish" in report.recommendation


def test_interpret_set_1_is_not_ready(tmp_path):
    _finished_run(tmp_path)
    report = ResultsInterpreter(str(tmp_path), set_index=1).analyze()
    assert report.recommendation.startswith("Not ready for synthesis")
    assert "has not been assessed" in report.recommendation


def test_interpret_set_0_of_a_clean_run_is_not_refused(tmp_path):
    _finished_run(tmp_path)
    report = ResultsInterpreter(str(tmp_path)).analyze()
    assert not report.recommendation.startswith("Not ready for synthesis")


def test_report_recommendation_refuses_a_blocked_run():
    metrics = PipelineMetrics(primer_count=6)
    clean = calculate_quality_grade(metrics)
    metrics.blocked_reason = "The last run in this directory did not finish."
    blocked = calculate_quality_grade(metrics)
    assert blocked.recommendation == "DO NOT ORDER"
    assert metrics.blocked_reason in blocked.recommendation_details
    assert blocked.grade == clean.grade


def test_unreadable_validation_report_is_not_ok(tmp_path):
    (tmp_path / "step4_improved_df_validation.json").write_text("{ truncated")
    issues, ok = _load_validation_report(tmp_path)
    assert ok is False
    assert issues and issues[0]["level"] == "error"
