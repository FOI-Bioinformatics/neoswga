"""Phase 17A — validator warnings (per_target_coverage_below_threshold,
blacklist_primer_in_set, duplicate_primers, ...) must surface in the
HTML report. Previously they lived only in step4_improved_df_validation.json
and a user looking at the rendered report could see Grade A with a
hidden warning. These tests lock in the visible rendering path.
"""

import json
from datetime import datetime
from pathlib import Path

import pytest

from neoswga.core.report.executive_summary import (
    create_executive_summary,
    render_executive_summary,
)
from neoswga.core.report.metrics import (
    PipelineMetrics,
    _load_validation_report,
    collect_pipeline_metrics,
)
from neoswga.core.report.quality import calculate_quality_grade
from neoswga.core.report.technical_report import (
    TechnicalReportData,
    render_technical_report,
)


def _metrics_with(issues):
    m = PipelineMetrics()
    m.validation_issues = list(issues)
    m.validation_ok = not any(i.get("level") == "error" for i in issues)
    return m


def test_load_validation_report_missing_file_yields_empty(tmp_path):
    """Missing validation.json is NOT an error — absence means 'no issues'."""
    issues, ok = _load_validation_report(tmp_path)
    assert issues == []
    assert ok is True


def test_load_validation_report_parses_issues(tmp_path):
    """The loader normalises the dict produced by
    OptimizationResult.validate()."""
    payload = {
        "optimizer": "hybrid",
        "num_primers": 6,
        "ok": False,
        "issues": [
            {
                "level": "warning",
                "code": "per_target_coverage_below_threshold",
                "detail": "1 target(s) below 0.50: {'pcDNA': 0.35}",
            },
            {
                "level": "error",
                "code": "blacklist_primer_in_set",
                "detail": "Primers in blacklist: ['AAAAAAAA']",
            },
        ],
    }
    (tmp_path / "step4_improved_df_validation.json").write_text(
        json.dumps(payload)
    )
    issues, ok = _load_validation_report(tmp_path)
    assert ok is False
    assert len(issues) == 2
    assert issues[0]["code"] == "per_target_coverage_below_threshold"
    assert issues[0]["level"] == "warning"
    assert issues[1]["code"] == "blacklist_primer_in_set"
    assert issues[1]["level"] == "error"


def test_load_validation_report_handles_malformed(tmp_path):
    """Malformed JSON must not crash report generation."""
    (tmp_path / "step4_improved_df_validation.json").write_text("{not json")
    issues, ok = _load_validation_report(tmp_path)
    assert issues == []
    assert ok is True  # best-effort: treat unreadable as no-issues


def test_executive_summary_without_issues_has_no_banner():
    """If the validator produced no issues, the banner div must be absent
    so Grade A results are not visually cluttered."""
    m = PipelineMetrics()  # validation_issues default = []
    q = calculate_quality_grade(m)
    html = render_executive_summary(create_executive_summary(m, q))
    assert '<div class="validation-banner' not in html


def test_executive_summary_renders_per_target_warning():
    """per_target_coverage_below_threshold must be user-visible: the
    issue code, the human label, and the warning CSS class."""
    issues = [{
        "level": "warning",
        "code": "per_target_coverage_below_threshold",
        "detail": "1 target(s) below 0.50: {'pcDNA': 0.35}",
    }]
    m = _metrics_with(issues)
    q = calculate_quality_grade(m)
    html = render_executive_summary(create_executive_summary(m, q))
    assert 'validation-banner level-warning' in html
    assert 'per_target_coverage_below_threshold' in html
    assert 'One or more targets below coverage threshold' in html
    assert 'pcDNA' in html
    assert 'Validation warnings' in html


def test_executive_summary_renders_error_banner():
    """An error-level issue (e.g. duplicate primers) must get the
    error-styled banner and a stronger heading, not a warning banner."""
    issues = [{
        "level": "error",
        "code": "duplicate_primers",
        "detail": "1 duplicate primer(s): ['ACGTACGT']",
    }]
    m = _metrics_with(issues)
    q = calculate_quality_grade(m)
    html = render_executive_summary(create_executive_summary(m, q))
    assert 'validation-banner level-error' in html
    assert 'Validation errors' in html
    assert 'Duplicate primers in set' in html


def test_executive_summary_mixed_issues_use_error_class():
    """When both warnings and errors are present, the banner adopts the
    error style (the higher-severity signal wins)."""
    issues = [
        {"level": "warning", "code": "set_size_mismatch", "detail": "5 vs 6"},
        {"level": "error", "code": "blacklist_primer_in_set",
         "detail": "blacklist AAAA"},
    ]
    m = _metrics_with(issues)
    q = calculate_quality_grade(m)
    html = render_executive_summary(create_executive_summary(m, q))
    assert 'level-error' in html
    assert 'level-warning' not in html.split('validation-banner')[1].split(
        '</div>'
    )[0]


def test_technical_report_also_surfaces_issues():
    """Both entry points (executive + technical) must render the banner;
    users pick one based on interactive/non-interactive mode."""
    issues = [{
        "level": "warning",
        "code": "per_target_coverage_below_threshold",
        "detail": "target_b: 0.40",
    }]
    m = _metrics_with(issues)
    q = calculate_quality_grade(m)
    data = TechnicalReportData(
        generated_at=datetime.now().isoformat(),
        metrics=m,
        quality=q,
    )
    html = render_technical_report(data)
    assert 'validation-banner level-warning' in html
    assert 'per_target_coverage_below_threshold' in html
    assert 'target_b' in html


def test_collect_pipeline_metrics_loads_validation(tmp_path):
    """End-to-end: dropping a validation JSON in results_dir makes the
    warnings show up on metrics.validation_issues."""
    # Minimal results_dir with just the validation file.
    payload = {
        "optimizer": "hybrid",
        "num_primers": 6,
        "ok": False,
        "issues": [{
            "level": "warning",
            "code": "coverage_below_threshold",
            "detail": "fg_coverage=0.40 < min_coverage=0.50",
        }],
    }
    (tmp_path / "step4_improved_df_validation.json").write_text(
        json.dumps(payload)
    )
    # Minimal step4_df so collect_pipeline_metrics loads primers (but our
    # banner test only needs validation_issues populated).
    (tmp_path / "step4_improved_df.csv").write_text(
        "primer,fg_count,bg_count,gini\nACGT,10,0,0.2\n"
    )
    metrics = collect_pipeline_metrics(str(tmp_path))
    assert len(metrics.validation_issues) == 1
    assert metrics.validation_issues[0]["code"] == "coverage_below_threshold"
    assert metrics.validation_ok is False
