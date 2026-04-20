"""Unit tests for OptimizationResult.validate() — the post-optimization
sanity checker that catches silent failures an optimizer can return
(duplicate primers, wrong set size, zero coverage, blacklist re-injection).
"""

import pytest

from neoswga.core.base_optimizer import (
    OptimizationResult, OptimizationStatus, PrimerSetMetrics,
)


def _make_result(primers, status=OptimizationStatus.SUCCESS, **metric_overrides):
    metrics = PrimerSetMetrics.empty()
    for k, v in metric_overrides.items():
        object.__setattr__(metrics, k, v) if hasattr(metrics, "__dataclass_fields__") else setattr(metrics, k, v)
    return OptimizationResult(
        primers=tuple(primers),
        score=1.0,
        status=status,
        metrics=metrics,
        iterations=1,
        optimizer_name="test",
    )


def test_validate_clean_result_reports_no_issues():
    r = _make_result(["ATCG", "GCTA", "TACG"], fg_coverage=0.9)
    report = r.validate(target_size=3, min_coverage=0.5)
    assert report["ok"] is True
    assert report["issues"] == []
    assert report["num_primers"] == 3


def test_validate_flags_duplicates_as_error():
    r = _make_result(["ATCG", "GCTA", "ATCG"], fg_coverage=0.9)
    report = r.validate(target_size=3)
    assert report["ok"] is False
    codes = [i["code"] for i in report["issues"]]
    assert "duplicate_primers" in codes
    err = [i for i in report["issues"] if i["code"] == "duplicate_primers"][0]
    assert err["level"] == "error"


def test_validate_flags_size_mismatch_as_error_on_success():
    r = _make_result(["ATCG", "GCTA"], fg_coverage=0.9)
    report = r.validate(target_size=4)
    codes = [i["code"] for i in report["issues"]]
    assert "set_size_mismatch" in codes
    issue = [i for i in report["issues"] if i["code"] == "set_size_mismatch"][0]
    assert issue["level"] == "error"


def test_validate_downgrades_size_mismatch_to_warning_on_partial():
    r = _make_result(["ATCG", "GCTA"], status=OptimizationStatus.PARTIAL,
                     fg_coverage=0.5)
    report = r.validate(target_size=4)
    issues = [i for i in report["issues"] if i["code"] == "set_size_mismatch"]
    assert issues and issues[0]["level"] == "warning"


def test_validate_flags_coverage_below_threshold():
    r = _make_result(["ATCG", "GCTA", "TACG"], fg_coverage=0.05)
    report = r.validate(target_size=3, min_coverage=0.3)
    codes = [i["code"] for i in report["issues"]]
    assert "coverage_below_threshold" in codes


def test_validate_flags_blacklist_reinjection():
    r = _make_result(["ATCG", "BLKL", "TACG"], fg_coverage=0.9)
    report = r.validate(target_size=3, forbidden_primers=["BLKL", "EVIL"])
    issues = [i for i in report["issues"] if i["code"] == "blacklist_primer_in_set"]
    assert issues and issues[0]["level"] == "error"


def test_validate_skips_coverage_check_on_error_status():
    r = _make_result([], status=OptimizationStatus.ERROR)
    report = r.validate(target_size=3, min_coverage=0.5)
    codes = [i["code"] for i in report["issues"]]
    assert "coverage_below_threshold" not in codes


def test_per_target_coverage_below_threshold_warning(monkeypatch):
    r = _make_result(["ATCG", "GCTA"], fg_coverage=0.9)
    # Inject per_target_coverage on metrics
    object.__setattr__(r.metrics, "per_target_coverage", {"target_a": 0.95, "target_b": 0.35})
    report = r.validate(
        target_size=2, min_coverage=0.5, min_per_target_coverage=0.50,
    )
    codes = [i["code"] for i in report["issues"]]
    assert "per_target_coverage_below_threshold" in codes
    issue = [i for i in report["issues"] if i["code"] == "per_target_coverage_below_threshold"][0]
    assert issue["level"] == "warning"
    assert "target_b" in issue["detail"]
