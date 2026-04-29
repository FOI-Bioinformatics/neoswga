"""Tests for surfacing extension_reach in coverage labelling.

The optimizer measures interval coverage at a specific per-primer reach
(default 3000 bp post-Phase 16). The report should record this value so
users can interpret the reported coverage instead of guessing whether
"95% coverage" means 95% bin-coverage at 10 kb, 95% interval coverage at
3 kb reach, or 95% theoretical reachability at 70 kb.
"""

import json
from pathlib import Path

import pytest

from neoswga.core.base_optimizer import PrimerSetMetrics


def _make_metrics(**overrides):
    defaults = dict(
        fg_coverage=0.95,
        bg_coverage=0.0,
        coverage_uniformity=0.3,
        total_fg_sites=100,
        total_bg_sites=0,
        selectivity_ratio=999.0,
        mean_tm=50.0,
        tm_range=(48.0, 52.0),
        dimer_risk_score=0.05,
        mean_gap=5000.0,
        max_gap=15000.0,
        gap_gini=0.3,
        gap_entropy=8.0,
        strand_alternation_score=0.5,
        strand_coverage_ratio=0.5,
    )
    defaults.update(overrides)
    return PrimerSetMetrics(**defaults)


def test_primer_set_metrics_to_dict_includes_extension_reach():
    metrics = _make_metrics(extension_reach=3000)
    assert metrics.to_dict()["extension_reach"] == 3000


def test_primer_set_metrics_default_extension_reach_matches_phase16():
    """Backward-compat: existing callers that do not pass extension_reach
    should get the Phase 16 default (3000 bp realistic per-primer reach).
    """
    metrics = _make_metrics()
    assert metrics.extension_reach == 3000


def test_report_metrics_loader_propagates_extension_reach(tmp_path):
    """When step4_improved_df_summary.json carries extension_reach, the
    report's CoverageMetrics must surface it.
    """
    from neoswga.core.report import metrics as report_metrics

    results = tmp_path
    summary = {
        "metrics": {
            "fg_coverage": 0.95,
            "extension_reach": 4000,
            "mean_gap": 1000.0,
            "max_gap": 5000.0,
            "gap_gini": 0.3,
            "gap_entropy": 7.0,
        }
    }
    (results / "step4_improved_df_summary.json").write_text(json.dumps(summary))

    loaded = report_metrics._load_optimizer_summary(results)
    assert loaded["metrics"]["extension_reach"] == 4000


def test_report_coverage_metrics_dataclass_has_extension_reach_field():
    """CoverageMetrics must carry extension_reach so report renderers can
    label the coverage value with its source reach.
    """
    from neoswga.core.report.metrics import CoverageMetrics

    cm = CoverageMetrics()
    assert hasattr(cm, "extension_reach")
    cm.extension_reach = 3000
    assert cm.extension_reach == 3000
