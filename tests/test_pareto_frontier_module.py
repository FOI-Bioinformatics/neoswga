"""Coverage for neoswga.core.pareto_frontier CLI/report helpers, driven by a
synthetic FrontierResult (no optimizer run needed)."""

import json

import pytest

from neoswga.core.pareto_frontier import (
    generate_frontier_json,
    generate_frontier_report,
    summarize_frontier_for_cli,
)
from neoswga.core.set_size_optimizer import FrontierResult, SetSizeMetrics


def _point(size, fg, bg, fg_sites, bg_sites):
    ratio = fg_sites / bg_sites if bg_sites else float(fg_sites)
    return SetSizeMetrics(
        set_size=size,
        fg_coverage=fg,
        bg_coverage=bg,
        fg_binding_sites=fg_sites,
        bg_binding_sites=bg_sites,
        fg_bg_ratio=ratio,
        primers=tuple(f"PRIMER{i}" for i in range(size)),
    )


@pytest.fixture
def frontier():
    pts = [
        _point(4, 0.55, 0.02, 400, 8),
        _point(6, 0.75, 0.03, 620, 12),
        _point(8, 0.88, 0.05, 810, 20),
    ]
    return FrontierResult(
        all_points=pts,
        pareto_points=pts,
        selected_point=pts[1],
        selection_explanation="balanced pick",
    )


def test_generate_frontier_json_is_valid_and_complete(frontier):
    out = generate_frontier_json(frontier, application="enrichment")
    data = json.loads(out) if isinstance(out, str) else out
    assert isinstance(data, dict)
    # the JSON must expose the pareto points and the selection
    blob = json.dumps(data)
    assert "selected" in blob or "selection" in blob
    assert "6" in blob  # the selected set size appears somewhere


def test_summarize_frontier_for_cli_returns_text(frontier):
    summary = summarize_frontier_for_cli(frontier, application="enrichment")
    assert isinstance(summary, str) and summary
    # the three set sizes should be referenced
    for size in ("4", "6", "8"):
        assert size in summary


def test_generate_frontier_report_text(frontier):
    report = generate_frontier_report(frontier, application="enrichment")
    assert isinstance(report, str) and report
    assert "Pareto Frontier" in report
    # the set sizes on the frontier are listed
    for size in ("4", "6", "8"):
        assert size in report


def test_helpers_handle_single_point_frontier():
    pts = [_point(6, 0.7, 0.03, 600, 15)]
    fr = FrontierResult(all_points=pts, pareto_points=pts, selected_point=pts[0])
    # Should not raise on a degenerate single-point frontier.
    assert summarize_frontier_for_cli(fr)
    assert generate_frontier_json(fr)
