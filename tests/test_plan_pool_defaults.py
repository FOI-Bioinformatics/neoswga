"""`plan-pool` takes length, coverage metric and method from params.json.

`--primer-length` defaulted to 12, `--coverage-metric` to "effective", and
`--method` fell through to background-aware or hybrid, so the params file's
k range, `coverage_metric` and `optimization_method` were ignored. A file with
`min_k = max_k = 10` failed with "No 12-mers" (2026-11-23 sweep, tier 2; the
Known Issue 8 shape).
"""

import argparse

import pytest

from neoswga.cli.plan_pool import _resolve_plan_options


def _args(**flags):
    base = {"primer_length": None, "coverage_metric": None, "method": None}
    base.update(flags)
    return argparse.Namespace(**base)


def test_a_single_length_comes_from_the_k_range():
    length, _, _ = _resolve_plan_options(_args(), {"min_k": 10, "max_k": 10}, has_background=True)
    assert length == 10


def test_a_k_range_without_the_flag_is_refused_by_name():
    with pytest.raises(ValueError, match="--primer-length"):
        _resolve_plan_options(_args(), {"min_k": 8, "max_k": 12}, has_background=True)


def test_the_flag_wins_over_the_file():
    length, metric, method = _resolve_plan_options(
        _args(primer_length=11, coverage_metric="raw", method="clique"),
        {"min_k": 10, "max_k": 10, "coverage_metric": "effective", "optimization_method": "hybrid"},
        has_background=True,
    )
    assert (length, metric, method) == (11, "raw", "clique")


def test_the_file_decides_metric_and_method():
    _, metric, method = _resolve_plan_options(
        _args(),
        {"min_k": 12, "max_k": 12, "coverage_metric": "raw", "optimization_method": "dominating-set"},
        has_background=True,
    )
    assert (metric, method) == ("raw", "dominating-set")


@pytest.mark.parametrize("background,expected", [(True, "background-aware"), (False, "hybrid")])
def test_without_either_the_old_defaults_hold(background, expected):
    _, metric, method = _resolve_plan_options(
        _args(), {"min_k": 12, "max_k": 12}, has_background=background
    )
    assert (metric, method) == ("effective", expected)


def test_a_method_plan_pool_cannot_run_falls_back():
    _, _, method = _resolve_plan_options(
        _args(), {"min_k": 12, "max_k": 12, "optimization_method": "ensemble"}, has_background=False
    )
    assert method == "hybrid"
