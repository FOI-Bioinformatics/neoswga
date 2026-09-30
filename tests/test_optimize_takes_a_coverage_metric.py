"""`optimize --coverage-metric` chooses the coverage a panel is judged on.

A configured design is judged on occupancy-weighted ("effective") coverage,
while the run prints the unweighted figure. `--target-coverage` is compared
against the first, so on the plasmid example at 100 bp reach a 12-primer panel
printed 46.7% against a 30% target and could not be reduced: the compared
figure was 19.7%. `plan-pool` has always let the user choose the metric;
`optimize` had no such option.

The default is unchanged. The flag is a None sentinel, so an absent flag
leaves the rule `objective_for_optimizer` already applied.
"""

from types import SimpleNamespace

import pytest

from neoswga.core.base_optimizer import OptimizerConfig
from neoswga.core.panel_refinement import objective_for_optimizer
from neoswga.core.pool_objective import PoolConstraints


def _optimizer(coverage_metric=None, conditions=object()):
    return SimpleNamespace(
        config=OptimizerConfig(coverage_metric=coverage_metric),
        conditions=conditions,
        compute_metrics=lambda primers: None,
        bg_prefixes=["host"],
        bg_seq_lengths=[1000],
    )


def test_a_configured_design_is_judged_on_effective_coverage_by_default():
    assert objective_for_optimizer(_optimizer()).constraints.coverage_metric == "effective"


def test_a_condition_free_evaluator_keeps_geometric_coverage():
    objective = objective_for_optimizer(_optimizer(conditions=None))

    assert objective.constraints.coverage_metric == "raw"


def test_the_configured_metric_is_the_one_the_objective_uses():
    objective = objective_for_optimizer(_optimizer("raw"))

    assert objective.constraints.coverage_metric == "raw"


def test_the_configured_metric_also_applies_when_limits_are_set():
    """Limits arrive as their own constraints, built with the default metric."""
    limits = PoolConstraints(max_background_sites=50)

    objective = objective_for_optimizer(_optimizer("raw"), limits)

    assert objective.constraints.coverage_metric == "raw"
    assert objective.constraints.max_background_sites == 50


def test_an_unknown_metric_is_refused_when_the_objective_is_built():
    """`PoolConstraints` holds the rule, so the config does not repeat it."""
    with pytest.raises(ValueError, match="coverage_metric"):
        objective_for_optimizer(_optimizer("measured"))


def test_the_keyword_reaches_the_optimizer_config():
    from neoswga.core.unified_optimizer import _build_optimizer_config

    config = _build_optimizer_config(6, False, 3000, False, {"coverage_metric": "raw"})

    assert config.coverage_metric == "raw"


def test_an_absent_keyword_leaves_the_metric_unset():
    from neoswga.core.unified_optimizer import _build_optimizer_config

    assert _build_optimizer_config(6, False, 3000, False, {}).coverage_metric is None


def _parse(argv):
    from neoswga.cli_unified import create_parser

    return create_parser().parse_args(argv)


def test_the_flag_is_a_sentinel_so_absence_changes_nothing():
    assert _parse(["optimize", "-j", "params.json"]).coverage_metric is None


def test_the_flag_travels_from_the_command_line_to_step_four():
    from neoswga.cli.pipeline import _step4_optimizer_kwargs

    args = _parse(["optimize", "-j", "params.json", "--coverage-metric", "raw"])

    assert _step4_optimizer_kwargs(args)["coverage_metric"] == "raw"


def test_params_json_can_set_it(tmp_path):
    """It resolves through `parameter`, like every other OptimizerConfig field."""
    from neoswga.core.search_control import resolve_search_settings

    assert resolve_search_settings({"coverage_metric": "raw"})["coverage_metric"] == "raw"
    assert resolve_search_settings({})["coverage_metric"] is None


def test_an_invalid_value_is_refused_at_load_time():
    """Before the position cache is built, so a typo costs nothing."""
    from neoswga.core.exceptions import InvalidDesignRequest
    from neoswga.core.search_control import resolve_search_settings

    with pytest.raises(InvalidDesignRequest, match="coverage_metric"):
        resolve_search_settings({"coverage_metric": "effektiv"})


def test_the_reach_metric_is_not_a_panel_metric():
    """`coverage.polymerase_extension_reach` takes a `coverage_metric` too, with
    values 'realistic' and 'processivity'. Same word, different question, so a
    value from one must not be accepted by the other."""
    from neoswga.core.exceptions import InvalidDesignRequest
    from neoswga.core.search_control import resolve_search_settings

    with pytest.raises(InvalidDesignRequest, match="coverage_metric"):
        resolve_search_settings({"coverage_metric": "realistic"})


def test_the_flag_beats_a_configured_value(monkeypatch):
    """A None argparse default is what makes this possible; see Known Issue 8."""
    from neoswga.core import parameter
    from neoswga.core.unified_optimizer import _build_optimizer_config

    monkeypatch.setattr(parameter, "coverage_metric", "effective", raising=False)

    assert (
        _build_optimizer_config(6, False, 3000, False, {"coverage_metric": "raw"}).coverage_metric
        == "raw"
    )


def test_a_configured_value_applies_when_no_flag_was_given(monkeypatch):
    from neoswga.core import parameter
    from neoswga.core.unified_optimizer import _build_optimizer_config

    monkeypatch.setattr(parameter, "coverage_metric", "raw", raising=False)

    assert _build_optimizer_config(6, False, 3000, False, {}).coverage_metric == "raw"


def test_the_target_help_names_the_coverage_it_is_compared_against(capsys):
    with pytest.raises(SystemExit):
        _parse(["optimize", "--help"])

    text = " ".join(capsys.readouterr().out.split())
    assert "--coverage-metric" in text
    assert "occupancy-weighted" in text
