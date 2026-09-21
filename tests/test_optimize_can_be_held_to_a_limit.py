"""`optimize` should be able to enforce a limit, not only `plan-pool`.

The audit found the constraint framework reaching one command:
`PoolConstraints` is enforced inside the search through `shortfall` and
repaired by a beam, and only `plan-pool` ever constructed one
(`docs/validation/pool_selection_audit_2026-09-18.md`). swga 1.0 has had two
hard constraints enforced inside its search since 2017.

Item 2 of the ranked list. Every limit is unset by default, so a run that asks
for nothing gets exactly the panel it got before: a new limit that changed a
delivered panel unasked would be the scoring change both wet-lab benchmarks
refuse.
"""

import json
import pathlib

import pytest

from neoswga.core.panel_acceptance import (
    LIMIT_KEYS,
    constraints_from_parameter,
    enforce_constraints,
)
from neoswga.core.pool_objective import PoolConstraints

SCHEMA = json.loads(
    (pathlib.Path("neoswga") / "core" / "schema" / "params.schema.json").read_text()
)


class _Params:
    """A stand-in for the `parameter` module, carrying only what is set."""

    def __init__(self, **limits):
        for key, value in limits.items():
            setattr(self, key, value)


class TestTheKeysAreDeclaredAndBound:
    @pytest.mark.parametrize("key", sorted(LIMIT_KEYS))
    def test_the_key_is_declared_in_the_schema(self, key):
        """Otherwise setting it produces the unknown-key warning and nothing."""
        assert key in SCHEMA["properties"]

    @pytest.mark.parametrize("key", sorted(LIMIT_KEYS))
    def test_the_key_binds_a_parameter_global(self, key):
        """`test_no_schema_key_is_inert` is the general ratchet; this names the
        six so a later edit cannot quietly drop one."""
        from neoswga.core import parameter

        assert hasattr(parameter, key)

    @pytest.mark.parametrize("key", sorted(LIMIT_KEYS))
    def test_the_key_defaults_to_unset(self, key):
        from neoswga.core import parameter

        assert getattr(parameter, key) is None

    @pytest.mark.parametrize("key", sorted(LIMIT_KEYS))
    def test_the_key_is_a_field_of_the_constraints(self, key):
        assert key in PoolConstraints().__dataclass_fields__


class TestBuildingConstraintsFromConfiguration:
    def test_no_limit_set_yields_no_constraints_at_all(self):
        """`None`, not an empty `PoolConstraints`. A run that asks for nothing
        must not acquire an objective, because building one is what would change
        the delivered panel."""
        assert constraints_from_parameter(_Params()) is None

    def test_an_explicit_none_is_still_nothing(self):
        assert constraints_from_parameter(_Params(max_worst_hole=None)) is None

    def test_one_limit_set_yields_constraints_carrying_it(self):
        constraints = constraints_from_parameter(_Params(max_worst_hole=20_000.0))

        assert constraints is not None
        assert constraints.max_worst_hole == 20_000.0
        assert constraints.max_evenness is None

    def test_every_key_reaches_the_field_of_the_same_name(self):
        """A transposed pair here would enforce the wrong limit silently."""
        values = {key: float(index + 2) for index, key in enumerate(sorted(LIMIT_KEYS))}

        constraints = constraints_from_parameter(_Params(**values))

        for key, value in values.items():
            assert getattr(constraints, key) == value, f"{key} did not reach its field"


class _Metrics:
    def __init__(self, max_gap=31_000.0, **kw):
        self.effective_fg_coverage = kw.get("effective_fg_coverage", 0.8)
        self.fg_coverage = 0.9
        self.selectivity_density = 45.0
        self.total_bg_sites = 120
        self.max_gap = max_gap
        self.mean_gap = 4_800.0
        self.gap_gini = 0.48
        self.bg_coverage = 0.02


class _Optimizer:
    def __init__(self, metrics):
        self._metrics = metrics

    def compute_metrics(self, primers):
        return self._metrics


PANEL = ["ACGTACGTAC", "TTGCATGCAT", "GGCCTTAAGG"]


class TestEnforcement:
    def test_a_panel_inside_every_limit_is_reported_as_meeting_them(self):
        report = enforce_constraints(
            PANEL,
            _Optimizer(_Metrics()),
            candidates=PANEL,
            constraints=PoolConstraints(max_worst_hole=40_000.0),
            config=None,
        )

        assert report.violations == ()
        assert report.primers == PANEL

    def test_a_panel_over_a_limit_names_the_limit_it_missed(self):
        report = enforce_constraints(
            PANEL,
            _Optimizer(_Metrics(max_gap=31_000.0)),
            candidates=PANEL,
            constraints=PoolConstraints(max_worst_hole=20_000.0),
            config=None,
        )

        assert "worst hole above maximum" in report.violations

    def test_the_panel_is_never_replaced_by_a_worse_one(self):
        """A repair that does not resolve the violation returns its input. The
        same rule `pool_planner._repair` follows: on a limit no panel can meet,
        chasing it would trade real coverage for a step toward a floor it never
        reaches."""
        report = enforce_constraints(
            PANEL,
            _Optimizer(_Metrics(max_gap=10_000_000.0)),
            candidates=PANEL,
            constraints=PoolConstraints(max_worst_hole=1.0),
            config=None,
        )

        assert report.primers == PANEL
        assert report.violations

    def test_the_report_formats_to_lines_naming_the_values(self):
        report = enforce_constraints(
            PANEL,
            _Optimizer(_Metrics(max_gap=31_000.0)),
            candidates=PANEL,
            constraints=PoolConstraints(max_worst_hole=20_000.0),
            config=None,
        )
        lines = "\n".join(report.lines())

        assert "worst hole above maximum" in lines
        assert "20" in lines


class TestTheOptimizePathReachesIt:
    def test_run_optimization_builds_constraints_from_configuration(self):
        """Asserted on the call, because both ends existing and the path not is
        exactly Known Issue 16."""
        import ast
        import inspect

        from neoswga.core import unified_optimizer

        tree = ast.parse(inspect.getsource(unified_optimizer))
        names = {
            node.func.id
            for node in ast.walk(tree)
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
        }

        assert "constraints_from_parameter" in names
        assert "apply_configured_limits" in names

    def test_the_acceptance_module_is_what_calls_the_objective(self):
        """`apply_configured_limits` is the seam, so the enforcement itself has
        to happen behind it rather than being a name nobody invokes."""
        import ast
        import inspect

        from neoswga.core import panel_acceptance

        tree = ast.parse(inspect.getsource(panel_acceptance))
        names = {
            node.func.id
            for node in ast.walk(tree)
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
        }

        assert "repair_result" in names
        assert "objective_for_optimizer" in names
