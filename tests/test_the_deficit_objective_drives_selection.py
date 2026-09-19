"""Measuring the deficit is not targeting it.

The previous increment gave expansion a per-base deficit and a way to score a
panel by how much of it a panel recovers, and selection went on ranking by the
ordinary global coverage. That is Known Issue 14's shape exactly: "a
background-aware stage that does not choose the panel changes nothing useful",
and the lesson recorded there is to check which stage produces the delivered
result before concluding that a measurement reaching the code reached the user.

So this pins the PATH, not the pieces. `DeficitObjective` composes
`PoolObjective` rather than replacing it -- same metrics, same violations, same
shortfall, with `coverage` redefined as recovered deficit -- which is what lets
`refine_by_swaps` drive it unchanged, since its score is
`(-shortfall, coverage, -total_bg_sites)`.
"""

import numpy as np
import pytest

from neoswga.core.deficit_objective import DeficitObjective


class _Metrics:
    def __init__(self, bg_sites=7):
        self.effective_fg_coverage = 0.5
        self.fg_coverage = 0.5
        self.total_bg_sites = bg_sites
        self.selectivity_density = 40.0


class _Inner:
    """A stand-in `PoolObjective`, recording what it was asked."""

    def __init__(self):
        self.asked = []

    def metrics(self, primers):
        self.asked.append(("metrics", tuple(primers)))
        return _Metrics()

    def coverage(self, primers):
        self.asked.append(("coverage", tuple(primers)))
        return 0.25

    def violations(self, primers):
        self.asked.append(("violations", tuple(primers)))
        return ("selectivity below minimum",)

    def shortfall(self, primers):
        self.asked.append(("shortfall", tuple(primers)))
        return 0.75


class _Cache:
    def __init__(self, positions):
        self._positions = positions

    def get_positions(self, prefix, primer, strand="both"):
        return np.asarray(self._positions.get((prefix, primer), []), dtype=np.int64)

    def get_record_starts(self, prefix):
        return []


def _objective(inner=None, positions=None, deficit=(400, 600), length=1000):
    weights = np.zeros(length)
    weights[deficit[0] : deficit[1]] = 1.0
    return DeficitObjective(
        inner or _Inner(),
        cache=_Cache(positions or {}),
        weights_by_prefix={"fg": weights},
        lengths_by_prefix={"fg": length},
        extension=50,
        circular=False,
    )


class TestItComposesRatherThanReplaces:
    """Everything except coverage must reach the wrapped objective untouched.

    Replacing it would drop the specificity floor, the dimer accounting and
    the constraint ordering that keeps a feasible panel ahead of an infeasible
    one, and expansion would quietly stop honouring any of them.
    """

    def test_violations_are_the_inner_objective_s(self):
        inner = _Inner()

        assert _objective(inner).violations(["A"]) == ("selectivity below minimum",)

    def test_shortfall_is_the_inner_objective_s(self):
        inner = _Inner()

        assert _objective(inner).shortfall(["A"]) == 0.75

    def test_metrics_are_the_inner_objective_s(self):
        inner = _Inner()

        assert _objective(inner).metrics(["A"]).total_bg_sites == 7

    def test_coverage_is_not(self):
        """The one thing it overrides. The inner objective's 0.25 must not
        survive, or nothing about selection has changed."""
        inner = _Inner()
        objective = _objective(inner, positions={("fg", "A"): [500]})

        assert objective.coverage(["A"]) != 0.25


class TestCoverageIsRecoveredDeficit:
    def test_a_panel_over_the_deficit_scores_high(self):
        objective = _objective(positions={("fg", "A"): [500]})

        assert objective.coverage(["A"]) == pytest.approx(100 / 200)

    def test_a_panel_missing_it_scores_zero(self):
        objective = _objective(positions={("fg", "A"): [50]})

        assert objective.coverage(["A"]) == 0.0

    def test_more_recovery_scores_higher(self):
        objective = _objective(positions={("fg", "NEAR"): [420], ("fg", "MID"): [500]})

        assert objective.coverage(["MID"]) >= objective.coverage(["NEAR"])

    def test_it_is_a_fraction(self):
        objective = _objective(positions={("fg", "A"): [450, 550]})

        value = objective.coverage(["A"])
        assert 0.0 <= value <= 1.0

    def test_no_deficit_anywhere_is_not_a_division_by_zero(self):
        weights = np.zeros(1000)
        objective = DeficitObjective(
            _Inner(),
            cache=_Cache({("fg", "A"): [500]}),
            weights_by_prefix={"fg": weights},
            lengths_by_prefix={"fg": 1000},
            extension=50,
            circular=False,
        )

        assert objective.coverage(["A"]) == 0.0


class TestTheSwapScoreCanDriveIt:
    def test_it_answers_everything_refine_by_swaps_asks(self):
        """`_score` is `(-shortfall, coverage, -metrics(panel).total_bg_sites)`.
        A missing one of those is an AttributeError mid-search."""
        objective = _objective(positions={("fg", "A"): [500]})

        assert objective.shortfall(["A"]) is not None
        assert objective.coverage(["A"]) is not None
        assert objective.metrics(["A"]).total_bg_sites is not None

    def test_a_feasible_panel_still_outranks_an_infeasible_one(self):
        """Redefining coverage must not let it outbid a constraint, which is
        the ordering `shortfall` exists to protect."""

        class _Feasible(_Inner):
            def shortfall(self, primers):
                return 0.0

        feasible = _objective(_Feasible(), positions={("fg", "A"): [50]})
        infeasible = _objective(_Inner(), positions={("fg", "B"): [500]})

        assert -feasible.shortfall(["A"]) > -infeasible.shortfall(["B"])


class TestItReachesTheExpansionPath:
    def test_expansion_attaches_it_to_the_optimizer(self):
        """Asserted on the call. `plan_pool` assigned an objective to a wrapper
        while the refinement read it off the delegate, and the path did not
        exist for months; `attach_search_config` is what fixed that."""
        import ast
        import inspect

        from neoswga.core import primer_expansion

        tree = ast.parse(inspect.getsource(primer_expansion))
        names = {
            node.func.id
            for node in ast.walk(tree)
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
        }

        assert "DeficitObjective" in names
        assert "attach_search_config" in names

    def test_expansion_asks_for_the_refinement_that_reads_an_objective(self):
        """`network`, the default Stage 2, does not read `pool_objective` at
        all. Attaching an objective to a run that cannot consult it is the
        defect, not the fix."""
        import inspect

        from neoswga.core.primer_expansion import PrimerExpander

        source = inspect.getsource(PrimerExpander._expand_hybrid)

        assert '"swap" if wants_deficit else "network"' in source

    def test_a_host_aware_run_keeps_the_refinement_that_carries_the_host_term(self):
        """The two Stage 2s carry different things and neither carries both.

        `network` holds the host term Known Issue 14 added, which is the whole
        of what `background-aware` buys, and
        `test_expansion_uses_the_background` fails if it is taken away. `swap`
        is the one that reads a `pool_objective`. So a host-aware expansion
        does not yet target the deficit, and this pins that the trade-off is
        made deliberately rather than by whichever was wired last.
        """
        import inspect

        from neoswga.core.primer_expansion import PrimerExpander

        source = inspect.getsource(PrimerExpander._expand_hybrid)

        assert "not background_pruning" in source, (
            "the deficit objective must not be attached to a host-aware run, "
            "whose Stage 2 cannot read it"
        )

    def test_the_limitation_is_said_out_loud(self):
        """A run that narrows the pool to the gaps and then ranks by something
        else must say so, or it reads as a gap-targeted design."""
        import inspect

        from neoswga.core.primer_expansion import PrimerExpander

        source = inspect.getsource(PrimerExpander._expand_hybrid)

        assert "logger.warning" in source
