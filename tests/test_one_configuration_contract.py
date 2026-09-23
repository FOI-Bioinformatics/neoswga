"""Two commands, one params file, one answer.

Finding F7 of the 2026-09-16 pipeline audit, the configuration half.
`plan-pool` resolves the reaction, the coverage reach and the dimer limits from
params.json. `expand-primers` resolved none of them: it ran at a hard-coded
3 kb reach with `conditions=None`, `tm_weight=0.0`, `dimer_penalty=0.0` and
`max_dimer_bp=None`.

So the same params file produced designs under different chemistry depending on
which command was run, and the one that exists to ADD to an existing panel was
the one running without the panel's chemistry.

`core/design_context.py` is the single resolution both use. The test that
matters is the one asserting they agree, because two resolutions that merely
look similar are how they drift.
"""

import pytest

from neoswga.core.design_context import DesignContext, design_context_from_params


def _params(**overrides):
    params = {
        "polymerase": "equiphi29",
        "reaction_temp": 42.0,
        "dmso_percent": 5.0,
        "max_dimer_bp": 4,
        "max_self_dimer_bp": 5,
        "min_tm": 37.0,
        "max_tm": 62.0,
        "fg_circular": True,
        "fg_prefixes": ["fg"],
        "fg_seq_lengths": [1000],
    }
    params.update(overrides)
    return params


class TestItResolvesWhatBothCommandsNeed:
    def test_the_reach_comes_from_the_polymerase(self):
        context = design_context_from_params(_params())

        assert context.coverage_reach > 0

    def test_an_explicit_reach_wins(self):
        context = design_context_from_params(_params(), coverage_reach_override=8000)

        assert context.coverage_reach == 8000

    def test_a_configured_reach_beats_the_polymerase_default(self):
        context = design_context_from_params(_params(coverage_reach=7000))

        assert context.coverage_reach == 7000

    def test_the_reaction_carries_the_additive(self):
        context = design_context_from_params(_params())

        assert context.conditions.dmso_percent == 5.0
        assert context.conditions.temp == 42.0

    def test_the_dimer_limits_come_from_the_file(self):
        context = design_context_from_params(_params())

        assert context.max_dimer_bp == 4
        assert context.max_self_dimer_bp == 5

    def test_the_tm_window_comes_from_the_file(self):
        context = design_context_from_params(_params())

        assert (context.min_tm, context.max_tm) == (37.0, 62.0)

    def test_circularity_comes_from_the_file(self):
        assert design_context_from_params(_params()).fg_circular is True

    def test_an_optional_stability_floor_is_carried(self):
        context = design_context_from_params(_params(max_dimer_dg=-4.0))

        assert context.max_dimer_dg == -4.0

    def test_it_is_frozen(self):
        context = design_context_from_params(_params())

        with pytest.raises(Exception):
            context.coverage_reach = 1


class TestTheTwoCommandsAgree:
    """The point of the module. Two resolutions that look similar drift."""

    def test_the_optimizer_config_carries_the_same_values(self):
        context = design_context_from_params(_params())
        config = context.optimizer_config()

        assert config.max_dimer_bp == context.max_dimer_bp
        assert config.max_self_dimer_bp == context.max_self_dimer_bp
        assert config.min_tm == context.min_tm
        assert config.max_tm == context.max_tm
        assert config.extension_reach == context.coverage_reach
        assert config.fg_circular == context.fg_circular

    def test_an_override_reaches_the_config(self):
        context = design_context_from_params(_params())

        config = context.optimizer_config(refinement_method="swap", verbose=True)

        assert config.refinement_method == "swap"
        assert config.verbose is True

    def test_the_expander_takes_the_same_context(self):
        """`PrimerExpander` defaulted its reach to 3 kb whatever the
        polymerase, so a design made at 8 kb was expanded at 3 kb."""
        import inspect

        from neoswga.core.primer_expansion import PrimerExpander

        assert "context" in inspect.signature(PrimerExpander.__init__).parameters

    def test_the_expander_uses_the_context_reach(self):
        from neoswga.core.primer_expansion import PrimerExpander

        context = design_context_from_params(_params(coverage_reach=7000))
        expander = PrimerExpander(
            position_cache=None, fg_prefixes=["fg"], fg_seq_lengths=[1000], context=context
        )

        assert expander.coverage_reach == 7000

    def test_the_expander_uses_the_context_chemistry(self):
        """It ran with `conditions=None`, so no additive, salt or temperature
        correction reached the primers it added to a panel designed with them.
        """
        from neoswga.core.primer_expansion import PrimerExpander

        context = design_context_from_params(_params())
        expander = PrimerExpander(
            position_cache=None, fg_prefixes=["fg"], fg_seq_lengths=[1000], context=context
        )

        assert expander.conditions is context.conditions

    def test_the_expander_uses_the_context_dimer_limit(self):
        from neoswga.core.primer_expansion import PrimerExpander

        context = design_context_from_params(_params())
        expander = PrimerExpander(
            position_cache=None, fg_prefixes=["fg"], fg_seq_lengths=[1000], context=context
        )

        assert expander.max_dimer_bp == 4

    def test_the_expander_uses_the_context_circularity(self):
        """It read `fg_circular` off nothing, so a circular target's windows
        did not wrap during expansion though they did during design."""
        from neoswga.core.primer_expansion import PrimerExpander

        context = design_context_from_params(_params())
        expander = PrimerExpander(
            position_cache=None, fg_prefixes=["fg"], fg_seq_lengths=[1000], context=context
        )

        assert expander.fg_circular is True


class TestWithoutAContextNothingChanges:
    def test_the_old_defaults_still_apply(self):
        """Three callers construct an expander without one, and a library
        caller may too."""
        from neoswga.core.primer_expansion import PrimerExpander

        expander = PrimerExpander(position_cache=None, fg_prefixes=["fg"], fg_seq_lengths=[1000])

        assert expander.coverage_reach == 3000
        assert expander.conditions is None

    def test_an_explicit_reach_still_wins_over_the_default(self):
        from neoswga.core.primer_expansion import PrimerExpander

        expander = PrimerExpander(
            position_cache=None,
            fg_prefixes=["fg"],
            fg_seq_lengths=[1000],
            coverage_reach=5000,
        )

        assert expander.coverage_reach == 5000


class TestItReachesTheLibraryEntryPoint:
    def test_expand_primer_set_builds_a_context(self):
        """It had `params` in hand the whole time and resolved nothing from it."""
        import ast
        import inspect

        from neoswga.core import primer_expansion

        tree = ast.parse(inspect.getsource(primer_expansion))
        names = {
            node.func.id
            for node in ast.walk(tree)
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
        }

        assert "design_context_from_params" in names

    def test_the_context_type_is_what_is_passed(self):
        context = design_context_from_params(_params())

        assert isinstance(context, DesignContext)
