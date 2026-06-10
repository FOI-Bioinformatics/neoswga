"""Guard against CLI/registry drift.

The `optimize --optimization-method` choices must equal the registered
optimizers plus the `ensemble` pseudo-method. This catches the class of bug
where the help/choices over-promise methods that were removed (greedy, milp,
genetic, moea, ...) or under-list new ones.
"""


def _optimize_method_choices():
    from neoswga.cli_unified import create_parser

    parser = create_parser()
    # Find the optimize subparser's --optimization-method action.
    subparsers_action = next(
        a for a in parser._actions if a.__class__.__name__ == "_SubParsersAction"
    )
    optimize = subparsers_action.choices["optimize"]
    action = next(a for a in optimize._actions if "--optimization-method" in a.option_strings)
    return set(action.choices)


def test_optimize_choices_match_registry_plus_ensemble():
    from neoswga.core import unified_optimizer as uo
    from neoswga.core.optimizer_factory import OptimizerFactory

    uo._ensure_optimizers_registered()
    registered = set(OptimizerFactory.list_optimizers().keys())

    choices = _optimize_method_choices()

    assert choices == registered | {"ensemble"}, (
        f"optimize --optimization-method choices {choices} must equal the "
        f"registered optimizers {registered} plus 'ensemble'."
    )


def test_ensemble_default_methods_are_all_registered():
    """The ensemble default fan-out must reference only registered methods."""
    from neoswga.core import unified_optimizer as uo
    from neoswga.core.optimizer_factory import OptimizerFactory

    uo._ensure_optimizers_registered()
    registered = set(OptimizerFactory.list_optimizers().keys())
    default_methods = {"hybrid", "dominating-set", "network", "background-aware"}
    assert default_methods.issubset(registered)


def test_method_guide_lists_all_registered_methods():
    """The printed method guide must mention every registered method + ensemble,
    so the user-facing guide does not drift from the actual choices."""
    import io
    from contextlib import redirect_stdout

    from neoswga.cli_unified import print_method_guide
    from neoswga.core import unified_optimizer as uo
    from neoswga.core.optimizer_factory import OptimizerFactory

    uo._ensure_optimizers_registered()
    buf = io.StringIO()
    with redirect_stdout(buf):
        print_method_guide()
    text = buf.getvalue()
    for method in set(OptimizerFactory.list_optimizers().keys()) | {"ensemble"}:
        assert method in text, f"method guide omits '{method}'"


def test_ensemble_combine_choices_documented():
    """--ensemble-combine choices must stay {best, union} (guards doc drift)."""
    choices = _ensemble_combine_choices()
    assert choices == {"best", "union"}


def _ensemble_combine_choices():
    from neoswga.cli_unified import create_parser

    parser = create_parser()
    subparsers_action = next(
        a for a in parser._actions if a.__class__.__name__ == "_SubParsersAction"
    )
    optimize = subparsers_action.choices["optimize"]
    action = next(a for a in optimize._actions if "--ensemble-combine" in a.option_strings)
    return set(action.choices)
