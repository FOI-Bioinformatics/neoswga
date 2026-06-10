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
        a for a in parser._actions
        if a.__class__.__name__ == "_SubParsersAction"
    )
    optimize = subparsers_action.choices["optimize"]
    action = next(
        a for a in optimize._actions if "--optimization-method" in a.option_strings
    )
    return set(action.choices)


def test_optimize_choices_match_registry_plus_ensemble():
    from neoswga.core.optimizer_factory import OptimizerFactory
    from neoswga.core import unified_optimizer as uo

    uo._ensure_optimizers_registered()
    registered = set(OptimizerFactory.list_optimizers().keys())

    choices = _optimize_method_choices()

    assert choices == registered | {"ensemble"}, (
        f"optimize --optimization-method choices {choices} must equal the "
        f"registered optimizers {registered} plus 'ensemble'."
    )


def test_ensemble_default_methods_are_all_registered():
    """The ensemble default fan-out must reference only registered methods."""
    from neoswga.core.optimizer_factory import OptimizerFactory
    from neoswga.core import unified_optimizer as uo

    uo._ensure_optimizers_registered()
    registered = set(OptimizerFactory.list_optimizers().keys())
    default_methods = {"hybrid", "dominating-set", "network", "background-aware"}
    assert default_methods.issubset(registered)
