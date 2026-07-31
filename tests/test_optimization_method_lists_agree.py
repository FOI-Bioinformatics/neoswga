"""Guard: the optimization-method lists must not drift apart.

`optimization_method` is enumerated in three independent places:

  1. ``OptimizerRegistry`` -- populated by the ``@OptimizerFactory.register()``
     decorators, plus the virtual ensemble names handled directly in
     ``unified_optimizer.run_optimization``. This is the real contract.
  2. ``param_validator.VALID_OPTIMIZATION_METHODS`` -- used by ``validate-params``.
  3. The ``optimization_method`` enum in ``core/schema/params.schema.json``.

These drifted before: the validator and the schema both still listed the retired
``greedy``/``genetic``/``milp``/``moea`` methods while omitting ``ensemble``, which
the wizard and the CLI both emit. The result was that a wizard-generated ensemble
config was rejected by the very validator the wizard tells you to run next.

The validator keeps a literal list on purpose -- it must stay importable without the
numpy/scipy optimizer stack -- so this test is what stops it going stale.
"""

import json
from pathlib import Path

import pytest

# The ensemble runner is dispatched by name in unified_optimizer.run_optimization
# rather than registered as an optimizer class.
VIRTUAL_METHODS = {"ensemble", "auto", "all"}


def _registered_methods():
    """Canonical names + aliases, from the registry the factory actually uses."""
    from neoswga.core.optimizer_factory import OptimizerRegistry

    # Importing the optimizer modules is what triggers their register decorators.
    import neoswga.core.background_aware_optimizer  # noqa: F401
    import neoswga.core.dominating_set_adapter  # noqa: F401
    import neoswga.core.hybrid_optimizer  # noqa: F401
    import neoswga.core.network_optimizer  # noqa: F401

    canonical = set(OptimizerRegistry._registry)
    aliases = set(OptimizerRegistry._aliases)
    assert canonical, "no optimizers registered - the register decorators stopped working"
    return canonical | aliases


def _schema_enum():
    schema_path = (
        Path(__file__).resolve().parents[1] / "neoswga" / "core" / "schema" / "params.schema.json"
    )
    schema = json.loads(schema_path.read_text())
    return set(schema["properties"]["optimization_method"]["enum"])


def test_validator_accepts_every_real_method():
    from neoswga.core.param_validator import VALID_OPTIMIZATION_METHODS

    expected = _registered_methods() | VIRTUAL_METHODS
    missing = expected - set(VALID_OPTIMIZATION_METHODS)
    assert (
        not missing
    ), f"param_validator would reject method(s) the pipeline accepts: {sorted(missing)}"


def test_validator_lists_nothing_the_pipeline_rejects():
    from neoswga.core.param_validator import VALID_OPTIMIZATION_METHODS

    allowed = _registered_methods() | VIRTUAL_METHODS
    phantom = set(VALID_OPTIMIZATION_METHODS) - allowed
    assert not phantom, f"param_validator advertises method(s) that do not exist: {sorted(phantom)}"


def test_schema_enum_matches_validator():
    from neoswga.core.param_validator import VALID_OPTIMIZATION_METHODS

    assert _schema_enum() == set(VALID_OPTIMIZATION_METHODS), (
        "params.schema.json optimization_method enum and "
        "param_validator.VALID_OPTIMIZATION_METHODS have diverged"
    )


@pytest.mark.parametrize("method", sorted(VIRTUAL_METHODS))
def test_ensemble_names_are_valid(method):
    """Regression: 'ensemble' was emitted by the wizard but rejected by the validator."""
    from neoswga.core.param_validator import VALID_OPTIMIZATION_METHODS

    assert method in VALID_OPTIMIZATION_METHODS
    assert method in _schema_enum()


def test_cli_choices_are_a_subset_of_valid_methods():
    """Everything argparse offers must survive validate-params."""
    from neoswga.cli_unified import create_parser
    from neoswga.core.param_validator import VALID_OPTIMIZATION_METHODS

    parser = create_parser()
    optimize = parser._subparsers._group_actions[0].choices["optimize"]
    choices = next(a.choices for a in optimize._actions if a.dest == "optimization_method")
    unsupported = set(choices) - set(VALID_OPTIMIZATION_METHODS)
    assert not unsupported, (
        f"CLI offers --optimization-method value(s) that validate-params rejects: "
        f"{sorted(unsupported)}"
    )
