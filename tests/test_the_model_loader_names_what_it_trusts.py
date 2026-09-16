"""The archive format is version-tolerant; its default trust list is not.

CLAUDE.md Known Issue 2 says the model ships in skops format precisely so a
minor sklearn upgrade does not require retraining. That is true of the FORMAT.
It is not true of which types skops trusts by default, and the two are easy to
conflate.

skops 0.15.0 stopped implicitly trusting `sklearn.tree._tree.Tree`. Every
model-loading test went red on continuous integration while passing locally on
0.14.0, because `pip install -e ".[dev]"` resolves `skops>=0.11,<1` to whatever
is newest and the lock file is not used for the test job.

The loader now states the types a forest legitimately needs and refuses anything
else. That is narrower than trusting the archive wholesale, and it does not
depend on a future release keeping today's defaults. Pinning skops would have
worked until the next release did the same thing.
"""

import pytest

from neoswga.core.rf_preprocessing import _TRUSTED_MODEL_TYPES, _load_skops_model


def _bundled_model_path():
    import pathlib

    import neoswga

    return str(pathlib.Path(neoswga.__file__).parent / "core/models/random_forest_filter.skops")


def test_the_bundled_model_loads():
    model = _load_skops_model(_bundled_model_path())

    assert type(model).__name__ == "RandomForestRegressor"


def test_the_loader_names_the_type_the_new_skops_stopped_trusting():
    """The specific type whose default trust changed."""
    assert "sklearn.tree._tree.Tree" in _TRUSTED_MODEL_TYPES


def test_the_trust_list_is_narrow():
    """Trusting the archive wholesale would defeat the point of the format."""
    assert len(_TRUSTED_MODEL_TYPES) <= 8
    assert all(
        name.startswith(("sklearn.", "numpy.")) for name in _TRUSTED_MODEL_TYPES
    ), "the loader trusts a type from outside sklearn and numpy"


def test_an_unexpected_type_is_reported_by_name():
    """A tampered archive must be refused, and say what it carried.

    Driven through the rule rather than through the loader: which types skops
    reports as untrusted depends on the installed version, so a test that loaded
    a file would exercise this on 0.15 and skip straight past it on 0.14.
    """
    from neoswga.core.rf_preprocessing import unexpected_model_types

    assert unexpected_model_types(["sklearn.tree._tree.Tree"]) == []
    assert unexpected_model_types(["posix.system"]) == ["posix.system"]
    assert unexpected_model_types(["sklearn.tree._tree.Tree", "builtins.eval"]) == ["builtins.eval"]


def test_the_loader_refuses_before_it_loads():
    """The order matters: naming the type after loading it is not a refusal."""
    import ast
    import inspect
    import textwrap

    tree = ast.parse(textwrap.dedent(inspect.getsource(_load_skops_model)))
    steps = [
        node.func.id if isinstance(node.func, ast.Name) else node.func.attr
        for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, (ast.Name, ast.Attribute))
    ]

    assert "unexpected_model_types" in steps
    assert "_skops_load" in steps
    assert steps.index("unexpected_model_types") < steps.index("_skops_load")


def test_the_loader_does_not_trust_blindly():
    """Guard the guard.

    `trusted=True` would load anything and every test above would still pass.
    The loader must pass a list it derived from the file and checked.
    """
    import inspect

    source = inspect.getsource(_load_skops_model)

    assert "trusted=True" not in source, "the loader trusts the archive wholesale"
    assert "get_untrusted_types" in source, "the loader does not inspect what it is trusting"
