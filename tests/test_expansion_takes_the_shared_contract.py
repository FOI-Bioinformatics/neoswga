"""Expansion goes through `run_panel_search` on every path a command takes.

`docs/validation/alternatives_through_the_contract_2026-09-28.md` said
`primer_expansion._expand_hybrid` "calls `optimizer.optimize` directly, with
expansion-specific arguments the contract does not currently carry", and listed
it as the last path not routed. That reads as a gap in production and is not
one.

`_expand_hybrid` routes to `_expand_configured` -- which builds through
`OptimizerFactory` and calls `run_panel_search` -- whenever the expander has a
design context. Every command that expands supplies one: `expand-primers` and
`iterate` both build their expander with `context=context`. The direct branch is
the fallback for an API caller with no context, and it exists for a reason
rather than by omission: `_expand_configured` reads the optimizer config off
that context, so without one there is nothing to build a factory optimizer from.

So the honest statement is that the contract covers every command path, and one
API shape falls back. These tests pin both halves, because "already routed" and
"cannot be routed here" are different claims and the record confused them.
"""

import ast
import inspect
import pathlib

import pytest

from neoswga.core import primer_expansion


def _calls_in(function_name, module=primer_expansion):
    tree = ast.parse(inspect.getsource(module))
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == function_name:
            return {
                getattr(call.func, "id", None) or getattr(call.func, "attr", None)
                for call in ast.walk(node)
                if isinstance(call, ast.Call)
            }
    raise AssertionError(f"{function_name} not found")


def test_the_configured_expansion_uses_the_shared_contract():
    assert "run_panel_search" in _calls_in("_expand_configured")


def test_a_context_sends_expansion_down_the_contract():
    """The branch, read from the source that decides it."""
    source = inspect.getsource(primer_expansion.PrimerExpander._expand_hybrid)
    guard = source.index("if self.context is not None:")
    routed = source.index("_expand_configured")
    assert guard < routed, "a context must route expansion to the configured path"


@pytest.mark.parametrize("module", ["neoswga/cli/iterate.py"])
def test_every_expanding_command_supplies_a_context(module):
    """`expand-primers` and `iterate` both build their expander with one.

    `cli/commands.py` also builds an expander, for `analyze-coverage`, and that
    one never expands -- it calls `identify_gaps`, which is read-only. So it is
    not covered here and does not need to be.
    """
    source = pathlib.Path(module).read_text()
    construction = source[source.index("PrimerExpander(") :]
    construction = construction[: construction.index(")\n")]
    assert "context=context" in construction, (
        f"{module} builds an expander without a design context, so expansion "
        f"there would bypass the shared contract"
    )


def test_the_direct_branch_needs_no_context_by_design():
    """Why the fallback exists, so it is not mistaken for an oversight.

    `_expand_configured` reads the optimizer config off the design context.
    With no context there is nothing to build a factory optimizer from, which
    is what the direct branch is for.
    """
    source = inspect.getsource(primer_expansion.PrimerExpander._expand_configured)
    assert "self.context.optimizer_config" in source, (
        "if the configured path stops reading the context, the direct branch "
        "loses its reason to exist and should be removed"
    )


def test_the_direct_branch_is_still_exercised():
    """It has one caller in the suite, so it is not dead code either."""
    covering = pathlib.Path("tests/test_expansion_uses_the_background.py")
    assert covering.exists(), (
        "the no-context expansion path lost its only coverage; either restore "
        "it or delete the branch"
    )
    assert (
        "context=" not in covering.read_text()
    ), "that file now passes a context, so nothing exercises the direct branch"
