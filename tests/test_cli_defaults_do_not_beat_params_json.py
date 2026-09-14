"""An absent CLI flag must not override a configured value.

Known Issue 8 records this defect class and names the fix: give the flag an
argparse default of `None`, so an explicit flag and an absent one can be told
apart. Three flags still carried a real default and so beat params.json on
every run.
"""

from types import SimpleNamespace

import pytest

from neoswga.cli_unified import create_parser


def _flag_default(command, dest):
    parser = create_parser()
    action = next(
        a
        for sub in parser._subparsers._group_actions
        for name, sp in sub.choices.items()
        if name == command
        for a in sp._actions
        if a.dest == dest
    )
    return action.default


@pytest.mark.parametrize(
    "command,dest",
    [
        ("filter", "gc_tolerance"),
        ("filter", "excl_threshold"),
        ("expand-primers", "optimization_method"),
        ("plan-pool", "swap_max_evaluations"),
    ],
)
def test_the_flag_uses_the_none_sentinel(command, dest):
    """A real default cannot be distinguished from the user asking for it."""
    assert _flag_default(command, dest) is None


def test_an_absent_gc_tolerance_leaves_the_configured_window_alone():
    from neoswga.cli.pipeline import _resolve_gc_window

    parameter = SimpleNamespace(genome_gc=0.19, gc_min=0.05, gc_max=0.40)
    _resolve_gc_window(SimpleNamespace(gc_tolerance=None, gc_min=None), parameter)
    assert (parameter.gc_min, parameter.gc_max) == (0.05, 0.40)


def test_an_explicit_gc_tolerance_keeps_the_extreme_at_release():
    """The block had its own formula, clamping where the real one releases.

    `adaptive_gc_window` drops the lower bound to zero below the extreme-AT
    threshold, because published AT-rich designs use primers of zero GC. The
    old block clamped at 0.20 and excluded exactly those.
    """
    from neoswga.cli.pipeline import _resolve_gc_window

    parameter = SimpleNamespace(genome_gc=0.19, gc_min=0.375, gc_max=0.625)
    _resolve_gc_window(SimpleNamespace(gc_tolerance=0.15, gc_min=None), parameter)
    assert parameter.gc_min == 0.0


def test_an_explicit_gc_min_still_wins_over_the_tolerance():
    from neoswga.cli.pipeline import _resolve_gc_window

    parameter = SimpleNamespace(genome_gc=0.19, gc_min=0.30, gc_max=0.60)
    _resolve_gc_window(SimpleNamespace(gc_tolerance=0.15, gc_min=0.30), parameter)
    assert (parameter.gc_min, parameter.gc_max) == (0.30, 0.60)


def test_expand_primers_routes_the_method_through_the_resolver():
    """Reading `args.optimization_method` directly bypasses params.json.

    The sentinel above is only half the fix: with a `None` default, a call site
    that still reads the attribute straight through passes `None` to the
    expander instead of the configured method.
    """
    import ast
    import pathlib

    source = pathlib.Path("neoswga/cli/iterate.py").read_text()
    tree = ast.parse(source)
    direct = [
        node
        for node in ast.walk(tree)
        if isinstance(node, ast.Attribute)
        and node.attr == "optimization_method"
        and isinstance(node.value, ast.Name)
        and node.value.id == "args"
    ]
    assert not direct, (
        "neoswga/cli/iterate.py reads args.optimization_method directly; "
        "call resolve_optimization_method(args) so params.json can reach it"
    )
    assert "resolve_optimization_method(args)" in source
