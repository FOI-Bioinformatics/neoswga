"""Guard: main()'s dispatch table must resolve every handler name.

The cli_unified split moved handlers into neoswga/cli/ submodules, re-exported
back into cli_unified for the dispatch dict. If a re-export is dropped, the dict
literal in main() raises NameError the moment any command runs. This test
exercises that path with a harmless no-arg command so a missing handler is caught
immediately (not only when that specific subcommand happens to be invoked).
"""

import sys

import pytest


def test_main_dispatch_table_fully_resolves(monkeypatch):
    from neoswga import cli_unified

    # 'show-presets' is a side-effect-free command that reaches the dispatch
    # dict construction in main(); building that dict references every handler
    # by name, so an unresolved (un-re-exported) handler raises NameError here.
    monkeypatch.setattr(sys, "argv", ["neoswga", "show-presets"])
    # The assertion is simply that building the dispatch dict + running a command
    # does not raise NameError (a dropped re-export). A clean return or
    # SystemExit(0) are both fine.
    try:
        cli_unified.main()
    except SystemExit as exc:
        assert exc.code in (0, None)


def test_all_cli_submodule_handlers_reexported():
    """Every run_* handler defined in a neoswga/cli/*.py submodule must be
    importable from cli_unified (so the dispatch table and external importers
    keep resolving them)."""
    import importlib
    import pkgutil

    from neoswga import cli as cli_pkg
    from neoswga import cli_unified

    missing = []
    for mod in pkgutil.iter_modules(cli_pkg.__path__):
        if mod.name.startswith("_"):
            continue
        m = importlib.import_module(f"neoswga.cli.{mod.name}")
        for attr in dir(m):
            if attr.startswith("run_") and callable(getattr(m, attr)):
                if not hasattr(cli_unified, attr):
                    missing.append(f"{mod.name}.{attr}")
    assert not missing, f"handlers not re-exported from cli_unified: {missing}"
