"""`show-presets` lists what `--preset` applies, and nothing it does not.

It printed `Optimization: hybrid` for four presets, from a display table of
optimizer hints that `--preset` never applies, and listed the remaining seven
presets through a second code path (2026-11-23 sweep, tier 3).
"""

import pytest

from neoswga.cli._common import load_preset_conditions
from neoswga.cli.analysis import show_presets
from neoswga.cli_unified import create_parser


def _preset_choices():
    parser = create_parser()
    sub = next(a for a in parser._actions if a.__class__.__name__ == "_SubParsersAction")
    return next(a for a in sub.choices["filter"]._actions if a.dest == "preset").choices


@pytest.fixture(scope="module")
def listing():
    import contextlib
    import io

    buffer = io.StringIO()
    with contextlib.redirect_stdout(buffer):
        show_presets()
    return buffer.getvalue()


def test_no_optimizer_hint_is_shown(listing):
    assert "Optimization" not in listing


@pytest.mark.parametrize("name", list(_preset_choices()))
def test_each_preset_shows_the_reaction_it_applies(listing, name):
    block = listing.split(f"\n{name}:\n", 1)[1].split("\n\n", 1)[0]
    applied = load_preset_conditions(name)
    assert f"Temperature: {applied['reaction_temp']}C" in block
    assert f"Polymerase: {applied['polymerase']}" in block
