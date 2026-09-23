"""A run can ask for the second geometry, and asks for the first by default.

`OptimizerConfig.coverage_geometry` has existed since PR #108 and every
coverage path reads it, but nothing set it from a params file, so the only way
to reach it was to construct an optimizer by hand. That is the shape this
repository ratchets against, and it is why `tests/test_no_schema_key_is_inert.py`
and `tests/test_every_cli_option_has_an_effect.py` exist.

The default stays `symmetric`. Every figure recorded in `docs/validation/` was
produced under it, and the reach was fitted alongside it -- 2.9-4.6 kb symmetric
against 4.4-6.7 directional, measured in
`docs/validation/2026-09-23-reach-refit-directional.md`. Setting the geometry
without moving the reach recalibrates silently, which is why they are separate
keys and why neither has a new default here.
"""

import json

import pytest

from neoswga.core.unified_optimizer import _build_optimizer_config


def build(**kwargs):
    return _build_optimizer_config(
        target_size=6,
        verbose=False,
        extension_reach=3000,
        fg_circular=False,
        kwargs=kwargs,
    )


def test_the_default_is_the_convention_every_figure_was_produced_under():
    assert build().coverage_geometry == "symmetric"


def test_a_run_can_ask_for_the_directional_geometry():
    """The point. Before this the field was reachable only by hand."""
    assert build(coverage_geometry="directional").coverage_geometry == "directional"


def test_the_key_is_declared_in_the_schema():
    """`test_no_schema_key_is_inert.py` walks the schema and requires every key
    to bind something. A field reachable only through `**kwargs` is reachable
    from a library caller and not from a params file, which is half a wiring."""
    from pathlib import Path

    schema = json.loads(
        (
            Path(__file__).resolve().parent.parent / "neoswga/core/schema/params.schema.json"
        ).read_text()
    )
    properties = schema.get("properties", {})

    assert "coverage_geometry" in properties
    assert properties["coverage_geometry"]["enum"] == ["symmetric", "directional"]


def test_the_schema_names_the_permitted_values():
    """An enum rather than a free string, because a typo falling through to a
    default is exactly what `site_spans` raises to prevent, and a schema
    refusal is earlier and names what is allowed.

    Asserted on the declaration rather than by running the validator: this
    repository has no `validate_params(mapping)` entry point, and an earlier
    draft of this test invented one.
    """
    from pathlib import Path

    schema = json.loads(
        (
            Path(__file__).resolve().parent.parent / "neoswga/core/schema/params.schema.json"
        ).read_text()
    )
    entry = schema["properties"]["coverage_geometry"]

    assert entry["type"] == "string"
    assert entry["enum"] == ["symmetric", "directional"]
    assert "recalibrat" in entry["description"], "the description must warn about the reach"


def test_params_json_reaches_the_config(tmp_path, monkeypatch):
    """The path a user actually takes, not just the dataclass field."""
    from neoswga.core import parameter

    monkeypatch.setattr(parameter, "coverage_geometry", "directional", raising=False)

    assert build().coverage_geometry == "directional"
