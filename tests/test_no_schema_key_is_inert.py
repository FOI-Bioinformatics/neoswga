"""A documented configuration key must reach something, or not be documented.

CLAUDE.md Known Issue 8 describes the class: a key that is declared in
`params.schema.json`, accepted by the validator, rendered into the parameter
reference, and read by nothing. `additionalProperties` was true when it was
written, so none of them warned.

Six survived the audit of 2026-09-14. Five never bound a module global at all,
so every reader took its fallback; one bound a global that no scoring code read.
They are closed here in the two ways available: `mismatch_penalty` is wired to
the consumer already written for it, and the rest are removed from the schema so
setting one produces the unknown-key warning rather than silence.

This test is the ratchet. A key added to the schema from now on must either
bind a parameter global or be listed below as deliberately parameter-free.
"""

import json
import pathlib

import pytest

SCHEMA = json.loads(
    (pathlib.Path("neoswga") / "core" / "schema" / "params.schema.json").read_text()
)

# Keys that legitimately never become a `parameter` global: they are consumed
# during loading, name files rather than settings, or are metadata.
NOT_PARAMETER_GLOBALS = {
    "schema_version",
    "fg_genomes",
    "bg_genomes",
    "fg_prefixes",
    "bg_prefixes",
    "fg_seq_lengths",
    "bg_seq_lengths",
    "data_dir",
    "src_dir",
    "json_file",
    # Read from the loaded `data` inside `get_params` itself, where it gates
    # whether the adaptive GC window is computed at all. It has effect without
    # becoming a global, which is the distinction this list exists to draw.
    "adaptive_gc",
    # Consumed by `design_request.resolve_design_request`, not by `parameter`.
    # The request is resolved from the params FILE before any global is
    # populated, which is the whole reason it exists, so these two cannot bind
    # a global by construction. Declared in the schema as of 2026-09-22 because
    # they were accepted through an internal whitelist while invisible to every
    # schema-driven validator. `fixed_total` is currently REFUSED, so
    # `total_primer_molar` has no effect and its description says so; that is a
    # documented refusal rather than the silent inertness this file guards.
    "concentration_mode",
    "total_primer_molar",
}

RETIRED = {"retries", "drop_iterations", "top_set_count", "selection_metric", "bl_penalty"}


def _keys_bound_by_loading():
    """Names `_apply_params_only_keys` declares global, read from the source.

    Those keys exist as module attributes only AFTER a config is loaded, so
    `hasattr` on a bare import cannot see them. Read statically so this holds
    without loading anything.
    """
    import ast
    import inspect

    from neoswga.core import parameter

    tree = ast.parse(inspect.getsource(parameter))
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == "_apply_params_only_keys":
            return {
                name
                for inner in ast.walk(node)
                if isinstance(inner, ast.Global)
                for name in inner.names
            }
    return set()


@pytest.mark.parametrize("key", sorted(SCHEMA["properties"]))
def test_every_schema_key_binds_a_parameter_global(key):
    """A documented key must bind SOMEWHERE in `parameter.py`.

    There are three binding sites and this used to check only one.
    `hasattr` on a bare import is true for a module-level global, and for a
    dataclass field or a load-assigned key it is true only once something has
    called `get_params` -- so under `pytest -n 8` this passed or failed by
    which worker drew the test. One run in three failed with 16 of these,
    all on the same worker, while the other two were clean.

    Checking all three is deterministic and no weaker: a key bound by none of
    them is exactly the inert key Known Issue 8 describes, and the five wired
    on 2026-09-14 bind through the third site rather than the first.
    """
    from neoswga.core import parameter
    from neoswga.core.parameter import PipelineParameters

    if key in NOT_PARAMETER_GLOBALS:
        return
    bound = (
        hasattr(parameter, key)
        or key in PipelineParameters.__dataclass_fields__
        or key in _keys_bound_by_loading()
    )
    assert bound, (
        f"{key} is declared in params.schema.json but binds no parameter "
        f"global, no PipelineParameters field and nothing in "
        f"_apply_params_only_keys, so every reader takes its fallback and the "
        f"key does nothing"
    )


@pytest.mark.parametrize("key", sorted(RETIRED))
def test_the_retired_keys_are_gone_from_the_schema(key):
    """Removed rather than wired: nothing implements what they name."""
    assert key not in SCHEMA["properties"]


def test_mismatch_penalty_reaches_its_consumer():
    """It had a consumer whose docstring called itself the first one.

    `occupancy.default_mismatch_penalty` read
    `getattr(parameter, "mismatch_penalty", None)`, which was always None, so it
    always returned the hardcoded 4.0.
    """
    from neoswga.core import occupancy, parameter

    original = getattr(parameter, "mismatch_penalty", None)
    try:
        parameter.mismatch_penalty = 7.5
        assert occupancy.default_mismatch_penalty() == 7.5
    finally:
        if original is None:
            delattr(parameter, "mismatch_penalty")
        else:
            parameter.mismatch_penalty = original
