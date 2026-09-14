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
}

RETIRED = {"retries", "drop_iterations", "top_set_count", "selection_metric", "bl_penalty"}


@pytest.mark.parametrize("key", sorted(SCHEMA["properties"]))
def test_every_schema_key_binds_a_parameter_global(key):
    """`hasattr` is the check, because a missing global is exactly the defect.

    A reader written as `getattr(parameter, name, default)` cannot tell an
    unset global from a configured value equal to the default, which is how
    these stayed invisible.
    """
    from neoswga.core import parameter

    if key in NOT_PARAMETER_GLOBALS:
        return
    assert hasattr(parameter, key), (
        f"{key} is declared in params.schema.json but binds no parameter global, "
        f"so every reader takes its fallback and the key does nothing"
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
