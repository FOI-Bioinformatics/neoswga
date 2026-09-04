"""Every field `OptimizationConfig` declares must reach `run_optimization`.

`run_optimization_from_config` forwarded eight of the eighteen fields the
dataclass declares and silently dropped the rest, so a library caller writing

    OptimizationConfig(minimize_primers=True, target_coverage=0.9)

got neither. That also contradicted the comment inside `run_optimization` which
justifies popping `minimize_primers` out of kwargs on the grounds that it
"lives on the separate OptimizationConfig used by run_optimization_from_config".
It does live there. It was not forwarded.

The remaining unforwarded fields are the ones with no counterpart in
`run_optimization`'s signature -- the file paths, and `quiet` -- and they are
listed here explicitly so that adding a nineteenth field without deciding what
to do with it fails rather than joining them silently.
"""

import dataclasses

import pytest

from neoswga.core.unified_optimizer import OptimizationConfig

#: Fields `run_optimization` has no parameter for. `data_dir`, `step3_file` and
#: `output_file` are resolved from the `parameter` module, and `quiet` is the
#: inverse of `verbose`, which IS forwarded.
NOT_FORWARDED = {"method", "data_dir", "step3_file", "output_file", "quiet"}


def test_every_declared_field_is_either_forwarded_or_listed():
    """The list above has to stay in step with the dataclass."""
    declared = {f.name for f in dataclasses.fields(OptimizationConfig)}
    stale = NOT_FORWARDED - declared
    assert not stale, f"NOT_FORWARDED names fields that no longer exist: {stale}"


def test_the_config_reaches_run_optimization(monkeypatch):
    from neoswga.core import unified_optimizer as uo

    seen = {}
    monkeypatch.setattr(uo, "run_optimization", lambda **kwargs: seen.update(kwargs))

    config = OptimizationConfig(
        method="network",
        target_set_size=9,
        max_iterations=17,
        fg_prefixes=["fg"],
        fg_seq_lengths=[10_000],
        bg_prefixes=["bg"],
        bg_seq_lengths=[10_000],
        use_position_cache=False,
        use_background_filter=False,
        verbose=False,
        uniformity_weight=0.4,
        minimize_primers=True,
        target_coverage=0.9,
    )
    uo.run_optimization_from_config(config)

    assert seen["minimize_primers"] is True
    assert seen["target_coverage"] == 0.9
    assert seen["uniformity_weight"] == 0.4
    assert seen["use_cache"] is False
    assert seen["use_background_filter"] is False
    assert seen["max_iterations"] == 17
    assert seen["target_size"] == 9


@pytest.mark.parametrize(
    "field",
    sorted({f.name for f in dataclasses.fields(OptimizationConfig)} - NOT_FORWARDED),
)
def test_each_forwardable_field_actually_travels(monkeypatch, field):
    """One named failing case per dropped field, rather than one lump.

    Parameterised so a field added to the dataclass and forgotten here produces
    a test named after it.
    """
    from neoswga.core import unified_optimizer as uo

    seen = {}
    monkeypatch.setattr(uo, "run_optimization", lambda **kwargs: seen.update(kwargs))
    uo.run_optimization_from_config(OptimizationConfig())

    # The parameter names differ from the field names in three places.
    alias = {
        "target_set_size": "target_size",
        "use_position_cache": "use_cache",
    }
    assert alias.get(field, field) in seen, (
        f"OptimizationConfig.{field} is declared and never forwarded to "
        "run_optimization; either pass it or add it to NOT_FORWARDED with a "
        "reason"
    )
