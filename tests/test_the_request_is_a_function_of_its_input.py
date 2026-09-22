"""A resolved design request must depend on its input and nothing else.

`resolve_design_request` promises in its own docstring that "nothing here reads
a `parameter` module global". It does, indirectly: it builds its chemistry
through `build_reaction_conditions`, whose `_resolve` falls through to
`getattr(_parameter, name, None)` for any field the mapping does not carry.

So a params file that omits `betaine_m` picks up whatever a previous command in
the same process left on the module. Measured 22 September 2026 on one mapping:

    clean process      betaine_m 0.0, hash e59a16eb...
    after another run  betaine_m 1.5, hash c396f1cc...

Same input, different chemistry, different identity. This is not only a
provenance problem: the resolved reaction that the design is filtered and
scored under changes with process history.

It is the same shape as the cached-conditions leak fixed in `filter` on
2026-09-21, where a preset run left conditions behind and the next design wrote
its candidate inventory under the wrong reaction fingerprint. That one reached
CI. This one cannot, because a request resolved in a fresh process is correct
and the suite mostly resolves in fresh processes.

The fix is a resolver that does not consult the module at all. Constructor
defaults fill what the mapping omits, which is a documented, stable answer
rather than a historical one.
"""

import pytest

from neoswga.core import parameter
from neoswga.core.design_request import resolve_design_request

MAPPING = {
    "fg_genomes": ["a.fna"],
    "fg_prefixes": ["a"],
    "polymerase": "phi29",
    "reaction_temp": 30.0,
    "min_k": 12,
    "max_k": 12,
    "num_primers": 12,
}


@pytest.fixture
def leftover_chemistry(monkeypatch):
    """What an earlier command in the same process leaves on the module."""
    monkeypatch.setattr(parameter, "betaine_m", 1.5, raising=False)
    monkeypatch.setattr(parameter, "dmso_percent", 7.5, raising=False)


def test_the_same_mapping_resolves_to_the_same_chemistry(leftover_chemistry):
    """The headline. A request must not inherit a previous run's reaction."""
    resolved = resolve_design_request(dict(MAPPING))

    assert resolved.conditions.betaine_m == 0.0, (
        "the request picked up betaine from the parameter module; its chemistry "
        "depends on what ran before it in this process"
    )
    assert resolved.conditions.dmso_percent == 0.0


def test_the_same_mapping_resolves_to_the_same_hash(leftover_chemistry):
    """Identity follows chemistry, so the hash moves for the same reason."""
    polluted = resolve_design_request(dict(MAPPING))

    # A value the module cannot supply, resolved the same way both times.
    clean_equivalent = resolve_design_request({**MAPPING, "betaine_m": 0.0, "dmso_percent": 0.0})

    assert polluted.request_hash == clean_equivalent.request_hash


def test_an_explicit_value_in_the_mapping_still_wins(leftover_chemistry):
    """Not reading the module must not mean ignoring the user."""
    resolved = resolve_design_request({**MAPPING, "betaine_m": 1.0})

    assert resolved.conditions.betaine_m == 1.0


def test_an_omitted_field_takes_the_documented_default(leftover_chemistry):
    """Absence resolves to the constructor's default, not to process history."""
    resolved = resolve_design_request(dict(MAPPING))

    assert resolved.conditions.na_conc == 50.0


# ---------------------------------------------------------------------------
# Explicit zero is a value
# ---------------------------------------------------------------------------


def test_an_explicit_zero_panel_size_is_refused():
    """`target_set_size` used an `or` chain, so 0 silently became 6.

    `_resolve_reach` already rejects an explicit 0 and its docstring explains
    why: `override or params.get(...)` treated a configured 0 as absent, and
    every coverage figure was then reported at a reach nobody asked for. The
    same sentinel-versus-value confusion lived one field away.
    """
    from neoswga.core.exceptions import InvalidDesignRequest

    with pytest.raises(InvalidDesignRequest) as excinfo:
        resolve_design_request({**MAPPING, "target_set_size": 0})

    assert "target_set_size" in str(excinfo.value)


def test_a_negative_panel_size_is_refused():
    from neoswga.core.exceptions import InvalidDesignRequest

    with pytest.raises(InvalidDesignRequest):
        resolve_design_request({**MAPPING, "num_primers": -1})


def test_a_real_panel_size_still_resolves():
    """Guard the guard: a resolver that refuses everything passes the two above."""
    assert resolve_design_request({**MAPPING, "target_set_size": 24}).target_size == 24


# ---------------------------------------------------------------------------
# A declared mode with no execution path
# ---------------------------------------------------------------------------


def test_a_fixed_total_concentration_is_refused():
    """Declared, validated, hashed, and then it changes nothing.

    Measured 2026-09-22: under `fixed_total` with 12 uM across 12 oligos, the
    reaction still reports `primer_conc` 5e-07 and the effective Tm of a
    12-mer is identical to the last digit under both modes. `total_primer_molar`
    is not a `ReactionConditions` field, so it is never forwarded;
    `concentrations_molar` computes the allocation and has no production caller.

    Accepting it is the silent-substitution failure this project's contract
    exists to remove: a user asking for 96 oligos to share 12 uM gets each one
    evaluated at 0.5 uM, a 24-fold error in a quantity that moves Tm about ten
    degrees across that panel-size range.

    Refused rather than warned because the wrong answer is silent and
    consequential. `--use-gpu` warns instead, and that is the right choice
    there: it degrades speed, not correctness.
    """
    from neoswga.core.exceptions import InvalidDesignRequest

    with pytest.raises(InvalidDesignRequest) as excinfo:
        resolve_design_request(
            {**MAPPING, "concentration_mode": "fixed_total", "total_primer_molar": 12e-6}
        )

    message = str(excinfo.value)
    assert "fixed_total" in message
    assert "per_oligo" in message, "the message must name what to use instead"


def test_per_oligo_concentration_still_works():
    resolved = resolve_design_request(
        {**MAPPING, "concentration_mode": "per_oligo", "primer_conc": 0.25e-6}
    )

    assert resolved.conditions.primer_conc == 0.25e-6


def test_the_concentration_keys_are_declared_in_the_schema():
    """They were accepted through a whitelist while absent from the schema, so
    no schema-driven validator or wizard could see them."""
    import json
    import pathlib

    schema = json.loads(
        (
            pathlib.Path(__file__).resolve().parent.parent
            / "neoswga"
            / "core"
            / "schema"
            / "params.schema.json"
        ).read_text()
    )

    for key in ("concentration_mode", "total_primer_molar"):
        assert key in schema["properties"], f"{key} is accepted but undeclared"


# ---------------------------------------------------------------------------
# Frozen means frozen, including what it points at
# ---------------------------------------------------------------------------


def test_mutating_the_conditions_does_not_move_the_hash():
    """`DesignRequest` is a frozen dataclass whose docstring says "It is
    frozen, including its nested content". It was not.

    `conditions` holds a plain mutable class, and the hash folds it in by
    calling `fingerprint()` at hash time. So setting `.temp` on it afterwards
    silently re-identified a record whose whole purpose is to say what a saved
    result was produced under.

    Freezing the CLASS is not the fix: `optimize_conditions_for_primers` and
    `recommend_conditions` both mutate conditions objects in place, on objects
    they build themselves. The fix is that the request captures its identity
    once, at construction.
    """
    resolved = resolve_design_request(dict(MAPPING))
    before = resolved.request_hash

    resolved.conditions.temp = 45.0

    assert resolved.request_hash == before, (
        "mutating the conditions changed the identity of a frozen request"
    )


def test_two_requests_with_different_chemistry_still_differ():
    """Guard the guard: pinning the hash at construction must not make every
    request hash alike."""
    plain = resolve_design_request(dict(MAPPING))
    with_betaine = resolve_design_request({**MAPPING, "betaine_m": 1.0})

    assert plain.request_hash != with_betaine.request_hash


# ---------------------------------------------------------------------------
# The identity must cover what changes the answer
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "key,value",
    [
        ("bg_circular", True),
        ("iterations", 3),
        ("refinement_method", "swap"),
        ("stage1_objective_width", 64),
        ("swap_max_evaluations", 250),
        ("swap_max_seconds", 2.5),
    ],
)
def test_a_setting_that_changes_the_search_changes_the_hash(key, value):
    """These six alter what the search does and left the identity untouched.

    They are declared in the schema and consumed elsewhere, so
    `resolve_design_request` accepted them without storing them: the acceptance
    check runs against the schema, not against the fields the request keeps.
    A run therefore recorded a hash that several settings could not move, which
    defeats the one thing the hash is for.

    `bg_circular` is the plainest case, since `fg_circular` was already a field
    and its background twin was not.
    """
    baseline = resolve_design_request(dict(MAPPING)).request_hash

    assert resolve_design_request({**MAPPING, key: value}).request_hash != baseline, (
        f"{key}={value} changes the search and not the recorded request"
    )


def test_an_unrelated_comment_key_does_not_change_the_hash():
    """Guard the guard: a hash that moves for everything identifies nothing.

    A leading underscore marks a comment and is accepted by the resolver.
    """
    baseline = resolve_design_request(dict(MAPPING)).request_hash

    assert resolve_design_request({**MAPPING, "_note": "for the lab book"}).request_hash == baseline
