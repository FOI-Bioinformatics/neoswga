"""`background_aggregate` says what "the background" is, and must be recorded.

Phase 7a. A design that judges each host on its own length is a different
design from one that pools the hosts, so the key has to reach the resolved
request, the request's hash and the run manifest before anything acts on it.
Two keys added earlier today were missing from exactly those two places, which
is why they are asserted here rather than when a gate starts reading the value.

**Updated by Phase 7b.** 7a plumbed the key and nothing read it, so the one
behavioural claim available then was that nothing moved. The step-2 background
gate now reads it, and
`tests/test_the_background_gate_reads_the_aggregate.py` is where that
behaviour is asserted; the non-effect test here was replaced by
`test_the_two_readings_differ_only_in_this_field_of_the_request`, which keeps
the part that is still true and still useful. Step 4's per-host figures remain
unwired, and that half of the 7a position is now stated alone at
`tests/test_design_options_have_effect.py::test_background_aggregate_moves_no_step4_figure_before_7c`.
"""

import dataclasses
import json

import pytest

from neoswga.core import background_aggregate as ba

BASE = {
    "fg_genomes": ["a.fna"],
    "fg_prefixes": ["a"],
    "bg_genomes": ["h.fna", "g.fna"],
    "bg_prefixes": ["h", "g"],
    "bg_seq_lengths": [3_348_624, 143_726_002],
    "polymerase": "phi29",
    "reaction_temp": 30.0,
    "min_k": 12,
    "max_k": 12,
    "num_primers": 6,
}


def _request(**overrides):
    from neoswga.core.design_request import resolve_design_request

    return resolve_design_request({**BASE, **overrides})


# ---------------------------------------------------------------------------
# The vocabulary
# ---------------------------------------------------------------------------


def test_the_two_values_are_the_two_the_schema_declares():
    """One vocabulary, in one place. A schema that accepts a third value the
    resolver refuses would refuse it after the user had already been told the
    file was valid."""
    import pathlib

    from neoswga.core import parameter

    schema = json.loads(
        (pathlib.Path(parameter.__file__).parent / "schema" / "params.schema.json").read_text()
    )

    assert tuple(schema["properties"]["background_aggregate"]["enum"]) == ba.BACKGROUND_AGGREGATES
    assert ba.BACKGROUND_AGGREGATES == (ba.SUM, ba.WORST_CASE)


def test_mean_is_not_offered():
    """`host_profile.aggregate_loads` accepts it, because it serves a panel of
    modelled profiles; the key does not, because an average over hosts is
    neither every host at once nor the worst one, and no use for it was shown.
    """
    from neoswga.core.host_profile import aggregate_loads

    assert aggregate_loads([2.0, 3.0], mode="mean") == 2.5
    assert "mean" not in ba.BACKGROUND_AGGREGATES

    with pytest.raises(Exception) as caught:
        ba.resolve_background_aggregate("mean")

    assert "sum" in str(caught.value) and "worst-case" in str(caught.value)


def test_an_unrecognised_value_is_refused_and_the_message_names_both(monkeypatch):
    """Not a fallback to `sum`. Pooling hosts for a user who asked for a worst
    case reports a specificity the panel does not have, which is worse than
    refusing."""
    from neoswga.core.exceptions import InvalidDesignRequest

    with pytest.raises(InvalidDesignRequest) as caught:
        ba.resolve_background_aggregate("max")

    message = str(caught.value)
    assert "background_aggregate" in message
    assert "sum" in message and "worst-case" in message


# ---------------------------------------------------------------------------
# The seam
# ---------------------------------------------------------------------------


def test_the_resolver_reads_the_parameter_global(monkeypatch):
    """The route for code running inside a step, where `get_params` has already
    populated the globals and no request object is in hand."""
    from neoswga.core import parameter

    monkeypatch.setattr(parameter, "background_aggregate", None, raising=False)
    assert ba.resolve_background_aggregate() == ba.SUM

    monkeypatch.setattr(parameter, "background_aggregate", "worst-case", raising=False)
    assert ba.resolve_background_aggregate() == ba.WORST_CASE


def test_an_explicit_argument_beats_the_global(monkeypatch):
    """So a caller that resolved the run's configuration once does not re-read
    a mutable global per evaluation. That is the defect
    `occupancy.weighted_site_load`'s `penalty` argument exists for."""
    from neoswga.core import parameter

    monkeypatch.setattr(parameter, "background_aggregate", "worst-case", raising=False)

    assert ba.resolve_background_aggregate("sum") == ba.SUM


def test_the_key_reaches_the_resolved_request():
    assert _request().background_aggregate == ba.SUM
    assert _request(background_aggregate="sum").background_aggregate == ba.SUM
    assert _request(background_aggregate="worst-case").background_aggregate == ba.WORST_CASE


def test_a_bad_value_is_refused_while_resolving_the_request():
    """Before any index is opened. The request is resolved from the params FILE
    at the start of every design command, which is where a value that cannot be
    applied should be named."""
    from neoswga.core.exceptions import InvalidDesignRequest

    with pytest.raises(InvalidDesignRequest) as caught:
        _request(background_aggregate="worst case")

    assert "background_aggregate" in str(caught.value)


# ---------------------------------------------------------------------------
# Identity
# ---------------------------------------------------------------------------


def test_two_designs_differing_only_in_this_key_have_different_hashes():
    """The defect this file was written against. A key outside the hashed
    fields let two designs that deliver different panels record the same
    identity."""
    pooled = _request(background_aggregate="sum")
    worst = _request(background_aggregate="worst-case")

    assert pooled.request_hash != worst.request_hash


def test_absent_and_sum_are_one_design():
    """They deliver byte-identical designs, and the hash is the identity of a
    design, so they must compare equal. Follows `candidate_retention`, which
    normalises absent to `all_qc`, rather than `optimization_method`, which
    keeps None because a CLI flag has to beat the file. There is no flag here.
    """
    assert _request().request_hash == _request(background_aggregate="sum").request_hash


def test_the_field_is_hashed_rather_than_merely_present():
    """`_HASHED_FIELDS` is derived from the dataclass fields, so this holds by
    construction today. Pinned because the derivation has an exclusion list,
    and a field added to it would stop moving the hash with nothing saying so.
    """
    from neoswga.core.design_request import _HASHED_FIELDS

    assert "background_aggregate" in _HASHED_FIELDS
    assert "background_aggregate" in _request().to_dict()


# ---------------------------------------------------------------------------
# The manifest
# ---------------------------------------------------------------------------


def test_the_manifest_records_which_reading_ran(monkeypatch):
    """A manifest that does not say what the background was describes numbers
    it cannot account for. The RESOLVED value, never None."""
    from neoswga.cli._common import _effective_conditions
    from neoswga.core import parameter

    monkeypatch.setattr(parameter, "background_aggregate", None, raising=False)
    assert _effective_conditions(parameter)["background_aggregate"] == ba.SUM

    monkeypatch.setattr(parameter, "background_aggregate", "worst-case", raising=False)
    assert _effective_conditions(parameter)["background_aggregate"] == ba.WORST_CASE


def test_it_survives_a_write_and_a_read_of_the_manifest(tmp_path, monkeypatch):
    from neoswga.cli._common import _effective_conditions
    from neoswga.core import parameter
    from neoswga.core.run_manifest import read_effective_conditions, write_manifest

    monkeypatch.setattr(parameter, "data_dir", str(tmp_path), raising=False)
    monkeypatch.setattr(parameter, "background_aggregate", "worst-case", raising=False)
    write_manifest(
        step="filter",
        data_dir=str(tmp_path),
        effective_conditions=_effective_conditions(parameter),
    )

    recorded = read_effective_conditions(str(tmp_path), step="filter")

    assert recorded["background_aggregate"] == ba.WORST_CASE


def test_a_change_between_steps_is_condition_drift(tmp_path, monkeypatch):
    """`filter` pooling the hosts and `optimize` taking the worst of them would
    select a pool under one reading of the background and score it under
    another. The drift check already exists; this makes it see the key."""
    from neoswga.cli._common import _effective_conditions, warn_on_condition_drift
    from neoswga.core import parameter
    from neoswga.core.run_manifest import write_manifest

    monkeypatch.setattr(parameter, "data_dir", str(tmp_path), raising=False)
    monkeypatch.setattr(parameter, "background_aggregate", "worst-case", raising=False)
    write_manifest(
        step="filter",
        data_dir=str(tmp_path),
        effective_conditions=_effective_conditions(parameter),
    )

    monkeypatch.setattr(parameter, "background_aggregate", None, raising=False)

    assert "background_aggregate" in warn_on_condition_drift(parameter, reference_step="filter")


# ---------------------------------------------------------------------------
# What the key now moves, and what is still waiting on 7c
# ---------------------------------------------------------------------------


def test_the_two_readings_differ_only_in_this_field_of_the_request():
    """The request is the same design either way except for the key itself.

    This was `test_background_aggregate_does_not_yet_move_a_figure` in 7a,
    where it carried the stronger claim that nothing read the key at all. Phase
    7b retired that claim: the step-2 gate reads it, and
    `tests/test_the_background_gate_reads_the_aggregate.py` is where that is
    asserted. What remains true, and worth pinning on its own, is that the two
    readings change nothing in the RESOLVED REQUEST but the one field -- so any
    difference a stage produces has to come from the stage reading the key, not
    from a second field moving with it.

    Checked on the request rather than by grep: every field a stage reads is
    here, and the only one that differs is the key itself.
    """
    pooled = _request(background_aggregate="sum")
    worst = _request(background_aggregate="worst-case")

    # `conditions` is a freshly built object per resolution and compares by
    # identity; its fingerprint is what the request hashes, so compare that.
    differing = {
        name
        for name in (f.name for f in dataclasses.fields(pooled))
        if name != "conditions" and getattr(pooled, name) != getattr(worst, name)
    }

    assert differing == {"background_aggregate"}
    assert pooled.conditions_fingerprint == worst.conditions_fingerprint


def test_the_step_2_gate_is_the_one_stage_that_reads_it_so_far():
    """7b wired step 2; 7c wires step 4's per-host figures.

    Stated as a positive assertion about the gate rather than as "nothing reads
    it", which is what 7a had to say. The step-4 half is still a non-effect
    claim and lives at
    `tests/test_design_options_have_effect.py::test_background_aggregate_moves_no_step4_figure_before_7c`.
    """
    from neoswga.core import background_gate

    panel = background_gate.host_panel(["small", "large"], [40_000, 2_000_000], ba.WORST_CASE)
    counts = {"small": {"ACGTACGTAC": 1}, "large": {"ACGTACGTAC": 1}}

    # 2 sites over the pooled 2,040,000 bp is 9.8e-7 and passes; 1 site over
    # the small host's own 40,000 bp is 2.5e-5 and does not.
    assert background_gate.pooled_verdict(2, 2_040_000, 5e-6) is True
    assert background_gate.worst_case_verdict(counts, "ACGTACGTAC", panel, 5e-6) is False
