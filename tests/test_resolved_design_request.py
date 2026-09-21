"""One immutable request, resolved once, with every default's source recorded.

Task 2 of the 2026-09-21 valid-design plan. The defect it addresses is not that
any single reader takes a wrong default; it is that the same params file is
resolved independently by several commands, each with its own fallbacks, so
"what chemistry was this designed under" has no single answer to read back.

Three properties are pinned here:

1. **Resolution happens once and is frozen.** A request cannot be edited after
   it is built, and neither can anything nested inside it, so a later stage
   cannot change the chemistry or a threshold the search already used.
2. **A setting is applied or refused, never quietly replaced.** That includes
   an explicit zero, which `or`-style defaulting turns into the default. The
   plan calls this out by name because `design_context_from_params` did exactly
   that to `coverage_reach`.
3. **Equivalent requests hash equal.** The hash is what ties a saved result to
   the configuration that produced it, so it has to be canonical rather than
   dependent on key order or on how a caller spelled a path.
"""

import dataclasses
import json

import pytest

from neoswga.core.design_request import (
    DesignRequest,
    ConcentrationPolicy,
    resolve_design_request,
)
from neoswga.core.exceptions import InvalidDesignRequest, UnsupportedModelError


def base_params(**overrides):
    """A params mapping that resolves cleanly, so a test can break one thing."""
    params = {
        "fg_prefixes": ["target"],
        "fg_seq_lengths": [100000],
        "bg_prefixes": [],
        "bg_seq_lengths": [],
        "polymerase": "phi29",
        "reaction_temp": 30.0,
        "min_k": 12,
        "max_k": 12,
        "num_primers": 12,
        "max_dimer_bp": 3,
        "data_dir": "results",
    }
    params.update(overrides)
    return params


# ---------------------------------------------------------------------------
# 1. The request is immutable
# ---------------------------------------------------------------------------


def test_a_resolved_request_cannot_be_edited():
    request = resolve_design_request(base_params())

    with pytest.raises(dataclasses.FrozenInstanceError):
        request.coverage_reach = 10000


def test_nested_content_is_immutable_too():
    """A frozen dataclass holding a list is not a frozen request.

    The stage that matters is the one that appends to a list it was handed.
    """
    request = resolve_design_request(base_params())

    for value in (
        request.fg_prefixes,
        request.bg_prefixes,
        request.fg_seq_lengths,
        request.fixed_oligos,
        request.excluded_oligos,
        request.primer_lengths,
    ):
        assert isinstance(value, tuple), type(value)


def test_runtime_caches_are_not_part_of_the_request():
    """Configuration, not state. A request must be safe to hash and to log."""
    names = {field.name for field in dataclasses.fields(DesignRequest)}

    assert not names & {"position_cache", "cache", "candidates", "optimizer", "inventory"}


# ---------------------------------------------------------------------------
# 2. Settings are applied or refused
# ---------------------------------------------------------------------------


def test_explicit_zero_reach_is_not_replaced_by_a_default():
    with pytest.raises(InvalidDesignRequest, match="coverage_reach"):
        resolve_design_request({"coverage_reach": 0})


def test_an_explicit_reach_survives_resolution():
    assert resolve_design_request(base_params(coverage_reach=8000)).coverage_reach == 8000


def test_an_absent_reach_takes_the_polymerase_reach_and_says_so():
    request = resolve_design_request(base_params())

    assert request.coverage_reach == 3000
    assert request.default_sources["coverage_reach"] == "polymerase:phi29"


def test_a_supplied_value_is_recorded_as_supplied():
    request = resolve_design_request(base_params(coverage_reach=8000))

    assert request.default_sources["coverage_reach"] == "request"


@pytest.mark.parametrize("bad", [float("nan"), float("inf"), float("-inf")])
def test_a_non_finite_value_is_refused(bad):
    with pytest.raises(InvalidDesignRequest, match="finite"):
        resolve_design_request(base_params(reaction_temp=bad))


def test_an_unknown_key_is_refused_rather_than_ignored():
    """`max_bg_freqency` was accepted in silence and the default applied."""
    with pytest.raises(InvalidDesignRequest, match="max_bg_freqency"):
        resolve_design_request(base_params(max_bg_freqency=1e-6))


def test_a_retired_key_names_what_replaced_it():
    with pytest.raises(InvalidDesignRequest, match="candidate_retention"):
        resolve_design_request(base_params(candidate_retention="legacy"))


@pytest.mark.parametrize(
    "key", ["total_search_evaluations", "total_search_seconds", "max_frontier_refills"]
)
def test_a_negative_budget_is_refused(key):
    with pytest.raises(InvalidDesignRequest, match=key):
        resolve_design_request(base_params(**{key: -1}))


def test_an_oligo_cannot_be_both_fixed_and_excluded():
    with pytest.raises(InvalidDesignRequest, match="ACGTACGTACGT"):
        resolve_design_request(
            base_params(
                fixed_oligos=["ACGTACGTACGT"],
                excluded_oligos=["ACGTACGTACGT", "TTTTTTTTTTTT"],
            )
        )


def test_an_unknown_polymerase_is_an_unsupported_model_not_a_bad_field():
    """The two are different remedies: fix the spelling, or model the enzyme."""
    with pytest.raises(UnsupportedModelError):
        resolve_design_request(base_params(polymerase="taq-in-an-isothermal-reaction"))


def test_a_specificity_limit_without_a_background_is_refused():
    """Reported as satisfied is the failure mode; it cannot be measured at all."""
    with pytest.raises(InvalidDesignRequest, match="min_selectivity_density"):
        resolve_design_request(base_params(min_selectivity_density=60.0))


def test_the_same_limit_with_a_background_resolves():
    request = resolve_design_request(
        base_params(
            bg_prefixes=["host"],
            bg_seq_lengths=[3000000],
            min_selectivity_density=60.0,
        )
    )

    assert request.constraints is not None


def test_a_missing_required_field_names_the_field():
    with pytest.raises(InvalidDesignRequest, match="fg_prefixes"):
        resolve_design_request({"polymerase": "phi29"})


# ---------------------------------------------------------------------------
# 3. The hash identifies the configuration
# ---------------------------------------------------------------------------


def test_key_order_does_not_change_the_hash():
    forward = base_params()
    reversed_order = dict(reversed(list(forward.items())))

    assert (
        resolve_design_request(forward).request_hash
        == resolve_design_request(reversed_order).request_hash
    )


def test_a_changed_setting_changes_the_hash():
    baseline = resolve_design_request(base_params()).request_hash

    assert resolve_design_request(base_params(coverage_reach=8000)).request_hash != baseline
    # Stays inside phi29's supported 20-40 C: an out-of-range temperature is
    # refused before it can reach a hash, which is a different test.
    assert resolve_design_request(base_params(reaction_temp=35.0)).request_hash != baseline


def test_the_hash_is_stable_across_processes():
    """Built from a canonical JSON form, not from `hash()`, which is salted."""
    request = resolve_design_request(base_params())

    assert len(request.request_hash) == 64
    assert request.request_hash == resolve_design_request(base_params()).request_hash


def test_the_request_serializes_for_the_result_file():
    payload = resolve_design_request(base_params()).to_dict()

    assert json.loads(json.dumps(payload)) == payload
    assert payload["request_hash"] == resolve_design_request(base_params()).request_hash


# ---------------------------------------------------------------------------
# 4. Concentration policy is explicit
# ---------------------------------------------------------------------------


def test_the_default_policy_is_recorded_rather_than_assumed():
    request = resolve_design_request(base_params())

    assert request.concentration_policy.mode in {"per_oligo", "fixed_total"}
    assert request.concentration_policy.molar > 0


def assert_concentration_conserved(request, oligos):
    values = request.concentrations_molar(tuple(oligos))
    assert len(values) == len(oligos)
    assert sum(values) == pytest.approx(request.total_primer_molar)
    assert all(value > 0 for value in values)


@pytest.mark.parametrize("count", [2, 4])
def test_a_fixed_total_is_shared_out_and_conserved(count):
    request = resolve_design_request(
        base_params(concentration_mode="fixed_total", total_primer_molar=4e-6)
    )

    assert_concentration_conserved(request, ["ACGTACGTACGT"[:12]] * 0 + [f"A{'CGT' * 3}{i:02d}" for i in range(count)])


def test_fixed_per_oligo_keeps_each_concentration_while_the_total_moves():
    request = resolve_design_request(
        base_params(concentration_mode="per_oligo", primer_conc=2e-6)
    )

    two = request.concentrations_molar(("AAAACCCCGGGG", "TTTTGGGGCCCC"))
    four = request.concentrations_molar(
        ("AAAACCCCGGGG", "TTTTGGGGCCCC", "ACACACACACAC", "GTGTGTGTGTGT")
    )

    assert set(two) == {2e-6}
    assert set(four) == {2e-6}
    assert sum(four) == pytest.approx(2 * sum(two))


def test_an_unknown_concentration_mode_is_refused():
    with pytest.raises(InvalidDesignRequest, match="concentration_mode"):
        resolve_design_request(base_params(concentration_mode="whatever"))


def test_the_policy_is_immutable():
    policy = resolve_design_request(base_params()).concentration_policy

    assert isinstance(policy, ConcentrationPolicy)
    with pytest.raises(dataclasses.FrozenInstanceError):
        policy.molar = 1e-6
