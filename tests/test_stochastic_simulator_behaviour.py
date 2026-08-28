"""Behavioural tests for the Gillespie stochastic simulator.

This is the "validation" tier of the simulation stack: it models individual
molecular events (binding, extension, dissociation, resource depletion) to check
what the network-based heuristics predict. That makes its kinetic constants
load-bearing in a way the heuristics' are not -- a rate that is wrong by a
factor of three changes every simulated timing.

It sat at 25% coverage, and one of its constants had indeed drifted: the schema
v2 correction moved phi29's extension rate from 150 to 53 nt/s everywhere except
here.
"""

import math

import pytest

from neoswga.core.registry import views
from neoswga.core.stochastic_simulator import ReactionParameters, ReactionState

REGISTRY_PHI29 = views.mechanistic_enzyme_params()["phi29"]


# ----------------------------------------------------------------------
# Kinetic constants must not drift from the registry
# ----------------------------------------------------------------------


def test_default_extension_rate_matches_the_registry():
    """A stale literal here made the simulator run phi29 ~2.8x too fast.

    `from_conditions` always reads the registry, so only direct construction
    was affected -- which is exactly the path a test or a notebook takes.
    """
    assert ReactionParameters().extension_rate == REGISTRY_PHI29["extension_rate"]


def test_default_processivity_matches_the_registry():
    assert ReactionParameters().max_extension == REGISTRY_PHI29["processivity"]


def test_extension_rate_is_the_corrected_literature_value():
    """53 nt/s (Blanco 1989), not the 150 the module previously carried."""
    assert ReactionParameters().extension_rate == pytest.approx(53, abs=1)


def test_conditions_path_and_default_path_agree_for_phi29():
    """The two ways of getting parameters must not disagree."""
    from neoswga.core.reaction_conditions import ReactionConditions

    from_conditions = ReactionParameters.from_conditions(
        ReactionConditions(temp=30.0, polymerase="phi29"),
        use_mechanistic_model=False,
    )
    assert from_conditions.extension_rate == ReactionParameters().extension_rate
    assert from_conditions.max_extension == ReactionParameters().max_extension


# ----------------------------------------------------------------------
# The kinetic model itself
# ----------------------------------------------------------------------


def test_extension_time_scales_with_distance():
    params = ReactionParameters()
    assert params.extension_time(2_000) == pytest.approx(2 * params.extension_time(1_000))


def test_extension_time_uses_the_configured_rate():
    fast = ReactionParameters(extension_rate=100.0)
    slow = ReactionParameters(extension_rate=50.0)
    assert slow.extension_time(1_000) == pytest.approx(2 * fast.extension_time(1_000))


def test_extension_probability_decays_with_distance():
    """A processive polymerase still falls off; longer runs are less likely."""
    params = ReactionParameters()
    p1 = params.extension_probability(1_000)
    p10 = params.extension_probability(10_000)
    p50 = params.extension_probability(50_000)

    assert 1.0 >= p1 > p10 > p50 > 0.0


def test_extension_probability_matches_the_documented_geometric_model():
    """The docstring states the numbers; they should be reproducible.

    Mean run length = 100 / -ln(processivity_step), which for the phi29 default
    is about 70 kb.
    """
    params = ReactionParameters(processivity_step=0.99857)

    assert params.extension_probability(1_000) == pytest.approx(0.986, abs=0.002)
    assert params.extension_probability(70_000) == pytest.approx(0.368, abs=0.005)

    mean_run = 100 / -math.log(params.processivity_step)
    assert mean_run == pytest.approx(70_000, rel=0.02)


def test_zero_distance_always_completes():
    assert ReactionParameters().extension_probability(0) == pytest.approx(1.0)


# ----------------------------------------------------------------------
# Low-processivity polymerases are flagged
# ----------------------------------------------------------------------


@pytest.mark.parametrize("polymerase", ["bst", "klenow"])
def test_low_processivity_polymerases_warn(polymerase):
    """bst and klenow are poor fits for whole-genome amplification.

    Choosing one silently would produce a plausible-looking simulation of a
    reaction that cannot work.
    """
    from neoswga.core.reaction_conditions import ReactionConditions

    temps = {"bst": 63.0, "klenow": 37.0}
    conditions = ReactionConditions(temp=temps[polymerase], polymerase=polymerase)

    with pytest.warns(UserWarning, match="processivity"):
        ReactionParameters.from_conditions(conditions, use_mechanistic_model=False)


def test_phi29_does_not_warn():
    import warnings

    from neoswga.core.reaction_conditions import ReactionConditions

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        ReactionParameters.from_conditions(
            ReactionConditions(temp=30.0, polymerase="phi29"),
            use_mechanistic_model=False,
        )


def test_each_polymerase_gets_its_own_kinetics():
    """Parameters that do not vary with the enzyme would make the choice inert."""
    import warnings

    from neoswga.core.reaction_conditions import ReactionConditions

    rates = {}
    for name, temp in (("phi29", 30.0), ("equiphi29", 42.0)):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            params = ReactionParameters.from_conditions(
                ReactionConditions(temp=temp, polymerase=name), use_mechanistic_model=False
            )
        rates[name] = params.extension_rate

    assert rates["phi29"] != rates["equiphi29"]


# ----------------------------------------------------------------------
# ReactionState
# ----------------------------------------------------------------------


def test_reaction_state_initialises_its_binding_maps():
    """They default to None and are filled in `__post_init__`; a caller that
    indexes them straight away must not hit a TypeError."""
    state = ReactionState()

    assert state.fg_bound is not None and state.bg_bound is not None
    state.fg_bound[42] += 1
    assert state.fg_bound[42] == 1
    assert state.bg_bound[7] == 0  # defaultdict, not KeyError


def test_reaction_state_defaults_describe_a_dilute_target():
    """1000 target copies against 999000 background -- 0.1% starting material,
    the regime SWGA exists for."""
    state = ReactionState()
    total = state.fg_molecules + state.bg_molecules
    assert state.fg_molecules / total == pytest.approx(0.001, rel=0.01)
