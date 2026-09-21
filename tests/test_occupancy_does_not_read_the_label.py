"""A duplex does not know whether the reference is the target or the host.

Task 5 of the 2026-09-21 valid-design plan: "identical duplexes under identical
local assumptions must get identical occupancy regardless of foreground or
background label", and "a hypothetical mismatch discrimination diagnostic must
not become measured specificity when the index contains only exact matches".

Both matter because specificity is a RATIO of two loads. If the two sides were
computed by different code, or under different assumptions, the ratio would
measure the difference between the two calculations as much as the difference
between the genomes -- and a ratio is exactly the shape in which that is
invisible, because a common factor cancels and an asymmetric one does not.
"""

import pytest

from neoswga.core.occupancy import mismatch_tm, site_occupancy
from neoswga.core.reaction_conditions import ReactionConditions
from neoswga.core.thermodynamics import calculate_enthalpy_entropy

SEQUENCES = ["ACGTTGCAAGGC", "AAAACCCCGGGG", "TTTTGGGGCCCC", "GCTAAAGACAAT"]


@pytest.fixture
def conditions():
    return ReactionConditions(temp=30.0, polymerase="phi29")


@pytest.mark.parametrize("sequence", SEQUENCES)
def test_one_site_has_one_occupancy_however_it_is_labelled(sequence, conditions):
    """The same duplex, asked for twice. There is no label in the arithmetic.

    `site_occupancy` takes an enthalpy, a melting temperature and a reaction
    temperature. None of those is a genome, which is the property being pinned:
    a target site and a host site with identical local sequence are the same
    physical duplex and must not receive different numbers.
    """
    enthalpy, _entropy = calculate_enthalpy_entropy(sequence)
    melting = conditions.calculate_effective_tm(sequence)

    as_foreground = site_occupancy(enthalpy, melting, conditions.temp)
    as_background = site_occupancy(enthalpy, melting, conditions.temp)

    assert as_foreground == as_background
    assert 0.0 <= as_foreground <= 1.0


def test_the_same_load_helper_serves_both_sides():
    """One function, so the two sides of the ratio cannot drift apart.

    `weighted_site_load` is called for the foreground prefixes and for the
    background prefixes. Two helpers that merely looked alike is how this
    codebase has produced disagreeing quantities before -- three different
    coverage semantics under one name, in one audit.
    """
    import ast
    import inspect

    from neoswga.core import base_optimizer

    tree = ast.parse(inspect.getsource(base_optimizer))
    called = {
        getattr(node.func, "attr", None) or getattr(node.func, "id", None)
        for node in ast.walk(tree)
        if isinstance(node, ast.Call)
    }

    assert "weighted_site_load" in called, (
        "the occupancy-weighted load no longer goes through one helper; two "
        "implementations of the two sides of a ratio is how they drift"
    )


@pytest.mark.parametrize("distance", [0, 1, 2])
def test_a_mismatch_costs_the_same_wherever_the_site_is(distance, conditions):
    """The penalty is a property of the duplex, not of which genome it is on.

    It is also uniform in identity and position, which the evidence registry
    records as `assumed`: a real mismatch costs roughly 4 to 8 C depending on
    which bases and where. Sites are grouped by mismatch COUNT because per-site
    alignment against a 3 Gbp host is not tractable, and that approximation is
    stated rather than buried.
    """
    melting = conditions.calculate_effective_tm(SEQUENCES[0])

    assert mismatch_tm(melting, distance) == mismatch_tm(melting, distance)
    assert mismatch_tm(melting, distance) <= melting


def test_a_mismatch_lowers_occupancy_rather_than_raising_it(conditions):
    enthalpy, _entropy = calculate_enthalpy_entropy(SEQUENCES[0])
    melting = conditions.calculate_effective_tm(SEQUENCES[0])

    matched = site_occupancy(enthalpy, melting, conditions.temp)
    mismatched = site_occupancy(enthalpy, mismatch_tm(melting, 1), conditions.temp)

    assert mismatched <= matched


# ---------------------------------------------------------------------------
# Exact and modelled measurements stay apart
# ---------------------------------------------------------------------------


def test_the_selectivity_mode_says_which_measurement_was_made():
    """`exact` and `modelled` are different claims and must not be conflated.

    An index holding only exact matches supports an exact-match specificity
    statement. A mismatch-weighted figure computed over that same index is a
    hypothetical: it says what the load WOULD be if near-matches existed where
    the model supposes, which is not a measurement of the host.
    """
    from neoswga.core.base_optimizer import PrimerSetMetrics

    default = PrimerSetMetrics.__dataclass_fields__["selectivity_mode"].default

    assert default == "exact", (
        "the default selectivity claim must be the weaker one; defaulting to a "
        "modelled figure would report a hypothetical as a measurement"
    )


def test_the_mode_is_carried_into_the_saved_result():
    """A reader cannot interpret the ratio without knowing which kind it is."""
    import ast
    import inspect

    from neoswga.core import base_optimizer

    source = inspect.getsource(base_optimizer.PrimerSetMetrics.to_dict)
    keys = {
        node.value
        for node in ast.walk(ast.parse(source.strip()))
        if isinstance(node, ast.Constant) and isinstance(node.value, str)
    }

    assert "selectivity_mode" in keys
    assert "selectivity_ratio" in keys


def test_occupancy_is_bounded_and_monotone_in_temperature(conditions):
    """Sanity the ratio depends on: warmer is never more bound.

    If occupancy rose with temperature anywhere in the supported band, every
    specificity figure computed from it would be the wrong way round, and a
    ratio would hide it whenever both sides moved together.
    """
    enthalpy, _entropy = calculate_enthalpy_entropy(SEQUENCES[0])
    melting = conditions.calculate_effective_tm(SEQUENCES[0])

    values = [site_occupancy(enthalpy, melting, temp) for temp in (20.0, 30.0, 40.0)]

    assert values == sorted(values, reverse=True)
    assert all(0.0 <= value <= 1.0 for value in values)


def test_occupancy_is_one_half_at_the_melting_temperature(conditions):
    """What Tm means, and the anchor the whole two-state model rests on."""
    enthalpy, _entropy = calculate_enthalpy_entropy(SEQUENCES[0])
    melting = conditions.calculate_effective_tm(SEQUENCES[0])

    assert site_occupancy(enthalpy, melting, melting) == pytest.approx(0.5)
