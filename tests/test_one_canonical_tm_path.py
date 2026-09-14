"""Three ways the same quantity was computed differently.

Audit finding F6. The project has one melting temperature and had three paths to
it, disagreeing in ways no test could see:

- `ReactionConditions.from_additives` forwarded ten of the eleven additive
  fields. A 1 M propanediol input arrived as 0 M.
- `MechanisticModel._calculate_effective_tm` takes a `tm_correction` argument
  and, on its primary path, returns `conditions.calculate_effective_tm(primer)`
  and discards it. Interaction corrections of -1 C and -10 C produced the same
  effective Tm.

The fix is not to make the mechanistic model add its whole correction on top of
the canonical one -- both compute the full additive effect, so that would count
every additive twice. Only the CHANGE an interaction makes belongs to the
mechanistic model; the canonical calculation already carries the rest.
"""

import dataclasses

import pytest

from neoswga.core.additives import AdditiveConcentrations
from neoswga.core.reaction_conditions import ReactionConditions


def _probe_value(field_name):
    """A distinctive value inside every field's own validated range.

    `AdditiveConcentrations` range-checks each field, and the ranges differ by
    two orders of magnitude (TMAC tops out at 0.1 M, betaine at 2.5 M), so one
    shared probe value is rejected for some fields. These are small enough to
    be valid everywhere and non-zero, which is all the round trip needs.
    """
    if field_name.endswith("_ml"):
        return 120.0
    if field_name.endswith("_percent"):
        return 3.0
    if field_name.endswith("_m"):
        return 0.05
    raise AssertionError(f"no probe value for {field_name}")


@pytest.mark.parametrize("field", [f.name for f in dataclasses.fields(AdditiveConcentrations)])
def test_every_additive_field_survives_conversion(field):
    """A field-by-field round trip, so the next omission fails here.

    Parametrised over the dataclass rather than a hand-written list: a new
    additive is covered the day it is added, which is what the propanediol
    omission needed and did not have.
    """
    value = _probe_value(field)
    additives = AdditiveConcentrations(**{field: value})

    conditions = ReactionConditions.from_additives(additives, temp=30.0)

    assert getattr(conditions, field) == pytest.approx(
        value
    ), f"{field} was dropped by ReactionConditions.from_additives"


def test_an_interaction_correction_changes_the_effective_tm():
    """The argument was accepted and ignored on the path that runs."""
    from neoswga.core.mechanistic_model import MechanisticModel

    conditions = ReactionConditions(temp=30.0, polymerase="phi29", dmso_percent=5.0)
    model = MechanisticModel(conditions)
    primer = "ACGTACGTACGT"
    gc = model._primer_gc(primer)

    baseline = model._calculate_effective_tm(primer, gc, 0.0)
    shifted = model._calculate_effective_tm(primer, gc, 0.0, interaction_delta=-5.0)

    assert shifted == pytest.approx(baseline - 5.0)


def test_the_canonical_correction_is_not_applied_twice():
    """Guard the fix against the obvious wrong version.

    Adding the mechanistic model's whole `tm_correction` on top of
    `conditions.calculate_effective_tm` would double every additive. With no
    interaction in play the two must agree exactly.
    """
    from neoswga.core.mechanistic_model import MechanisticModel

    conditions = ReactionConditions(temp=30.0, polymerase="phi29", dmso_percent=10.0)
    model = MechanisticModel(conditions)
    primer = "ACGTACGTACGT"
    gc = model._primer_gc(primer)
    correction = model._calculate_tm_correction(gc, len(primer))

    assert correction != 0.0  # DMSO is in play, so the guard is not vacuous
    assert model._calculate_effective_tm(primer, gc, correction) == pytest.approx(
        conditions.calculate_effective_tm(primer)
    )
