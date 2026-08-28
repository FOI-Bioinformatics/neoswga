"""An absent reagent contributes nothing.

`AdditiveConcentrations._gc_normalization_correction` modelled the
GC-equalising effect of betaine and TMAC with sigmoids that were not anchored at
zero:

    sigmoid(0.0, midpoint=1.5, steepness=1.5) = 0.095
    sigmoid(0.0, midpoint=1.0, steepness=2.0) = 0.119

so a reaction with no isostabilising agent at all was modelled as having ~20% of
its GC dependence removed. An empty `AdditiveConcentrations()` shifted a GC-rich
12-mer's Tm by -1.95 C.

`ReactionConditions` never reached that path with zero additives, so no shipped
Tm carried the error, but `AdditiveConcentrations` is public and its own module
docstring demonstrates calling it directly.

Sigmoids are a reasonable shape for a saturating dose response and a bad one for
a term that must vanish at zero unless the offset is subtracted. The equalisation
model is now a ramp to the concentration each source reports for full
equalisation, which is zero at zero by construction.
"""

import pytest

from neoswga.core.additives import AdditiveConcentrations, equalization_factor
from neoswga.core.reaction_conditions import ReactionConditions

GC_RICH, AT_RICH = 0.9, 0.1


@pytest.mark.parametrize("gc_content", [0.1, 0.3, 0.5, 0.7, 0.9])
@pytest.mark.parametrize("primer_length", [8, 12, 18])
def test_no_additives_means_no_tm_correction(gc_content, primer_length):
    empty = AdditiveConcentrations()

    assert empty.calculate_tm_correction(gc_content, primer_length) == pytest.approx(0.0)


def test_the_equalisation_factor_is_zero_at_zero():
    assert equalization_factor() == 0.0
    assert equalization_factor(betaine=0.0, tmac=0.0) == 0.0


def test_reaction_conditions_agrees_on_both_paths():
    conditions = ReactionConditions(temp=30.0, polymerase="phi29")

    for gc in (0.1, 0.5, 0.9):
        assert conditions.calculate_tm_correction(gc, 12, use_arrhenius=True) == pytest.approx(0.0)
        assert conditions.calculate_tm_correction(gc, 12, use_arrhenius=False) == pytest.approx(0.0)


def test_an_equaliser_still_moves_both_directions():
    """The term has to survive being anchored: it pulls a GC-rich primer down
    and an AT-rich primer up, toward the 50%-GC baseline."""
    betaine = AdditiveConcentrations(betaine_m=2.0)

    assert betaine.calculate_tm_correction(GC_RICH, 12) < 0
    assert betaine.calculate_tm_correction(AT_RICH, 12) > 0
    assert betaine.calculate_tm_correction(0.5, 12) == pytest.approx(
        AdditiveConcentrations(betaine_m=2.0)._betaine_uniform_correction()
    )
