"""Two equalising additives cannot each equalise the whole GC deviation.

Betaine and TMAC both work by moving a duplex's stability toward what it would
be at 50% GC. `_gc_adjustment` models that as

    -2 * primer_length * (gc - 0.5) * equalization

where `equalization` is a sigmoid in concentration running 0 to 1. The `2 * L`
scale is the Wallace slope, so the bracket is the entire GC-dependent part of
the melting temperature: at `equalization = 1` the primer melts as if it were
50% GC, and there is nothing left to equalise.

`calculate_total_correction` summed the per-additive corrections, so betaine and
TMAC each contributed a full-magnitude term and the two added. Measured at
GC 0.75 on a 12-mer, with betaine 2.5 M and TMAC 0.1 M, the combined GC term was
-5.792 C against a physical ceiling of -6.0 -- and it stayed inside that ceiling
only because the concentration caps happen to keep the two sigmoids summing to
just under 1. Raising either cap drives the model past the 50%-GC baseline and
out the other side, which is not a thing an isostabilising agent does.

The legacy path (`ReactionConditions._calculate_gc_normalization`) already
combined them as `1 - (1-b)(1-t)`: the fraction of the deviation *not* removed by
either. That is the correct composition and it is what both paths now use.

Urea is deliberately excluded. Its `gc_preference` is a preferential
destabilisation of GC pairs, not equalisation toward a midpoint, so two such
agents genuinely do add.
"""

import pytest

from neoswga.core.additives import ArrheniusTmCorrector
from neoswga.core.mechanistic_params import ADDITIVE_TM_PARAMS

PRIMER_LENGTH = 12


def _gc_term(corrector, additives, gc_content):
    """The GC-dependent part: total at this GC minus total at 50% GC."""
    at_gc = corrector.calculate_total_correction(additives, gc_content, PRIMER_LENGTH)
    at_half = corrector.calculate_total_correction(additives, 0.5, PRIMER_LENGTH)
    return at_gc - at_half


@pytest.fixture
def corrector():
    return ArrheniusTmCorrector(37.0)


@pytest.mark.parametrize("gc_content", [0.25, 0.35, 0.65, 0.75, 0.9])
def test_the_gc_term_never_exceeds_the_wallace_bound(corrector, gc_content):
    """Full equalisation removes the whole GC deviation and no more."""
    ceiling = 2.0 * PRIMER_LENGTH * abs(gc_content - 0.5)
    saturating = {
        "betaine": ADDITIVE_TM_PARAMS["betaine"]["max_concentration"],
        "tmac": ADDITIVE_TM_PARAMS["tmac"]["max_concentration"],
    }

    assert abs(_gc_term(corrector, saturating, gc_content)) <= ceiling + 1e-9


def test_two_equalisers_compose_rather_than_add(corrector):
    """The regression. Summing the two terms overstates the effect; combining
    the fractions they each remove does not."""
    betaine, tmac = 2.0, 0.05

    both = _gc_term(corrector, {"betaine": betaine, "tmac": tmac}, 0.75)
    b_only = _gc_term(corrector, {"betaine": betaine}, 0.75)
    t_only = _gc_term(corrector, {"tmac": tmac}, 0.75)

    assert abs(both) < abs(b_only) + abs(t_only), "the two GC terms are being summed"
    assert abs(both) >= max(abs(b_only), abs(t_only)), "adding an equaliser must not undo one"

    # 1 - (1-b)(1-t), read back through the shared Wallace scale.
    ceiling = 2.0 * PRIMER_LENGTH * 0.25
    fb, ft = abs(b_only) / ceiling, abs(t_only) / ceiling
    assert abs(both) / ceiling == pytest.approx(1 - (1 - fb) * (1 - ft), abs=1e-6)


def test_a_single_equaliser_follows_the_measured_dose_response():
    """The dose response is a ramp to the concentration where equalisation is
    complete, which is the quantity the primary sources report: 5.2 M betaine
    (Rees 1993), 3.0 M TMAC (Melchior & von Hippel 1973).

    The Arrhenius path used a sigmoid centred at `gc_equalization_conc / 3`
    instead -- 1.73 M for betaine, with no source for either the centre or the
    steepness. At the 0.5-1.5 M betaine protocols actually use it returns 43%
    more equalisation than the ramp, so the two paths disagreed by 1.8 C on the
    same reaction. Recorded as known-issues #8 before this.
    """
    from neoswga.core.additives import FULL_EQUALIZATION_M, equalization_factor

    assert equalization_factor(betaine=1.5) == pytest.approx(1.5 / 5.2)
    assert equalization_factor(tmac=0.05) == pytest.approx(0.05 / 3.0)
    assert FULL_EQUALIZATION_M == {"betaine": 5.2, "tmac": 3.0}

    # Saturating and beyond both mean "fully equalised", never more.
    assert equalization_factor(betaine=5.2) == pytest.approx(1.0)
    assert equalization_factor(betaine=100.0) == pytest.approx(1.0)


def test_the_uniform_terms_still_add(corrector):
    """Only the equalisation composes. The concentration-linear part of each
    additive is genuinely additive and must stay so."""
    at_half = corrector.calculate_total_correction(
        {"betaine": 1.0, "dmso": 5.0}, 0.5, PRIMER_LENGTH
    )
    betaine_only = corrector.calculate_total_correction({"betaine": 1.0}, 0.5, PRIMER_LENGTH)
    dmso_only = corrector.calculate_total_correction({"dmso": 5.0}, 0.5, PRIMER_LENGTH)

    assert at_half == pytest.approx(betaine_only + dmso_only, abs=1e-9)


def test_ureas_gc_preference_is_not_treated_as_equalisation(corrector):
    """It destabilises GC pairs rather than moving toward a midpoint, so it is
    not part of the composition and it still adds."""
    both = corrector.calculate_total_correction({"urea": 1.0, "betaine": 1.0}, 0.75, PRIMER_LENGTH)
    urea_only = corrector.calculate_total_correction({"urea": 1.0}, 0.75, PRIMER_LENGTH)
    betaine_only = corrector.calculate_total_correction({"betaine": 1.0}, 0.75, PRIMER_LENGTH)

    assert both == pytest.approx(urea_only + betaine_only, abs=1e-9)


def test_the_two_gc_models_agree(corrector):
    """`ReactionConditions` has an Arrhenius path and a legacy path. They
    disagreed on both the composition and the betaine midpoint, so the same
    reaction gave two different Tm values depending on which was taken."""
    from neoswga.core.reaction_conditions import ReactionConditions

    conditions = ReactionConditions(
        temp=37.0, polymerase="klenow", betaine_m=1.5, tmac_m=0.05, mg_conc=10.0
    )
    arrhenius = conditions.calculate_tm_correction(0.75, PRIMER_LENGTH, use_arrhenius=True)
    legacy = conditions.calculate_tm_correction(0.75, PRIMER_LENGTH, use_arrhenius=False)

    assert arrhenius == pytest.approx(legacy, abs=0.5)
