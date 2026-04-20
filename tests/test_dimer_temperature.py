"""Dimer severity must respond to reaction temperature.

A primer-pair dimer stable at 30 C (phi29) can melt at 42-45 C (equiphi29).
The StructurePrediction machinery takes ReactionConditions.temp and is
supposed to reflect this in the `severity` it returns. Without an explicit
regression test, future refactors could silently decouple severity from
temperature — optimizers would then treat every polymerase identically at
the dimer-scoring step, which is wrong.
"""

import pytest

from neoswga.core.reaction_conditions import ReactionConditions
from neoswga.core.secondary_structure import StructurePrediction


# A deliberately complementary pair; forms a strong dimer at low temp.
HOT_DIMER_SEQ_1 = "ATCGATCGAT"
HOT_DIMER_SEQ_2 = "ATCGATCGAT"  # palindrome-like — binds to its own reverse complement


def _severity(temp: float, polymerase: str = "phi29", betaine_m: float = 0.0):
    cond = ReactionConditions(
        temp=temp, polymerase=polymerase, betaine_m=betaine_m,
    )
    pred = StructurePrediction(cond)
    result = pred.predict_heterodimer(HOT_DIMER_SEQ_1, HOT_DIMER_SEQ_2)
    return result["severity"], result.get("tm", 0.0), result["forms_dimer"]


def test_higher_temperature_reduces_dimer_severity():
    """Raising temp from 30 C to 45 C should lower severity of the same pair."""
    sev_30, tm_30, stable_30 = _severity(temp=30.0, polymerase="phi29")
    sev_45, tm_45, stable_45 = _severity(temp=45.0, polymerase="equiphi29")

    # The pair is engineered to form at 30 C.
    assert stable_30, (
        f"Test precondition broken: chosen pair should form dimer at 30 C "
        f"(severity={sev_30}, tm={tm_30:.1f})"
    )

    # The two severities must differ (not frozen) and severity at 45 C
    # must be no greater than at 30 C. A strict drop is the happy path;
    # we allow equality only when both drop to 0 (neither forms a dimer).
    assert sev_45 <= sev_30 + 1e-6, (
        f"Severity should not increase as temperature rises: "
        f"30 C -> {sev_30:.3f}, 45 C -> {sev_45:.3f}"
    )
    # Require at least a 20% relative drop OR the 45 C pair must be
    # reported as no longer forming a stable dimer.
    if stable_45 and sev_45 > 0:
        rel_drop = (sev_30 - sev_45) / sev_30 if sev_30 else 0.0
        assert rel_drop >= 0.20, (
            f"Raising temp from 30 C to 45 C should drop severity by >= 20%; "
            f"got {rel_drop*100:.1f}% (30 C sev={sev_30:.3f}, 45 C sev={sev_45:.3f})"
        )


def test_tm_estimate_is_temperature_aware():
    """The Tm estimate itself is a property of the dimer, not the reaction
    temperature — but the "stable_at_temp" flag should flip when reaction
    temp crosses dimer Tm."""
    sev_30, tm_30, stable_30 = _severity(temp=30.0, polymerase="phi29")
    # If tm_30 is known, pick a temperature well above it and expect stable_at_temp=False
    if tm_30 > 35.0:
        hot_temp = tm_30 + 5.0
        sev_hot, tm_hot, stable_hot = _severity(
            temp=min(hot_temp, 55.0),  # cap at bst range to stay in polymerase validation
            polymerase="bst" if hot_temp > 50 else "equiphi29",
        )
        assert not stable_hot or sev_hot < sev_30, (
            f"Above dimer Tm, pair should be reported unstable or with lower severity"
        )


def test_betaine_does_not_increase_severity_at_constant_temp():
    """Betaine lowers Tm, so at fixed reaction temperature it should not
    push a pair into higher severity. This is a weak but useful guardrail."""
    sev_plain, _, _ = _severity(temp=30.0, polymerase="phi29", betaine_m=0.0)
    sev_betaine, _, _ = _severity(temp=30.0, polymerase="phi29", betaine_m=1.5)
    # Tolerate tiny floating noise
    assert sev_betaine <= sev_plain + 1e-3, (
        f"Betaine 1.5 M should not increase dimer severity at fixed temp: "
        f"plain={sev_plain:.3f}, betaine={sev_betaine:.3f}"
    )
