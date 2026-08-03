"""Per-site binding occupancy, so reaction conditions can affect specificity.

Selectivity in this project is a count of exact k-mer matches
(``base_optimizer.py``): every site contributes 1.0 regardless of how stable
that duplex is under the actual conditions, and no additive or temperature
appears anywhere in the arithmetic. Additives move Tm; Tm is not in the formula.
So DMSO cannot change predicted specificity by any amount, however it is set.

This module supplies the missing term. For a duplex melting at ``tm`` under the
configured conditions, at reaction temperature ``temp``, the two-state fraction
bound is

    theta = 1 / (1 + exp( (dH / R) * (1/T - 1/Tm) ))

with both temperatures in Kelvin. At ``T == Tm`` it is 0.5, which is what Tm
means, and it approaches 1 well below and 0 well above.

Nothing here is new thermodynamics. ``dH`` comes from
``thermodynamics.calculate_enthalpy_entropy`` and the additive- and
salt-corrected Tm from ``ReactionConditions.calculate_effective_tm``; occupancy
is layered onto corrections that already existed and were only ever applied
per-primer.

Why it delivers condition-dependent specificity: a foreground site matches
perfectly while a background site typically carries a mismatch or two, so it
melts lower. Raising stringency -- a hotter reaction, or an additive that lowers
Tm -- moves the weaker duplex down the sigmoid faster than the stronger one, and
the ratio between them grows. That ratio is the discrimination an additive buys.

The effect is bounded at both ends, which is worth knowing before tuning a
protocol around it. Far below Tm everything is saturated and there is nothing to
discriminate; far above, nothing is bound and there is nothing to amplify. The
useful window is near Tm, which is the argument for equiphi29 at 42 C over
phi29 at 30 C when specificity is the binding constraint.
"""

import math
from typing import Optional

from neoswga.core.thermodynamics import R

# Gas constant is in cal/(mol*K) while the nearest-neighbour tables report
# enthalpy in kcal/mol, so dH crosses that factor before entering the exponent.
_KCAL_TO_CAL = 1000.0

# Clamp on the exponent. exp() overflows around 709, and any |x| past ~40 has
# already saturated theta to within floating-point noise of 0 or 1, so this
# changes no answer -- it only stops a far-from-Tm site raising OverflowError.
_MAX_EXPONENT = 40.0

# Fallback destabilisation per mismatch, in Celsius, when params.json does not
# say. A single internal mismatch costs roughly 4-8 C depending on identity and
# position; 4 is the conservative end and matches the long-standing (unused)
# default in `utility.mismatch_penalty`.
DEFAULT_MISMATCH_PENALTY_C = 4.0

_ABSOLUTE_ZERO_C = -273.15


def default_mismatch_penalty() -> float:
    """Per-mismatch Tm penalty in Celsius, from configuration.

    ``mismatch_penalty`` has been a params.json key with no reader: three
    occurrences in the package (``utility.py``, and two lines of ``parameter.py``
    declaring it), none of them a use. The one place a penalty is actually
    applied hardcodes 4 because the call site never passes the configured value.
    This is its first consumer, and its intended meaning.
    """
    from neoswga.core import parameter

    value = getattr(parameter, "mismatch_penalty", None)
    if value is None:
        return DEFAULT_MISMATCH_PENALTY_C
    return float(value)


def mismatch_tm(tm: float, n_mismatches: int, penalty: Optional[float] = None) -> float:
    """Melting temperature of a duplex carrying ``n_mismatches``.

    A linear penalty is a simplification: real destabilisation depends on which
    bases are involved and where they sit, and a mismatch near the 3' end costs
    far more than one in the middle. Sites are grouped by mismatch count rather
    than modelled individually because that is what makes background counting
    tractable at all -- the alternative is a per-site alignment against a 3 Gbp
    genome. The linear form is stated here so the approximation is visible
    rather than buried.
    """
    if penalty is None:
        penalty = default_mismatch_penalty()
    return tm - n_mismatches * penalty


def site_occupancy(dh_kcal: float, tm: float, temp: float) -> float:
    """Fraction of the time a site with melting temperature ``tm`` is bound.

    Args:
        dh_kcal: Duplex enthalpy in kcal/mol, negative for hybridisation, as
            returned by ``thermodynamics.calculate_enthalpy_entropy``.
        tm: Melting temperature in Celsius, additive-corrected.
        temp: Reaction temperature in Celsius.

    Returns:
        Occupancy in [0, 1]; exactly 0.5 when ``temp == tm``.
    """
    temp_k = temp - _ABSOLUTE_ZERO_C
    tm_k = tm - _ABSOLUTE_ZERO_C

    # A duplex at or below absolute zero is not a physical state; treating it as
    # fully bound keeps a nonsense input from propagating a nan into a score.
    if temp_k <= 0 or tm_k <= 0:
        return 1.0 if temp_k <= tm_k else 0.0

    # No enthalpy means no transition to model, so the site carries no
    # information about stability. 0.5 is the honest answer, and it keeps a
    # degenerate sequence from dividing by zero.
    if dh_kcal == 0:
        return 0.5

    exponent = (dh_kcal * _KCAL_TO_CAL / R) * (1.0 / temp_k - 1.0 / tm_k)
    exponent = max(-_MAX_EXPONENT, min(_MAX_EXPONENT, exponent))

    return 1.0 / (1.0 + math.exp(exponent))
