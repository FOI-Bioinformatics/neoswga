"""A mismatch model in which where the mismatch sits changes the answer.

The shipped model is uniform: ``occupancy.mismatch_tm`` subtracts one number
per mismatch, independent of identity and position, and its registry record
(``mismatch_penalty``) has status ``assumed``. This module is the alternative,
reachable through the ``mismatch_model`` key and OFF by default.

**It is an extrapolation, and the code says so.** The duplex term is built on
``thermodynamics.DELTA_G_MISMATCH``, the SantaLucia (1998) nearest-neighbour
internal-mismatch set. Those values were measured at 37 C in 1 M Na+ on duplex
DNA; SWGA runs near 30 C in a magnesium buffer on 8-12mers. The table is a
measurement, this use of it is not, and the registry records the two separately
(``mismatch_duplex_delta_g``). Nothing here is a measurement of mismatch
discrimination, and a load it produces is reported with its own mode string so
it can never be read as the uniform one.

**Two effects, kept apart.** They have different evidence and are switched
separately, so a later measurement can overturn one without the other.

1. *Duplex stability*, graded. Position enters through the flanking stacks
   rather than through a position index: a mismatch next to the end has one
   neighbouring stack inside the duplex, one in the middle has two, and the
   identity of both flanks is in the table. This is the term
   ``position-dependent`` turns on.
2. *3'-terminal extension*, close to binary and independent of (1). A
   strand-displacing polymerase extends poorly from a mismatched 3' terminus
   even when the duplex is stable, so such a site is not a weak site but
   largely a non-site. This is the term ``position-dependent-3prime`` adds. Its
   evidence is thin -- the Innis (1988) citation in
   ``three_prime_stability.py`` is a docstring, and no measurement for a
   strand-displacing polymerase was found -- so it is a separate value of the
   key rather than part of the first.

**Its own walk, deliberately.**
``thermodynamics.compute_free_energy_for_two_strings`` abandons the walk once
the running total passes ``penalty * 10``, which is right for the dimer screen
it was written for and wrong for a whole-duplex dG, where an early stop
silently truncates the sum. It also applies no nearest-neighbour symmetry, so
every mismatch's 3'-side doublet falls through to the flat ``penalty``. Both are
fixed here and neither is changed there: that function's behaviour for its one
existing caller (``rf_preprocessing``, the retired random-forest feature path)
is untouched.

**The symmetry.** ``DELTA_G_MISMATCH`` holds the 64 doublets whose 5' position
is Watson-Crick paired. A nearest-neighbour doublet is unchanged by rotating
the duplex 180 degrees, so the doublet ``5'-CG-3'/3'-TC-5'`` -- absent as
written -- is the tabulated ``5'-CT-3'/3'-GC-5'``. Applying that rotation
covers every doublet with at least one Watson-Crick position, which is every
doublet a single mismatch can produce. A doublet with two mismatches is in
neither form, because the source does not parameterise tandem mismatches, and
falls back to ``penalty``; at the default depth of one mismatch that case does
not arise.

**The table's domain is INTERNAL positions, and a terminus is outside it.**
``DELTA_G_MISMATCH`` is an internal single-mismatch set. The first and last
base of a duplex have one flanking stack rather than two, and the applicable
parameters there are dangling-end and terminal-mismatch sets (Bommarito et al.
2000) which this repository does not carry. Applying the internal set at a
terminus, while also dropping that end's initiation correction, measured as a
net STABILISING result for 1 of 36 one-mismatch neighbours of a 12-mer: a
mismatched host site weighted as well as a perfect one. So a neighbour with a
mismatch at either terminus does not use the table at all and takes the uniform
assumed penalty instead -- the same number the shipped model charges it. See
``mismatch_delta_tm``. The 3'-terminal extension factor of the ``-3prime``
variant still applies on top, because that is an extension effect and not a
duplex one.

**From dG to occupancy.** ``occupancy.site_occupancy`` takes an enthalpy and a
melting temperature, not a free energy, and the melting temperature it is given
already carries the salt and additive corrections. So the duplex term enters as
a Tm SHIFT on the perfect duplex rather than as a replacement for it:

    dTm = Tm * ddG(37 C) / dH

which is the first-order perturbation of ``Tm = dH / (dS + R ln C)`` under
``dH -> dH + ddH`` and ``dS -> dS + ddS``, evaluated where the mismatch's cost
is read at the table's reference temperature rather than at Tm. It replaces
``- distance * penalty`` and nothing else, which is what makes the uniform model
the special case: force dTm to that constant and the two loads agree.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

from neoswga.core.mismatch_sites import NeighbourReading, NeighbourSite, neighbour_sites
from neoswga.core.thermodynamics import (
    DELTA_G_MISMATCH,
    NN_INIT_CORRECTIONS,
)
from neoswga.core.utility import complement as _complement

__all__ = [
    "MISMATCH_MODELS",
    "UNIFORM",
    "POSITION_DEPENDENT",
    "POSITION_DEPENDENT_THREE_PRIME",
    "THREE_PRIME_EXTENSION_FACTOR",
    "THREE_PRIME_WINDOW",
    "TABLE_REFERENCE_BUFFER",
    "TABLE_REFERENCE_K",
    "duplex_delta_g",
    "extrapolation_notice",
    "log_extrapolation_once",
    "mismatch_delta_tm",
    "position_dependent_site_load",
    "resolve_mismatch_model",
    "site_load_mode",
]

#: Today's behaviour, and the default.
UNIFORM = "uniform"
#: The duplex term alone.
POSITION_DEPENDENT = "position-dependent"
#: The duplex term plus the 3'-terminal extension rule.
POSITION_DEPENDENT_THREE_PRIME = "position-dependent-3prime"

#: Every accepted value of the ``mismatch_model`` key, weakest claim first.
MISMATCH_MODELS = (UNIFORM, POSITION_DEPENDENT, POSITION_DEPENDENT_THREE_PRIME)

#: Fallback for a doublet the source does not parameterise, in kcal/mol. The
#: same number and the same meaning as the ``penalty`` argument of
#: ``thermodynamics.compute_free_energy_for_two_strings``.
UNPARAMETERISED_DOUBLET_PENALTY = 4.0

#: Fallback uniform Tm penalty, in Celsius, for a position the table does not
#: cover. Only a default: every path that has the run's resolved
#: ``mismatch_penalty`` passes it instead.
DEFAULT_UNIFORM_PENALTY_C = 4.0

#: The reference temperature of ``DELTA_G_MISMATCH``, in Kelvin. 37 C.
TABLE_REFERENCE_K = 310.15

#: The buffer ``DELTA_G_MISMATCH`` was measured in. Named so the extrapolation
#: notice can state what it is extrapolating FROM rather than only how far.
TABLE_REFERENCE_BUFFER = "1 M NaCl"

#: Palindrome symmetry correction, in kcal/mol, as the existing walk applies it.
_SYMMETRY_CORRECTION = 0.43

_ABSOLUTE_ZERO_C = -273.15

#: How close to the 3' end a mismatch must sit for the extension rule to apply.
#: 1 is the terminal base alone, which is the default and the only span the
#: cited observation is about.
THREE_PRIME_WINDOW = 1

#: What is left of a site whose 3' terminus is mismatched. Not a measurement:
#: it states "largely a non-site" as a number so the claim is visible and can
#: be overturned. Its registry record is ``three_prime_mismatch_extension``
#: with status ``assumed``.
THREE_PRIME_EXTENSION_FACTOR = 0.05


def _doublet(x_pair: str, y_pair: str, penalty: float) -> float:
    """One nearest-neighbour doublet's dG, with the rotation fallback.

    ``x_pair`` is two bases 5' to 3'; ``y_pair`` is the two bases facing them,
    3' to 5'. Rotating the duplex 180 degrees exchanges the strands and
    reverses both, which is the same physical stack, so the rotated key is
    looked up when the written one is absent.
    """
    direct = DELTA_G_MISMATCH.get(f"{x_pair}/{y_pair}")
    if direct is not None:
        return float(direct)
    rotated = DELTA_G_MISMATCH.get(f"{y_pair[::-1]}/{x_pair[::-1]}")
    if rotated is not None:
        return float(rotated)
    return float(penalty)


def duplex_delta_g(
    primer: str,
    neighbour: str,
    *,
    penalty: float = UNPARAMETERISED_DOUBLET_PENALTY,
) -> float:
    """Free energy of the duplex ``primer`` forms with the site ``neighbour``.

    ``neighbour`` is a genomic k-mer written in the primer's own orientation
    and index order, as ``mismatch_sites`` reports it, so the strand the primer
    actually anneals to is its complement read in the same order. With
    ``neighbour == primer`` this is the perfect duplex, and it equals what
    ``thermodynamics.compute_free_energy_for_two_strings(primer,
    complement(primer))`` returns: a perfect duplex reaches neither the early
    stop nor the rotation.

    Returns dG in kcal/mol at the table's reference state, more negative being
    more stable. The whole walk runs; there is no early stop.
    """
    x = primer.upper()
    y = _complement(neighbour.upper())
    if not x or not y:
        return 0.0

    delta_g = 0.0

    # Terminal corrections, applied at an end only when it is Watson-Crick
    # paired. A mismatched terminus takes none, which is the existing walk's
    # behaviour and is the conservative direction: the correction is positive.
    for base_x, base_y in ((x[0], y[0]), (x[-1], y[-1])):
        if _complement(base_x) == base_y:
            correction = NN_INIT_CORRECTIONS.get(f"{base_x}/{base_y}")
            if correction is not None:
                delta_g += float(correction)

    limit = min(len(x), len(y))
    for i in range(1, limit):
        delta_g += _doublet(x[i - 1 : i + 1], y[i - 1 : i + 1], penalty)

    # Reproduced from `thermodynamics._compute_free_energy_for_two_strings_impl`
    # exactly, which is what makes the perfect-duplex equality hold. Note that
    # this condition is a string palindrome and not self-complementarity: for a
    # perfect duplex it reduces to `x == reverse(x)`, so it fires for ACGTTGCA
    # and not for ACGCGT. That is inherited, not verified here, and it is not
    # introduced by this module.
    if _complement(x) == y[::-1]:
        delta_g += _SYMMETRY_CORRECTION

    return delta_g


def mismatch_delta_tm(
    primer: str,
    reading: NeighbourReading,
    *,
    dh_kcal: float,
    tm: float,
    uniform_penalty: float = DEFAULT_UNIFORM_PENALTY_C,
    penalty: float = UNPARAMETERISED_DOUBLET_PENALTY,
) -> float:
    """Melting-temperature shift of one mismatched duplex, in Celsius.

    Negative, destabilising. First order in the nearest-neighbour cost of the
    mismatch, read at the table's 37 C reference rather than at Tm; see the
    module docstring for the derivation and for what that approximation is.

    **A mismatch at either terminus of the duplex does not use the table.**
    ``DELTA_G_MISMATCH`` is an INTERNAL single-mismatch set and a terminal
    position is outside its domain: that position has one flanking stack rather
    than two, and dangling-end and terminal-mismatch parameters (Bommarito et
    al. 2000) are the applicable set, which this repository does not carry.
    Applying the internal set there, while also dropping that end's initiation
    correction, measured as a NET STABILISING result for 1 of 36 neighbours of
    a 12-mer -- a mismatched host site weighted as well as a perfect one. So a
    neighbour with any mismatch at a terminus falls back to the uniform assumed
    penalty, ``- distance * uniform_penalty``, which is what the shipped model
    charges it. The 3'-terminal extension factor of the ``-3prime`` variant
    still applies on top, because that is a separate effect and not a duplex
    one.

    Returns 0.0 for a perfect duplex, so the exact-match site is weighted by
    the same Tm the uniform path gives it.

    A duplex with no enthalpy has no transition to perturb and no shift can be
    derived from it, so the shift is 0.0 rather than an arbitrary number; the
    occupancy of such a site is already reported as 0.5 by
    ``occupancy.site_occupancy`` for the same reason.
    """
    if not reading.pairs:
        return 0.0

    if reading.is_at_a_duplex_terminus:
        return -float(len(reading.pairs)) * float(uniform_penalty)

    if not dh_kcal:
        return 0.0

    perfect = duplex_delta_g(primer, primer, penalty=penalty)
    mismatched = duplex_delta_g(primer, reading.neighbour, penalty=penalty)
    tm_k = tm - _ABSOLUTE_ZERO_C
    if tm_k <= 0:
        return 0.0
    shift = tm_k * (mismatched - perfect) / dh_kcal
    # Kept as a guard, not as a model. With the terminal positions handled
    # above this is not expected to bind at all, and
    # `test_the_clamp_never_binds_at_an_internal_position` is what makes that a
    # measured statement rather than a hope. If it ever does bind, a mismatch
    # would otherwise raise a site's occupancy above its perfect-match value,
    # which is a mismatch turning a rejection into a pass.
    return min(0.0, shift)


def _extension_factor(reading: NeighbourReading, window: int) -> float:
    """What is left of a site after the 3'-terminal extension rule.

    1.0 unless a mismatch falls inside ``window`` bases of the 3' end. At the
    default window of 1 that is exactly `NeighbourReading.is_three_prime_terminal`,
    which is asked rather than re-derived so the two cannot disagree about what
    "at the 3' terminus" means.
    """
    if not reading.offsets_from_three_prime:
        return 1.0
    span = max(1, int(window))
    terminal = (
        reading.is_three_prime_terminal
        if span == 1
        else min(reading.offsets_from_three_prime) < span
    )
    return THREE_PRIME_EXTENSION_FACTOR if terminal else 1.0


def _group_weight(
    primer: str,
    group: NeighbourSite,
    *,
    dh_kcal: float,
    tm: float,
    temp: float,
    penalty: float,
    duplex_penalty: float,
    three_prime_rule: bool,
    three_prime_window: int,
    uniform_delta_tm: bool,
) -> float:
    """Occupancy of one site group, per site.

    A group has more than one reading only for a near-palindromic primer, where
    a neighbour and its reverse complement are both neighbours and a site
    stored under that canonical form could be either duplex. The more stable
    reading is taken: a site binds in whichever orientation binds better, and
    the alternative is to pick by sort order, which is arbitrary.
    """
    from neoswga.core.occupancy import site_occupancy

    best = 0.0
    for reading in group.readings:
        if uniform_delta_tm:
            shift = -float(group.distance) * float(penalty)
        else:
            shift = mismatch_delta_tm(
                primer,
                reading,
                dh_kcal=dh_kcal,
                tm=tm,
                uniform_penalty=penalty,
                penalty=duplex_penalty,
            )
        weight = site_occupancy(dh_kcal, tm + shift, temp)
        if three_prime_rule:
            weight *= _extension_factor(reading, three_prime_window)
        best = max(best, weight)
    return best


def position_dependent_site_load(
    primers: Sequence[str],
    prefixes: Sequence[str],
    conditions: Any,
    max_mismatches: int = 1,
    penalty: float | None = None,
    *,
    three_prime_rule: bool = False,
    three_prime_window: int = THREE_PRIME_WINDOW,
    duplex_penalty: float = UNPARAMETERISED_DOUBLET_PENALTY,
    uniform_delta_tm: bool = False,
) -> float:
    """Occupancy-weighted binding load with a per-neighbour mismatch cost.

        load = SUM_primers SUM_groups count * theta(dH, Tm + dTm(group), T)

    The same quantity ``occupancy.weighted_site_load`` computes, differing only
    in where ``dTm`` comes from. ``uniform_delta_tm=True`` forces it back to
    ``- distance * penalty``, which is the reduction the two paths are checked
    against: they then differ only in the order the sum is taken.

    Args:
        primers: The panel.
        prefixes: K-mer table prefixes, summed over.
        conditions: Reaction conditions; supplies ``temp`` and the
            additive-corrected Tm.
        max_mismatches: Highest mismatch distance counted.
        penalty: Uniform per-mismatch Tm penalty, used only by the reduction
            and for a caller that wants the run's resolved value threaded
            through. Resolved from configuration when None, as the uniform path
            resolves it.
        three_prime_rule: Apply the 3'-terminal extension factor. Separate from
            the duplex term on purpose; its evidence is weaker.
        three_prime_window: How many 3'-terminal bases the rule covers.
        duplex_penalty: Fallback for a doublet the source does not
            parameterise, in kcal/mol.
        uniform_delta_tm: Weight every mismatch alike. For the reduction test.

    Raises:
        FileNotFoundError, OSError: when a k-mer table is missing, exactly as
            ``occupancy.weighted_site_load`` raises, so a caller falls back to
            exact counting deliberately rather than receiving a quietly
            different number.
    """
    from neoswga.core.mismatch_counts import load_kmer_counts
    from neoswga.core.occupancy import default_mismatch_penalty
    from neoswga.core.thermodynamics import calculate_enthalpy_entropy

    if penalty is None:
        penalty = default_mismatch_penalty()

    temp = conditions.temp
    total = 0.0
    tables_by_k: dict[int, list[dict[str, int]]] = {}
    for primer in primers:
        k = len(primer)
        if k not in tables_by_k:
            tables_by_k[k] = [load_kmer_counts(prefix, k) for prefix in prefixes]
        dh, _ds = calculate_enthalpy_entropy(primer)
        tm = conditions.calculate_effective_tm(primer)
        for group in neighbour_sites(primer, prefixes, max_mismatches, tables=tables_by_k[k]):
            if not group.count:
                continue
            total += group.count * _group_weight(
                primer,
                group,
                dh_kcal=dh,
                tm=tm,
                temp=temp,
                penalty=penalty,
                duplex_penalty=duplex_penalty,
                three_prime_rule=three_prime_rule,
                three_prime_window=three_prime_window,
                uniform_delta_tm=uniform_delta_tm,
            )
    return total


def resolve_mismatch_model(value: Any = None) -> str:
    """Which mismatch model to apply.

    ``value`` wins when given, so a caller that already resolved the run's
    configuration does not re-read a module global per evaluation -- the defect
    ``occupancy.weighted_site_load``'s ``penalty`` argument exists for. With no
    argument the ``mismatch_model`` key is read, and an absent key is
    ``uniform``: the shipped behaviour, unchanged.

    An unrecognised value raises rather than falling back to ``uniform``.
    Silently running the old model for a user who asked for the new one is the
    inert-option defect of Known Issue 8 with an extra step.
    """
    from neoswga.core.exceptions import InvalidDesignRequest

    if value is None:
        from neoswga.core import parameter

        value = getattr(parameter, "mismatch_model", None)
    if value is None:
        return UNIFORM
    text = str(value)
    if text not in MISMATCH_MODELS:
        raise InvalidDesignRequest(
            "mismatch_model",
            f"must be one of {', '.join(MISMATCH_MODELS)}",
            text,
        )
    return text


def extrapolation_notice(model: str | None, temp: float | None) -> str:
    """How far out of its measured domain this model is being used, as prose.

    Empty under ``uniform``, and under no model at all: the uniform penalty has
    its own record and makes no claim about a temperature.

    No refusal threshold exists, and none is invented here. The registry's rule
    is that a range is enforced where one is recorded and nowhere else, and
    neither new record states one -- so a design at 63 C runs. What this
    function does is make the distance visible: the degrees between the
    reaction and ``TABLE_REFERENCE_K``, and the buffer the table was measured
    in against the one the reaction uses. It is carried in the assessment's
    evidence, written into the run manifest, and logged once per run.
    """
    if model is None or model == UNIFORM:
        return ""

    reference_c = TABLE_REFERENCE_K + _ABSOLUTE_ZERO_C
    buffer_part = (
        f"measured in {TABLE_REFERENCE_BUFFER}, applied in a magnesium "
        f"isothermal buffer with additives"
    )
    if temp is None:
        return (
            f"mismatch free energies extrapolated from {reference_c:.0f} C to a "
            f"reaction whose temperature is not recorded; {buffer_part}"
        )
    distance = float(temp) - reference_c
    return (
        f"mismatch free energies extrapolated {abs(distance):.1f} C "
        f"{'above' if distance > 0 else 'below'} their {reference_c:.0f} C "
        f"reference (reaction {float(temp):.1f} C); {buffer_part}"
    )


#: Runs this process has already logged the notice for, so a search that
#: evaluates thousands of panels logs once rather than once per evaluation.
_LOGGED: set[str] = set()


def log_extrapolation_once(model: str | None, temp: float | None) -> str:
    """Warn once per (model, temperature) that the model is out of domain.

    Returns the notice, so a caller that also has to record it does not build
    it twice. Warning level, not info: a number produced outside the domain its
    source covers is exactly what this project's rules say must not pass
    silently.
    """
    notice = extrapolation_notice(model, temp)
    if not notice:
        return ""
    key = f"{model}@{temp}"
    if key not in _LOGGED:
        _LOGGED.add(key)
        import logging

        logging.getLogger(__name__).warning(
            f"mismatch_model={model}: {notice}. No threshold refuses this; see "
            f"the mismatch_duplex_delta_g record in "
            f"core/registry/model_evidence.json."
        )
    return notice


def site_load_mode(model: str | None) -> str:
    """The ``selectivity_mode`` string a load under ``model`` is reported as.

    A position-dependent load must never be reported as the uniform one. The
    uniform path's own string stays ``occupancy``, so no existing output moves,
    and so does an unconfigured one: None is the shipped model, not an unknown.
    """
    if model is None or model == UNIFORM:
        return "occupancy"
    return f"occupancy-{model}"
