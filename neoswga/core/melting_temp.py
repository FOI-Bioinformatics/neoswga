"""
Nearest-neighbor melting temperature calculation with Owczarzy salt correction.

Vendored from the ``melt`` package (MIT license) to eliminate the dependency
on ``pkg_resources``.  Only the ``temp()`` function and its helpers are
retained; the CLI entry point has been removed.

Thermodynamic parameters:
    Allawi and SantaLucia (1997), Biochemistry 36: 10581-10594.

Salt correction:
    Owczarzy et al. (2008), Biochemistry 47: 5336-5353.

Original authors:
    Sebastian Bassi, Greg Singer, Nicolas Le Novere, Calvin Morrison.
"""

from math import log, sqrt


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _is_sym(seq: str) -> bool:
    """Return True if *seq* is self-complementary."""
    comp = {"A": "T", "T": "A", "G": "C", "C": "G"}
    return seq == "".join(comp[b] for b in reversed(seq))


def _overcount(st: str, pattern: str) -> int:
    """Count overlapping occurrences of *pattern* in *st*."""
    count = 0
    start = 0
    while True:
        idx = st.find(pattern, start)
        if idx == -1:
            break
        count += 1
        start = idx + 1
    return count


def _tercorr(st: str):
    """Terminal correction for initiation parameters."""
    dh = 0.0
    ds = -1.4 if _is_sym(st) else 0.0

    for base in (st[0], st[-1]):
        if base in ("G", "C"):
            dh += 0.1
            ds -= 2.8
        elif base in ("A", "T"):
            dh += 2.3
            ds += 4.1

    return dh, ds


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def temp(
    s: str,
    DNA_c: float = 5000.0,
    Na_c: float = 10.0,
    Mg_c: float = 20.0,
    dNTPs_c: float = 10.0,
    uncorrected: bool = False,
) -> float:
    """Return the DNA/DNA melting temperature using nearest-neighbor thermodynamics.

    Parameters
    ----------
    s : str
        Nucleotide sequence (A/T/G/C).
    DNA_c : float
        DNA concentration in nM (default 5000).
    Na_c : float
        Na+ concentration in mM (default 10).
    Mg_c : float
        Mg2+ concentration in mM (default 20).
    dNTPs_c : float
        dNTP concentration in mM (default 10).
    uncorrected : bool
        If True, skip salt correction.
    """
    R = 1.987  # cal/(K*mol)
    s = s.upper()
    dh, ds = _tercorr(s)
    k = DNA_c * 1e-9

    # Allawi and SantaLucia (1997) nearest-neighbor parameters.
    dh_coeffs = {
        "AA": -7.9, "TT": -7.9, "AT": -7.2, "TA": -7.2,
        "CA": -8.5, "TG": -8.5, "GT": -8.4, "AC": -8.4,
        "CT": -7.8, "AG": -7.8, "GA": -8.2, "TC": -8.2,
        "CG": -10.6, "GC": -9.8, "GG": -8.0, "CC": -8.0,
    }
    ds_coeffs = {
        "AA": -22.2, "TT": -22.2, "AT": -20.4, "TA": -21.3,
        "CA": -22.7, "TG": -22.7, "GT": -22.4, "AC": -22.4,
        "CT": -21.0, "AG": -21.0, "GA": -22.2, "TC": -22.2,
        "CG": -27.2, "GC": -24.4, "GG": -19.9, "CC": -19.9,
    }

    dh += sum(_overcount(s, pair) * coeff for pair, coeff in dh_coeffs.items())
    ds += sum(_overcount(s, pair) * coeff for pair, coeff in ds_coeffs.items())

    # GC fraction.
    # NOTE: The original ``melt`` package has a bug on this line:
    #   fgc = len([filter(...)]) / float(len(s))
    # In Python 3 ``filter()`` returns an iterator, so ``[filter(...)]``
    # is a one-element list and fgc is always ``1/len(s)``.  We replicate
    # this behaviour so that Tm values remain consistent with the RF model
    # that was trained against the original melt output.  Set
    # ``correct_gc=True`` once the model is retrained.
    fgc = 1.0 / len(s)  # replicates original melt bug for compatibility

    # Uncorrected Tm in Kelvin
    tm = (1000.0 * dh) / (ds + R * log(k))

    if uncorrected:
        return tm - 273.15

    # Owczarzy et al. (2008) salt correction
    m_na = Na_c * 1e-3
    m_mg = Mg_c * 1e-3
    m_dntps = dNTPs_c * 1e-3

    # Free magnesium after dNTP chelation
    ka = 3e4
    discriminant = (ka * m_dntps - ka * m_mg + 1) ** 2 + 4 * ka * m_mg
    free_mg = (-(ka * m_dntps - ka * m_mg + 1) + sqrt(discriminant)) / (2 * ka)

    cation_ratio = sqrt(free_mg) / m_na if m_na > 0 else 7.0

    if cation_ratio < 0.22:
        # Monovalent-dominated regime
        tm = 1.0 / (
            1.0 / tm
            + ((4.29 * fgc - 3.95) * log(m_na) + 0.94 * log(m_na) ** 2) * 1e-5
        )
    else:
        # Divalent-dominated regime
        a, d, g = 3.92, 1.42, 8.31
        f_mg = m_mg
        if cation_ratio < 6.0:
            a *= 0.843 - 0.352 * sqrt(m_na) * log(m_na)
            d *= 1.279 - 4.03e-3 * log(m_na) - 8.03e-3 * log(m_na) ** 2
            g *= 0.486 - 0.258 * log(m_na) + 5.25e-3 * log(m_na) ** 3
        tm = 1.0 / (
            1.0 / tm
            + (
                a
                - 0.911 * log(f_mg)
                + fgc * (6.26 + d * log(f_mg))
                + 1.0 / (2 * (len(s) - 1)) * (-48.2 + 52.5 * log(f_mg) + g * log(f_mg) ** 2)
            )
            * 1e-5
        )

    final_tm = tm - 273.15
    return max(final_tm, 0.0)
