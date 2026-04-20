"""Numeric regression tests for scientific constants used in thermodynamic
and mechanistic calculations.

These lock each published value to primary literature. Someone changing a
value accidentally or intentionally must update both the constant and this
test file (plus `docs/SCIENCE_CITATIONS.md`). For empirically-tuned values
we assert a reasonable range rather than an exact number, so the test
catches gross drift without requiring lockstep updates for small tuning.
"""

import math

import pytest


# ----------------------------------------------------------------------
# SantaLucia 1998 NN parameters (Table 1, PNAS 95:1460-1465)
# ----------------------------------------------------------------------

# Keys are "XY/X'Y'" where X'Y' is the reverse-complement (SantaLucia 1998 Table 1).
SANTALUCIA_NN_DH = {
    "AA/TT": -7.9, "AT/TA": -7.2, "TA/AT": -7.2,
    "CA/GT": -8.5, "GT/CA": -8.4, "CT/GA": -7.8, "GA/CT": -8.2,
    "CG/GC": -10.6, "GC/CG": -9.8, "GG/CC": -8.0,
}
SANTALUCIA_NN_DS = {
    "AA/TT": -22.2, "AT/TA": -20.4, "TA/AT": -21.3,
    "CA/GT": -22.7, "GT/CA": -22.4, "CT/GA": -21.0, "GA/CT": -22.2,
    "CG/GC": -27.2, "GC/CG": -24.4, "GG/CC": -19.9,
}


def test_santalucia_nn_enthalpy_matches_table_1():
    from neoswga.core import thermodynamics as td
    nn = td.ENTHALPY_NN
    for pair, expected in SANTALUCIA_NN_DH.items():
        assert abs(nn[pair] - expected) < 1e-6, (
            f"NN dH[{pair}] drifted from SantaLucia 1998 Table 1: "
            f"got {nn[pair]}, expected {expected}"
        )


def test_santalucia_nn_entropy_matches_table_1():
    from neoswga.core import thermodynamics as td
    nn = td.ENTROPY_NN
    for pair, expected in SANTALUCIA_NN_DS.items():
        assert abs(nn[pair] - expected) < 1e-6, (
            f"NN dS[{pair}] drifted from SantaLucia 1998 Table 1: "
            f"got {nn[pair]}, expected {expected}"
        )


# ----------------------------------------------------------------------
# Owczarzy 2004 salt correction coefficient
# ----------------------------------------------------------------------

def test_owczarzy_entropy_salt_coefficient():
    """Owczarzy 2004 coefficient 0.368 must match Biochemistry 43:3537."""
    import inspect
    from neoswga.core import thermodynamics as td
    source = inspect.getsource(td)
    assert "0.368" in source, (
        "Owczarzy 2004 entropy salt coefficient 0.368 should appear in "
        "thermodynamics.py"
    )


# ----------------------------------------------------------------------
# Literature-backed additive Tm coefficients at 37 C
# ----------------------------------------------------------------------

ADDITIVE_COEFFS_37C = {
    "dmso_coef": (-0.55, "Chester & Marshak 1993 / Varadaraj 1994"),
    "formamide_coef": (-0.65, "Blake & Delcourt 1996"),
    "trehalose_coef": (-3.0, "Spiess 2004 midpoint"),
    "ethanol_coef": (-0.4, "Cheng 1994"),
    "betaine_uniform_coef": (-1.2, "Rees 1993 / Henke 1997 avg"),
    "urea_coef": (-2.5, "Lesnick & Bhalla 1995 short oligos"),
    "tmac_uniform_coef": (-0.5, "Melchior 1973 uniform component"),
}


@pytest.mark.parametrize("key,expected_value", [
    (k, v[0]) for k, v in ADDITIVE_COEFFS_37C.items()
])
def test_additive_coefficient_at_37c(key, expected_value):
    """Each additive's reference Tm coefficient must match its citation."""
    from neoswga.core.mechanistic_params import MECHANISTIC_MODEL_PARAMS
    tm_params = MECHANISTIC_MODEL_PARAMS["tm"]
    actual = tm_params[key]
    assert abs(actual - expected_value) < 1e-6, (
        f"{key} = {actual}, expected {expected_value}. "
        f"Update both the constant and docs/SCIENCE_CITATIONS.md if this "
        f"change is intentional."
    )


# ----------------------------------------------------------------------
# Polymerase characteristics — primary literature
# ----------------------------------------------------------------------

def test_phi29_processivity_matches_blanco_1989():
    """phi29 processivity 70,000 bp per Blanco 1989 JBC 264:8935."""
    from neoswga.core.reaction_conditions import POLYMERASE_CHARACTERISTICS
    assert POLYMERASE_CHARACTERISTICS["phi29"]["processivity"] == 70_000


def test_equiphi29_processivity_vendor_data():
    from neoswga.core.reaction_conditions import POLYMERASE_CHARACTERISTICS
    assert POLYMERASE_CHARACTERISTICS["equiphi29"]["processivity"] == 80_000


def test_bst_processivity_matches_notomi_2000():
    from neoswga.core.reaction_conditions import POLYMERASE_CHARACTERISTICS
    assert POLYMERASE_CHARACTERISTICS["bst"]["processivity"] == 2_000


def test_klenow_processivity_matches_bambara_1978():
    from neoswga.core.reaction_conditions import POLYMERASE_CHARACTERISTICS
    assert POLYMERASE_CHARACTERISTICS["klenow"]["processivity"] == 10_000


@pytest.mark.parametrize("polymerase,expected_temp", [
    ("phi29", 30.0),
    ("equiphi29", 42.0),
    ("bst", 63.0),
    ("klenow", 37.0),
])
def test_polymerase_optimal_temperatures(polymerase, expected_temp):
    from neoswga.core.reaction_conditions import POLYMERASE_CHARACTERISTICS
    assert POLYMERASE_CHARACTERISTICS[polymerase]["optimal_temp"] == expected_temp


# ----------------------------------------------------------------------
# Empirical tuning parameters — bounded sanity checks
# ----------------------------------------------------------------------

def test_phi29_dmso_threshold_in_sensible_range():
    """phi29 DMSO threshold is empirical; assert it stays within the
    practitioner-accepted 4-6% window."""
    from neoswga.core.mechanistic_params import MECHANISTIC_MODEL_PARAMS
    threshold = MECHANISTIC_MODEL_PARAMS["enzyme"]["phi29"]["dmso_threshold"]
    assert 4.0 <= threshold <= 6.0


def test_betaine_peak_in_operating_range():
    """Betaine peak enhancement is empirical; assert 0.5-2 M operating range."""
    from neoswga.core.mechanistic_params import MECHANISTIC_MODEL_PARAMS
    peak = MECHANISTIC_MODEL_PARAMS["enzyme"]["betaine_peak"]
    assert 0.5 <= peak <= 2.0


def test_mg_optimal_consistent_with_rahman_2014():
    """Mg2+ optimum should be in the 2-3 mM range for balanced GC targets,
    per Rahman et al. (2014) PLoS One 9:e112515."""
    from neoswga.core.mechanistic_params import MECHANISTIC_MODEL_PARAMS
    mg_opt = MECHANISTIC_MODEL_PARAMS["enzyme"]["mg_optimal"]
    assert 2.0 <= mg_opt <= 3.0


# ----------------------------------------------------------------------
# Dose-response sanity: applied additive shifts Tm by the expected amount
# ----------------------------------------------------------------------

def test_dmso_5pct_shifts_tm_by_roughly_minus_2_75():
    """5% DMSO on a 50% GC primer should lower Tm by about -2.75 C
    (5 * 0.55). Tolerance widened to account for Arrhenius correction
    at non-37 C temperatures."""
    from neoswga.core.reaction_conditions import ReactionConditions

    cond = ReactionConditions(temp=37.0, polymerase="phi29", dmso_percent=5.0)
    shift = cond.calculate_tm_correction(gc_content=0.5, primer_length=12)
    # Expected ~-2.75 C at 37 C reference.
    assert -3.2 <= shift <= -2.3, f"5% DMSO Tm shift out of range: got {shift}"


def test_betaine_1m_shifts_tm_by_about_minus_1_2():
    """1 M betaine at 50 % GC should lower Tm by ~-1.2 C."""
    from neoswga.core.reaction_conditions import ReactionConditions
    cond = ReactionConditions(temp=37.0, polymerase="phi29", betaine_m=1.0)
    shift = cond.calculate_tm_correction(gc_content=0.5, primer_length=12)
    assert -1.7 <= shift <= -0.7, f"1 M betaine Tm shift out of range: got {shift}"


def test_formamide_10pct_shifts_tm_by_about_minus_6_5():
    """10 % formamide should lower Tm by ~-6.5 C (10 * 0.65)."""
    from neoswga.core.reaction_conditions import ReactionConditions
    cond = ReactionConditions(temp=37.0, polymerase="phi29", formamide_percent=10.0)
    shift = cond.calculate_tm_correction(gc_content=0.5, primer_length=12)
    assert -7.5 <= shift <= -5.5, f"10 % formamide Tm shift out of range: got {shift}"
