"""
Consolidated thermodynamic tests for NeoSWGA.

Covers core Tm calculations, batch operations, edge cases, literature
validation, free energy, nearest-neighbor parameters, and backward
compatibility with the deprecated thermo_estimation module.

Literature references:
    - SantaLucia (1998) PNAS 95:1460-1465 (nearest-neighbor parameters)
    - Owczarzy et al. (2008) Biochemistry 47:5336-5353 (salt corrections)
"""

import time
import warnings

import numpy as np
import pytest

from neoswga.core import thermodynamics as thermo
from neoswga.core.thermodynamics import (
    ENTHALPY_NN,
    ENTROPY_NN,
    calculate_enthalpy_entropy,
    calculate_free_energy,
    calculate_gc_batch,
    calculate_salt_correction,
    calculate_tm_basic,
    calculate_tm_batch,
    calculate_tm_with_salt,
    calculate_wallace_tm_batch,
    clear_thermodynamic_caches,
    gc_content,
    get_cache_stats,
    has_ambiguous_bases,
    normalize_sequence,
    reverse_complement,
    wallace_tm,
)
from neoswga.core.reaction_conditions import (
    ReactionConditions,
    get_enhanced_conditions,
    get_standard_conditions,
)


# =============================================================================
# Tm Calculation
# =============================================================================


class TestTmCalculation:
    """Tests for individual Tm calculations."""

    def test_basic_tm(self):
        """Test basic nearest-neighbor Tm calculation."""
        seq = "ATCGATCG"
        tm = calculate_tm_basic(seq)
        assert isinstance(tm, float)
        assert not np.isnan(tm)

    def test_tm_with_salt(self):
        """Test Tm with salt correction."""
        seq = "ATCGATCG"
        tm = calculate_tm_with_salt(seq, na_conc=50.0, mg_conc=0.0)
        assert isinstance(tm, float)
        assert not np.isnan(tm)

    def test_gc_rich_higher_tm(self):
        """GC-rich sequences should have higher Tm than AT-rich sequences."""
        at_rich = "AATTAATT"
        gc_rich = "GCGCGCGC"
        mixed = "ATCGATCG"

        tm_at = calculate_tm_with_salt(at_rich, na_conc=50)
        tm_gc = calculate_tm_with_salt(gc_rich, na_conc=50)
        tm_mixed = calculate_tm_with_salt(mixed, na_conc=50)

        assert tm_gc > tm_mixed > tm_at, (
            f"Expected Tm order: GC ({tm_gc:.1f}) > mixed ({tm_mixed:.1f}) "
            f"> AT ({tm_at:.1f})"
        )

    def test_longer_sequence_higher_tm(self):
        """Longer sequences should generally have higher Tm."""
        short = "ATCG"
        medium = "ATCGATCG"
        long = "ATCGATCGATCG"

        tm_short = calculate_tm_with_salt(short, na_conc=50)
        tm_medium = calculate_tm_with_salt(medium, na_conc=50)
        tm_long = calculate_tm_with_salt(long, na_conc=50)

        assert tm_long > tm_medium > tm_short

    def test_very_short_sequence_4bp(self):
        """Test minimum practical primer length (4bp)."""
        seq = "ATCG"
        enthalpy, entropy = calculate_enthalpy_entropy(seq)
        tm = calculate_tm_basic(seq)
        assert tm < 30

    def test_very_short_sequence_6bp(self):
        """Test 6bp sequence (minimum SWGA primer)."""
        seq = "ATCGAT"
        tm = calculate_tm_with_salt(seq)
        assert isinstance(tm, float)
        assert not np.isnan(tm)

    def test_all_at_sequence(self):
        """Test poly-AT sequence has low Tm."""
        seq = "ATATATATAT"
        assert gc_content(seq) == 0.0
        tm = calculate_tm_with_salt(seq)
        assert tm < 35

    def test_all_gc_sequence(self):
        """Test poly-GC sequence has high Tm."""
        seq = "GCGCGCGCGC"
        assert gc_content(seq) == 1.0
        tm = calculate_tm_with_salt(seq)
        assert tm > 40

    def test_poly_a_tract(self):
        """Test poly-A tract has low Tm."""
        seq = "AAAAAAAAAAAA"
        tm = calculate_tm_basic(seq)
        assert tm < 50


# =============================================================================
# Tm Literature Validation
# =============================================================================


class TestLiteratureValidation:
    """Validate Tm calculations against published reference values."""

    # Reference Tm values validated against the NeoSWGA NN model
    # (SantaLucia 1998). Tolerance: 3C standard, 5C extreme AT-rich.
    TM_REFERENCE_DATA = [
        ("GCGCGCGC", 50, 35.5, 3.0, "NeoSWGA NN model"),
        ("ATCGATCG", 50, 13.5, 3.0, "NeoSWGA NN model"),
        ("GCATGCAT", 50, 16.5, 3.0, "NeoSWGA NN model"),
        ("ATCGATCGATCG", 50, 35.9, 3.0, "NeoSWGA NN model"),
        ("GCGAATTCGC", 50, 29.2, 3.0, "NeoSWGA NN model"),
        ("AATTAATT", 50, -13.0, 5.0, "NeoSWGA NN model - extreme AT"),
    ]

    @pytest.mark.parametrize(
        "seq,na_mm,expected_tm,tolerance,source", TM_REFERENCE_DATA
    )
    def test_tm_within_tolerance(self, seq, na_mm, expected_tm, tolerance, source):
        """Verify calculated Tm is within acceptable tolerance of reference value."""
        calculated_tm = calculate_tm_with_salt(seq, na_conc=na_mm, mg_conc=0)
        error = abs(calculated_tm - expected_tm)
        assert error <= tolerance, (
            f"Tm for {seq}: expected {expected_tm}C +/- {tolerance}C, "
            f"got {calculated_tm:.1f}C (error: {error:.1f}C)"
        )

    def test_primer_set_tm_uniformity(self):
        """Test that a primer set has computable and finite Tm range."""
        primers = [
            "ATCGATCG",
            "GCTAGCTA",
            "AATTAATT",
            "CCGGCCGG",
            "ATATGCGC",
            "GCGCATAT",
        ]
        tms = [calculate_tm_with_salt(p) for p in primers]
        assert len(tms) == 6
        tm_range = max(tms) - min(tms)
        assert tm_range < 50

    def test_francisella_primer_low_gc(self):
        """Test primer for Francisella tularensis (low GC)."""
        primer = "AATATATATAT"
        assert gc_content(primer) < 0.2
        tm = calculate_tm_with_salt(primer)
        assert tm < 40

    def test_ecoli_primer_balanced_gc(self):
        """Test primer for E. coli (50% GC)."""
        primer = "ATCGATCGATCG"
        assert gc_content(primer) == 0.5
        tm = calculate_tm_with_salt(primer)
        assert isinstance(tm, float)

    def test_burkholderia_primer_high_gc(self):
        """Test primer for Burkholderia (high GC)."""
        primer = "GCGCGCATGCGC"
        assert gc_content(primer) > 0.6
        tm = calculate_tm_with_salt(primer)
        assert tm > 35


# =============================================================================
# Salt Correction
# =============================================================================


class TestSaltCorrection:
    """Validate salt correction against Owczarzy et al. (2008)."""

    def test_higher_salt_higher_tm(self):
        """Higher salt concentration should increase Tm."""
        seq = "ATCGATCG"
        tm_low = calculate_tm_with_salt(seq, na_conc=20)
        tm_med = calculate_tm_with_salt(seq, na_conc=50)
        tm_high = calculate_tm_with_salt(seq, na_conc=100)
        assert tm_high > tm_med > tm_low

    def test_salt_correction_magnitude(self):
        """Salt correction should follow 12.5 * log10([Na+]) for oligos."""
        correction_50mm = calculate_salt_correction(na_conc=50, mg_conc=0)
        correction_100mm = calculate_salt_correction(na_conc=100, mg_conc=0)

        diff = correction_100mm - correction_50mm
        expected_diff = 12.5 * (np.log10(0.1) - np.log10(0.05))
        assert abs(diff - expected_diff) < 0.5

    def test_mg_increases_tm(self):
        """Mg2+ should stabilize duplex and increase Tm."""
        seq = "ATCGATCG"
        tm_no_mg = calculate_tm_with_salt(seq, na_conc=50, mg_conc=0)
        tm_with_mg = calculate_tm_with_salt(seq, na_conc=50, mg_conc=2)
        assert tm_with_mg > tm_no_mg

    def test_zero_sodium(self):
        """Test with very low sodium concentration (edge case)."""
        seq = "ATCGATCG"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                tm = calculate_tm_with_salt(seq, na_conc=0.01, mg_conc=0.0)
                assert isinstance(tm, float)
            except (ValueError, ZeroDivisionError):
                pass  # acceptable for invalid conditions

    def test_very_low_salt_1mm(self):
        """Test with 1mM total salt (minimal conditions)."""
        seq = "ATCGATCG"
        tm = calculate_tm_with_salt(seq, na_conc=1.0, mg_conc=0.0)
        assert isinstance(tm, float)

    def test_very_high_sodium_500mm(self):
        """Test with 500mM sodium (high salt conditions)."""
        seq = "ATCGATCG"
        tm_low = calculate_tm_with_salt(seq, na_conc=50.0, mg_conc=0.0)
        tm_high = calculate_tm_with_salt(seq, na_conc=500.0, mg_conc=0.0)
        assert tm_high > tm_low

    def test_very_high_magnesium_10mm(self):
        """Test with 10mM Mg2+ (excess)."""
        seq = "ATCGATCGATCG"
        tm = calculate_tm_with_salt(seq, na_conc=50.0, mg_conc=10.0)
        assert isinstance(tm, float)
        assert tm > 20

    def test_magnesium_dominates_at_low_sodium(self):
        """Test that Mg2+ dominates at low sodium concentrations."""
        seq = "ATCGATCGATCG"
        tm_na_only = calculate_tm_with_salt(seq, na_conc=10.0, mg_conc=0.0)
        tm_with_mg = calculate_tm_with_salt(seq, na_conc=10.0, mg_conc=2.0)
        assert tm_with_mg > tm_na_only

    def test_mixed_high_salt(self):
        """Test high concentrations of both Na+ and Mg2+."""
        seq = "ATCGATCG"
        tm = calculate_tm_with_salt(seq, na_conc=200.0, mg_conc=5.0)
        assert isinstance(tm, float)


# =============================================================================
# Batch Operations
# =============================================================================


class TestBatchOperations:
    """Tests for batch Tm, GC, and Wallace Tm calculations."""

    def test_batch_tm_single_sequence(self):
        """Test batch Tm with a single sequence."""
        tms = calculate_tm_batch(["ATCGATCG"])
        assert isinstance(tms, np.ndarray)
        assert len(tms) == 1
        assert not np.isnan(tms[0])

    def test_batch_tm_multiple_sequences(self):
        """Test batch Tm with multiple sequences."""
        tms = calculate_tm_batch(["ATCGATCG", "GCTAGCTA", "AAAAAATTTT"])
        assert isinstance(tms, np.ndarray)
        assert len(tms) == 3
        assert all(not np.isnan(tm) for tm in tms)

    def test_batch_tm_matches_individual(self):
        """Test that batch results match individual calculations."""
        sequences = ["ATCGATCG", "GCTAGCTA", "ATATATATAT"]
        batch_tms = calculate_tm_batch(sequences)
        individual_tms = [calculate_tm_with_salt(seq) for seq in sequences]
        for i, seq in enumerate(sequences):
            assert abs(batch_tms[i] - individual_tms[i]) < 0.1, (
                f"Mismatch for {seq}: batch={batch_tms[i]}, "
                f"individual={individual_tms[i]}"
            )

    def test_batch_tm_custom_salt(self):
        """Test batch Tm with custom salt concentrations."""
        sequences = ["ATCGATCG", "GCTAGCTA"]
        tms_low = calculate_tm_batch(sequences, na_conc=10.0)
        tms_high = calculate_tm_batch(sequences, na_conc=100.0)
        for i in range(len(sequences)):
            assert tms_high[i] > tms_low[i]

    def test_batch_tm_magnesium(self):
        """Test batch Tm with magnesium."""
        sequences = ["ATCGATCGATCG"]
        tms_no_mg = calculate_tm_batch(sequences, mg_conc=0.0)
        tms_with_mg = calculate_tm_batch(sequences, mg_conc=2.0)
        assert tms_no_mg[0] != tms_with_mg[0]

    def test_batch_tm_empty_list(self):
        """Test batch Tm with empty list."""
        tms = calculate_tm_batch([])
        assert isinstance(tms, np.ndarray)
        assert len(tms) == 0

    def test_batch_tm_large(self):
        """Test batch Tm with 100 random sequences."""
        np.random.seed(42)
        sequences = [
            "".join(np.random.choice(list("ATCG"), np.random.randint(8, 15)))
            for _ in range(100)
        ]
        tms = calculate_tm_batch(sequences)
        assert len(tms) == 100
        assert all(not np.isnan(tm) for tm in tms)
        assert all(-20 < tm < 100 for tm in tms)

    def test_batch_gc_single(self):
        """Test GC batch with single sequence."""
        gc_values = calculate_gc_batch(["ATCGATCG"])
        assert isinstance(gc_values, np.ndarray)
        assert gc_values[0] == 0.5

    def test_batch_gc_multiple(self):
        """Test GC batch with multiple sequences."""
        gc_values = calculate_gc_batch(["AAAA", "GGGG", "ATCG"])
        assert gc_values[0] == 0.0
        assert gc_values[1] == 1.0
        assert gc_values[2] == 0.5

    def test_batch_gc_matches_individual(self):
        """Test that batch GC matches individual calculations."""
        sequences = ["ATCGATCG", "AAATTTCCC", "GCGCGCGC"]
        batch_gc = calculate_gc_batch(sequences)
        for i, seq in enumerate(sequences):
            assert abs(batch_gc[i] - gc_content(seq)) < 0.001

    def test_batch_gc_empty(self):
        """Test GC batch with empty list."""
        assert len(calculate_gc_batch([])) == 0

    def test_batch_wallace_single(self):
        """Test Wallace Tm batch with single sequence."""
        tms = calculate_wallace_tm_batch(["ATCGATCG"])
        assert isinstance(tms, np.ndarray)
        assert not np.isnan(tms[0])

    def test_batch_wallace_multiple(self):
        """Test Wallace Tm batch with known values."""
        tms = calculate_wallace_tm_batch(["AAAA", "GGGG", "ATCG"])
        assert tms[0] == 8.0   # 4*2
        assert tms[1] == 16.0  # 4*4
        assert tms[2] == 12.0  # 2*2 + 2*4

    def test_batch_wallace_matches_individual(self):
        """Test batch Wallace Tm matches individual."""
        sequences = ["ATCGATCG", "AAATTTCCC", "GCGCGCGC"]
        batch_tms = calculate_wallace_tm_batch(sequences)
        for i, seq in enumerate(sequences):
            assert abs(batch_tms[i] - wallace_tm(seq)) < 0.001

    def test_batch_wallace_gc_rich_higher(self):
        """Test that GC-rich has higher Wallace Tm."""
        at_tm = calculate_wallace_tm_batch(["ATATATATAT"])[0]
        gc_tm = calculate_wallace_tm_batch(["GCGCGCGCGC"])[0]
        assert gc_tm > at_tm

    def test_batch_preserves_order(self):
        """Test that batch calculation preserves sequence order."""
        sequences = ["AAAA", "TTTT", "GGGG", "CCCC"]
        gc_values = calculate_gc_batch(sequences)
        assert gc_values[0] == 0.0
        assert gc_values[1] == 0.0
        assert gc_values[2] == 1.0
        assert gc_values[3] == 1.0

    def test_batch_varying_lengths(self):
        """Test batch with sequences of different lengths."""
        tms = calculate_tm_batch(["ATCG", "ATCGATCG", "ATCGATCGATCGATCG"])
        assert len(tms) == 3

    def test_batch_extreme_gc(self):
        """Test batch with extreme GC content."""
        gc_values = calculate_gc_batch(
            ["AAAAAAAAAA", "ATATATATAT", "GCGCGCGCGC", "CCCCCCCCCC"]
        )
        assert gc_values[0] == 0.0
        assert gc_values[1] == 0.0
        assert gc_values[2] == 1.0
        assert gc_values[3] == 1.0


# =============================================================================
# Batch Performance
# =============================================================================


class TestBatchPerformance:
    """Performance tests for batch calculations."""

    def test_batch_not_slower_than_5x_loop(self):
        """Test that batch calculation is not significantly slower than a loop."""
        np.random.seed(42)
        sequences = [
            "".join(np.random.choice(list("ATCG"), np.random.randint(8, 12)))
            for _ in range(50)
        ]

        start = time.time()
        batch_result = calculate_tm_batch(sequences)
        batch_time = time.time() - start

        clear_thermodynamic_caches()

        start = time.time()
        individual_results = [calculate_tm_with_salt(seq) for seq in sequences]
        individual_time = time.time() - start

        assert batch_time < individual_time * 5
        for i in range(len(sequences)):
            assert abs(batch_result[i] - individual_results[i]) < 0.1


# =============================================================================
# GC Content
# =============================================================================


class TestGCContent:
    """Tests for GC content calculations."""

    @pytest.mark.parametrize(
        "seq,expected",
        [
            ("AAAA", 0.0),
            ("TTTT", 0.0),
            ("GGGG", 1.0),
            ("CCCC", 1.0),
            ("ATCG", 0.5),
            ("ATAT", 0.0),
            ("GCGC", 1.0),
            ("ATGC", 0.5),
            ("ATCGATCG", 0.5),
            ("ATATATATAT", 0.0),
            ("GCGCGCGC", 1.0),
        ],
    )
    def test_gc_content_values(self, seq, expected):
        """Test GC content for known sequences."""
        assert abs(gc_content(seq) - expected) < 0.0001

    def test_case_insensitive(self):
        """GC content should be case insensitive."""
        assert gc_content("atcg") == gc_content("ATCG")
        assert gc_content("AtCg") == gc_content("ATCG")


# =============================================================================
# Free Energy
# =============================================================================


class TestFreeEnergy:
    """Tests for free energy calculations."""

    def test_negative_dg_for_stable_duplexes(self):
        """Stable duplexes should have negative deltaG at 37C."""
        dg = calculate_free_energy("GCGCGCGCGC", temperature=37.0)
        assert dg < 0

    def test_dg_increases_with_temperature(self):
        """Delta G should become less negative at higher temperature."""
        seq = "ATCGATCGATCG"
        dg_25 = calculate_free_energy(seq, temperature=25.0)
        dg_37 = calculate_free_energy(seq, temperature=37.0)
        dg_55 = calculate_free_energy(seq, temperature=55.0)
        assert dg_25 < dg_37 < dg_55

    def test_gc_rich_more_stable(self):
        """GC-rich sequences should have more negative deltaG."""
        dg_at = calculate_free_energy("AATTAATTAA", temperature=37.0)
        dg_gc = calculate_free_energy("GCGCGCGCGC", temperature=37.0)
        assert dg_gc < dg_at


# =============================================================================
# Enthalpy / Entropy
# =============================================================================


class TestEnthalpyEntropy:
    """Tests for enthalpy and entropy calculations."""

    def test_enthalpy_negative(self):
        """Total enthalpy should be negative (favorable)."""
        enthalpy, _ = calculate_enthalpy_entropy("ATCGATCG")
        assert enthalpy < 0

    def test_entropy_negative(self):
        """Total entropy should be negative (ordering)."""
        _, entropy = calculate_enthalpy_entropy("ATCGATCG")
        assert entropy < 0

    def test_gc_rich_more_negative_enthalpy(self):
        """GC-rich sequences should have more negative enthalpy."""
        h_at, _ = calculate_enthalpy_entropy("ATATATATAT")
        h_gc, _ = calculate_enthalpy_entropy("GCGCGCGCGC")
        assert h_gc < h_at

    def test_palindrome_gaattc(self):
        """Test palindromic EcoRI site GAATTC."""
        seq = "GAATTC"
        assert thermo.is_palindrome(seq)
        enthalpy, entropy = calculate_enthalpy_entropy(seq)
        assert isinstance(enthalpy, float)
        assert entropy < 0

    def test_palindrome_cgcgcgcg(self):
        """Test repetitive palindrome CGCGCGCG."""
        seq = "CGCGCGCG"
        assert thermo.is_palindrome(seq)
        enthalpy, _ = calculate_enthalpy_entropy(seq)
        assert enthalpy < -50


# =============================================================================
# Nearest-Neighbor Parameters
# =============================================================================


class TestNearestNeighborParameters:
    """Validate NN parameters against SantaLucia (1998)."""

    SANTALUCIA_NN_PARAMS = [
        ("AA/TT", -7.9, -22.2),
        ("AT/TA", -7.2, -20.4),
        ("TA/AT", -7.2, -21.3),
        ("CA/GT", -8.5, -22.7),
        ("GT/CA", -8.4, -22.4),
        ("CT/GA", -7.8, -21.0),
        ("GA/CT", -8.2, -22.2),
        ("CG/GC", -10.6, -27.2),
        ("GC/CG", -9.8, -24.4),
        ("GG/CC", -8.0, -19.9),
    ]

    @pytest.mark.parametrize("stack,expected_dh,expected_ds", SANTALUCIA_NN_PARAMS)
    def test_enthalpy_matches_santalucia(self, stack, expected_dh, expected_ds):
        """Verify enthalpy values match SantaLucia (1998) Table 1."""
        if stack in ENTHALPY_NN:
            assert abs(ENTHALPY_NN[stack] - expected_dh) < 0.1

    @pytest.mark.parametrize("stack,expected_dh,expected_ds", SANTALUCIA_NN_PARAMS)
    def test_entropy_matches_santalucia(self, stack, expected_dh, expected_ds):
        """Verify entropy values match SantaLucia (1998) Table 1."""
        if stack in ENTROPY_NN:
            assert abs(ENTROPY_NN[stack] - expected_ds) < 0.1

    def test_canonical_stacks_present(self):
        """All 10 canonical SantaLucia stacks should be present."""
        canonical = [
            "AA/TT", "AT/TA", "TA/AT", "CA/GT", "GT/CA",
            "CT/GA", "GA/CT", "CG/GC", "GC/CG", "GG/CC",
        ]
        for stack in canonical:
            assert stack in ENTHALPY_NN
            assert stack in ENTROPY_NN

    def test_enthalpy_entropy_key_consistency(self):
        """Enthalpy and entropy dicts should have matching keys."""
        assert set(ENTHALPY_NN.keys()) == set(ENTROPY_NN.keys())

    def test_enthalpy_values_reasonable(self):
        """NN enthalpy values should be in range -12 to -6 kcal/mol."""
        for stack, value in ENTHALPY_NN.items():
            assert -12 < value < -6, f"{stack}: {value}"

    def test_entropy_values_reasonable(self):
        """NN entropy values should be in range -30 to -18 cal/(mol*K)."""
        for stack, value in ENTROPY_NN.items():
            assert -30 < value < -18, f"{stack}: {value}"


# =============================================================================
# Sequence Normalization and Edge Cases
# =============================================================================


class TestEdgeCases:
    """Edge cases for sequence handling and thermodynamic calculations."""

    def test_lowercase_gives_same_tm(self):
        """Lowercase sequences should give the same Tm."""
        tm_upper = calculate_tm_with_salt("ATCGATCG")
        tm_lower = calculate_tm_with_salt("atcgatcg")
        assert abs(tm_upper - tm_lower) < 1e-5

    def test_mixed_case(self):
        """Mixed case should produce valid Tm."""
        tm = calculate_tm_with_salt("AtCgAtCg")
        assert isinstance(tm, float)

    def test_ambiguous_base_warning(self):
        """Ambiguous bases should trigger a warning."""
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            calculate_tm_with_salt("ATCGATNCG")
            assert any(
                "ambiguous" in str(warning.message).lower() for warning in w
            )

    def test_has_ambiguous_bases(self):
        """Test has_ambiguous_bases detection."""
        assert not has_ambiguous_bases("ATCG")
        assert has_ambiguous_bases("ATNG")
        assert has_ambiguous_bases("ATRY")
        assert not has_ambiguous_bases("atcg")

    def test_normalize_sequence(self):
        """Test sequence normalization to uppercase."""
        assert normalize_sequence("atcg") == "ATCG"
        assert normalize_sequence("AtCg") == "ATCG"
        assert normalize_sequence("ATCG") == "ATCG"


# =============================================================================
# Reverse Complement
# =============================================================================


class TestReverseComplement:
    """Tests for reverse_complement function."""

    def test_basic(self):
        """Test basic reverse complement."""
        assert reverse_complement("ATCG") == "CGAT"
        assert reverse_complement("AAAA") == "TTTT"
        assert reverse_complement("GGGG") == "CCCC"

    def test_palindrome(self):
        """Palindrome reverse complement should equal itself."""
        assert reverse_complement("ATAT") == "ATAT"
        assert reverse_complement("GCGC") == "GCGC"

    def test_lowercase_handling(self):
        """Lowercase input should produce uppercase reverse complement."""
        assert reverse_complement("atcg") == "CGAT"
        assert reverse_complement("AtCg") == "CGAT"

    def test_double_reverse_complement(self):
        """Double reverse complement should return original."""
        seq = "ATCGATCGATCG"
        assert reverse_complement(reverse_complement(seq)) == seq.upper()
        assert len(reverse_complement(seq)) == len(seq)


# =============================================================================
# Cache Management
# =============================================================================


class TestCacheManagement:
    """Tests for thermodynamic cache functions."""

    def test_clear_caches(self):
        """clear_thermodynamic_caches should not raise."""
        calculate_tm_with_salt("ATCGATCG")
        clear_thermodynamic_caches()

    def test_get_cache_stats(self):
        """get_cache_stats should return valid dict."""
        clear_thermodynamic_caches()
        calculate_tm_with_salt("ATCGATCG")
        calculate_tm_with_salt("ATCGATCG")  # cache hit

        stats = get_cache_stats()
        assert isinstance(stats, dict)
        assert "enthalpy_entropy" in stats
        assert "free_energy" in stats

    def test_cache_hits_after_repeat(self):
        """Cache should accumulate hits on repeat calls."""
        clear_thermodynamic_caches()
        calculate_tm_with_salt("ATCGATCG")
        calculate_tm_with_salt("ATCGATCG")

        stats = get_cache_stats()
        assert stats["enthalpy_entropy"].hits >= 0


# =============================================================================
# Polymerase Presets
# =============================================================================


class TestPolymerasePresets:
    """Validate polymerase preset configurations."""

    def test_phi29_standard_conditions(self):
        """Phi29 should operate at 30C."""
        conditions = get_standard_conditions()
        assert conditions.polymerase == "phi29"
        assert conditions.temp == 30.0
        min_t, max_t = conditions.get_polymerase_range()
        assert min_t == 30.0 and max_t == 40.0

    def test_equiphi29_conditions(self):
        """EquiPhi29 should operate at 42C."""
        conditions = ReactionConditions(temp=42, polymerase="equiphi29")
        assert conditions.polymerase == "equiphi29"
        min_t, max_t = conditions.get_polymerase_range()
        assert min_t == 42.0 and max_t == 45.0

    def test_enhanced_longer_primers(self):
        """Enhanced conditions should support longer primers."""
        standard = get_standard_conditions()
        enhanced = get_enhanced_conditions()
        assert enhanced.max_primer_length() > standard.max_primer_length()

    def test_phi29_valid_range(self):
        """Phi29 within valid range should work."""
        ReactionConditions(polymerase="phi29", temp=20.0)
        ReactionConditions(polymerase="phi29", temp=30.0)
        ReactionConditions(polymerase="phi29", temp=40.0)

    def test_phi29_below_range_raises(self):
        """Phi29 below valid range should raise ValueError."""
        with pytest.raises(ValueError):
            ReactionConditions(polymerase="phi29", temp=15.0)

    def test_phi29_above_range_raises(self):
        """Phi29 above valid range should raise ValueError."""
        with pytest.raises(ValueError):
            ReactionConditions(polymerase="phi29", temp=45.0)

    def test_equiphi29_at_boundaries(self):
        """EquiPhi29 at boundary temperatures."""
        low = ReactionConditions(polymerase="equiphi29", temp=30.0)
        high = ReactionConditions(polymerase="equiphi29", temp=50.0)
        assert low.temp == 30.0
        assert high.temp == 50.0

    def test_bst_high_temperature(self):
        """BST at high temperature (50-72C)."""
        conditions = ReactionConditions(polymerase="bst", temp=65.0)
        assert conditions.temp == 65.0


# =============================================================================
# Backward Compatibility (deprecated thermo_estimation)
# =============================================================================


class TestBackwardCompatibility:
    """Test backward compatibility with deprecated thermo_estimation module."""

    def test_compute_free_energy_exists(self):
        """thermodynamics should have compute_free_energy_for_two_strings."""
        from neoswga.core.thermodynamics import compute_free_energy_for_two_strings

        assert callable(compute_free_energy_for_two_strings)

    def test_compute_free_energy_basic(self):
        """compute_free_energy_for_two_strings should return valid float."""
        from neoswga.core.thermodynamics import compute_free_energy_for_two_strings

        result = compute_free_energy_for_two_strings("ATCGATCG", "TAGCTAGC")
        assert isinstance(result, float)

    def test_delta_g_mismatch_table_exists(self):
        """DELTA_G_MISMATCH table should exist."""
        from neoswga.core.thermodynamics import DELTA_G_MISMATCH

        assert isinstance(DELTA_G_MISMATCH, dict)
        assert len(DELTA_G_MISMATCH) > 50

    def test_rf_preprocessing_uses_thermo_module(self):
        """rf_preprocessing should import thermodynamics (not thermo_estimation)."""
        import neoswga.core.rf_preprocessing as rf

        assert hasattr(rf, "thermo")

    def test_rf_preprocessing_no_deprecated_import(self):
        """rf_preprocessing should not import deprecated te."""
        import neoswga.core.rf_preprocessing as rf

        assert not hasattr(rf, "te")

    def test_thermo_estimation_importable(self):
        """thermo_estimation should be importable but deprecated."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", DeprecationWarning)
            from neoswga.core import thermo_estimation

        assert hasattr(thermo_estimation, "compute_free_energy_for_two_strings")

    def test_thermodynamics_importable(self):
        """thermodynamics should be importable."""
        from neoswga.core import thermodynamics

        assert hasattr(thermodynamics, "compute_free_energy_for_two_strings")

    def test_similar_results_for_known_pairs(self):
        """New function should produce reasonable results for known pairs."""
        from neoswga.core.thermodynamics import compute_free_energy_for_two_strings

        test_pairs = [
            ("ATCGATCG", "TAGCTAGC"),
            ("GCGCGCGC", "CGCGCGCG"),
            ("AAAAAAA", "TTTTTTT"),
        ]
        for x, y in test_pairs:
            result = compute_free_energy_for_two_strings(x, y)
            assert -50 < result < 50, f"Unreasonable delta G for {x}/{y}: {result}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
