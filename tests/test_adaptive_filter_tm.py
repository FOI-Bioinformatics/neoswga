"""Test that adaptive_filters uses the canonical thermodynamics module for Tm."""
import pytest


def test_thermodynamic_filter_tm_matches_canonical():
    """ThermodynamicFilter Tm must agree with thermodynamics.calculate_tm_with_salt."""
    from neoswga.core.adaptive_filters import ThermodynamicFilter
    from neoswga.core.thermodynamics import calculate_tm_with_salt

    tf = ThermodynamicFilter(min_tm=15.0, max_tm=60.0, na_conc=50.0)

    test_seqs = [
        'ATCGATCGATCG',
        'GCGCGCGCGC',
        'AAAATTTTCCCC',
        'ATATATATATAT',
        'GCTAGCTAGCTA',
    ]

    for seq in test_seqs:
        adaptive_tm = tf._calculate_tm(seq)
        canonical_tm = calculate_tm_with_salt(seq, na_conc=50.0)
        assert abs(adaptive_tm - canonical_tm) < 0.01, \
            f"Tm mismatch for {seq}: adaptive={adaptive_tm:.1f}, canonical={canonical_tm:.1f}"


def test_thermodynamic_filter_accepts_valid_primers():
    """Primers with Tm in range should pass the filter."""
    from neoswga.core.adaptive_filters import ThermodynamicFilter

    tf = ThermodynamicFilter(min_tm=15.0, max_tm=60.0, na_conc=50.0)
    assert tf.passes_filter('ATCGATCGATCG')


def test_thermodynamic_filter_rejects_extreme_primers():
    """Very AT-rich short primers should have low Tm and be rejected with narrow window."""
    from neoswga.core.adaptive_filters import ThermodynamicFilter

    tf = ThermodynamicFilter(min_tm=50.0, max_tm=70.0, na_conc=50.0)
    result = tf.passes_filter('AAAAAA')
    assert result is False
