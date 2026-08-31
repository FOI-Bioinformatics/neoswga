"""Behavioural tests for the adaptive filtering pipeline.

`adaptive_filters` is the path that lets extreme-GC organisms work at all: a
fixed 37.5-62.5% primer GC window rejects everything for Wolbachia (35% GC) or
Burkholderia (67% GC). Its filters were changed earlier in this audit so that
strongly AT-rich targets are not excluded by a positive lower bound -- published
sets against *P. falciparum* and *Borrelia* use primers with zero GC.

The tests below check the filters do what their names claim across the GC range,
and that the widening is confined to compositionally extreme targets rather than
quietly loosening the normal case.
"""

import pytest

from neoswga.core.adaptive_filters import (
    AdaptiveGCFilter,
    GCClampFilter,
    RepeatFilter,
    ThermodynamicFilter,
)
from neoswga.core.parameter import EXTREME_AT_GENOME_GC, EXTREME_GC_GENOME_GC

# ----------------------------------------------------------------------
# AdaptiveGCFilter
# ----------------------------------------------------------------------


@pytest.mark.parametrize(
    "genome_gc,label",
    [(0.35, "Wolbachia"), (0.50, "E. coli"), (0.67, "Burkholderia")],
)
def test_window_is_centred_on_the_genome_for_normal_targets(genome_gc, label):
    """The whole point: a primer matching the target's composition must pass."""
    f = AdaptiveGCFilter(genome_gc, tolerance=0.15)
    assert f.gc_min <= genome_gc <= f.gc_max, f"{label}: window excludes its own genome GC"


@pytest.mark.parametrize("genome_gc", [0.05, 0.19, 0.28])
def test_at_rich_targets_admit_zero_gc_primers(genome_gc):
    """Published AT-rich sets use primers with no G or C at all.

    A window of genome_gc +/- tolerance cannot reach zero at any tolerance a
    normal target would want, so on extreme targets the lower bound is released.
    """
    assert genome_gc < EXTREME_AT_GENOME_GC
    f = AdaptiveGCFilter(genome_gc)
    assert f.gc_min == 0.0
    assert f.passes("AAAAAAAAAA")
    assert f.passes("ATATATATAT")


@pytest.mark.parametrize("genome_gc", [0.75, 0.85])
def test_gc_rich_targets_admit_all_gc_primers(genome_gc):
    """The mirror image, released on the other side."""
    assert genome_gc > EXTREME_GC_GENOME_GC
    f = AdaptiveGCFilter(genome_gc)
    assert f.gc_max == 1.0
    assert f.passes("GCGCGCGCGC")


def test_widening_does_not_leak_into_normal_targets():
    """The change must be confined; a 50% GC target keeps its old window."""
    f = AdaptiveGCFilter(0.50, tolerance=0.15)
    assert (f.gc_min, f.gc_max) == (0.35, 0.65)
    assert not f.passes("AAAAAAAAAA")
    assert not f.passes("GCGCGCGCGC")


def test_at_rich_target_still_rejects_gc_rich_primers():
    """Releasing the lower bound must not release the upper one as well."""
    f = AdaptiveGCFilter(0.19)
    assert not f.passes("GCGCGCGCGC"), "AT-rich target accepted a 100% GC primer"


def test_rejection_is_explained():
    """`explain_rejection` is what tells a user why their pool came back empty."""
    f = AdaptiveGCFilter(0.50)
    assert "PASSES" in f.explain_rejection("ACGTACGTAC")
    assert "low" in f.explain_rejection("AAAAAAAAAA")
    assert "high" in f.explain_rejection("GCGCGCGCGC")


# ----------------------------------------------------------------------
# GCClampFilter
# ----------------------------------------------------------------------


def test_clamp_is_waived_only_for_at_rich_targets():
    """A zero-GC primer cannot satisfy a GC clamp at any threshold.

    The clamp is PCR primer-design practice; it does not transfer to SWGA
    against an AT-rich genome, where AT-rich primers are the point.
    """
    at_rich = GCClampFilter(genome_gc=0.19)
    normal = GCClampFilter(genome_gc=0.50)

    assert at_rich.passes("AAAAAAAAAA")
    assert not normal.passes("AAAAAAAAAA")


def test_clamp_default_is_unchanged_when_target_gc_is_unknown():
    """Callers that pass no genome GC keep the previous behaviour."""
    assert GCClampFilter().min_gc == 1
    assert not GCClampFilter().passes("AAAAAAAAAA")


def test_clamp_still_rejects_an_over_strong_3_prime_end():
    """Too much GC in the last five bases promotes mispriming."""
    f = GCClampFilter(genome_gc=0.50)
    assert not f.passes("AAAAAGCGCG"), "four G/C in the last five should be rejected"
    assert f.passes("AAAAAAAAGC")


def test_short_primers_are_not_clamp_checked():
    """With fewer than five bases the rule has no meaning."""
    assert GCClampFilter().passes("ACGT")


# ----------------------------------------------------------------------
# RepeatFilter
# ----------------------------------------------------------------------


def test_repeat_filter_rejects_long_homopolymers():
    f = RepeatFilter(max_homopolymer=5)
    assert not f.passes("AAAAAAACGT")
    assert f.passes("AAAACGTACG")


def test_repeat_filter_rejects_dinucleotide_runs():
    f = RepeatFilter(max_dinuc_repeat=5)
    assert not f.passes("ATATATATATAT")


def test_repeat_filter_is_more_permissive_than_pcr_rules():
    """SWGA tolerates more repetition than PCR; a moderately repetitive primer
    that a PCR designer would reject should survive here."""
    f = RepeatFilter()
    assert f.passes("ATATACGTAC")


# ----------------------------------------------------------------------
# ThermodynamicFilter
# ----------------------------------------------------------------------


def test_tm_window_tracks_the_reaction_temperature():
    """SWGA is isothermal: the target Tm sits just above the reaction temp,
    which is why PCR's 55-65 C rules do not apply."""
    cold = ThermodynamicFilter(reaction_temp=30.0)
    hot = ThermodynamicFilter(reaction_temp=63.0)
    assert hot.tm_min > cold.tm_min
    assert cold.tm_min < 30.0 + 15


def test_explicit_tm_bounds_win():
    f = ThermodynamicFilter(reaction_temp=30.0, min_tm=10.0, max_tm=90.0)
    assert (f.tm_min, f.tm_max) == (10.0, 90.0)


def test_degenerate_sequence_does_not_raise():
    """A sequence too short to have a nearest-neighbour stack must not explode.

    It warns and returns -273.15 (0 K) rather than raising, so the filter
    rejects it as far below any Tm window. Pinned because the value is
    surprising enough that someone might "tidy" it into an exception.
    """
    f = ThermodynamicFilter(reaction_temp=30.0)
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert f._calculate_tm("") < f.tm_min


def test_filter_alias_matches_passes():
    f = ThermodynamicFilter(reaction_temp=30.0, min_tm=0.0, max_tm=100.0)
    primer = "ACGTACGTACGT"
    assert f.passes_filter(primer) == f.passes(primer)


# ----------------------------------------------------------------------
# Case handling, which the filters inherit from thermodynamics
# ----------------------------------------------------------------------


@pytest.mark.parametrize("primer", ["ACGTACGTACGT", "GAATTC", "TTGACCATGA", "GCATTACGGT"])
def test_melting_temperature_is_case_insensitive(primer):
    """Lowercase must not change the answer.

    `reverse_complement` uppercases its output, so the self-complementarity
    check compared a raw sequence against an uppercased one and returned False
    for any lowercase palindrome. Those sequences then lost the SantaLucia
    symmetry correction, shifting Tm by about 1.5 C -- silently, and precisely
    for the self-complementary primers whose thermodynamics matter most.

    Lowercase reaches the pipeline from soft-masked FASTA and hand-typed
    primer lists.
    """
    import warnings

    from neoswga.core.thermodynamics import calculate_tm_with_salt

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert calculate_tm_with_salt(primer) == pytest.approx(
            calculate_tm_with_salt(primer.lower())
        )


def test_palindrome_detection_is_case_insensitive():
    from neoswga.core.thermodynamics import is_palindrome

    assert is_palindrome("GAATTC")
    assert is_palindrome("gaattc")
    assert is_palindrome("GaAtTc")
    assert not is_palindrome("AAAAAA")
