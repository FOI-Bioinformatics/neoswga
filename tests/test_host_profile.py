"""Background specificity without a background genome.

Metagenomic SWGA has no host to design against -- the background varies between
samples or is unknown. Host-free mode drops the background, so every design
reports `MAX_SELECTIVITY` and coverage becomes the only rankable quantity.
That is not host-free design; it is design with specificity switched off.

A compositional profile restores it. Expected k-mer counts follow from a Markov
transition table (16 numbers at order 1, 64 at order 2) rather than from a
genome. Predicting the real occupancy-weighted chr21 background load of 400
random 12-mers:

    profile fitted to        order 0   order 1   order 2   order 3
    chr21 (the actual host)    0.491     0.881     0.888     0.891
    Prevotella (foreign)       0.491     0.701     0.650     0.605

Two facts drive the design. A **wrong-organism** profile still reaches rho 0.70,
so generic and metagenomic designs are possible without knowing the host. And
higher order helps only when the profile is the real host -- on a foreign one it
makes things worse, because the extra structure it captures belongs to the wrong
genome. `recommended_order` encodes that rather than fixing one value.

Recorded so it is not rebuilt: fitting the null to the *target* and treating
over-representation as specificity does not work (rho 0.03-0.26 against real
host load, unhelpful sign). Over-representation inside the target measures
target repeats, not host absence.
"""

import math

import pytest

from neoswga.core.host_profile import (
    ORDER_FOR_FOREIGN_HOST,
    ORDER_FOR_OWN_HOST,
    HostProfile,
    aggregate_loads,
    expected_site_load,
    fit_profile,
    recommended_order,
)

UNIFORM = "ACGT" * 2000


# ----------------------------------------------------------------------
# Fitting
# ----------------------------------------------------------------------


def test_a_uniform_sequence_gives_uniform_base_frequencies():
    profile = fit_profile([UNIFORM], order=0)

    for base in "ACGT":
        assert profile.base_freq[base] == pytest.approx(0.25)


def test_it_learns_composition_bias():
    """The order-0 signal: a GC-poor background makes GC-rich primers rare."""
    profile = fit_profile(["AT" * 5000], order=0)

    assert profile.base_freq["A"] == pytest.approx(0.5)
    assert profile.base_freq["G"] == 0.0


def test_it_learns_dinucleotide_bias():
    """The order-1 signal, and the reason a profile beats GC content: vertebrate
    genomes are ~5x CpG-depleted, so a primer containing CG is far rarer than
    base composition alone predicts."""
    # A sequence where C is never followed by G.
    profile = fit_profile(["CACACACATATA" * 500], order=1)

    assert profile.transitions["C"]["G"] == 0.0
    assert profile.transitions["C"]["A"] > 0.9


def test_non_acgt_breaks_the_context_rather_than_being_counted():
    """A transition observed across an N-run is not a transition. Counting it
    would smear real structure toward uniform, which is exactly the signal the
    profile exists to capture."""
    profile = fit_profile(["AAAA" + "N" * 50 + "GGGG"], order=1)

    assert profile.transitions.get("A", {}).get("G", 0.0) == 0.0
    assert "N" not in profile.base_freq


def test_it_streams_across_chunk_boundaries():
    """A multi-GB genome is read in chunks; a context spanning a boundary must
    still be counted or the fit silently loses one transition per line of the
    FASTA."""
    whole = fit_profile(["ACGTACGTACGT" * 100], order=1)
    chunked = fit_profile(["ACGTACGT" * 100, "ACGT" * 100], order=1)

    assert chunked.transitions["A"]["C"] == pytest.approx(whole.transitions["A"]["C"], abs=0.05)


def test_an_empty_sequence_is_an_error_not_a_uniform_profile():
    """Returning a uniform profile would silently produce a background model
    fitted to nothing, indistinguishable from one fitted to a real genome."""
    with pytest.raises(ValueError, match="No ACGT"):
        fit_profile(["NNNN"], order=1)


def test_a_negative_order_is_rejected():
    with pytest.raises(ValueError, match="order"):
        fit_profile([UNIFORM], order=-1)


# ----------------------------------------------------------------------
# Expected counts
# ----------------------------------------------------------------------


def test_expected_count_matches_the_uniform_analytic_value():
    """Under a uniform order-0 model a 12-mer has probability 4^-12 per strand,
    and both strands are counted because jellyfish -C stores canonical k-mers."""
    profile = fit_profile([UNIFORM], order=0)
    length = 1_000_000

    expected = profile.expected_count("ACGTACGTACGT", length)

    assert expected == pytest.approx(2 * (length - 11) * 4.0**-12, rel=1e-6)


def test_a_palindrome_is_not_counted_twice():
    """A k-mer equal to its own reverse complement is one canonical entry, not
    two. Doubling the forward probability would overstate it."""
    profile = fit_profile([UNIFORM], order=0)
    length = 1_000_000

    palindrome = profile.expected_count("ACGCGCGCGCGT", length)
    non_palindrome = profile.expected_count("AAAAAAAAAAAA", length)

    assert palindrome == pytest.approx((length - 11) * 2 * 4.0**-12, rel=1e-6)
    assert non_palindrome == pytest.approx(2 * (length - 11) * 4.0**-12, rel=1e-6)


def test_rare_composition_gives_a_lower_expected_count():
    """The whole point: a primer made of bases the background lacks should be
    predicted rare there."""
    at_rich = fit_profile(["AT" * 20000 + "GC" * 500], order=0)
    length = 1_000_000

    assert at_rich.expected_count("GCGCGCGCGCGC", length) < at_rich.expected_count(
        "ATATATATATAT", length
    )


def test_probabilities_are_computed_in_log_space():
    """A 12-mer's probability is ~1e-8 and a set's product underflows float64.
    Checked directly because an underflow would silently make every long primer
    look infinitely specific."""
    profile = fit_profile([UNIFORM], order=1)

    assert profile.log_probability("A" * 40) < -50
    assert math.isfinite(profile.log_probability("A" * 40))


def test_a_zero_length_genome_expects_nothing():
    profile = fit_profile([UNIFORM], order=0)

    assert profile.expected_count("ACGTACGTACGT", 0) == 0.0


# ----------------------------------------------------------------------
# Model order follows provenance
# ----------------------------------------------------------------------


def test_order_is_lower_for_a_foreign_profile():
    """Measured, not assumed. Predicting real chr21 load from a Prevotella
    profile: order 1 gives rho 0.701, order 2 gives 0.650, order 3 gives 0.605.
    On the actual host the ordering reverses (0.881 / 0.888 / 0.891). What the
    extra order captures is the profile organism's own structure."""
    assert recommended_order(is_actual_host=False) == ORDER_FOR_FOREIGN_HOST
    assert recommended_order(is_actual_host=True) == ORDER_FOR_OWN_HOST
    assert ORDER_FOR_FOREIGN_HOST < ORDER_FOR_OWN_HOST


# ----------------------------------------------------------------------
# The modelled load mirrors the real one
# ----------------------------------------------------------------------


def test_modelled_load_has_the_same_shape_as_the_real_one():
    """`expected_site_load` mirrors `occupancy.weighted_site_load` so a modelled
    and a real background are one quantity measured two ways, and comparable."""
    from neoswga.core.reaction_conditions import ReactionConditions

    profile = fit_profile([UNIFORM], order=0)
    conditions = ReactionConditions(temp=42.0, polymerase="equiphi29")

    load = expected_site_load(["ACGTACGTACGT"], profile, 1_000_000, conditions, 1)

    assert load > 0
    assert math.isfinite(load)


def test_more_mismatch_classes_mean_more_load():
    """Near-matches add background, never remove it."""
    from neoswga.core.reaction_conditions import ReactionConditions

    profile = fit_profile([UNIFORM], order=0)
    conditions = ReactionConditions(temp=42.0, polymerase="equiphi29")

    exact = expected_site_load(["ACGTACGTACGT"], profile, 1_000_000, conditions, 0)
    with_mm = expected_site_load(["ACGTACGTACGT"], profile, 1_000_000, conditions, 1)

    assert with_mm > exact


def test_mismatched_variants_get_their_own_expected_count():
    """A variant has a different composition and so a different abundance.
    Reusing the perfect match's expectation for all 3k neighbours would ignore
    that, which matters most in a compositionally skewed background -- exactly
    where the model is meant to help."""
    from neoswga.core.reaction_conditions import ReactionConditions

    at_rich = fit_profile(["AT" * 20000 + "GC" * 500], order=0)
    conditions = ReactionConditions(temp=42.0, polymerase="equiphi29")

    # Its 1-mismatch neighbours differ in how many G/C they carry, so a
    # composition-blind implementation would return the same load for both.
    gc_primer = expected_site_load(["GCGCGCGCGCGC"], at_rich, 1_000_000, conditions, 1)
    at_primer = expected_site_load(["ATATATATATAT"], at_rich, 1_000_000, conditions, 1)

    assert at_primer > gc_primer * 10


# ----------------------------------------------------------------------
# Panels: what "the background" means when the host is unknown
# ----------------------------------------------------------------------


def test_worst_case_is_the_default_for_an_unknown_host():
    """A set clean against one host and filthy against another is not generic.
    Averaging hides that; the maximum is what a design must survive."""
    assert aggregate_loads([1.0, 50.0], mode="worst-case") == 50.0


def test_sum_is_available_for_backgrounds_present_together():
    """A defined community or co-infection is a different question, and there
    every background really is present at once."""
    assert aggregate_loads([1.0, 50.0], mode="sum") == 51.0


def test_the_existing_behaviour_is_sum_so_the_choice_must_be_explicit():
    """`weighted_site_load` sums across prefixes. Applying that to a panel of
    ALTERNATIVE hosts would penalise a primer for background it will never
    meet, which is why the mode is named rather than inferred."""
    assert aggregate_loads([2.0, 3.0], mode="sum") == 5.0
    assert aggregate_loads([2.0, 3.0], mode="worst-case") == 3.0
    assert aggregate_loads([2.0, 3.0], mode="mean") == 2.5


def test_an_empty_panel_is_zero_rather_than_an_error():
    assert aggregate_loads([]) == 0.0


def test_an_unknown_mode_is_rejected():
    """Silently defaulting would pick a semantics the user did not choose."""
    with pytest.raises(ValueError, match="aggregation mode"):
        aggregate_loads([1.0], mode="average")


# ----------------------------------------------------------------------
# What the model is and is not evidence for
# ----------------------------------------------------------------------


def test_a_modelled_background_is_labelled_as_modelled():
    """A predicted background must never be read as a measured one. The mode is
    "modelled", distinct from both "occupancy" (real counts, occupancy-weighted)
    and "exact" (real counts, unweighted)."""
    import numpy as np

    h5py = pytest.importorskip("h5py")
    import os
    import tempfile

    from neoswga.core.base_optimizer import OptimizerConfig
    from neoswga.core.dominating_set_adapter import DominatingSetAdapter
    from neoswga.core.position_cache import PositionCache
    from neoswga.core.reaction_conditions import ReactionConditions
    from neoswga.core.thermodynamics import reverse_complement

    primer = "GCGCGGCCGC"
    directory = tempfile.mkdtemp()
    prefix = os.path.join(directory, "fg")
    with h5py.File(f"{prefix}_10mer_positions.h5", "w") as handle:
        for key in sorted({primer, reverse_complement(primer)}):
            handle.create_dataset(key, data=np.array([100, 5000], dtype=np.int32))
    with open(f"{prefix}_10mer_all.txt", "w") as handle:
        handle.write(f"{min(primer, reverse_complement(primer))} 2\n")

    cache = PositionCache([prefix], [primer])
    conditions = ReactionConditions(temp=42.0, polymerase="equiphi29")

    def build(**kw):
        return DominatingSetAdapter(
            position_cache=cache,
            fg_prefixes=[prefix],
            fg_seq_lengths=[20000],
            bg_prefixes=[],
            bg_seq_lengths=[],
            config=OptimizerConfig(target_set_size=1, extension_reach=1000),
            conditions=conditions,
            **kw,
        ).compute_metrics([primer])

    # Without a profile there is no background term at all, and selectivity is
    # the "nothing detected" sentinel -- which cannot rank two designs.
    assert build().selectivity_mode == "exact"

    profile = fit_profile([UNIFORM], order=0)
    modelled = build(background_profile=profile, background_profile_lengths=[3_000_000])
    assert modelled.selectivity_mode == "modelled"
    assert 0 < modelled.selectivity_ratio < 1e6


def test_the_optimizer_actually_receives_the_profile():
    """Every optimizer subclass forwards **kwargs to BaseOptimizer. They did
    not: each listed its own arguments and dropped the rest, which made a
    compositional background silently inert on all five methods while still
    appearing to be configured."""
    import inspect

    from neoswga.core.background_aware_optimizer import BackgroundAwareBaseOptimizer
    from neoswga.core.clique_optimizer import CliqueOptimizer
    from neoswga.core.dominating_set_adapter import DominatingSetAdapter
    from neoswga.core.hybrid_optimizer import HybridBaseOptimizer
    from neoswga.core.network_optimizer import NetworkBaseOptimizer

    for cls in (
        DominatingSetAdapter,
        CliqueOptimizer,
        NetworkBaseOptimizer,
        HybridBaseOptimizer,
        BackgroundAwareBaseOptimizer,
    ):
        source = inspect.getsource(cls.__init__)
        assert "**kwargs" in source, f"{cls.__name__} does not accept extra kwargs"
        assert "**kwargs,\n        )" in source or "**kwargs)" in source, (
            f"{cls.__name__} accepts **kwargs but does not forward them to super(); "
            f"a background profile passed to it would be silently ignored"
        )
