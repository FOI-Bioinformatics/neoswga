"""The retraining script must compute features the way serving computes them.

`scripts/retrain_rf_model.py` and `neoswga/core/rf_preprocessing.py` both fill a
feature named `melting_tm`, and they used two different functions to do it.

Serving (`rf_preprocessing.py:399`) calls `melting_temp.temp`, nearest-neighbour
thermodynamics with an Owczarzy salt correction. Training used a local
`compute_melting_temp` that applies the Wallace rule below 14 bp and a GC
regression at or above it. Measured on 400 real 12-mers from the E. coli pool:

    trained on (Wallace):  mean 37.01  sd 1.32  range 34.0-38.0
    served (NN + salt):    mean 60.79  sd 3.27  range 49.7-67.5
    Spearman(trained, served) = 0.732

The offset is not the worst of it. The training formula switches branch at
length 14, and the two branches do not meet:

    len 13    32.0 .. 44.0
    len 14    67.1 .. 87.6     <- +23.1 C of pure formula change

So no training primer of any length ever took a value between 44 and 67 C. That
band is empty by construction. Meanwhile 99% of a real 12-mer pool lands inside
it (396 of 400), because the nearest-neighbour function is continuous and puts a
12-mer at 49.7-67.5 C.

The consequence is not a shifted feature, it is an inverted one. In training,
the empty band is exactly the gap separating short primers from long ones, so a
value above it means "14 bp or longer". At serving, 99% of 12-mers claim that
value. The model reads a 12-mer as a long primer while `sequence.length`, its
second-largest feature at 25.6% of importance, says 12. Those two features carry
41.9% of total importance between them and, on a fixed-k pool, they disagree.

This file pins the two paths together so the feature cannot drift apart again.
It does NOT assert the shipped model is good; the synthetic-label question is
recorded separately in `tests/test_step3_order_is_deterministic.py`.
"""

import importlib.util
import pathlib

import pytest

_SCRIPT = pathlib.Path(__file__).resolve().parents[1] / "scripts" / "retrain_rf_model.py"


@pytest.fixture(scope="module")
def retrain():
    """Load the retraining script as a module without running it."""
    if not _SCRIPT.exists():
        pytest.skip("retrain_rf_model.py is not present")
    spec = importlib.util.spec_from_file_location("_retrain_under_test", _SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# A pool spanning the lengths the training generator actually draws from.
_SEQS = [
    "ACGTAC",
    "ACGTACGT",
    "ACGTACGTAC",
    "AAACCCGGGT",
    "ACGTACGTACGT",
    "CACCGACGACGA",
    "ACGTACGTACGTA",
    "ACGTACGTACGTAC",
    "ACGTACGTACGTACGT",
    "ACGTACGTACGTACGTAC",
]


# ---------------------------------------------------------------------------
# The property that was missing
# ---------------------------------------------------------------------------


def test_training_and_serving_agree_on_melting_tm(retrain):
    """One feature name must mean one function."""
    from neoswga.core.melting_temp import temp as serving_tm

    disagreements = []
    for seq in _SEQS:
        trained = retrain.compute_melting_temp(seq)
        served = serving_tm(seq)
        if abs(trained - served) > 0.5:
            disagreements.append(
                f"  {seq} ({len(seq)} bp): trained {trained:.1f}, served {served:.1f}"
            )

    assert not disagreements, (
        "the retraining script and the serving code fill `melting_tm` with "
        "different functions:\n" + "\n".join(disagreements)
    )


def test_the_extracted_feature_row_agrees_too(retrain):
    """Guard the whole extractor, not just the helper.

    `extract_features` is what actually builds a training row, so pinning
    `compute_melting_temp` alone would not catch a caller that recomputed it.
    """
    from neoswga.core.melting_temp import temp as serving_tm

    for seq in _SEQS:
        row = retrain.extract_features(seq)
        assert row["melting_tm"] == pytest.approx(serving_tm(seq), abs=0.5), (
            f"extract_features({seq!r}) put {row['melting_tm']:.1f} in `melting_tm` "
            f"where serving supplies {serving_tm(seq):.1f}"
        )


def test_no_discontinuity_across_the_length_14_boundary(retrain):
    """The training formula jumped 23.1 C between a 13-mer and a 14-mer.

    That jump created a band of Tm values no training primer could occupy, and
    99% of a real 12-mer pool lands in it.
    """
    thirteen = retrain.compute_melting_temp("ACGTACGTACGTA")
    fourteen = retrain.compute_melting_temp("ACGTACGTACGTAC")

    assert abs(fourteen - thirteen) < 12.0, (
        f"a 13-mer scores {thirteen:.1f} and a 14-mer {fourteen:.1f}, a "
        f"{fourteen - thirteen:+.1f} C step from adding one base; the feature is "
        "reporting which formula ran, not a melting temperature"
    )


def test_real_primers_land_where_training_data_lived(retrain):
    """The end-to-end statement: a 12-mer pool must not fall in a band the
    training distribution left empty."""
    lengths = [6, 8, 10, 12, 13, 14, 16, 18]
    trained = [retrain.compute_melting_temp(("ACGT" * 6)[:n]) for n in lengths]

    # With one continuous function, Tm rises with length. With two branches it
    # leaps. Sorting by length must not reorder the values wildly.
    assert trained == sorted(trained), (
        f"melting_tm is not monotonic in primer length across {lengths}: "
        f"{[round(t, 1) for t in trained]}"
    )


# ---------------------------------------------------------------------------
# The label, which is computed FROM the feature
# ---------------------------------------------------------------------------


def test_the_label_discriminates_on_a_realistic_pool(retrain):
    """`compute_target_score` bands `melting_tm` at 30-42 C, which was tuned to
    the Wallace scale.

    On the served scale a real 12-mer pool sits at 49.7-67.5 C, so the band
    would never fire and every primer would take the same Tm contribution. A
    label term that is constant across the pool teaches the model nothing.
    """
    seqs = [
        "ACGTACGTACGT",
        "CACCGACGACGA",
        "AAAATTTTAAAT",
        "GGGCCCGGGCCC",
        "ACGCGCGTACGT",
        "TTTATTTATTTA",
    ]
    scores = {s: retrain.compute_target_score(retrain.extract_features(s)) for s in seqs}

    assert len(set(scores.values())) > 1, (
        f"every primer in the pool got the same target score ({scores}); the "
        "label cannot rank them"
    )


def test_the_tm_term_is_not_constant_across_the_pool(retrain):
    """Isolate the Tm contribution specifically.

    The whole score can vary through GC content while the Tm term contributes a
    flat offset, which is what happened on the Wallace scale: every real primer
    fell inside 30-42 and collected the same +3.0.
    """
    base = retrain.extract_features("ACGTACGTACGT")

    contributions = set()
    for tm in (35.0, 50.0, 60.0, 70.0, 80.0):
        row = dict(base)
        row["melting_tm"] = tm
        contributions.add(retrain.compute_target_score(row))

    assert len(contributions) > 1, (
        "varying melting_tm across the served range 35-80 C did not change the "
        f"target score at all (always {contributions}); the label's Tm bands do "
        "not cover the scale the feature is computed on"
    )


# ---------------------------------------------------------------------------
# The serving side, where the tempting "fix" is the wrong one
# ---------------------------------------------------------------------------


def test_the_feature_uses_the_training_defaults_not_the_configured_chemistry():
    """`rf_preprocessing` must NOT pass a run's salt concentrations here.

    The bare `_melting_temp(primer)` call looks like an oversight: it takes
    10 mM Na and 20 mM Mg while a typical run configures 50 and 10, which puts
    the feature about 7.9 C off that reaction. Feeding the configured chemistry
    in is the obvious repair and it is the wrong one -- those defaults are the
    scale the bundled model was fitted on, so honouring the reaction here would
    move the served value off the trained distribution.

    The reaction's real Tm is computed elsewhere, and correctly:
    `filter.filter_extra` gates on `ReactionConditions.calculate_effective_tm`
    and `NetworkOptimizer._get_primer_tm` uses `calculate_tm_with_salt`.
    """
    from neoswga.core.melting_temp import temp as serving_tm
    from neoswga.core.rf_preprocessing import get_features

    primer = "CACCGACGACGA"
    # `get_features` returns [primer, target, *values in `base_features` order].
    from neoswga.core.rf_preprocessing import base_features

    row = get_features(primer)
    served_value = row[2 + base_features.index("melting_tm")]

    assert served_value == pytest.approx(serving_tm(primer), abs=1e-9), (
        "the RF feature no longer matches `melting_temp.temp` at its defaults; "
        "if this changed to honour the configured reaction, the bundled model "
        "is being served a feature it was never fitted on"
    )

    # And specifically not the configured-chemistry value.
    at_reaction = serving_tm(primer, Na_c=50.0, Mg_c=10.0)
    assert abs(served_value - at_reaction) > 1.0, (
        "the feature is being computed at a run's configured salt; that is "
        "train/serve skew, not a fix"
    )
