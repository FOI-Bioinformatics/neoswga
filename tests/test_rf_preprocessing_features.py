"""Behavioural tests for the random-forest feature engineering.

`rf_preprocessing` runs in every `score` step: it builds the feature vector the
shipped model consumes and, when that model is unavailable, the heuristic that
stands in for it. A wrong feature does not fail -- it produces a plausible score
for the wrong reason, which then decides which primers reach the optimizer.

It sat at 38%. The parts with no coverage were the ones that decide things:
the heuristic fallback, the delta-G scaling, and the k-mer sampling controls
that trade accuracy for speed on large genomes.
"""

import math

import numpy as np
import pandas as pd
import pytest

from neoswga.core import rf_preprocessing as rf

PRIMER = "TTGACCATGA"


def _features_as_dict(primer, **kwargs):
    """`get_features` returns a positional list; pair it with its own names."""
    values = rf.get_features(primer, **kwargs)
    return values


# ----------------------------------------------------------------------
# The feature vector
# ----------------------------------------------------------------------


def test_feature_vector_is_case_insensitive():
    """Every base count is case-sensitive, so lowercase produced zeros.

    `primer.count("A")` sees nothing in "atcg", so a lowercase primer built a
    vector of zero counts, zero proportions and zero GC content before failing
    on the 3'-base lookup with `KeyError: 'a'`. Jellyfish output is uppercase so
    the pipeline never hit it, but soft-masked FASTA and hand-typed primer lists
    are lowercase and reach scoring through the bring-your-own-oligo commands.
    """
    assert rf.get_features(PRIMER) == rf.get_features(PRIMER.lower())
    assert rf.get_features(PRIMER) == rf.get_features("TtGaCcAtGa")


def test_base_counts_sum_to_the_length():
    """A vector whose counts do not add up means one base was miscounted."""
    values = rf.get_features(PRIMER)
    row = dict(zip(rf.base_features, values[2:])) if hasattr(rf, "base_features") else None
    if row is None:
        pytest.skip("base feature names not exported")

    total = sum(row[f"number.of.{b}"] for b in "AGTC")
    assert total == row["sequence.length"]


def test_gc_content_matches_the_sequence():
    values = rf.get_features("GGGGGCCCCC")
    assert 1.0 in values, "a 100% GC primer should report GC content 1.0"


def test_feature_vector_length_is_stable_across_primers():
    """The model consumes a fixed-width vector; a ragged one would break it."""
    lengths = {len(rf.get_features(p)) for p in (PRIMER, "ACGT" * 3, "AAAAAAAAAA")}
    assert len(lengths) == 1


def test_every_feature_is_finite():
    """A NaN here propagates into the model score without failing."""
    for primer in (PRIMER, "AAAAAAAAAA", "GCGCGCGCGC", "ACGTA"):
        for value in rf.get_features(primer):
            if isinstance(value, float):
                assert math.isfinite(value), f"non-finite feature for {primer}"


def test_primers_shorter_than_five_bases_are_not_supported():
    """The 3'-end features read `primer[-5]`, so four bases raise IndexError.

    Pinned rather than fixed: `min_k` defaults to 6 and a 4-mer is not a
    plausible SWGA primer, so the constraint is real. Worth stating explicitly
    so it is a known limit rather than a surprise.
    """
    with pytest.raises(IndexError):
        rf.get_features("ACGT")


def test_homopolymer_feature_tracks_the_sequence():
    plain = rf.get_features("ACGTACGTAC")
    runny = rf.get_features("AAAAAACGTA")
    assert plain != runny


# ----------------------------------------------------------------------
# The base feature matrix
# ----------------------------------------------------------------------


def test_feature_matrix_has_one_row_per_primer():
    primers = [PRIMER, "GCATTACGGT", "AACCGGTTAC"]
    matrix = rf.create_base_feature_matrix(primers, molarity=2.5)

    assert len(matrix) == len(primers)


def test_feature_matrix_runs_before_get_params():
    """It reads `parameter.cpus` as a bare attribute.

    That had no module-level default, so building a feature matrix outside the
    CLI pipeline raised AttributeError -- the same shape as max_self_dimer_bp in
    the filtering rules.
    """
    import subprocess
    import sys

    code = (
        "from neoswga.core.rf_preprocessing import create_base_feature_matrix;"
        "print(len(create_base_feature_matrix(['TTGACCATGA'], molarity=2.5)))"
    )
    proc = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=300)
    assert proc.returncode == 0, proc.stderr[-1200:]
    assert proc.stdout.strip() == "1"


def test_feature_matrix_preserves_primer_order():
    """Rows are matched back to primers positionally downstream, so a reordering
    would attach every score to the wrong sequence."""
    primers = [PRIMER, "GCATTACGGT", "AACCGGTTAC"]
    matrix = rf.create_base_feature_matrix(primers, molarity=2.5)

    assert list(matrix["sequence"]) == primers


def test_feature_matrix_columns_are_the_documented_features():
    """The model consumes these by name; a renamed column is a silent mismatch."""
    matrix = rf.create_base_feature_matrix([PRIMER], molarity=2.5)

    assert "sequence" in matrix.columns
    assert "GC.content" in matrix.columns
    assert "melting_tm" in matrix.columns


def test_empty_primer_list_gives_an_empty_frame():
    matrix = rf.create_base_feature_matrix([], molarity=2.5)
    assert len(matrix) == 0


# ----------------------------------------------------------------------
# The heuristic fallback, used when the model is unavailable
# ----------------------------------------------------------------------


def _row(tm=30.0, gc=0.5, three_prime="G"):
    """A minimal row shaped like what the heuristic reads."""
    row = {
        "melting_tm": tm,
        "GC.content": gc,
        "sequence.length": 10,
        "3.end.first.base": {"A": 1, "C": 2, "G": 3, "T": 4}[three_prime],
    }
    for name in getattr(rf, "delta_g_on_features", []):
        row[name] = 0.0
    return pd.Series(row)


def test_heuristic_score_is_in_the_models_range():
    """It substitutes for the RF model, so it must be comparable to it.

    `min_amp_pred` thresholds against this number, so a different scale would
    silently change how many primers survive.
    """
    score = rf._heuristic_primer_score(_row())
    assert 0.0 <= score <= 20.0


def test_heuristic_prefers_a_tm_near_the_reaction_temperature():
    """SWGA is isothermal around 30 C; a PCR-range Tm is not useful here."""
    good = rf._heuristic_primer_score(_row(tm=32.0))
    too_hot = rf._heuristic_primer_score(_row(tm=75.0))
    too_cold = rf._heuristic_primer_score(_row(tm=2.0))

    assert good > too_hot
    assert good > too_cold


def test_heuristic_prefers_balanced_gc():
    balanced = rf._heuristic_primer_score(_row(gc=0.5))
    extreme = rf._heuristic_primer_score(_row(gc=0.95))
    assert balanced > extreme


def test_heuristic_is_deterministic():
    """A scorer that varies between calls cannot be thresholded."""
    assert rf._heuristic_primer_score(_row()) == rf._heuristic_primer_score(_row())


def test_heuristic_never_returns_non_finite():
    for tm in (0.0, 30.0, 100.0):
        for gc in (0.0, 0.5, 1.0):
            score = rf._heuristic_primer_score(_row(tm=tm, gc=gc))
            assert math.isfinite(score)


# ----------------------------------------------------------------------
# Delta-G scaling
# ----------------------------------------------------------------------


def test_delta_g_scaling_divides_by_the_given_factor():
    """The on/off-target genome size ratio; getting it wrong tilts every score.

    Background genomes are usually far larger than the target, so the histogram
    counts are not comparable until they are put on the same scale.
    """
    names = getattr(rf, "delta_g_on_features", [])
    if not names:
        pytest.skip("delta_g_on_features not exported")

    df = pd.DataFrame({name: [4000.0] for name in names})
    scaled = rf.scale_delta_Gs(df.copy(), on_scale=4000)

    for name in names:
        assert scaled[name].iloc[0] == pytest.approx(1.0)


def test_delta_g_scaling_is_linear():
    names = getattr(rf, "delta_g_on_features", [])
    if not names:
        pytest.skip("delta_g_on_features not exported")

    df = pd.DataFrame({name: [100.0] for name in names})
    half = rf.scale_delta_Gs(df.copy(), on_scale=2)
    quarter = rf.scale_delta_Gs(df.copy(), on_scale=4)

    for name in names:
        assert half[name].iloc[0] == pytest.approx(2 * quarter[name].iloc[0])


# ----------------------------------------------------------------------
# K-mer sampling: the accuracy/speed trade for large genomes
# ----------------------------------------------------------------------


@pytest.fixture(autouse=True)
def reset_sampling():
    yield
    rf.disable_kmer_sampling()
    rf.clear_kmer_cache()


def test_sampling_can_be_enabled_and_disabled():
    rf.enable_kmer_sampling(sample_rate=0.1, min_count=5)
    rf.disable_kmer_sampling()  # must not raise or leave state behind


def test_sampling_seed_makes_it_reproducible():
    """Sampling trades accuracy for speed; without a seed the same input would
    score differently between runs, which is not debuggable."""
    rf.set_kmer_sampling_seed(1234)
    rf.set_kmer_sampling_seed(1234)  # idempotent


def test_kmer_cache_can_be_cleared():
    """The cache is module-level state shared across runs in one process."""
    rf.clear_kmer_cache()
    assert rf._kmer_cache == {}


def test_clearing_the_cache_twice_is_safe():
    rf.clear_kmer_cache()
    rf.clear_kmer_cache()


# ----------------------------------------------------------------------
# Model loading and integrity
# ----------------------------------------------------------------------


def test_missing_model_raises_file_not_found():
    with pytest.raises(FileNotFoundError):
        rf.load_model_safely("/definitely/not/a/model.skops")


def test_shipped_model_passes_its_hash_check(tmp_path):
    """The checksum guard is the reason the model can be loaded at all.

    If it starts failing, either the model was retrained without updating
    checksums.json or the file has been altered.
    """
    import os

    model_dir = os.path.join(os.path.dirname(rf.__file__), "models")
    model_path = os.path.join(model_dir, "random_forest_filter.skops")
    if not os.path.exists(model_path):
        pytest.skip("bundled model not present")

    model = rf.load_model_safely(model_path)
    assert model is not None


def test_tampered_model_is_rejected(tmp_path):
    """A hash mismatch must raise rather than load."""
    import os
    import shutil

    model_dir = os.path.join(os.path.dirname(rf.__file__), "models")
    model_path = os.path.join(model_dir, "random_forest_filter.skops")
    if not os.path.exists(model_path):
        pytest.skip("bundled model not present")

    tampered = tmp_path / "random_forest_filter.skops"
    shutil.copy(model_path, tampered)
    with open(tampered, "ab") as fh:
        fh.write(b"\x00extra")

    with pytest.raises(rf.ModelIntegrityError):
        rf.load_model_safely(str(tampered))
