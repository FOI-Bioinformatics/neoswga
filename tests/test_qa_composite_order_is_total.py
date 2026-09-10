"""`score --enable-qa` must not reintroduce the unstable sort.

Audit finding D1d. `pipeline_qa_integration.py:718` sorted on composite_score
alone with the default quicksort and no secondary key. Measured on the
449-candidate E. coli pool, two input orders gave 65 differing positions, the
first at rank 14 where two primers swap, and 100 of the 449 rows share a
composite value with another row.

That is the same defect `order_step3_rows` exists to remove, on a different
path. The optimizer is order-sensitive: on the E. coli pool at target size 24,
dominating-set returned a set with a Jaccard of 0.600 against the as-written
order when the candidates were reversed.
"""

import pandas as pd
import pytest


def _tied_pool(n=200):
    """A pool with heavy composite ties, which is the real shape: 100 of the
    449 E. coli rows share a composite value with another row."""
    return pd.DataFrame(
        {
            "primer": [f"{i:04d}ACGTAC" for i in range(n)],
            # Five distinct values across 200 rows: ties everywhere.
            "composite_score": [round(0.90 + 0.01 * (i % 5), 2) for i in range(n)],
        }
    )


def test_the_pool_really_is_mostly_ties():
    """Guard the guard: with distinct values the sort is total already."""
    df = _tied_pool()
    assert df["composite_score"].duplicated(keep=False).sum() > 0.9 * len(df)


def test_the_order_does_not_depend_on_how_the_rows_arrived():
    """The regression."""
    from neoswga.core.pipeline_qa_integration import order_by_composite_score

    df = _tied_pool()
    a = order_by_composite_score(df.sample(frac=1, random_state=0))["primer"].tolist()
    b = order_by_composite_score(df.sample(frac=1, random_state=1))["primer"].tolist()
    assert a == b


def test_the_composite_still_leads():
    from neoswga.core.pipeline_qa_integration import order_by_composite_score

    ordered = order_by_composite_score(_tied_pool())
    assert ordered["composite_score"].is_monotonic_decreasing


def test_ties_are_broken_by_the_primer_sequence():
    from neoswga.core.pipeline_qa_integration import order_by_composite_score

    ordered = order_by_composite_score(_tied_pool())
    for _score, group in ordered.groupby("composite_score", sort=False):
        primers = group["primer"].tolist()
        assert primers == sorted(primers)


def test_the_qa_naming_of_the_primer_column_also_works():
    """`primer_column` prefers 'primer' and falls back to 'seq'."""
    from neoswga.core.pipeline_qa_integration import order_by_composite_score

    df = _tied_pool().rename(columns={"primer": "seq"})
    ordered = order_by_composite_score(df)
    for _score, group in ordered.groupby("composite_score", sort=False):
        assert group["seq"].tolist() == sorted(group["seq"].tolist())


def test_an_empty_frame_is_not_an_error():
    from neoswga.core.pipeline_qa_integration import order_by_composite_score

    assert len(order_by_composite_score(_tied_pool().iloc[:0])) == 0


# ---------------------------------------------------------------------------
# 55% of the declared weight could not move
# ---------------------------------------------------------------------------


def _scorer():
    from neoswga.core.integrated_quality_scorer import create_quality_scorer

    return create_quality_scorer("moderate")


def test_the_declared_weights_still_sum_to_one():
    """Guard the guard. The renormalization is only meaningful against a
    normalized declaration."""
    assert sum(_scorer().weights.values()) == pytest.approx(1.0)


def test_dimer_and_strand_are_the_larger_half_of_the_declaration():
    """The size of the problem: 0.35 + 0.20 of a declared 1.0."""
    weights = _scorer().weights
    assert weights["dimer"] + weights["strand_bias"] == pytest.approx(0.55)


def test_a_composite_over_a_subset_renormalizes_to_that_subset():
    scorer = _scorer()
    scores = {
        "three_prime": 0.5,
        "thermodynamics": 0.5,
        "complexity": 0.5,
        "dimer": 1.0,
        "strand_bias": 1.0,
    }
    measured = ("three_prime", "thermodynamics", "complexity")
    assert scorer._composite(scores, measured) == pytest.approx(0.5)


def test_a_composite_over_every_component_matches_the_plain_weighted_sum():
    scorer = _scorer()
    scores = {
        "three_prime": 0.4,
        "thermodynamics": 0.6,
        "complexity": 0.8,
        "dimer": 0.2,
        "strand_bias": 1.0,
    }
    expected = sum(scorer.weights[name] * value for name, value in scores.items())
    assert scorer._composite(scores, scores.keys()) == pytest.approx(expected)


def test_scoring_without_binding_sites_excludes_strand_and_dimer():
    """The regression. Without binding sites, strand and dimer are constants,
    so they must not carry 55% of the weight."""
    scorer = _scorer()
    score = scorer.score_primer("ACGTACGTACGT")
    assert set(score.measured_components) == {"three_prime", "thermodynamics", "complexity"}


def test_the_unmeasured_components_no_longer_prop_the_score_up():
    """A primer that is poor on every component that CAN be measured must not
    score 0.55 + 0.45 * poor."""
    scorer = _scorer()
    score = scorer.score_primer("ACGTACGTACGT")
    measured_only = sum(
        scorer.weights[name] * getattr(score, _FIELD[name]) for name in score.measured_components
    ) / sum(scorer.weights[name] for name in score.measured_components)
    assert score.overall_score == pytest.approx(measured_only)


_FIELD = {
    "three_prime": "three_prime_score",
    "thermodynamics": "thermo_score",
    "complexity": "complexity_score",
    "dimer": "dimer_score",
    "strand_bias": "strand_bias_score",
}


def test_a_set_analysis_measures_dimer_and_says_so():
    """`analyze_primer_set` computes a real dimer score, so the composite must
    include it again."""
    scorer = _scorer()
    primers = ["ACGTACGTACGT", "TTGACCATGACC", "GGCATTACGATC"]
    primer_scores, _set_score = scorer.analyze_primer_set(primers, verbose=False)
    for score in primer_scores:
        assert "dimer" in score.measured_components


def test_an_unmeasurable_composite_is_not_a_division_by_zero():
    scorer = _scorer()
    assert scorer._composite({"dimer": 1.0}, ()) == 0.0
