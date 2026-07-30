"""Guard: the `--enable-qa` filtering path must actually run.

`pipeline_qa_integration` consumed `dimer_network_analyzer` through an API that does
not exist, and nothing exercised it:

  * ``analyze_primer_set`` returns ``Dict[str, PrimerDimerProfile]``, but the code did
    ``for profile in profiles`` -- iterating a dict yields its *keys*, i.e. primer
    strings, so the very next attribute access failed.
  * ``PrimerDimerProfile`` has ``num_interactions``/``is_hub``; the code read
    ``profile.degree``.
  * ``DimerNetworkMetrics`` has no ``mean_degree``; the verbose branch read it.

Any one of these raised ``AttributeError`` on the first primer, so `--enable-qa` was
dead on arrival. These tests exercise the real path rather than asserting on source.
"""

import pandas as pd
import pytest

from neoswga.core.pipeline_qa_integration import (
    QAAwareOptimizer,
    QAFilterConfig,
    apply_post_step2_qa_filter,
)


@pytest.fixture
def step2_df():
    return pd.DataFrame(
        {
            "primer": [
                "ATCGATCGATCG",
                "GCTAGCTAGCTA",
                "AAATTTAAATTT",
                "GCGCGCGCGCGC",
                "ACGTACGTACGT",
                "TTTTAAAATTTT",
            ],
            "fg_freq": [1e-4] * 6,
            "bg_freq": [1e-7] * 6,
            "gini": [0.3] * 6,
            "tm": [30.0] * 6,
        }
    )


def test_dimer_profile_fields_are_what_the_qa_code_expects():
    """Pin the contract the QA integration relies on."""
    from neoswga.core.dimer_network_analyzer import DimerNetworkMetrics, PrimerDimerProfile

    profile_fields = set(PrimerDimerProfile.__dataclass_fields__)
    assert "num_interactions" in profile_fields
    assert "primer" in profile_fields
    assert "degree" not in profile_fields, (
        "PrimerDimerProfile gained a `degree` field - reconcile pipeline_qa_integration"
    )
    assert "mean_degree" not in set(DimerNetworkMetrics.__dataclass_fields__)


def test_analyze_primer_set_returns_a_mapping():
    """The QA code must iterate .values(), not the dict itself."""
    from neoswga.core.dimer_network_analyzer import create_dimer_network_analyzer

    analyzer = create_dimer_network_analyzer("moderate", None)
    _metrics, profiles, _matrix = analyzer.analyze_primer_set(
        ["ATCGATCGATCG", "GCTAGCTAGCTA"]
    )
    assert isinstance(profiles, dict)
    for key, value in profiles.items():
        assert isinstance(key, str), "iterating this dict directly yields primer strings"
        assert hasattr(value, "num_interactions")


def test_qa_filter_runs_with_dimer_network_enabled(step2_df):
    """Regression: this raised AttributeError before the fix."""
    config = QAFilterConfig(enable_dimer_network=True)
    result = apply_post_step2_qa_filter(step2_df, config=config, verbose=False)

    assert result is not None
    assert hasattr(result, "filtered_df")
    assert len(result.filtered_df) <= len(step2_df)


def test_qa_filter_verbose_branch_runs(step2_df, caplog):
    """The verbose branch read the non-existent metrics.mean_degree."""
    config = QAFilterConfig(enable_dimer_network=True)
    result = apply_post_step2_qa_filter(step2_df, config=config, verbose=True)
    assert result is not None


def test_qa_filter_marks_hub_primers_rather_than_crashing(step2_df):
    """A very low hub threshold must filter, not raise."""
    config = QAFilterConfig(enable_dimer_network=True, max_hub_degree=0)
    result = apply_post_step2_qa_filter(step2_df, config=config, verbose=False)
    assert result is not None


def test_qa_aware_optimizer_prefilter_runs(step2_df):
    """The second broken site lived in QAAwareOptimizer's candidate prefilter."""
    from neoswga.core.reaction_conditions import get_standard_conditions

    candidates = step2_df["primer"].tolist()
    optimizer = QAAwareOptimizer(
        config=QAFilterConfig(enable_dimer_network=True),
        conditions=get_standard_conditions(),
    )

    # verbose=True is deliberate: the mean_degree crash lived in that branch.
    filtered, metadata = optimizer.prefilter_candidates(candidates, verbose=True)

    assert isinstance(filtered, list)
    assert set(filtered) <= set(candidates)
    assert metadata["original_count"] == len(candidates)


def test_qa_aware_optimizer_prefilter_handles_empty_input():
    from neoswga.core.reaction_conditions import get_standard_conditions

    optimizer = QAAwareOptimizer(
        config=QAFilterConfig(enable_dimer_network=True),
        conditions=get_standard_conditions(),
    )
    filtered, metadata = optimizer.prefilter_candidates([], verbose=True)
    assert filtered == []
    assert metadata["original_count"] == 0
