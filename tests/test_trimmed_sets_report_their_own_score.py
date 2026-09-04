"""A trimmed set must not carry the untrimmed set's score, and 0 means 0.

Two small defects in the `--minimize-primers` path.

`_minimize_primer_count` returned `replace(result, primers=..., metrics=...)`
and left `score` untouched, so a set trimmed from four primers to one still
reported the four-primer score -- written to the `score` column of
step4_improved_df.csv and to the summary JSON. The number described a set that
no longer existed.

And `float(kwargs.pop("target_coverage", 0.70) or 0.70)` treated an explicit
zero as "unset", so `--target-coverage 0`, which reasonably means "trim as far
as you can", silently became the 0.70 default.
"""

from dataclasses import replace

import pytest

from neoswga.core.base_optimizer import OptimizationResult, OptimizationStatus, PrimerSetMetrics

PRIMERS = ["AAACCCGGGT", "ACCCGGGTTT", "AAGGCCTTAC", "GTAAGGCCTT"]


class _Optimizer:
    """Coverage falls only when the pool drops below two primers, so the
    minimiser has an unambiguous trim to make."""

    def compute_metrics(self, primers):
        coverage = 0.95 if len(primers) >= 2 else 0.10
        return replace(PrimerSetMetrics.empty(), fg_coverage=coverage)


def _result(score=42.0):
    return OptimizationResult(
        primers=tuple(PRIMERS),
        score=score,
        metrics=replace(PrimerSetMetrics.empty(), fg_coverage=0.99),
        status=OptimizationStatus.SUCCESS,
        optimizer_name="hybrid",
        iterations=1,
    )


def test_a_trimmed_set_does_not_report_the_untrimmed_score():
    from neoswga.core.unified_optimizer import _minimize_primer_count

    original = _result(score=42.0)
    trimmed = _minimize_primer_count(
        result=original, optimizer=_Optimizer(), target_coverage=0.70, verbose=False
    )

    assert len(trimmed.primers) < len(original.primers), "nothing was trimmed"
    assert trimmed.score != original.score, (
        f"the trimmed set still reports {trimmed.score}, the score of the "
        f"{len(original.primers)}-primer set it no longer is"
    )
    assert trimmed.score == pytest.approx(trimmed.metrics.normalized_score())


def test_an_untrimmed_set_is_returned_unchanged():
    """The guard in the other direction: if nothing is trimmed, nothing moves."""
    from neoswga.core.unified_optimizer import _minimize_primer_count

    class _NeverEnough:
        def compute_metrics(self, primers):
            return replace(PrimerSetMetrics.empty(), fg_coverage=0.10)

    original = _result(score=42.0)
    same = _minimize_primer_count(
        result=original, optimizer=_NeverEnough(), target_coverage=0.70, verbose=False
    )

    assert same.primers == original.primers
    assert same.score == original.score


@pytest.mark.parametrize(
    "supplied,expected",
    [(None, 0.70), (0, 0.0), (0.0, 0.0), (0.9, 0.9)],
)
def test_target_coverage_zero_is_not_treated_as_unset(supplied, expected):
    """Only a missing value may fall back to the default."""
    from neoswga.core.unified_optimizer import _resolve_target_coverage

    kwargs = {} if supplied is None else {"target_coverage": supplied}
    assert _resolve_target_coverage(kwargs) == expected
    assert "target_coverage" not in kwargs, "the key must be popped, not left to travel twice"
