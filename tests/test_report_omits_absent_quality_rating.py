"""The report omits the per-primer quality rating rather than faking one.

Finding D0. Retiring the amplification model on 2026-09-05 left `amp_pred` at
0.0 for every primer, and three rendering sites substitute 0.5 when it is
falsy. `int(0.5 * 5 + 0.5)` is 3, so every primer in every report rendered the
same three stars. A per-primer column that cannot vary is worse than no column;
`a2a58c0` fixed that once already, from the other end.
"""

import pytest

from neoswga.core.report import executive_summary, technical_report
from neoswga.core.report.metrics import PrimerMetrics, amp_pred_is_available


def _primer(sequence: str, amp_pred: float = 0.0) -> PrimerMetrics:
    return PrimerMetrics(
        sequence=sequence,
        length=len(sequence),
        gc_content=0.5,
        tm=42.0,
        fg_freq=1e-5,
        bg_freq=1e-7,
        fg_sites=20,
        bg_sites=2,
        gini=0.4,
        specificity=10.0,
        amp_pred=amp_pred,
        strand_ratio=1.0,
    )


def test_availability_helper_reports_absence_and_presence():
    absent = [_primer("AAACCCGGGTTT"), _primer("ACCCGGGTTTAA")]
    present = [_primer("AAACCCGGGTTT"), _primer("ACCCGGGTTTAA", amp_pred=0.7)]
    assert amp_pred_is_available(absent) is False
    assert amp_pred_is_available(present) is True
    assert amp_pred_is_available([]) is False


def test_primer_row_has_no_quality_cell_when_no_primer_carries_a_score():
    row = executive_summary._format_primer_row(1, _primer("AAACCCGGGTTT"), show_quality=False)
    assert "quality-stars" not in row
    assert "★" not in row
    assert "☆" not in row
    assert row.count("<td") == 6


def test_primer_row_keeps_the_quality_cell_when_a_score_exists():
    row = executive_summary._format_primer_row(
        1, _primer("AAACCCGGGTTT", amp_pred=0.8), show_quality=True
    )
    assert "quality-stars" in row
    assert row.count("<td") == 7


def test_contribution_score_renormalises_when_no_score_exists():
    """Dropping the 0.2 amplification term must not shrink every primer's score."""
    p = _primer("AAACCCGGGTTT")
    p.specificity = 100.0  # specificity_score 1.0
    p.gini = 0.0  # uniformity_score 1.0
    p.strand_ratio = 1.0  # strand_score 1.0
    contribution = technical_report._contribution_score(p, amp_available=False)
    assert contribution == pytest.approx(1.0)


def test_contribution_score_keeps_four_terms_when_a_score_exists():
    p = _primer("AAACCCGGGTTT", amp_pred=1.0)
    p.specificity = 100.0
    p.gini = 0.0
    p.strand_ratio = 1.0
    contribution = technical_report._contribution_score(p, amp_available=True)
    assert contribution == pytest.approx(1.0)
