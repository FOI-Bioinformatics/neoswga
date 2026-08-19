"""`report` and `interpret` graded the same result differently, and both wrongly.

Two independent threshold tables existed. Enrichment "excellent" was 500x in
`report/quality.py` and 200x in `results_interpreter.py`, so the two commands
could hand a user different verdicts on one run. Uniformity disagreed more
sharply: one graded `1 - Gini` against 0.85, the other Gini against 0.30.

Both were also miscalibrated against the repository's own benchmark suite. Every
published set -- including `Prev06`, the best in the suite at 120-fold
enrichment and 72% of reads on target -- graded "poor" or "critical" on
uniformity. Nothing that has ever been published reached "good".

A scale on which no achievable design scores well is a constant with extra
steps. These tests pin the calibration against the eleven published sets that
carry both a Gini and an outcome, and pin the two commands to one source.
"""

import json
from pathlib import Path

import pytest

from neoswga.core import quality_thresholds as T

DATA = Path(__file__).parent / "validation" / "data"


def _published_sets():
    """Every benchmark set carrying a published Gini and an outcome."""
    out = []
    for name in ("clarke_2017_wolbachia", "dwivedi_yu_2023_prevotella"):
        for key, entry in json.loads((DATA / f"{name}.json").read_text())["sets"].items():
            gini = entry.get("published", {}).get("fg_gini")
            if gini is not None:
                out.append((key, gini, entry["outcome"]))
    return out


def _grade(value, thresholds, lower_is_better=False):
    for label in T.RATING_ORDER:
        if (value <= thresholds[label]) if lower_is_better else (value >= thresholds[label]):
            return label
    return "critical"


# ----------------------------------------------------------------------
# The scale spans what is achievable
# ----------------------------------------------------------------------


def test_the_best_published_set_is_not_graded_poor():
    """`Prev06`: 120-fold enrichment, 72% of reads on target, and the report
    called its binding evenness poor."""
    grade = _grade(0.552, T.GINI, lower_is_better=True)

    assert grade in ("excellent", "good"), grade


def test_the_grades_span_the_published_distribution():
    """Not every set should land in one bucket, in either direction."""
    grades = {_grade(gini, T.GINI, lower_is_better=True) for _, gini, _ in _published_sets()}

    assert len(grades) >= 3, grades


def test_the_worst_published_set_grades_worst():
    """The Leichty set has the worst evenness of the eleven (0.712) and Clarke
    reported it did not improve sequencing efficiency."""
    assert _grade(0.712, T.GINI, lower_is_better=True) in ("poor", "critical")


@pytest.mark.parametrize("fold,expected", [(120.0, "excellent"), (96.0, "good"), (2.07, "poor")])
def test_enrichment_thresholds_sit_inside_measured_outcomes(fold, expected):
    """The old bars -- 500x and 200x -- were above every result ever published.
    The best measured is 120x; 2x designs failed."""
    assert _grade(fold, T.ENRICHMENT) == expected


# ----------------------------------------------------------------------
# One source
# ----------------------------------------------------------------------


def test_report_and_interpret_share_the_enrichment_scale():
    from neoswga.core.report.quality import ENRICHMENT_THRESHOLDS as report_side
    from neoswga.core.results_interpreter import ENRICHMENT_THRESHOLDS as interpret_side

    assert set(report_side.values()) == set(interpret_side.values())


def test_report_and_interpret_share_the_uniformity_scale():
    """One reports `1 - Gini`, the other Gini, so the check is that they invert
    to each other rather than that they are equal."""
    from neoswga.core.report.quality import UNIFORMITY_THRESHOLDS as report_side
    from neoswga.core.results_interpreter import UNIFORMITY_THRESHOLDS as interpret_side

    interpret_by_label = {rating.value.lower(): v for rating, v in interpret_side.items()}
    for label, score in report_side.items():
        assert score == pytest.approx(1.0 - interpret_by_label[label]), label


def test_neither_module_defines_its_own_numbers():
    """A hand-written copy is how the two came to disagree by 2.5x."""
    import ast

    for module in ("core/report/quality.py", "core/results_interpreter.py"):
        source = (Path(__file__).parent.parent / "neoswga" / module).read_text()
        tree = ast.parse(source)
        for node in ast.walk(tree):
            if not isinstance(node, ast.Assign):
                continue
            names = {t.id for t in node.targets if isinstance(t, ast.Name)}
            if not names & {
                "ENRICHMENT_THRESHOLDS",
                "UNIFORMITY_THRESHOLDS",
                "COVERAGE_THRESHOLDS",
            }:
                continue
            assert not isinstance(node.value, ast.Dict) or not any(
                isinstance(v, ast.Constant) and isinstance(v.value, float)
                for v in node.value.values
            ), f"{module} restates threshold values; derive them from quality_thresholds"


def test_the_gini_grade_is_documented_as_descriptive():
    """It does not separate winners from losers: `TmL/Even` (effective) and
    `TmH/Even` (ineffective) have identical Gini. Anyone reading the grade as a
    prediction is being misled, so the module has to say so."""
    assert "not predictive" in (T.__doc__ or "")
