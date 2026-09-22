"""The panel assessment must judge the panel, not only measure it.

`evaluate_panel` is described as the authoritative record of whether a panel is
acceptable. It measured carefully and enforced almost nothing: violations came
solely from an optional objective, so with none supplied a one-primer panel
qualified against a request for twelve.

    evaluate_panel(request_for_12, ["ACGT..."], metrics).qualified  -> True

It also reported geometric `fg_coverage` while the live acceptance path selects
on EFFECTIVE coverage, so the record answered a different question from the
search that produced it.

Three things are added here, each already enforced somewhere else in the
codebase, which is the point: they are being brought into one record rather
than invented.

- **Requested size**, enforced today in `optimization_service`'s `assess`
  closure and again in `base_optimizer.validate`.
- **Delivered-panel dimers**, enforced today during the search in
  `optimization_service.panel_violations` and reported again, at warning level
  only, after it.
- **Effective coverage**, which the objective path already measures and the
  record omitted.

What is deliberately NOT added is candidate composition and QC. In this
codebase those are admission rules applied during `filter`, not properties of a
delivered panel, and no panel-level path re-checks them. Faulting the
assessment for that would be a category error.

This increment completes the record. Making it the SOLE gate, and deleting the
paths it subsumes, is the next one: several of those live on the search's
stopping rule, so removing them moves delivered panels and needs measuring.
"""

import pytest

from neoswga.core.design_request import resolve_design_request
from neoswga.core.panel_evaluation import evaluate_panel

MAPPING = {
    "fg_genomes": ["a.fna"],
    "fg_prefixes": ["a"],
    "polymerase": "phi29",
    "reaction_temp": 30.0,
    "min_k": 12,
    "max_k": 12,
    "num_primers": 12,
    "max_dimer_bp": 3,
}

CLEAN = ["ACCACAGATAGC", "GTTGTAGATGGA", "ATCAGCAGACCA"]
# 10 of 12 bases complementary, the pair this project's dimer tests use.
DIMERISING = ["TGTACTGCCAAG", "CTTGGCAGTACA"]


class Metrics:
    """The evaluator's own measurements, as the assessment receives them."""

    def __init__(self, **kwargs):
        self.fg_coverage = kwargs.pop("fg_coverage", 0.5)
        self.bg_coverage = kwargs.pop("bg_coverage", 0.01)
        self.total_fg_sites = kwargs.pop("total_fg_sites", 100)
        self.total_bg_sites = kwargs.pop("total_bg_sites", 5)
        self.mean_gap = kwargs.pop("mean_gap", 1000.0)
        self.max_gap = kwargs.pop("max_gap", 5000.0)
        self.gap_gini = kwargs.pop("gap_gini", 0.3)
        self.mean_tm = kwargs.pop("mean_tm", 35.0)
        self.selectivity_ratio = kwargs.pop("selectivity_ratio", 20.0)
        self.per_target_coverage = kwargs.pop("per_target_coverage", {})
        self.effective_fg_coverage = kwargs.pop("effective_fg_coverage", 0.45)
        for key, value in kwargs.items():
            setattr(self, key, value)


def request_for(size=12, **overrides):
    return resolve_design_request({**MAPPING, "num_primers": size, **overrides})


# ---------------------------------------------------------------------------
# Requested size
# ---------------------------------------------------------------------------


def test_a_short_panel_does_not_qualify():
    """The headline. One primer against a request for twelve qualified."""
    assessment = evaluate_panel(request_for(12), CLEAN[:1], Metrics())

    assert not assessment.qualified
    assert any("size" in v.lower() for v in assessment.violations), assessment.violations


def test_a_panel_of_the_requested_size_qualifies():
    """Guard the guard: an assessment that refuses everything is no better."""
    assessment = evaluate_panel(request_for(3), CLEAN, Metrics())

    assert assessment.qualified, assessment.violations


def test_the_violation_names_both_counts():
    """A reader has to see what was asked for and what arrived."""
    assessment = evaluate_panel(request_for(12), CLEAN, Metrics())

    joined = " ".join(assessment.violations)
    assert "12" in joined and "3" in joined, joined


def test_a_larger_panel_than_requested_is_not_a_violation():
    """`num_primers` is a request and the delivered panel is never larger, so
    this cannot arise from selection. Asserted so that a future change which
    makes it possible is a deliberate one rather than a silent rejection."""
    assessment = evaluate_panel(request_for(2), CLEAN, Metrics())

    assert assessment.qualified, assessment.violations


# ---------------------------------------------------------------------------
# Delivered-panel dimers
# ---------------------------------------------------------------------------


def test_a_dimerising_pair_does_not_qualify():
    """Enforced hard during the search and, after it, only at warning level."""
    assessment = evaluate_panel(request_for(2), DIMERISING, Metrics())

    assert not assessment.qualified
    assert any("dimer" in v.lower() for v in assessment.violations), assessment.violations


def test_the_dimer_violation_names_the_pair():
    assessment = evaluate_panel(request_for(2), DIMERISING, Metrics())

    joined = " ".join(assessment.violations)
    assert DIMERISING[0] in joined or DIMERISING[1] in joined, joined


def test_a_clean_panel_carries_no_dimer_violation():
    assessment = evaluate_panel(request_for(3), CLEAN, Metrics())

    assert not any("dimer" in v.lower() for v in assessment.violations)


# ---------------------------------------------------------------------------
# Effective coverage
# ---------------------------------------------------------------------------


def test_effective_coverage_is_recorded():
    """The record reported the geometric figure while the acceptance path
    selects on the effective one, so the two answered different questions."""
    assessment = evaluate_panel(request_for(3), CLEAN, Metrics(effective_fg_coverage=0.45))

    assert "effective_fg_coverage" in assessment.metrics
    assert assessment.metrics["effective_fg_coverage"].value == pytest.approx(0.45)


def test_an_unmeasured_effective_coverage_says_so_rather_than_reading_zero():
    """Absence and zero must not share a representation. Without reaction
    conditions there is no temperature at which to evaluate occupancy, so this
    is genuinely unavailable rather than measured as none."""
    metrics = Metrics()
    metrics.effective_fg_coverage = None

    assessment = evaluate_panel(request_for(3), CLEAN, metrics)

    measurement = assessment.metrics["effective_fg_coverage"]
    assert measurement.value is None
    assert measurement.unavailable


def test_geometric_coverage_is_still_reported_and_still_labelled():
    """Both are kept. Which one a reader is looking at must be unambiguous."""
    assessment = evaluate_panel(request_for(3), CLEAN, Metrics())

    assert assessment.metrics["fg_coverage"].value == pytest.approx(0.5)
    assert "reach" in assessment.metrics["fg_coverage"].basis


# ---------------------------------------------------------------------------
# What is deliberately unchanged
# ---------------------------------------------------------------------------


def test_a_non_finite_required_metric_still_fails_the_run():
    """NaN compares False against every threshold, so a panel carrying one
    passes no limit and fails none. That refusal predates this work."""
    from neoswga.core.exceptions import ModelEvaluationError

    with pytest.raises(ModelEvaluationError):
        evaluate_panel(request_for(3), CLEAN, Metrics(fg_coverage=float("nan")))


def test_configured_limits_still_come_from_the_objective():
    """Panel limits are the objective's job and stay there. This increment
    adds the three rules that were enforced elsewhere, not a fourth owner of
    the configured limits."""

    class Objective:
        def violations(self, primers):
            return ["selectivity density below minimum"]

    assessment = evaluate_panel(request_for(3), CLEAN, Metrics(), objective=Objective())

    assert not assessment.qualified
    assert "selectivity density below minimum" in assessment.violations
