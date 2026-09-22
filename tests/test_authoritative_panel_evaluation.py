"""One acceptance record, so a metric and a verdict cannot come from different arithmetic.

Task 5 of the 2026-09-21 valid-design plan. A panel is measured in several
places today and assembled differently at each: the optimizer holds
`PrimerSetMetrics`, the acceptance path builds an `AcceptanceReport`, the
summary JSON is written by a third piece of code and the report renders a
fourth. Nothing makes them agree, and a hybrid run already prints two coverage
figures that differ by construction.

What is pinned here is the shape of the one record, and three rules it carries:

- a required quantity that is NaN or infinite fails the run rather than being
  reported, because a non-finite number passes no threshold and fails none;
- a verified zero background is represented as a zero denominator, not divided
  by and not serialised as infinity;
- an unavailable quantity and a measured zero have different representations,
  because a panel that binds the host nowhere and a panel whose host index was
  never opened otherwise read the same.
"""

import json
import math

import pytest

from neoswga.core.design_request import resolve_design_request
from neoswga.core.exceptions import ModelEvaluationError
from neoswga.core.panel_evaluation import Measurement, PanelAssessment, evaluate_panel

# A pair whose worst complementary run is 3 bp, at the default
# `max_dimer_bp` of 3 and so within it. The previous fixture,
# ("AAAACCCCGGGG", "TTTTGGGGCCCC"), shares a 4 bp run and broke it.
#
# That did not matter while `evaluate_panel` enforced nothing but measurement
# validity. It does now: as of 2026-09-22 the assessment also checks the
# requested size and the delivered-panel dimer limit, both of which were
# already enforced elsewhere and neither of which reached this record. A
# fixture that violates the rule under test would make these tests assert the
# wrong thing.
PRIMERS = ("ACCACAGATAGC", "GTTGTAGATGGA")


class Metrics:
    """A stand-in for `PrimerSetMetrics`, so the record's own rules are tested.

    Deliberately not the production dataclass: this file is about what
    `evaluate_panel` does with measurements, and building a real one would
    couple these assertions to an evaluator that is not under test.
    """

    def __init__(self, **values):
        defaults = {
            "fg_coverage": 0.62,
            "bg_coverage": 0.04,
            "total_fg_sites": 240,
            "total_bg_sites": 31,
            "selectivity_ratio": 7.7,
            "mean_gap": 1800.0,
            "max_gap": 9100.0,
            "gap_gini": 0.41,
            "mean_tm": 33.5,
            "per_target_coverage": {},
        }
        defaults.update(values)
        for key, value in defaults.items():
            setattr(self, key, value)


def request(**overrides):
    params = {
        "fg_prefixes": ["target"],
        "fg_genomes": ["target.fna"],
        "polymerase": "phi29",
        "reaction_temp": 30.0,
        "min_k": 8,
        "max_k": 12,
        # Was 12 while the fixture delivered 2 and nothing checked. The
        # assessment now enforces the requested size, so the fixture has to
        # ask for what it hands over.
        "num_primers": len(PRIMERS),
    }
    params.update(overrides)
    return resolve_design_request(params)


class Objective:
    def __init__(self, violations=()):
        self._violations = tuple(violations)

    def violations(self, primers):
        return self._violations


# ---------------------------------------------------------------------------
# The record
# ---------------------------------------------------------------------------


def test_the_assessment_names_its_panel_and_its_request():
    assessment = evaluate_panel(request(), PRIMERS, Metrics())

    assert assessment.primers == PRIMERS
    assert assessment.request_hash == request().request_hash
    assert isinstance(assessment, PanelAssessment)


def test_every_metric_carries_its_units_and_what_it_was_computed_against():
    """A coverage figure without its reach carries almost no information.

    One saved 26-oligo panel reads 41.3 percent at 1 kb and 93.5 at 5 kb.
    """
    assessment = evaluate_panel(request(coverage_reach=5000), PRIMERS, Metrics())

    coverage = assessment.metrics["fg_coverage"]
    assert coverage.units == "fraction of target bases"
    assert "5000" in coverage.basis
    assert "denominator" in coverage.basis


def test_the_record_is_immutable():
    import dataclasses

    assessment = evaluate_panel(request(), PRIMERS, Metrics())
    with pytest.raises(dataclasses.FrozenInstanceError):
        assessment.qualified = False


def test_it_serialises_without_nan_or_infinity():
    payload = evaluate_panel(request(), PRIMERS, Metrics()).as_dict()
    text = json.dumps(payload)

    assert json.loads(text) == payload
    assert "Infinity" not in text and "NaN" not in text


def test_the_model_versions_travel_with_the_result():
    """A saved assessment names the code that produced it."""
    assessment = evaluate_panel(request(), PRIMERS, Metrics())

    assert any(name == "neoswga" for name, _version in assessment.model_versions)


# ---------------------------------------------------------------------------
# Non-finite quantities fail
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("bad", [float("nan"), float("inf"), float("-inf")])
def test_a_non_finite_required_quantity_fails_the_run(bad):
    """NaN compares False against every threshold.

    A panel carrying one passes no limit and fails no limit, so a constrained
    design reports a clean result computed from a number arithmetic cannot use.
    """
    with pytest.raises(ModelEvaluationError, match="fg_coverage"):
        evaluate_panel(request(), PRIMERS, Metrics(fg_coverage=bad))


def test_a_missing_required_quantity_fails_rather_than_reading_zero():
    with pytest.raises(ModelEvaluationError, match="fg_coverage"):
        evaluate_panel(request(), PRIMERS, Metrics(fg_coverage=None))


def test_a_non_finite_optional_quantity_also_fails():
    """Optional means it may be absent, not that it may be nonsense."""
    with pytest.raises(ModelEvaluationError, match="max_gap"):
        evaluate_panel(request(), PRIMERS, Metrics(max_gap=float("inf")))


# ---------------------------------------------------------------------------
# Unavailable is not zero
# ---------------------------------------------------------------------------


def test_an_unavailable_quantity_says_so_rather_than_reading_zero():
    assessment = evaluate_panel(request(), PRIMERS, Metrics(bg_coverage=None))
    measurement = assessment.metrics["bg_coverage"]

    assert measurement.value is None
    assert measurement.unavailable
    assert assessment.value("bg_coverage") is None


def test_a_measured_zero_is_a_value_and_not_an_absence():
    assessment = evaluate_panel(request(), PRIMERS, Metrics(bg_coverage=0.0))
    measurement = assessment.metrics["bg_coverage"]

    assert measurement.value == 0.0
    assert not measurement.unavailable


def test_a_measurement_cannot_be_both_present_and_absent():
    with pytest.raises(ValueError):
        Measurement("x", 1.0, "unit", unavailable="also missing")
    with pytest.raises(ValueError):
        Measurement("x", None, "unit")


# ---------------------------------------------------------------------------
# A zero background is represented, not divided by
# ---------------------------------------------------------------------------


def test_no_host_sites_gives_an_undefined_ratio_rather_than_infinity():
    """Infinity is not a number JSON carries, and not a claim anyone has evidence for.

    A panel with no exact host match is a real and good measurement. What it is
    not is an enrichment estimate: mismatched binding is not counted, so the
    honest record is a zero denominator with the site count beside it.
    """
    assessment = evaluate_panel(
        request(), PRIMERS, Metrics(total_bg_sites=0, selectivity_ratio=float("inf"))
    )

    assert assessment.zero_background is True
    ratio = assessment.metrics["selectivity_ratio"]
    assert ratio.value is None
    assert ratio.unavailable == "zero denominator"
    assert assessment.metrics["total_bg_sites"].value == 0
    assert "mismatched binding is not counted" in " ".join(assessment.notes)
    json.dumps(assessment.as_dict())


def test_a_non_zero_background_keeps_its_ratio():
    assessment = evaluate_panel(request(), PRIMERS, Metrics())

    assert assessment.zero_background is False
    assert assessment.metrics["selectivity_ratio"].value == pytest.approx(7.7)


# ---------------------------------------------------------------------------
# Qualification
# ---------------------------------------------------------------------------


def test_no_configured_limits_means_no_configured_limit_can_fail():
    """Setting no limit must leave the delivered panel byte-identical.

    Narrowed 2026-09-22. The promise is about CONFIGURED limits, and still
    holds. It never was that an assessment cannot fail at all: a panel short of
    its requested size, or carrying a pair above `max_dimer_bp`, fails whether
    or not a limit is configured, because neither is configurable.
    """
    assessment = evaluate_panel(request(), PRIMERS, Metrics(), objective=None)

    assert assessment.violations == ()
    assert assessment.qualified is True


def test_a_violated_limit_is_named_and_disqualifies():
    assessment = evaluate_panel(
        request(), PRIMERS, Metrics(), objective=Objective(["selectivity below minimum"])
    )

    assert assessment.violations == ("selectivity below minimum",)
    assert assessment.qualified is False


def test_qualification_is_exactly_the_absence_of_violations():
    """One boolean, derived, so a record cannot say acceptable and list failures."""
    for violations in ((), ("a",), ("a", "b")):
        assessment = evaluate_panel(request(), PRIMERS, Metrics(), objective=Objective(violations))
        assert assessment.qualified == (not violations)


# ---------------------------------------------------------------------------
# What the numbers are, and are not
# ---------------------------------------------------------------------------


def test_the_record_does_not_call_coverage_a_recovery_probability():
    assessment = evaluate_panel(request(), PRIMERS, Metrics())
    notes = " ".join(assessment.notes).lower()

    assert "proxy" in notes
    assert "not a predicted sequencing breadth" in notes
    assert "recovery probability" not in notes


def test_the_evidence_status_of_the_reach_travels_with_the_number():
    """The reach is an assumption, and the record says so where it is read.

    A coverage figure is only interpretable beside the status of the reach it
    was computed at, and that status is `assumed`: a design-density convention
    from sets with wet-lab success, never measured in this repository.
    """
    assessment = evaluate_panel(request(), PRIMERS, Metrics())

    assert "assumed" in assessment.evidence["coverage_reach"]
    assert "never measured" in assessment.evidence["coverage_reach"]


def test_per_target_results_are_carried_individually():
    """Aggregate coverage hides a starved target.

    A panel covering one target 0.9 and another 0.1 beats a balanced 0.5/0.5
    panel on the mean, and nothing in selection balances across targets.
    """
    assessment = evaluate_panel(
        request(),
        PRIMERS,
        Metrics(per_target_coverage={"chrA": 0.9, "chrB": 0.1}),
    )

    assert assessment.per_target["chrA"].value == pytest.approx(0.9)
    assert assessment.per_target["chrB"].value == pytest.approx(0.1)


def test_exact_match_counts_say_that_they_are_exact_match_counts():
    """A hypothetical mismatch diagnostic must not read as measured specificity."""
    assessment = evaluate_panel(request(), PRIMERS, Metrics())

    assert assessment.metrics["total_bg_sites"].basis == "exact matches only"
    assert assessment.metrics["total_fg_sites"].basis == "exact matches only"


def test_the_same_metrics_under_the_same_request_give_the_same_record():
    """One record, so a metric and a verdict cannot be computed separately."""
    first = evaluate_panel(request(), PRIMERS, Metrics()).as_dict()
    second = evaluate_panel(request(), PRIMERS, Metrics()).as_dict()

    assert first == second
    assert not math.isnan(first["metrics"]["fg_coverage"]["value"])


def test_the_record_replaces_the_large_finite_sentinel():
    """`base_optimizer` reports 1e6 for a zero background, deliberately.

    That sentinel is better than infinity: it is finite and JSON carries it.
    Its own docstring says it means "no background binding was detected", not
    "measured this well" -- which is an acknowledgement that a reader cannot
    tell the two apart from the number alone. A downstream consumer comparing
    selectivity across designs sees 1e6 and a genuinely excellent 900 as two
    points on one scale.

    The assessment says the quantity is undefined and why, and keeps the site
    count that makes it interpretable. The sentinel is left in place because
    changing it moves every saved summary; this test records the difference so
    the choice between them is visible.
    """
    from neoswga.core.base_optimizer import MAX_SELECTIVITY, _selectivity_from_loads

    assert _selectivity_from_loads(120.0, 0.0) == MAX_SELECTIVITY

    assessment = evaluate_panel(
        request(), PRIMERS, Metrics(total_bg_sites=0, selectivity_ratio=MAX_SELECTIVITY)
    )
    ratio = assessment.metrics["selectivity_ratio"]

    assert ratio.value is None
    assert ratio.unavailable == "zero denominator"
    assert "no exact-match site" in ratio.basis
