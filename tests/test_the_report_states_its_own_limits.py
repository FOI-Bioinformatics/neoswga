"""A report that shows a coverage percentage must say what it is.

Nothing in this project has ever measured sequencing recovery. Coverage is a
geometric proxy at a declared reach, its reach is an assumed constant, and some
quantities cannot be measured at all on a given run -- a selectivity ratio has
no value when the host carries no site to divide by. All of that lives in the
saved panel assessment, and until 2026-09-27 neither report said any of it: a
reader met a number and nothing qualifying it.

`render_claim_limits` is shared by both reports, in the manner of
`render_validation_banner`, and returns "" for a directory written before the
assessment existed. An empty string is the honest output there: an older record
cannot state caveats it never stored, and assembling them from elsewhere would
attribute them to a measurement that did not make them.
"""

import pytest

from neoswga.core.report.utils import render_claim_limits


def _assessment(**overrides):
    base = {
        "metrics": {
            "fg_coverage": {
                "value": 0.63,
                "units": "fraction of target bases",
                "basis": "reach 3000 bp; denominator is total target length",
            },
        },
        "notes": [
            "Coverage is a geometric proxy at the declared reach. It is not a "
            "predicted sequencing breadth and has not been calibrated against one."
        ],
        "blocking_violations": [],
        "advisory_violations": [],
        "evidence": {"coverage_reach": "assumed (never measured)"},
    }
    base.update(overrides)
    return base


def test_the_reach_a_figure_was_computed_at_is_stated():
    html = render_claim_limits(_assessment())
    assert "reach 3000 bp" in html
    assert "not comparable" in html, (
        "two coverage figures at different reaches are not comparable and the "
        "report must say so, because the number alone looks comparable"
    )


def test_the_proxy_is_named_as_a_proxy():
    html = render_claim_limits(_assessment())
    assert "geometric proxy" in html
    assert "not a predicted sequencing breadth" in html


def test_an_unmeasured_quantity_is_named_with_its_reason():
    """Absent must not be readable as zero."""
    html = render_claim_limits(
        _assessment(
            metrics={
                "selectivity_ratio": {
                    "value": None,
                    "units": "dimensionless",
                    "basis": "undefined",
                    "unavailable": "zero denominator",
                }
            }
        )
    )
    assert "selectivity_ratio was not measured" in html
    assert "zero denominator" in html
    assert "absent above rather than zero" in html


def test_an_advisory_violation_is_not_presented_as_a_defect():
    html = render_claim_limits(
        _assessment(advisory_violations=["panel size 4 is below the requested 6"])
    )
    assert "below the requested 6" in html
    assert "documented outcome, not a defect" in html


def test_a_blocking_violation_is_presented_as_unmet():
    html = render_claim_limits(_assessment(blocking_violations=["selectivity below minimum"]))
    assert "Unmet requirement: selectivity below minimum" in html


def test_the_evidence_status_of_a_constant_is_shown():
    html = render_claim_limits(_assessment())
    assert "coverage_reach" in html and "assumed" in html


def test_nothing_is_rendered_without_an_assessment():
    assert render_claim_limits(None) == ""
    assert render_claim_limits({}) == ""
    assert render_claim_limits("not a record") == ""


def test_content_is_escaped():
    html = render_claim_limits(_assessment(advisory_violations=["<script>alert(1)</script>"]))
    assert "<script>" not in html
    assert "&lt;script&gt;" in html


@pytest.mark.parametrize("report", ["technical", "executive"])
def test_both_reports_splice_the_section(report):
    """A helper nothing calls is the defect class this repository calls out."""
    import pathlib

    module = {
        "technical": "neoswga/core/report/technical_report.py",
        "executive": "neoswga/core/report/executive_summary.py",
    }[report]
    source = pathlib.Path(module).read_text()
    assert "render_claim_limits(" in source, f"{report} report does not render the limits"
    assert "claim_limits_css" in source, f"{report} report does not carry the styling"
