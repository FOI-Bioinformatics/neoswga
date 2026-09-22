"""A design must say how much of the available candidate pool it examined.

`search_frontiers` stops at the first frontier that satisfies the configured
limits. On a real design that is about 2,000 candidates of the 365,073 the
inventory holds -- half a per cent -- and until 2026-09-22 nothing recorded it.
Both numbers were known at the moment the pool was opened.

That silence is not neutral. A reader comparing two coverage figures has no way
to tell a search over everything from a search over a shortlist, and the
project's own documentation described the frontier as reaching "all hard-QC
survivors", which four refills at doubling cannot do from a 2,000 shortlist.

It is a deliberate default and the measurement supports it: widening the
frontier raises coverage 0.3-1.5% relative and costs 9-26% of selectivity
density, because `max_primer` cuts on `bg_count / fg_count` ascending so the
candidates a refill reaches bind the host more. Reporting it does not change
it; it stops the default being invisible.

The rule that matters here is the one this project keeps relearning: a figure
that is absent must not render as a favourable value. A run over a named
candidate list has no universe behind it, and a directory written before this
field existed has no record -- and "0% of the pool searched" is a different and
alarming claim from "not recorded".
"""

import json

import pytest

from neoswga.core.candidate_source import CandidateFrontier, ListCandidateSource, describe_reach
from neoswga.core.report.funnel_section import render_candidate_reach as _render_candidate_reach

# ---------------------------------------------------------------------------
# The measurement
# ---------------------------------------------------------------------------


class Source:
    """A source whose frontier is smaller than its universe."""

    kind = "inventory"

    def __init__(self, universe, frontier):
        self._universe, self._frontier = universe, frontier

    def initial(self, limit=None):
        return [f"P{i}" for i in range(self._frontier)]

    def describe(self):
        return {"kind": self.kind, "universe": self._universe, "examined": self._frontier}


def test_a_partial_search_reports_both_numbers():
    reach = describe_reach(CandidateFrontier(Source(365073, 2000)))

    assert reach["examined"] == 2000
    assert reach["universe"] == 365073
    assert reach["complete"] is False
    assert reach["fraction"] == pytest.approx(0.005478, abs=1e-6)


def test_a_search_over_everything_says_so():
    """A named candidate list is the ordinary case: the user chose the pool."""
    reach = describe_reach(CandidateFrontier(ListCandidateSource(["AAA", "CCC", "GGG"])))

    assert reach["complete"] is True
    assert reach["examined"] == reach["universe"] == 3


def test_a_plain_list_gives_no_answer_rather_than_zero():
    """Absence, not a reach of nothing. The two must not render alike."""
    assert describe_reach(["AAA", "CCC"]) is None


def test_a_source_that_cannot_describe_itself_does_not_fail_the_run():
    """Reporting how much was searched must never lose the design."""

    class Broken(Source):
        def describe(self):
            raise RuntimeError("no")

    assert describe_reach(CandidateFrontier(Broken(10, 2))) is None


# ---------------------------------------------------------------------------
# Where a reader meets it
# ---------------------------------------------------------------------------


def test_the_report_states_a_partial_search():
    rendered = _render_candidate_reach(
        {"universe": 365073, "examined": 2000, "complete": False, "fraction": 0.005478}
    )

    assert "2,000" in rendered and "365,073" in rendered
    assert "0.5%" in rendered
    assert "deliberate default" in rendered, "a reader must not read this as a broken run"


def test_the_report_renders_nothing_when_the_figure_is_absent():
    """The half that matters. A missing measurement must not become 0%."""
    for absent in (None, {}, {"universe": 0, "examined": 0}):
        assert _render_candidate_reach(absent) == "", absent


def test_the_report_does_not_apologise_for_a_complete_search():
    rendered = _render_candidate_reach(
        {"universe": 500, "examined": 500, "complete": True, "fraction": 1.0}
    )

    assert "all 500" in rendered
    assert "deliberate default" not in rendered


# ---------------------------------------------------------------------------
# The saved result carries it
# ---------------------------------------------------------------------------


def test_the_summary_records_the_reach(tmp_path):
    from neoswga.core.base_optimizer import (
        OptimizationResult,
        OptimizationStatus,
        PrimerSetMetrics,
    )
    from neoswga.core.step4_output import save_results

    result = OptimizationResult(
        primers=("ACGTACGTACGT",),
        score=1.0,
        status=OptimizationStatus.SUCCESS,
        metrics=PrimerSetMetrics.empty(),
        iterations=1,
        optimizer_name="test",
    )
    out = tmp_path / "step4_improved_df.csv"
    reach = {"kind": "inventory", "universe": 999, "examined": 10, "complete": False}
    save_results(result, str(out), candidate_reach=reach)

    summary = json.loads((tmp_path / "step4_improved_df_summary.json").read_text())

    assert summary["candidate_reach"] == reach


def test_the_summary_omits_the_reach_when_there_is_none(tmp_path):
    from neoswga.core.base_optimizer import (
        OptimizationResult,
        OptimizationStatus,
        PrimerSetMetrics,
    )
    from neoswga.core.step4_output import save_results

    result = OptimizationResult(
        primers=("ACGTACGTACGT",),
        score=1.0,
        status=OptimizationStatus.SUCCESS,
        metrics=PrimerSetMetrics.empty(),
        iterations=1,
        optimizer_name="test",
    )
    out = tmp_path / "step4_improved_df.csv"
    save_results(result, str(out), candidate_reach=None)

    summary = json.loads((tmp_path / "step4_improved_df_summary.json").read_text())

    assert "candidate_reach" not in summary, (
        "an absent reach was written as a value; a reader cannot then tell "
        "'not recorded' from a measurement"
    )
