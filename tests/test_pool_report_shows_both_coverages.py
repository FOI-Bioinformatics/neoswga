"""A report must not let one coverage number stand for two different things.

Task 7 of the condition-aware pool design plan.

The report printed a single coverage figure and named the metric beside it.
Two quantities are in play and they answer different questions:

geometric coverage
    The union of extension windows around exact binding sites. How much of the
    target a panel could reach.

occupancy-weighted coverage
    The same windows, each weighted by how much of the time its site is actually
    bound at the reaction temperature. How much of the target a panel plausibly
    reaches under this chemistry.

Showing one invites the reader to treat a modelling choice as a result. Showing
both, side by side, makes the gap between them visible -- and that gap is the
whole contribution of the additive and temperature model.

It also has to be said plainly that a higher modelled density is not
demonstrated enrichment. Every number in this report is computed; none is a
laboratory measurement, and the evidence matrix records which of them have no
compatible observation at all.
"""

import json
import pathlib

import pytest

from neoswga.core.pool_plan_report import write_pool_plan

A, C = "AAAAAAAAAAAA", "CCCCCCCCCCCC"


def _plan(**overrides):
    plan = {
        "rows": [
            {
                "requested_size": 2,
                "size": 2,
                "primers": [A, C],
                "eligible": True,
                "status": "evaluated",
                "coverage": 0.62,
                "raw_coverage": 0.88,
                "effective_coverage": 0.62,
                "selectivity_density": 40.0,
                "background_sites": 12,
            }
        ],
        "recommendations": [{"target_coverage": 0.6, "row_index": 0, "size": 2}],
        "primer_length": 12,
        "coverage_metric": "effective",
        "extension_reach": 3000,
        "background_assessed": True,
        "min_selectivity_density": 10.0,
        "max_background_sites": None,
        "interpretation": "Smallest qualifying panels found among evaluated sizes.",
    }
    plan.update(overrides)
    return plan


def test_both_coverage_figures_appear_for_every_evaluated_panel(tmp_path):
    page = write_pool_plan(_plan(), tmp_path / "r").read_text()

    assert "88.0%" in page, "the geometric figure is missing"
    assert "62.0%" in page, "the occupancy-weighted figure is missing"


def test_the_two_columns_are_labelled_distinctly(tmp_path):
    page = write_pool_plan(_plan(), tmp_path / "r").read_text()

    assert "Geometric" in page
    assert "Occupancy-weighted" in page


def test_the_report_states_that_a_modelled_gain_is_not_measured_enrichment(tmp_path):
    page = write_pool_plan(_plan(), tmp_path / "r").read_text()

    assert "not demonstrated enrichment" in page or "not a demonstrated" in page


def test_the_resolved_chemistry_is_recorded_when_supplied(tmp_path):
    plan = _plan(
        design_parameters={"polymerase": "phi29", "reaction_temp": 30.0, "primer_conc": 2e-6},
        condition_fingerprint="tm-2026-09-14:deadbeef",
    )
    page = write_pool_plan(plan, tmp_path / "r").read_text()

    assert "tm-2026-09-14:deadbeef" in page, "the reaction identity is not recorded"


def test_the_per_oligo_and_total_concentration_are_both_stated(tmp_path):
    """They are different experiments, and the model assumes the first."""
    plan = _plan(design_parameters={"primer_conc": 5e-7})
    page = write_pool_plan(plan, tmp_path / "r").read_text()

    assert "per oligo" in page.lower()
    assert "total" in page.lower()


def test_search_accounting_is_reported_when_present(tmp_path):
    plan = _plan(
        design_counts={
            "counted": 5190,
            "assessed": 5190,
            "hard_qc_passed": 5190,
            "shortlisted": 37,
            "indexed": 5190,
            "examined": 37,
        },
        stop_reason="budget_exhausted",
    )
    page = write_pool_plan(plan, tmp_path / "r").read_text()

    assert "5,190" in page or "5190" in page
    assert "budget_exhausted" in page


def test_an_older_plan_without_the_new_fields_still_renders(tmp_path):
    """Historical reports must keep working; none of this is required."""
    page = write_pool_plan(_plan(), tmp_path / "r").read_text()

    assert "Evaluated panels" in page


def test_a_saved_plan_regenerates_without_rerunning_the_design(tmp_path):
    """The JSON the report writes is enough to rebuild the report."""
    first = tmp_path / "first"
    write_pool_plan(_plan(), first)
    saved = json.loads((first / "pool_plan.json").read_text())

    second = tmp_path / "second"
    page = write_pool_plan(saved, second).read_text()

    assert "Evaluated panels" in page
    assert (second / "pool_sizes.csv").is_file()
