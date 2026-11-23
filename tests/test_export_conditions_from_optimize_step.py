"""`export` and `report` take the reaction from the step that chose the panel.

The manifest is append-only. A `filter` re-run after a design appends an entry
with its own conditions, and the unfiltered read returned that entry, so the
optimize step's panel was Tm-corrected under the filter step's buffer
(2026-11-23 sweep, tier 1 item 11).
"""

import json

from neoswga.core.run_manifest import write_manifest


def _results_with_two_entries(tmp_path):
    results = tmp_path / "results"
    results.mkdir()
    (results / "step4_improved_df.csv").write_text("primer\nATGCATGCATGC\n")
    (results / "params.json").write_text(json.dumps({"polymerase": "phi29"}))
    write_manifest(
        "optimize", str(results), effective_conditions={"polymerase": "phi29", "betaine_m": 1.0}
    )
    write_manifest(
        "filter", str(results), effective_conditions={"polymerase": "phi29", "betaine_m": 1.5}
    )
    return results


def test_export_conditions_ignore_a_later_filter_entry(tmp_path):
    from neoswga.core.export import PrimerExporter

    results = _results_with_two_entries(tmp_path)
    exporter = PrimerExporter.from_results_dir(str(results))
    assert exporter.additives.get("betaine_m") == 1.0


def test_report_conditions_ignore_a_later_filter_entry(tmp_path):
    from neoswga.core.report.metrics import _recorded_conditions

    results = _results_with_two_entries(tmp_path)
    assert _recorded_conditions(results)["betaine_m"] == 1.0
