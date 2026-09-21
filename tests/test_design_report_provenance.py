"""Every exported artifact names the same panel, and a failed run exports none.

Task 8 of the 2026-09-21 valid-design plan: "prove JSON, HTML, CSV and
recommended FASTA identify the same panel and qualification", and "test failed
final validation and stale output directories so neither can leave an
apparently current recommendation".

An output directory is the only thing a later command sees. It holds a CSV of
primers, a summary JSON of metrics, and whatever the last run left behind, and
none of those carries a timestamp anyone compares. So a directory where the
most recent run FAILED looks exactly like one where it succeeded, minus
whatever the failure prevented being written -- and the previous run's outputs
are still sitting there.

Task 1 added `design_failure.json` for exactly this, and nothing read it. A
record written and never consulted is the Known Issue 8 class in its
artifact form: the evidence exists, the check does not.
"""

import csv
import json

import pytest

from neoswga.core.design_result import RunState, recommendation_allowed

PANEL = ["AAAACCCCGGGG", "TTTTGGGGCCCC", "ACGTACGTACGT"]


def build_results_dir(tmp_path, primers=PANEL, failure=None):
    """A results directory as the pipeline leaves one."""
    directory = tmp_path / "results"
    directory.mkdir(exist_ok=True)

    with open(directory / "step4_improved_df.csv", "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["primer", "set_index", "fg_count", "bg_count"])
        for primer in primers:
            writer.writerow([primer, 0, 80, 10])

    (directory / "step4_improved_df_summary.json").write_text(
        json.dumps(
            {
                "primers": list(primers),
                "metrics": {"fg_coverage": 0.62, "total_bg_sites": 31, "extension_reach": 3000},
            }
        )
    )
    (directory / "params.json").write_text(
        json.dumps(
            {
                "fg_genomes": ["target.fna"],
                "fg_prefixes": ["target"],
                "fg_seq_lengths": [100000],
                "bg_prefixes": [],
                "polymerase": "phi29",
                "reaction_temp": 30.0,
                "min_k": 12,
                "max_k": 12,
            }
        )
    )
    if failure is not None:
        (directory / "design_failure.json").write_text(json.dumps(failure))
    return directory


def failure_record(run_state=RunState.FAILED):
    return {
        "run_state": run_state,
        "termination": "error",
        "stage": "optimize",
        "error_type": "ReferenceDataError",
        "message": "position index built from a different reference",
        "recommendation_written": False,
        "qualified": False,
    }


# ---------------------------------------------------------------------------
# The artifacts agree with each other
# ---------------------------------------------------------------------------


def test_the_csv_and_the_summary_name_the_same_panel(tmp_path):
    directory = build_results_dir(tmp_path)

    with open(directory / "step4_improved_df.csv") as handle:
        from_csv = [row["primer"] for row in csv.DictReader(handle)]
    from_summary = json.loads((directory / "step4_improved_df_summary.json").read_text())["primers"]

    assert from_csv == from_summary


def test_the_exported_fasta_holds_exactly_the_delivered_panel(tmp_path):
    """A recommendation that omits or adds an oligo is a different pool."""
    from neoswga.core.export import PrimerExporter

    directory = build_results_dir(tmp_path)
    exporter = PrimerExporter.from_results_dir(str(directory))
    path = tmp_path / "order.fasta"
    exporter.export_fasta(str(path), prefix="test", include_metadata=False)

    sequences = [
        line.strip()
        for line in path.read_text().splitlines()
        if line.strip() and not line.startswith(">")
    ]

    assert sequences == PANEL


def test_the_report_and_the_summary_name_the_same_panel(tmp_path):
    from neoswga.core.report.metrics import collect_pipeline_metrics

    directory = build_results_dir(tmp_path)
    metrics = collect_pipeline_metrics(str(directory))

    assert [primer.sequence for primer in metrics.primers] == PANEL


# ---------------------------------------------------------------------------
# A failed run leaves nothing exportable
# ---------------------------------------------------------------------------


def test_a_recorded_failure_makes_the_directory_not_recommendable(tmp_path):
    """The gate, at the level of the record.

    `recommendation_allowed` already refuses this. What follows is whether
    anything asks it.
    """
    record = failure_record()

    assert recommendation_allowed(record["run_state"], qualified=True) is False


@pytest.mark.parametrize("state", [RunState.FAILED, RunState.INTERRUPTED])
def test_export_refuses_a_directory_whose_last_run_did_not_finish(tmp_path, state):
    """The check that was missing.

    `design_failure.json` was written by Task 1 and read by nothing, so a
    directory holding last week's `step4_improved_df.csv` beside this week's
    failure exported the old panel as though it were current. The CSV is real;
    it is just not the result of the run the user last asked for.
    """
    from neoswga.core.export import export_is_blocked

    directory = build_results_dir(tmp_path, failure=failure_record(state))
    blocked = export_is_blocked(str(directory))

    assert blocked, "a failed run must not leave an exportable recommendation"
    assert "optimize" in blocked or "did not finish" in blocked


def test_export_allows_a_directory_with_no_failure_record(tmp_path):
    from neoswga.core.export import export_is_blocked

    assert export_is_blocked(str(build_results_dir(tmp_path))) is None


def test_a_successful_rerun_clears_the_stale_failure_record(tmp_path):
    """Otherwise the first failure would block every later success.

    A record that is never cleared is as wrong as one that is never read: the
    user fixes the problem, the run succeeds, and the export still refuses.
    """
    from neoswga.cli._failure import clear_failure_artifact
    from neoswga.core.export import export_is_blocked

    directory = build_results_dir(tmp_path, failure=failure_record())
    assert export_is_blocked(str(directory))

    clear_failure_artifact(str(directory))

    assert not (directory / "design_failure.json").exists()
    assert export_is_blocked(str(directory)) is None


def test_clearing_is_safe_when_there_is_nothing_to_clear(tmp_path):
    from neoswga.cli._failure import clear_failure_artifact

    directory = build_results_dir(tmp_path)
    clear_failure_artifact(str(directory))
    clear_failure_artifact(str(directory))


def test_a_finished_run_that_found_nothing_is_also_refused(tmp_path):
    """Finished is not qualified. A run can complete having found no panel."""
    from neoswga.core.export import export_is_blocked

    record = failure_record(RunState.FINISHED)
    record["termination"] = "candidates_exhausted"
    record["qualified"] = False
    directory = build_results_dir(tmp_path, failure=record)

    assert export_is_blocked(str(directory))


def test_an_unreadable_failure_record_blocks_rather_than_passes(tmp_path):
    """A record that cannot be parsed is unknown, and unknown is not success.

    Defaulting to "export anyway" would make a corrupted artifact the most
    permissive state available, which is the shape of every silent-zero this
    project has recorded.
    """
    from neoswga.core.export import export_is_blocked

    directory = build_results_dir(tmp_path)
    (directory / "design_failure.json").write_text("{not json")

    assert export_is_blocked(str(directory))
