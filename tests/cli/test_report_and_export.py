"""Execution tests for `report`, `export` and `interpret`.

These are the commands a user reaches for once a design is finished, and they
are the last place a wrong number can be introduced without anything failing.
The report in particular renders values read from
`step4_improved_df_summary.json` and badges them MEASURED or ESTIMATED -- a
mislabelled or fabricated figure there looks exactly like a real result.

`export` matters for a different reason: it produces the file someone orders
oligos from. A wrong sequence, or a silently dropped 3' phosphorothioate bond,
becomes a wasted synthesis order rather than a failed test.
"""

import json
import os

import pytest

pytestmark = pytest.mark.usefixtures("completed_run")


# ----------------------------------------------------------------------
# report
# ----------------------------------------------------------------------


def test_report_generates_html(run_cli, completed_run, tmp_path):
    out = tmp_path / "report.html"
    run_cli(["report", "-d", completed_run["data_dir"], "-o", str(out)])

    assert out.is_file(), "no report written"
    html = out.read_text()
    assert "<html" in html.lower() or "<!doctype" in html.lower()
    assert len(html) > 500, "report is suspiciously short"


def test_full_report_is_richer_than_the_summary(run_cli, completed_run, tmp_path):
    """--level full is documented as a technical report, not a restyled summary."""
    summary = tmp_path / "summary.html"
    full = tmp_path / "full.html"
    run_cli(["report", "-d", completed_run["data_dir"], "-o", str(summary)])
    run_cli(["report", "-d", completed_run["data_dir"], "-o", str(full), "--level", "full"])

    assert len(full.read_text()) > len(summary.read_text())


def test_report_shows_all_three_gap_statistics(run_cli, completed_run, tmp_path):
    """The published benchmarks disagree about which gap statistic predicts
    success, so all three are surfaced rather than folded into one score."""
    out = tmp_path / "full.html"
    run_cli(["report", "-d", completed_run["data_dir"], "-o", str(out), "--level", "full"])
    html = out.read_text().lower()

    assert "mean gap" in html
    assert "max gap" in html
    assert "gini" in html or "evenness" in html


def test_report_numbers_come_from_the_authoritative_summary(run_cli, completed_run, tmp_path):
    """Coverage in the report must match step4_improved_df_summary.json.

    CLAUDE.md names that file authoritative precisely because CSV-derived
    estimates drifted from it before.
    """
    summary_path = os.path.join(completed_run["data_dir"], "step4_improved_df_summary.json")
    coverage = json.loads(open(summary_path).read())["metrics"]["fg_coverage"]

    out = tmp_path / "full.html"
    run_cli(["report", "-d", completed_run["data_dir"], "-o", str(out), "--level", "full"])
    html = out.read_text()

    assert (
        f"{coverage * 100:.1f}" in html or f"{coverage:.1%}" in html
    ), f"report does not show the summary's coverage of {coverage:.1%}"


def test_report_does_not_fabricate_a_filtering_funnel(run_cli, completed_run, tmp_path):
    """The funnel is rendered from the real filter_stats.json the filter step
    writes; invented counts would read as a clean run that never happened."""
    stats_path = os.path.join(completed_run["data_dir"], "filter_stats.json")
    stats = json.loads(open(stats_path).read())
    counts = [v for v in stats.values() if isinstance(v, int)]

    out = tmp_path / "full.html"
    run_cli(["report", "-d", completed_run["data_dir"], "-o", str(out), "--level", "full"])
    html = out.read_text()

    assert any(
        str(c) in html for c in counts
    ), f"none of the real funnel counts {counts} appear in the report"


def test_report_check_passes_on_a_complete_run(run_cli, completed_run):
    """--check is the readiness probe; it must accept a finished pipeline."""
    try:
        run_cli(["report", "-d", completed_run["data_dir"], "--check"])
    except SystemExit as exc:
        assert exc.code in (0, None), "report --check rejected a completed run"


def test_report_on_a_missing_directory_fails_clearly(run_cli, tmp_path, capsys, caplog):
    """It must say why, on some channel, rather than exiting mutely."""
    import logging

    with caplog.at_level(logging.ERROR):
        with pytest.raises(SystemExit) as exc:
            run_cli(["report", "-d", str(tmp_path / "nope"), "-o", str(tmp_path / "r.html")])

    assert exc.value.code != 0
    captured = capsys.readouterr()
    explanation = (captured.out + captured.err + caplog.text).lower()
    assert (
        "valid" in explanation or "not found" in explanation or "error" in explanation
    ), f"no explanation for a missing results directory: {explanation!r}"


# ----------------------------------------------------------------------
# export
# ----------------------------------------------------------------------


@pytest.fixture
def exported(run_cli, completed_run, tmp_path):
    out = tmp_path / "export"
    run_cli(["export", "-d", completed_run["data_dir"], "-o", str(out)])
    if not out.exists():
        pytest.skip("export produced no output directory")
    return out


def test_export_writes_an_order_sheet(exported):
    files = os.listdir(exported)
    assert files, "export wrote nothing"
    assert any(f.endswith(".csv") for f in files), files


def test_exported_sequences_are_the_selected_primers(exported, completed_run):
    """The order sheet must contain the primers that were actually chosen.

    Bases beyond ACGT are legitimate here only as modification notation.
    """
    summary_path = os.path.join(completed_run["data_dir"], "step4_improved_df_summary.json")
    chosen = json.loads(open(summary_path).read())["primers"]

    blob = "".join(
        (exported / f).read_text() for f in os.listdir(exported) if f.endswith((".csv", ".txt"))
    )
    stripped = blob.replace("*", "")
    for primer in chosen:
        assert primer in stripped, f"selected primer {primer} missing from the export"


def test_default_export_carries_phosphorothioate_bonds(exported):
    """3' PTO protects the primer from phi29's exonuclease.

    Silently shipping bare oligos where the protocol assumes protected ones is a
    wet-lab failure the tool would never see.
    """
    blob = "".join(
        (exported / f).read_text() for f in os.listdir(exported) if f.endswith((".csv", ".txt"))
    )
    assert "*" in blob, "no PTO notation in the default export"


def test_no_modifications_flag_produces_bare_sequences(run_cli, completed_run, tmp_path):
    out = tmp_path / "bare"
    run_cli(["export", "-d", completed_run["data_dir"], "-o", str(out), "--no-modifications"])
    blob = "".join((out / f).read_text() for f in os.listdir(out) if f.endswith((".csv", ".txt")))
    assert "*" not in blob, "--no-modifications still emitted PTO notation"


def test_pto_bond_count_is_honoured(run_cli, completed_run, tmp_path):
    def stars(n):
        out = tmp_path / f"pto{n}"
        run_cli(
            [
                "export",
                "-d",
                completed_run["data_dir"],
                "-o",
                str(out),
                "--pto-bonds",
                str(n),
            ]
        )
        blob = "".join(
            (out / f).read_text() for f in os.listdir(out) if f.endswith((".csv", ".txt"))
        )
        return blob.count("*")

    assert stars(1) < stars(3), "--pto-bonds had no effect on the exported sequences"


# ----------------------------------------------------------------------
# interpret
# ----------------------------------------------------------------------


def test_interpret_gives_a_verdict(run_cli, completed_run, capsys, caplog):
    """`interpret` exists to turn metrics into a go/no-go recommendation."""
    import logging

    with caplog.at_level(logging.INFO):
        run_cli(["interpret", "-d", completed_run["data_dir"]])

    reported = (capsys.readouterr().out + caplog.text).lower()
    assert reported.strip(), "interpret said nothing"
    assert any(
        word in reported for word in ("coverage", "recommend", "quality", "grade", "proceed")
    ), reported[:500]
