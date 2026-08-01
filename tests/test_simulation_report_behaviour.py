"""Tests for the simulation HTML reports.

Two reports live here and they consume different result shapes:
`generate_html_report` takes a `SwgaSimulator` `SimulationResult`, while
`generate_replication_report` takes the dict `simulate_primer_set` returns. Only
the second had any coverage, and the first was untestable until the simulator's
position loader was fixed.

A report is the last thing between a computed number and a decision to order
oligos, so what matters is that the figures it shows are the figures that were
computed. Both docstrings claim to render only real simulation data; these tests
check that against the source values rather than just asserting HTML came out.
"""

import re

import pytest

h5py = pytest.importorskip("h5py")

import numpy as np

from neoswga.core.swga_simulator import SwgaSimulator
from neoswga.core.thermodynamics import reverse_complement

PRIMERS = ["TTGACCATGA", "GCATTACGGT", "AACCGGTTAC"]


def _occurrences(sequence, pattern):
    out, i = [], sequence.find(pattern)
    while i != -1:
        out.append(i)
        i = sequence.find(pattern, i + 1)
    return out


@pytest.fixture(scope="module")
def sequence():
    import random

    rng = random.Random(31337)
    seq = list("".join(rng.choice("ACGT") for _ in range(20_000)))
    for offset, primer in enumerate(PRIMERS):
        for pos in (1_000, 7_000, 13_000):
            start = pos + offset * 100
            seq[start : start + len(primer)] = list(primer)
    return "".join(seq)


@pytest.fixture
def simulator(tmp_path, sequence):
    fasta = tmp_path / "target.fasta"
    fasta.write_text(
        ">target\n" + "\n".join(sequence[i : i + 70] for i in range(0, len(sequence), 70)) + "\n"
    )
    h5 = tmp_path / "positions.h5"
    with h5py.File(h5, "w") as f:
        for primer in PRIMERS:
            for key in {primer, reverse_complement(primer)}:
                pos = _occurrences(sequence, key)
                if pos:
                    f.create_dataset(key, data=np.array(pos, dtype=np.int32))

    return SwgaSimulator(
        primers=PRIMERS,
        fg_genome=str(fasta),
        bg_genome=str(fasta),
        fg_positions_h5=str(h5),
        bg_positions_h5=str(h5),
        bin_size=1_000,
    )


@pytest.fixture
def result(simulator):
    return simulator.simulate_fast()


# ----------------------------------------------------------------------
# generate_html_report -- the SwgaSimulator path
# ----------------------------------------------------------------------


def test_html_report_is_written_and_self_contained(result, tmp_path):
    """It must not depend on a stylesheet or script the user does not have."""
    from neoswga.core.simulation_report import generate_html_report

    out = tmp_path / "report.html"
    generate_html_report(result, str(out))

    html = out.read_text()
    assert "<html" in html.lower()
    assert len(html) > 1_000

    external = re.findall(r'(?:src|href)\s*=\s*["\'](https?://[^"\']+)', html)
    assert not external, f"report pulls in external resources: {external}"


def test_report_shows_the_computed_coverage(result, tmp_path):
    """The figure in the report must be the figure that was computed."""
    from neoswga.core.simulation_report import generate_html_report

    out = tmp_path / "report.html"
    generate_html_report(result, str(out))
    html = out.read_text()

    assert f"{result.target_coverage * 100:.1f}" in html


def test_report_shows_the_recommendation(result, tmp_path):
    from neoswga.core.simulation_report import generate_html_report

    out = tmp_path / "report.html"
    generate_html_report(result, str(out))

    assert result.recommendation in out.read_text()


def test_report_names_the_simulation_mode(result, tmp_path):
    """`fast`, `detailed` and `validation` differ in what they can support.

    A report that does not say which one produced it invites a fast heuristic
    being read as a kinetic validation.
    """
    from neoswga.core.simulation_report import generate_html_report

    out = tmp_path / "report.html"
    generate_html_report(result, str(out))

    assert result.mode in out.read_text().lower()


def test_report_works_without_an_analysis(result, tmp_path):
    from neoswga.core.simulation_report import generate_html_report

    out = tmp_path / "no_analysis.html"
    generate_html_report(result, str(out), analysis=None)
    assert out.is_file()


def test_report_works_with_an_analysis(simulator, result, tmp_path):
    from neoswga.core.simulation_analysis import SimulationAnalyzer
    from neoswga.core.simulation_report import generate_html_report

    analysis = SimulationAnalyzer(
        result=result,
        fg_positions=simulator.fg_positions,
        bg_positions=simulator.bg_positions,
        fg_length=simulator.fg_length,
        bg_length=simulator.bg_length,
        bin_size=1_000,
    ).analyze()

    out = tmp_path / "with_analysis.html"
    generate_html_report(result, str(out), analysis=analysis)

    html = out.read_text()
    assert out.is_file()
    assert len(html) > 1_000


def test_report_handles_a_result_with_no_gaps(result, tmp_path):
    """No gaps is the good outcome; the gap section must cope with an empty list."""
    import dataclasses

    from neoswga.core.simulation_report import generate_html_report

    clean = dataclasses.replace(result, target_gaps=[])
    out = tmp_path / "no_gaps.html"
    generate_html_report(clean, str(out))
    assert out.is_file()


def test_report_handles_a_failed_result(result, tmp_path):
    """Zero coverage and no enrichment is exactly when a report is consulted."""
    import dataclasses

    from neoswga.core.simulation_report import generate_html_report

    failed = dataclasses.replace(
        result,
        target_coverage=0.0,
        target_uniformity=0.0,
        enrichment=0.0,
        composite_score=0.0,
        recommendation="POOR",
        target_gaps=[],
    )

    out = tmp_path / "failed.html"
    generate_html_report(failed, str(out))

    html = out.read_text()
    assert "POOR" in html
    assert "nan" not in html.lower() and "inf" not in html.lower().replace("info", "")


@pytest.mark.parametrize("enrichment", [1.0, 50.0, 500.0, 5_000.0])
def test_enrichment_commentary_covers_every_band(result, tmp_path, enrichment):
    """The report branches on enrichment (<10, <100, >1000); each branch must render."""
    import dataclasses

    from neoswga.core.simulation_report import generate_html_report

    variant = dataclasses.replace(result, enrichment=enrichment)
    out = tmp_path / f"enrich_{enrichment:.0f}.html"
    generate_html_report(variant, str(out))
    assert out.is_file()


# ----------------------------------------------------------------------
# generate_replication_report -- the simulate_primer_set path
# ----------------------------------------------------------------------


def _replication_results(n_replicates=3):
    class _Replicate:
        final_coverage_fraction = 0.42
        amplification_fold = 12.0
        coverage = np.concatenate([np.ones(5_000), np.zeros(5_000)])

    return {
        "mean_coverage": 0.42,
        "std_coverage": 0.03,
        "mean_amplification": 12.0,
        "std_amplification": 1.5,
        "mean_forks_created": 7.0,
        "mean_fork_travel": 9_824.0,
        "all_results": [_Replicate() for _ in range(n_replicates)],
    }


def test_replication_report_shows_the_primers_and_metrics(tmp_path):
    from neoswga.core.simulation_report import generate_replication_report

    out = tmp_path / "replication.html"
    generate_replication_report(_replication_results(), PRIMERS, str(out))

    html = out.read_text()
    for primer in PRIMERS:
        assert primer in html, f"{primer} missing from the report"
    assert "42" in html, "the mean coverage does not appear"


def test_replication_report_reports_the_replicate_count(tmp_path):
    """A single-replicate run and a ten-replicate run carry different weight."""
    from neoswga.core.simulation_report import generate_replication_report

    out = tmp_path / "reps.html"
    generate_replication_report(_replication_results(n_replicates=7), PRIMERS, str(out))
    assert "7" in out.read_text()


def test_replication_report_handles_no_replicates(tmp_path):
    from neoswga.core.simulation_report import generate_replication_report

    out = tmp_path / "empty.html"
    generate_replication_report({"all_results": []}, PRIMERS, str(out))

    html = out.read_text()
    assert out.is_file()
    assert "nan" not in html.lower()


def test_replication_report_is_self_contained(tmp_path):
    from neoswga.core.simulation_report import generate_replication_report

    out = tmp_path / "replication.html"
    generate_replication_report(_replication_results(), PRIMERS, str(out))

    external = re.findall(r'(?:src|href)\s*=\s*["\'](https?://[^"\']+)', out.read_text())
    assert not external, f"report pulls in external resources: {external}"
