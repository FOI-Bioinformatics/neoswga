"""Tests for the simulation plotting module.

`simulation_plots` sat at 25%. It renders the figure a user actually looks at
after a simulation, and plotting code fails in a characteristic way: it works on
the happy path and raises on the degenerate one -- an empty gap list, a primer
with no binding sites, a single replicate, a genome shorter than the bin size.
Those are exactly the cases a poor primer set produces, so the failure lands on
the run the user most needs to see.

`generate_plots` takes a `SwgaSimulator`, which only became constructible once
its position loader was fixed; that is why this module had no execution coverage.
"""

import os

import pytest

plt = pytest.importorskip("matplotlib.pyplot")
h5py = pytest.importorskip("h5py")

import numpy as np

from neoswga.core.swga_simulator import SwgaSimulator
from neoswga.core.thermodynamics import reverse_complement

PRIMERS = ["TTGACCATGA", "GCATTACGGT", "AACCGGTTAC"]


@pytest.fixture(autouse=True)
def close_figures():
    yield
    plt.close("all")


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
# The whole figure
# ----------------------------------------------------------------------


def test_generate_plots_writes_a_real_png(simulator, result, tmp_path):
    from neoswga.core.simulation_plots import generate_plots

    out = tmp_path / "plots.png"
    generate_plots(result, simulator, output_file=str(out))

    assert out.is_file()
    assert out.stat().st_size > 5_000, "PNG is too small to contain a figure"
    assert out.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n", "not a PNG"


def test_generate_plots_without_analysis(simulator, result, tmp_path):
    """The analysis argument is optional and changes the grid layout."""
    from neoswga.core.simulation_plots import generate_plots

    out = tmp_path / "no_analysis.png"
    generate_plots(result, simulator, analysis=None, output_file=str(out))
    assert out.is_file()


def test_generate_plots_with_analysis(simulator, result, tmp_path):
    from neoswga.core.simulation_analysis import SimulationAnalyzer
    from neoswga.core.simulation_plots import generate_plots

    analyzer = SimulationAnalyzer(
        result=result,
        fg_positions=simulator.fg_positions,
        bg_positions=simulator.bg_positions,
        fg_length=simulator.fg_length,
        bg_length=simulator.bg_length,
        bin_size=1_000,
    )
    analysis = analyzer.analyze()

    out = tmp_path / "with_analysis.png"
    generate_plots(result, simulator, analysis=analysis, output_file=str(out))

    assert out.is_file()
    assert out.stat().st_size > 5_000


# ----------------------------------------------------------------------
# Degenerate inputs -- what a poor primer set produces
# ----------------------------------------------------------------------


def test_coverage_heatmap_survives_a_primer_with_no_sites(simulator):
    """A set containing a dud primer must still plot."""
    from neoswga.core.simulation_plots import plot_coverage_heatmap

    positions = dict(simulator.fg_positions)
    positions["CGCGCGCGCG"] = {
        "+": np.array([]),
        "-": np.array([]),
        "forward": np.array([]),
        "reverse": np.array([]),
    }

    _fig, ax = plt.subplots()
    plot_coverage_heatmap(ax, positions, simulator.fg_length)


def test_coverage_heatmap_survives_no_positions_at_all(simulator):
    """The zero-coverage case, which is exactly when a user wants the plot."""
    from neoswga.core.simulation_plots import plot_coverage_heatmap

    empty = {
        primer: {
            "+": np.array([]),
            "-": np.array([]),
            "forward": np.array([]),
            "reverse": np.array([]),
        }
        for primer in PRIMERS
    }

    _fig, ax = plt.subplots()
    plot_coverage_heatmap(ax, empty, simulator.fg_length)


def test_primer_distribution_survives_empty_input(simulator):
    from neoswga.core.simulation_plots import plot_primer_distribution

    _fig, ax = plt.subplots()
    plot_primer_distribution(ax, {}, simulator.fg_length)


def test_gap_panel_survives_an_empty_gap_list():
    """No gaps is a good result; it must not break the figure that reports it."""
    from neoswga.core.simulation_plots import plot_gaps

    _fig, ax = plt.subplots()
    plot_gaps(ax, [])


def test_gap_panel_renders_labelled_gaps():
    """The shape `simulation_analysis` produces, which is what the figure uses."""
    from neoswga.core.simulation_plots import plot_gaps

    gaps = [
        {"start": 1_000, "end": 4_000, "length": 3_000, "severity": "low"},
        {"start": 12_000, "end": 160_000, "length": 148_000, "severity": "critical"},
    ]
    _fig, ax = plt.subplots()
    plot_gaps(ax, gaps)


def test_gap_panel_accepts_unlabelled_gaps():
    """`SwgaSimulator` emits gaps without a severity key.

    Handing those to plot_gaps -- the obvious thing to try, since both are gap
    lists in the same project -- raised a bare KeyError from inside the plotting
    code. They are now bucketed by length on simulation_analysis's own
    thresholds.
    """
    from neoswga.core.simulation_plots import plot_gaps

    gaps = [
        {"start": 1_000, "end": 4_000, "length": 3_000, "mean_coverage": 0.0},
        {"start": 20_000, "end": 180_000, "length": 160_000, "mean_coverage": 0.0},
    ]
    _fig, ax = plt.subplots()
    plot_gaps(ax, gaps)


def test_severity_bucketing_matches_the_analyser_thresholds():
    from neoswga.core.simulation_plots import _severity_from_size

    assert _severity_from_size(150_000) == "critical"
    assert _severity_from_size(60_000) == "high"
    assert _severity_from_size(30_000) == "medium"
    assert _severity_from_size(5_000) == "low"


def test_primer_contribution_panel_survives_empty_input():
    from neoswga.core.simulation_plots import plot_primer_contributions

    _fig, ax = plt.subplots()
    plot_primer_contributions(ax, [])


def test_metrics_panel_renders_the_result(result):
    from neoswga.core.simulation_plots import plot_metrics_summary

    _fig, ax = plt.subplots()
    plot_metrics_summary(ax, result)


def test_enrichment_panel_renders(simulator, result):
    from neoswga.core.simulation_plots import plot_enrichment

    _fig, ax = plt.subplots()
    plot_enrichment(ax, result, simulator)


# ----------------------------------------------------------------------
# The replication-summary figure, which reads a different result shape
# ----------------------------------------------------------------------


def test_replication_summary_consumes_simulate_primer_set_output(tmp_path):
    """It reads the dict `replication_simulator.simulate_primer_set` returns.

    That shape (mean_coverage, mean_amplification, all_results) is a contract
    between two modules and nothing else checks it.
    """
    from neoswga.core.simulation_plots import plot_replication_summary

    class _Replicate:
        final_coverage_fraction = 0.42
        amplification_factor = 12.0
        coverage = np.concatenate([np.ones(5_000), np.zeros(5_000)])

    results = {
        "mean_coverage": 0.42,
        "mean_amplification": 12.0,
        "all_results": [_Replicate(), _Replicate()],
    }

    out = tmp_path / "replication.png"
    plot_replication_summary(results, output_file=str(out))

    assert out.is_file()
    assert out.stat().st_size > 5_000


def test_replication_summary_survives_no_replicates(tmp_path):
    """A run that produced nothing still needs a figure saying so."""
    from neoswga.core.simulation_plots import plot_replication_summary

    out = tmp_path / "empty.png"
    plot_replication_summary({"all_results": []}, output_file=str(out))
    assert out.is_file()
