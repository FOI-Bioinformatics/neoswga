"""Phase 1 (production-readiness v2): crash + garbage-data fixes.

- design-oligos command + optimal_oligo_generator placeholder module removed.
- simulate --visualize/--report wired to dict-compatible functions (were crashing).
- simulation_report dict report + swga_simulator.simulate_validation honest label.
"""

import numpy as np
import pytest

# ---------------------------------------------------------------------------
# design-oligos removal
# ---------------------------------------------------------------------------


def test_design_oligos_command_removed():
    from neoswga.cli_unified import create_parser

    parser = create_parser()
    sub = next(a for a in parser._actions if a.__class__.__name__ == "_SubParsersAction")
    assert "design-oligos" not in sub.choices


def test_optimal_oligo_generator_module_deleted():
    with pytest.raises(ModuleNotFoundError):
        __import__("neoswga.core.optimal_oligo_generator")


# ---------------------------------------------------------------------------
# Replication-simulation report / plots accept the simulate_primer_set dict
# ---------------------------------------------------------------------------


class _RepResult:
    """Stand-in for a replication_simulator.SimulationResult."""

    def __init__(self, cov, amp, forks, length=200):
        self.final_coverage_fraction = cov
        self.amplification_fold = amp
        self.num_forks_created = forks
        self.coverage = np.arange(length) < int(length * cov)


def _results_dict():
    reps = [_RepResult(0.8, 12.0, 30), _RepResult(0.9, 15.0, 35)]
    return {
        "mean_coverage": 0.85,
        "std_coverage": 0.05,
        "mean_amplification": 13.5,
        "std_amplification": 1.5,
        "mean_forks_created": 32,
        "mean_fork_travel": 2500,
        "mean_displaced_strands": 10,
        "mean_bases_synthesized": 1_000_000,
        "all_results": reps,
    }


def test_generate_replication_report_writes_real_html(tmp_path):
    from neoswga.core.simulation_report import generate_replication_report

    out = tmp_path / "sim.html"
    generate_replication_report(_results_dict(), ["ATCG", "GGGG"], str(out))
    html = out.read_text()
    assert "Replication Simulation Report" in html
    assert "85.0%" in html  # mean coverage, measured
    assert "N/A" not in html  # no fabricated placeholder
    assert "Replicate" in html


def test_plot_replication_summary_runs(tmp_path):
    pytest.importorskip("matplotlib")
    from neoswga.core.simulation_plots import plot_replication_summary

    out = tmp_path / "sim.png"
    plot_replication_summary(_results_dict(), str(out))
    assert out.exists() and out.stat().st_size > 0


def test_plot_replication_summary_handles_empty(tmp_path):
    pytest.importorskip("matplotlib")
    from neoswga.core.simulation_plots import plot_replication_summary

    out = tmp_path / "empty.png"
    plot_replication_summary({"all_results": []}, str(out))
    assert out.exists()


# ---------------------------------------------------------------------------
# simulate_validation honest labelling
# ---------------------------------------------------------------------------


def test_simulate_validation_labels_fallback(monkeypatch):
    from neoswga.core import swga_simulator as ss

    class _FakeResult:
        mode = "detailed"
        confidence = 0.9

    sim = ss.SwgaSimulator.__new__(ss.SwgaSimulator)  # bypass heavy __init__
    monkeypatch.setattr(
        sim, "simulate_detailed", lambda num_replicates=10: _FakeResult(), raising=False
    )
    out = sim.simulate_validation()
    assert "validation fallback" in out.mode
    assert out.confidence <= 0.7
