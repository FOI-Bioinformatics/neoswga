"""Tests for simulation_analysis (previously dead module, now wired into
SwgaSimulator.analyze()). Uses synthetic positions so no genome I/O is needed.
"""

from types import SimpleNamespace

import numpy as np
import pytest

from neoswga.core.simulation_analysis import (
    ComprehensiveAnalysis,
    SimulationAnalyzer,
    print_analysis_report,
)


def _result(uniformity=0.8, enrichment=50.0, specificity=0.7):
    # SimulationAnalyzer reads only these three scalar fields off the result.
    return SimpleNamespace(
        target_uniformity=uniformity,
        enrichment=enrichment,
        specificity_score=specificity,
    )


def _positions(sites_by_primer):
    """Build {primer: {'+': array, '-': array}} from a dict of forward sites."""
    return {
        primer: {"+": np.array(fwd, dtype=int), "-": np.array(rev, dtype=int)}
        for primer, (fwd, rev) in sites_by_primer.items()
    }


def test_analyzer_produces_comprehensive_analysis():
    fg = _positions(
        {
            "AAAA": ([1000, 25000, 80000], [5000, 30000]),
            "CCCC": ([60000, 95000], [12000]),
        }
    )
    bg = _positions({"AAAA": ([200], []), "CCCC": ([], [])})

    analyzer = SimulationAnalyzer(
        result=_result(),
        fg_positions=fg,
        bg_positions=bg,
        fg_length=100_000,
        bg_length=1_000_000,
        bin_size=10_000,
    )
    analysis = analyzer.analyze()

    assert isinstance(analysis, ComprehensiveAnalysis)
    # Coverage is fraction of 10 bins (100kb / 10kb) that hold >=1 site.
    assert 0.0 < analysis.coverage.overall_coverage <= 1.0
    assert analysis.coverage.bins_total == 10
    # Specificity: 8 fg sites vs 1 bg site -> high enrichment density.
    assert analysis.specificity.target_sites == 8
    assert analysis.specificity.background_sites == 1
    # One contribution row per primer, sorted by contribution score desc.
    assert len(analysis.primer_contributions) == 2
    scores = [c.contribution_score for c in analysis.primer_contributions]
    assert scores == sorted(scores, reverse=True)
    assert 0.0 <= analysis.quality_score <= 1.0
    assert isinstance(analysis.recommendations, list)
    assert isinstance(analysis.issues, list)


def test_analyzer_flags_low_coverage_and_gaps():
    # All sites in the first bin -> huge downstream gap, low coverage.
    fg = _positions({"AAAA": ([100, 200, 300], [150])})
    bg = _positions({"AAAA": ([], [])})
    analyzer = SimulationAnalyzer(
        result=_result(uniformity=0.2, enrichment=3.0, specificity=0.2),
        fg_positions=fg,
        bg_positions=bg,
        fg_length=500_000,
        bg_length=1_000_000,
        bin_size=10_000,
    )
    analysis = analyzer.analyze()
    # 1 covered bin of 50 -> low coverage; a critical (>100kb) gap is impossible
    # here (single bin), but coverage and issues should reflect the sparsity.
    assert analysis.coverage.overall_coverage < 0.1
    assert any("coverage" in r.lower() for r in analysis.recommendations)
    assert any("CRITICAL" in issue for issue in analysis.issues)


def test_print_analysis_report_runs(capsys):
    fg = _positions({"AAAA": ([1000, 50000], [25000])})
    bg = _positions({"AAAA": ([], [])})
    analysis = SimulationAnalyzer(
        result=_result(),
        fg_positions=fg,
        bg_positions=bg,
        fg_length=100_000,
        bg_length=1_000_000,
    ).analyze()
    print_analysis_report(analysis)
    out = capsys.readouterr().out
    assert "COMPREHENSIVE SIMULATION ANALYSIS" in out


def test_swga_simulator_analyze_wires_the_analyzer():
    """SwgaSimulator.analyze() must build the analyzer from its own state."""
    from neoswga.core.swga_simulator import SwgaSimulator

    sim = SwgaSimulator.__new__(SwgaSimulator)  # bypass genome/HDF5 loading
    sim.fg_positions = _positions({"AAAA": ([1000, 40000], [20000])})
    sim.bg_positions = _positions({"AAAA": ([100], [])})
    sim.fg_length = 100_000
    sim.bg_length = 1_000_000
    sim.bin_size = 10_000

    # analyze(result=...) must NOT re-run a simulation when a result is given.
    analysis = sim.analyze(result=_result())
    assert isinstance(analysis, ComprehensiveAnalysis)
    assert analysis.specificity.target_sites == 3
