"""Phase 1 hardening regression tests.

Each test exercises a code path that previously raised at runtime:
- dimer_network_analyzer: empty primer list built DimerNetworkMetrics with
  nonexistent field names -> TypeError.
- amplicon_network: warnings.warn used without `import warnings` -> NameError.
- swga_simulator.simulate_detailed: stale Phi29Simulator API (wrong import,
  kwargs, run signature, and result attribute). We validate the corrected
  Phi29Simulator API surface that the fix now depends on.
"""

import warnings

import numpy as np
import pytest


def test_empty_primer_set_returns_metrics_without_typeerror():
    """analyze_primer_set([]) must return valid metrics, not raise TypeError."""
    from neoswga.core.dimer_network_analyzer import (
        DimerNetworkAnalyzer,
        DimerNetworkMetrics,
    )

    analyzer = DimerNetworkAnalyzer()
    metrics, profiles, matrix = analyzer.analyze_primer_set([])

    assert isinstance(metrics, DimerNetworkMetrics)
    assert metrics.num_primers == 0
    assert metrics.total_interactions == 0
    assert metrics.num_hub_primers == 0
    assert metrics.severity_distribution == {}
    # An empty primer set is degenerate, not a passing result.
    assert metrics.passes is False
    assert profiles == {}
    assert matrix.size == 0


def test_amplicon_network_missing_positions_warns_not_nameerror(tmp_path):
    """A missing HDF5 prefix should emit a warning, not NameError."""
    from neoswga.core.amplicon_network import AmpliconNetwork

    network = AmpliconNetwork(primers=["ATCGATCG"], genome_length=1000)
    missing_prefix = str(tmp_path / "does_not_exist")

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        # Must not raise NameError; should record a warning and leave the
        # primer with empty position lists.
        network.load_positions_from_hdf5(missing_prefix)

    assert any("Could not load positions" in str(w.message) for w in caught)
    assert network.positions_forward["ATCGATCG"] == []
    assert network.positions_reverse["ATCGATCG"] == []


def test_phi29_simulator_api_matches_swga_simulator_fix(simple_genome, phi29_conditions):
    """The corrected swga_simulator.simulate_detailed relies on this exact
    Phi29Simulator API: positional genome_length, run(verbose=), .coverage."""
    from neoswga.core.replication_simulator import Phi29Simulator, SimulationConfig

    primers = ["ATCG"]
    positions = {"ATCG": {"forward": np.array([0, 4, 8]), "reverse": np.array([2, 6])}}
    config = SimulationConfig(duration=5.0, time_step=1.0)

    sim = Phi29Simulator(
        primers=primers,
        primer_positions=positions,
        genome_length=len(simple_genome),
        genome_sequence=simple_genome,
        conditions=phi29_conditions,
        config=config,
    )
    result = sim.run(verbose=False)

    # Attribute used by the swga_simulator fix is `.coverage` (not coverage_array).
    assert hasattr(result, "coverage")
    assert not hasattr(result, "coverage_array")
    assert isinstance(result.coverage, np.ndarray)
