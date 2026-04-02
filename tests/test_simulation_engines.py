"""
Tests for the simulation engines in NeoSWGA.

Covers the agent-based phi29 replication simulator, including
fork dynamics, strand displacement, and basic simulation runs.
"""

import random

import numpy as np
import pytest

from neoswga.core.replication_simulator import (
    DisplacedStrand,
    ForkState,
    Phi29Simulator,
    ReplicationFork,
    SimulationConfig,
    SimulationResult,
)
from neoswga.core.reaction_conditions import ReactionConditions, get_standard_conditions


# ---------------------------------------------------------------------------
# ForkState enum
# ---------------------------------------------------------------------------

class TestForkState:
    """Verify ForkState enum values."""

    def test_active_state(self):
        assert ForkState.ACTIVE.value == "active"

    def test_paused_state(self):
        assert ForkState.PAUSED.value == "paused"

    def test_terminated_state(self):
        assert ForkState.TERMINATED.value == "terminated"

    def test_collided_state(self):
        assert ForkState.COLLIDED.value == "collided"

    def test_displaced_state(self):
        assert ForkState.DISPLACED.value == "displaced"

    def test_all_states_present(self):
        expected = {"ACTIVE", "PAUSED", "TERMINATED", "COLLIDED", "DISPLACED"}
        assert set(ForkState.__members__.keys()) == expected


# ---------------------------------------------------------------------------
# ReplicationFork dataclass
# ---------------------------------------------------------------------------

class TestReplicationFork:
    """Verify ReplicationFork construction and default state."""

    def test_basic_creation(self):
        fork = ReplicationFork(
            fork_id=0,
            primer="ATCGATCG",
            start_position=100,
            current_position=100,
            direction="forward",
            strand="top",
        )
        assert fork.fork_id == 0
        assert fork.primer == "ATCGATCG"
        assert fork.start_position == 100
        assert fork.current_position == 100
        assert fork.direction == "forward"
        assert fork.strand == "top"

    def test_default_state_is_active(self):
        fork = ReplicationFork(
            fork_id=1,
            primer="GCTAGCTA",
            start_position=0,
            current_position=0,
            direction="reverse",
            strand="bottom",
        )
        assert fork.state == ForkState.ACTIVE

    def test_default_speed(self):
        fork = ReplicationFork(
            fork_id=2,
            primer="AAAA",
            start_position=0,
            current_position=0,
            direction="forward",
            strand="top",
        )
        assert fork.speed == 167.0

    def test_defaults(self):
        fork = ReplicationFork(
            fork_id=0,
            primer="ATCG",
            start_position=0,
            current_position=0,
            direction="forward",
            strand="top",
        )
        assert fork.birth_time == 0.0
        assert fork.termination_time is None
        assert fork.bases_synthesized == 0
        assert fork.template_id == 0
        assert fork.displaced_by is None


# ---------------------------------------------------------------------------
# DisplacedStrand dataclass
# ---------------------------------------------------------------------------

class TestDisplacedStrand:
    """Verify DisplacedStrand construction and properties."""

    def test_length_property(self):
        strand = DisplacedStrand(
            strand_id=0,
            start_pos=100,
            end_pos=600,
            strand="top",
            template_id=0,
            creation_time=1.0,
        )
        assert strand.length == 500

    def test_contains_position(self):
        strand = DisplacedStrand(
            strand_id=0,
            start_pos=100,
            end_pos=600,
            strand="top",
            template_id=0,
            creation_time=1.0,
        )
        assert strand.contains_position(100)
        assert strand.contains_position(300)
        assert strand.contains_position(599)
        assert not strand.contains_position(600)
        assert not strand.contains_position(99)

    def test_defaults(self):
        strand = DisplacedStrand(
            strand_id=0,
            start_pos=0,
            end_pos=100,
            strand="bottom",
            template_id=0,
            creation_time=0.0,
        )
        assert strand.is_available is True
        assert strand.copy_number == 1


# ---------------------------------------------------------------------------
# SimulationConfig defaults
# ---------------------------------------------------------------------------

class TestSimulationConfig:
    """Verify SimulationConfig default values."""

    def test_default_duration(self):
        cfg = SimulationConfig()
        assert cfg.duration == 28800.0

    def test_default_time_step(self):
        cfg = SimulationConfig()
        assert cfg.time_step == 1.0

    def test_default_polymerase_type(self):
        cfg = SimulationConfig()
        assert cfg.polymerase_type == "phi29"

    def test_default_extension_rate(self):
        cfg = SimulationConfig()
        assert cfg.extension_rate == 150.0

    def test_default_processivity_limit(self):
        cfg = SimulationConfig()
        assert cfg.processivity_limit == 70000

    def test_mechanistic_model_defaults(self):
        cfg = SimulationConfig()
        assert cfg.use_mechanistic_model is True
        assert cfg.probabilistic_processivity is True

    def test_strand_displacement_defaults(self):
        cfg = SimulationConfig()
        assert cfg.enable_strand_displacement is True
        assert cfg.displacement_probability == 0.95
        assert cfg.displaced_strand_binding_boost == 3.0

    def test_custom_values(self):
        cfg = SimulationConfig(duration=10, time_step=0.5, extension_rate=100)
        assert cfg.duration == 10
        assert cfg.time_step == 0.5
        assert cfg.extension_rate == 100


# ---------------------------------------------------------------------------
# Helper to build a small simulator
# ---------------------------------------------------------------------------

def _make_small_genome(length=1000):
    """Generate a simple repeating genome sequence."""
    random.seed(42)
    bases = "ATCG"
    return "".join(random.choice(bases) for _ in range(length))


def _make_simulator(
    genome_length=1000,
    primers=None,
    primer_positions=None,
    duration=10,
    time_step=1.0,
    use_mechanistic=False,
    enable_displacement=False,
):
    """Build a Phi29Simulator with short-run settings for testing."""
    genome_seq = _make_small_genome(genome_length)
    conditions = get_standard_conditions()

    if primers is None:
        primers = ["ATCGATCG"]
    if primer_positions is None:
        primer_positions = {
            p: {"forward": [100, 500], "reverse": [300, 700]}
            for p in primers
        }

    config = SimulationConfig(
        duration=duration,
        time_step=time_step,
        use_mechanistic_model=use_mechanistic,
        enable_strand_displacement=enable_displacement,
        probabilistic_processivity=False,
        use_temperature_effects=False,
        record_interval=1.0,
    )

    sim = Phi29Simulator(
        primers=primers,
        primer_positions=primer_positions,
        genome_length=genome_length,
        genome_sequence=genome_seq,
        conditions=conditions,
        config=config,
    )
    return sim


# ---------------------------------------------------------------------------
# Phi29Simulator initialization
# ---------------------------------------------------------------------------

class TestPhi29SimulatorInit:
    """Verify Phi29Simulator initialization."""

    def test_basic_init(self):
        sim = _make_simulator()
        assert sim.genome_length == 1000
        assert len(sim.primers) == 1
        assert sim.current_time == 0.0

    def test_coverage_array_initialized(self):
        sim = _make_simulator()
        assert sim.coverage.shape == (1000,)
        assert sim.coverage.dtype == bool
        assert sim.coverage.sum() == 0

    def test_fork_counters_start_at_zero(self):
        sim = _make_simulator()
        assert sim.forks_created == 0
        assert sim.forks_terminated == 0
        assert sim.fork_id_counter == 0

    def test_conditions_stored(self):
        sim = _make_simulator()
        assert sim.conditions is not None
        assert sim.conditions.polymerase == "phi29"

    def test_template_gc_calculated(self):
        sim = _make_simulator()
        assert 0.0 < sim.template_gc < 1.0


# ---------------------------------------------------------------------------
# Simulation run
# ---------------------------------------------------------------------------

class TestSimulationRun:
    """Test short simulation runs and result structure."""

    def test_short_run_returns_result(self):
        """A short simulation should return a SimulationResult."""
        random.seed(42)
        np.random.seed(42)
        sim = _make_simulator(duration=10, time_step=1.0)
        result = sim.run(verbose=False)
        assert isinstance(result, SimulationResult)

    def test_result_has_coverage_array(self):
        random.seed(42)
        np.random.seed(42)
        sim = _make_simulator(duration=10)
        result = sim.run(verbose=False)
        assert isinstance(result.coverage, np.ndarray)
        assert result.coverage.shape == (1000,)

    def test_result_has_coverage_fraction(self):
        random.seed(42)
        np.random.seed(42)
        sim = _make_simulator(duration=10)
        result = sim.run(verbose=False)
        assert 0.0 <= result.final_coverage_fraction <= 1.0

    def test_result_has_fork_counts(self):
        random.seed(42)
        np.random.seed(42)
        sim = _make_simulator(duration=10)
        result = sim.run(verbose=False)
        assert result.num_forks_created >= 0
        assert result.num_forks_terminated >= 0

    def test_result_has_amplification_fields(self):
        random.seed(42)
        np.random.seed(42)
        sim = _make_simulator(duration=10)
        result = sim.run(verbose=False)
        assert result.amplification_fold >= 1.0
        assert result.total_bases_synthesized >= 0

    def test_result_has_coverage_over_time(self):
        random.seed(42)
        np.random.seed(42)
        sim = _make_simulator(duration=10, time_step=1.0)
        result = sim.run(verbose=False)
        assert isinstance(result.coverage_over_time, list)

    def test_result_has_fork_history(self):
        random.seed(42)
        np.random.seed(42)
        sim = _make_simulator(duration=10)
        result = sim.run(verbose=False)
        assert isinstance(result.fork_history, list)

    def test_result_has_copy_number_array(self):
        random.seed(42)
        np.random.seed(42)
        sim = _make_simulator(duration=10)
        result = sim.run(verbose=False)
        assert result.copy_number_array is not None
        assert result.copy_number_array.shape == (1000,)


# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------

class TestEdgeCases:
    """Edge case scenarios for the simulator."""

    def test_no_primer_positions(self):
        """No binding sites should yield zero or near-zero coverage."""
        random.seed(42)
        np.random.seed(42)
        sim = _make_simulator(
            primers=["ATCGATCG"],
            primer_positions={"ATCGATCG": {"forward": [], "reverse": []}},
            duration=10,
        )
        result = sim.run(verbose=False)
        assert isinstance(result, SimulationResult)
        assert result.num_forks_created == 0
        assert result.final_coverage_fraction == 0.0

    def test_single_primer_single_position(self):
        """Single primer at one position should still produce a valid result."""
        random.seed(42)
        np.random.seed(42)
        sim = _make_simulator(
            primers=["ATCGATCG"],
            primer_positions={"ATCGATCG": {"forward": [500], "reverse": []}},
            duration=10,
        )
        result = sim.run(verbose=False)
        assert isinstance(result, SimulationResult)
        assert result.final_coverage_fraction >= 0.0

    def test_very_small_genome(self):
        """A very small genome (100 bp) should run without error."""
        random.seed(42)
        np.random.seed(42)
        sim = _make_simulator(
            genome_length=100,
            primers=["ATCG"],
            primer_positions={"ATCG": {"forward": [10], "reverse": [50]}},
            duration=5,
        )
        result = sim.run(verbose=False)
        assert isinstance(result, SimulationResult)
        assert result.coverage.shape == (100,)


# ---------------------------------------------------------------------------
# Determinism with fixed seed
# ---------------------------------------------------------------------------

class TestDeterminism:
    """Verify that simulation results are reproducible with a fixed seed."""

    def test_reproducible_coverage(self):
        """Two runs with the same seed should yield identical coverage."""
        random.seed(123)
        np.random.seed(123)
        sim1 = _make_simulator(duration=10)
        result1 = sim1.run(verbose=False)

        random.seed(123)
        np.random.seed(123)
        sim2 = _make_simulator(duration=10)
        result2 = sim2.run(verbose=False)

        np.testing.assert_array_equal(result1.coverage, result2.coverage)

    def test_reproducible_fork_count(self):
        """Fork counts should match across seeded runs."""
        random.seed(123)
        np.random.seed(123)
        sim1 = _make_simulator(duration=10)
        result1 = sim1.run(verbose=False)

        random.seed(123)
        np.random.seed(123)
        sim2 = _make_simulator(duration=10)
        result2 = sim2.run(verbose=False)

        assert result1.num_forks_created == result2.num_forks_created
        assert result1.num_forks_terminated == result2.num_forks_terminated
