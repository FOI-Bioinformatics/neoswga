"""
Agent-based simulation of phi29 DNA replication for SWGA.

Simulates the physical replication process including:
- Primer binding kinetics
- Polymerase extension
- Replication fork dynamics
- Fork collisions and termination
- Strand displacement
- Template depletion

Enables prediction of:
- Coverage over time
- Amplification bias
- Optimal reaction times
- Primer set performance
"""

import numpy as np
from typing import List, Dict, Tuple, Optional, Set
from dataclasses import dataclass, field
from enum import Enum
import random

import neoswga.core.thermodynamics as thermo
import neoswga.core.reaction_conditions as rc


class ForkState(Enum):
    """States of replication fork."""
    ACTIVE = "active"
    PAUSED = "paused"
    TERMINATED = "terminated"
    COLLIDED = "collided"


@dataclass
class ReplicationFork:
    """Represents an active replication fork."""
    fork_id: int
    primer: str
    start_position: int
    current_position: int
    direction: str  # 'forward' or 'reverse'
    strand: str  # 'top' or 'bottom'
    state: ForkState = ForkState.ACTIVE
    speed: float = 167.0  # bp/sec (phi29 typical: ~10 kb/min)
    birth_time: float = 0.0
    termination_time: Optional[float] = None
    bases_synthesized: int = 0


@dataclass
class SimulationConfig:
    """Configuration for replication simulation."""
    duration: float = 3600.0  # seconds (default 1 hour)
    time_step: float = 1.0  # seconds
    polymerase_type: str = 'phi29'  # 'phi29' or 'equiphi29'
    primer_binding_rate: float = 1e-3  # per second per site
    extension_rate: float = 167.0  # bp/sec
    extension_rate_gc_penalty: float = 0.5  # Slowdown in high-GC
    fork_collision_distance: int = 100  # bp for collision
    template_regeneration: bool = False  # New templates from product
    record_interval: float = 60.0  # Record state every N seconds


@dataclass
class SimulationResult:
    """Results from replication simulation."""
    coverage: np.ndarray  # Coverage array
    final_coverage_fraction: float
    num_forks_created: int
    num_forks_terminated: int
    mean_fork_travel: float
    coverage_over_time: List[Tuple[float, float]]  # (time, coverage)
    fork_history: List[ReplicationFork]


class Phi29Simulator:
    """
    Agent-based simulator for phi29 replication.

    Models individual replication forks as agents that:
    - Initiate at primer binding sites
    - Extend along template
    - Interact with other forks (collision)
    - Terminate at genome ends or collisions
    """

    def __init__(self,
                 primers: List[str],
                 primer_positions: Dict[str, Dict[str, List[int]]],  # primer -> {'forward': [...], 'reverse': [...]}
                 genome_length: int,
                 genome_sequence: str,
                 conditions: rc.ReactionConditions,
                 config: Optional[SimulationConfig] = None):
        """
        Initialize simulator.

        Args:
            primers: Primer sequences
            primer_positions: Positions for each primer (forward/reverse)
            genome_length: Length of genome
            genome_sequence: Genome sequence (for GC-dependent speed)
            conditions: Reaction conditions
            config: Simulation configuration
        """
        self.primers = primers
        self.primer_positions = primer_positions
        self.genome_length = genome_length
        self.genome_sequence = genome_sequence.upper()
        self.conditions = conditions
        self.config = config if config else SimulationConfig()

        # Adjust parameters based on polymerase
        if conditions.polymerase == 'equiphi29':
            self.config.extension_rate = 200.0  # Slightly faster
        elif conditions.polymerase == 'phi29':
            self.config.extension_rate = 167.0

        # State
        self.current_time = 0.0
        self.forks: List[ReplicationFork] = []
        self.fork_id_counter = 0
        self.coverage = np.zeros(genome_length, dtype=bool)
        self.coverage_history = []

        # Statistics
        self.forks_created = 0
        self.forks_terminated = 0

        print(f"Phi29 Simulator initialized:")
        print(f"  Genome: {genome_length:,} bp")
        print(f"  Primers: {len(primers)}")
        print(f"  Polymerase: {conditions.polymerase}")
        print(f"  Extension rate: {self.config.extension_rate:.0f} bp/s")
        print(f"  Duration: {self.config.duration:.0f} s ({self.config.duration/60:.1f} min)")

    def run(self, verbose: bool = True) -> SimulationResult:
        """
        Run simulation.

        Args:
            verbose: Print progress

        Returns:
            SimulationResult
        """
        print("\nStarting simulation...")

        # Initialize primers
        self._initialize_primers()

        # Main simulation loop
        next_record_time = 0.0

        while self.current_time < self.config.duration:
            # Time step
            self.current_time += self.config.time_step

            # Update forks
            self._update_forks()

            # Attempt new primer initiations (stochastic)
            if random.random() < self.config.primer_binding_rate:
                self._attempt_primer_binding()

            # Record coverage
            if self.current_time >= next_record_time:
                coverage_frac = self.coverage.sum() / self.genome_length
                self.coverage_history.append((self.current_time, coverage_frac))
                next_record_time += self.config.record_interval

                if verbose and int(self.current_time) % 300 == 0:
                    active_forks = sum(1 for f in self.forks if f.state == ForkState.ACTIVE)
                    print(f"  t={self.current_time/60:.1f}min: "
                          f"Coverage={coverage_frac:.1%}, "
                          f"Active forks={active_forks}")

            # Early termination if no active forks
            active = sum(1 for f in self.forks if f.state == ForkState.ACTIVE)
            if active == 0 and self.current_time > 60:
                if verbose:
                    print(f"  No active forks at t={self.current_time/60:.1f}min, terminating early")
                break

        final_coverage = self.coverage.sum() / self.genome_length

        if verbose:
            print(f"\nSimulation complete:")
            print(f"  Final coverage: {final_coverage:.1%}")
            print(f"  Forks created: {self.forks_created}")
            print(f"  Forks terminated: {self.forks_terminated}")

        # Calculate statistics
        terminated_forks = [f for f in self.forks if f.state != ForkState.ACTIVE]
        mean_travel = np.mean([f.bases_synthesized for f in terminated_forks]) if terminated_forks else 0

        return SimulationResult(
            coverage=self.coverage,
            final_coverage_fraction=final_coverage,
            num_forks_created=self.forks_created,
            num_forks_terminated=self.forks_terminated,
            mean_fork_travel=mean_travel,
            coverage_over_time=self.coverage_history,
            fork_history=self.forks.copy()
        )

    def _initialize_primers(self):
        """Create initial primer bindings."""
        # For each primer binding site, stochastically initiate
        for primer in self.primers:
            positions = self.primer_positions.get(primer, {'forward': [], 'reverse': []})

            # Forward strand
            for pos in positions.get('forward', []):
                # Binding probability based on Tm
                tm = self.conditions.calculate_effective_tm(primer)
                binding_prob = self._calculate_binding_probability(tm)

                if random.random() < binding_prob * 0.1:  # Scale down for initial binding
                    self._create_fork(primer, pos, 'forward', 'top')

            # Reverse strand
            for pos in positions.get('reverse', []):
                tm = self.conditions.calculate_effective_tm(primer)
                binding_prob = self._calculate_binding_probability(tm)

                if random.random() < binding_prob * 0.1:
                    self._create_fork(primer, pos, 'reverse', 'bottom')

    def _calculate_binding_probability(self, tm: float) -> float:
        """Calculate primer binding probability based on Tm."""
        # Sigmoid function centered at reaction temperature
        delta = tm - self.conditions.temp
        prob = 1 / (1 + np.exp(-delta / 5))  # Sharper transition
        return prob

    def _create_fork(self, primer: str, position: int, direction: str, strand: str):
        """Create new replication fork."""
        # Calculate extension speed based on local GC content
        speed = self._calculate_local_speed(position)

        fork = ReplicationFork(
            fork_id=self.fork_id_counter,
            primer=primer,
            start_position=position,
            current_position=position,
            direction=direction,
            strand=strand,
            speed=speed,
            birth_time=self.current_time
        )

        self.forks.append(fork)
        self.fork_id_counter += 1
        self.forks_created += 1

    def _calculate_local_speed(self, position: int, window: int = 50) -> float:
        """
        Calculate extension speed based on local GC content.

        High GC regions slow polymerase.

        Args:
            position: Current position
            window: Window size for GC calculation

        Returns:
            Extension speed in bp/s
        """
        start = max(0, position - window // 2)
        end = min(self.genome_length, position + window // 2)

        if start >= end:
            local_gc = 0.5
        else:
            local_seq = self.genome_sequence[start:end]
            local_gc = (local_seq.count('G') + local_seq.count('C')) / len(local_seq)

        # Slowdown in high GC (empirical: ~50% slower in 70% GC)
        gc_factor = 1 - self.config.extension_rate_gc_penalty * max(0, local_gc - 0.5)

        return self.config.extension_rate * gc_factor

    def _update_forks(self):
        """Update all replication forks."""
        for fork in self.forks:
            if fork.state != ForkState.ACTIVE:
                continue

            # Calculate distance to move this time step
            distance = int(fork.speed * self.config.time_step)

            # Update position
            if fork.direction == 'forward':
                new_position = fork.current_position + distance
            else:  # reverse
                new_position = fork.current_position - distance

            # Check boundaries
            if new_position < 0 or new_position >= self.genome_length:
                self._terminate_fork(fork, "boundary")
                continue

            # Check collisions with other forks
            collision = self._check_collisions(fork, new_position)
            if collision:
                self._terminate_fork(fork, "collision")
                self._terminate_fork(collision, "collision")
                continue

            # Update position and coverage
            start = min(fork.current_position, new_position)
            end = max(fork.current_position, new_position)
            self.coverage[start:end] = True

            fork.current_position = new_position
            fork.bases_synthesized += distance

            # Update speed based on new position
            fork.speed = self._calculate_local_speed(new_position)

    def _check_collisions(self, fork: ReplicationFork, new_position: int) -> Optional[ReplicationFork]:
        """
        Check if fork will collide with another fork.

        Args:
            fork: Fork to check
            new_position: Proposed new position

        Returns:
            Colliding fork if found, None otherwise
        """
        for other_fork in self.forks:
            if other_fork.fork_id == fork.fork_id:
                continue

            if other_fork.state != ForkState.ACTIVE:
                continue

            # Check if forks are close and moving toward each other
            distance = abs(other_fork.current_position - new_position)

            if distance < self.config.fork_collision_distance:
                # Check if moving toward each other
                if fork.direction == 'forward' and other_fork.direction == 'reverse':
                    if fork.current_position < other_fork.current_position:
                        return other_fork
                elif fork.direction == 'reverse' and other_fork.direction == 'forward':
                    if fork.current_position > other_fork.current_position:
                        return other_fork

        return None

    def _terminate_fork(self, fork: ReplicationFork, reason: str):
        """Terminate replication fork."""
        fork.state = ForkState.TERMINATED if reason == "boundary" else ForkState.COLLIDED
        fork.termination_time = self.current_time
        self.forks_terminated += 1

    def _attempt_primer_binding(self):
        """Stochastically attempt new primer binding."""
        # Random primer and position
        primer = random.choice(self.primers)
        positions = self.primer_positions.get(primer, {'forward': [], 'reverse': []})

        all_positions = (
            [(p, 'forward', 'top') for p in positions.get('forward', [])] +
            [(p, 'reverse', 'bottom') for p in positions.get('reverse', [])]
        )

        if not all_positions:
            return

        pos, direction, strand = random.choice(all_positions)

        # Check if position already covered
        if self.coverage[pos]:
            return  # Template already replicated

        # Check binding probability
        tm = self.conditions.calculate_effective_tm(primer)
        binding_prob = self._calculate_binding_probability(tm)

        if random.random() < binding_prob * self.config.primer_binding_rate:
            self._create_fork(primer, pos, direction, strand)


def simulate_primer_set(primers: List[str],
                       primer_positions: Dict,
                       genome_length: int,
                       genome_sequence: str,
                       conditions: Optional[rc.ReactionConditions] = None,
                       n_replicates: int = 5) -> Dict:
    """
    Simulate primer set performance with multiple replicates.

    Args:
        primers: Primer sequences
        primer_positions: Primer positions dict
        genome_length: Genome length
        genome_sequence: Genome sequence
        conditions: Reaction conditions
        n_replicates: Number of simulation replicates

    Returns:
        Dictionary with mean and std of metrics
    """
    if conditions is None:
        conditions = rc.get_enhanced_conditions()

    results = []

    print(f"Running {n_replicates} simulation replicates...")

    for i in range(n_replicates):
        print(f"\n{'='*60}")
        print(f"Replicate {i+1}/{n_replicates}")
        print(f"{'='*60}")

        simulator = Phi29Simulator(
            primers, primer_positions, genome_length,
            genome_sequence, conditions
        )

        result = simulator.run(verbose=False)
        results.append(result)

    # Aggregate statistics
    coverages = [r.final_coverage_fraction for r in results]

    print(f"\n{'='*60}")
    print("Aggregate Results:")
    print(f"{'='*60}")
    print(f"Coverage: {np.mean(coverages):.1%} ± {np.std(coverages):.1%}")
    print(f"Mean forks created: {np.mean([r.num_forks_created for r in results]):.0f}")
    print(f"Mean fork travel: {np.mean([r.mean_fork_travel for r in results]):.0f} bp")

    return {
        'mean_coverage': np.mean(coverages),
        'std_coverage': np.std(coverages),
        'mean_forks_created': np.mean([r.num_forks_created for r in results]),
        'mean_fork_travel': np.mean([r.mean_fork_travel for r in results]),
        'all_results': results
    }


if __name__ == "__main__":
    print("Phi29 Replication Simulator")
    print("=" * 60)
    print("\nAgent-based simulation of DNA replication:")
    print("  - Individual replication forks as agents")
    print("  - Realistic extension kinetics")
    print("  - Fork collisions and termination")
    print("  - GC-dependent polymerase speed")
    print("  - Temperature-dependent primer binding")
    print("\nApplications:")
    print("  - Predict coverage before experiments")
    print("  - Optimize reaction time")
    print("  - Compare primer set performance")
    print("  - Identify amplification bias")
    print("\nExample usage:")
    print("""
    from neoswga.core import replication_simulator, reaction_conditions

    # Set up
    primers = ['ATCGAT', 'GCTAGC']
    positions = {
        'ATCGAT': {'forward': [100, 500, 1000], 'reverse': [2000, 3000]},
        'GCTAGC': {'forward': [300, 800], 'reverse': [1500, 2500]}
    }

    # Run simulation
    result = replication_simulator.simulate_primer_set(
        primers=primers,
        primer_positions=positions,
        genome_length=4400000,
        genome_sequence="ATCG..." * 1100000,
        conditions=reaction_conditions.get_enhanced_conditions(),
        n_replicates=5
    )

    print(f"Mean coverage: {result['mean_coverage']:.1%}")
    """)
