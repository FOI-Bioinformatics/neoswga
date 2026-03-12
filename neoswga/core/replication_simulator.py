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

import math
import numpy as np
from typing import List, Dict, Tuple, Optional, Set, TYPE_CHECKING
from dataclasses import dataclass, field
from enum import Enum
import random
import warnings

from neoswga.core import thermodynamics as thermo
from neoswga.core import reaction_conditions as rc
from neoswga.core.mechanistic_params import get_polymerase_params, MECHANISTIC_MODEL_PARAMS

if TYPE_CHECKING:
    from neoswga.core.mechanistic_model import MechanisticModel, MechanisticEffects


class ForkState(Enum):
    """States of replication fork."""
    ACTIVE = "active"
    PAUSED = "paused"
    TERMINATED = "terminated"
    COLLIDED = "collided"
    DISPLACED = "displaced"  # Fork displaced by another fork


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
    template_id: int = 0  # Which template this fork is on (0 = original)
    displaced_by: Optional[int] = None  # Fork ID that displaced this one


@dataclass
class DisplacedStrand:
    """
    Represents a displaced single-stranded DNA segment.

    When phi29 polymerase encounters a previously synthesized strand,
    it displaces it, creating a new single-stranded template. This
    displaced strand can then serve as a template for new synthesis,
    leading to exponential (hyperbranching) amplification.

    Attributes:
        strand_id: Unique identifier
        start_pos: Start position on original genome
        end_pos: End position on original genome
        strand: 'top' or 'bottom' (which strand of original genome)
        template_id: Parent template this was displaced from
        creation_time: When this strand was displaced
        is_available: Whether this strand can be used as template
        copy_number: How many copies of this strand exist
    """
    strand_id: int
    start_pos: int
    end_pos: int
    strand: str  # 'top' or 'bottom'
    template_id: int  # Which template this came from
    creation_time: float
    is_available: bool = True
    copy_number: int = 1

    @property
    def length(self) -> int:
        """Length of this displaced strand."""
        return self.end_pos - self.start_pos

    def contains_position(self, pos: int) -> bool:
        """Check if position is within this strand."""
        return self.start_pos <= pos < self.end_pos


@dataclass
class SimulationConfig:
    """Configuration for replication simulation.

    Defaults calibrated to literature values for phi29:
        - Processivity: 70,000 bp (Blanco 1989)
        - Extension rate: 100-170 nt/s (Esteban 1993, Blanco 1989)
        - Duration: 8 hours (typical SWGA protocol)

    Mechanistic model integration:
        When use_mechanistic_model=True, extension rates, processivity,
        and binding kinetics are modified by MechanisticEffects factors:
        - speed_factor: Modifies extension_rate
        - processivity_factor: Modifies processivity_limit
        - kon_factor: Modifies primer_binding_rate
        - accessibility_factor: Modifies binding probability for GC-rich regions

    Probabilistic processivity:
        When probabilistic_processivity=True, uses geometric distribution
        for polymerase dissociation instead of hard limit. This gives
        more realistic behavior where some forks terminate early and
        others extend beyond the nominal processivity.

    Strand displacement / Hyperbranching:
        When enable_strand_displacement=True, the simulator models phi29's
        ability to displace existing strands when it encounters them. This
        creates new single-stranded templates that can be primed, leading
        to exponential (hyperbranching) amplification.

        Key parameters:
        - displacement_probability: Probability of displacing vs terminating
        - displaced_strand_binding_boost: ssDNA binds primers more efficiently
        - max_amplification_fold: Resource limit on total amplification
    """
    duration: float = 28800.0  # seconds (default 8 hours, typical SWGA)
    time_step: float = 1.0  # seconds
    polymerase_type: str = 'phi29'  # 'phi29' or 'equiphi29'
    primer_binding_rate: float = 1e-3  # per second per site
    extension_rate: float = 150.0  # bp/sec (median of 100-170 range)
    extension_rate_gc_penalty: float = 0.5  # Slowdown in high-GC
    fork_collision_distance: int = 100  # bp for collision
    template_regeneration: bool = True  # Legacy flag (use enable_strand_displacement)
    record_interval: float = 60.0  # Record state every N seconds
    processivity_limit: int = 70000  # bp before polymerase dissociates (phi29)

    # Mechanistic model integration
    use_mechanistic_model: bool = True  # Apply MechanisticEffects to simulation
    probabilistic_processivity: bool = True  # Use probabilistic termination

    # Temperature effects
    use_temperature_effects: bool = True  # Apply temp deviation penalty

    # Strand displacement / Hyperbranching
    enable_strand_displacement: bool = True  # Model strand displacement

    # Phi29 displaces downstream strands with high efficiency owing to its
    # intrinsic helicase activity (Blanco et al. 1989 JBC 264:8935-8940).
    # Paez et al. (2004) PNAS 101:6261-6266 observed near-complete strand
    # displacement in MDA reactions.  A value of 0.95 reflects this high
    # but not absolute displacement efficiency.
    displacement_probability: float = 0.95

    # Single-stranded DNA binds short primers more readily than dsDNA because
    # no strand separation is needed.  Kool (1998) Annu Rev Biophys Biomol
    # Struct 27:1-44 measured a ~2-5x binding advantage for ssDNA templates.
    # We use 3.0 as a moderate estimate within that range.
    displaced_strand_binding_boost: float = 3.0

    min_displaced_strand_length: int = 500  # Don't track very short displacements
    max_displaced_strands: int = 10000  # Memory limit on tracked strands

    # Dean et al. (2002) PNAS 99:5261-5266 reported 10^4-10^6 fold
    # amplification in MDA reactions.  We use the upper bound as the
    # resource-depletion cap.
    max_amplification_fold: float = 1e6

    max_concurrent_forks: int = 1000  # Polymerase availability limit


@dataclass
class SimulationResult:
    """Results from replication simulation.

    Includes both coverage metrics (fraction of genome covered) and
    amplification metrics (total DNA synthesized / hyperbranching).
    """
    coverage: np.ndarray  # Coverage array (bool per position)
    final_coverage_fraction: float
    num_forks_created: int
    num_forks_terminated: int
    mean_fork_travel: float
    coverage_over_time: List[Tuple[float, float]]  # (time, coverage)
    fork_history: List[ReplicationFork]

    # Amplification metrics (for hyperbranching)
    amplification_fold: float = 1.0  # Total DNA / original genome
    total_bases_synthesized: int = 0  # Total bases synthesized
    num_displaced_strands: int = 0  # Number of displaced strand templates
    amplification_over_time: List[Tuple[float, float]] = field(default_factory=list)  # (time, fold)
    displaced_strand_history: List['DisplacedStrand'] = field(default_factory=list)

    # Copy number per region (for bias analysis)
    copy_number_array: Optional[np.ndarray] = None  # Copies per position


class Phi29Simulator:
    """
    Agent-based simulator for phi29 replication.

    Models individual replication forks as agents that:
    - Initiate at primer binding sites
    - Extend along template
    - Interact with other forks (collision)
    - Terminate at genome ends or collisions

    Mechanistic model integration:
        When config.use_mechanistic_model=True, the simulator uses
        MechanisticModel to calculate effects of reaction conditions
        (additives, temperature) on:
        - Extension rate (speed_factor)
        - Processivity (processivity_factor)
        - Primer binding rate (kon_factor)
        - Template accessibility (accessibility_factor)

    Probabilistic processivity:
        When config.probabilistic_processivity=True, polymerase
        dissociation follows a geometric distribution instead of
        a hard limit, giving more realistic fork length distributions.
    """

    def __init__(self,
                 primers: List[str],
                 primer_positions: Dict[str, Dict[str, List[int]]],  # primer -> {'forward': [...], 'reverse': [...]}
                 genome_length: int,
                 genome_sequence: str,
                 conditions: rc.ReactionConditions,
                 config: Optional[SimulationConfig] = None,
                 mechanistic_model: Optional['MechanisticModel'] = None):
        """
        Initialize simulator.

        Args:
            primers: Primer sequences
            primer_positions: Positions for each primer (forward/reverse)
            genome_length: Length of genome
            genome_sequence: Genome sequence (for GC-dependent speed)
            conditions: Reaction conditions
            config: Simulation configuration
            mechanistic_model: Optional pre-initialized MechanisticModel.
                              If not provided and config.use_mechanistic_model=True,
                              one will be created from conditions.
        """
        self.primers = primers
        self.primer_positions = primer_positions
        self.genome_length = genome_length
        self.genome_sequence = genome_sequence.upper()
        self.conditions = conditions
        self.config = config if config else SimulationConfig()

        # Calculate template GC content
        self.template_gc = self._calculate_template_gc()

        # Initialize mechanistic model if requested
        self.mechanistic_model = mechanistic_model
        self.primer_effects: Dict[str, 'MechanisticEffects'] = {}

        if self.config.use_mechanistic_model:
            self._initialize_mechanistic_model()

        # Get base polymerase parameters
        poly_params = get_polymerase_params(conditions.polymerase)
        base_extension_rate = poly_params['extension_rate']
        base_processivity = poly_params['processivity']
        optimal_temp = poly_params['optimal_temp']

        # Apply mechanistic effects to base parameters
        if self.config.use_mechanistic_model and self.mechanistic_model is not None:
            # Get enzyme activity from model (primer-independent)
            enzyme_params = self.mechanistic_model.get_enzyme_parameters()
            speed_factor = enzyme_params['speed_factor']
            processivity_factor = enzyme_params['processivity_factor']

            # Apply temperature effects if enabled
            if self.config.use_temperature_effects:
                temp_coef = MECHANISTIC_MODEL_PARAMS['enzyme']['temp_activity_coef']
                temp_dev = abs(conditions.temp - optimal_temp)
                temp_penalty = 1.0 - min(0.5, temp_coef * temp_dev)
                speed_factor *= temp_penalty

            self.config.extension_rate = base_extension_rate * speed_factor
            self.config.processivity_limit = int(base_processivity * processivity_factor)

            # Store factors for later use
            self._speed_factor = speed_factor
            self._processivity_factor = processivity_factor
        else:
            # Use base parameters without mechanistic modifications
            self.config.extension_rate = base_extension_rate
            self.config.processivity_limit = base_processivity
            self._speed_factor = 1.0
            self._processivity_factor = 1.0

        # Warn about low-processivity polymerases
        if conditions.polymerase in ('bst', 'klenow'):
            warnings.warn(
                f"{conditions.polymerase} has low processivity ({base_processivity:,} bp) "
                f"and is not well-suited for whole-genome amplification. "
                f"Consider phi29 (70kb) or equiphi29 (80kb) for SWGA.",
                UserWarning
            )

        # State
        self.current_time = 0.0
        self.forks: List[ReplicationFork] = []
        self.fork_id_counter = 0
        self.coverage = np.zeros(genome_length, dtype=bool)
        self.coverage_history = []

        # Strand displacement / Hyperbranching state
        self.displaced_strands: List[DisplacedStrand] = []
        self.displaced_strand_counter = 0
        self.copy_number = np.ones(genome_length, dtype=np.float32)  # Start with 1 copy
        self.total_bases_synthesized = 0
        self.amplification_history = []

        # Statistics
        self.forks_created = 0
        self.forks_terminated = 0
        self.forks_displaced = 0
        self.displacement_events = 0

        print(f"Phi29 Simulator initialized:")
        print(f"  Genome: {genome_length:,} bp (GC: {self.template_gc:.1%})")
        print(f"  Primers: {len(primers)}")
        print(f"  Polymerase: {conditions.polymerase}")
        print(f"  Extension rate: {self.config.extension_rate:.0f} bp/s")
        print(f"  Processivity: {self.config.processivity_limit:,} bp")
        print(f"  Duration: {self.config.duration:.0f} s ({self.config.duration/3600:.1f} hours)")
        if self.config.use_mechanistic_model:
            print(f"  Mechanistic model: enabled (speed={self._speed_factor:.2f}, proc={self._processivity_factor:.2f})")
        if self.config.probabilistic_processivity:
            print(f"  Processivity mode: probabilistic (geometric distribution)")
        if self.config.enable_strand_displacement:
            print(f"  Strand displacement: enabled (hyperbranching)")
            print(f"    Displacement probability: {self.config.displacement_probability:.0%}")
            print(f"    ssDNA binding boost: {self.config.displaced_strand_binding_boost:.1f}x")
            print(f"    Max amplification: {self.config.max_amplification_fold:.0e}x")

    def _calculate_template_gc(self) -> float:
        """Calculate overall GC content of template."""
        if not self.genome_sequence:
            return 0.5
        gc = self.genome_sequence.count('G') + self.genome_sequence.count('C')
        return gc / len(self.genome_sequence)

    def _initialize_mechanistic_model(self) -> None:
        """Initialize mechanistic model and pre-calculate primer effects."""
        if self.mechanistic_model is None:
            try:
                from neoswga.core.mechanistic_model import MechanisticModel
                self.mechanistic_model = MechanisticModel(self.conditions)
            except ImportError:
                warnings.warn(
                    "Could not import MechanisticModel. "
                    "Falling back to non-mechanistic simulation.",
                    UserWarning
                )
                self.config.use_mechanistic_model = False
                return

        # Pre-calculate effects for each primer
        for primer in self.primers:
            self.primer_effects[primer] = self.mechanistic_model.calculate_effects(
                primer, self.template_gc
            )

    def run(self, verbose: bool = True) -> SimulationResult:
        """
        Run simulation with strand displacement / hyperbranching.

        Args:
            verbose: Print progress

        Returns:
            SimulationResult with coverage and amplification metrics
        """
        print("\nStarting simulation...")
        if self.config.enable_strand_displacement:
            print("  Hyperbranching mode: enabled")

        # Initialize primers
        self._initialize_primers()

        # Main simulation loop
        next_record_time = 0.0
        resource_limited = False

        while self.current_time < self.config.duration:
            # Time step
            self.current_time += self.config.time_step

            # Update forks
            self._update_forks()

            # Attempt new primer initiations (stochastic)
            if random.random() < self.config.primer_binding_rate:
                self._attempt_primer_binding()

            # Additional binding attempts on displaced strands (hyperbranching)
            if self.config.enable_strand_displacement and self.displaced_strands:
                if random.random() < self.config.primer_binding_rate * 2:  # Higher rate for ssDNA
                    self._attempt_displaced_strand_binding()

            # Record metrics
            if self.current_time >= next_record_time:
                coverage_frac = self.coverage.sum() / self.genome_length
                self.coverage_history.append((self.current_time, coverage_frac))

                amplification = self._calculate_amplification_fold()
                self.amplification_history.append((self.current_time, amplification))

                next_record_time += self.config.record_interval

                if verbose and int(self.current_time) % 300 == 0:
                    active_forks = sum(1 for f in self.forks if f.state == ForkState.ACTIVE)
                    msg = (f"  t={self.current_time/60:.1f}min: "
                           f"Coverage={coverage_frac:.1%}, "
                           f"Active forks={active_forks}")
                    if self.config.enable_strand_displacement:
                        msg += f", Amplification={amplification:.1f}x"
                        msg += f", Displaced strands={len(self.displaced_strands)}"
                    print(msg)

            # Check for resource exhaustion
            if self._check_resource_limits():
                if verbose:
                    print(f"  Resource limit reached at t={self.current_time/60:.1f}min")
                    print(f"    Amplification: {self._calculate_amplification_fold():.1f}x")
                resource_limited = True
                break

            # Early termination if no active forks
            active = sum(1 for f in self.forks if f.state == ForkState.ACTIVE)
            if active == 0 and self.current_time > 60:
                if verbose:
                    print(f"  No active forks at t={self.current_time/60:.1f}min, terminating early")
                break

        final_coverage = self.coverage.sum() / self.genome_length
        final_amplification = self._calculate_amplification_fold()

        if verbose:
            print(f"\nSimulation complete:")
            print(f"  Final coverage: {final_coverage:.1%}")
            print(f"  Forks created: {self.forks_created}")
            print(f"  Forks terminated: {self.forks_terminated}")
            if self.config.enable_strand_displacement:
                print(f"  Displacement events: {self.displacement_events}")
                print(f"  Displaced strands: {len(self.displaced_strands)}")
                print(f"  Total bases synthesized: {self.total_bases_synthesized:,}")
                print(f"  Final amplification: {final_amplification:.1f}x")
                if resource_limited:
                    print(f"  Note: Simulation ended due to resource limits")

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
            fork_history=self.forks.copy(),
            # Hyperbranching metrics
            amplification_fold=final_amplification,
            total_bases_synthesized=self.total_bases_synthesized,
            num_displaced_strands=len(self.displaced_strands),
            amplification_over_time=self.amplification_history,
            displaced_strand_history=self.displaced_strands.copy(),
            copy_number_array=self.copy_number.copy(),
        )

    def _initialize_primers(self):
        """Create initial primer bindings."""
        # For each primer binding site, stochastically initiate
        for primer in self.primers:
            positions = self.primer_positions.get(primer, {'forward': [], 'reverse': []})

            # Forward strand
            for pos in positions.get('forward', []):
                # Binding probability based on Tm and mechanistic effects
                tm = self.conditions.calculate_effective_tm(primer)
                binding_prob = self._calculate_binding_probability(tm, primer)

                if random.random() < binding_prob * 0.1:  # Scale down for initial binding
                    self._create_fork(primer, pos, 'forward', 'top')

            # Reverse strand
            for pos in positions.get('reverse', []):
                tm = self.conditions.calculate_effective_tm(primer)
                binding_prob = self._calculate_binding_probability(tm, primer)

                if random.random() < binding_prob * 0.1:
                    self._create_fork(primer, pos, 'reverse', 'bottom')

    def _calculate_binding_probability(self, tm: float, primer: Optional[str] = None) -> float:
        """
        Calculate primer binding probability based on Tm and mechanistic effects.

        Args:
            tm: Effective melting temperature
            primer: Optional primer sequence (for mechanistic model lookup)

        Returns:
            Binding probability (0-1)
        """
        # Base probability from Tm vs reaction temperature
        delta = tm - self.conditions.temp
        base_prob = 1 / (1 + np.exp(-delta / 5))  # Sharper transition

        # Apply mechanistic effects if available
        if (self.config.use_mechanistic_model and
                primer is not None and
                primer in self.primer_effects):
            effects = self.primer_effects[primer]
            # Modify by kon_factor and accessibility_factor
            base_prob *= effects.kon_factor * effects.accessibility_factor

        return min(1.0, max(0.0, base_prob))

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
        """Update all replication forks, including strand displacement."""
        for fork in self.forks:
            if fork.state != ForkState.ACTIVE:
                continue

            # Check resource limits before processing
            if self._check_resource_limits():
                self._terminate_fork(fork, "resources")
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
                # With strand displacement, collision can lead to displacement
                if self.config.enable_strand_displacement:
                    if random.random() < self.config.displacement_probability:
                        # Displacing fork continues, other fork is displaced
                        collision.displaced_by = fork.fork_id
                        self._terminate_fork(collision, "displaced")

                        # Create displaced strand from the collided fork's synthesis
                        if collision.bases_synthesized >= self.config.min_displaced_strand_length:
                            coll_start = min(collision.start_position, collision.current_position)
                            coll_end = max(collision.start_position, collision.current_position)
                            self._create_displaced_strand(
                                coll_start, coll_end,
                                collision.strand,
                                collision.template_id
                            )
                    else:
                        # Traditional collision - both terminate
                        self._terminate_fork(fork, "collision")
                        self._terminate_fork(collision, "collision")
                        continue
                else:
                    self._terminate_fork(fork, "collision")
                    self._terminate_fork(collision, "collision")
                    continue

            # Update position and coverage
            start = min(fork.current_position, new_position)
            end = max(fork.current_position, new_position)

            # Check for strand displacement opportunity
            # (extending through previously replicated region)
            if self._check_displacement_opportunity(fork, start, end):
                # Create a displaced strand for the region being overwritten
                opposite_strand = 'bottom' if fork.strand == 'top' else 'top'
                self._create_displaced_strand(start, end, opposite_strand, fork.template_id)

            # Update coverage and copy number
            self.coverage[start:end] = True
            self.copy_number[start:end] += 1.0

            # Track total synthesis
            self.total_bases_synthesized += distance

            fork.current_position = new_position
            fork.bases_synthesized += distance

            # Check processivity - probabilistic or deterministic
            if self._should_terminate_processivity(fork, distance):
                self._terminate_fork(fork, "processivity")
                continue

            # Update speed based on new position
            fork.speed = self._calculate_local_speed(new_position)

    def _should_terminate_processivity(self, fork: ReplicationFork, step_size: int) -> bool:
        """
        Determine if fork should terminate due to processivity.

        Uses either probabilistic (geometric distribution) or deterministic
        (hard limit) model based on config.

        Probabilistic model:
            For each step, probability of continuing is:
            P(continue) = exp(-step_size / processivity)

            This gives a geometric distribution with mean equal to
            the nominal processivity.

            For phi29 with processivity=70,000 bp and step_size=150 bp:
            P(continue per step) = exp(-150/70000) = 0.998
            After 70,000 bp: P(still running) = 0.998^467 = 0.37 (37%)

        Args:
            fork: The replication fork
            step_size: Distance traveled this step

        Returns:
            True if fork should terminate, False otherwise
        """
        if self.config.probabilistic_processivity:
            # Probabilistic termination using geometric model
            # Probability of continuing this step
            continue_prob = math.exp(-step_size / self.config.processivity_limit)
            return random.random() > continue_prob
        else:
            # Deterministic hard limit
            return fork.bases_synthesized >= self.config.processivity_limit

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
        """Terminate replication fork.

        Args:
            fork: Fork to terminate
            reason: Termination reason ('boundary', 'collision', 'processivity', 'displaced', 'resources')
        """
        if reason == "collision":
            fork.state = ForkState.COLLIDED
        elif reason == "displaced":
            fork.state = ForkState.DISPLACED
            self.forks_displaced += 1
        else:
            fork.state = ForkState.TERMINATED
        fork.termination_time = self.current_time
        self.forks_terminated += 1

    # =========================================================================
    # Strand Displacement / Hyperbranching Methods
    # =========================================================================

    def _check_displacement_opportunity(
        self,
        fork: ReplicationFork,
        start_pos: int,
        end_pos: int
    ) -> bool:
        """
        Check if extending through a region creates a displacement event.

        A displacement occurs when:
        1. The fork extends through a region that was already replicated
        2. The region has copy_number > 1 (has been synthesized before)

        Args:
            fork: The extending fork
            start_pos: Start of the region being replicated
            end_pos: End of the region being replicated

        Returns:
            True if displacement should occur
        """
        if not self.config.enable_strand_displacement:
            return False

        # Check if any part of this region has been replicated before
        region_copy_number = self.copy_number[start_pos:end_pos]
        if np.any(region_copy_number > 1):
            # Probabilistic displacement
            return random.random() < self.config.displacement_probability

        return False

    def _create_displaced_strand(
        self,
        start_pos: int,
        end_pos: int,
        strand: str,
        parent_template_id: int = 0
    ) -> Optional[DisplacedStrand]:
        """
        Create a new displaced strand that can serve as template.

        Args:
            start_pos: Start position on genome
            end_pos: End position on genome
            strand: 'top' or 'bottom'
            parent_template_id: Which template this was displaced from

        Returns:
            DisplacedStrand object, or None if too short or at limit
        """
        length = end_pos - start_pos

        # Don't track very short displacements
        if length < self.config.min_displaced_strand_length:
            return None

        # Check memory limit
        if len(self.displaced_strands) >= self.config.max_displaced_strands:
            return None

        displaced = DisplacedStrand(
            strand_id=self.displaced_strand_counter,
            start_pos=start_pos,
            end_pos=end_pos,
            strand=strand,
            template_id=parent_template_id,
            creation_time=self.current_time,
        )

        self.displaced_strands.append(displaced)
        self.displaced_strand_counter += 1
        self.displacement_events += 1

        return displaced

    def _get_available_displaced_strands(self, position: int) -> List[DisplacedStrand]:
        """
        Get displaced strands that cover a given position and are available.

        Args:
            position: Position to check

        Returns:
            List of available displaced strands at this position
        """
        return [
            ds for ds in self.displaced_strands
            if ds.is_available and ds.contains_position(position)
        ]

    def _calculate_amplification_fold(self) -> float:
        """
        Calculate current amplification fold.

        Returns:
            Total DNA synthesized / original genome length
        """
        return self.total_bases_synthesized / self.genome_length + 1.0

    def _check_resource_limits(self) -> bool:
        """
        Check if resource limits have been reached.

        Returns:
            True if resources exhausted (should stop simulation)
        """
        # Check amplification limit (dNTP exhaustion)
        if self._calculate_amplification_fold() >= self.config.max_amplification_fold:
            return True

        # Check concurrent fork limit (polymerase availability)
        active_forks = sum(1 for f in self.forks if f.state == ForkState.ACTIVE)
        if active_forks >= self.config.max_concurrent_forks:
            return True

        return False

    def _attempt_primer_binding(self):
        """
        Stochastically attempt new primer binding.

        With strand displacement enabled, primers can bind to:
        1. Original template (double-stranded, position not yet covered)
        2. Displaced strands (single-stranded, higher binding affinity)

        Single-stranded displaced DNA binds primers more efficiently
        (no need to displace the complementary strand), controlled by
        config.displaced_strand_binding_boost.
        """
        # Check resource limits
        if self._check_resource_limits():
            return

        # Random primer
        primer = random.choice(self.primers)
        positions = self.primer_positions.get(primer, {'forward': [], 'reverse': []})

        # Calculate binding probability
        tm = self.conditions.calculate_effective_tm(primer)
        base_binding_prob = self._calculate_binding_probability(tm, primer)

        # Apply mechanistic kon_factor to binding rate if available
        effective_binding_rate = self.config.primer_binding_rate
        if (self.config.use_mechanistic_model and
                primer in self.primer_effects):
            effective_binding_rate *= self.primer_effects[primer].kon_factor

        # Try binding to original template
        all_positions = (
            [(p, 'forward', 'top') for p in positions.get('forward', [])] +
            [(p, 'reverse', 'bottom') for p in positions.get('reverse', [])]
        )

        if all_positions:
            pos, direction, strand = random.choice(all_positions)

            # Original template binding (only if not yet covered)
            if not self.coverage[pos]:
                if random.random() < base_binding_prob * effective_binding_rate:
                    self._create_fork(primer, pos, direction, strand)
                    return

        # Try binding to displaced strands (if enabled)
        if self.config.enable_strand_displacement and self.displaced_strands:
            # Get available displaced strands
            available_strands = [ds for ds in self.displaced_strands if ds.is_available]

            if available_strands:
                # Pick a random displaced strand
                ds = random.choice(available_strands)

                # Pick a random position within the strand
                if ds.length > len(primer):
                    pos = random.randint(ds.start_pos, ds.end_pos - len(primer))

                    # Enhanced binding to single-stranded DNA
                    ssdna_binding_prob = base_binding_prob * self.config.displaced_strand_binding_boost
                    ssdna_binding_prob = min(1.0, ssdna_binding_prob)

                    if random.random() < ssdna_binding_prob * effective_binding_rate:
                        # Determine direction based on strand
                        # On displaced top strand, new synthesis goes reverse (5'->3' on complement)
                        # On displaced bottom strand, new synthesis goes forward
                        if ds.strand == 'top':
                            direction = 'reverse'
                            new_strand = 'bottom'
                        else:
                            direction = 'forward'
                            new_strand = 'top'

                        # Create fork on the displaced strand
                        fork = self._create_fork(primer, pos, direction, new_strand)
                        if fork is not None:
                            # Mark as on a displaced template
                            self.forks[-1].template_id = ds.strand_id + 1  # +1 to distinguish from original

    def _attempt_displaced_strand_binding(self):
        """
        Additional binding attempts specifically for displaced strands.

        Called during simulation to model hyperbranching where displaced
        strands quickly become templates for new synthesis.
        """
        if not self.config.enable_strand_displacement:
            return

        if not self.displaced_strands:
            return

        # Get available displaced strands
        available = [ds for ds in self.displaced_strands if ds.is_available]
        if not available:
            return

        # Try each primer on displaced strands
        for primer in self.primers:
            positions = self.primer_positions.get(primer, {'forward': [], 'reverse': []})

            tm = self.conditions.calculate_effective_tm(primer)
            binding_prob = self._calculate_binding_probability(tm, primer)
            binding_prob *= self.config.displaced_strand_binding_boost
            binding_prob = min(1.0, binding_prob)

            effective_rate = self.config.primer_binding_rate * 0.1  # Lower rate for displaced
            if (self.config.use_mechanistic_model and primer in self.primer_effects):
                effective_rate *= self.primer_effects[primer].kon_factor

            for ds in available:
                # Check if primer can bind to this displaced strand
                # (simplified: check if primer positions overlap with strand)
                for pos in positions.get('forward', []) + positions.get('reverse', []):
                    if ds.contains_position(pos):
                        if random.random() < binding_prob * effective_rate:
                            direction = 'forward' if ds.strand == 'bottom' else 'reverse'
                            new_strand = 'top' if ds.strand == 'bottom' else 'bottom'
                            self._create_fork(primer, pos, direction, new_strand)
                            return  # One binding per call


def simulate_primer_set(primers: List[str],
                       primer_positions: Dict,
                       genome_length: int,
                       genome_sequence: str,
                       conditions: Optional[rc.ReactionConditions] = None,
                       n_replicates: int = 5,
                       config: Optional[SimulationConfig] = None,
                       mechanistic_model: Optional['MechanisticModel'] = None) -> Dict:
    """
    Simulate primer set performance with multiple replicates.

    Args:
        primers: Primer sequences
        primer_positions: Primer positions dict
        genome_length: Genome length
        genome_sequence: Genome sequence
        conditions: Reaction conditions
        n_replicates: Number of simulation replicates
        config: Optional SimulationConfig (if not provided, uses defaults)
        mechanistic_model: Optional pre-initialized MechanisticModel

    Returns:
        Dictionary with mean and std of metrics
    """
    if conditions is None:
        conditions = rc.get_enhanced_conditions()

    # Create mechanistic model once and reuse
    if config is None:
        config = SimulationConfig()

    if mechanistic_model is None and config.use_mechanistic_model:
        try:
            from neoswga.core.mechanistic_model import MechanisticModel
            mechanistic_model = MechanisticModel(conditions)
        except ImportError:
            pass

    results = []

    print(f"Running {n_replicates} simulation replicates...")

    for i in range(n_replicates):
        print(f"\n{'='*60}")
        print(f"Replicate {i+1}/{n_replicates}")
        print(f"{'='*60}")

        simulator = Phi29Simulator(
            primers, primer_positions, genome_length,
            genome_sequence, conditions,
            config=config,
            mechanistic_model=mechanistic_model
        )

        result = simulator.run(verbose=False)
        results.append(result)

    # Aggregate statistics
    coverages = [r.final_coverage_fraction for r in results]
    amplifications = [r.amplification_fold for r in results]

    print(f"\n{'='*60}")
    print("Aggregate Results:")
    print(f"{'='*60}")
    print(f"Coverage: {np.mean(coverages):.1%} +/- {np.std(coverages):.1%}")
    print(f"Mean forks created: {np.mean([r.num_forks_created for r in results]):.0f}")
    print(f"Mean fork travel: {np.mean([r.mean_fork_travel for r in results]):.0f} bp")

    # Hyperbranching metrics
    if config.enable_strand_displacement:
        print(f"Amplification: {np.mean(amplifications):.1f}x +/- {np.std(amplifications):.1f}x")
        print(f"Mean displaced strands: {np.mean([r.num_displaced_strands for r in results]):.0f}")
        print(f"Mean bases synthesized: {np.mean([r.total_bases_synthesized for r in results]):,.0f}")

    return {
        'mean_coverage': np.mean(coverages),
        'std_coverage': np.std(coverages),
        'mean_forks_created': np.mean([r.num_forks_created for r in results]),
        'mean_fork_travel': np.mean([r.mean_fork_travel for r in results]),
        # Hyperbranching metrics
        'mean_amplification': np.mean(amplifications),
        'std_amplification': np.std(amplifications),
        'mean_displaced_strands': np.mean([r.num_displaced_strands for r in results]),
        'mean_bases_synthesized': np.mean([r.total_bases_synthesized for r in results]),
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
