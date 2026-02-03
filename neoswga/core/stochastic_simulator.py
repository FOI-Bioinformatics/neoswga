"""
Stochastic simulation of SWGA amplification dynamics.

Uses Gillespie algorithm to model actual molecular reactions:
- Primer binding (stochastic)
- Extension (distance-dependent)
- Rebinding (network effect)
- Resource depletion (dNTP, polymerase)

This VALIDATES network-based predictions against actual kinetics.

Mechanistic model integration:
    When use_mechanistic_model=True, the simulator uses MechanisticModel
    to modify kinetic parameters based on reaction conditions (additives,
    temperature, polymerase):
    - kon: Modified by kon_factor (binding kinetics)
    - koff: Modified by koff_factor (dissociation kinetics)
    - extension_rate: Modified by speed_factor
    - processivity_step: Modified by processivity_factor
"""

import math
import numpy as np
import logging
import warnings
from typing import List, Dict, Tuple, Optional, TYPE_CHECKING
from dataclasses import dataclass, field
from collections import defaultdict
import time

from neoswga.core.mechanistic_params import get_polymerase_params, MECHANISTIC_MODEL_PARAMS

if TYPE_CHECKING:
    from neoswga.core.mechanistic_model import MechanisticModel, MechanisticEffects
    from neoswga.core.reaction_conditions import ReactionConditions

logger = logging.getLogger(__name__)


@dataclass
class ReactionState:
    """Current state of amplification reaction"""
    time: float = 0.0  # seconds

    # Molecule counts
    fg_molecules: int = 1000  # Target genome copies
    bg_molecules: int = 999000  # Background genome copies

    # Bound primers (per binding site)
    fg_bound: Dict[int, int] = None  # site -> count
    bg_bound: Dict[int, int] = None

    # Extended products
    fg_products: int = 0
    bg_products: int = 0

    # Resources
    primer_conc: float = 1e-6  # M (1 µM)
    dNTP: float = 200e-6  # M (200 µM)
    polymerase: float = 10e-9  # M (10 nM Phi29)

    def __post_init__(self):
        if self.fg_bound is None:
            self.fg_bound = defaultdict(int)
        if self.bg_bound is None:
            self.bg_bound = defaultdict(int)


@dataclass
class ReactionParameters:
    """Kinetic parameters for SWGA.

    Processivity model:
        Literature phi29 processivity is ~70,000 bp (Blanco 1989).
        For stochastic model, we use per-step probability:
            processivity_step = exp(-step_size / processivity)
        With step_size=100bp and processivity=70,000bp:
            processivity_step = exp(-100/70000) = 0.99857
        This gives mean run length = step_size / -ln(processivity_step) = 70,000 bp

    Mechanistic model integration:
        When created via from_conditions(), these parameters are modified
        by MechanisticEffects:
        - k_on: multiplied by kon_factor
        - k_off: multiplied by koff_factor
        - extension_rate: multiplied by speed_factor
        - processivity_step: adjusted by processivity_factor
    """
    # Binding kinetics
    k_on: float = 1e6  # M^-1 s^-1 (primer binding)
    k_off: float = 1e-2  # s^-1 (dissociation)

    # Extension kinetics
    extension_rate: float = 150.0  # nt/s (Phi29 median, literature 100-170)
    processivity_step: float = 0.99857  # Per 100bp step (phi29 default)
    max_extension: int = 70000  # Maximum extension length (bp)

    # Resource consumption
    dNTP_per_nt: float = 4.0  # Average dNTPs per nucleotide
    polymerase_binding: float = 1e7  # M^-1 s^-1

    # Network effects
    rebinding_rate: float = 10.0  # Enhanced binding to amplicons

    # Temperature-dependent
    temperature: float = 30.0  # C

    # Mechanistic model flags
    use_mechanistic_model: bool = True

    def extension_time(self, distance: int) -> float:
        """Time to extend given distance."""
        return distance / self.extension_rate

    def extension_probability(self, distance: int) -> float:
        """Probability of completing extension to given distance.

        Uses geometric model where processivity_step is the probability
        of continuing extension per 100bp step.

        For phi29 (processivity_step=0.99857):
            - At 1 kb: P = 0.99857^10 = 0.986 (98.6% complete)
            - At 10 kb: P = 0.99857^100 = 0.867 (86.7% complete)
            - At 50 kb: P = 0.99857^500 = 0.489 (48.9% complete)
            - At 70 kb: P = 0.99857^700 = 0.368 (36.8% complete)

        Mean run length = 100 / -ln(0.99857) = ~70,000 bp
        """
        steps = distance / 100  # 100 bp per step
        return self.processivity_step ** steps

    @classmethod
    def from_conditions(
        cls,
        conditions: 'ReactionConditions',
        template_gc: float = 0.5,
        use_mechanistic_model: bool = True,
    ) -> 'ReactionParameters':
        """
        Create ReactionParameters from ReactionConditions with mechanistic model.

        Args:
            conditions: ReactionConditions with temp, additives, polymerase
            template_gc: Template genome GC content (0-1)
            use_mechanistic_model: Whether to apply mechanistic effects

        Returns:
            ReactionParameters with literature-based values modified by
            mechanistic model effects if enabled.
        """
        # Get base polymerase parameters
        poly_params = get_polymerase_params(conditions.polymerase)
        base_extension_rate = poly_params['extension_rate']
        base_processivity = poly_params['processivity']
        base_processivity_step = poly_params['processivity_step']
        optimal_temp = poly_params['optimal_temp']

        # Warn about low-processivity polymerases
        if conditions.polymerase in ('bst', 'klenow'):
            warnings.warn(
                f"{conditions.polymerase} has low processivity ({base_processivity:,} bp) "
                f"and is not well-suited for whole-genome amplification. "
                f"Consider phi29 (70kb) or equiphi29 (80kb) for SWGA.",
                UserWarning
            )

        # Start with base parameters
        params = cls(
            extension_rate=base_extension_rate,
            processivity_step=base_processivity_step,
            max_extension=base_processivity,
            temperature=conditions.temp,
            use_mechanistic_model=use_mechanistic_model,
        )

        if not use_mechanistic_model:
            return params

        # Apply mechanistic model effects
        try:
            from neoswga.core.mechanistic_model import MechanisticModel

            model = MechanisticModel(conditions)

            # Get enzyme parameters (primer-independent)
            enzyme_params = model.get_enzyme_parameters()
            speed_factor = enzyme_params['speed_factor']
            processivity_factor = enzyme_params['processivity_factor']

            # Apply temperature deviation penalty
            temp_coef = MECHANISTIC_MODEL_PARAMS['enzyme']['temp_activity_coef']
            temp_dev = abs(conditions.temp - optimal_temp)
            temp_penalty = 1.0 - min(0.5, temp_coef * temp_dev)
            speed_factor *= temp_penalty

            # Modify extension rate
            params.extension_rate = base_extension_rate * speed_factor

            # Modify processivity
            # processivity_step = exp(-100 / effective_processivity)
            effective_processivity = base_processivity * processivity_factor
            params.processivity_step = math.exp(-100.0 / effective_processivity)
            params.max_extension = int(effective_processivity)

            # Get kinetics parameters from model (use a representative primer)
            # We use a balanced 10-mer as reference
            ref_primer = 'ATCGATCGAT'
            effects = model.calculate_effects(ref_primer, template_gc)

            # Modify kon and koff
            params.k_on = 1e6 * effects.kon_factor
            params.k_off = 1e-2 * effects.koff_factor

            logger.info(
                f"Mechanistic model applied: "
                f"speed={speed_factor:.2f}, proc={processivity_factor:.2f}, "
                f"kon={effects.kon_factor:.2f}, koff={effects.koff_factor:.2f}"
            )

        except ImportError:
            warnings.warn(
                "Could not import MechanisticModel. "
                "Using base parameters without mechanistic corrections.",
                UserWarning
            )

        return params


class GillespieSimulator:
    """
    Gillespie algorithm for stochastic simulation of SWGA.

    Reactions modeled:
    1. P + T -> PT        (primer binding to template)
    2. PT -> P + T        (primer dissociation)
    3. PT + Pol -> EXT    (extension initiation)
    4. EXT -> 2T          (completion of extension)

    Mechanistic model integration:
        When conditions are provided (or params created via from_conditions),
        kinetic parameters are modified by MechanisticEffects:
        - k_on, k_off: Binding kinetics
        - extension_rate: Polymerase speed
        - processivity_step: Continuation probability
    """

    def __init__(self, primers: List[str], fg_network, bg_network,
                 params: Optional[ReactionParameters] = None,
                 conditions: Optional['ReactionConditions'] = None,
                 template_gc: float = 0.5):
        """
        Initialize simulator.

        Args:
            primers: List of primer sequences
            fg_network: AmplificationNetwork for target
            bg_network: AmplificationNetwork for background
            params: Reaction parameters (or use defaults). If not provided
                   and conditions is provided, will create from conditions.
            conditions: Optional ReactionConditions for mechanistic model.
                       If provided and params is None, will create params
                       using ReactionParameters.from_conditions().
            template_gc: Template genome GC content (0-1), used for
                        mechanistic model calculations.
        """
        self.primers = primers
        self.fg_network = fg_network
        self.bg_network = bg_network
        self.template_gc = template_gc

        # Create params from conditions if needed
        if params is not None:
            self.params = params
        elif conditions is not None:
            self.params = ReactionParameters.from_conditions(
                conditions, template_gc=template_gc
            )
        else:
            self.params = ReactionParameters()

        # Extract binding sites
        self.fg_sites = list(fg_network.binding_sites)
        self.bg_sites = list(bg_network.binding_sites)

        # Build extension network (site -> reachable sites)
        self.fg_extensions = self._build_extension_network(fg_network)
        self.bg_extensions = self._build_extension_network(bg_network)

        logger.info(
            f"GillespieSimulator initialized: "
            f"kon={self.params.k_on:.2e}, koff={self.params.k_off:.2e}, "
            f"extension={self.params.extension_rate:.0f} nt/s, "
            f"processivity_step={self.params.processivity_step:.5f}"
        )

    def _build_extension_network(self, network) -> Dict[int, List[Tuple[int, int]]]:
        """
        Build extension network: site -> [(target_site, distance), ...]
        """
        extensions = defaultdict(list)

        for site1, site2, data in network.graph.edges(data=True):
            site1_idx = self.fg_sites.index(site1) if site1 in self.fg_sites else self.bg_sites.index(site1)
            site2_idx = self.fg_sites.index(site2) if site2 in self.fg_sites else self.bg_sites.index(site2)
            distance = data.get('distance', 0)

            extensions[site1_idx].append((site2_idx, distance))

        return dict(extensions)

    def simulate(self, max_time: float = 28800.0, sample_interval: float = 60.0) -> List[Dict]:
        """
        Run Gillespie simulation.

        Args:
            max_time: Maximum simulation time (seconds). Default 8 hours (28800s)
                      to match typical SWGA protocol duration.
            sample_interval: Sampling interval for recording state

        Returns:
            List of state snapshots
        """
        state = ReactionState()
        history = [self._snapshot(state)]
        next_sample = sample_interval

        logger.info(f"Starting Gillespie simulation: {max_time/3600:.1f} hours")

        iteration = 0
        while state.time < max_time:
            iteration += 1

            # Calculate reaction propensities
            reactions = self._get_reactions(state)

            if not reactions:
                logger.warning("No reactions possible - simulation terminated")
                break

            total_propensity = sum(r[1] for r in reactions)

            if total_propensity == 0:
                logger.warning("Zero total propensity - simulation terminated")
                break

            # Time to next reaction (exponential distribution)
            dt = np.random.exponential(1.0 / total_propensity)
            state.time += dt

            # Select which reaction occurs
            rand = np.random.random() * total_propensity
            cumsum = 0
            selected_reaction = None

            for reaction_type, propensity, *args in reactions:
                cumsum += propensity
                if rand < cumsum:
                    selected_reaction = (reaction_type, *args)
                    break

            if selected_reaction:
                self._execute_reaction(state, selected_reaction)

            # Sample state
            if state.time >= next_sample:
                history.append(self._snapshot(state))
                next_sample += sample_interval

                if iteration % 1000 == 0:
                    logger.debug(f"Time: {state.time/60:.1f} min, "
                               f"FG products: {state.fg_products}, "
                               f"BG products: {state.bg_products}")

        logger.info(f"Simulation complete: {iteration} reactions in {state.time/60:.1f} min")

        return history

    def _get_reactions(self, state: ReactionState) -> List[Tuple]:
        """
        Calculate all possible reactions and their propensities.

        Returns:
            List of (reaction_type, propensity, *args)
        """
        reactions = []

        # 1. Primer binding to target sites
        for i, site in enumerate(self.fg_sites):
            # Binding rate depends on primer concentration and template availability
            templates_available = state.fg_molecules + state.fg_products
            if templates_available > 0:
                propensity = self.params.k_on * state.primer_conc * templates_available
                reactions.append(('fg_bind', propensity, i))

        # 2. Primer binding to background sites
        for i, site in enumerate(self.bg_sites):
            templates_available = state.bg_molecules + state.bg_products
            if templates_available > 0:
                propensity = self.params.k_on * state.primer_conc * templates_available * 0.1  # Lower affinity
                reactions.append(('bg_bind', propensity, i))

        # 3. Extension from bound primers (target)
        for site_idx, count in state.fg_bound.items():
            if count > 0 and site_idx in self.fg_extensions:
                # Can extend to any connected site
                for target_site, distance in self.fg_extensions[site_idx]:
                    if distance <= self.params.max_extension:
                        # Extension rate depends on polymerase availability
                        propensity = self.params.polymerase_binding * self.params.polymerase * count
                        propensity *= self.params.extension_probability(distance)
                        reactions.append(('fg_extend', propensity, site_idx, target_site, distance))

        # 4. Extension from bound primers (background)
        for site_idx, count in state.bg_bound.items():
            if count > 0 and site_idx in self.bg_extensions:
                for target_site, distance in self.bg_extensions[site_idx]:
                    if distance <= self.params.max_extension:
                        propensity = self.params.polymerase_binding * self.params.polymerase * count * 0.5
                        propensity *= self.params.extension_probability(distance)
                        reactions.append(('bg_extend', propensity, site_idx, target_site, distance))

        # 5. Primer dissociation
        for site_idx, count in state.fg_bound.items():
            if count > 0:
                propensity = self.params.k_off * count
                reactions.append(('fg_unbind', propensity, site_idx))

        for site_idx, count in state.bg_bound.items():
            if count > 0:
                propensity = self.params.k_off * count
                reactions.append(('bg_unbind', propensity, site_idx))

        return reactions

    def _execute_reaction(self, state: ReactionState, reaction: Tuple):
        """Execute selected reaction, update state"""
        reaction_type = reaction[0]

        if reaction_type == 'fg_bind':
            site_idx = reaction[1]
            state.fg_bound[site_idx] += 1

        elif reaction_type == 'bg_bind':
            site_idx = reaction[1]
            state.bg_bound[site_idx] += 1

        elif reaction_type == 'fg_extend':
            site_idx, target_site, distance = reaction[1:]

            # Consume resources
            dNTPs_needed = distance * self.params.dNTP_per_nt
            if state.dNTP >= dNTPs_needed * 1e-9:  # Check sufficient dNTPs
                state.dNTP -= dNTPs_needed * 1e-9

                # Create new amplicon
                state.fg_products += 1

                # Unbind primer
                state.fg_bound[site_idx] -= 1

                # New amplicon can template more synthesis (network effect!)
                # Enhanced binding to products
                state.fg_bound[target_site] += 1  # Immediate rebinding

        elif reaction_type == 'bg_extend':
            site_idx, target_site, distance = reaction[1:]

            dNTPs_needed = distance * self.params.dNTP_per_nt
            if state.dNTP >= dNTPs_needed * 1e-9:
                state.dNTP -= dNTPs_needed * 1e-9
                state.bg_products += 1
                state.bg_bound[site_idx] -= 1

        elif reaction_type == 'fg_unbind':
            site_idx = reaction[1]
            state.fg_bound[site_idx] = max(0, state.fg_bound[site_idx] - 1)

        elif reaction_type == 'bg_unbind':
            site_idx = reaction[1]
            state.bg_bound[site_idx] = max(0, state.bg_bound[site_idx] - 1)

    def _snapshot(self, state: ReactionState) -> Dict:
        """Create snapshot of current state"""
        total_fg = state.fg_molecules + state.fg_products
        total_bg = state.bg_molecules + state.bg_products

        enrichment = total_fg / total_bg if total_bg > 0 else 0
        fg_fold = total_fg / state.fg_molecules if state.fg_molecules > 0 else 0
        bg_fold = total_bg / state.bg_molecules if state.bg_molecules > 0 else 0

        return {
            'time': state.time,
            'fg_molecules': state.fg_molecules,
            'fg_products': state.fg_products,
            'total_fg': total_fg,
            'bg_molecules': state.bg_molecules,
            'bg_products': state.bg_products,
            'total_bg': total_bg,
            'enrichment': enrichment,
            'fg_amplification': fg_fold,
            'bg_amplification': bg_fold,
            'dNTP': state.dNTP,
            'primer_conc': state.primer_conc,
        }


def validate_network_predictions(
    primers: List[str],
    fg_network,
    bg_network,
    conditions: Optional['ReactionConditions'] = None,
    template_gc: float = 0.5,
) -> Dict:
    """
    Validate network-based predictions with stochastic simulation.

    Compares:
    - Predicted amplification (from network connectivity)
    - Simulated amplification (from Gillespie)

    Args:
        primers: List of primer sequences
        fg_network: AmplificationNetwork for target
        bg_network: AmplificationNetwork for background
        conditions: Optional ReactionConditions for mechanistic model
        template_gc: Template genome GC content (0-1)

    Returns:
        Validation results
    """
    logger.info("Validating network predictions with stochastic simulation...")

    # Network-based prediction
    fg_largest = fg_network.largest_component_size()
    bg_avg = bg_network.average_component_size()

    predicted_fg_fold = fg_network.predict_amplification_fold()
    predicted_bg_fold = bg_network.predict_amplification_fold()
    predicted_enrichment = predicted_fg_fold / (predicted_bg_fold + 1)

    logger.info(f"Network prediction: {predicted_enrichment:.0f}x enrichment")

    # Stochastic simulation (8 hours typical SWGA protocol)
    simulator = GillespieSimulator(
        primers, fg_network, bg_network,
        conditions=conditions,
        template_gc=template_gc
    )
    history = simulator.simulate(max_time=28800.0, sample_interval=60.0)  # 8 hours

    final_state = history[-1]
    simulated_enrichment = final_state['enrichment']
    simulated_fg_fold = final_state['fg_amplification']
    simulated_bg_fold = final_state['bg_amplification']

    logger.info(f"Simulation result: {simulated_enrichment:.0f}x enrichment")

    # Compare
    prediction_error = abs(predicted_enrichment - simulated_enrichment) / max(simulated_enrichment, 1e-10)

    return {
        'predicted': {
            'enrichment': predicted_enrichment,
            'fg_amplification': predicted_fg_fold,
            'bg_amplification': predicted_bg_fold,
            'fg_largest_component': fg_largest,
            'bg_avg_component': bg_avg,
        },
        'simulated': {
            'enrichment': simulated_enrichment,
            'fg_amplification': simulated_fg_fold,
            'bg_amplification': simulated_bg_fold,
        },
        'validation': {
            'prediction_error': prediction_error,
            'prediction_accurate': prediction_error < 0.5,  # Within 50%
        },
        'history': history,
        'mechanistic_model_used': conditions is not None,
    }


def plot_simulation_results(history: List[Dict], output_path: Optional[str] = None):
    """
    Plot simulation results.

    Args:
        history: Simulation history
        output_path: Path to save plot (or None to display)
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning("matplotlib not available, skipping plot")
        return

    times = [h['time'] / 3600 for h in history]  # Hours
    fg_amp = [h['fg_amplification'] for h in history]
    bg_amp = [h['bg_amplification'] for h in history]
    enrichment = [h['enrichment'] for h in history]

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Target amplification
    axes[0, 0].semilogy(times, fg_amp, 'b-', linewidth=2)
    axes[0, 0].set_xlabel('Time (hours)')
    axes[0, 0].set_ylabel('Target amplification (fold)')
    axes[0, 0].set_title('Target Genome Amplification')
    axes[0, 0].grid(True, alpha=0.3)

    # Background amplification
    axes[0, 1].semilogy(times, bg_amp, 'r-', linewidth=2)
    axes[0, 1].set_xlabel('Time (hours)')
    axes[0, 1].set_ylabel('Background amplification (fold)')
    axes[0, 1].set_title('Background Genome Amplification')
    axes[0, 1].grid(True, alpha=0.3)

    # Enrichment over time
    axes[1, 0].semilogy(times, enrichment, 'g-', linewidth=2)
    axes[1, 0].set_xlabel('Time (hours)')
    axes[1, 0].set_ylabel('Enrichment (target/background)')
    axes[1, 0].set_title('Enrichment Over Time')
    axes[1, 0].grid(True, alpha=0.3)

    # Resource depletion
    dNTP = [h['dNTP'] * 1e6 for h in history]  # µM
    axes[1, 1].plot(times, dNTP, 'purple', linewidth=2)
    axes[1, 1].set_xlabel('Time (hours)')
    axes[1, 1].set_ylabel('dNTP concentration (µM)')
    axes[1, 1].set_title('Resource Depletion')
    axes[1, 1].grid(True, alpha=0.3)

    plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=300)
        logger.info(f"Plot saved to {output_path}")
    else:
        plt.show()


if __name__ == "__main__":
    import sys
    logging.basicConfig(level=logging.INFO)

    print("Stochastic SWGA Simulator")
    print("Uses Gillespie algorithm to model amplification dynamics")
    print("\nRequires:")
    print("  - Primers")
    print("  - AmplificationNetwork for target and background")
    print("\nUsage:")
    print("  from neoswga.core.network_optimizer import AmplificationNetwork")
    print("  from neoswga.core.stochastic_simulator import GillespieSimulator")
    print("  ")
    print("  simulator = GillespieSimulator(primers, fg_network, bg_network)")
    print("  history = simulator.simulate(max_time=3600.0)")
    print("  plot_simulation_results(history)")
