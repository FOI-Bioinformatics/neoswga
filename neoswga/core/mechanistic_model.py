"""
Mechanistic model for additive effects on SWGA.

Four-pathway model that predicts how reaction conditions affect
SWGA amplification efficiency:

1. Tm modification (primer-template stability)
   - How additives change primer melting temperature
   - GC normalization effects

2. Secondary structure accessibility (template melting)
   - How additives melt template secondary structure
   - Temperature and GC effects

3. Enzyme activity (polymerase performance)
   - Processivity, speed, and stability modifiers
   - DMSO inhibition, betaine enhancement

4. Binding kinetics (association/dissociation rates)
   - kon and koff modifiers
   - SSB effects

This model integrates with the network optimizer to weight primers
by their predicted amplification efficiency under specific conditions.

Usage:
    from neoswga.core.mechanistic_model import MechanisticModel
    from neoswga.core.reaction_conditions import ReactionConditions

    conditions = ReactionConditions(temp=42.0, dmso_percent=5.0, betaine_m=1.0)
    model = MechanisticModel(conditions)

    effects = model.calculate_effects(primer='ATCGATCG', template_gc=0.5)
    print(f"Amplification factor: {effects.predicted_amplification_factor:.2f}")
"""

import math
from dataclasses import dataclass
from typing import Dict, Optional, TYPE_CHECKING

from neoswga.core.mechanistic_params import MECHANISTIC_MODEL_PARAMS, get_polymerase_params

if TYPE_CHECKING:
    from neoswga.core.reaction_conditions import ReactionConditions


@dataclass(frozen=True)
class MechanisticEffects:
    """
    Combined effects from mechanistic model calculation.

    This dataclass contains the outputs from all four pathways
    as well as combined metrics for primer evaluation.

    Attributes:
        tm_correction: Total Tm correction from additives (C)
        effective_tm: Effective melting temperature (C)

        accessibility_factor: Template accessibility (0-1)

        processivity_factor: Polymerase processivity modifier (0-1.2)
        speed_factor: Extension rate modifier (0-1.1)
        stability_factor: Enzyme stability modifier (0-1.3)

        kon_factor: Association rate modifier (0.1-3.0)
        koff_factor: Dissociation rate modifier (0.2-2.0)

        effective_binding_rate: Combined binding efficiency (0-1)
        effective_extension_rate: Combined extension efficiency (0-1)
        predicted_amplification_factor: Overall amplification prediction (0-1)
    """
    # Pathway 1: Tm
    tm_correction: float
    effective_tm: float

    # Pathway 2: Accessibility
    accessibility_factor: float

    # Pathway 3: Enzyme
    processivity_factor: float
    speed_factor: float
    stability_factor: float

    # Pathway 4: Kinetics
    kon_factor: float
    koff_factor: float

    # Combined metrics
    effective_binding_rate: float
    effective_extension_rate: float
    predicted_amplification_factor: float


class MechanisticModel:
    """
    Calculate mechanistic effects of reaction conditions on SWGA.

    This class implements a four-pathway model to predict how
    reaction conditions (temperature, additives, polymerase) affect
    primer amplification efficiency.

    The model is initialized with ReactionConditions and can then
    calculate effects for any primer sequence.

    Attributes:
        conditions: The reaction conditions being modeled
        params: Literature-based model parameters
    """

    def __init__(self, conditions: 'ReactionConditions'):
        """
        Initialize mechanistic model with reaction conditions.

        Args:
            conditions: ReactionConditions object with temp, additives, polymerase
        """
        self.conditions = conditions
        self.params = MECHANISTIC_MODEL_PARAMS

        # Get polymerase-specific parameters
        self._poly_params = get_polymerase_params(conditions.polymerase)

        # Pre-calculate condition-dependent enzyme activity (doesn't depend on primer)
        self._enzyme_activity = self._calculate_enzyme_activity()

    def calculate_effects(
        self,
        primer: str,
        template_gc: float
    ) -> MechanisticEffects:
        """
        Calculate all mechanistic effects for a primer.

        Args:
            primer: Primer DNA sequence
            template_gc: GC content of the template genome (0-1)

        Returns:
            MechanisticEffects with all pathway outputs and combined metrics
        """
        # Pathway 1: Tm modification
        primer_gc = self._primer_gc(primer)
        primer_length = len(primer)
        tm_correction = self._calculate_tm_correction(primer_gc, primer_length)
        effective_tm = self._calculate_effective_tm(primer, primer_gc, tm_correction)

        # Pathway 2: Template accessibility
        accessibility = self._calculate_accessibility(template_gc)

        # Pathway 3: Enzyme activity (pre-calculated)
        enzyme = self._enzyme_activity

        # Pathway 4: Binding kinetics
        kinetics = self._calculate_kinetics(effective_tm)

        # Combined binding rate
        base_binding = self._tm_to_binding_prob(effective_tm)
        effective_binding = base_binding * accessibility * kinetics['kon_factor']
        effective_binding = max(0.0, min(1.0, effective_binding))

        # Combined extension rate
        effective_extension = enzyme['speed_factor'] * enzyme['processivity_factor']

        # Overall amplification factor
        # Product of binding efficiency, extension efficiency, and accessibility
        amplification = effective_binding * effective_extension * accessibility
        amplification = max(0.0, min(1.0, amplification))

        return MechanisticEffects(
            tm_correction=tm_correction,
            effective_tm=effective_tm,
            accessibility_factor=accessibility,
            processivity_factor=enzyme['processivity_factor'],
            speed_factor=enzyme['speed_factor'],
            stability_factor=enzyme['stability_factor'],
            kon_factor=kinetics['kon_factor'],
            koff_factor=kinetics['koff_factor'],
            effective_binding_rate=effective_binding,
            effective_extension_rate=effective_extension,
            predicted_amplification_factor=amplification,
        )

    # =========================================================================
    # Pathway 1: Tm Modification
    # =========================================================================

    def _calculate_tm_correction(self, gc_content: float, primer_length: int) -> float:
        """
        Calculate total Tm correction from additives.

        Includes both uniform corrections and GC-dependent normalization.

        Args:
            gc_content: Primer GC content (0-1)
            primer_length: Primer length in bp

        Returns:
            Total Tm correction in degrees Celsius
        """
        p = self.params['tm']
        c = self.conditions

        correction = 0.0

        # Uniform corrections (GC-independent)
        correction += p['dmso_coef'] * c.dmso_percent
        correction += p['formamide_coef'] * c.formamide_percent
        correction += p['trehalose_coef'] * c.trehalose_m
        correction += p['ethanol_coef'] * c.ethanol_percent
        correction += p['betaine_uniform_coef'] * c.betaine_m
        correction += p['urea_coef'] * c.urea_m
        correction += p['tmac_uniform_coef'] * c.tmac_m

        # GC-dependent normalization using sigmoid model
        gc_correction = self._calculate_gc_normalization(gc_content, primer_length)
        correction += gc_correction

        return correction

    def _calculate_gc_normalization(
        self,
        gc_content: float,
        primer_length: int
    ) -> float:
        """
        GC-normalization with sigmoid dose-response.

        Replaces linear model with more accurate sigmoid for betaine and TMAC.
        At saturation, these additives completely equalize AT and GC Tm.

        Args:
            gc_content: Primer GC content (0-1)
            primer_length: Primer length in bp

        Returns:
            GC normalization correction (negative for GC-rich, positive for AT-rich)
        """
        p = self.params['tm']
        c = self.conditions

        gc_deviation = gc_content - 0.5

        # Sigmoid factors (not linear!)
        betaine_factor = self._sigmoid(
            c.betaine_m,
            p['betaine_gc_midpoint'],
            p['betaine_gc_steepness']
        )

        # TMAC has stronger effect per M, scale for practical range
        tmac_factor = self._sigmoid(
            c.tmac_m * 10,  # Scale for 0-0.1M practical range
            p['tmac_gc_midpoint'],
            p['tmac_gc_steepness']
        )

        # Multiplicative combination (each independently contributes)
        equalization = 1.0 - (1.0 - betaine_factor) * (1.0 - tmac_factor)

        # Length-dependent scaling from Wallace rule
        scale = 2.0 * primer_length

        return -scale * gc_deviation * equalization

    def _calculate_effective_tm(
        self,
        primer: str,
        gc_content: float,
        tm_correction: float
    ) -> float:
        """
        Calculate effective Tm for primer under these conditions.

        Uses the ReactionConditions method if available, otherwise
        estimates from sequence.

        Args:
            primer: Primer sequence
            gc_content: Primer GC content
            tm_correction: Pre-calculated Tm correction

        Returns:
            Effective Tm in degrees Celsius
        """
        # Try to use ReactionConditions method
        try:
            return self.conditions.calculate_effective_tm(primer)
        except (AttributeError, Exception):
            # Fall back to simple estimation
            # Wallace rule: Tm = 2*(AT) + 4*(GC)
            length = len(primer)
            gc_count = int(gc_content * length)
            at_count = length - gc_count
            base_tm = 2 * at_count + 4 * gc_count
            return base_tm + tm_correction

    # =========================================================================
    # Pathway 2: Secondary Structure Accessibility
    # =========================================================================

    def _calculate_accessibility(self, template_gc: float) -> float:
        """
        Calculate template accessibility based on secondary structure.

        High GC templates have more secondary structure, reducing accessibility.
        DMSO and betaine melt structure, improving accessibility.

        Args:
            template_gc: Template genome GC content (0-1)

        Returns:
            Accessibility factor (0.1-1.0)
        """
        p = self.params['structure']
        c = self.conditions

        # Base accessibility from GC content
        # GC > 50% reduces accessibility
        if template_gc <= 0.5:
            base = 1.0
        else:
            base = 1.0 - p['gc_structure_coef'] * (template_gc - 0.5)
            base = max(0.3, base)  # Floor at 30%

        # DMSO melts structure (saturating effect)
        dmso_melt = p['dmso_max_effect'] * (
            1.0 - math.exp(-p['dmso_melt_rate'] * c.dmso_percent)
        )

        # Betaine also helps melt structure
        betaine_melt = p['betaine_max_effect'] * (
            1.0 - math.exp(-p['betaine_melt_rate'] * c.betaine_m)
        )

        # Temperature effect (higher temp melts more structure)
        temp_effect = p['temp_structure_coef'] * max(0, c.temp - 25.0)

        # Combined: each independently reduces remaining structure
        structure_remaining = (1 - base)
        structure_remaining *= (1 - dmso_melt)
        structure_remaining *= (1 - betaine_melt)
        structure_remaining *= (1 - min(0.5, temp_effect))  # Cap temp effect

        accessibility = 1.0 - max(0, structure_remaining)

        return max(0.1, min(1.0, accessibility))

    # =========================================================================
    # Pathway 3: Enzyme Activity
    # =========================================================================

    def _calculate_enzyme_activity(self) -> Dict[str, float]:
        """
        Calculate enzyme activity modifiers.

        This includes processivity, speed, and stability factors
        based on additives and temperature.

        Returns:
            Dictionary with processivity_factor, speed_factor, stability_factor
        """
        p = self.params['enzyme']
        c = self.conditions
        poly = self._poly_params

        processivity = 1.0
        speed = 1.0
        stability = 1.0

        # DMSO effect (threshold then steep drop)
        dmso = c.dmso_percent
        if dmso <= poly['dmso_threshold']:
            dmso_effect = poly['dmso_mild_coef'] * dmso
        else:
            dmso_effect = (
                poly['dmso_mild_coef'] * poly['dmso_threshold']
                + poly['dmso_steep_coef'] * (dmso - poly['dmso_threshold'])
            )
        processivity *= (1 - dmso_effect)
        speed *= (1 - 0.5 * dmso_effect)

        # Betaine effect (enhancement then inhibition)
        betaine = c.betaine_m
        if betaine <= p['betaine_peak']:
            # Enhancement phase
            betaine_boost = p['betaine_enhancement'] * betaine / p['betaine_peak']
            stability *= (1 + betaine_boost)
            processivity *= (1 + 0.5 * betaine_boost)
        elif betaine <= p['betaine_inhibition_start']:
            # Plateau phase
            stability *= (1 + p['betaine_enhancement'])
        else:
            # Inhibition phase
            excess = betaine - p['betaine_inhibition_start']
            inhibition = p['betaine_inhibition_coef'] * excess
            processivity *= (1 - inhibition)
            # Stability still enhanced but less
            stability *= (1 + p['betaine_enhancement'] * 0.5)

        # Formamide (always inhibitory)
        formamide_effect = p['formamide_coef'] * c.formamide_percent
        processivity *= (1 - formamide_effect)
        stability *= (1 - 0.5 * formamide_effect)

        # Mg2+ effect
        mg = c.mg_conc
        if mg < p['mg_low_threshold']:
            # Low Mg reduces activity
            mg_factor = mg / p['mg_low_threshold'] if mg > 0 else 0.2
        elif mg > p['mg_high_threshold']:
            # High Mg can inhibit
            mg_factor = 1.0 - 0.1 * (mg - p['mg_high_threshold'])
        else:
            mg_factor = 1.0
        processivity *= max(0.2, mg_factor)

        # Temperature deviation from optimal
        temp_dev = abs(c.temp - poly['optimal_temp'])
        temp_penalty = p['temp_activity_coef'] * temp_dev
        speed *= (1 - min(0.5, temp_penalty))

        # Glycerol
        glycerol = c.glycerol_percent
        stability *= (1 + p['glycerol_stability'] * glycerol)
        speed *= (1 - p['glycerol_speed_penalty'] * glycerol)

        # DMSO-Mg interaction
        interactions = self.params['interactions']
        if dmso > interactions['dmso_mg_threshold'] and mg < 3.0:
            chelation = (
                interactions['dmso_mg_chelation']
                * (dmso - interactions['dmso_mg_threshold'])
                * (3.0 - mg) / 3.0
            )
            processivity *= (1 - chelation)

        return {
            'processivity_factor': max(0.1, min(1.2, processivity)),
            'speed_factor': max(0.3, min(1.1, speed)),
            'stability_factor': max(0.2, min(1.3, stability)),
        }

    # =========================================================================
    # Pathway 4: Binding Kinetics
    # =========================================================================

    def _calculate_kinetics(self, effective_tm: float) -> Dict[str, float]:
        """
        Calculate binding kinetics modifiers.

        Optimal binding occurs when Tm is ~5-10C above reaction temp.
        Additives can boost kon (faster binding) but may increase koff.

        Args:
            effective_tm: Primer effective Tm (C)

        Returns:
            Dictionary with kon_factor and koff_factor
        """
        p = self.params['kinetics']
        c = self.conditions

        # Optimal Tm is above reaction temp
        delta = effective_tm - c.temp
        optimal = p['optimal_delta_t']
        width = p['delta_t_width']

        # kon: Gaussian centered on optimal delta
        kon_tm_factor = math.exp(-0.5 * ((delta - optimal) / width) ** 2)

        # Additive effects on kon
        kon_factor = kon_tm_factor
        kon_factor *= (1 + p['betaine_kon_boost'] * min(c.betaine_m, 2.0))
        kon_factor *= (1 + p['dmso_kon_boost'] * min(c.dmso_percent, 5.0))

        # SSB dramatically increases kon
        if c.ssb:
            kon_factor *= p['ssb_kon_multiplier']

        # koff: depends on Tm vs reaction temp
        if delta > 0:
            # Tm above reaction temp = more stable binding = lower koff
            koff_factor = math.exp(-0.1 * delta)
        else:
            # Tm below reaction temp = less stable = higher koff
            koff_factor = math.exp(-0.2 * delta)  # More sensitive to low Tm

        # Additives can increase koff (destabilizing)
        koff_factor *= (1 + p['betaine_koff_penalty'] * c.betaine_m)
        koff_factor *= (1 + p['dmso_koff_penalty'] * c.dmso_percent)

        return {
            'kon_factor': max(0.1, min(3.0, kon_factor)),
            'koff_factor': max(0.2, min(2.0, koff_factor)),
        }

    # =========================================================================
    # Helper Methods
    # =========================================================================

    def _tm_to_binding_prob(self, effective_tm: float) -> float:
        """
        Convert effective Tm to binding probability.

        Uses Gaussian centered on optimal Tm (reaction temp + optimal_delta_t).

        Args:
            effective_tm: Primer effective Tm (C)

        Returns:
            Binding probability (0-1)
        """
        p = self.params['kinetics']
        optimal = p['optimal_delta_t']
        width = p['delta_t_width']

        delta = effective_tm - self.conditions.temp
        prob = math.exp(-0.5 * ((delta - optimal) / width) ** 2)

        return max(0.0, min(1.0, prob))

    @staticmethod
    def _primer_gc(primer: str) -> float:
        """
        Calculate primer GC content.

        Args:
            primer: DNA sequence

        Returns:
            GC content as fraction (0-1)
        """
        if not primer:
            return 0.5
        gc = primer.upper().count('G') + primer.upper().count('C')
        return gc / len(primer)

    @staticmethod
    def _sigmoid(x: float, midpoint: float, steepness: float) -> float:
        """
        Sigmoidal dose-response function.

        Args:
            x: Input value
            midpoint: x value at sigmoid center (50% response)
            steepness: Sigmoid steepness (higher = steeper)

        Returns:
            Value between 0 and 1
        """
        return 1.0 / (1.0 + math.exp(-steepness * (x - midpoint)))

    def get_enzyme_parameters(self) -> Dict[str, float]:
        """
        Get the pre-calculated enzyme activity parameters.

        Returns:
            Dictionary with processivity_factor, speed_factor, stability_factor
        """
        return self._enzyme_activity.copy()

    def __repr__(self) -> str:
        """String representation."""
        return (
            f"MechanisticModel("
            f"temp={self.conditions.temp}C, "
            f"polymerase={self.conditions.polymerase})"
        )
