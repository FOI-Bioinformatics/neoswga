"""
Additive concentrations for SWGA reaction conditions.

Extracted from ReactionConditions to follow Single Responsibility Principle.
This module handles additive-specific logic:
- Tm corrections from additives
- GC-normalization effects
- Additive interaction modeling

The AdditiveConcentrations class is immutable (frozen dataclass) to prevent
accidental modification during pipeline execution.

Usage:
    # Create with factory methods
    additives = AdditiveConcentrations.for_standard_phi29()
    additives = AdditiveConcentrations.for_enhanced_equiphi29()
    additives = AdditiveConcentrations.for_extreme_gc()

    # Or with explicit values
    additives = AdditiveConcentrations(
        dmso_percent=5.0,
        betaine_m=1.0,
    )

    # Calculate Tm correction
    correction = additives.calculate_tm_correction(gc_content=0.6, primer_length=12)

References:
    - DMSO: Varadaraj & Skinner (1994) Gene 140:1-5
    - Betaine: Henke et al. (1997) NAR 25:3957-3958
    - Trehalose: Spiess et al. (2004) Biotechniques 36:732-736
    - Formamide: Blake & Delcourt (1996) NAR 24:2095-2103
    - TMAC: Melchior & von Hippel (1973) PNAS 70:298-302
"""

from dataclasses import dataclass
from typing import Dict, Any, Optional, Tuple, List
import math


# =============================================================================
# Arrhenius-Based Tm Correction Calculator
# =============================================================================

class ArrheniusTmCorrector:
    """
    Calculate temperature-dependent Tm corrections using Arrhenius kinetics.

    The Tm correction from an additive varies with reaction temperature.
    This class models this using the Arrhenius equation:

        dTm(T) = dTm_ref * exp(-Ea/R * (1/T - 1/T_ref))

    Where:
        dTm_ref: Reference Tm correction at T_ref (typically 37C)
        Ea: Activation energy for the additive's destabilizing effect
        R: Gas constant (8.314 J/mol*K)
        T, T_ref: Temperatures in Kelvin

    Literature basis:
        - Chester & Marshak (1993) measured DMSO effects at multiple temps
        - Rees et al. (1993) measured betaine effects at 25C, 37C, 65C
        - Hutton (1977) measured urea effects at multiple temps

    Example:
        >>> corrector = ArrheniusTmCorrector(42.0)  # EquiPhi29 at 42C
        >>> correction = corrector.calculate_correction(
        ...     'dmso', concentration=5.0, gc_content=0.5
        ... )
        >>> print(f"DMSO correction at 42C: {correction:.1f}C")
        DMSO correction at 42C: -2.9C
    """

    R = 8.314  # Gas constant in J/mol*K

    def __init__(self, reaction_temp_celsius: float):
        """
        Initialize the Arrhenius Tm corrector.

        Args:
            reaction_temp_celsius: Reaction temperature in Celsius
        """
        self.reaction_temp_celsius = reaction_temp_celsius
        self.reaction_temp_k = reaction_temp_celsius + 273.15
        self._params = None  # Lazy-loaded

    @property
    def params(self) -> Dict[str, Dict[str, Any]]:
        """Lazy-load additive parameters."""
        if self._params is None:
            from neoswga.core.mechanistic_params import ADDITIVE_TM_PARAMS
            self._params = ADDITIVE_TM_PARAMS
        return self._params

    def calculate_correction(
        self,
        additive: str,
        concentration: float,
        gc_content: float = 0.5,
        primer_length: int = 10
    ) -> float:
        """
        Calculate Tm correction for an additive at the reaction temperature.

        Args:
            additive: Additive name (dmso, betaine, formamide, etc.)
            concentration: Additive concentration (% for DMSO/formamide/ethanol,
                          M for betaine/trehalose/urea/tmac)
            gc_content: GC fraction of primer (0-1), used for GC-dependent
                       additives like betaine, TMAC, and urea
            primer_length: Primer length in bp, used for GC normalization scaling

        Returns:
            Tm correction in degrees Celsius (negative = Tm lowering)

        Raises:
            ValueError: If additive is not recognized
        """
        additive_lower = additive.lower()
        if additive_lower not in self.params:
            available = list(self.params.keys())
            raise ValueError(
                f"Unknown additive '{additive}'. "
                f"Available: {', '.join(available)}"
            )

        params = self.params[additive_lower]

        # Arrhenius temperature scaling factor
        temp_factor = self._arrhenius_factor(
            params['activation_energy'],
            params['ref_temp']
        )

        # Effective coefficient at this temperature
        effective_coef = params['ref_coef'] * temp_factor

        # Base correction from concentration
        correction = effective_coef * concentration

        # GC-dependent adjustment if applicable
        if params.get('gc_dependent', False):
            gc_adjustment = self._gc_adjustment(
                additive_lower, params, gc_content, concentration, primer_length
            )
            correction += gc_adjustment

        return correction

    def _arrhenius_factor(self, activation_energy: float, ref_temp_k: float) -> float:
        """
        Calculate Arrhenius temperature scaling factor.

        Args:
            activation_energy: Activation energy in J/mol
            ref_temp_k: Reference temperature in Kelvin

        Returns:
            Scaling factor (>1 at higher temps, <1 at lower temps)
        """
        exponent = -activation_energy / self.R * (
            1.0 / self.reaction_temp_k - 1.0 / ref_temp_k
        )
        return math.exp(exponent)

    def _gc_adjustment(
        self,
        additive: str,
        params: Dict[str, Any],
        gc_content: float,
        concentration: float,
        primer_length: int
    ) -> float:
        """
        Calculate GC-dependent adjustment for an additive.

        For betaine and TMAC: GC equalization effect moves Tm toward 50% GC.
        For urea: Preferentially destabilizes GC base pairs.

        Args:
            additive: Additive name
            params: Additive parameters dictionary
            gc_content: GC content fraction (0-1)
            concentration: Additive concentration
            primer_length: Primer length in bp

        Returns:
            Additional Tm correction for GC effect
        """
        gc_deviation = gc_content - 0.5

        if additive in ('betaine', 'tmac'):
            # GC equalization effect
            eq_conc = params.get('gc_equalization_conc', 5.0)
            # Sigmoid equalization factor
            equalization = self._sigmoid(concentration, eq_conc / 3.0, 1.5)

            # Length-dependent scaling (Wallace rule: dTm/d(gc) = 2*length)
            scale = 2.0 * primer_length

            # Correction moves Tm towards 50% GC
            return -scale * gc_deviation * equalization

        elif additive == 'urea':
            # Preferential GC destabilization
            gc_preference = params.get('gc_preference', 1.3)

            # Extra destabilization for GC-rich sequences
            if gc_content > 0.5:
                extra_factor = (gc_preference - 1.0) * (gc_content - 0.5) * 2.0
                # Return additional correction (already negative from base)
                return params['ref_coef'] * concentration * extra_factor

        return 0.0

    @staticmethod
    def _sigmoid(x: float, midpoint: float, steepness: float) -> float:
        """Sigmoidal dose-response function."""
        return 1.0 / (1.0 + math.exp(-steepness * (x - midpoint)))

    def calculate_total_correction(
        self,
        additives: Dict[str, float],
        gc_content: float = 0.5,
        primer_length: int = 10
    ) -> float:
        """
        Calculate total Tm correction from multiple additives.

        Args:
            additives: Dictionary mapping additive names to concentrations
                      e.g., {'dmso': 5.0, 'betaine': 1.0}
            gc_content: GC content fraction (0-1)
            primer_length: Primer length in bp

        Returns:
            Total Tm correction in degrees Celsius
        """
        total = 0.0
        for additive, concentration in additives.items():
            if concentration > 0:
                total += self.calculate_correction(
                    additive, concentration, gc_content, primer_length
                )
        return total

    def get_temperature_sensitivity(self, additive: str) -> float:
        """
        Get temperature sensitivity for an additive.

        Higher activation energy means more temperature-dependent effect.

        Args:
            additive: Additive name

        Returns:
            Activation energy in J/mol (higher = more temp sensitive)
        """
        params = self.params.get(additive.lower())
        if params is None:
            raise ValueError(f"Unknown additive: {additive}")
        return params['activation_energy']

    def compare_temperatures(
        self,
        additive: str,
        concentration: float,
        temp1_celsius: float,
        temp2_celsius: float,
        gc_content: float = 0.5
    ) -> Tuple[float, float, float]:
        """
        Compare Tm correction at two different temperatures.

        Useful for understanding how additive effects change between
        phi29 (30C) and equiphi29 (42C) conditions.

        Args:
            additive: Additive name
            concentration: Additive concentration
            temp1_celsius: First temperature (C)
            temp2_celsius: Second temperature (C)
            gc_content: GC content fraction

        Returns:
            (correction_at_temp1, correction_at_temp2, ratio)
        """
        # Calculate at temp1
        corrector1 = ArrheniusTmCorrector(temp1_celsius)
        corr1 = corrector1.calculate_correction(additive, concentration, gc_content)

        # Calculate at temp2
        corrector2 = ArrheniusTmCorrector(temp2_celsius)
        corr2 = corrector2.calculate_correction(additive, concentration, gc_content)

        # Ratio (avoiding division by zero)
        ratio = corr2 / corr1 if abs(corr1) > 0.01 else 1.0

        return (corr1, corr2, ratio)


@dataclass(frozen=True)
class AdditiveConcentrations:
    """
    Immutable container for reaction additive concentrations.

    All concentrations are validated at construction time.
    Use factory methods for common configurations.

    Attributes:
        dmso_percent: DMSO concentration (0-10%)
        betaine_m: Betaine concentration (0-2.5 M)
        trehalose_m: Trehalose concentration (0-1.0 M)
        formamide_percent: Formamide concentration (0-10%)
        glycerol_percent: Glycerol concentration (0-15%)
        bsa_ug_ml: BSA concentration (0-400 ug/mL)
        peg_percent: PEG concentration (0-15%)
        ethanol_percent: Ethanol concentration (0-5%)
        urea_m: Urea concentration (0-2.0 M)
        tmac_m: TMAC concentration (0-0.1 M)
    """
    dmso_percent: float = 0.0
    betaine_m: float = 0.0
    trehalose_m: float = 0.0
    formamide_percent: float = 0.0
    glycerol_percent: float = 0.0
    bsa_ug_ml: float = 0.0
    peg_percent: float = 0.0
    ethanol_percent: float = 0.0
    urea_m: float = 0.0
    tmac_m: float = 0.0

    def __post_init__(self):
        """Validate all concentrations are within valid ranges."""
        self._validate_range('dmso_percent', self.dmso_percent, 0, 10)
        self._validate_range('betaine_m', self.betaine_m, 0, 2.5)
        self._validate_range('trehalose_m', self.trehalose_m, 0, 1.0)
        self._validate_range('formamide_percent', self.formamide_percent, 0, 10)
        self._validate_range('glycerol_percent', self.glycerol_percent, 0, 15)
        self._validate_range('bsa_ug_ml', self.bsa_ug_ml, 0, 400)
        self._validate_range('peg_percent', self.peg_percent, 0, 15)
        self._validate_range('ethanol_percent', self.ethanol_percent, 0, 5)
        self._validate_range('urea_m', self.urea_m, 0, 2.0)
        self._validate_range('tmac_m', self.tmac_m, 0, 0.1)

    def _validate_range(self, name: str, value: float, min_val: float, max_val: float):
        """Validate a value is within range."""
        if value < min_val or value > max_val:
            raise ValueError(
                f"{name} = {value} is outside valid range [{min_val}, {max_val}]"
            )

    # =========================================================================
    # Tm Correction Calculations
    # =========================================================================

    def calculate_tm_correction(
        self,
        gc_content: float = 0.5,
        primer_length: int = 10,
        reaction_temp_celsius: Optional[float] = None
    ) -> float:
        """
        Calculate total Tm correction from all additives.

        The correction is additive across all components:
        - Uniform corrections (DMSO, trehalose, formamide, etc.)
        - GC-dependent corrections (betaine, TMAC)

        When reaction_temp_celsius is provided, uses Arrhenius-based
        temperature-dependent corrections for more accurate predictions.

        Args:
            gc_content: GC fraction of primer (0-1)
            primer_length: Primer length in bp
            reaction_temp_celsius: Optional reaction temperature for
                Arrhenius-based corrections. If None, uses fixed coefficients
                (backward compatible).

        Returns:
            Total Tm correction in degrees Celsius (negative = Tm lowering)

        Example:
            >>> additives = AdditiveConcentrations(dmso_percent=5.0, betaine_m=1.0)
            >>> # Fixed coefficients (backward compatible)
            >>> additives.calculate_tm_correction(gc_content=0.6, primer_length=12)
            -4.8
            >>> # Arrhenius-based at 42C (equiphi29)
            >>> additives.calculate_tm_correction(gc_content=0.6, primer_length=12,
            ...                                    reaction_temp_celsius=42.0)
            -5.1
        """
        if reaction_temp_celsius is not None:
            return self._calculate_tm_correction_arrhenius(
                gc_content, primer_length, reaction_temp_celsius
            )

        # Legacy fixed-coefficient calculation
        correction = 0.0

        # Uniform corrections (GC-independent)
        correction += self._dmso_correction()
        correction += self._trehalose_correction()
        correction += self._formamide_correction()
        correction += self._ethanol_correction()
        correction += self._urea_correction()
        correction += self._tmac_uniform_correction()
        correction += self._betaine_uniform_correction()

        # GC-dependent corrections
        correction += self._gc_normalization_correction(gc_content, primer_length)

        return correction

    def _calculate_tm_correction_arrhenius(
        self,
        gc_content: float,
        primer_length: int,
        reaction_temp_celsius: float
    ) -> float:
        """
        Calculate Tm correction using Arrhenius-based temperature dependence.

        Uses literature-based activation energies to model how additive
        effects change with temperature.

        Args:
            gc_content: GC fraction (0-1)
            primer_length: Primer length in bp
            reaction_temp_celsius: Reaction temperature in Celsius

        Returns:
            Total Tm correction in degrees Celsius
        """
        corrector = ArrheniusTmCorrector(reaction_temp_celsius)

        # Build additives dictionary
        additives_dict = {}
        if self.dmso_percent > 0:
            additives_dict['dmso'] = self.dmso_percent
        if self.betaine_m > 0:
            additives_dict['betaine'] = self.betaine_m
        if self.trehalose_m > 0:
            additives_dict['trehalose'] = self.trehalose_m
        if self.formamide_percent > 0:
            additives_dict['formamide'] = self.formamide_percent
        if self.ethanol_percent > 0:
            additives_dict['ethanol'] = self.ethanol_percent
        if self.urea_m > 0:
            additives_dict['urea'] = self.urea_m
        if self.tmac_m > 0:
            additives_dict['tmac'] = self.tmac_m

        return corrector.calculate_total_correction(
            additives_dict, gc_content, primer_length
        )

    def _dmso_correction(self) -> float:
        """
        DMSO Tm correction.

        DMSO lowers Tm by approximately 0.55C per percent.
        Mechanism: Destabilizes AT base pairs, reduces secondary structure.

        Literature: Chester & Marshak (1993) measured -0.5 to -0.6C/%
        at 37C. We use -0.55C/% as the reference value.

        Reference: Varadaraj & Skinner (1994) Gene; Chester & Marshak (1993)
        """
        return -0.55 * self.dmso_percent

    def _betaine_uniform_correction(self) -> float:
        """
        Betaine uniform Tm correction (GC-independent component).

        Betaine has a uniform Tm reduction of ~1.2C per M,
        in addition to its GC-dependent effect.

        Literature: Rees et al. (1993) measured -1.0 to -1.5C/M
        at 37C. We use -1.2C/M as the reference value.

        Reference: Rees et al. (1993) Biochemistry; Henke et al. (1997) NAR
        """
        return -1.2 * self.betaine_m

    def _trehalose_correction(self) -> float:
        """
        Trehalose Tm correction.

        Trehalose lowers Tm by approximately 3C per M.
        Mechanism: Stabilizes proteins, modifies water structure.

        Literature: Spiess et al. (2004) measured ~2-4C/M depending
        on conditions. We use -3.0C/M as a recalibrated value.

        Reference: Spiess et al. (2004) Biotechniques 36:732-736
        """
        return -3.0 * self.trehalose_m

    def _formamide_correction(self) -> float:
        """
        Formamide Tm correction.

        Formamide lowers Tm by approximately 0.65C per percent.
        Mechanism: Destabilizes hydrogen bonding.

        Reference: Blake & Delcourt (1996) NAR 24:2095-2103
        """
        return -0.65 * self.formamide_percent

    def _ethanol_correction(self) -> float:
        """
        Ethanol Tm correction.

        Ethanol lowers Tm by approximately 0.4C per percent.
        Mechanism: Reduces secondary structure formation.

        Literature: Cheng et al. (1994) measured ~0.3-0.5C/%.
        We use -0.4C/% as the recalibrated value.

        Reference: Cheng et al. (1994) PNAS 91:5695-5699
        """
        return -0.4 * self.ethanol_percent

    def _urea_correction(self) -> float:
        """
        Urea Tm correction.

        Urea lowers Tm by approximately 2.5C per M.
        Mechanism: Denatures GC-rich regions by disrupting stacking.

        Literature: Lesnick & Bhalla (1995) measured -2.0 to -3.0C/M.
        Previous value of -5.0C/M was too high.

        Reference: Lesnick & Bhalla (1995) NAR; Hutton (1977) NAR
        """
        return -2.5 * self.urea_m

    def _tmac_uniform_correction(self) -> float:
        """
        TMAC uniform Tm correction (GC-independent component).

        At low concentrations (0.01-0.1M), TMAC has a small uniform
        Tm reduction of ~0.5C per M. The primary effect of TMAC is
        GC-dependent (handled separately).

        Literature: Melchior & von Hippel (1973) showed that TMAC's
        main effect is GC equalization, with minimal uniform reduction.

        Reference: Melchior & von Hippel (1973) PNAS 70:298-302
        """
        return -0.5 * self.tmac_m

    def _gc_normalization_correction(
        self,
        gc_content: float,
        primer_length: int
    ) -> float:
        """
        Calculate GC-normalization effect from TMAC and betaine.

        Uses sigmoid dose-response model for improved accuracy over linear.
        TMAC and betaine equalize AT and GC contributions to Tm,
        moving the Tm towards the 50% GC baseline.

        The sigmoid model better captures the saturation behavior observed
        experimentally: effect increases rapidly at low concentrations then
        plateaus as binding sites become saturated.

        Args:
            gc_content: GC fraction (0-1)
            primer_length: Primer length in bp

        Returns:
            Tm correction (negative for GC-rich, positive for AT-rich)

        References:
            - Rees et al. (1993): Betaine full equalization at 5.2M
            - Melchior & von Hippel (1973): TMAC full equalization at 3M
        """
        # Import sigmoid parameters from mechanistic model
        try:
            from neoswga.core.mechanistic_params import MECHANISTIC_MODEL_PARAMS
            params = MECHANISTIC_MODEL_PARAMS['tm']
            betaine_midpoint = params['betaine_gc_midpoint']
            betaine_steepness = params['betaine_gc_steepness']
            tmac_midpoint = params['tmac_gc_midpoint']
            tmac_steepness = params['tmac_gc_steepness']
        except ImportError:
            # Fallback to default values if mechanistic_params not available
            betaine_midpoint = 1.5
            betaine_steepness = 1.5
            tmac_midpoint = 1.0
            tmac_steepness = 2.0

        # Deviation from balanced GC
        gc_deviation = gc_content - 0.5

        # Sigmoid equalization factors (not linear!)
        # Sigmoid: 1 / (1 + exp(-steepness * (x - midpoint)))
        betaine_factor = self._sigmoid(
            self.betaine_m, betaine_midpoint, betaine_steepness
        )

        # TMAC has stronger effect per M, scale for practical 0-0.1M range
        tmac_factor = self._sigmoid(
            self.tmac_m * 10, tmac_midpoint, tmac_steepness
        )

        # Combined equalization (multiplicative model)
        # Each additive independently contributes to GC equalization
        equalization_factor = 1.0 - (1.0 - betaine_factor) * (1.0 - tmac_factor)

        # Length-dependent scaling from Wallace rule
        # dTm/d(gc_fraction) = 2 * length
        scale_factor = 2.0 * primer_length

        # Correction moves Tm towards balanced
        # GC-rich: negative correction (reduce Tm)
        # AT-rich: positive correction (increase Tm)
        return -scale_factor * gc_deviation * equalization_factor

    @staticmethod
    def _sigmoid(x: float, midpoint: float, steepness: float) -> float:
        """
        Sigmoidal dose-response function.

        Args:
            x: Input concentration value
            midpoint: x value at sigmoid center (50% response)
            steepness: Sigmoid steepness (higher = sharper transition)

        Returns:
            Value between 0 and 1
        """
        return 1.0 / (1.0 + math.exp(-steepness * (x - midpoint)))

    # =========================================================================
    # Primer Length Recommendations
    # =========================================================================

    def max_supported_primer_length(self, polymerase: str = 'phi29') -> int:
        """
        Calculate maximum supported primer length with these additives.

        Additives like betaine and DMSO enable longer primers by:
        - Reducing GC-dependence of Tm
        - Melting secondary structure
        - Improving primer annealing specificity

        Args:
            polymerase: Polymerase type ('phi29', 'equiphi29', 'bst', 'klenow')

        Returns:
            Maximum recommended primer length (6-18)
        """
        base_max = 12  # Default swga2 maximum

        # Betaine enables longer primers
        if self.betaine_m >= 0.5:
            base_max += 2
        if self.betaine_m >= 1.0:
            base_max += 1

        # DMSO helps with secondary structure
        if self.dmso_percent >= 3.0:
            base_max += 1
        if self.dmso_percent >= 5.0:
            base_max += 1

        # EquiPhi29 at higher temp can handle longer primers
        if polymerase == 'equiphi29':
            base_max += 1

        # Apply caps based on additive support level
        supported_max = base_max

        # 16bp+: Need moderate additive support
        if supported_max >= 16:
            if not (self.betaine_m >= 1.0 or self.dmso_percent >= 5.0):
                supported_max = 15

        # 17bp+: Need strong combined support
        if supported_max >= 17:
            if not (self.betaine_m >= 1.5 and self.dmso_percent >= 3.0):
                supported_max = 16

        # 18bp+: Need extreme combined support
        if supported_max >= 18:
            if not (self.betaine_m >= 2.0 and self.dmso_percent >= 5.0):
                supported_max = 17

        return min(supported_max, 18)

    def gc_content_range(self) -> Tuple[float, float]:
        """
        Recommended GC content range with these additives.

        Betaine widens the acceptable GC range by equalizing
        AT and GC contributions to Tm.

        Returns:
            (min_gc, max_gc) as fractions (0-1)
        """
        if self.betaine_m >= 1.0:
            return (0.30, 0.70)
        elif self.betaine_m >= 0.5:
            return (0.35, 0.65)
        else:
            return (0.375, 0.625)

    # =========================================================================
    # Factory Methods
    # =========================================================================

    @classmethod
    def none(cls) -> 'AdditiveConcentrations':
        """No additives (all zeros)."""
        return cls()

    @classmethod
    def for_standard_phi29(cls) -> 'AdditiveConcentrations':
        """
        Standard phi29 conditions with no additives.

        Use for: Standard SWGA at 30C with 6-12bp primers.
        """
        return cls()

    @classmethod
    def for_enhanced_equiphi29(cls) -> 'AdditiveConcentrations':
        """
        Enhanced conditions for EquiPhi29 with longer primers.

        Use for: EquiPhi29 at 42C with 12-15bp primers.
        Enables: Higher specificity, reduced off-target.

        Based on: Henke et al. (1997) - betaine + DMSO combination
        """
        return cls(
            dmso_percent=5.0,
            betaine_m=1.0,
        )

    @classmethod
    def for_high_gc(cls) -> 'AdditiveConcentrations':
        """
        High-GC genome conditions.

        Use for: Genomes with >60% GC content.
        Enables: Better coverage of GC-rich regions.

        Based on: High betaine to equalize AT/GC stability
        """
        return cls(
            dmso_percent=5.0,
            betaine_m=2.0,
        )

    @classmethod
    def for_extreme_gc(cls) -> 'AdditiveConcentrations':
        """
        Extreme GC genome conditions (>70% or <30% GC).

        Use for: Mycobacterium, Streptomyces, or AT-rich Plasmodium.
        Enables: Full GC-normalization with TMAC.

        Based on: Melchior & von Hippel (1973) - TMAC isostabilization
        """
        return cls(
            dmso_percent=5.0,
            betaine_m=2.0,
            urea_m=0.5,
            tmac_m=0.05,
        )

    @classmethod
    def for_long_primers(cls) -> 'AdditiveConcentrations':
        """
        Conditions for maximum primer length (16-18bp).

        Use for: Maximum specificity applications.
        Enables: Up to 18bp primers with good annealing.

        Based on: Musso et al. (2006) - extreme additive conditions
        """
        return cls(
            dmso_percent=7.0,
            betaine_m=2.0,
            trehalose_m=0.3,
        )

    @classmethod
    def for_crude_samples(cls) -> 'AdditiveConcentrations':
        """
        Conditions for crude or inhibitor-containing samples.

        Use for: Blood, soil, crude lysates.
        Enables: Tolerance to PCR inhibitors.
        """
        return cls(
            betaine_m=1.0,
            glycerol_percent=10.0,
            bsa_ug_ml=400.0,
            peg_percent=5.0,
        )

    @classmethod
    def q_solution_equivalent(cls) -> 'AdditiveConcentrations':
        """
        Qiagen Q-Solution equivalent conditions.

        Approximates the proprietary Q-Solution formulation.
        """
        return cls(
            betaine_m=1.5,
            glycerol_percent=10.0,
            bsa_ug_ml=200.0,
        )

    # =========================================================================
    # Serialization
    # =========================================================================

    def to_dict(self) -> Dict[str, float]:
        """Convert to dictionary for JSON serialization."""
        return {
            'dmso_percent': self.dmso_percent,
            'betaine_m': self.betaine_m,
            'trehalose_m': self.trehalose_m,
            'formamide_percent': self.formamide_percent,
            'glycerol_percent': self.glycerol_percent,
            'bsa_ug_ml': self.bsa_ug_ml,
            'peg_percent': self.peg_percent,
            'ethanol_percent': self.ethanol_percent,
            'urea_m': self.urea_m,
            'tmac_m': self.tmac_m,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, float]) -> 'AdditiveConcentrations':
        """Create from dictionary."""
        return cls(**{k: v for k, v in data.items() if k in cls.__dataclass_fields__})

    def __repr__(self) -> str:
        """Compact string representation showing non-zero values."""
        parts = []
        if self.dmso_percent > 0:
            parts.append(f"DMSO={self.dmso_percent}%")
        if self.betaine_m > 0:
            parts.append(f"betaine={self.betaine_m}M")
        if self.trehalose_m > 0:
            parts.append(f"trehalose={self.trehalose_m}M")
        if self.formamide_percent > 0:
            parts.append(f"formamide={self.formamide_percent}%")
        if self.glycerol_percent > 0:
            parts.append(f"glycerol={self.glycerol_percent}%")
        if self.bsa_ug_ml > 0:
            parts.append(f"BSA={self.bsa_ug_ml}ug/mL")
        if self.peg_percent > 0:
            parts.append(f"PEG={self.peg_percent}%")
        if self.ethanol_percent > 0:
            parts.append(f"ethanol={self.ethanol_percent}%")
        if self.urea_m > 0:
            parts.append(f"urea={self.urea_m}M")
        if self.tmac_m > 0:
            parts.append(f"TMAC={self.tmac_m}M")

        if not parts:
            return "AdditiveConcentrations(none)"
        return f"AdditiveConcentrations({', '.join(parts)})"


# =============================================================================
# Additive Interaction Modeling
# =============================================================================

def estimate_combined_effect(
    additives: AdditiveConcentrations,
    gc_content: float,
    primer_length: int
) -> Dict[str, Any]:
    """
    Estimate combined effects of additives on reaction.

    Provides a summary of how additives affect the reaction,
    including synergistic and antagonistic effects.

    Args:
        additives: Additive concentrations
        gc_content: Primer GC content (0-1)
        primer_length: Primer length in bp

    Returns:
        Dict with effect estimates
    """
    tm_correction = additives.calculate_tm_correction(gc_content, primer_length)
    max_length = additives.max_supported_primer_length()
    gc_range = additives.gc_content_range()

    # Estimate polymerase activity retention
    # Some additives (DMSO, betaine at low conc) are well-tolerated
    # Others (high DMSO, formamide) can inhibit
    activity_factor = 1.0
    if additives.dmso_percent > 7:
        activity_factor *= 0.9
    if additives.formamide_percent > 5:
        activity_factor *= 0.85
    if additives.urea_m > 1.0:
        activity_factor *= 0.8
    if additives.betaine_m > 0 and additives.betaine_m <= 1.5:
        activity_factor *= 1.1  # Slight enhancement

    return {
        'tm_correction': tm_correction,
        'max_primer_length': max_length,
        'gc_range': gc_range,
        'polymerase_activity_factor': activity_factor,
        'secondary_structure_reduction': min(1.0, additives.dmso_percent / 10.0),
        'gc_normalization': min(1.0, additives.betaine_m / 2.0 + additives.tmac_m * 10),
    }
