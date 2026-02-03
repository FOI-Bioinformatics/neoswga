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
from typing import Dict, Any, Optional, Tuple
import math


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
        primer_length: int = 10
    ) -> float:
        """
        Calculate total Tm correction from all additives.

        The correction is additive across all components:
        - Uniform corrections (DMSO, trehalose, formamide, etc.)
        - GC-dependent corrections (betaine, TMAC)

        Args:
            gc_content: GC fraction of primer (0-1)
            primer_length: Primer length in bp

        Returns:
            Total Tm correction in degrees Celsius (negative = Tm lowering)

        Example:
            >>> additives = AdditiveConcentrations(dmso_percent=5.0, betaine_m=1.0)
            >>> additives.calculate_tm_correction(gc_content=0.6, primer_length=12)
            -4.2  # DMSO lowers Tm, betaine normalizes GC
        """
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

    def _dmso_correction(self) -> float:
        """
        DMSO Tm correction.

        DMSO lowers Tm by approximately 0.6C per percent.
        Mechanism: Destabilizes AT base pairs, reduces secondary structure.

        Reference: Varadaraj & Skinner (1994) Gene 140:1-5
        """
        return -0.6 * self.dmso_percent

    def _betaine_uniform_correction(self) -> float:
        """
        Betaine uniform Tm correction (GC-independent component).

        Betaine has a small uniform Tm reduction of ~0.5C per M,
        in addition to its GC-dependent effect.

        Reference: Rees et al. (1993) Biochemistry
        """
        return -0.5 * self.betaine_m

    def _trehalose_correction(self) -> float:
        """
        Trehalose Tm correction.

        Trehalose lowers Tm by approximately 5C per M.
        Mechanism: Stabilizes proteins, modifies water structure.

        Reference: Spiess et al. (2004) Biotechniques 36:732-736
        """
        return -5.0 * self.trehalose_m

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

        Ethanol lowers Tm by approximately 0.5C per percent.
        Mechanism: Reduces secondary structure formation.

        Reference: Cheng et al. (1994) PNAS 91:5695-5699
        """
        return -0.5 * self.ethanol_percent

    def _urea_correction(self) -> float:
        """
        Urea Tm correction.

        Urea lowers Tm by approximately 5C per M.
        Mechanism: Denatures GC-rich regions by disrupting stacking.

        Reference: Hutton (1977) NAR 4:3537-3555
        """
        return -5.0 * self.urea_m

    def _tmac_uniform_correction(self) -> float:
        """
        TMAC uniform Tm correction (GC-independent component).

        At low concentrations (0.01-0.1M), TMAC has a small uniform
        Tm reduction of ~1C per 0.1M.

        Reference: Melchior & von Hippel (1973) PNAS 70:298-302
        """
        return -10.0 * self.tmac_m  # -1C per 0.1M

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
