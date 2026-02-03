"""
Reaction condition modeling for selective whole genome amplification.

Handles the effects of thermodynamic additives and buffer conditions on:
- Melting temperature
- Primer binding
- Polymerase activity
- Optimal primer length selection

Supports additives:
- DMSO (dimethyl sulfoxide)
- Betaine (trimethylglycine)
- Trehalose (disaccharide)
- Formamide
- Glycerol (enzyme stabilizer)
- BSA (bovine serum albumin, inhibitor neutralizer)
- PEG (polyethylene glycol, molecular crowding agent)
- SSB (single-stranded binding proteins)
- Ethanol (secondary structure reducer)
- Urea (GC denaturant)
- TMAC (tetramethylammonium chloride, AT/GC equalizer)

Supports polymerases:
- Phi29 (30C, high processivity ~70kb)
- EquiPhi29 (42-45C, thermostable variant)
- Bst 2.0/3.0 (60-65C, LAMP-compatible)
- Klenow exo- (37C, moderate processivity)

References:
- Henke et al. (1997) Nucleic Acids Res 25:3957-3958 (betaine)
- Varadaraj & Skinner (1994) Gene 140:1-5 (DMSO)
- Spiess et al. (2004) Biotechniques 36:732-736 (trehalose)
- Rahman et al. (2014) PLoS One 9:e112515 (Mg2+ optimization)
- Owczarzy et al. (2008) Biochemistry 47:5336-5353 (salt corrections)
- Melchior & von Hippel (1973) PNAS 70:298-302 (TMAC AT/GC equalization)
- Cheng et al. (1994) PNAS 91:5695-5699 (ethanol effects)
- Notomi et al. (2000) Nucleic Acids Res 28:e63 (Bst polymerase)
"""

import numpy as np
from typing import Dict, Tuple, Optional
import neoswga.core.thermodynamics as thermo
from neoswga.core.additives import AdditiveConcentrations


# ========================================
# Polymerase Characteristics Database
# ========================================
# Comprehensive database of polymerase properties for SWGA applications.
# All values are based on published literature.

POLYMERASE_CHARACTERISTICS = {
    'phi29': {
        'name': 'Phi29 DNA Polymerase',
        'temp_range': (30.0, 40.0),
        'optimal_temp': 30.0,
        'processivity': 70000,  # ~70kb, Blanco et al. (1989)
        'strand_displacement': True,
        'exonuclease': '3to5',  # Proofreading
        'error_rate': 1e-6,     # Very high fidelity
        'requires_primer': True,
        'description': 'High-fidelity, high-processivity polymerase for MDA/SWGA'
    },
    'equiphi29': {
        'name': 'EquiPhi29 DNA Polymerase',
        'temp_range': (42.0, 45.0),
        'optimal_temp': 42.0,
        'processivity': 80000,  # ~80kb at elevated temp
        'strand_displacement': True,
        'exonuclease': '3to5',
        'error_rate': 1e-6,
        'requires_primer': True,
        'description': 'Thermostable phi29 variant for higher specificity'
    },
    'bst': {
        'name': 'Bst 2.0/3.0 DNA Polymerase',
        'temp_range': (60.0, 65.0),
        'optimal_temp': 63.0,
        'processivity': 2000,   # ~1-2kb, Notomi et al. (2000)
        'strand_displacement': True,
        'exonuclease': 'none',  # Large fragment lacks exo
        'error_rate': 1e-4,     # Lower fidelity
        'requires_primer': True,
        'description': 'LAMP-compatible, high-temperature isothermal amplification'
    },
    'klenow': {
        'name': 'Klenow Fragment (exo-)',
        'temp_range': (25.0, 40.0),
        'optimal_temp': 37.0,
        'processivity': 10000,  # ~10kb, Bambara et al. (1978)
        'strand_displacement': True,  # Moderate
        'exonuclease': 'none',  # exo- variant
        'error_rate': 1e-4,
        'requires_primer': True,
        'description': 'Budget-friendly alternative with moderate processivity'
    }
}


def get_polymerase_processivity(polymerase: str) -> int:
    """
    Get processivity for a polymerase by name.

    Args:
        polymerase: Polymerase name ('phi29', 'equiphi29', 'bst', 'klenow')

    Returns:
        Maximum extension length in base pairs

    Example:
        >>> get_polymerase_processivity('phi29')
        70000
        >>> get_polymerase_processivity('bst')
        2000
    """
    poly = polymerase.lower()
    if poly not in POLYMERASE_CHARACTERISTICS:
        raise ValueError(f"Unknown polymerase: {polymerase}")
    return POLYMERASE_CHARACTERISTICS[poly]['processivity']


def list_polymerases() -> Dict[str, str]:
    """
    List all supported polymerases with descriptions.

    Returns:
        Dictionary mapping polymerase name to description
    """
    return {
        name: info['description']
        for name, info in POLYMERASE_CHARACTERISTICS.items()
    }


class ReactionConditions:
    """
    Comprehensive reaction condition model for SWGA.

    Attributes:
        temp: Reaction temperature in degrees Celsius
        dmso_percent: DMSO concentration as percentage (0-10)
        betaine_m: Betaine concentration in M (0-2.5)
        trehalose_m: Trehalose concentration in M (0-1.0)
        formamide_percent: Formamide concentration as percentage (0-10)
        glycerol_percent: Glycerol concentration as percentage (0-15)
        bsa_ug_ml: BSA concentration in micrograms/mL (0-400)
        peg_percent: PEG concentration as percentage (0-15)
        ethanol_percent: Ethanol concentration as percentage (0-5)
        urea_m: Urea concentration in M (0-2.0)
        tmac_m: TMAC concentration in M (0-0.1)
        na_conc: Na+ concentration in mM
        mg_conc: Mg2+ concentration in mM
        ssb: Whether SSB protein is present
        polymerase: Polymerase type ('phi29', 'equiphi29', 'bst', 'klenow')
    """

    def __init__(self,
                 temp: float = 30.0,
                 dmso_percent: float = 0.0,
                 betaine_m: float = 0.0,
                 trehalose_m: float = 0.0,
                 formamide_percent: float = 0.0,
                 glycerol_percent: float = 0.0,
                 bsa_ug_ml: float = 0.0,
                 peg_percent: float = 0.0,
                 ethanol_percent: float = 0.0,
                 urea_m: float = 0.0,
                 tmac_m: float = 0.0,
                 na_conc: float = 50.0,
                 mg_conc: float = 0.0,
                 ssb: bool = False,
                 polymerase: str = 'phi29'):
        """
        Initialize reaction conditions.

        Args:
            temp: Reaction temperature (C)
            dmso_percent: DMSO concentration (%)
            betaine_m: Betaine concentration (M)
            trehalose_m: Trehalose concentration (M)
            formamide_percent: Formamide concentration (%)
            glycerol_percent: Glycerol concentration (%)
            bsa_ug_ml: BSA concentration (ug/mL)
            peg_percent: PEG concentration (%)
            ethanol_percent: Ethanol concentration (%), reduces secondary structure
            urea_m: Urea concentration (M), denatures GC-rich regions
            tmac_m: TMAC concentration (M), equalizes AT/GC Tm
            na_conc: Na+ concentration (mM)
            mg_conc: Mg2+ concentration (mM)
            ssb: SSB protein present
            polymerase: 'phi29', 'equiphi29', 'bst', or 'klenow'
        """
        self.temp = temp
        self.dmso_percent = dmso_percent
        self.betaine_m = betaine_m
        self.trehalose_m = trehalose_m
        self.formamide_percent = formamide_percent
        self.glycerol_percent = glycerol_percent
        self.bsa_ug_ml = bsa_ug_ml
        self.peg_percent = peg_percent
        self.ethanol_percent = ethanol_percent
        self.urea_m = urea_m
        self.tmac_m = tmac_m
        self.na_conc = na_conc
        self.mg_conc = mg_conc
        self.ssb = ssb
        self.polymerase = polymerase.lower()

        # Validate inputs
        self._validate()

    def _validate(self):
        """Validate reaction condition parameters."""
        # Temperature range depends on polymerase
        polymerase_temp_ranges = {
            'phi29': (20, 40),
            'equiphi29': (30, 50),
            'bst': (50, 72),
            'klenow': (25, 40)
        }

        if self.polymerase not in polymerase_temp_ranges:
            raise ValueError(
                f"Unknown polymerase '{self.polymerase}'. "
                f"Use: {', '.join(polymerase_temp_ranges.keys())}"
            )

        min_temp, max_temp = polymerase_temp_ranges[self.polymerase]
        if self.temp < min_temp or self.temp > max_temp:
            raise ValueError(
                f"Temperature {self.temp}C outside valid range for {self.polymerase} "
                f"({min_temp}-{max_temp}C)"
            )

        if self.dmso_percent < 0 or self.dmso_percent > 10:
            raise ValueError(f"DMSO {self.dmso_percent}% outside valid range (0-10%)")

        if self.betaine_m < 0 or self.betaine_m > 2.5:
            raise ValueError(f"Betaine {self.betaine_m}M outside valid range (0-2.5M)")

        if self.trehalose_m < 0 or self.trehalose_m > 1.0:
            raise ValueError(f"Trehalose {self.trehalose_m}M outside valid range (0-1.0M)")

        if self.glycerol_percent < 0 or self.glycerol_percent > 15:
            raise ValueError(f"Glycerol {self.glycerol_percent}% outside valid range (0-15%)")

        if self.bsa_ug_ml < 0 or self.bsa_ug_ml > 400:
            raise ValueError(f"BSA {self.bsa_ug_ml} ug/mL outside valid range (0-400 ug/mL)")

        if self.peg_percent < 0 or self.peg_percent > 15:
            raise ValueError(f"PEG {self.peg_percent}% outside valid range (0-15%)")

        # Validate new additives
        if self.ethanol_percent < 0 or self.ethanol_percent > 5:
            raise ValueError(f"Ethanol {self.ethanol_percent}% outside valid range (0-5%)")

        if self.urea_m < 0 or self.urea_m > 2.0:
            raise ValueError(f"Urea {self.urea_m}M outside valid range (0-2.0M)")

        if self.tmac_m < 0 or self.tmac_m > 0.1:
            raise ValueError(f"TMAC {self.tmac_m}M outside valid range (0-0.1M)")

        if self.formamide_percent < 0 or self.formamide_percent > 10:
            raise ValueError(f"Formamide {self.formamide_percent}% outside valid range (0-10%)")

    def calculate_tm_correction(self, gc_content: float = 0.5,
                                 primer_length: int = 10) -> float:
        """
        Calculate total Tm correction from all additives.

        Args:
            gc_content: GC content of the sequence as fraction (0-1).
                        Used for GC-dependent corrections (TMAC, betaine).
                        Default 0.5 (balanced) for backwards compatibility.
            primer_length: Length of the primer in base pairs.
                          Used for length-dependent GC correction scaling.
                          Default 10 for backwards compatibility.

        Returns:
            Total Tm shift in degrees Celsius (negative = Tm lowering)

        References:
            - DMSO: Varadaraj & Skinner (1994) Gene 140:1-5
            - Betaine: Rees et al. (1993) Biochemistry; Henke et al. (1997) NAR 25:3957
            - Trehalose: Spiess et al. (2004) Biotechniques 36:732-736
            - Formamide: Blake & Delcourt (1996) NAR 24:2095-2103
            - Ethanol: Cheng et al. (1994) PNAS 91:5695-5699
            - Urea: Hutton (1977) NAR 4:3537-3555
            - TMAC: Melchior & von Hippel (1973) PNAS 70:298-302
        """
        correction = 0.0

        # DMSO: -0.6C per percent
        # Empirical data: 10% DMSO lowers Tm by ~6C
        correction -= 0.6 * self.dmso_percent

        # Betaine: GC-dependent effect
        # At 5.2M betaine: AT and GC base pairs have equal stability (Rees 1993)
        # Also has a small uniform Tm reduction (-0.5C per M)
        # The GC-dependent effect normalizes Tm towards 50% GC behavior
        correction -= 0.5 * self.betaine_m  # Uniform component

        # Trehalose: -5C per molar
        # Empirical: 0.4M trehalose lowers Tm by ~2C
        correction -= 5.0 * self.trehalose_m

        # Formamide: -0.65C per percent
        # Empirical: 10% formamide lowers Tm by ~6.5C
        correction -= 0.65 * self.formamide_percent

        # Ethanol: -0.5C per percent
        # Reduces secondary structure formation
        # Empirical: 5% ethanol lowers Tm by ~2.5C
        correction -= 0.5 * self.ethanol_percent

        # Urea: -0.5C per 0.1M (or -5C per M)
        # Denatures GC-rich regions, useful for extreme GC genomes
        correction -= 5.0 * self.urea_m

        # TMAC: small uniform reduction at low concentrations
        # At 0.01-0.1M: ~-1C per 0.1M (GC-dependent effect handled separately)
        correction -= 1.0 * self.tmac_m * 10  # -1C per 0.1M

        # GC-dependent corrections for TMAC and betaine
        # These additives equalize AT and GC contributions to Tm
        gc_correction = self._calculate_gc_normalization(gc_content, primer_length)
        correction += gc_correction

        return correction

    def _calculate_gc_normalization(self, gc_content: float,
                                     primer_length: int = 10) -> float:
        """
        Calculate the GC-normalization effect of TMAC and betaine.

        The Wallace rule gives Tm = 2*(AT) + 4*(GC) per base.
        At 50% GC, this gives Tm = 3*length.
        GC-rich sequences have Tm > 3*length, AT-rich sequences have Tm < 3*length.

        TMAC and betaine equalize AT and GC contributions, moving Tm towards
        the 50% GC baseline. At full isostabilizing concentration, Tm becomes
        independent of GC content.

        Args:
            gc_content: GC content as fraction (0-1)
            primer_length: Length of the primer in base pairs (default 10)

        Returns:
            Tm correction to apply (negative for GC-rich, positive for AT-rich)

        Scientific basis:
            - Melchior & von Hippel (1973): 3M TMAC gives full GC-independence
            - Rees et al. (1993): 5.2M betaine gives full GC-independence
            - At lower concentrations: partial equalization proportional to conc.
            - Wallace rule: dTm/d(gc_fraction) = 2*length (length-dependent)
        """
        # Calculate deviation from balanced (50%) GC content
        # GC-rich (>50%): positive deviation, normally higher Tm
        # AT-rich (<50%): negative deviation, normally lower Tm
        gc_deviation = gc_content - 0.5

        # Equalization factors (fraction of full effect, capped at 1.0)
        # TMAC: full equalization at 3M (Melchior & von Hippel 1973)
        tmac_factor = min(1.0, self.tmac_m / 3.0)

        # Betaine: full equalization at 5.2M (Rees et al. 1993)
        betaine_factor = min(1.0, self.betaine_m / 5.2)

        # Combined equalization using multiplicative model
        # Each additive independently contributes to GC equalization.
        # Combined effect: 1 - (1 - tmac_factor) * (1 - betaine_factor)
        # Example: 50% TMAC effect + 50% betaine effect = 75% combined (not 50%)
        # This prevents double-counting while allowing synergy.
        equalization_factor = 1.0 - (1.0 - tmac_factor) * (1.0 - betaine_factor)

        # The GC-dependent Tm difference from Wallace rule:
        # For a primer of length L with GC fraction f:
        #   Tm = 2*(AT pairs) + 4*(GC pairs)
        #      = 2*L*(1-f) + 4*L*f
        #      = 2L + 2Lf
        #   dTm/df = 2L
        #
        # So the GC-dependent component scales linearly with primer length.
        # For a 10bp primer: scale = 20
        # For a 6bp primer: scale = 12
        # For an 18bp primer: scale = 36
        #
        # This length-dependent scaling is important for accurate corrections
        # across the range of SWGA primer lengths (6-18bp).
        scale_factor = 2.0 * primer_length

        # We apply a correction that moves Tm towards balanced
        # GC-rich sequences get negative correction (reduce Tm)
        # AT-rich sequences get positive correction (increase Tm)
        gc_correction = -scale_factor * gc_deviation * equalization_factor

        return gc_correction

    def adjust_tm(self, tm_base: float, gc_content: float = 0.5,
                  primer_length: int = 10) -> float:
        """
        Apply all corrections to base Tm.

        Args:
            tm_base: Base melting temperature (with salt corrections)
            gc_content: GC content of sequence (0-1) for GC-dependent corrections.
                        Default 0.5 for backwards compatibility.
            primer_length: Length of the primer in base pairs.
                          Default 10 for backwards compatibility.

        Returns:
            Effective Tm accounting for all additives
        """
        correction = self.calculate_tm_correction(gc_content, primer_length)
        return tm_base + correction

    def calculate_effective_tm(self, seq: str, primer_conc: float = 0.5e-6) -> float:
        """
        Calculate effective Tm for sequence under these conditions.

        Applies all additive corrections including GC-dependent effects
        from TMAC and betaine. Uses actual primer length for accurate
        length-dependent GC correction.

        Args:
            seq: DNA sequence
            primer_conc: Primer concentration in M

        Returns:
            Effective Tm in degrees Celsius
        """
        # Calculate base Tm with salt
        tm_base = thermo.calculate_tm_with_salt(
            seq, self.na_conc, self.mg_conc, primer_conc
        )

        # Calculate GC content for GC-dependent corrections
        gc = thermo.gc_content(seq)

        # Use actual primer length for length-dependent GC correction
        primer_length = len(seq)

        # Apply additive corrections (including GC-dependent effects)
        tm_effective = self.adjust_tm(tm_base, gc, primer_length)

        return tm_effective

    def get_polymerase_range(self) -> Tuple[float, float]:
        """
        Get optimal temperature range for the polymerase.

        Returns:
            (min_temp, max_temp) in degrees Celsius
        """
        return POLYMERASE_CHARACTERISTICS[self.polymerase]['temp_range']

    def get_processivity(self) -> int:
        """
        Get processivity (maximum extension length) for the polymerase.

        Returns:
            Maximum extension length in base pairs

        Literature references:
            phi29: Blanco et al. (1989) JBC 264:8935 - ~70kb
            equiphi29: NEB data - ~80kb at elevated temp
            bst: Notomi et al. (2000) NAR 28:e63 - ~1-2kb
            klenow: Bambara et al. (1978) JBC 253:413 - ~10kb
        """
        return POLYMERASE_CHARACTERISTICS[self.polymerase]['processivity']

    def get_strand_displacement(self) -> bool:
        """
        Check if polymerase has strand displacement activity.

        Strand displacement is essential for SWGA/MDA.

        Returns:
            True if polymerase has strand displacement activity
        """
        return POLYMERASE_CHARACTERISTICS[self.polymerase]['strand_displacement']

    def get_exonuclease_activity(self) -> str:
        """
        Get exonuclease activity profile.

        Returns:
            Exonuclease activity: 'none', '3to5', '5to3', or 'both'
        """
        return POLYMERASE_CHARACTERISTICS[self.polymerase]['exonuclease']

    def get_fidelity(self) -> float:
        """
        Get polymerase fidelity (error rate).

        Returns:
            Error rate per nucleotide incorporated
        """
        return POLYMERASE_CHARACTERISTICS[self.polymerase]['error_rate']

    def effective_annealing_temp(self) -> float:
        """
        Calculate effective annealing temperature accounting for SSB.

        SSB prevents re-annealing, effectively lowering the
        temperature needed for primer binding.

        Returns:
            Effective annealing temperature in degrees Celsius
        """
        effective = self.temp

        # SSB reduces effective annealing temp by ~5°C
        if self.ssb:
            effective -= 5.0

        return effective

    def is_temp_in_polymerase_range(self) -> bool:
        """
        Check if reaction temperature is within polymerase optimal range.

        Returns:
            True if temperature is optimal
        """
        min_temp, max_temp = self.get_polymerase_range()
        return min_temp <= self.temp <= max_temp

    def max_primer_length(self) -> int:
        """
        Calculate maximum recommended primer length for these conditions.

        Longer primers provide exponentially better specificity:
        - 15bp: 4^15 = 1.07 billion combinations (baseline)
        - 16bp: 4^16 = 4.29 billion (4× improvement)
        - 17bp: 4^17 = 17.2 billion (16× improvement)
        - 18bp: 4^18 = 68.7 billion (64× improvement)

        SWGA (isothermal 30°C) can handle longer primers better than PCR (55-60°C)
        because lower reaction temperature = larger ΔT = more stable binding.

        Literature support:
        - Henke et al. (1997): 18bp primers work with 1M betaine at 55°C (PCR)
        - Musso et al. (2006): 20bp primers work with 2M betaine + 5% DMSO at 60°C (PCR)
        - If it works in PCR, it works better in SWGA due to lower temperature

        Returns:
            Maximum k-mer length (6-18)
        """
        base_max = 12  # Default swga2 maximum

        # Betaine allows longer primers by reducing GC-dependence
        if self.betaine_m >= 0.5:
            base_max += 2
        if self.betaine_m >= 1.0:
            base_max += 1

        # DMSO helps with secondary structure, enables longer primers
        if self.dmso_percent >= 3.0:
            base_max += 1
        if self.dmso_percent >= 5.0:
            base_max += 1

        # EquiPhi29 at higher temp can handle longer primers
        if self.polymerase == 'equiphi29' and self.temp >= 42:
            base_max += 1

        # Graduated cap based on additive support
        # Start with the incrementally calculated maximum
        supported_max = base_max

        # Apply safety caps for longer primers

        # For 16bp+: Need at least moderate additive support
        if supported_max >= 16:
            if not (self.betaine_m >= 1.0 or self.dmso_percent >= 5.0):
                supported_max = 15  # Insufficient support for 16bp

        # For 17bp+: Need strong combined support (Henke 1997 level)
        if supported_max >= 17:
            if not (self.betaine_m >= 1.5 and self.dmso_percent >= 3.0):
                supported_max = 16  # Insufficient support for 17bp

        # For 18bp+: Need extreme combined support (Musso 2006 level)
        if supported_max >= 18:
            if not (self.betaine_m >= 2.0 and self.dmso_percent >= 5.0):
                supported_max = 17  # Insufficient support for 18bp

        # Absolute maximum cap at 18bp (literature-validated limit)
        return min(supported_max, 18)

    def min_primer_length(self) -> int:
        """
        Calculate minimum recommended primer length.

        Returns:
            Minimum k-mer length (typically 6)
        """
        return 6  # Standard minimum for specificity

    def gc_content_range(self) -> Tuple[float, float]:
        """
        Calculate recommended GC content range.

        Betaine reduces GC-dependence, allowing wider range.

        Returns:
            (min_gc, max_gc) as fractions (0-1)
        """
        if self.betaine_m >= 1.0:
            # Betaine equalizes AT vs GC
            return (0.30, 0.70)
        else:
            # Standard range
            return (0.375, 0.625)

    def optimize_mg_concentration(self, genome_gc: float) -> float:
        """
        Optimize Mg2+ concentration based on genome GC content.

        Literature evidence (Rahman et al. 2014, Owczarzy et al. 2008):
        - AT-rich genomes (<35% GC): Higher Mg2+ stabilizes AT-rich duplexes
        - GC-rich genomes (>65% GC): Lower Mg2+ reduces non-specific binding
        - Balanced genomes: Standard 2.0 mM

        Args:
            genome_gc: Genome GC content (0-1)

        Returns:
            Optimal Mg2+ concentration (mM)

        Example:
            >>> conditions = ReactionConditions()
            >>> conditions.optimize_mg_concentration(0.33)  # Francisella
            2.5
            >>> conditions.optimize_mg_concentration(0.67)  # Burkholderia
            1.5
        """
        if genome_gc < 0 or genome_gc > 1:
            raise ValueError(f"GC content {genome_gc} outside valid range (0-1)")
        if genome_gc < 0.35:
            # AT-rich: need extra stabilization
            return 2.5  # mM
        elif genome_gc > 0.65:
            # GC-rich: reduce to minimize background
            return 1.5  # mM
        else:
            # Balanced
            return 2.0  # mM

    def estimate_coverage_improvement(self) -> float:
        """
        Estimate expected genome coverage improvement from conditions.

        Based on:
        - Longer primers (better specificity)
        - Wider GC range (more candidates)
        - Better polymerase performance

        Returns:
            Expected improvement factor (1.0 = no improvement)
        """
        improvement = 1.0

        # Longer primers increase specificity (fewer off-targets)
        max_length = self.max_primer_length()
        if max_length > 12:
            improvement *= 1.05 * (max_length - 12)  # 5% per extra base

        # Betaine improves evenness
        if self.betaine_m > 0:
            improvement *= 1.0 + 0.1 * self.betaine_m  # Up to 25% at 2.5M

        # EquiPhi29 gives better yield
        if self.polymerase == 'equiphi29':
            improvement *= 1.15

        return improvement

    def __repr__(self) -> str:
        """String representation of reaction conditions."""
        parts = [
            f"Temp={self.temp}C",
            f"Polymerase={self.polymerase}"
        ]

        if self.dmso_percent > 0:
            parts.append(f"DMSO={self.dmso_percent}%")
        if self.betaine_m > 0:
            parts.append(f"Betaine={self.betaine_m}M")
        if self.trehalose_m > 0:
            parts.append(f"Trehalose={self.trehalose_m}M")
        if self.formamide_percent > 0:
            parts.append(f"Formamide={self.formamide_percent}%")
        if self.glycerol_percent > 0:
            parts.append(f"Glycerol={self.glycerol_percent}%")
        if self.bsa_ug_ml > 0:
            parts.append(f"BSA={self.bsa_ug_ml}ug/mL")
        if self.peg_percent > 0:
            parts.append(f"PEG={self.peg_percent}%")
        if self.ethanol_percent > 0:
            parts.append(f"Ethanol={self.ethanol_percent}%")
        if self.urea_m > 0:
            parts.append(f"Urea={self.urea_m}M")
        if self.tmac_m > 0:
            parts.append(f"TMAC={self.tmac_m}M")
        if self.ssb:
            parts.append("SSB=True")

        parts.append(f"Na+={self.na_conc}mM")
        if self.mg_conc > 0:
            parts.append(f"Mg2+={self.mg_conc}mM")

        return f"ReactionConditions({', '.join(parts)})"

    def to_dict(self) -> Dict:
        """Convert to dictionary for JSON serialization."""
        return {
            'temp': self.temp,
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
            'na_conc': self.na_conc,
            'mg_conc': self.mg_conc,
            'ssb': self.ssb,
            'polymerase': self.polymerase
        }

    @classmethod
    def from_dict(cls, data: Dict) -> 'ReactionConditions':
        """Create from dictionary."""
        return cls(**data)

    @property
    def additives(self) -> AdditiveConcentrations:
        """
        Get additive concentrations as immutable AdditiveConcentrations object.

        Returns:
            AdditiveConcentrations with all additive values from this instance.

        Example:
            >>> conditions = get_enhanced_conditions()
            >>> additives = conditions.additives
            >>> correction = additives.calculate_tm_correction(gc_content=0.6)
        """
        return AdditiveConcentrations(
            dmso_percent=self.dmso_percent,
            betaine_m=self.betaine_m,
            trehalose_m=self.trehalose_m,
            formamide_percent=self.formamide_percent,
            glycerol_percent=self.glycerol_percent,
            bsa_ug_ml=self.bsa_ug_ml,
            peg_percent=self.peg_percent,
            ethanol_percent=self.ethanol_percent,
            urea_m=self.urea_m,
            tmac_m=self.tmac_m,
        )

    @classmethod
    def from_additives(
        cls,
        additives: AdditiveConcentrations,
        temp: float = 30.0,
        polymerase: str = 'phi29',
        na_conc: float = 50.0,
        mg_conc: float = 0.0,
        ssb: bool = False
    ) -> 'ReactionConditions':
        """
        Create ReactionConditions from AdditiveConcentrations.

        Useful when you have configured additives separately and want to
        combine with temperature and polymerase settings.

        Args:
            additives: AdditiveConcentrations with additive values
            temp: Reaction temperature in Celsius
            polymerase: Polymerase type
            na_conc: Na+ concentration in mM
            mg_conc: Mg2+ concentration in mM
            ssb: Whether SSB protein is present

        Returns:
            ReactionConditions with all settings

        Example:
            >>> additives = AdditiveConcentrations.for_enhanced_equiphi29()
            >>> conditions = ReactionConditions.from_additives(
            ...     additives, temp=42.0, polymerase='equiphi29'
            ... )
        """
        return cls(
            temp=temp,
            dmso_percent=additives.dmso_percent,
            betaine_m=additives.betaine_m,
            trehalose_m=additives.trehalose_m,
            formamide_percent=additives.formamide_percent,
            glycerol_percent=additives.glycerol_percent,
            bsa_ug_ml=additives.bsa_ug_ml,
            peg_percent=additives.peg_percent,
            ethanol_percent=additives.ethanol_percent,
            urea_m=additives.urea_m,
            tmac_m=additives.tmac_m,
            na_conc=na_conc,
            mg_conc=mg_conc,
            ssb=ssb,
            polymerase=polymerase,
        )


# ========================================
# Preset Conditions
# ========================================

def get_standard_conditions() -> ReactionConditions:
    """Standard phi29 SWGA conditions (original swga2)."""
    return ReactionConditions(
        temp=30.0,
        polymerase='phi29',
        na_conc=50.0
    )


def get_equiphi_conditions() -> ReactionConditions:
    """EquiPhi29 optimized conditions."""
    return ReactionConditions(
        temp=42.0,
        polymerase='equiphi29',
        na_conc=50.0
    )


def get_enhanced_conditions() -> ReactionConditions:
    """
    Enhanced conditions with additives for longer primers.

    Uses betaine and DMSO to enable 13-15mer primers.
    """
    return ReactionConditions(
        temp=42.0,
        dmso_percent=5.0,
        betaine_m=1.0,
        polymerase='equiphi29',
        na_conc=50.0
    )


def get_high_gc_conditions() -> ReactionConditions:
    """
    Optimized for high-GC genomes.

    Uses high betaine to equalize AT/GC stability.
    """
    return ReactionConditions(
        temp=42.0,
        dmso_percent=5.0,
        betaine_m=2.0,
        polymerase='equiphi29',
        na_conc=50.0
    )


def get_low_temp_conditions() -> ReactionConditions:
    """
    Low temperature conditions with SSB.

    For sensitive samples or thermolabile templates.
    """
    return ReactionConditions(
        temp=30.0,
        dmso_percent=3.0,
        betaine_m=1.0,
        ssb=True,
        polymerase='phi29',
        na_conc=50.0
    )


def get_q_solution_equivalent() -> ReactionConditions:
    """
    Qiagen Q-Solution equivalent conditions.

    Q-Solution contains proprietary betaine mix plus stabilizers.
    Approximation: 1.5M betaine plus glycerol stabilizer.

    Use for: Difficult templates, high GC content, secondary structure.
    Expected improvement: 10-20% better amplification of GC-rich regions.
    """
    return ReactionConditions(
        temp=42.0,
        betaine_m=1.5,
        glycerol_percent=10.0,
        bsa_ug_ml=200.0,
        polymerase='equiphi29',
        na_conc=50.0
    )


def get_gc_melt_conditions() -> ReactionConditions:
    """
    GC-Melt equivalent for extreme GC genomes.

    Combines DMSO, betaine, and trehalose for maximum GC tolerance.

    Use for: >70% GC genomes, strong secondary structure.
    Expected improvement: 20-30% better coverage of GC-rich regions.
    """
    return ReactionConditions(
        temp=45.0,
        dmso_percent=5.0,
        betaine_m=2.0,
        trehalose_m=0.5,
        mg_conc=1.5,
        polymerase='equiphi29',
        na_conc=50.0
    )


def get_crude_sample_conditions() -> ReactionConditions:
    """
    Optimized for crude or inhibitor-containing samples.

    High BSA, moderate glycerol, standard additives.

    Use for: Crude lysates, blood, soil, inhibitor-rich samples.
    Expected improvement: 20-30% better performance vs standard conditions.
    """
    return ReactionConditions(
        temp=30.0,
        betaine_m=1.0,
        glycerol_percent=10.0,
        bsa_ug_ml=400.0,
        peg_percent=5.0,
        polymerase='phi29',
        na_conc=50.0
    )


def get_bst_conditions() -> ReactionConditions:
    """
    Bst polymerase conditions for high-temperature isothermal amplification.

    Bst 2.0/3.0 operates at 60-65C with strand displacement capability.
    Lower processivity (~1-2kb) than phi29 but higher temperature tolerance.

    Use for:
        - Samples with phi29 inhibitors
        - Integration with LAMP workflows
        - Higher specificity due to elevated temperature

    Reference: Notomi et al. (2000) Nucleic Acids Res 28:e63
    """
    return ReactionConditions(
        temp=63.0,
        betaine_m=0.8,  # Moderate betaine for high-temp stability
        mg_conc=8.0,    # Bst requires higher Mg2+ (typically 6-10mM)
        polymerase='bst',
        na_conc=50.0
    )


def get_klenow_conditions() -> ReactionConditions:
    """
    Klenow exo- polymerase conditions.

    Klenow fragment (exonuclease deficient) has moderate processivity (~10kb)
    and operates at 37C. Lower cost alternative to phi29.

    Use for:
        - Budget-constrained applications
        - 37C incubation (standard lab equipment)
        - Moderate amplification requirements

    Note: Lower processivity than phi29 may result in more fragmented products.
    """
    return ReactionConditions(
        temp=37.0,
        betaine_m=1.0,
        dmso_percent=5.0,
        mg_conc=10.0,   # Klenow requires higher Mg2+ (typically 10mM)
        polymerase='klenow',
        na_conc=50.0
    )


def get_extreme_gc_conditions() -> ReactionConditions:
    """
    Conditions for extreme GC genomes (>70% or <30% GC).

    Uses TMAC to equalize AT/GC Tm, combined with urea for denaturation.

    Use for:
        - Very high GC genomes (>70%): Mycobacterium, Streptomyces
        - Very low GC genomes (<30%): Some Firmicutes, Plasmodium

    Reference: Melchior & von Hippel (1973) PNAS 70:298-302
    """
    return ReactionConditions(
        temp=42.0,
        dmso_percent=5.0,
        betaine_m=2.0,
        urea_m=0.5,     # Moderate urea for denaturation
        tmac_m=0.05,    # Low TMAC for AT/GC equalization
        polymerase='equiphi29',
        na_conc=50.0
    )


# ========================================
# Condition Optimization
# ========================================

def optimize_conditions_for_primers(primers: list,
                                   target_tm_range: Tuple[float, float] = (30, 45),
                                   polymerase: str = 'phi29') -> ReactionConditions:
    """
    Optimize reaction conditions to bring primer Tm into target range.

    Args:
        primers: List of primer sequences
        target_tm_range: Desired (min, max) Tm in °C
        polymerase: Polymerase to use

    Returns:
        Optimized ReactionConditions
    """
    # Calculate current Tm distribution
    tms = [thermo.calculate_tm_with_salt(p, na_conc=50) for p in primers]
    mean_tm = np.mean(tms)
    target_mean = np.mean(target_tm_range)

    # Calculate needed Tm shift
    tm_shift_needed = target_mean - mean_tm

    # Optimize additives to achieve shift
    conditions = ReactionConditions(polymerase=polymerase)

    if tm_shift_needed < -15:
        # Need large Tm reduction: use betaine + DMSO + formamide
        conditions.betaine_m = 2.0  # -4.6°C
        conditions.dmso_percent = 10.0  # -6.0°C
        conditions.formamide_percent = 5.0  # -3.25°C
    elif tm_shift_needed < -8:
        # Moderate reduction: betaine + DMSO
        conditions.betaine_m = 1.5  # -3.45°C
        conditions.dmso_percent = 5.0  # -3.0°C
    elif tm_shift_needed < -4:
        # Small reduction: DMSO only
        conditions.dmso_percent = abs(tm_shift_needed) / 0.6
    elif tm_shift_needed > 5:
        # Need higher Tm: use EquiPhi29 at higher temp
        conditions.polymerase = 'equiphi29'
        conditions.temp = 42.0

    return conditions


def recommend_conditions(genome_seq: str, target_k: Optional[int] = None) -> Dict:
    """
    Analyze genome and recommend optimal reaction conditions.

    Args:
        genome_seq: Genome sequence string
        target_k: Target k-mer length (optional, will recommend if not provided)

    Returns:
        Dictionary with recommended conditions and analysis
    """
    # Calculate genome properties
    genome_length = len(genome_seq)
    gc_content = (genome_seq.count('G') + genome_seq.count('C')) / genome_length

    # Recommend conditions based on GC content
    if gc_content < 0.30:
        # AT-rich genome (e.g., Francisella 32% GC)
        conditions = ReactionConditions(
            temp=37.0,
            dmso_percent=3.0,
            betaine_m=0.8,
            mg_conc=2.5,  # Extra Mg2+ for AT-rich
            polymerase='phi29',
            na_conc=50.0
        )
        optimal_k = 10
        rationale = "AT-rich genome: moderate additives + extra Mg2+ for stability"

    elif gc_content > 0.65:
        # GC-rich genome (e.g., Burkholderia 67% GC)
        conditions = get_high_gc_conditions()
        optimal_k = 12
        rationale = "GC-rich genome: high betaine + DMSO to equalize AT/GC stability"

    elif 0.48 < gc_content < 0.52:
        # Balanced genome (e.g., E. coli 50% GC)
        conditions = get_standard_conditions()
        optimal_k = 11
        rationale = "Balanced GC: standard phi29 conditions sufficient"

    else:
        # Moderately AT or GC-biased
        conditions = get_enhanced_conditions()
        optimal_k = 11
        rationale = "Moderate GC bias: enhanced conditions recommended"

    # Override if user specified target k
    if target_k is not None:
        optimal_k = target_k
        # Validate and adjust conditions if needed
        max_supported = conditions.max_primer_length()
        if target_k > max_supported:
            # Need stronger additives
            if target_k >= 15:
                conditions.betaine_m = max(2.0, conditions.betaine_m)
                conditions.dmso_percent = max(7.0, conditions.dmso_percent)
                conditions.polymerase = 'equiphi29'
                conditions.temp = 42.0
            elif target_k >= 13:
                conditions.betaine_m = max(1.5, conditions.betaine_m)
                conditions.dmso_percent = max(5.0, conditions.dmso_percent)

    # Prepare recommendations
    return {
        'optimal_k': optimal_k,
        'genome_gc': gc_content,
        'genome_length': genome_length,
        'temperature': conditions.temp,
        'polymerase': conditions.polymerase,
        'dmso_percent': conditions.dmso_percent,
        'betaine_m': conditions.betaine_m,
        'mg_conc': conditions.optimize_mg_concentration(gc_content),
        'na_conc': conditions.na_conc,
        'ssb': conditions.ssb,
        'max_primer_length': conditions.max_primer_length(),
        'rationale': rationale,
        'conditions': conditions
    }


if __name__ == "__main__":
    # Test reaction conditions
    print("Standard conditions:")
    std = get_standard_conditions()
    print(std)
    print(f"Max primer length: {std.max_primer_length()}")
    print()

    print("Enhanced conditions:")
    enh = get_enhanced_conditions()
    print(enh)
    print(f"Max primer length: {enh.max_primer_length()}")
    print(f"Expected coverage improvement: {enh.estimate_coverage_improvement():.2f}x")
    print()

    # Test Tm calculation
    test_seq = "ATCGATCGATCGATCG"  # 16-mer
    print(f"Test sequence: {test_seq}")
    print(f"Base Tm: {thermo.calculate_tm_with_salt(test_seq):.1f}°C")
    print(f"Tm with standard conditions: {std.calculate_effective_tm(test_seq):.1f}°C")
    print(f"Tm with enhanced conditions: {enh.calculate_effective_tm(test_seq):.1f}°C")
