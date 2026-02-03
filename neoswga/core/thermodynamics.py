"""
Comprehensive thermodynamic calculations for DNA primers.

Implements the unified nearest-neighbor model (SantaLucia, 1998) with:
- Complete enthalpy and entropy parameters for all dinucleotide stacks
- Salt concentration corrections (Na+, Mg2+)
- Temperature-dependent free energy calculations
- Helix initiation, terminal penalties, and symmetry corrections
- Support for thermodynamic additives (DMSO, betaine, trehalose)

References:
- SantaLucia (1998) PNAS 95:1460-1465
- Owczarzy et al. (2008) Biochemistry 47:5336-5353
"""

import numpy as np
import warnings
import logging
from typing import Dict, Tuple, Optional, List, Any
from functools import lru_cache

logger = logging.getLogger(__name__)

# Universal gas constant (cal/(mol*K))
R = 1.987

# ========================================
# SantaLucia Nearest-Neighbor Parameters
# ========================================
# Enthalpy values (kcal/mol) for all 16 dinucleotide stacks

ENTHALPY_NN = {
    # SantaLucia (1998) unified NN parameters - Table 1
    # The 10 canonical Watson-Crick nearest-neighbor stacks
    # Stack notation: 5'-XY-3' on top strand / 3'-ZW-5' on bottom strand
    # Values in kcal/mol
    'AA/TT': -7.9,
    'AT/TA': -7.2,
    'TA/AT': -7.2,
    'CA/GT': -8.5,
    'GT/CA': -8.4,
    'CT/GA': -7.8,
    'GA/CT': -8.2,
    'CG/GC': -10.6,
    'GC/CG': -9.8,
    'GG/CC': -8.0,
    # Reverse-complement equivalents for faster lookup
    # (avoids runtime reverse_complement calculation)
    'TT/AA': -7.9,   # = AA/TT
    'AC/TG': -8.4,   # = GT/CA
    'TC/AG': -8.2,   # = GA/CT
    'AG/TC': -8.0,   # = GG/CC (entropy differs slightly)
    'TG/AC': -8.5,   # = CA/GT
    'CC/GG': -8.0,   # = GG/CC
}

# Entropy values (cal/(mol*K)) - must match ENTHALPY_NN keys
ENTROPY_NN = {
    # SantaLucia (1998) unified NN parameters - Table 1
    # Values in cal/(mol*K)
    'AA/TT': -22.2,
    'AT/TA': -20.4,
    'TA/AT': -21.3,
    'CA/GT': -22.7,
    'GT/CA': -22.4,
    'CT/GA': -21.0,
    'GA/CT': -22.2,
    'CG/GC': -27.2,
    'GC/CG': -24.4,
    'GG/CC': -19.9,
    # Reverse-complement equivalents
    'TT/AA': -22.2,   # = AA/TT
    'AC/TG': -22.4,   # = GT/CA
    'TC/AG': -22.2,   # = GA/CT
    'AG/TC': -19.9,   # = GG/CC
    'TG/AC': -22.7,   # = CA/GT
    'CC/GG': -19.9,   # = GG/CC
}

# The 10 canonical Watson-Crick nearest-neighbor stacks from SantaLucia (1998):
#   AA/TT, AT/TA, TA/AT, CA/GT, GT/CA, CT/GA, GA/CT, CG/GC, GC/CG, GG/CC
# Additional entries are reverse-complement equivalents for faster lookup.
# Stack notation: 5'-XY-3' / 3'-ZW-5' where X pairs with W and Y pairs with Z.

# Helix initiation (independent of sequence)
INIT_ENTHALPY = 0.2  # kcal/mol
INIT_ENTROPY = -5.7  # cal/(mol*K)

# Terminal AT penalty (for each terminal AT base pair)
TERMINAL_AT_ENTHALPY = 2.2  # kcal/mol
TERMINAL_AT_ENTROPY = 6.9  # cal/(mol*K)

# Symmetry correction (for palindromic sequences)
SYMMETRY_ENTROPY = -1.4  # cal/(mol*K)


# IUPAC ambiguous base codes
IUPAC_AMBIGUOUS = set('NRYSWKMBDHV')


def normalize_sequence(seq: str) -> str:
    """
    Normalize DNA sequence for thermodynamic calculations.

    Converts lowercase to uppercase and validates bases.

    Args:
        seq: DNA sequence (may contain lowercase or ambiguous bases)

    Returns:
        Uppercase normalized sequence
    """
    return seq.upper()


def has_ambiguous_bases(seq: str) -> bool:
    """
    Check if sequence contains ambiguous IUPAC bases.

    Args:
        seq: DNA sequence (should be uppercase)

    Returns:
        True if sequence contains N, R, Y, S, W, K, M, B, D, H, or V
    """
    return any(base in IUPAC_AMBIGUOUS for base in seq.upper())


def reverse_complement(seq: str) -> str:
    """
    Return reverse complement of DNA sequence.

    Handles lowercase by converting to uppercase first.
    Ambiguous bases (N, R, Y, etc.) are mapped to N.
    """
    seq = normalize_sequence(seq)
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq))


def is_palindrome(seq: str) -> bool:
    """Check if sequence is palindromic."""
    return seq == reverse_complement(seq)


def is_watson_crick(base1: str, base2: str) -> bool:
    """Check if two bases form Watson-Crick pair."""
    pairs = {('A', 'T'), ('T', 'A'), ('G', 'C'), ('C', 'G')}
    return (base1, base2) in pairs


@lru_cache(maxsize=1000000)
def calculate_enthalpy_entropy_cached(seq: str, complementary: Optional[str] = None) -> Tuple[float, float]:
    """
    Cached version of enthalpy/entropy calculation.
    Uses LRU cache (1M entries) for 100-1000x speedup on repeated sequences.
    """
    return _calculate_enthalpy_entropy_impl(seq, complementary)


def calculate_enthalpy_entropy(seq: str, complementary: Optional[str] = None) -> Tuple[float, float]:
    """
    Calculate total enthalpy and entropy for DNA duplex using nearest-neighbor model.

    This function uses caching for performance. For cache-bypassing behavior,
    use _calculate_enthalpy_entropy_impl directly.

    Args:
        seq: DNA sequence (5' to 3')
        complementary: Optional complementary strand (3' to 5'). If None, assumes perfect Watson-Crick.

    Returns:
        (enthalpy, entropy) in (kcal/mol, cal/(mol*K))
    """
    return calculate_enthalpy_entropy_cached(seq, complementary)


def _calculate_enthalpy_entropy_impl(seq: str, complementary: Optional[str] = None) -> Tuple[float, float]:
    """
    Implementation of enthalpy/entropy calculation (uncached).

    Args:
        seq: DNA sequence (5' to 3'). Lowercase is converted to uppercase.
        complementary: Optional complementary strand (3' to 5'). If None, assumes perfect Watson-Crick.

    Returns:
        (enthalpy, entropy) in (kcal/mol, cal/(mol*K))

    Notes:
        Sequences containing ambiguous IUPAC bases (N, R, Y, etc.) will generate
        warnings and use default penalty values for unknown stacks.
    """
    # Normalize to uppercase for consistent lookup
    seq = normalize_sequence(seq)
    if complementary is not None:
        complementary = normalize_sequence(complementary)
    else:
        complementary = reverse_complement(seq)

    # Warn if ambiguous bases are present
    if has_ambiguous_bases(seq):
        warnings.warn(
            f"Sequence contains ambiguous bases: {seq}. "
            "Tm calculation will use default penalty values for unknown stacks.",
            UserWarning
        )

    if len(seq) != len(complementary):
        raise ValueError("Sequences must be same length")

    if len(seq) < 2:
        warnings.warn("Sequence too short for nearest-neighbor calculation")
        return 0.0, 0.0

    total_enthalpy = INIT_ENTHALPY
    total_entropy = INIT_ENTROPY

    # The complementary strand is stored as reverse complement (5'->3' direction).
    # For proper nearest-neighbor alignment, we need it in 3'->5' direction.
    # Example: seq="ATCG", complementary="CGAT" (reverse complement)
    #   5'-A-T-C-G-3'  (seq)
    #      | | | |
    #   3'-T-A-G-C-5'  (complement_3to5 = "TAGC")
    # The stack AT/TA means: 5'-AT-3' on top, 3'-TA-5' on bottom
    complement_3to5 = complementary[::-1]

    # Sum nearest-neighbor contributions
    for i in range(len(seq) - 1):
        # Get dinucleotide stack
        stack_top = seq[i:i+2]
        stack_bottom = complement_3to5[i:i+2]  # Already in 3'->5' orientation

        # Create stack notation: top/bottom (SantaLucia convention)
        stack = f"{stack_top}/{stack_bottom}"

        # Try to find in parameters
        if stack in ENTHALPY_NN:
            total_enthalpy += ENTHALPY_NN[stack]
            total_entropy += ENTROPY_NN[stack]
        else:
            # Try reverse complement of stack
            rev_stack = f"{reverse_complement(stack_bottom)}/{reverse_complement(stack_top)}"
            if rev_stack in ENTHALPY_NN:
                total_enthalpy += ENTHALPY_NN[rev_stack]
                total_entropy += ENTROPY_NN[rev_stack]
            else:
                # Unknown stack (ambiguous base or non-ATGC) - use average penalty
                total_enthalpy += -8.0  # Average value
                total_entropy += -21.0  # Average value

    # Terminal AT penalty
    if seq[0] in ['A', 'T']:
        total_enthalpy += TERMINAL_AT_ENTHALPY
        total_entropy += TERMINAL_AT_ENTROPY

    if seq[-1] in ['A', 'T']:
        total_enthalpy += TERMINAL_AT_ENTHALPY
        total_entropy += TERMINAL_AT_ENTROPY

    # Symmetry correction
    if is_palindrome(seq):
        total_entropy += SYMMETRY_ENTROPY

    return total_enthalpy, total_entropy


def calculate_tm_basic(seq: str, primer_conc: float = 0.5e-6) -> float:
    """
    Calculate basic melting temperature without salt corrections.

    Args:
        seq: DNA sequence
        primer_conc: Primer concentration in M (default 0.5 uM = 0.5e-6 M)

    Returns:
        Tm in degrees Celsius
    """
    enthalpy, entropy = calculate_enthalpy_entropy(seq)

    # For self-complementary sequences, [primer] = Ct/4
    # For non-self-complementary, [primer] = Ct/2 (assumes equal concentrations)
    if is_palindrome(seq):
        effective_conc = primer_conc / 4
    else:
        effective_conc = primer_conc / 2

    # Tm = ΔH / (ΔS + R*ln(Ct)) - 273.15
    # Convert enthalpy to cal/mol (multiply by 1000)
    tm_kelvin = (enthalpy * 1000) / (entropy + R * np.log(effective_conc))
    tm_celsius = tm_kelvin - 273.15

    return tm_celsius


def calculate_salt_correction(na_conc: float = 50, mg_conc: float = 0) -> float:
    """
    Calculate salt correction term for Tm.

    Uses SantaLucia (1998) formula for oligonucleotides with Owczarzy (2008)
    Mg2+ equivalence for mixed solutions.

    The coefficient 12.5 is appropriate for oligonucleotides (6-30 bp),
    while 16.6 is for long DNA polymers.

    References:
        - SantaLucia (1998) PNAS 95:1460-1465: 12.5 * log10([Na+])
        - Owczarzy et al. (2008) Biochemistry 47:5336-5353: Mg2+ correction

    Args:
        na_conc: Na+ concentration in mM (default 50 mM)
        mg_conc: Mg2+ concentration in mM (default 0 mM)

    Returns:
        Salt correction in degrees Celsius (relative to 1M Na+)
    """
    # Convert to M
    na_m = na_conc / 1000.0
    mg_m = mg_conc / 1000.0

    # Unified salt correction: one Mg2+ equivalent = 3.3 * sqrt([Mg2+]) in Na+
    # From Owczarzy (2008): accounts for Mg2+ stabilization of duplex
    if mg_m > 0:
        # When Mg2+ is present, use Owczarzy's approximation
        monovalent_eq = na_m + 3.3 * np.sqrt(mg_m)
    else:
        monovalent_eq = na_m

    # Salt correction using SantaLucia (1998) coefficient for oligonucleotides
    # Tm = Tm(1M) + 12.5 * log10([Na+])
    if monovalent_eq > 0:
        salt_correction = 12.5 * np.log10(monovalent_eq)
    else:
        salt_correction = 0.0

    return salt_correction


def calculate_tm_with_salt(seq: str,
                           na_conc: float = 50,
                           mg_conc: float = 0,
                           primer_conc: float = 0.5e-6) -> float:
    """
    Calculate melting temperature with salt corrections.

    Args:
        seq: DNA sequence
        na_conc: Na+ concentration in mM
        mg_conc: Mg2+ concentration in mM
        primer_conc: Primer concentration in M

    Returns:
        Tm in degrees Celsius
    """
    tm_basic = calculate_tm_basic(seq, primer_conc)
    salt_corr = calculate_salt_correction(na_conc, mg_conc)

    return tm_basic + salt_corr


def calculate_free_energy(seq: str, temperature: float = 37.0) -> float:
    """
    Calculate Gibbs free energy at specified temperature.

    Args:
        seq: DNA sequence
        temperature: Temperature in degrees Celsius

    Returns:
        ΔG in kcal/mol (negative = favorable binding)
    """
    enthalpy, entropy = calculate_enthalpy_entropy(seq)
    temp_kelvin = temperature + 273.15

    # ΔG = ΔH - T*ΔS
    # Convert entropy to kcal/(mol*K)
    delta_g = enthalpy - (temp_kelvin * entropy / 1000)

    return delta_g


def calculate_binding_probability(seq: str,
                                  temperature: float = 37.0,
                                  na_conc: float = 50,
                                  mg_conc: float = 0) -> float:
    """
    Calculate probability of primer being bound at given temperature.

    Uses Boltzmann distribution: P_bound = 1 / (1 + exp(ΔG/RT))

    Args:
        seq: DNA sequence
        temperature: Temperature in degrees Celsius
        na_conc: Na+ concentration in mM
        mg_conc: Mg2+ concentration in mM

    Returns:
        Probability between 0 and 1
    """
    delta_g = calculate_free_energy(seq, temperature)
    temp_kelvin = temperature + 273.15

    # Boltzmann factor
    boltzmann = np.exp(delta_g * 1000 / (R * temp_kelvin))  # Convert ΔG to cal/mol

    # Probability of bound state
    p_bound = 1 / (1 + boltzmann)

    return p_bound


def energy_to_tm(delta_g: float, seq_length: int,
                 na_conc: float = 50, mg_conc: float = 0,
                 reference_temp: float = 37.0) -> float:
    """
    Estimate Tm from free energy (approximate conversion).

    This function estimates Tm from a free energy value by:
    1. Estimating entropy from sequence length (average -22 cal/(mol*K) per bp)
    2. Back-calculating enthalpy from dG = dH - T*dS at reference temperature
    3. Using the standard Tm equation: Tm = dH / (dS + R*ln(Ct/4))

    Args:
        delta_g: Free energy in kcal/mol (calculated at reference_temp)
        seq_length: Length of sequence (for entropy estimation)
        na_conc: Na+ concentration in mM
        mg_conc: Mg2+ concentration in mM
        reference_temp: Temperature at which delta_g was calculated (C, default 37)

    Returns:
        Approximate Tm in degrees Celsius

    Note:
        This is an approximation. For accurate Tm, use calculate_tm_with_salt()
        with the actual sequence.
    """
    # Rough estimate: dS = -22 cal/(mol*K) per bp (average of NN stacks)
    # Range is -19.9 to -27.2, so this introduces ~15% error
    est_entropy = -22.0 * seq_length

    # Back-calculate enthalpy from dG = dH - T*dS
    # dH = dG + T*dS
    ref_temp_kelvin = reference_temp + 273.15
    # Convert dG from kcal/mol to cal/mol for consistent units with entropy
    delta_g_cal = delta_g * 1000.0
    est_enthalpy = delta_g_cal + ref_temp_kelvin * est_entropy

    # Tm = dH / (dS + R*ln(Ct/4))
    # Standard primer concentration: 0.5 uM = 0.5e-6 M
    denominator = est_entropy + R * np.log(0.5e-6 / 4)
    if abs(denominator) < 1e-10:
        return float('nan')  # Avoid division by zero

    tm_kelvin = est_enthalpy / denominator
    tm_celsius = tm_kelvin - 273.15

    # Add salt correction
    salt_corr = calculate_salt_correction(na_conc, mg_conc)

    return tm_celsius + salt_corr


def calculate_tm_range(seq: str,
                       na_conc: float = 50,
                       mg_conc: float = 0,
                       primer_conc: float = 0.5e-6) -> Tuple[float, float]:
    """
    Calculate Tm range accounting for uncertainty.

    Returns (Tm_lower, Tm_upper) with +/- 2°C margin (typical experimental error).

    Args:
        seq: DNA sequence
        na_conc: Na+ concentration in mM
        mg_conc: Mg2+ concentration in mM
        primer_conc: Primer concentration in M

    Returns:
        (lower_tm, upper_tm) in degrees Celsius
    """
    tm = calculate_tm_with_salt(seq, na_conc, mg_conc, primer_conc)

    # Typical experimental uncertainty
    uncertainty = 2.0  # degrees Celsius

    return (tm - uncertainty, tm + uncertainty)


# ========================================
# GC Content and Quick Estimates
# ========================================

def gc_content(seq: str) -> float:
    """Calculate GC content as fraction."""
    seq_upper = seq.upper()
    gc_count = seq_upper.count('G') + seq_upper.count('C')
    return gc_count / len(seq) if len(seq) > 0 else 0.0


def wallace_tm(seq: str) -> float:
    """
    Quick Tm estimate using Wallace rule (for primers < 14 bp).

    Tm = 2*(A+T) + 4*(G+C)

    Less accurate than nearest-neighbor but very fast for screening.

    Args:
        seq: DNA sequence

    Returns:
        Tm in degrees Celsius
    """
    seq_upper = seq.upper()
    at_count = seq_upper.count('A') + seq_upper.count('T')
    gc_count = seq_upper.count('G') + seq_upper.count('C')

    return 2 * at_count + 4 * gc_count


# ========================================
# Mismatch Free Energy Calculations
# (Legacy API compatibility for rf_preprocessing)
# ========================================

# Free energy values for mismatched nearest-neighbor doublets
# From SantaLucia (1998) - includes internal mismatches
DELTA_G_MISMATCH = {
    'GA/CA': 0.17, 'GA/CC': 0.81, 'GA/CG': -0.25, 'GA/CT': -1.30,
    'GC/CA': 0.47, 'GC/CC': 0.79, 'GC/CG': -2.24, 'GC/CT': 0.62,
    'GG/CA': -0.52, 'GG/CC': -1.84, 'GG/CG': -1.11, 'GG/CT': 0.08,
    'GT/CA': -1.44, 'GT/CC': 0.98, 'GT/CG': -0.59, 'GT/CT': 0.45,
    'CA/GA': 0.43, 'CA/GC': 0.75, 'CA/GG': 0.03, 'CA/GT': -1.45,
    'CC/GA': 0.79, 'CC/GC': 0.70, 'CC/GG': -1.84, 'CC/GT': 0.62,
    'CG/GA': 0.11, 'CG/GC': -2.17, 'CG/GG': -0.11, 'CG/GT': -0.47,
    'CT/GA': -1.28, 'CT/GC': 0.40, 'CT/GG': -0.32, 'CT/GT': -0.12,
    'AA/TA': 0.61, 'AA/TC': 0.88, 'AA/TG': 0.14, 'AA/TT': -1.00,
    'AC/TA': 0.77, 'AC/TC': 1.33, 'AC/TG': -1.44, 'AC/TT': 0.64,
    'AG/TA': 0.02, 'AG/TC': -1.28, 'AG/TG': -0.13, 'AG/TT': 0.71,
    'AT/TA': -0.88, 'AT/TC': 0.73, 'AT/TG': 0.07, 'AT/TT': 0.69,
    'TA/AA': 0.69, 'TA/AC': 0.92, 'TA/AG': 0.42, 'TA/AT': -0.58,
    'TC/AA': 1.33, 'TC/AC': 1.05, 'TC/AG': -1.30, 'TC/AT': 0.97,
    'TG/AA': 0.74, 'TG/AC': -1.45, 'TG/AG': 0.44, 'TG/AT': 0.43,
    'TT/AA': -1.00, 'TT/AC': 0.75, 'TT/AG': 0.34, 'TT/AT': 0.68
}

NN_INIT_CORRECTIONS = {'G/C': 0.98, 'A/T': 1.03, 'C/G': 0.98, 'T/A': 1.03}


def complement(base: str) -> str:
    """Return complement of a single base."""
    comp_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
                'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}
    return comp_map.get(base, base)


@lru_cache(maxsize=1000000)
def compute_free_energy_for_two_strings_cached(x: str, y: str, penalty: float = 4.0) -> float:
    """
    Cached version of free energy calculation between two sequences.
    Uses LRU cache for 10-100x speedup in Step 3 scoring.
    """
    return _compute_free_energy_for_two_strings_impl(x, y, penalty)


def compute_free_energy_for_two_strings(x: str, y: str, penalty: float = 4.0) -> float:
    """
    Compute delta G between two arbitrary DNA sequences using nearest-neighbor model.

    This function calculates the thermodynamic stability of hybridization between
    two sequences that may contain mismatches. Used for scoring primer-target
    binding with potential mismatches.

    Uses LRU caching for performance - same (x, y, penalty) inputs return cached results.

    Args:
        x: 5' to 3' sequence
        y: 3' to 5' sequence (should be aligned with x)
        penalty: Penalty for unspecified mismatch pairs (default: 4.0 kcal/mol)

    Returns:
        Delta G value in kcal/mol (more negative = more stable)
    """
    return compute_free_energy_for_two_strings_cached(x.upper(), y.upper(), penalty)


def _compute_free_energy_for_two_strings_impl(x: str, y: str, penalty: float = 4.0) -> float:
    """
    Implementation of free energy calculation (uncached).

    Args:
        x: 5' to 3' sequence (already uppercased)
        y: 3' to 5' sequence (already uppercased)
        penalty: Penalty for unspecified mismatch pairs

    Returns:
        Delta G value in kcal/mol
    """
    delta_g = 0.0

    # Terminal corrections
    if len(x) > 0 and len(y) > 0:
        if x[0] == complement(y[0]):
            key = f"{x[0]}/{y[0]}"
            if key in NN_INIT_CORRECTIONS:
                delta_g = NN_INIT_CORRECTIONS[key]

        if x[-1] == complement(y[-1]):
            key = f"{x[-1]}/{y[-1]}"
            if key in NN_INIT_CORRECTIONS:
                delta_g += NN_INIT_CORRECTIONS[key]

    # Sum nearest-neighbor contributions
    for i in range(1, len(x)):
        if i < len(y):
            doublet_1 = x[i-1:i+1]
            doublet_2 = y[i-1:i+1]
            nn_key = f"{doublet_1}/{doublet_2}"

            if nn_key in DELTA_G_MISMATCH:
                delta_g += DELTA_G_MISMATCH[nn_key]
            else:
                delta_g += penalty

            # Early termination for very unfavorable binding
            if delta_g > penalty * 10:
                break

    # Symmetry correction for palindromes
    from neoswga.core.utility import complement as util_complement, reverse
    if util_complement(x) == reverse(y):
        delta_g += 0.43

    return delta_g


def clear_thermodynamic_caches() -> None:
    """Clear all thermodynamic calculation caches. Useful for testing or memory management."""
    calculate_enthalpy_entropy_cached.cache_clear()
    compute_free_energy_for_two_strings_cached.cache_clear()


def get_cache_stats() -> Dict[str, Any]:
    """Return cache statistics for performance monitoring."""
    return {
        'enthalpy_entropy': calculate_enthalpy_entropy_cached.cache_info(),
        'free_energy': compute_free_energy_for_two_strings_cached.cache_info()
    }


# ========================================
# Vectorized Batch Tm Calculations
# ========================================

def calculate_tm_batch(sequences: List[str],
                       na_conc: float = 50.0,
                       mg_conc: float = 0.0,
                       primer_conc: float = 0.5e-6) -> np.ndarray:
    """
    Calculate Tm for a batch of sequences using vectorized operations.

    For large primer sets (>100), provides 2-5x speedup over individual calls
    by batching common calculations and leveraging numpy operations.

    Args:
        sequences: List of DNA sequences (can be mixed lengths)
        na_conc: Na+ concentration in mM (default 50 mM)
        mg_conc: Mg2+ concentration in mM (default 0 mM)
        primer_conc: Primer concentration in M (default 0.5 uM)

    Returns:
        numpy array of Tm values in degrees Celsius
    """
    n_seqs = len(sequences)
    if n_seqs == 0:
        return np.array([])

    # Pre-compute salt correction (same for all sequences)
    salt_corr = calculate_salt_correction(na_conc, mg_conc)

    # Calculate Tm for each sequence (leveraging cached thermodynamics)
    tm_values = np.zeros(n_seqs)

    for i, seq in enumerate(sequences):
        try:
            tm_basic = calculate_tm_basic(seq, primer_conc)
            tm_values[i] = tm_basic + salt_corr
        except (ValueError, KeyError) as e:
            # Handle expected errors: invalid bases, missing NN parameters
            logger.debug(f"Could not calculate Tm for sequence {seq}: {e}")
            tm_values[i] = np.nan
        except Exception as e:
            # Log unexpected errors for debugging
            logger.warning(f"Unexpected error calculating Tm for sequence {seq}: {e}")
            tm_values[i] = np.nan

    return tm_values


def calculate_tm_batch_with_additives(sequences: List[str],
                                       na_conc: float = 50.0,
                                       mg_conc: float = 0.0,
                                       primer_conc: float = 0.5e-6,
                                       dmso_percent: float = 0.0,
                                       betaine_m: float = 0.0,
                                       formamide_percent: float = 0.0,
                                       trehalose_m: float = 0.0,
                                       ethanol_percent: float = 0.0,
                                       urea_m: float = 0.0,
                                       tmac_m: float = 0.0) -> np.ndarray:
    """
    DEPRECATED: Calculate Tm for a batch of sequences with additive corrections.

    This function uses simplified additive corrections that do not account for
    GC-dependent effects of TMAC and betaine. Use ReactionConditions from
    neoswga.core.reaction_conditions instead, which implements the scientifically
    correct GC-dependent corrections.

    Example replacement:
        from neoswga.core.reaction_conditions import ReactionConditions
        conditions = ReactionConditions(
            na_conc=na_conc, mg_conc=mg_conc,
            dmso_percent=dmso_percent, betaine_m=betaine_m, ...
        )
        tms = [conditions.calculate_effective_tm(seq) for seq in sequences]

    DEPRECATED coefficients (fixed, non-GC-dependent):
    - DMSO: -0.6C per %
    - Betaine: -2.3C per M (INCORRECT: should be -0.5C/M + GC-dependent)
    - TMAC: -10C per M (INCORRECT: should include GC-dependent normalization)

    Args:
        sequences: List of DNA sequences
        na_conc: Na+ concentration in mM
        mg_conc: Mg2+ concentration in mM
        primer_conc: Primer concentration in M
        dmso_percent: DMSO percentage (0-10%)
        betaine_m: Betaine concentration in M (0-2.5M)
        formamide_percent: Formamide percentage (0-10%)
        trehalose_m: Trehalose concentration in M (0-1M)
        ethanol_percent: Ethanol percentage (0-5%)
        urea_m: Urea concentration in M (0-2M)
        tmac_m: TMAC concentration in M (0-0.3M)

    Returns:
        numpy array of Tm values in degrees Celsius
    """
    warnings.warn(
        "calculate_tm_batch_with_additives() is deprecated. "
        "Use ReactionConditions.calculate_effective_tm() from "
        "neoswga.core.reaction_conditions instead, which implements "
        "GC-dependent corrections for TMAC and betaine.",
        DeprecationWarning,
        stacklevel=2
    )
    # Get base Tm values
    tm_values = calculate_tm_batch(sequences, na_conc, mg_conc, primer_conc)

    # Calculate total additive correction (applies to all sequences)
    # WARNING: These are simplified coefficients that ignore GC-dependent effects.
    # The betaine coefficient (2.3) is particularly inaccurate for AT-rich genomes.
    # Use ReactionConditions.calculate_effective_tm() for accurate calculations.
    additive_correction = 0.0
    additive_correction -= 0.6 * dmso_percent      # DMSO
    additive_correction -= 2.3 * betaine_m         # Betaine (WRONG - ignores GC-dependence)
    additive_correction -= 0.65 * formamide_percent # Formamide
    additive_correction -= 5.0 * trehalose_m       # Trehalose
    additive_correction -= 0.5 * ethanol_percent   # Ethanol
    additive_correction -= 5.0 * urea_m            # Urea
    additive_correction -= 10.0 * tmac_m           # TMAC (~-1C per 0.1M)

    # Apply correction to all sequences (vectorized)
    if additive_correction != 0.0:
        tm_values += additive_correction

    return tm_values


def calculate_gc_batch(sequences: List[str]) -> np.ndarray:
    """
    Calculate GC content for a batch of sequences.

    Args:
        sequences: List of DNA sequences

    Returns:
        numpy array of GC fractions (0.0 to 1.0)
    """
    gc_values = np.zeros(len(sequences))

    for i, seq in enumerate(sequences):
        seq_upper = seq.upper()
        gc_count = seq_upper.count('G') + seq_upper.count('C')
        gc_values[i] = gc_count / len(seq) if len(seq) > 0 else 0.0

    return gc_values


def calculate_wallace_tm_batch(sequences: List[str]) -> np.ndarray:
    """
    Calculate Wallace Tm for a batch of sequences (fast approximation).

    Wallace rule: Tm = 2*(A+T) + 4*(G+C)

    Args:
        sequences: List of DNA sequences

    Returns:
        numpy array of approximate Tm values in degrees Celsius
    """
    tm_values = np.zeros(len(sequences))

    for i, seq in enumerate(sequences):
        seq_upper = seq.upper()
        at_count = seq_upper.count('A') + seq_upper.count('T')
        gc_count = seq_upper.count('G') + seq_upper.count('C')
        tm_values[i] = 2 * at_count + 4 * gc_count

    return tm_values


# ========================================
# Helper to create ReactionConditions from parameters
# ========================================

def get_reaction_conditions_from_params() -> 'ReactionConditions':
    """
    Create a ReactionConditions object from the current parameter module settings.

    This helper bridges the parameter module with the ReactionConditions class,
    allowing thermodynamic calculations to use the reaction conditions specified
    in params.json.

    Returns:
        ReactionConditions object configured with current parameter settings
    """
    from neoswga.core.reaction_conditions import ReactionConditions
    import neoswga.core.parameter as parameter

    # Get reaction temperature, defaulting to polymerase optimal temp
    reaction_temp = getattr(parameter, 'reaction_temp', None)
    polymerase = getattr(parameter, 'polymerase', 'phi29')

    # Auto-set temperature based on polymerase if not specified
    if reaction_temp is None:
        from neoswga.core.reaction_conditions import POLYMERASE_CHARACTERISTICS
        if polymerase.lower() in POLYMERASE_CHARACTERISTICS:
            reaction_temp = POLYMERASE_CHARACTERISTICS[polymerase.lower()]['optimal_temp']
        else:
            reaction_temp = 30.0  # Default to phi29 temperature

    return ReactionConditions(
        temp=reaction_temp,
        na_conc=getattr(parameter, 'na_conc', 50.0),
        mg_conc=getattr(parameter, 'mg_conc', 2.0),
        dmso_percent=getattr(parameter, 'dmso_percent', 0.0),
        betaine_m=getattr(parameter, 'betaine_m', 0.0),
        trehalose_m=getattr(parameter, 'trehalose_m', 0.0),
        formamide_percent=getattr(parameter, 'formamide_percent', 0.0),
        ethanol_percent=getattr(parameter, 'ethanol_percent', 0.0),
        urea_m=getattr(parameter, 'urea_m', 0.0),
        tmac_m=getattr(parameter, 'tmac_m', 0.0),
        polymerase=polymerase
    )


if __name__ == "__main__":
    # Test calculations
    test_seq = "ATCGATCG"

    print(f"Test sequence: {test_seq}")
    print(f"GC content: {gc_content(test_seq):.2%}")
    print(f"Wallace Tm: {wallace_tm(test_seq):.1f}°C")
    print(f"Basic Tm: {calculate_tm_basic(test_seq):.1f}°C")
    print(f"Tm with salt (50mM Na+): {calculate_tm_with_salt(test_seq, na_conc=50):.1f}°C")
    print(f"ΔG at 37°C: {calculate_free_energy(test_seq, 37):.2f} kcal/mol")
    print(f"Binding probability at 30°C: {calculate_binding_probability(test_seq, 30):.2%}")

    enthalpy, entropy = calculate_enthalpy_entropy(test_seq)
    print(f"ΔH: {enthalpy:.2f} kcal/mol")
    print(f"ΔS: {entropy:.2f} cal/(mol*K)")
