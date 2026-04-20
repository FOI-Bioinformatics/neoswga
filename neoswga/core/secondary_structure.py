"""
Secondary structure prediction for primers.

Implements thermodynamics-based prediction of:
- Heterodimers (primer-primer interactions)
- Homodimers (self-complementary structures)
- Hairpins (stem-loop structures)
- 3' end extension risk

Uses dynamic programming approach based on Zuker algorithm
with nearest-neighbor thermodynamic parameters.

References:
- Zuker (2003) Nucleic Acids Res 31:3406-3415
- Markham & Zuker (2008) Methods Mol Biol 453:3-31
"""

import numpy as np
from typing import List, Tuple, Dict, Optional
import warnings
from neoswga.core import thermodynamics as thermo
from neoswga.core import reaction_conditions as rc


# ========================================
# Loop and Bulge Penalties
# ========================================
# Loop initiation penalties from Serra & Turner (1995) Methods Enzymol
# 259:242-261, refined by Mathews et al. (1999) J Mol Biol 288:911-940.
# For loops larger than 6 nt, Jacobson-Stockmayer extrapolation is used:
#   G_loop(n) = G_loop(6) + 1.75 * RT * ln(n/6)
# (Jacobson & Stockmayer (1950) J Chem Phys 18:1600-1606).
#
# Bulge penalties follow Mathews et al. (1999) and Freier et al. (1986)
# PNAS 83:9373-9377 for single-nucleotide bulges.

def loop_penalty(size: int) -> float:
    """
    Calculate thermodynamic penalty for loops.

    Empirical values for 3-6 nt from Serra & Turner (1995) and
    Mathews et al. (1999). Larger loops use Jacobson-Stockmayer
    extrapolation.

    Args:
        size: Loop size in nucleotides

    Returns:
        Penalty in kcal/mol
    """
    if size < 3:
        return 1000.0  # Physically impossible
    elif size <= 6:
        # Small loops: empirical values (Serra & Turner 1995, Table 1)
        penalties = {3: 5.6, 4: 5.0, 5: 5.2, 6: 5.4}
        return penalties[size]
    else:
        # Larger loops: Jacobson-Stockmayer logarithmic extrapolation
        # G_loop(n) = G_loop(6) + 1.75 * R * T * ln(n/6)
        # R = 1.987 cal/(mol*K), T = 310.15 K (37C reference), convert to kcal/mol
        R_kcal = 1.987e-3  # kcal/(mol*K)
        T_ref = 310.15     # 37C in Kelvin (standard reference temperature)
        return 5.4 + 1.75 * R_kcal * T_ref * np.log(size / 6.0)


def bulge_penalty(size: int) -> float:
    """
    Calculate thermodynamic penalty for bulges.

    Single-nucleotide bulge value from Freier et al. (1986) PNAS 83:9373-9377.
    Multi-nucleotide bulges from Mathews et al. (1999) J Mol Biol 288:911-940.

    Args:
        size: Bulge size in nucleotides

    Returns:
        Penalty in kcal/mol
    """
    if size == 0:
        return 0.0
    elif size == 1:
        return 3.3  # Single nucleotide bulge (Freier 1986)
    elif size <= 6:
        # Small bulges (Mathews 1999)
        return 3.3 + 1.7 * (size - 1)
    else:
        # Larger bulges: linear extrapolation beyond 6 nt
        return 3.3 + 1.7 * 5 + 2.0 * (size - 6)


# ========================================
# Dynamic Programming for Structure Prediction
# ========================================

class StructurePrediction:
    """
    Predict secondary structures using dynamic programming.
    """

    def __init__(self, conditions: Optional[rc.ReactionConditions] = None):
        """
        Initialize structure predictor.

        Args:
            conditions: Reaction conditions for temperature-dependent calculations
        """
        self.conditions = conditions if conditions else rc.get_standard_conditions()

    def is_complementary(self, base1: str, base2: str) -> bool:
        """Check if two bases are Watson-Crick complementary."""
        pairs = {('A', 'T'), ('T', 'A'), ('G', 'C'), ('C', 'G')}
        return (base1.upper(), base2.upper()) in pairs

    def predict_heterodimer(self, seq1: str, seq2: str) -> Dict:
        """
        Predict heterodimer structure and stability.

        Args:
            seq1: First primer sequence
            seq2: Second primer sequence

        Returns:
            Dictionary with:
                - energy: Minimum free energy (kcal/mol)
                - tm: Estimated Tm of dimer
                - forms_dimer: Boolean if stable at reaction temp
                - binding_region: (start1, end1, start2, end2) of main binding
                - severity: Score 0-1 indicating dimer severity
        """
        n, m = len(seq1), len(seq2)

        # Initialize DP table for minimum free energy
        # E[i][j] = min energy for seq1[0:i] binding to seq2[0:j]
        E = np.full((n + 1, m + 1), np.inf)
        E[0][0] = 0.0

        # Traceback table
        parent = {}

        seq2_rc = thermo.reverse_complement(seq2)

        # Allow alignment to start anywhere in either strand (not just at 0, 0).
        # Without this, a dimer that binds in the middle of either primer is
        # unreachable in the DP: E[0][*] and E[*][0] remain inf, so options 2/3
        # cannot propagate bulges from those edges. Setting E[i][0]=E[0][j]=0
        # lets the DP "start" the alignment at any offset in either primer.
        for i in range(n + 1):
            E[i][0] = 0.0
        for j in range(m + 1):
            E[0][j] = 0.0

        # Fill DP table
        for i in range(n + 1):
            for j in range(m + 1):
                if i == 0 and j == 0:
                    continue

                current_min = E[i][j]

                # Option 1: Extend alignment with base pair. seq2_rc is already
                # the reverse complement of seq2, so a Watson-Crick pair between
                # seq1[i-1] and seq2[len-j] (antiparallel binding) becomes an
                # identity check against seq2_rc[j-1]. The earlier is_complementary
                # check was a bug: it required e.g. C-G (complementary) against
                # seq2_rc which already maps that G back to C, so the match
                # never fired and every pair reported severity=0.
                if i > 0 and j > 0:
                    if seq1[i-1].upper() == seq2_rc[j-1].upper():
                        # Check if we can form a stack
                        if i >= 2 and j >= 2:
                            # Get stacking energy
                            stack1 = seq1[i-2:i]
                            stack2 = seq2_rc[j-2:j]

                            # Calculate ΔG for this stack
                            try:
                                enthalpy, entropy = thermo.calculate_enthalpy_entropy(stack1 + stack2)
                                T = self.conditions.temp + 273.15
                                stack_energy = enthalpy - (T * entropy / 1000)
                            except (KeyError, ValueError, TypeError):
                                # Unknown dinucleotide stack, use average energy
                                stack_energy = -2.0  # Average stack energy

                            if E[i-1][j-1] + stack_energy < current_min:
                                current_min = E[i-1][j-1] + stack_energy
                                parent[(i, j)] = (i-1, j-1, 'pair')
                        else:
                            # First pair - initiation energy
                            if E[i-1][j-1] < current_min:
                                current_min = E[i-1][j-1]
                                parent[(i, j)] = (i-1, j-1, 'init')

                # Option 2: Bulge in seq1
                if i > 0:
                    for k in range(1, min(i, 7)):  # Max bulge size 6
                        penalty = bulge_penalty(k)
                        if E[i-k][j] + penalty < current_min:
                            current_min = E[i-k][j] + penalty
                            parent[(i, j)] = (i-k, j, f'bulge1_{k}')

                # Option 3: Bulge in seq2
                if j > 0:
                    for k in range(1, min(j, 7)):
                        penalty = bulge_penalty(k)
                        if E[i][j-k] + penalty < current_min:
                            current_min = E[i][j-k] + penalty
                            parent[(i, j)] = (i, j-k, f'bulge2_{k}')

                # Option 4: Internal loop
                if i > 0 and j > 0:
                    for k in range(1, min(i, j, 7)):
                        penalty = loop_penalty(k * 2)
                        if E[i-k][j-k] + penalty < current_min:
                            current_min = E[i-k][j-k] + penalty
                            parent[(i, j)] = (i-k, j-k, f'loop_{k}')

                E[i][j] = current_min

        # Find minimum energy
        min_energy = E[n][m]

        # Check if dimer forms
        # Dimer stable if ΔG < -9 kcal/mol at reaction temperature
        forms_dimer = min_energy < -9.0

        # Estimate Tm of dimer
        # Very rough estimate
        estimated_tm = thermo.energy_to_tm(
            min_energy,
            min(n, m) // 2,  # Approximate binding length
            self.conditions.na_conc,
            self.conditions.mg_conc
        )

        # Check if stable at reaction temperature
        stable_at_temp = estimated_tm > self.conditions.temp

        # Calculate severity score (0-1)
        # Based on: energy, Tm margin, and binding length
        if not forms_dimer:
            severity = 0.0
        else:
            # Energy component (0-1, worse at lower energy)
            energy_score = max(0, min(1, (-min_energy - 9) / 15))  # Scale -9 to -24

            # Tm margin component
            tm_margin = estimated_tm - self.conditions.temp
            tm_score = max(0, min(1, tm_margin / 20))  # Scale 0-20°C

            # Binding length estimate
            binding_len = self._estimate_binding_length(parent, n, m)
            length_score = min(1, binding_len / 10)  # Scale 0-10 bp

            # Combined severity
            severity = (energy_score * 0.4 + tm_score * 0.4 + length_score * 0.2)

        return {
            'energy': min_energy,
            'tm': estimated_tm,
            'forms_dimer': stable_at_temp,
            'severity': severity,
            'binding_length': self._estimate_binding_length(parent, n, m)
        }

    def _estimate_binding_length(self, parent: Dict, n: int, m: int) -> int:
        """
        Estimate number of bound base pairs from traceback.

        Args:
            parent: Traceback dictionary
            n: Length of seq1
            m: Length of seq2

        Returns:
            Number of base pairs
        """
        count = 0
        pos = (n, m)

        while pos in parent:
            prev_i, prev_j, move_type = parent[pos]
            if 'pair' in move_type or 'init' in move_type:
                count += 1
            pos = (prev_i, prev_j)

            # Safety check
            if count > max(n, m):
                break

        return count

    def check_3prime_extension_risk(self, seq1: str, seq2: str) -> float:
        """
        Assess risk that 3' ends will allow primer extension.

        3' end complementarity is most critical for PCR extension.

        Args:
            seq1: First primer sequence
            seq2: Second primer sequence

        Returns:
            Risk score 0-1 (1 = high risk of extension)
        """
        # Check last 5 bases (most critical for extension)
        terminal_len = 5
        term1 = seq1[-terminal_len:]
        term2 = thermo.reverse_complement(seq2[-terminal_len:])

        # Count complementary bases
        matches = sum(1 for a, b in zip(term1, term2) if self.is_complementary(a, b))

        # Check consecutive matches at very end (most critical)
        consecutive_3prime = 0
        for i in range(min(len(term1), len(term2))):
            if self.is_complementary(term1[-(i+1)], term2[-(i+1)]):
                consecutive_3prime += 1
            else:
                break

        # Calculate risk
        # High risk if: 3+ matches in last 5, or 2+ consecutive at 3' end
        risk = 0.0

        if matches >= 4:
            risk = 1.0  # Very high risk
        elif matches == 3:
            risk = 0.6
        elif matches == 2:
            risk = 0.3
        elif matches == 1:
            risk = 0.1

        # Boost risk if consecutive at 3' end
        if consecutive_3prime >= 3:
            risk = max(risk, 0.9)
        elif consecutive_3prime >= 2:
            risk = max(risk, 0.5)

        return risk

    def predict_hairpin(self, seq: str) -> List[Dict]:
        """
        Predict all possible hairpin structures.

        Args:
            seq: DNA sequence

        Returns:
            List of hairpins, each with:
                - position: Start of hairpin
                - loop_size: Size of loop (nt)
                - stem_length: Length of stem (bp)
                - energy: ΔG (kcal/mol)
                - tm: Estimated Tm
                - stable: Whether stable at reaction temp
        """
        hairpins = []
        n = len(seq)

        # Minimum hairpin: 3bp stem + 3nt loop
        min_stem = 3
        min_loop = 3
        max_loop = 8  # Typical hairpin loops

        for loop_start in range(min_stem, n - min_stem - min_loop):
            for loop_size in range(min_loop, max_loop + 1):
                if loop_start + loop_size >= n:
                    break

                loop_end = loop_start + loop_size

                # Try different stem lengths
                max_stem = min(loop_start, n - loop_end)

                for stem_len in range(min_stem, max_stem + 1):
                    # Check if stem regions are complementary
                    stem_5prime = seq[loop_start - stem_len:loop_start]
                    stem_3prime = seq[loop_end:loop_end + stem_len]

                    if self._is_stem(stem_5prime, stem_3prime):
                        # Calculate hairpin energy
                        energy = self._calculate_hairpin_energy(
                            stem_5prime, stem_3prime, loop_size
                        )

                        # Estimate Tm
                        tm = thermo.energy_to_tm(
                            energy, stem_len,
                            self.conditions.na_conc,
                            self.conditions.mg_conc
                        )

                        # Apply additive correction with GC-dependent effects
                        # Use stem GC content and length for hairpin Tm correction
                        stem_gc = thermo.gc_content(stem_5prime + stem_3prime)
                        tm_effective = self.conditions.adjust_tm(tm, stem_gc, stem_len)

                        # Check if stable
                        stable = tm_effective > self.conditions.temp

                        hairpins.append({
                            'position': loop_start - stem_len,
                            'loop_size': loop_size,
                            'stem_length': stem_len,
                            'energy': energy,
                            'tm': tm_effective,
                            'stable': stable
                        })

        # Sort by energy (most stable first)
        hairpins.sort(key=lambda x: x['energy'])

        return hairpins

    def _is_stem(self, seq5: str, seq3: str) -> bool:
        """
        Check if two sequences can form a stem.

        Requires at least 70% complementarity.

        Args:
            seq5: 5' side of stem
            seq3: 3' side of stem

        Returns:
            True if sequences can form stem
        """
        if len(seq5) != len(seq3):
            return False

        seq3_rc = thermo.reverse_complement(seq3)
        matches = sum(1 for a, b in zip(seq5, seq3_rc) if a == b)

        return matches / len(seq5) >= 0.7

    def _calculate_hairpin_energy(self, stem5: str, stem3: str, loop_size: int) -> float:
        """
        Calculate free energy of hairpin.

        Args:
            stem5: 5' side of stem
            stem3: 3' side of stem
            loop_size: Size of loop

        Returns:
            ΔG in kcal/mol
        """
        # Stem stability
        stem_seq = stem5 + thermo.reverse_complement(stem3)
        try:
            enthalpy, entropy = thermo.calculate_enthalpy_entropy(stem_seq)
            T = self.conditions.temp + 273.15
            stem_energy = enthalpy - (T * entropy / 1000)
        except (KeyError, ValueError, TypeError):
            # Fallback: rough estimate for unknown stacks
            stem_energy = -2.0 * len(stem5)  # ~2 kcal/mol per bp

        # Loop penalty
        loop_cost = loop_penalty(loop_size)

        return stem_energy + loop_cost


# ========================================
# High-Level Interface Functions
# ========================================

def check_heterodimer(seq1: str, seq2: str,
                      conditions: Optional[rc.ReactionConditions] = None) -> Dict:
    """
    Check if two primers form heterodimers.

    Args:
        seq1: First primer
        seq2: Second primer
        conditions: Reaction conditions

    Returns:
        Dimer prediction dictionary
    """
    predictor = StructurePrediction(conditions)
    return predictor.predict_heterodimer(seq1, seq2)


def check_homodimer(seq: str,
                    conditions: Optional[rc.ReactionConditions] = None) -> Dict:
    """
    Check if primer forms homodimers (self-complementary).

    Args:
        seq: Primer sequence
        conditions: Reaction conditions

    Returns:
        Dimer prediction dictionary
    """
    predictor = StructurePrediction(conditions)
    return predictor.predict_heterodimer(seq, seq)


def check_hairpins(seq: str,
                  conditions: Optional[rc.ReactionConditions] = None) -> List[Dict]:
    """
    Find all hairpin structures in sequence.

    Args:
        seq: DNA sequence
        conditions: Reaction conditions

    Returns:
        List of hairpin predictions
    """
    predictor = StructurePrediction(conditions)
    return predictor.predict_hairpin(seq)


def calculate_dimer_matrix(primers: List[str],
                          conditions: Optional[rc.ReactionConditions] = None) -> np.ndarray:
    """
    Calculate pairwise dimer severity matrix.

    Args:
        primers: List of primer sequences
        conditions: Reaction conditions

    Returns:
        n x n matrix of dimer severity scores (0-1)
    """
    n = len(primers)
    matrix = np.zeros((n, n))

    predictor = StructurePrediction(conditions)

    for i in range(n):
        for j in range(i, n):
            if i == j:
                # Self-dimer (homodimer)
                result = predictor.predict_heterodimer(primers[i], primers[i])
            else:
                # Heterodimer
                result = predictor.predict_heterodimer(primers[i], primers[j])

            severity = result['severity']
            matrix[i][j] = severity
            matrix[j][i] = severity

    return matrix


def filter_primers_by_structure(primers: List[str],
                                max_hairpin_tm: float = 35.0,
                                max_self_dimer_severity: float = 0.3,
                                conditions: Optional[rc.ReactionConditions] = None) -> List[str]:
    """
    Filter primers that have problematic secondary structures.

    Args:
        primers: List of primer sequences
        max_hairpin_tm: Maximum allowed hairpin Tm (°C)
        max_self_dimer_severity: Maximum allowed self-dimer severity (0-1)
        conditions: Reaction conditions

    Returns:
        Filtered list of primers
    """
    predictor = StructurePrediction(conditions)
    filtered = []

    for primer in primers:
        # Check hairpins
        hairpins = predictor.predict_hairpin(primer)
        if hairpins:
            max_hp_tm = max(h['tm'] for h in hairpins)
            if max_hp_tm > max_hairpin_tm:
                continue  # Skip primers with stable hairpins

        # Check self-dimers
        self_dimer = predictor.predict_heterodimer(primer, primer)
        if self_dimer['severity'] > max_self_dimer_severity:
            continue  # Skip primers with strong self-dimers

        filtered.append(primer)

    return filtered


if __name__ == "__main__":
    # Test secondary structure prediction
    print("=== Secondary Structure Prediction Tests ===\n")

    # Test hairpin
    seq_hairpin = "GCGATCGCAAAAGCGATCGC"  # Palindromic ends
    print(f"Sequence: {seq_hairpin}")

    hairpins = check_hairpins(seq_hairpin)
    print(f"Found {len(hairpins)} hairpins")
    if hairpins:
        best = hairpins[0]
        print(f"Best hairpin: stem={best['stem_length']}bp, loop={best['loop_size']}nt")
        print(f"  ΔG={best['energy']:.2f} kcal/mol, Tm={best['tm']:.1f}°C, stable={best['stable']}")

    print()

    # Test heterodimer
    primer1 = "ATCGATCGAT"
    primer2 = "ATCGATCGAT"
    print(f"Primer 1: {primer1}")
    print(f"Primer 2: {primer2}")

    dimer = check_heterodimer(primer1, primer2)
    print(f"Heterodimer: ΔG={dimer['energy']:.2f} kcal/mol")
    print(f"  Tm={dimer['tm']:.1f}°C, forms_dimer={dimer['forms_dimer']}")
    print(f"  Severity={dimer['severity']:.2f}, binding_length={dimer['binding_length']}bp")

    print()

    # Test 3' extension risk
    predictor = StructurePrediction()
    risk = predictor.check_3prime_extension_risk(primer1, primer2)
    print(f"3' extension risk: {risk:.2f}")
