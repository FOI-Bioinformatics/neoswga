"""
Advanced feature engineering for machine learning-based primer prediction.

Engineers 95+ sophisticated features across multiple categories:
- Base composition and sequence statistics (20 features)
- Thermodynamic properties (15 features)
- Secondary structure propensity (12 features)
- Positional and spatial features (18 features)
- Sequence complexity and information theory (10 features)
- K-mer motif frequencies (20 features)
- Binding landscape features (15 features)
- Interaction and context features (10 features)

Total: 120+ features for random forest and deep learning models.
"""

import numpy as np
import pandas as pd
from typing import List, Dict, Tuple, Optional
from collections import Counter
import warnings
from scipy.stats import entropy, skew, kurtosis
from scipy.signal import find_peaks
import gzip

from neoswga.core import thermodynamics as thermo
from neoswga.core import reaction_conditions as rc
from neoswga.core import secondary_structure as ss


class AdvancedFeatureEngineer:
    """
    Comprehensive feature engineering for primer prediction.

    Generates 120+ features per primer for ML models.
    """

    def __init__(self,
                 genome_sequence: str,
                 conditions: rc.ReactionConditions,
                 primer_positions: Optional[Dict] = None):
        """
        Initialize feature engineer.

        Args:
            genome_sequence: Reference genome sequence
            conditions: Reaction conditions
            primer_positions: Optional pre-loaded positions
        """
        self.genome = genome_sequence.upper()
        self.genome_length = len(genome_sequence)
        self.conditions = conditions
        self.primer_positions = primer_positions or {}

        # Precompute genome statistics for efficiency
        self._precompute_genome_stats()

    def _precompute_genome_stats(self):
        """Precompute genome-level statistics."""
        self.genome_gc = thermo.gc_content(self.genome)
        self.genome_dinuc_freqs = self._calculate_dinucleotide_frequencies(self.genome)

        # K-mer backgrounds (2-4mers)
        self.genome_2mer_freqs = self._kmer_frequencies(self.genome, 2)
        self.genome_3mer_freqs = self._kmer_frequencies(self.genome, 3)
        self.genome_4mer_freqs = self._kmer_frequencies(self.genome, 4)

    def engineer_features(self, primers: List[str]) -> pd.DataFrame:
        """
        Engineer all features for list of primers.

        Args:
            primers: List of primer sequences

        Returns:
            DataFrame with one row per primer, 120+ feature columns
        """
        features_list = []

        for primer in primers:
            features = {}

            # Category 1: Base composition (20 features)
            features.update(self._base_composition_features(primer))

            # Category 2: Thermodynamic properties (15 features)
            features.update(self._thermodynamic_features(primer))

            # Category 3: Secondary structure (12 features)
            features.update(self._secondary_structure_features(primer))

            # Category 4: Positional and spatial (18 features)
            features.update(self._positional_features(primer))

            # Category 5: Sequence complexity (10 features)
            features.update(self._complexity_features(primer))

            # Category 6: K-mer motifs (20 features)
            features.update(self._motif_features(primer))

            # Category 7: Binding landscape (15 features)
            features.update(self._binding_landscape_features(primer))

            # Category 8: Context and interactions (10 features)
            features.update(self._context_features(primer))

            features['primer'] = primer
            features_list.append(features)

        df = pd.DataFrame(features_list)

        # Reorder columns: primer first, then features
        cols = ['primer'] + [c for c in df.columns if c != 'primer']
        return df[cols]

    # ========================================
    # Category 1: Base Composition Features
    # ========================================

    def _base_composition_features(self, primer: str) -> Dict:
        """Base composition and nucleotide statistics."""
        n = len(primer)

        # Individual base counts and frequencies
        a_count = primer.count('A')
        t_count = primer.count('T')
        g_count = primer.count('G')
        c_count = primer.count('C')

        gc_count = g_count + c_count
        at_count = a_count + t_count

        # Purine/pyrimidine
        purine = a_count + g_count  # A, G
        pyrimidine = t_count + c_count  # T, C

        # Strong/weak bonds
        strong = gc_count  # 3 H-bonds
        weak = at_count  # 2 H-bonds

        # Longest homopolymer runs
        longest_a = self._longest_run(primer, 'A')
        longest_t = self._longest_run(primer, 'T')
        longest_g = self._longest_run(primer, 'G')
        longest_c = self._longest_run(primer, 'C')
        longest_any = max(longest_a, longest_t, longest_g, longest_c)

        # GC content in windows
        gc_5prime = (primer[:5].count('G') + primer[:5].count('C')) / 5 if n >= 5 else 0
        gc_3prime = (primer[-5:].count('G') + primer[-5:].count('C')) / 5 if n >= 5 else 0
        gc_middle = (primer[n//4:3*n//4].count('G') + primer[n//4:3*n//4].count('C')) / (n//2) if n >= 8 else 0

        # GC skew
        gc_skew = (g_count - c_count) / (g_count + c_count) if (g_count + c_count) > 0 else 0
        at_skew = (a_count - t_count) / (a_count + t_count) if (a_count + t_count) > 0 else 0

        return {
            'length': n,
            'a_count': a_count,
            't_count': t_count,
            'g_count': g_count,
            'c_count': c_count,
            'gc_content': gc_count / n,
            'at_content': at_count / n,
            'purine_content': purine / n,
            'pyrimidine_content': pyrimidine / n,
            'strong_bonds_content': strong / n,
            'weak_bonds_content': weak / n,
            'longest_homopolymer': longest_any,
            'longest_a_run': longest_a,
            'longest_g_run': longest_g,
            'gc_content_5prime': gc_5prime,
            'gc_content_3prime': gc_3prime,
            'gc_content_middle': gc_middle,
            'gc_skew': gc_skew,
            'at_skew': at_skew,
            'gc_clamp': self._count_gc_clamp(primer),
            'purine_pyrimidine_ratio': purine / pyrimidine if pyrimidine > 0 else np.inf
        }

    def _longest_run(self, seq: str, base: str) -> int:
        """Find longest run of a specific base."""
        max_run = 0
        current_run = 0
        for b in seq:
            if b == base:
                current_run += 1
                max_run = max(max_run, current_run)
            else:
                current_run = 0
        return max_run

    def _count_gc_clamp(self, primer: str) -> int:
        """Count G/C bases in last 5 positions."""
        if len(primer) < 5:
            return sum(1 for b in primer if b in 'GC')
        return sum(1 for b in primer[-5:] if b in 'GC')

    # ========================================
    # Category 2: Thermodynamic Features
    # ========================================

    def _thermodynamic_features(self, primer: str) -> Dict:
        """Thermodynamic stability and binding features."""
        # Basic thermodynamics
        tm_base = thermo.calculate_tm_with_salt(primer, self.conditions.na_conc, self.conditions.mg_conc)
        tm_effective = self.conditions.calculate_effective_tm(primer)
        dg_37 = thermo.calculate_free_energy(primer, 37.0)
        dg_reaction = thermo.calculate_free_energy(primer, self.conditions.temp)

        enthalpy, entropy = thermo.calculate_enthalpy_entropy(primer)

        # Binding probability at different temperatures
        prob_reaction = thermo.calculate_binding_probability(primer, self.conditions.temp)
        prob_room = thermo.calculate_binding_probability(primer, 25.0)
        prob_cold = thermo.calculate_binding_probability(primer, 4.0)

        # Tm margin from reaction temperature
        tm_margin = tm_effective - self.conditions.temp

        # Temperature sensitivity (derivative)
        dg_25 = thermo.calculate_free_energy(primer, 25.0)
        dg_45 = thermo.calculate_free_energy(primer, 45.0)
        tm_sensitivity = (dg_45 - dg_25) / 20.0  # ΔΔG per °C

        # Polymerase compatibility
        polymerase_range = self.conditions.get_polymerase_range()
        in_polymerase_range = 1 if polymerase_range[0] <= tm_effective <= polymerase_range[1] else 0

        return {
            'tm_base': tm_base,
            'tm_effective': tm_effective,
            'tm_margin': tm_margin,
            'dg_37c': dg_37,
            'dg_reaction_temp': dg_reaction,
            'enthalpy': enthalpy,
            'entropy': entropy,
            'binding_prob_reaction': prob_reaction,
            'binding_prob_room': prob_room,
            'binding_prob_cold': prob_cold,
            'tm_sensitivity': tm_sensitivity,
            'in_polymerase_range': in_polymerase_range,
            'enthalpy_per_bp': enthalpy / len(primer),
            'entropy_per_bp': entropy / len(primer),
            'dg_per_bp': dg_reaction / len(primer)
        }

    # ========================================
    # Category 3: Secondary Structure Features
    # ========================================

    def _secondary_structure_features(self, primer: str) -> Dict:
        """Secondary structure propensity features."""
        # Hairpins
        hairpins = ss.check_hairpins(primer, self.conditions)

        if hairpins:
            max_hairpin_dg = min(h['energy'] for h in hairpins)
            max_hairpin_tm = max(h['tm'] for h in hairpins)
            max_hairpin_stem = max(h['stem_length'] for h in hairpins)
            num_hairpins = len(hairpins)
            hairpin_stable = 1 if max_hairpin_tm > self.conditions.temp + 5 else 0
        else:
            max_hairpin_dg = 0.0
            max_hairpin_tm = 0.0
            max_hairpin_stem = 0
            num_hairpins = 0
            hairpin_stable = 0

        # Self-dimer
        self_dimer = ss.check_homodimer(primer, self.conditions)

        # 3' end structure propensity
        terminal_3prime = primer[-5:] if len(primer) >= 5 else primer
        terminal_3prime_gc = thermo.gc_content(terminal_3prime)

        # 5' end structure propensity
        terminal_5prime = primer[:5] if len(primer) >= 5 else primer
        terminal_5prime_gc = thermo.gc_content(terminal_5prime)

        return {
            'num_hairpins': num_hairpins,
            'max_hairpin_dg': max_hairpin_dg,
            'max_hairpin_tm': max_hairpin_tm,
            'max_hairpin_stem_length': max_hairpin_stem,
            'hairpin_stable': hairpin_stable,
            'self_dimer_severity': self_dimer['severity'],
            'self_dimer_dg': self_dimer['energy'],
            'self_dimer_tm': self_dimer['tm'],
            'terminal_3prime_gc': terminal_3prime_gc,
            'terminal_5prime_gc': terminal_5prime_gc,
            'terminal_gc_asymmetry': abs(terminal_3prime_gc - terminal_5prime_gc),
            'palindrome_score': self._palindrome_score(primer)
        }

    def _palindrome_score(self, primer: str) -> float:
        """Calculate how palindromic the sequence is (0-1)."""
        rc_primer = thermo.reverse_complement(primer)
        matches = sum(1 for a, b in zip(primer, rc_primer) if a == b)
        return matches / len(primer)

    # ========================================
    # Category 4: Positional and Spatial Features
    # ========================================

    def _positional_features(self, primer: str) -> Dict:
        """Positional and spatial binding features."""
        positions = self.primer_positions.get(primer, {'forward': [], 'reverse': []})

        fwd_pos = positions.get('forward', [])
        rev_pos = positions.get('reverse', [])

        total_sites = len(fwd_pos) + len(rev_pos)

        if total_sites == 0:
            return self._empty_positional_features()

        # Binding frequency
        binding_freq = total_sites / self.genome_length

        # Strand bias
        strand_bias = abs(len(fwd_pos) - len(rev_pos)) / total_sites

        # Gap statistics (forward strand)
        if len(fwd_pos) > 1:
            fwd_gaps = np.diff(sorted(fwd_pos))
            mean_fwd_gap = np.mean(fwd_gaps)
            std_fwd_gap = np.std(fwd_gaps)
            min_fwd_gap = np.min(fwd_gaps)
            max_fwd_gap = np.max(fwd_gaps)
            median_fwd_gap = np.median(fwd_gaps)

            # Gini index
            gini_fwd = self._calculate_gini(fwd_gaps)

            # Clustering coefficient
            clustering_fwd = self._calculate_clustering(fwd_gaps, mean_fwd_gap)
        else:
            mean_fwd_gap = std_fwd_gap = min_fwd_gap = max_fwd_gap = median_fwd_gap = 0
            gini_fwd = 0
            clustering_fwd = 0

        # Nearest opposite-strand distance
        if len(fwd_pos) > 0 and len(rev_pos) > 0:
            nearest_opposite = self._mean_nearest_opposite_distance(fwd_pos, rev_pos)
        else:
            nearest_opposite = self.genome_length

        # Positional entropy (how spread out)
        position_entropy = self._positional_entropy(fwd_pos + rev_pos)

        # Coverage (what fraction of genome within max_gap of a site)
        coverage_fraction = self._estimate_coverage(fwd_pos + rev_pos, max_gap=3000)

        return {
            'total_binding_sites': total_sites,
            'binding_frequency': binding_freq,
            'strand_bias': strand_bias,
            'mean_forward_gap': mean_fwd_gap,
            'std_forward_gap': std_fwd_gap,
            'min_forward_gap': min_fwd_gap,
            'max_forward_gap': max_fwd_gap,
            'median_forward_gap': median_fwd_gap,
            'gini_index': gini_fwd,
            'clustering_coefficient': clustering_fwd,
            'mean_nearest_opposite_distance': nearest_opposite,
            'positional_entropy': position_entropy,
            'coverage_fraction': coverage_fraction,
            'sites_per_100kb': (total_sites / self.genome_length) * 100000,
            'forward_site_count': len(fwd_pos),
            'reverse_site_count': len(rev_pos),
            'gap_cv': std_fwd_gap / mean_fwd_gap if mean_fwd_gap > 0 else 0,
            'gap_skewness': self._safe_skew(fwd_gaps) if len(fwd_pos) > 2 else 0
        }

    @staticmethod
    def _safe_skew(values):
        """Compute skewness, returning 0 for degenerate data."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            result = skew(values)
        return float(np.nan_to_num(result))

    def _empty_positional_features(self) -> Dict:
        """Return zero-filled positional features."""
        return {
            'total_binding_sites': 0, 'binding_frequency': 0, 'strand_bias': 0,
            'mean_forward_gap': 0, 'std_forward_gap': 0, 'min_forward_gap': 0,
            'max_forward_gap': 0, 'median_forward_gap': 0, 'gini_index': 0,
            'clustering_coefficient': 0, 'mean_nearest_opposite_distance': 0,
            'positional_entropy': 0, 'coverage_fraction': 0, 'sites_per_100kb': 0,
            'forward_site_count': 0, 'reverse_site_count': 0, 'gap_cv': 0,
            'gap_skewness': 0
        }

    def _calculate_gini(self, gaps: np.ndarray) -> float:
        """Calculate Gini coefficient of gap distribution."""
        if len(gaps) == 0:
            return 0.0
        sorted_gaps = np.sort(gaps)
        n = len(gaps)
        index = np.arange(1, n + 1)
        return ((np.sum((2 * index - n - 1) * sorted_gaps)) / (n * np.sum(sorted_gaps)))

    def _calculate_clustering(self, gaps: np.ndarray, mean_gap: float) -> float:
        """Calculate clustering coefficient (fraction of small gaps)."""
        if len(gaps) == 0 or mean_gap == 0:
            return 0.0
        small_gaps = np.sum(gaps < mean_gap / 2)
        return small_gaps / len(gaps)

    def _mean_nearest_opposite_distance(self, fwd_pos: List[int], rev_pos: List[int]) -> float:
        """Calculate mean distance from forward to nearest reverse site."""
        if not fwd_pos or not rev_pos:
            return self.genome_length

        distances = []
        for fp in fwd_pos:
            min_dist = min(abs(fp - rp) for rp in rev_pos)
            distances.append(min_dist)

        return np.mean(distances)

    def _positional_entropy(self, positions: List[int]) -> float:
        """Calculate entropy of positional distribution."""
        if len(positions) < 2:
            return 0.0

        # Bin genome into windows
        bins = 100
        hist, _ = np.histogram(positions, bins=bins, range=(0, self.genome_length))
        hist = hist / np.sum(hist)  # Normalize

        return entropy(hist + 1e-10)  # Add small constant to avoid log(0)

    def _estimate_coverage(self, positions: List[int], max_gap: int = 3000) -> float:
        """Estimate coverage fraction within max_gap of binding sites."""
        if not positions:
            return 0.0

        covered = set()
        for pos in positions:
            start = max(0, pos - max_gap)
            end = min(self.genome_length, pos + max_gap)
            covered.update(range(start, end))

        return len(covered) / self.genome_length

    # ========================================
    # Category 5: Sequence Complexity Features
    # ========================================

    def _complexity_features(self, primer: str) -> Dict:
        """Information-theoretic complexity features."""
        # Shannon entropy
        shannon_ent = self._shannon_entropy(primer)

        # Linguistic complexity
        ling_complex = self._linguistic_complexity(primer)

        # Compressibility (gzip compression ratio)
        compressibility = self._compressibility(primer)

        # Dinucleotide diversity
        dinuc_diversity = self._dinucleotide_diversity(primer)

        # Repeat content
        tandem_repeats = self._count_tandem_repeats(primer)
        dinuc_repeats = self._count_dinucleotide_repeats(primer)

        # Uniqueness vs genome
        genome_similarity = self._genome_similarity(primer)

        return {
            'shannon_entropy': shannon_ent,
            'linguistic_complexity': ling_complex,
            'compressibility': compressibility,
            'dinucleotide_diversity': dinuc_diversity,
            'tandem_repeat_count': tandem_repeats,
            'dinucleotide_repeat_count': dinuc_repeats,
            'genome_similarity': genome_similarity,
            'normalized_entropy': shannon_ent / np.log2(4),  # Normalized to [0,1]
            'complexity_score': (shannon_ent + ling_complex) / 2,
            'uniqueness_score': 1 - genome_similarity
        }

    def _shannon_entropy(self, seq: str) -> float:
        """Calculate Shannon entropy of sequence."""
        counts = Counter(seq)
        probs = np.array([counts[b] / len(seq) for b in 'ATGC'])
        return -np.sum(probs * np.log2(probs + 1e-10))

    def _linguistic_complexity(self, seq: str) -> float:
        """Calculate linguistic complexity (vocabulary/possible)."""
        n = len(seq)
        substrings = set()
        for i in range(n):
            for j in range(i + 1, n + 1):
                substrings.add(seq[i:j])

        max_possible = n * (n + 1) / 2
        return len(substrings) / max_possible

    def _compressibility(self, seq: str) -> float:
        """Calculate compression ratio (1 = not compressible)."""
        compressed = gzip.compress(seq.encode())
        return len(compressed) / len(seq)

    def _dinucleotide_diversity(self, seq: str) -> float:
        """Calculate diversity of dinucleotides (0-1)."""
        if len(seq) < 2:
            return 0.0
        dinucs = [seq[i:i+2] for i in range(len(seq) - 1)]
        return len(set(dinucs)) / 16  # 16 possible dinucleotides

    def _count_tandem_repeats(self, seq: str) -> int:
        """Count tandem repeats (AA, TT, GG, CC, etc.)."""
        count = 0
        for i in range(len(seq) - 1):
            if seq[i] == seq[i + 1]:
                count += 1
        return count

    def _count_dinucleotide_repeats(self, seq: str, max_count: int = 10) -> int:
        """Count long dinucleotide repeats (ATATAT...)."""
        if len(seq) < 4:
            return 0

        max_repeat = 0
        for i in range(len(seq) - 3):
            dinuc = seq[i:i+2]
            repeat_len = 2
            pos = i + 2
            while pos + 1 < len(seq) and seq[pos:pos+2] == dinuc:
                repeat_len += 2
                pos += 2
            max_repeat = max(max_repeat, repeat_len)

        return min(max_repeat, max_count)

    def _genome_similarity(self, primer: str) -> float:
        """Calculate how similar primer is to genome background."""
        # Compare dinucleotide frequencies
        primer_dinucs = self._calculate_dinucleotide_frequencies(primer)

        # Euclidean distance
        dist = np.sqrt(sum((primer_dinucs.get(dn, 0) - self.genome_dinuc_freqs.get(dn, 0))**2
                          for dn in set(primer_dinucs.keys()) | set(self.genome_dinuc_freqs.keys())))

        # Normalize to [0, 1] (0 = very different, 1 = very similar)
        return np.exp(-dist)

    def _calculate_dinucleotide_frequencies(self, seq: str) -> Dict[str, float]:
        """Calculate normalized dinucleotide frequencies."""
        if len(seq) < 2:
            return {}

        dinucs = [seq[i:i+2] for i in range(len(seq) - 1)]
        counts = Counter(dinucs)
        total = sum(counts.values())

        return {dn: count / total for dn, count in counts.items()}

    # ========================================
    # Category 6: K-mer Motif Features
    # ========================================

    def _motif_features(self, primer: str) -> Dict:
        """K-mer motif frequency features."""
        # 2-mer frequencies
        twoer_freqs = self._kmer_frequencies(primer, 2)

        # 3-mer frequencies for key motifs
        threemer_freqs = self._kmer_frequencies(primer, 3)

        # 4-mer frequencies for key motifs
        fourmer_freqs = self._kmer_frequencies(primer, 4)

        # Important motifs
        features = {}

        # 2-mers
        for kmer in ['AA', 'TT', 'GG', 'CC', 'AT', 'TA', 'GC', 'CG']:
            features[f'freq_2mer_{kmer}'] = twoer_freqs.get(kmer, 0)

        # 3-mers (select important ones)
        for kmer in ['AAA', 'TTT', 'GGG', 'CCC', 'ATC', 'TAG', 'GCA', 'CGT']:
            features[f'freq_3mer_{kmer}'] = threemer_freqs.get(kmer, 0)

        # 4-mers (select important ones)
        for kmer in ['AAAA', 'TTTT', 'ATCG', 'CGAT']:
            features[f'freq_4mer_{kmer}'] = fourmer_freqs.get(kmer, 0)

        return features

    def _kmer_frequencies(self, seq: str, k: int) -> Dict[str, float]:
        """Calculate k-mer frequencies."""
        if len(seq) < k:
            return {}

        kmers = [seq[i:i+k] for i in range(len(seq) - k + 1)]
        counts = Counter(kmers)
        total = sum(counts.values())

        return {kmer: count / total for kmer, count in counts.items()}

    # ========================================
    # Category 7: Binding Landscape Features
    # ========================================

    def _binding_landscape_features(self, primer: str) -> Dict:
        """Binding affinity landscape features."""
        # Calculate ΔG against all possible k-mers in genome
        k = len(primer)

        # Sample genome k-mers (too expensive to do all)
        sample_size = min(10000, self.genome_length - k + 1)
        sample_positions = np.random.choice(self.genome_length - k + 1, sample_size, replace=False)

        dg_values = []
        for pos in sample_positions:
            target = self.genome[pos:pos+k]
            if len(target) == k:
                try:
                    dg = thermo.calculate_free_energy(primer + target, self.conditions.temp)
                    dg_values.append(dg)
                except (KeyError, ValueError, TypeError):
                    # Skip positions with invalid sequences
                    pass

        if len(dg_values) == 0:
            return self._empty_binding_landscape_features()

        dg_array = np.array(dg_values)

        # Histogram-based features (binned ΔG distribution)
        bins = [-20, -15, -12, -10, -8, -6, -4, -2, 0, 2]
        hist, _ = np.histogram(dg_array, bins=bins)
        hist_total = np.sum(hist)
        hist_norm = hist / hist_total if hist_total > 0 else np.zeros_like(hist, dtype=float)

        features = {}
        for i, (low, high) in enumerate(zip(bins[:-1], bins[1:])):
            features[f'binding_hist_{low}to{high}'] = hist_norm[i]

        # Statistical features
        features['binding_mean_dg'] = np.mean(dg_array)
        features['binding_std_dg'] = np.std(dg_array)
        features['binding_min_dg'] = np.min(dg_array)
        features['binding_median_dg'] = np.median(dg_array)
        features['binding_strong_fraction'] = np.sum(dg_array < -10) / len(dg_array)
        features['binding_weak_fraction'] = np.sum(dg_array > -5) / len(dg_array)

        return features

    def _empty_binding_landscape_features(self) -> Dict:
        """Return zero-filled binding landscape features."""
        features = {}
        bins = [-20, -15, -12, -10, -8, -6, -4, -2, 0, 2]
        for i in range(len(bins) - 1):
            features[f'binding_hist_{bins[i]}to{bins[i+1]}'] = 0
        features.update({
            'binding_mean_dg': 0, 'binding_std_dg': 0, 'binding_min_dg': 0,
            'binding_median_dg': 0, 'binding_strong_fraction': 0, 'binding_weak_fraction': 0
        })
        return features

    # ========================================
    # Category 8: Context and Interaction Features
    # ========================================

    def _context_features(self, primer: str) -> Dict:
        """Contextual and interaction features."""
        # Flanking GC content (if positions known)
        positions = self.primer_positions.get(primer, {'forward': [], 'reverse': []})

        if positions.get('forward'):
            # Sample first position
            pos = positions['forward'][0]
            flanking_gc = self._flanking_gc_content(pos, window=50)
        else:
            flanking_gc = self.genome_gc

        # GC deviation from genome
        gc_deviation = abs(thermo.gc_content(primer) - self.genome_gc)

        # Length relative to genome average k-mer
        length_ratio = len(primer) / 10  # Assuming 10bp average

        # Tm relative to genome average
        genome_sample_tm = self._sample_genome_tm(k=len(primer))
        tm_deviation = self.conditions.calculate_effective_tm(primer) - genome_sample_tm

        return {
            'flanking_gc_content': flanking_gc,
            'gc_deviation_from_genome': gc_deviation,
            'length_ratio': length_ratio,
            'tm_deviation_from_genome': tm_deviation,
            'gc_content_rank': self._gc_content_rank(primer),
            'tm_rank': self._tm_rank(primer),
            'complexity_rank': self._complexity_rank(primer),
            'binding_strength_rank': self._binding_strength_rank(primer),
            'overall_quality_score': self._overall_quality_score(primer),
            'interaction_risk_score': self._interaction_risk_score(primer)
        }

    def _flanking_gc_content(self, position: int, window: int = 50) -> float:
        """Calculate GC content in flanking region."""
        start = max(0, position - window)
        end = min(self.genome_length, position + window)

        flanking = self.genome[start:end]
        return thermo.gc_content(flanking)

    def _sample_genome_tm(self, k: int, n_samples: int = 100) -> float:
        """Sample average Tm of k-mers from genome."""
        tms = []
        for _ in range(n_samples):
            pos = np.random.randint(0, self.genome_length - k + 1)
            kmer = self.genome[pos:pos+k]
            if len(kmer) == k and 'N' not in kmer:
                tm = self.conditions.calculate_effective_tm(kmer)
                tms.append(tm)

        return np.mean(tms) if tms else 37.0

    def _gc_content_rank(self, primer: str) -> float:
        """Rank GC content (0=low, 1=high)."""
        gc = thermo.gc_content(primer)
        return min(1.0, max(0.0, (gc - 0.2) / 0.6))  # Map 20-80% to 0-1

    def _tm_rank(self, primer: str) -> float:
        """Rank Tm optimality (0=poor, 1=optimal)."""
        tm = self.conditions.calculate_effective_tm(primer)
        tm_range = self.conditions.get_polymerase_range()
        midpoint = np.mean(tm_range)

        # Gaussian around midpoint
        return np.exp(-((tm - midpoint) / 10)**2)

    def _complexity_rank(self, primer: str) -> float:
        """Rank sequence complexity (0=low, 1=high)."""
        return self._linguistic_complexity(primer)

    def _binding_strength_rank(self, primer: str) -> float:
        """Rank binding strength (0=weak, 1=strong)."""
        dg = thermo.calculate_free_energy(primer, self.conditions.temp)
        # Strong binding: ΔG < -15, Weak: ΔG > -5
        return min(1.0, max(0.0, (-dg - 5) / 10))

    def _overall_quality_score(self, primer: str) -> float:
        """Composite quality score."""
        gc_score = 1 - abs(thermo.gc_content(primer) - 0.5) / 0.5
        tm_score = self._tm_rank(primer)
        complexity_score = self._complexity_rank(primer)

        return (gc_score * 0.3 + tm_score * 0.4 + complexity_score * 0.3)

    def _interaction_risk_score(self, primer: str) -> float:
        """Risk of unwanted interactions."""
        self_dimer = ss.check_homodimer(primer, self.conditions)
        hairpins = ss.check_hairpins(primer, self.conditions)

        dimer_risk = self_dimer['severity']
        hairpin_risk = 1 if hairpins and max(h['tm'] for h in hairpins) > self.conditions.temp + 5 else 0

        return (dimer_risk * 0.6 + hairpin_risk * 0.4)


def engineer_features_for_primers(primers: List[str],
                                  genome_sequence: str,
                                  conditions: rc.ReactionConditions,
                                  primer_positions: Optional[Dict] = None) -> pd.DataFrame:
    """
    High-level function to engineer features for primer list.

    Args:
        primers: List of primer sequences
        genome_sequence: Reference genome
        conditions: Reaction conditions
        primer_positions: Optional pre-loaded positions

    Returns:
        DataFrame with 120+ features per primer
    """
    engineer = AdvancedFeatureEngineer(genome_sequence, conditions, primer_positions)
    return engineer.engineer_features(primers)


if __name__ == "__main__":
    print("Advanced Feature Engineering Module")
    print("=" * 60)
    print("\n120+ sophisticated features for machine learning:")
    print("  1. Base composition (20 features)")
    print("  2. Thermodynamics (15 features)")
    print("  3. Secondary structure (12 features)")
    print("  4. Positional/spatial (18 features)")
    print("  5. Sequence complexity (10 features)")
    print("  6. K-mer motifs (20 features)")
    print("  7. Binding landscape (15 features)")
    print("  8. Context/interactions (10 features)")
    print("\nExample usage:")
    print("""
    from neoswga.core import advanced_features, reaction_conditions

    conditions = reaction_conditions.get_enhanced_conditions()

    features_df = advanced_features.engineer_features_for_primers(
        primers=['ATCGATCGAT', 'GCTAGCTAGC'],
        genome_sequence=genome_seq,
        conditions=conditions,
        primer_positions=positions_dict
    )

    print(f"Features shape: {features_df.shape}")
    # (2, 121) - 2 primers × 121 features (120 + primer sequence)
    """)
