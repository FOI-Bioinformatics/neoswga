"""
Adaptive primer length search with thermodynamic optimization.

Intelligently determines optimal k-mer length (6-15) based on:
- Reaction conditions (additives, polymerase, temperature)
- Genome complexity and composition
- Target specificity requirements
- Coverage goals

Uses advanced filtering with thermodynamic models and secondary structure prediction.
"""

import numpy as np
import pandas as pd
from typing import List, Dict, Tuple, Optional, Set
from collections import defaultdict
import warnings
from dataclasses import dataclass

import neoswga.core.thermodynamics as thermo
import neoswga.core.reaction_conditions as rc
import neoswga.core.secondary_structure as ss
import neoswga.core.utility


@dataclass
class SearchResult:
    """Results from adaptive primer length search."""
    optimal_k: int
    primers: List[str]
    metrics: Dict
    all_k_results: Dict[int, Dict]


class AdaptivePrimerSearch:
    """
    Adaptive search for optimal primer length and candidates.

    Uses multi-objective optimization to balance:
    - Specificity (fg/bg ratio)
    - Coverage (number of binding sites)
    - Thermodynamic favorability
    - Secondary structure quality
    """

    def __init__(self,
                 fg_genomes: List[str],
                 bg_genomes: List[str],
                 fg_prefixes: List[str],
                 bg_prefixes: List[str],
                 conditions: rc.ReactionConditions,
                 min_specificity: float = 100.0,
                 min_coverage_sites: int = 100,
                 max_primers_per_k: int = 1000):
        """
        Initialize adaptive search.

        Args:
            fg_genomes: Paths to foreground genome files
            bg_genomes: Paths to background genome files
            fg_prefixes: Prefixes for foreground kmer files
            bg_prefixes: Prefixes for background kmer files
            conditions: Reaction conditions
            min_specificity: Minimum fg/bg binding ratio
            min_coverage_sites: Minimum binding sites in foreground
            max_primers_per_k: Max primers to keep per k-mer length
        """
        self.fg_genomes = fg_genomes
        self.bg_genomes = bg_genomes
        self.fg_prefixes = fg_prefixes
        self.bg_prefixes = bg_prefixes
        self.conditions = conditions
        self.min_specificity = min_specificity
        self.min_coverage_sites = min_coverage_sites
        self.max_primers_per_k = max_primers_per_k

        # Calculate genome lengths
        self.fg_lengths = neoswga.core.utility.get_all_seq_lengths(fg_genomes)
        self.bg_lengths = neoswga.core.utility.get_all_seq_lengths(bg_genomes)
        self.fg_total_length = sum(self.fg_lengths)
        self.bg_total_length = sum(self.bg_lengths)

        # Get optimal k range from conditions
        self.min_k = conditions.min_primer_length()
        self.max_k = conditions.max_primer_length()

        # Get polymerase temperature range
        self.temp_range = conditions.get_polymerase_range()

        print(f"Adaptive search initialized:")
        print(f"  K-mer range: {self.min_k}-{self.max_k}")
        print(f"  Target Tm range: {self.temp_range[0]:.1f}-{self.temp_range[1]:.1f}°C")
        print(f"  Reaction temp: {conditions.temp:.1f}°C")
        print(f"  Foreground: {self.fg_total_length:,} bp")
        print(f"  Background: {self.bg_total_length:,} bp")

    def search(self, verbose: bool = True) -> SearchResult:
        """
        Perform adaptive search across k-mer lengths.

        Returns:
            SearchResult with optimal k and primer candidates
        """
        results_by_k = {}

        for k in range(self.min_k, self.max_k + 1):
            if verbose:
                print(f"\n{'='*60}")
                print(f"Searching k={k}...")
                print(f"{'='*60}")

            result = self._search_single_k(k, verbose)
            results_by_k[k] = result

            if verbose:
                print(f"\nK={k} summary:")
                print(f"  Candidates: {result['num_candidates']}")
                print(f"  Mean specificity: {result['mean_specificity']:.1f}")
                print(f"  Mean Tm: {result['mean_tm']:.1f}°C")
                print(f"  Score: {result['composite_score']:.3f}")

        # Select optimal k
        optimal_k = self._select_optimal_k(results_by_k)
        optimal_result = results_by_k[optimal_k]

        if verbose:
            print(f"\n{'='*60}")
            print(f"OPTIMAL K = {optimal_k}")
            print(f"{'='*60}")
            print(f"  Primers: {len(optimal_result['primers'])}")
            print(f"  Mean specificity: {optimal_result['mean_specificity']:.1f}")
            print(f"  Expected coverage: {optimal_result['expected_coverage']:.1%}")

        return SearchResult(
            optimal_k=optimal_k,
            primers=optimal_result['primers'],
            metrics=optimal_result,
            all_k_results=results_by_k
        )

    def _search_single_k(self, k: int, verbose: bool = False) -> Dict:
        """
        Search for optimal primers at specific k-mer length.

        Args:
            k: K-mer length
            verbose: Print progress

        Returns:
            Dictionary with primers and metrics
        """
        # Step 1: Get all k-mers from foreground
        if verbose:
            print(f"Step 1: Loading {k}-mers from foreground...")

        all_kmers = self._get_kmers_from_foreground(k)

        if verbose:
            print(f"  Found {len(all_kmers):,} unique {k}-mers")

        # Step 2: Thermodynamic filtering
        if verbose:
            print(f"Step 2: Thermodynamic filtering...")

        tm_filtered = self._filter_by_tm(all_kmers, verbose)

        if verbose:
            print(f"  Passed Tm filter: {len(tm_filtered):,} ({len(tm_filtered)/len(all_kmers):.1%})")

        # Step 3: Secondary structure filtering
        if verbose:
            print(f"Step 3: Secondary structure filtering...")

        struct_filtered = self._filter_by_structure(tm_filtered, verbose)

        if verbose:
            print(f"  Passed structure filter: {len(struct_filtered):,} ({len(struct_filtered)/len(tm_filtered):.1%})")

        # Step 4: Specificity filtering
        if verbose:
            print(f"Step 4: Specificity filtering...")

        spec_filtered = self._filter_by_specificity(struct_filtered, verbose)

        if verbose:
            print(f"  Passed specificity filter: {len(spec_filtered):,}")

        # Step 5: Rank and select top primers
        if verbose:
            print(f"Step 5: Ranking primers...")

        ranked_primers = self._rank_primers(spec_filtered, k)
        top_primers = ranked_primers[:self.max_primers_per_k]

        # Calculate metrics
        metrics = self._calculate_metrics(top_primers, k)
        metrics['primers'] = [p['seq'] for p in top_primers]
        metrics['num_candidates'] = len(top_primers)

        return metrics

    def _get_kmers_from_foreground(self, k: int) -> Set[str]:
        """Get all k-mers from foreground genomes."""
        kmers = set()

        for prefix in self.fg_prefixes:
            kmer_file = f"{prefix}_{k}mer_all.txt"
            try:
                with open(kmer_file, 'r') as f:
                    for line in f:
                        parts = line.strip().split()
                        if len(parts) >= 2:
                            kmer = parts[0]
                            if len(kmer) == k:
                                kmers.add(kmer.upper())
            except FileNotFoundError:
                warnings.warn(f"K-mer file not found: {kmer_file}")

        return kmers

    def _filter_by_tm(self, kmers: Set[str], verbose: bool = False) -> List[str]:
        """
        Filter primers by melting temperature.

        Keeps primers with Tm in polymerase optimal range.
        """
        passed = []
        tm_min, tm_max = self.temp_range

        for kmer in kmers:
            tm = self.conditions.calculate_effective_tm(kmer)

            # Allow 3°C margin on either side
            if (tm_min - 3) <= tm <= (tm_max + 3):
                passed.append(kmer)

        return passed

    def _filter_by_structure(self, kmers: List[str], verbose: bool = False) -> List[str]:
        """
        Filter primers by secondary structure.

        Removes primers with:
        - Stable hairpins
        - Strong self-dimers
        """
        passed = []

        for kmer in kmers:
            # Check hairpins
            hairpins = ss.check_hairpins(kmer, self.conditions)
            if hairpins:
                max_hairpin_tm = max(h['tm'] for h in hairpins)
                # Reject if hairpin stable above reaction temp + 5°C margin
                if max_hairpin_tm > self.conditions.temp + 5:
                    continue

            # Check self-dimers
            self_dimer = ss.check_homodimer(kmer, self.conditions)
            # Reject if self-dimer severity > 0.3
            if self_dimer['severity'] > 0.3:
                continue

            passed.append(kmer)

        return passed

    def _filter_by_specificity(self, kmers: List[str], verbose: bool = False) -> List[Dict]:
        """
        Filter primers by foreground/background specificity.

        Returns list of dicts with primer and binding counts.
        """
        passed = []

        for kmer in kmers:
            # Count occurrences in foreground and background
            fg_count = self._count_kmer_occurrences(kmer, self.fg_genomes)
            bg_count = self._count_kmer_occurrences(kmer, self.bg_genomes)

            # Normalize by genome length
            fg_freq = fg_count / self.fg_total_length
            bg_freq = bg_count / self.bg_total_length if bg_count > 0 else 0

            # Calculate specificity ratio
            if bg_freq > 0:
                specificity = fg_freq / bg_freq
            else:
                specificity = np.inf

            # Check minimum requirements
            if fg_count >= self.min_coverage_sites and specificity >= self.min_specificity:
                passed.append({
                    'seq': kmer,
                    'fg_count': fg_count,
                    'bg_count': bg_count,
                    'fg_freq': fg_freq,
                    'bg_freq': bg_freq,
                    'specificity': specificity
                })

        return passed

    def _count_kmer_occurrences(self, kmer: str, genome_files: List[str]) -> int:
        """Count occurrences of k-mer in genome files."""
        count = 0
        kmer_rc = thermo.reverse_complement(kmer)

        for genome_file in genome_files:
            try:
                with open(genome_file, 'r') as f:
                    seq = ""
                    for line in f:
                        if not line.startswith('>'):
                            seq += line.strip().upper()

                    # Count both forward and reverse complement
                    count += seq.count(kmer)
                    count += seq.count(kmer_rc)
            except FileNotFoundError:
                warnings.warn(f"Genome file not found: {genome_file}")

        return count

    def _rank_primers(self, primers: List[Dict], k: int) -> List[Dict]:
        """
        Rank primers by composite score.

        Considers:
        - Specificity (fg/bg ratio)
        - Coverage (number of binding sites)
        - Thermodynamic favorability
        - GC content balance
        """
        for primer_data in primers:
            seq = primer_data['seq']

            # Specificity score (normalized log scale)
            spec_score = np.log10(min(primer_data['specificity'], 10000)) / 4  # 0-1

            # Coverage score (sigmoid)
            coverage_score = 1 / (1 + np.exp(-(primer_data['fg_count'] - 200) / 50))

            # Tm optimality score (Gaussian around midpoint)
            tm = self.conditions.calculate_effective_tm(seq)
            tm_midpoint = np.mean(self.temp_range)
            tm_score = np.exp(-((tm - tm_midpoint) / 5)**2)

            # GC content balance (optimal 40-60%)
            gc = thermo.gc_content(seq)
            gc_score = 1 - abs(gc - 0.5) / 0.5  # 1 at 50%, 0 at 0/100%

            # Free energy score (more negative = better)
            dg = thermo.calculate_free_energy(seq, self.conditions.temp)
            dg_score = 1 / (1 + np.exp((dg + 10) / 3))  # Sigmoid around -10

            # Composite score (weighted average)
            composite = (
                0.40 * spec_score +
                0.25 * coverage_score +
                0.20 * tm_score +
                0.10 * gc_score +
                0.05 * dg_score
            )

            primer_data['tm'] = tm
            primer_data['gc'] = gc
            primer_data['dg'] = dg
            primer_data['scores'] = {
                'specificity': spec_score,
                'coverage': coverage_score,
                'tm': tm_score,
                'gc': gc_score,
                'dg': dg_score
            }
            primer_data['composite_score'] = composite

        # Sort by composite score
        primers.sort(key=lambda x: x['composite_score'], reverse=True)

        return primers

    def _calculate_metrics(self, primers: List[Dict], k: int) -> Dict:
        """Calculate aggregate metrics for primer set."""
        if not primers:
            return {
                'num_candidates': 0,
                'mean_specificity': 0,
                'mean_tm': 0,
                'mean_gc': 0,
                'expected_coverage': 0,
                'composite_score': 0
            }

        return {
            'k': k,
            'num_candidates': len(primers),
            'mean_specificity': np.mean([p['specificity'] for p in primers]),
            'median_specificity': np.median([p['specificity'] for p in primers]),
            'mean_tm': np.mean([p['tm'] for p in primers]),
            'std_tm': np.std([p['tm'] for p in primers]),
            'mean_gc': np.mean([p['gc'] for p in primers]),
            'std_gc': np.std([p['gc'] for p in primers]),
            'total_fg_sites': sum(p['fg_count'] for p in primers),
            'expected_coverage': self._estimate_coverage(primers),
            'composite_score': np.mean([p['composite_score'] for p in primers])
        }

    def _estimate_coverage(self, primers: List[Dict]) -> float:
        """
        Estimate genome coverage from primer binding sites.

        Assumes random binding and calculates expected coverage
        using Poisson statistics.
        """
        if not primers:
            return 0.0

        # Average binding density (sites per bp)
        total_sites = sum(p['fg_count'] for p in primers)
        density = total_sites / self.fg_total_length

        # Expected coverage using Poisson: 1 - exp(-lambda)
        # where lambda = density * amplicon_length
        assumed_amplicon_length = 3000  # Typical SWGA amplicon
        lambda_param = density * assumed_amplicon_length

        coverage = 1 - np.exp(-lambda_param)

        return coverage

    def _select_optimal_k(self, results_by_k: Dict[int, Dict]) -> int:
        """
        Select optimal k based on multi-objective criteria.

        Prioritizes:
        1. Sufficient candidates (>100)
        2. High specificity
        3. Good coverage
        4. Longer k for better specificity (tie-breaker)
        """
        valid_ks = []

        for k, result in results_by_k.items():
            # Must have sufficient candidates
            if result['num_candidates'] >= 100:
                valid_ks.append((k, result['composite_score']))

        if not valid_ks:
            # Fallback: pick k with most candidates
            k_by_count = sorted(results_by_k.items(),
                              key=lambda x: x[1]['num_candidates'],
                              reverse=True)
            return k_by_count[0][0]

        # Select k with highest composite score
        # If tie, prefer longer k
        valid_ks.sort(key=lambda x: (x[1], x[0]), reverse=True)

        return valid_ks[0][0]


def adaptive_search_pipeline(fg_genomes: List[str],
                             bg_genomes: List[str],
                             fg_prefixes: List[str],
                             bg_prefixes: List[str],
                             conditions: Optional[rc.ReactionConditions] = None,
                             **kwargs) -> SearchResult:
    """
    High-level function for adaptive primer search.

    Args:
        fg_genomes: Foreground genome files
        bg_genomes: Background genome files
        fg_prefixes: Foreground kmer prefixes
        bg_prefixes: Background kmer prefixes
        conditions: Reaction conditions (default: enhanced)
        **kwargs: Additional search parameters

    Returns:
        SearchResult with optimal k and primers
    """
    if conditions is None:
        conditions = rc.get_enhanced_conditions()

    searcher = AdaptivePrimerSearch(
        fg_genomes, bg_genomes,
        fg_prefixes, bg_prefixes,
        conditions,
        **kwargs
    )

    return searcher.search(verbose=True)


if __name__ == "__main__":
    # Example usage
    print("Adaptive Primer Search Module")
    print("=" * 60)
    print("\nThis module provides intelligent k-mer length selection")
    print("based on reaction conditions and thermodynamic optimization.")
    print("\nKey features:")
    print("  - Adaptive k-mer length (6-15 based on conditions)")
    print("  - Thermodynamic filtering (Tm, ΔG)")
    print("  - Secondary structure filtering (hairpins, dimers)")
    print("  - Multi-objective ranking (specificity, coverage, Tm)")
    print("\nExample usage:")
    print("""
    from neoswga.core import adaptive_search, reaction_conditions

    # Set up enhanced conditions
    conditions = reaction_conditions.get_enhanced_conditions()

    # Run adaptive search
    result = adaptive_search.adaptive_search_pipeline(
        fg_genomes=['target.fasta'],
        bg_genomes=['offtarget.fasta'],
        fg_prefixes=['target_kmers'],
        bg_prefixes=['offtarget_kmers'],
        conditions=conditions
    )

    print(f"Optimal k: {result.optimal_k}")
    print(f"Primers found: {len(result.primers)}")
    """)
