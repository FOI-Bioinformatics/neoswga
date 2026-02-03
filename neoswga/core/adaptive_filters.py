"""
Adaptive filtering that adjusts to genome composition.

Fixes critical bug: Current fixed GC thresholds (37.5-62.5%) fail for
GC-extreme organisms like Francisella (33% GC) and Burkholderia (67% GC).
"""

import numpy as np
from typing import List, Dict, Tuple, Optional
from dataclasses import dataclass
import logging
from collections import Counter
import os

logger = logging.getLogger(__name__)


@dataclass
class FilterThresholds:
    """Adaptive thresholds based on genome properties"""
    gc_min: float
    gc_max: float
    tm_min: float
    tm_max: float
    max_repeat_len: int
    max_self_complementarity: int

    def __str__(self):
        return (f"GC: [{self.gc_min:.3f}, {self.gc_max:.3f}], "
                f"Tm: [{self.tm_min:.1f}, {self.tm_max:.1f}]°C")


class AdaptiveGCFilter:
    """
    GC content filter that adapts to genome composition.

    Current bug: filter.py uses fixed thresholds (37.5-62.5%)
    This rejects ALL primers for extreme GC organisms.

    Fix: Adapt thresholds to genome GC content.
    """

    def __init__(self, genome_gc: float, tolerance: float = 0.15):
        """
        Initialize adaptive filter.

        Args:
            genome_gc: GC content of target genome (0.0-1.0)
            tolerance: Allowed deviation from genome GC (default ±15%)

        Examples:
            Francisella (33% GC): Accept primers 18-48% GC
            E. coli (50% GC): Accept primers 35-65% GC
            Burkholderia (67% GC): Accept primers 52-82% GC
        """
        self.genome_gc = genome_gc
        self.tolerance = tolerance

        # Compute adaptive bounds
        self.gc_min = max(0.20, genome_gc - tolerance)
        self.gc_max = min(0.80, genome_gc + tolerance)

        logger.info(f"Adaptive GC filter: genome={genome_gc:.3f}, "
                   f"range=[{self.gc_min:.3f}, {self.gc_max:.3f}]")

    def passes(self, primer: str) -> bool:
        """Check if primer passes GC filter"""
        gc = self._gc_content(primer)
        return self.gc_min <= gc <= self.gc_max

    def _gc_content(self, seq: str) -> float:
        """Compute GC fraction"""
        gc_count = seq.count('G') + seq.count('C')
        return gc_count / len(seq) if len(seq) > 0 else 0.0

    def explain_rejection(self, primer: str) -> str:
        """Explain why primer was rejected"""
        gc = self._gc_content(primer)
        if gc < self.gc_min:
            return f"GC too low: {gc:.3f} < {self.gc_min:.3f}"
        elif gc > self.gc_max:
            return f"GC too high: {gc:.3f} > {self.gc_max:.3f}"
        else:
            return "PASSES"


class ThermodynamicFilter:
    """
    Filter based on thermodynamic properties appropriate for SWGA.

    Key insight: SWGA uses isothermal amplification (30°C), NOT PCR.
    Different rules apply!
    """

    def __init__(self, reaction_temp: float = 30.0, na_conc: float = 50.0,
                 target_tm_offset: float = 7.5):
        """
        Initialize thermodynamic filter.

        Args:
            reaction_temp: Reaction temperature (°C)
            na_conc: Sodium concentration (mM)
            target_tm_offset: Target Tm = reaction_temp + offset
        """
        self.reaction_temp = reaction_temp
        self.na_conc = na_conc
        self.target_tm = reaction_temp + target_tm_offset

        # Tm range: ±5°C around target
        self.tm_min = self.target_tm - 5
        self.tm_max = self.target_tm + 5

        logger.info(f"Thermodynamic filter: target Tm={self.target_tm:.1f}°C, "
                   f"range=[{self.tm_min:.1f}, {self.tm_max:.1f}]")

    def passes(self, primer: str) -> bool:
        """Check if primer Tm is in acceptable range"""
        tm = self._calculate_tm(primer)
        return self.tm_min <= tm <= self.tm_max

    def _calculate_tm(self, seq: str) -> float:
        """
        Calculate melting temperature with salt correction.

        Uses nearest-neighbor thermodynamics + SantaLucia parameters.
        """
        # For short primers, use simple approximation
        if len(seq) <= 13:
            # Wallace rule with GC adjustment
            gc_count = seq.count('G') + seq.count('C')
            at_count = len(seq) - gc_count
            tm_basic = 2 * at_count + 4 * gc_count

            # Salt correction (von Ahsen et al.)
            salt_factor = 16.6 * np.log10(self.na_conc / 1000.0)  # mM to M

            return tm_basic + salt_factor
        else:
            # Use nearest-neighbor (more accurate for longer primers)
            return self._nn_tm(seq)

    def _nn_tm(self, seq: str) -> float:
        """
        Nearest-neighbor Tm calculation.

        Simplified version - for production, use Bio.SeqUtils.MeltingTemp
        """
        # SantaLucia unified NN parameters (kcal/mol)
        nn_table = {
            'AA': (-7.9, -22.2), 'TT': (-7.9, -22.2),
            'AT': (-7.2, -20.4), 'TA': (-7.2, -21.3),
            'CA': (-8.5, -22.7), 'TG': (-8.5, -22.7),
            'GT': (-8.4, -22.4), 'AC': (-8.4, -22.4),
            'CT': (-7.8, -21.0), 'AG': (-7.8, -21.0),
            'GA': (-8.2, -22.2), 'TC': (-8.2, -22.2),
            'CG': (-10.6, -27.2), 'GC': (-9.8, -24.4),
            'GG': (-8.0, -19.9), 'CC': (-8.0, -19.9),
        }

        dH = 0.0  # Enthalpy (kcal/mol)
        dS = 0.0  # Entropy (cal/mol·K)

        # Sum nearest-neighbor contributions
        for i in range(len(seq) - 1):
            dinuc = seq[i:i+2]
            if dinuc in nn_table:
                h, s = nn_table[dinuc]
                dH += h
                dS += s

        # Initiation parameters
        dH += 0.2  # Helix initiation
        dS += -5.7

        # Salt correction for entropy (SantaLucia 1998)
        dS += 0.368 * (len(seq) - 1) * np.log(self.na_conc / 1000.0)

        # Tm = dH / (dS + R * ln(C/4))
        # Simplified: Tm ≈ dH / dS - 273.15
        if dS != 0:
            tm = (dH * 1000) / dS - 273.15  # Convert to Celsius
            return tm
        else:
            return 0.0

    def explain_rejection(self, primer: str) -> str:
        """Explain why primer was rejected"""
        tm = self._calculate_tm(primer)
        if tm < self.tm_min:
            return f"Tm too low: {tm:.1f}°C < {self.tm_min:.1f}°C"
        elif tm > self.tm_max:
            return f"Tm too high: {tm:.1f}°C > {self.tm_max:.1f}°C"
        else:
            return "PASSES"


class RepeatFilter:
    """
    Filter primers with excessive repeats (homopolymers, dinucleotide repeats).

    More permissive than PCR rules (SWGA is more tolerant).
    """

    def __init__(self, max_homopolymer: int = 5, max_dinuc_repeat: int = 5):
        """
        Initialize repeat filter.

        Args:
            max_homopolymer: Max consecutive identical bases (default 5)
            max_dinuc_repeat: Max dinucleotide repeats (default 5, i.e., 10 bp)
        """
        self.max_homopolymer = max_homopolymer
        self.max_dinuc_repeat = max_dinuc_repeat

    def passes(self, primer: str) -> bool:
        """Check if primer passes repeat filters"""
        # Check homopolymers
        for base in 'ATCG':
            if base * (self.max_homopolymer + 1) in primer:
                return False

        # Check dinucleotide repeats
        if len(primer) >= 10:
            for i in range(len(primer) - 9):
                dinuc = primer[i:i+2]
                if dinuc * self.max_dinuc_repeat in primer:
                    return False

        return True

    def explain_rejection(self, primer: str) -> str:
        """Explain why primer was rejected"""
        # Check homopolymers
        for base in 'ATCG':
            repeat = base * (self.max_homopolymer + 1)
            if repeat in primer:
                return f"Homopolymer: {repeat} found"

        # Check dinucleotide repeats
        for i in range(len(primer) - 9):
            dinuc = primer[i:i+2]
            repeat = dinuc * self.max_dinuc_repeat
            if repeat in primer:
                return f"Dinucleotide repeat: {repeat} found"

        return "PASSES"


class GCClampFilter:
    """
    Filter based on GC clamp at 3' end.

    GC clamp promotes specific binding, but excessive GC can cause problems.
    """

    def __init__(self, min_gc_in_last5: int = 1, max_gc_in_last5: int = 3):
        """
        Initialize GC clamp filter.

        Args:
            min_gc_in_last5: Minimum G/C in last 5 bases (default 1)
            max_gc_in_last5: Maximum G/C in last 5 bases (default 3)
        """
        self.min_gc = min_gc_in_last5
        self.max_gc = max_gc_in_last5

    def passes(self, primer: str) -> bool:
        """Check GC clamp"""
        if len(primer) < 5:
            return True

        last_five = primer[-5:]
        gc_count = last_five.count('G') + last_five.count('C')

        return self.min_gc <= gc_count <= self.max_gc

    def explain_rejection(self, primer: str) -> str:
        """Explain rejection"""
        last_five = primer[-5:]
        gc_count = last_five.count('G') + last_five.count('C')

        if gc_count < self.min_gc:
            return f"GC clamp too weak: {gc_count} < {self.min_gc}"
        elif gc_count > self.max_gc:
            return f"GC clamp too strong: {gc_count} > {self.max_gc}"
        else:
            return "PASSES"


class AdaptiveFilterPipeline:
    """
    Complete adaptive filtering pipeline.

    Automatically adjusts all filters to genome properties.
    """

    def __init__(self, target_genome: str,
                 reaction_temp: float = 30.0,
                 na_conc: float = 50.0,
                 gc_tolerance: float = 0.15):
        """
        Initialize adaptive pipeline.

        Args:
            target_genome: Target genome sequence (or path to FASTA)
            reaction_temp: Reaction temperature
            na_conc: Sodium concentration (mM)
            gc_tolerance: GC tolerance (default: ±15%)
        """
        # Analyze genome
        if os.path.isfile(target_genome):
            genome_seq = self._load_genome(target_genome)
        else:
            genome_seq = target_genome

        genome_gc = self._compute_gc(genome_seq)
        logger.info(f"Target genome GC: {genome_gc:.3f}")

        # Initialize adaptive filters
        self.gc_filter = AdaptiveGCFilter(genome_gc, tolerance=gc_tolerance)
        self.tm_filter = ThermodynamicFilter(reaction_temp, na_conc)
        self.repeat_filter = RepeatFilter(max_homopolymer=5, max_dinuc_repeat=5)
        self.clamp_filter = GCClampFilter(min_gc_in_last5=1, max_gc_in_last5=3)

        self.filters = {
            'gc': self.gc_filter,
            'tm': self.tm_filter,
            'repeat': self.repeat_filter,
            'clamp': self.clamp_filter,
        }

    def filter_primers(self, candidates: List[str], verbose: bool = False) -> List[str]:
        """
        Apply all filters.

        Args:
            candidates: List of candidate primers
            verbose: Print rejection reasons

        Returns:
            Filtered primer list
        """
        passed = []
        rejected_by = Counter()

        for primer in candidates:
            failed = []

            for name, filter_obj in self.filters.items():
                if not filter_obj.passes(primer):
                    failed.append(name)
                    rejected_by[name] += 1

            if not failed:
                passed.append(primer)
            elif verbose:
                reasons = [f"{name}: {self.filters[name].explain_rejection(primer)}"
                          for name in failed]
                logger.debug(f"{primer} rejected: {'; '.join(reasons)}")

        logger.info(f"Adaptive filtering: {len(candidates)} → {len(passed)} "
                   f"({100*len(passed)/len(candidates):.1f}% passed)")
        logger.info(f"Rejection reasons: {dict(rejected_by)}")

        return passed

    def get_thresholds(self) -> FilterThresholds:
        """Get current filter thresholds"""
        return FilterThresholds(
            gc_min=self.gc_filter.gc_min,
            gc_max=self.gc_filter.gc_max,
            tm_min=self.tm_filter.tm_min,
            tm_max=self.tm_filter.tm_max,
            max_repeat_len=self.repeat_filter.max_homopolymer,
            max_self_complementarity=3  # From clamp filter
        )

    def _load_genome(self, fasta_path: str) -> str:
        """Load genome from FASTA"""
        from Bio import SeqIO
        sequences = []
        for record in SeqIO.parse(fasta_path, "fasta"):
            sequences.append(str(record.seq).upper())
        return ''.join(sequences)

    def _compute_gc(self, seq: str) -> float:
        """Compute GC content"""
        gc = seq.count('G') + seq.count('C')
        return gc / len(seq) if len(seq) > 0 else 0.5


def compare_filters(primers: List[str], genome_seq: str):
    """
    Compare old (fixed) vs. new (adaptive) filters.

    Demonstrates how many primers are lost with fixed thresholds.
    """
    genome_gc = sum(1 for b in genome_seq if b in 'GC') / len(genome_seq)

    # Old filter (fixed thresholds)
    old_passed = []
    for primer in primers:
        primer_gc = sum(1 for b in primer if b in 'GC') / len(primer)
        if 0.375 <= primer_gc <= 0.625:  # Fixed thresholds!
            old_passed.append(primer)

    # New filter (adaptive)
    new_pipeline = AdaptiveFilterPipeline(genome_seq)
    new_passed = new_pipeline.filter_primers(primers)

    print(f"Genome GC: {genome_gc:.3f}")
    print(f"Old filter (fixed): {len(old_passed)}/{len(primers)} passed")
    print(f"New filter (adaptive): {len(new_passed)}/{len(primers)} passed")
    print(f"Improvement: {len(new_passed) - len(old_passed)} additional primers")

    # Show examples of newly accepted primers
    newly_accepted = set(new_passed) - set(old_passed)
    if newly_accepted:
        print(f"\nExample newly accepted primers:")
        for primer in list(newly_accepted)[:5]:
            primer_gc = sum(1 for b in primer if b in 'GC') / len(primer)
            print(f"  {primer} (GC={primer_gc:.3f})")


if __name__ == "__main__":
    import os
    import sys

    logging.basicConfig(level=logging.INFO)

    if len(sys.argv) > 1:
        # Test on provided genome
        genome_file = sys.argv[1]

        print("Testing adaptive filters...")
        pipeline = AdaptiveFilterPipeline(genome_file)

        # Generate test primers
        test_primers = [
            'ATCGATCG',  # 50% GC
            'AAAAAAAA',  # 0% GC (should fail repeat filter)
            'GCGCGCGC',  # 100% GC (may fail GC filter for AT-rich genomes)
            'ATATATAT',  # Dinuc repeat (should fail)
            'ATCGGGGG',  # Homopolymer (should fail)
            'ATCGATAA',  # No GC clamp (should fail)
        ]

        passed = pipeline.filter_primers(test_primers, verbose=True)
        print(f"\nPassed: {passed}")
        print(f"\nThresholds: {pipeline.get_thresholds()}")
    else:
        print("Usage: python adaptive_filters.py <genome.fasta>")


def run_step2_with_adaptive_gc(gc_tolerance: float = 0.15):
    """
    Run filter with adaptive GC filtering instead of fixed thresholds.

    This is a drop-in replacement that wraps the original filter
    but replaces the GC filter with an adaptive one.

    Args:
        gc_tolerance: GC tolerance (default: ±15%)
    """
    import neoswga.core.parameter as parameter
    import neoswga.core.filter
    import pandas as pd
    import os

    logger.info("Running filter with adaptive GC filtering")

    # Read parameters from JSON file to get genome paths
    # We need to read directly since get_params() hasn't been called yet
    import json
    import sys

    # Find JSON file from sys.argv or use params.json
    json_file = None
    for i, arg in enumerate(sys.argv):
        if arg in ['-j', '--json']:
            if i + 1 < len(sys.argv):
                json_file = sys.argv[i + 1]
                break

    if not json_file:
        json_file = 'params.json'

    if not os.path.exists(json_file):
        raise ValueError(f"JSON parameter file not found: {json_file}")

    with open(json_file, 'r') as f:
        params_data = json.load(f)

    # Get genome path
    fg_genomes = params_data.get('fg_genomes', [])
    if not fg_genomes or len(fg_genomes) == 0:
        raise ValueError("No target genome specified in params.json. Cannot compute adaptive GC thresholds.")

    fg_genome = fg_genomes[0]

    # Create adaptive filter
    logger.info(f"Computing adaptive GC thresholds from {fg_genome}")
    pipeline = AdaptiveFilterPipeline(fg_genome, gc_tolerance=gc_tolerance)
    thresholds = pipeline.get_thresholds()

    # Access dataclass attributes, not dictionary keys
    genome_gc = pipeline.gc_filter.genome_gc
    logger.info(f"Genome GC: {genome_gc:.3f}")
    logger.info(f"GC range: {thresholds.gc_min:.3f}-{thresholds.gc_max:.3f}")
    logger.info(f"Old fixed range: 0.375-0.625")

    # Monkey-patch filter_extra to use adaptive GC thresholds
    original_filter_extra = neoswga.core.filter.filter_extra
    gc_min = thresholds.gc_min
    gc_max = thresholds.gc_max

    def adaptive_filter_extra(primer):
        """Wraps filter_extra with adaptive GC thresholds"""
        import melting
        from collections import Counter
        import neoswga.core.parameter as parameter
        import neoswga.core.dimer as dimer

        primer_tm = melting.temp(primer)
        if primer_tm > 45 or primer_tm < 15:
            return False

        # Rule 5: No runs of 5+ nucleotides
        if 'GGGGG' in primer or 'CCCCC' in primer or 'AAAAA' in primer or 'TTTTT' in primer:
            return False

        # Rule 2: GC content (ADAPTIVE - replaces fixed 0.375-0.625)
        all_cnt = Counter(primer)
        gc_count = all_cnt.get('G', 0) + all_cnt.get('C', 0)
        GC_content = gc_count / float(len(primer))
        if GC_content <= gc_min or GC_content >= gc_max:
            return False

        # Rule 3: GC clamp (1-3 G/C in last 5 bases)
        last_five_count = Counter(primer[-5:])
        gc_last5 = last_five_count.get('G', 0) + last_five_count.get('C', 0)
        if gc_last5 > 3 or gc_last5 == 0:
            return False

        # Rule 1: No 3 G/C at 3' end
        last_three_count = Counter(primer[-3:])
        if last_three_count.get('G', 0) + last_three_count.get('C', 0) == 3:
            return False

        # Rule 4: Check dinucleotide repeats
        if len(primer) >= 10:
            more_than_5 = [nucleo for nucleo, count in all_cnt.items() if count >= 5]
            if len(more_than_5) > 1:
                # Check all dinucleotide repeat patterns
                repeats = ['ATATATATAT', 'TATATATATA', 'AGAGAGAGAG', 'GAGAGAGAGA',
                          'ACACACACAC', 'CACACACACA', 'TCTCTCTCTC', 'CTCTCTCTCT',
                          'GTGTGTGTGT', 'TGTGTGTGTG', 'CGCGCGCGCG', 'GCGCGCGCGC']
                if any(repeat in primer for repeat in repeats):
                    return False

        # Self-dimer check
        if dimer.is_dimer(primer, primer, parameter.max_self_dimer_bp):
            return False

        return True

    neoswga.core.filter.filter_extra = adaptive_filter_extra

    try:
        # Now run filter with patched filter
        logger.info("Running filter with adaptive filter...")

        # Temporarily remove --gc-tolerance from sys.argv since old pipeline doesn't know about it
        import sys
        original_argv = sys.argv.copy()
        if '--gc-tolerance' in sys.argv:
            idx = sys.argv.index('--gc-tolerance')
            sys.argv.pop(idx)  # Remove --gc-tolerance
            if idx < len(sys.argv):  # Remove the value too
                sys.argv.pop(idx)

        try:
            from neoswga.core import pipeline as core_pipeline
            core_pipeline.step2()
        finally:
            # Restore original sys.argv
            sys.argv = original_argv

        logger.info("Step2 complete with adaptive GC filtering!")

    finally:
        # Restore original filter
        neoswga.core.filter.filter_extra = original_filter_extra
