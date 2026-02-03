#!/usr/bin/env python3
"""
Diagnostic script to analyze terminal Tm distribution for different primer lengths
from Francisella genome.
"""

import sys
from pathlib import Path
from Bio import SeqIO
import numpy as np

sys.path.insert(0, str(Path.cwd()))

from neoswga.core.three_prime_stability import ThreePrimeStabilityAnalyzer
from neoswga.core.reaction_conditions import ReactionConditions

def load_genome(path, max_bp=500000):
    """Load genome sequence."""
    seq = ""
    with open(path) as f:
        for record in SeqIO.parse(f, "fasta"):
            seq += str(record.seq).upper()
            if len(seq) >= max_bp:
                break
    return seq[:max_bp]

def generate_kmers(seq, k, max_count=2000):
    """Generate k-mers from sequence."""
    kmers = []
    for i in range(min(len(seq) - k + 1, max_count * 10)):
        kmer = seq[i:i+k]
        if all(base in 'ATGC' for base in kmer):
            # Filter by GC (30-70%)
            gc = (kmer.count('G') + kmer.count('C')) / len(kmer) * 100
            if 30 <= gc <= 70:
                kmers.append(kmer)
        if len(kmers) >= max_count:
            break
    return kmers

def analyze_terminal_tm_distribution(length_range, genome_seq):
    """Analyze terminal Tm distribution for primers of given length."""
    from neoswga.core.three_prime_stability import create_three_prime_analyzer

    conditions = ReactionConditions(temp=30.0, mg_conc=3.0, na_conc=50.0)

    # Create analyzer using PROPER factory with SWGA mode (no salt correction for <13bp)
    # Use lenient mode with min_terminal_tm=0 to analyze all primers
    analyzer = create_three_prime_analyzer(
        stringency='lenient',
        conditions=conditions,
        swga_mode=True  # IMPORTANT: Use SWGA mode (no salt correction for short primers!)
    )
    # Override threshold to 0 to capture all primers
    analyzer.min_terminal_tm = 0

    terminal_tms = []

    for k in range(length_range[0], length_range[1] + 1):
        print(f"\nAnalyzing {k}-mers...")
        kmers = generate_kmers(genome_seq, k, max_count=500)
        print(f"  Generated {len(kmers)} candidates")

        for primer in kmers:
            result = analyzer.analyze_primer(primer)
            terminal_tms.append(result.terminal_tm)

    return np.array(terminal_tms)

def main():
    genome_path = "/Users/andreassjodin/Code/swga-dev/test/francisella_GCF_000008985.1_ASM898v1_genomic.fna"

    print("="*80)
    print("Terminal Tm Distribution Analysis - Francisella")
    print("="*80)

    # Load genome
    print("\nLoading genome...")
    genome_seq = load_genome(genome_path, max_bp=500000)
    gc_content = (genome_seq.count('G') + genome_seq.count('C')) / len(genome_seq) * 100
    print(f"Genome: {len(genome_seq):,} bp, GC = {gc_content:.1f}%")

    # Analyze different length ranges
    length_ranges = [
        (11, 13, "11-13bp (Config 4 - WORKS)"),
        (12, 15, "12-15bp (Configs 1,2,3 - WORK)"),
        (13, 16, "13-16bp (Config 5 - FAILS)"),
        (15, 18, "15-18bp (Configs 7,8 - FAIL)")
    ]

    for min_len, max_len, description in length_ranges:
        print(f"\n{'='*80}")
        print(f"{description}")
        print(f"{'='*80}")

        terminal_tms = analyze_terminal_tm_distribution((min_len, max_len), genome_seq)

        print(f"\nTerminal Tm Statistics:")
        print(f"  Count: {len(terminal_tms)}")
        print(f"  Mean:  {np.mean(terminal_tms):.2f}°C")
        print(f"  Std:   {np.std(terminal_tms):.2f}°C")
        print(f"  Min:   {np.min(terminal_tms):.2f}°C")
        print(f"  Max:   {np.max(terminal_tms):.2f}°C")
        print(f"  Median: {np.median(terminal_tms):.2f}°C")

        # Check pass rates at different thresholds
        print(f"\nPass rates at SWGA thresholds:")
        lenient_pass = np.sum(terminal_tms >= 10.0) / len(terminal_tms) * 100
        moderate_pass = np.sum(terminal_tms >= 14.0) / len(terminal_tms) * 100
        strict_pass = np.sum(terminal_tms >= 18.0) / len(terminal_tms) * 100

        print(f"  Lenient  (≥10°C): {lenient_pass:5.1f}%")
        print(f"  Moderate (≥14°C): {moderate_pass:5.1f}%")
        print(f"  Strict   (≥18°C): {strict_pass:5.1f}%")

        # Distribution by 2°C bins
        print(f"\nTerminal Tm Distribution:")
        bins = [(i, i+2) for i in range(6, 22, 2)]
        for low, high in bins:
            count = np.sum((terminal_tms >= low) & (terminal_tms < high))
            pct = count / len(terminal_tms) * 100
            print(f"  {low:2d}-{high:2d}°C: {count:4d} primers ({pct:5.1f}%)")

    print(f"\n{'='*80}")
    print("ANALYSIS COMPLETE")
    print(f"{'='*80}")

if __name__ == "__main__":
    main()
