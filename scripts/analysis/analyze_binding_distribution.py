#!/usr/bin/env python3
"""
Analyze spatial distribution of primer binding sites in Francisella genome.

This script examines:
1. Binding site clustering vs uniform distribution
2. GC content correlation with binding density
3. Coverage gaps and their genomic features
4. Comparison across primer configurations
"""

import json
from pathlib import Path
from typing import List, Dict, Tuple
import logging
from collections import defaultdict
from Bio import SeqIO

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)


def load_primers(config_dir: Path) -> List[str]:
    """Load primers from configuration directory"""
    primers_file = config_dir / "primers.txt"
    primers = []
    with open(primers_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            # Format: PRIMER\tSCORE
            parts = line.split('\t')
            if parts:
                primers.append(parts[0])
    return primers


def find_binding_sites(genome_file: Path, primers: List[str]) -> Dict[str, List[int]]:
    """Find all binding positions for each primer in genome"""
    binding_sites = defaultdict(list)

    # Load genome
    genome_seq = ""
    for record in SeqIO.parse(genome_file, "fasta"):
        genome_seq += str(record.seq).upper()

    # Find binding sites for each primer
    for primer in primers:
        primer_upper = primer.upper()
        # Forward strand
        pos = 0
        while True:
            pos = genome_seq.find(primer_upper, pos)
            if pos == -1:
                break
            binding_sites[primer].append(pos)
            pos += 1

        # Reverse complement
        rev_comp = str(primer_upper.translate(str.maketrans('ATCG', 'TAGC'))[::-1])
        pos = 0
        while True:
            pos = genome_seq.find(rev_comp, pos)
            if pos == -1:
                break
            binding_sites[primer].append(pos)
            pos += 1

    return binding_sites, genome_seq


def calculate_gc_content_windows(genome_seq: str, window_size: int = 10000) -> List[Tuple[int, float]]:
    """Calculate GC content in sliding windows"""
    gc_windows = []

    for i in range(0, len(genome_seq) - window_size, window_size):
        window = genome_seq[i:i+window_size]
        gc_count = window.count('G') + window.count('C')
        gc_percent = (gc_count / len(window)) * 100
        gc_windows.append((i, gc_percent))

    return gc_windows


def analyze_coverage_gaps(binding_positions: List[int], genome_length: int, bin_size: int = 10000) -> Dict:
    """Analyze regions with no binding sites"""
    # Create bins
    n_bins = (genome_length // bin_size) + 1
    bins_covered = set()

    for pos in binding_positions:
        bin_idx = pos // bin_size
        bins_covered.add(bin_idx)

    # Find gaps
    gaps = []
    gap_start = None
    for i in range(n_bins):
        if i not in bins_covered:
            if gap_start is None:
                gap_start = i
        else:
            if gap_start is not None:
                gap_length = (i - gap_start) * bin_size
                gaps.append({
                    'start_bin': gap_start,
                    'end_bin': i - 1,
                    'start_pos': gap_start * bin_size,
                    'end_pos': i * bin_size,
                    'length': gap_length
                })
                gap_start = None

    # Close final gap if exists
    if gap_start is not None:
        gaps.append({
            'start_bin': gap_start,
            'end_bin': n_bins - 1,
            'start_pos': gap_start * bin_size,
            'end_pos': genome_length,
            'length': (n_bins - gap_start) * bin_size
        })

    return {
        'total_bins': n_bins,
        'covered_bins': len(bins_covered),
        'coverage_percent': (len(bins_covered) / n_bins) * 100,
        'gaps': gaps,
        'n_gaps': len(gaps),
        'total_gap_length': sum(g['length'] for g in gaps),
        'largest_gap': max(gaps, key=lambda g: g['length']) if gaps else None
    }


def calculate_clustering_metric(positions: List[int], genome_length: int) -> float:
    """
    Calculate clustering coefficient.

    Returns ratio of observed nearest-neighbor distance to expected (uniform).
    Values < 1 indicate clustering, > 1 indicate overdispersion.
    """
    if len(positions) < 2:
        return 1.0

    sorted_pos = sorted(positions)

    # Calculate nearest neighbor distances
    nn_distances = []
    for i in range(len(sorted_pos) - 1):
        nn_distances.append(sorted_pos[i+1] - sorted_pos[i])

    # Add circular distance
    nn_distances.append((genome_length - sorted_pos[-1]) + sorted_pos[0])

    mean_nn_distance = sum(nn_distances) / len(nn_distances)
    expected_distance = genome_length / len(positions)

    return mean_nn_distance / expected_distance


def analyze_config(config_name: str, config_dir: Path, genome_file: Path, genome_seq: str) -> Dict:
    """Analyze binding distribution for a single configuration"""

    logger.info(f"\n{'='*80}")
    logger.info(f"ANALYZING: {config_name}")
    logger.info(f"{'='*80}")

    # Load primers
    primers = load_primers(config_dir)
    logger.info(f"Primers: {len(primers)}")

    # Find binding sites
    binding_sites, _ = find_binding_sites(genome_file, primers)

    # Collect all binding positions
    all_positions = []
    for positions in binding_sites.values():
        all_positions.extend(positions)

    total_sites = len(all_positions)
    genome_length = len(genome_seq)

    logger.info(f"Total binding sites: {total_sites}")
    logger.info(f"Sites per primer: {total_sites / len(primers):.2f}")
    logger.info(f"Binding frequency: {(total_sites / genome_length) * 1000:.2f} per kb")

    # Analyze coverage gaps
    gap_analysis = analyze_coverage_gaps(all_positions, genome_length, bin_size=10000)
    logger.info(f"\nCoverage Analysis:")
    logger.info(f"  Bins covered: {gap_analysis['covered_bins']}/{gap_analysis['total_bins']} ({gap_analysis['coverage_percent']:.1f}%)")
    logger.info(f"  Number of gaps: {gap_analysis['n_gaps']}")
    logger.info(f"  Total gap length: {gap_analysis['total_gap_length']:,} bp ({gap_analysis['total_gap_length']/genome_length*100:.1f}%)")

    if gap_analysis['largest_gap']:
        lg = gap_analysis['largest_gap']
        logger.info(f"  Largest gap: {lg['length']:,} bp (bins {lg['start_bin']}-{lg['end_bin']})")

    # Clustering analysis
    clustering = calculate_clustering_metric(all_positions, genome_length)
    logger.info(f"\nSpatial Distribution:")
    logger.info(f"  Clustering coefficient: {clustering:.3f}")
    if clustering < 0.8:
        logger.info(f"  → CLUSTERED (sites are not evenly distributed)")
    elif clustering > 1.2:
        logger.info(f"  → OVERDISPERSED (sites are too evenly spaced)")
    else:
        logger.info(f"  → RANDOM (sites are uniformly distributed)")

    # Primer length distribution
    primer_lengths = [len(p) for p in primers]
    logger.info(f"\nPrimer Characteristics:")
    logger.info(f"  Length range: {min(primer_lengths)}-{max(primer_lengths)} bp")
    logger.info(f"  Mean length: {sum(primer_lengths)/len(primer_lengths):.1f} bp")

    # GC content of primers
    primer_gc = []
    for p in primers:
        gc = (p.count('G') + p.count('C')) / len(p) * 100
        primer_gc.append(gc)
    logger.info(f"  Mean primer GC: {sum(primer_gc)/len(primer_gc):.1f}%")

    # Return results
    return {
        'config_name': config_name,
        'n_primers': len(primers),
        'total_binding_sites': total_sites,
        'sites_per_primer': total_sites / len(primers),
        'binding_frequency_per_kb': (total_sites / genome_length) * 1000,
        'coverage_percent': gap_analysis['coverage_percent'],
        'n_gaps': gap_analysis['n_gaps'],
        'total_gap_length': gap_analysis['total_gap_length'],
        'largest_gap_length': gap_analysis['largest_gap']['length'] if gap_analysis['largest_gap'] else 0,
        'clustering_coefficient': clustering,
        'mean_primer_length': sum(primer_lengths) / len(primer_lengths),
        'mean_primer_gc': sum(primer_gc) / len(primer_gc),
        'gaps': gap_analysis['gaps'][:10]  # Top 10 largest gaps
    }


def analyze_gc_correlation(genome_file: Path) -> None:
    """Analyze GC content variation across genome"""

    logger.info(f"\n{'='*80}")
    logger.info(f"GENOME GC CONTENT ANALYSIS")
    logger.info(f"{'='*80}")

    genome_seq = ""
    for record in SeqIO.parse(genome_file, "fasta"):
        genome_seq += str(record.seq).upper()

    # Overall GC
    gc_count = genome_seq.count('G') + genome_seq.count('C')
    overall_gc = (gc_count / len(genome_seq)) * 100
    logger.info(f"Whole genome GC: {overall_gc:.2f}%")

    # Window analysis
    gc_windows = calculate_gc_content_windows(genome_seq, window_size=10000)
    gc_values = [gc for _, gc in gc_windows]

    logger.info(f"\nGC content variation (10kb windows):")
    logger.info(f"  Mean: {sum(gc_values)/len(gc_values):.2f}%")
    logger.info(f"  Min: {min(gc_values):.2f}%")
    logger.info(f"  Max: {max(gc_values):.2f}%")
    logger.info(f"  Range: {max(gc_values) - min(gc_values):.2f}%")

    # Find GC-poor regions (potential coverage gaps)
    gc_poor_windows = [(pos, gc) for pos, gc in gc_windows if gc < 25]
    logger.info(f"\nGC-poor regions (<25% GC): {len(gc_poor_windows)} windows ({len(gc_poor_windows)/len(gc_windows)*100:.1f}%)")

    if gc_poor_windows:
        logger.info(f"Lowest GC regions:")
        for pos, gc in sorted(gc_poor_windows, key=lambda x: x[1])[:5]:
            logger.info(f"  Position {pos:,}: {gc:.1f}% GC")

    return genome_seq


def main():
    base_dir = Path(__file__).parent
    genome_file = Path("/Users/andreassjodin/Code/swga-dev/test/francisella_GCF_000008985.1_ASM898v1_genomic.fna")

    logger.info("="*80)
    logger.info("FRANCISELLA BINDING SITE DISTRIBUTION ANALYSIS")
    logger.info("="*80)

    # First analyze genome GC variation
    genome_seq = analyze_gc_correlation(genome_file)

    # Analyze each configuration
    configs = [
        ("Config1_Long_Lenient", base_dir / "francisella_results/Config1_Long_Lenient"),
        ("Config4_MediumLong_Moderate", base_dir / "francisella_results/Config4_MediumLong_Moderate"),
        ("Config5_VeryLong_Moderate", base_dir / "francisella_results/Config5_VeryLong_Moderate"),
        ("Config7_SuperLong_Optimized", base_dir / "francisella_results/Config7_SuperLong_Optimized"),
    ]

    all_results = []
    for config_name, config_dir in configs:
        result = analyze_config(config_name, config_dir, genome_file, genome_seq)
        all_results.append(result)

    # Summary comparison
    logger.info(f"\n{'='*80}")
    logger.info(f"COMPARATIVE SUMMARY")
    logger.info(f"{'='*80}")
    logger.info(f"\n{'Config':<30} {'Primers':>8} {'Sites':>8} {'Coverage':>10} {'Clustering':>12}")
    logger.info(f"{'-'*80}")

    for r in all_results:
        logger.info(f"{r['config_name']:<30} {r['n_primers']:>8} {r['total_binding_sites']:>8} "
                   f"{r['coverage_percent']:>9.1f}% {r['clustering_coefficient']:>12.3f}")

    # Save detailed results
    output_file = base_dir / "binding_distribution_analysis.json"
    with open(output_file, 'w') as f:
        json.dump(all_results, f, indent=2)

    logger.info(f"\n{'='*80}")
    logger.info(f"Analysis complete. Results saved to: {output_file}")
    logger.info(f"{'='*80}")


if __name__ == "__main__":
    main()
