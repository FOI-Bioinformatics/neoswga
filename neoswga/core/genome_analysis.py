"""
Genome Composition Analysis Utilities.

This module provides functions for analyzing genome composition,
particularly GC content calculation for genome-adaptive QA thresholds.

Author: NeoSWGA Development Team
Date: November 2025
Version: 3.5 - Genome-Adaptive QA System
"""

import logging
from pathlib import Path
from typing import Union, Dict, Optional
from Bio import SeqIO

logger = logging.getLogger(__name__)


def calculate_genome_gc(genome_path: Union[str, Path]) -> float:
    """
    Calculate GC content from genome FASTA file.

    Args:
        genome_path: Path to genome FASTA file

    Returns:
        GC content as fraction (0-1)

    Examples:
        >>> gc = calculate_genome_gc("francisella.fasta")
        >>> print(f"Francisella GC: {gc:.1%}")
        Francisella GC: 32.3%

    Raises:
        FileNotFoundError: If genome file doesn't exist
        ValueError: If file is empty or invalid FASTA
    """
    genome_path = Path(genome_path)

    if not genome_path.exists():
        raise FileNotFoundError(f"Genome file not found: {genome_path}")

    total_length = 0
    gc_count = 0

    try:
        with open(genome_path) as fh:
            for record in SeqIO.parse(fh, "fasta"):
                seq = str(record.seq).upper()
                total_length += len(seq)
                gc_count += seq.count('G') + seq.count('C')

        if total_length == 0:
            raise ValueError(f"Empty or invalid FASTA file: {genome_path}")

        gc_fraction = gc_count / total_length

        logger.info(
            f"Genome GC content: {gc_fraction:.2%} "
            f"({gc_count:,} GC / {total_length:,} bp) "
            f"[{genome_path.name}]"
        )

        return gc_fraction

    except Exception as e:
        logger.error(f"Error reading genome file {genome_path}: {e}")
        raise


def calculate_genome_stats(genome_path: Union[str, Path]) -> Dict[str, float]:
    """
    Calculate comprehensive genome composition statistics.

    Args:
        genome_path: Path to genome FASTA file

    Returns:
        Dictionary with composition statistics:
        - gc_content: Overall GC content (0-1)
        - at_content: Overall AT content (0-1)
        - length: Total genome length (bp)
        - n_contigs: Number of contigs/chromosomes
        - gc_std: Standard deviation of GC across 10kb windows

    Example:
        >>> stats = calculate_genome_stats("francisella.fasta")
        >>> print(f"GC: {stats['gc_content']:.1%}, Length: {stats['length']:,} bp")
        GC: 32.3%, Length: 1,892,775 bp
    """
    genome_path = Path(genome_path)

    if not genome_path.exists():
        raise FileNotFoundError(f"Genome file not found: {genome_path}")

    total_length = 0
    gc_count = 0
    n_contigs = 0
    window_gcs = []

    # Analyze genome
    for record in SeqIO.parse(genome_path, "fasta"):
        seq = str(record.seq).upper()
        contig_length = len(seq)
        total_length += contig_length
        gc_count += seq.count('G') + seq.count('C')
        n_contigs += 1

        # Calculate GC in 10kb windows for standard deviation
        window_size = 10000
        for i in range(0, contig_length - window_size, window_size):
            window = seq[i:i+window_size]
            window_gc = (window.count('G') + window.count('C')) / len(window)
            window_gcs.append(window_gc)

    if total_length == 0:
        raise ValueError(f"Empty or invalid FASTA file: {genome_path}")

    gc_fraction = gc_count / total_length
    at_fraction = 1.0 - gc_fraction

    # Calculate GC standard deviation
    import numpy as np
    gc_std = np.std(window_gcs) if window_gcs else 0.0

    stats = {
        'gc_content': gc_fraction,
        'at_content': at_fraction,
        'length': total_length,
        'n_contigs': n_contigs,
        'gc_std': gc_std
    }

    logger.info(
        f"Genome stats [{genome_path.name}]: "
        f"GC={gc_fraction:.2%}, Length={total_length:,} bp, "
        f"Contigs={n_contigs}, GC_std={gc_std:.4f}"
    )

    return stats


def get_gc_class(gc_content: float) -> str:
    """
    Classify genome by GC content.

    Args:
        gc_content: GC fraction (0-1)

    Returns:
        GC class: 'extreme_at', 'at_rich', 'balanced', 'gc_rich', or 'extreme_gc'

    Classification:
        - Extreme AT-rich: <25% GC (e.g., Plasmodium 19%)
        - AT-rich: 25-40% GC (e.g., Francisella 32%, Wolbachia 34%)
        - Balanced: 40-60% GC (e.g., E. coli 50%)
        - GC-rich: 60-70% GC (e.g., Mycobacterium 66%, Burkholderia 67%)
        - Extreme GC-rich: >70% GC (e.g., Streptomyces 72%)

    Examples:
        >>> get_gc_class(0.19)
        'extreme_at'
        >>> get_gc_class(0.32)
        'at_rich'
        >>> get_gc_class(0.50)
        'balanced'
        >>> get_gc_class(0.67)
        'gc_rich'
        >>> get_gc_class(0.72)
        'extreme_gc'
    """
    if gc_content < 0.25:
        return 'extreme_at'
    elif gc_content < 0.40:
        return 'at_rich'
    elif gc_content < 0.60:
        return 'balanced'
    elif gc_content < 0.70:
        return 'gc_rich'
    else:
        return 'extreme_gc'


def recommend_adaptive_qa(genome_gc: float) -> Dict[str, Union[bool, str]]:
    """
    Recommend whether to use genome-adaptive QA based on GC content.

    Args:
        genome_gc: Genome GC fraction (0-1)

    Returns:
        Dictionary with recommendations:
        - use_adaptive: Boolean, whether to use adaptive QA
        - reason: String explaining recommendation
        - gc_class: Genome GC classification
        - expected_improvement: Estimated coverage improvement

    Examples:
        >>> rec = recommend_adaptive_qa(0.32)  # Francisella
        >>> print(rec['use_adaptive'])
        True
        >>> print(rec['reason'])
        AT-rich genome: adaptive QA strongly recommended
        >>> print(rec['expected_improvement'])
        +100% coverage improvement expected

        >>> rec = recommend_adaptive_qa(0.50)  # E. coli
        >>> print(rec['use_adaptive'])
        False
        >>> print(rec['reason'])
        Balanced genome: adaptive QA optional
    """
    gc_class = get_gc_class(genome_gc)

    if gc_class == 'extreme_at':
        return {
            'use_adaptive': True,
            'reason': 'Extreme AT-rich genome: adaptive QA CRITICAL for coverage',
            'gc_class': gc_class,
            'expected_improvement': '+150-200% coverage improvement expected'
        }
    elif gc_class == 'at_rich':
        return {
            'use_adaptive': True,
            'reason': 'AT-rich genome: adaptive QA strongly recommended',
            'gc_class': gc_class,
            'expected_improvement': '+75-100% coverage improvement expected'
        }
    elif gc_class == 'balanced':
        return {
            'use_adaptive': False,
            'reason': 'Balanced genome: adaptive QA optional (minimal benefit)',
            'gc_class': gc_class,
            'expected_improvement': '+10-20% coverage improvement expected'
        }
    elif gc_class == 'gc_rich':
        return {
            'use_adaptive': True,
            'reason': 'GC-rich genome: adaptive QA strongly recommended',
            'gc_class': gc_class,
            'expected_improvement': '+70-100% coverage improvement expected'
        }
    else:  # extreme_gc
        return {
            'use_adaptive': True,
            'reason': 'Extreme GC-rich genome: adaptive QA CRITICAL for coverage',
            'gc_class': gc_class,
            'expected_improvement': '+120-150% coverage improvement expected'
        }


def analyze_genome_for_qa(genome_path: Union[str, Path]) -> Dict:
    """
    Comprehensive genome analysis for QA parameter selection.

    Combines GC calculation, classification, and QA recommendations.

    Args:
        genome_path: Path to genome FASTA file

    Returns:
        Dictionary with complete analysis:
        - gc_content: GC fraction
        - genome_stats: Full composition statistics
        - gc_class: GC classification
        - adaptive_qa_recommendation: Recommendation dictionary
        - suggested_stringency: Recommended QA stringency level

    Example:
        >>> analysis = analyze_genome_for_qa("francisella.fasta")
        >>> print(f"GC: {analysis['gc_content']:.1%}")
        GC: 32.3%
        >>> print(f"Class: {analysis['gc_class']}")
        Class: at_rich
        >>> print(f"Use adaptive: {analysis['adaptive_qa_recommendation']['use_adaptive']}")
        Use adaptive: True
    """
    genome_path = Path(genome_path)

    # Calculate comprehensive stats
    stats = calculate_genome_stats(genome_path)
    gc_content = stats['gc_content']

    # Get classification
    gc_class = get_gc_class(gc_content)

    # Get recommendation
    recommendation = recommend_adaptive_qa(gc_content)

    # Suggest stringency based on genome size and GC variability
    if stats['length'] < 3_000_000:  # Small genome
        if stats['gc_std'] < 0.05:  # Low GC variability
            suggested_stringency = 'moderate'
        else:
            suggested_stringency = 'lenient'
    else:  # Large genome
        if stats['gc_std'] < 0.05:
            suggested_stringency = 'strict'
        else:
            suggested_stringency = 'moderate'

    return {
        'genome_path': str(genome_path),
        'gc_content': gc_content,
        'genome_stats': stats,
        'gc_class': gc_class,
        'adaptive_qa_recommendation': recommendation,
        'suggested_stringency': suggested_stringency
    }


def analyze_genome(
    genome_path: Union[str, Path],
    window_size: int = 1000,
    output_dir: Optional[Union[str, Path]] = None
) -> Dict:
    """
    Comprehensive genome suitability analysis for SWGA.

    Analyzes genome composition and provides recommendations for
    primer design parameters. Optionally saves report to output directory.

    Args:
        genome_path: Path to genome FASTA file
        window_size: Window size for GC profiling (default: 1000bp)
        output_dir: Optional directory for output files

    Returns:
        Dictionary with complete analysis:
        - genome_info: Basic genome statistics
        - gc_analysis: GC content analysis
        - swga_suitability: Suitability assessment for SWGA
        - parameter_recommendations: Suggested pipeline parameters
        - warnings: List of potential issues

    Example:
        >>> results = analyze_genome("target.fasta", output_dir="analysis/")
        >>> print(f"Suitability: {results['swga_suitability']['rating']}")
        Suitability: GOOD
    """
    import json
    import numpy as np

    genome_path = Path(genome_path)
    if not genome_path.exists():
        raise FileNotFoundError(f"Genome file not found: {genome_path}")

    # Create output directory if specified
    if output_dir:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

    # Collect analysis results
    results = {
        'genome_path': str(genome_path),
        'genome_info': {},
        'gc_analysis': {},
        'gc_profile': [],
        'swga_suitability': {},
        'parameter_recommendations': {},
        'warnings': []
    }

    # Parse genome and collect statistics
    total_length = 0
    gc_count = 0
    n_count = 0
    contigs = []
    window_gcs = []

    for record in SeqIO.parse(genome_path, "fasta"):
        seq = str(record.seq).upper()
        contig_length = len(seq)
        contig_gc = (seq.count('G') + seq.count('C')) / contig_length if contig_length > 0 else 0
        contig_n = seq.count('N')

        contigs.append({
            'id': record.id,
            'length': contig_length,
            'gc': contig_gc,
            'n_count': contig_n
        })

        total_length += contig_length
        gc_count += seq.count('G') + seq.count('C')
        n_count += contig_n

        # Calculate GC in windows for profiling
        for i in range(0, contig_length - window_size, window_size):
            window = seq[i:i+window_size]
            window_gc = (window.count('G') + window.count('C')) / len(window)
            window_gcs.append(window_gc)

    if total_length == 0:
        raise ValueError(f"Empty or invalid FASTA file: {genome_path}")

    # Calculate overall statistics
    gc_content = gc_count / total_length
    gc_class = get_gc_class(gc_content)
    gc_std = np.std(window_gcs) if window_gcs else 0.0
    gc_min = min(window_gcs) if window_gcs else gc_content
    gc_max = max(window_gcs) if window_gcs else gc_content

    # Genome info
    results['genome_info'] = {
        'name': genome_path.name,
        'total_length': total_length,
        'n_contigs': len(contigs),
        'n50': _calculate_n50([c['length'] for c in contigs]),
        'largest_contig': max(c['length'] for c in contigs) if contigs else 0,
        'n_bases': n_count,
        'n_fraction': n_count / total_length,
        'contigs': contigs[:10] if len(contigs) > 10 else contigs  # Top 10 only
    }

    # GC analysis
    results['gc_analysis'] = {
        'overall_gc': gc_content,
        'gc_class': gc_class,
        'gc_class_description': _gc_class_description(gc_class),
        'gc_std': gc_std,
        'gc_min': gc_min,
        'gc_max': gc_max,
        'gc_range': gc_max - gc_min,
        'window_size': window_size,
        'n_windows': len(window_gcs)
    }

    # GC profile (subsampled for output)
    if len(window_gcs) > 100:
        step = len(window_gcs) // 100
        results['gc_profile'] = window_gcs[::step]
    else:
        results['gc_profile'] = window_gcs

    # SWGA suitability assessment
    suitability_score = _calculate_suitability_score(results)
    results['swga_suitability'] = {
        'score': suitability_score,
        'rating': _score_to_rating(suitability_score),
        'factors': _get_suitability_factors(results)
    }

    # Parameter recommendations
    qa_rec = recommend_adaptive_qa(gc_content)
    results['parameter_recommendations'] = {
        'polymerase': 'equiphi29' if gc_class in ('gc_rich', 'extreme_gc') else 'phi29',
        'reaction_temp': 42.0 if gc_class in ('gc_rich', 'extreme_gc') else 30.0,
        'primer_length': _recommend_primer_length(gc_class),
        'gc_bounds': {
            'min': max(0.0, gc_content - 0.15),
            'max': min(1.0, gc_content + 0.15)
        },
        'additives': _recommend_additives(gc_class),
        'use_adaptive_qa': qa_rec['use_adaptive'],
        'adaptive_qa_reason': qa_rec['reason']
    }

    # Warnings
    if gc_content < 0.25:
        results['warnings'].append(
            f"Extreme AT-rich genome ({gc_content:.1%} GC) - limited primer candidates expected"
        )
    elif gc_content > 0.70:
        results['warnings'].append(
            f"Extreme GC-rich genome ({gc_content:.1%} GC) - high additive concentrations required"
        )

    if gc_std > 0.10:
        results['warnings'].append(
            f"High GC variability (std={gc_std:.3f}) - consider using adaptive GC filtering"
        )

    if n_count / total_length > 0.01:
        results['warnings'].append(
            f"High N content ({n_count/total_length:.1%}) - assembly quality may affect results"
        )

    if len(contigs) > 100:
        results['warnings'].append(
            f"Fragmented assembly ({len(contigs)} contigs) - some primer positions may be missed"
        )

    # Save results if output directory specified
    if output_dir:
        # Save JSON report
        report_path = output_dir / 'genome_analysis.json'
        with open(report_path, 'w') as f:
            json.dump(results, f, indent=2, default=str)

        # Save text summary
        summary_path = output_dir / 'genome_analysis_summary.txt'
        _write_text_summary(results, summary_path)

        logger.info(f"Analysis saved to {output_dir}")

    return results


def _calculate_n50(lengths: list) -> int:
    """Calculate N50 from contig lengths."""
    if not lengths:
        return 0
    sorted_lengths = sorted(lengths, reverse=True)
    total = sum(sorted_lengths)
    cumsum = 0
    for length in sorted_lengths:
        cumsum += length
        if cumsum >= total / 2:
            return length
    return sorted_lengths[-1]


def _gc_class_description(gc_class: str) -> str:
    """Get human-readable description for GC class."""
    descriptions = {
        'extreme_at': 'Extreme AT-rich (like Plasmodium, <25% GC)',
        'at_rich': 'AT-rich (like Francisella, Wolbachia, 25-40% GC)',
        'balanced': 'Balanced (like E. coli, 40-60% GC)',
        'gc_rich': 'GC-rich (like Mycobacterium, Burkholderia, 60-70% GC)',
        'extreme_gc': 'Extreme GC-rich (like Streptomyces, >70% GC)'
    }
    return descriptions.get(gc_class, gc_class)


def _calculate_suitability_score(results: Dict) -> float:
    """Calculate SWGA suitability score (0-100)."""
    score = 100.0

    gc = results['gc_analysis']['overall_gc']
    gc_std = results['gc_analysis']['gc_std']
    n_fraction = results['genome_info']['n_fraction']
    n_contigs = results['genome_info']['n_contigs']

    # Penalize extreme GC
    if gc < 0.25 or gc > 0.75:
        score -= 30
    elif gc < 0.35 or gc > 0.65:
        score -= 10

    # Penalize high GC variability
    if gc_std > 0.15:
        score -= 20
    elif gc_std > 0.10:
        score -= 10

    # Penalize high N content
    if n_fraction > 0.05:
        score -= 20
    elif n_fraction > 0.01:
        score -= 10

    # Penalize fragmentation
    if n_contigs > 500:
        score -= 15
    elif n_contigs > 100:
        score -= 5

    return max(0, min(100, score))


def _score_to_rating(score: float) -> str:
    """Convert suitability score to rating."""
    if score >= 90:
        return 'EXCELLENT'
    elif score >= 75:
        return 'GOOD'
    elif score >= 60:
        return 'MODERATE'
    elif score >= 40:
        return 'CHALLENGING'
    else:
        return 'DIFFICULT'


def _get_suitability_factors(results: Dict) -> list:
    """Get list of factors affecting suitability."""
    factors = []
    gc = results['gc_analysis']['overall_gc']
    gc_class = results['gc_analysis']['gc_class']

    if gc_class == 'balanced':
        factors.append({'factor': 'GC content', 'impact': 'positive', 'note': 'Balanced GC is optimal'})
    elif gc_class in ('at_rich', 'gc_rich'):
        factors.append({'factor': 'GC content', 'impact': 'neutral', 'note': 'May need adjusted parameters'})
    else:
        factors.append({'factor': 'GC content', 'impact': 'negative', 'note': 'Extreme GC requires special handling'})

    if results['gc_analysis']['gc_std'] < 0.05:
        factors.append({'factor': 'GC uniformity', 'impact': 'positive', 'note': 'Low GC variability'})
    elif results['gc_analysis']['gc_std'] > 0.10:
        factors.append({'factor': 'GC uniformity', 'impact': 'negative', 'note': 'High GC variability'})

    if results['genome_info']['n_contigs'] == 1:
        factors.append({'factor': 'Assembly', 'impact': 'positive', 'note': 'Complete chromosome'})
    elif results['genome_info']['n_contigs'] > 100:
        factors.append({'factor': 'Assembly', 'impact': 'negative', 'note': 'Fragmented assembly'})

    return factors


def _recommend_primer_length(gc_class: str) -> Dict:
    """Recommend primer length range based on GC class."""
    if gc_class in ('extreme_at', 'at_rich'):
        return {'min': 8, 'max': 12, 'optimal': 10}
    elif gc_class == 'balanced':
        return {'min': 10, 'max': 14, 'optimal': 12}
    else:  # gc_rich, extreme_gc
        return {'min': 12, 'max': 16, 'optimal': 14}


def _recommend_additives(gc_class: str) -> Dict:
    """Recommend additive concentrations based on GC class."""
    if gc_class == 'extreme_at':
        return {'betaine_m': 0.5, 'dmso_percent': 0.0, 'mg_conc': 2.5}
    elif gc_class == 'at_rich':
        return {'betaine_m': 0.8, 'dmso_percent': 3.0, 'mg_conc': 2.0}
    elif gc_class == 'balanced':
        return {'betaine_m': 1.0, 'dmso_percent': 3.0, 'mg_conc': 1.5}
    elif gc_class == 'gc_rich':
        return {'betaine_m': 1.5, 'dmso_percent': 5.0, 'mg_conc': 1.0}
    else:  # extreme_gc
        return {'betaine_m': 2.0, 'dmso_percent': 5.0, 'mg_conc': 0.5}


def _write_text_summary(results: Dict, output_path: Path) -> None:
    """Write text summary report."""
    with open(output_path, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write("GENOME SUITABILITY ANALYSIS FOR SWGA\n")
        f.write("=" * 70 + "\n\n")

        f.write(f"Genome: {results['genome_info']['name']}\n")
        f.write(f"Length: {results['genome_info']['total_length']:,} bp\n")
        f.write(f"Contigs: {results['genome_info']['n_contigs']}\n")
        f.write(f"N50: {results['genome_info']['n50']:,} bp\n\n")

        f.write("--- GC Content Analysis ---\n")
        f.write(f"Overall GC: {results['gc_analysis']['overall_gc']:.1%}\n")
        f.write(f"GC Class: {results['gc_analysis']['gc_class_description']}\n")
        f.write(f"GC Std Dev: {results['gc_analysis']['gc_std']:.4f}\n")
        f.write(f"GC Range: {results['gc_analysis']['gc_min']:.1%} - {results['gc_analysis']['gc_max']:.1%}\n\n")

        f.write("--- SWGA Suitability ---\n")
        f.write(f"Score: {results['swga_suitability']['score']:.0f}/100\n")
        f.write(f"Rating: {results['swga_suitability']['rating']}\n\n")

        f.write("--- Recommended Parameters ---\n")
        rec = results['parameter_recommendations']
        f.write(f"Polymerase: {rec['polymerase']}\n")
        f.write(f"Temperature: {rec['reaction_temp']}C\n")
        f.write(f"Primer length: {rec['primer_length']['min']}-{rec['primer_length']['max']} bp\n")
        f.write(f"GC bounds: {rec['gc_bounds']['min']:.2f}-{rec['gc_bounds']['max']:.2f}\n")
        f.write(f"Betaine: {rec['additives']['betaine_m']}M\n")
        f.write(f"DMSO: {rec['additives']['dmso_percent']}%\n\n")

        if results['warnings']:
            f.write("--- Warnings ---\n")
            for warning in results['warnings']:
                f.write(f"! {warning}\n")
            f.write("\n")

        f.write("=" * 70 + "\n")


if __name__ == "__main__":
    # Example usage
    import sys

    if len(sys.argv) > 1:
        genome_file = sys.argv[1]
        output_dir = sys.argv[2] if len(sys.argv) > 2 else None
    else:
        print("Usage: python genome_analysis.py <genome.fasta> [output_dir]")
        print("\nExample genomes to test:")
        print("  - Francisella (32% GC, AT-rich)")
        print("  - E. coli (50% GC, balanced)")
        print("  - Burkholderia (67% GC, GC-rich)")
        raise SystemExit(1)

    print("="*80)
    print("GENOME SUITABILITY ANALYSIS FOR SWGA")
    print("="*80)

    results = analyze_genome(genome_file, output_dir=output_dir)

    print(f"\nGenome: {results['genome_info']['name']}")
    print(f"Length: {results['genome_info']['total_length']:,} bp")
    print(f"Contigs: {results['genome_info']['n_contigs']}")
    print(f"GC content: {results['gc_analysis']['overall_gc']:.2%}")
    print(f"GC class: {results['gc_analysis']['gc_class_description']}")

    print(f"\nSWGA Suitability: {results['swga_suitability']['rating']} ({results['swga_suitability']['score']:.0f}/100)")

    if results['warnings']:
        print("\nWarnings:")
        for w in results['warnings']:
            print(f"  ! {w}")

    print("\n" + "="*80)
