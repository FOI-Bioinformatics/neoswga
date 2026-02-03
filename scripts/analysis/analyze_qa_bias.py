#!/usr/bin/env python3
"""
Analyze QA filter bias to understand why primers are GC-rich.

This script examines:
1. GC content distribution before vs after QA
2. Terminal GC content (3' end 5bp)
3. Which specific QA components cause GC bias
4. Correlation between GC content and QA scores
"""

import sys
from pathlib import Path
import logging
from typing import List, Dict, Tuple
import json

# Add neoswga to path
sys.path.insert(0, str(Path.cwd()))

from neoswga.core.three_prime_stability import ThreePrimeStabilityChecker
from neoswga.core.integrated_quality_scorer import IntegratedQualityScorer

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)


def analyze_candidates(candidates: List[str], config_name: str, mode: str, stringency: str) -> Dict:
    """Analyze GC distribution and QA scores for candidate primers"""

    logger.info(f"\n{'='*80}")
    logger.info(f"ANALYZING: {config_name}")
    logger.info(f"Mode: {mode}, Stringency: {stringency}")
    logger.info(f"{'='*80}")

    # Calculate GC metrics for all candidates
    gc_contents = []
    terminal_gc_contents = []  # Last 5bp

    for primer in candidates:
        seq = primer.upper()
        gc = (seq.count('G') + seq.count('C')) / len(seq) * 100
        gc_contents.append(gc)

        # Terminal 5bp GC (3' end)
        terminal = seq[-5:] if len(seq) >= 5 else seq
        term_gc = (terminal.count('G') + terminal.count('C')) / len(terminal) * 100
        terminal_gc_contents.append(term_gc)

    logger.info(f"\nCandidate Pool ({len(candidates)} primers):")
    logger.info(f"  GC content: {sum(gc_contents)/len(gc_contents):.1f}% mean")
    logger.info(f"  GC range: {min(gc_contents):.1f}% - {max(gc_contents):.1f}%")
    logger.info(f"  Terminal GC: {sum(terminal_gc_contents)/len(terminal_gc_contents):.1f}% mean")

    # Initialize QA components
    three_prime_checker = ThreePrimeStabilityChecker(mode=mode, stringency=stringency)
    integrated_scorer = IntegratedQualityScorer(mode=mode)

    # Score all candidates
    scores = []
    three_prime_scores = []
    passed = []
    failed = []

    for i, primer in enumerate(candidates):
        # Get integrated QA score
        qa_result = integrated_scorer.score_primer(primer)
        scores.append(qa_result.composite_score)

        # Get 3' stability score
        three_prime_result = three_prime_checker.evaluate(primer)
        three_prime_scores.append(three_prime_result.score)

        # Check if passed
        if qa_result.passed:
            passed.append((primer, qa_result.composite_score, gc_contents[i], terminal_gc_contents[i]))
        else:
            failed.append((primer, qa_result.composite_score, gc_contents[i], terminal_gc_contents[i]))

    logger.info(f"\nQA Results:")
    logger.info(f"  Passed: {len(passed)}/{len(candidates)} ({len(passed)/len(candidates)*100:.1f}%)")
    logger.info(f"  Failed: {len(failed)}/{len(candidates)} ({len(failed)/len(candidates)*100:.1f}%)")

    # Analyze passed primers
    if passed:
        passed_gc = [gc for _, _, gc, _ in passed]
        passed_term_gc = [tgc for _, _, _, tgc in passed]
        passed_scores = [score for _, score, _, _ in passed]

        logger.info(f"\nPassed Primers ({len(passed)}):")
        logger.info(f"  GC content: {sum(passed_gc)/len(passed_gc):.1f}% mean")
        logger.info(f"  GC range: {min(passed_gc):.1f}% - {max(passed_gc):.1f}%")
        logger.info(f"  Terminal GC: {sum(passed_term_gc)/len(passed_term_gc):.1f}% mean")
        logger.info(f"  Mean QA score: {sum(passed_scores)/len(passed_scores):.3f}")

        # GC bias
        gc_bias = (sum(passed_gc)/len(passed_gc)) - (sum(gc_contents)/len(gc_contents))
        term_gc_bias = (sum(passed_term_gc)/len(passed_term_gc)) - (sum(terminal_gc_contents)/len(terminal_gc_contents))

        logger.info(f"\nGC Bias:")
        logger.info(f"  Overall GC bias: +{gc_bias:.1f}% (passed primers are more GC-rich)")
        logger.info(f"  Terminal GC bias: +{term_gc_bias:.1f}% (3' ends are more GC-rich)")

    # Analyze failed primers
    if failed:
        failed_gc = [gc for _, _, gc, _ in failed]
        failed_scores = [score for _, score, _, _ in failed]

        logger.info(f"\nFailed Primers ({len(failed)}):")
        logger.info(f"  GC content: {sum(failed_gc)/len(failed_gc):.1f}% mean")
        logger.info(f"  GC range: {min(failed_gc):.1f}% - {max(failed_gc):.1f}%")
        logger.info(f"  Mean QA score: {sum(failed_scores)/len(failed_scores):.3f}")

    # Correlation analysis
    if len(scores) > 1:
        # Calculate correlation between GC and QA score
        mean_gc = sum(gc_contents) / len(gc_contents)
        mean_score = sum(scores) / len(scores)

        numerator = sum((gc - mean_gc) * (score - mean_score)
                       for gc, score in zip(gc_contents, scores))
        denom_gc = sum((gc - mean_gc) ** 2 for gc in gc_contents) ** 0.5
        denom_score = sum((score - mean_score) ** 2 for score in scores) ** 0.5

        if denom_gc > 0 and denom_score > 0:
            correlation = numerator / (denom_gc * denom_score)
            logger.info(f"\nCorrelation Analysis:")
            logger.info(f"  GC vs QA score correlation: {correlation:.3f}")
            if correlation > 0.3:
                logger.info(f"  → STRONG POSITIVE: Higher GC primers get higher QA scores")
            elif correlation > 0.1:
                logger.info(f"  → WEAK POSITIVE: Slight bias toward GC-rich primers")
            elif correlation < -0.1:
                logger.info(f"  → NEGATIVE: QA favors AT-rich primers")
            else:
                logger.info(f"  → NO CORRELATION: GC content independent of QA score")

    # Detailed component analysis for a sample of failed primers
    if failed and len(failed) >= 10:
        logger.info(f"\nDetailed Failure Analysis (Sample of 10 failed primers):")
        logger.info(f"{'Primer':<15} {'GC%':>6} {'TmGC%':>7} {'QA':>6} {'Reason'}")
        logger.info(f"{'-'*80}")

        for primer, score, gc, term_gc in failed[:10]:
            qa_result = integrated_scorer.score_primer(primer)
            three_prime_result = three_prime_checker.evaluate(primer)

            reason = "Unknown"
            if not three_prime_result.passed:
                reason = "3'stability"
            elif qa_result.dimer_score < 0.4:
                reason = "Dimer"
            elif qa_result.strand_bias_score < 0.4:
                reason = "StrandBias"
            elif qa_result.complexity_score < 0.4:
                reason = "Complexity"

            logger.info(f"{primer:<15} {gc:>6.1f} {term_gc:>7.1f} {score:>6.3f} {reason}")

    return {
        'config_name': config_name,
        'n_candidates': len(candidates),
        'n_passed': len(passed),
        'pass_rate': len(passed) / len(candidates),
        'candidate_mean_gc': sum(gc_contents) / len(gc_contents),
        'candidate_mean_terminal_gc': sum(terminal_gc_contents) / len(terminal_gc_contents),
        'passed_mean_gc': sum(passed_gc) / len(passed_gc) if passed else 0,
        'passed_mean_terminal_gc': sum(passed_term_gc) / len(passed_term_gc) if passed else 0,
        'gc_bias': (sum(passed_gc)/len(passed_gc) - sum(gc_contents)/len(gc_contents)) if passed else 0,
        'terminal_gc_bias': (sum(passed_term_gc)/len(passed_term_gc) - sum(terminal_gc_contents)/len(terminal_gc_contents)) if passed else 0
    }


def load_candidates_from_file(config_dir: Path) -> List[str]:
    """Load candidate primers before QA (from all_candidates.txt if exists)"""
    candidates_file = config_dir / "all_candidates.txt"

    if not candidates_file.exists():
        logger.warning(f"  No all_candidates.txt found, using primers.txt instead")
        return load_passed_from_file(config_dir)

    candidates = []
    with open(candidates_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if parts:
                candidates.append(parts[0])

    return candidates


def load_passed_from_file(config_dir: Path) -> List[str]:
    """Load passed primers from primers.txt"""
    primers_file = config_dir / "primers.txt"
    primers = []
    with open(primers_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if parts:
                primers.append(parts[0])
    return primers


def generate_synthetic_candidates(genome_gc: float, n_candidates: int = 1000, length: int = 12) -> List[str]:
    """Generate synthetic primers with GC matching genome"""
    import random

    bases_gc = ['G', 'C']
    bases_at = ['A', 'T']

    candidates = []
    for _ in range(n_candidates):
        seq = []
        for _ in range(length):
            if random.random() < (genome_gc / 100):
                seq.append(random.choice(bases_gc))
            else:
                seq.append(random.choice(bases_at))
        candidates.append(''.join(seq))

    return candidates


def main():
    base_dir = Path(__file__).parent

    logger.info("="*80)
    logger.info("QA FILTER GC BIAS ANALYSIS")
    logger.info("="*80)

    # First, generate synthetic candidates matching genome GC
    logger.info(f"\n{'='*80}")
    logger.info(f"SYNTHETIC TEST: Genome-matched GC candidates")
    logger.info(f"{'='*80}")
    logger.info(f"\nGenerating 1000 synthetic 12bp primers with 32% GC (matching Francisella)...")

    synthetic = generate_synthetic_candidates(genome_gc=32, n_candidates=1000, length=12)
    synthetic_result = analyze_candidates(synthetic, "Synthetic_32GC", mode='swga', stringency='moderate')

    # Now analyze real Francisella configs
    configs = [
        ("Config1_Long_Lenient", base_dir / "francisella_results/Config1_Long_Lenient", 'swga', 'lenient'),
        ("Config4_MediumLong_Moderate", base_dir / "francisella_results/Config4_MediumLong_Moderate", 'swga', 'moderate'),
    ]

    all_results = [synthetic_result]

    for config_name, config_dir, mode, stringency in configs:
        # Load primers (ideally would load candidates, but using passed for now)
        primers = load_passed_from_file(config_dir)

        # Generate synthetic candidates matching this primer length
        primer_length = int(sum(len(p) for p in primers) / len(primers))
        candidates = generate_synthetic_candidates(genome_gc=32, n_candidates=len(primers)*20, length=primer_length)

        result = analyze_candidates(candidates, config_name, mode, stringency)
        all_results.append(result)

    # Summary
    logger.info(f"\n{'='*80}")
    logger.info(f"SUMMARY: GC BIAS IN QA FILTERS")
    logger.info(f"{'='*80}")
    logger.info(f"\n{'Config':<30} {'PassRate':>10} {'InputGC':>10} {'PassGC':>10} {'Bias':>8}")
    logger.info(f"{'-'*80}")

    for r in all_results:
        logger.info(f"{r['config_name']:<30} {r['pass_rate']*100:>9.1f}% "
                   f"{r['candidate_mean_gc']:>9.1f}% {r['passed_mean_gc']:>9.1f}% "
                   f"{r['gc_bias']:>+7.1f}%")

    logger.info(f"\n{'='*80}")
    logger.info(f"CONCLUSION:")
    logger.info(f"{'='*80}")

    avg_bias = sum(r['gc_bias'] for r in all_results) / len(all_results)
    avg_term_bias = sum(r['terminal_gc_bias'] for r in all_results) / len(all_results)

    logger.info(f"Average GC bias: +{avg_bias:.1f}%")
    logger.info(f"Average terminal GC bias: +{avg_term_bias:.1f}%")

    if avg_bias > 5:
        logger.info(f"\nQA filters strongly favor GC-rich primers!")
        logger.info(f"This explains why 32% GC genome has 39-40% GC primers.")
    if avg_term_bias > 10:
        logger.info(f"\nTerminal stability requirements strongly favor GC-rich 3' ends!")
        logger.info(f"This is likely the primary cause of GC bias.")

    # Save results
    output_file = base_dir / "qa_bias_analysis.json"
    with open(output_file, 'w') as f:
        json.dump(all_results, f, indent=2)

    logger.info(f"\nDetailed results saved to: {output_file}")


if __name__ == "__main__":
    main()
