#!/usr/bin/env python3
"""
SWGA Simulation Analysis and Metrics

Provides detailed analysis of simulation results including:
- Coverage uniformity analysis
- Gap identification and characterization
- GC-bias detection
- Primer contribution analysis
- Recommendations for improvement

Usage:
    from neoswga.core.simulation_analysis import SimulationAnalyzer

    analyzer = SimulationAnalyzer(simulation_result)
    analysis = analyzer.analyze()
    recommendations = analyzer.generate_recommendations()
"""

import numpy as np
from typing import List, Dict, Tuple, Optional
from dataclasses import dataclass
import logging

logger = logging.getLogger(__name__)


@dataclass
class CoverageAnalysis:
    """Detailed coverage analysis"""
    overall_coverage: float
    uniformity: float
    gaps: List[Dict]
    largest_gap: Optional[Dict]
    mean_inter_site_distance: float
    max_inter_site_distance: int
    bins_covered: int
    bins_total: int


@dataclass
class SpecificityAnalysis:
    """Specificity and enrichment analysis"""
    enrichment: float
    target_sites: int
    background_sites: int
    specificity_score: float
    target_density: float  # Sites per Mbp
    background_density: float


@dataclass
class PrimerContribution:
    """Individual primer contribution to overall performance"""
    primer: str
    target_sites: int
    background_sites: int
    specificity: float
    bins_covered: int
    unique_bins: int  # Bins only this primer covers
    contribution_score: float


@dataclass
class ComprehensiveAnalysis:
    """Complete simulation analysis"""
    coverage: CoverageAnalysis
    specificity: SpecificityAnalysis
    primer_contributions: List[PrimerContribution]
    quality_score: float
    recommendations: List[str]
    issues: List[str]


class SimulationAnalyzer:
    """
    Analyzes SWGA simulation results to identify strengths, weaknesses, and improvements.
    """

    def __init__(self, result, fg_positions: Dict, bg_positions: Dict,
                 fg_length: int, bg_length: int, bin_size: int = 10000):
        """
        Initialize analyzer.

        Args:
            result: SimulationResult from SwgaSimulator
            fg_positions: Dict of primer positions for target genome
            bg_positions: Dict of primer positions for background genome
            fg_length: Target genome length
            bg_length: Background genome length
            bin_size: Bin size for coverage analysis (default 10kb)
        """
        self.result = result
        self.fg_positions = fg_positions
        self.bg_positions = bg_positions
        self.fg_length = fg_length
        self.bg_length = bg_length
        self.bin_size = bin_size

    def analyze(self) -> ComprehensiveAnalysis:
        """
        Perform comprehensive analysis.

        Returns:
            ComprehensiveAnalysis with detailed metrics
        """
        logger.info("Performing comprehensive simulation analysis...")

        # Coverage analysis
        coverage = self._analyze_coverage()

        # Specificity analysis
        specificity = self._analyze_specificity()

        # Primer contributions
        primer_contributions = self._analyze_primer_contributions()

        # Overall quality score
        quality_score = self._calculate_quality_score(coverage, specificity)

        # Generate recommendations
        recommendations = self._generate_recommendations(coverage, specificity, primer_contributions)

        # Identify issues
        issues = self._identify_issues(coverage, specificity, primer_contributions)

        return ComprehensiveAnalysis(
            coverage=coverage,
            specificity=specificity,
            primer_contributions=primer_contributions,
            quality_score=quality_score,
            recommendations=recommendations,
            issues=issues
        )

    def _analyze_coverage(self) -> CoverageAnalysis:
        """Analyze coverage patterns"""
        # Get all binding positions
        all_positions = []
        for primer, pos_data in self.fg_positions.items():
            all_positions.extend(pos_data['+'])
            all_positions.extend(pos_data['-'])

        all_positions = sorted(all_positions)

        # Calculate inter-site distances
        if len(all_positions) > 1:
            distances = [all_positions[i+1] - all_positions[i]
                        for i in range(len(all_positions)-1)]
            mean_distance = np.mean(distances)
            max_distance = max(distances)
        else:
            mean_distance = 0
            max_distance = 0

        # Bin coverage
        n_bins = (self.fg_length + self.bin_size - 1) // self.bin_size
        covered_bins = set()

        for pos in all_positions:
            bin_idx = int(pos) // self.bin_size
            if 0 <= bin_idx < n_bins:
                covered_bins.add(bin_idx)

        # Identify gaps
        gaps = self._identify_detailed_gaps(covered_bins, n_bins)
        largest_gap = max(gaps, key=lambda g: g['length']) if gaps else None

        return CoverageAnalysis(
            overall_coverage=len(covered_bins) / n_bins if n_bins > 0 else 0,
            uniformity=self.result.target_uniformity,
            gaps=gaps,
            largest_gap=largest_gap,
            mean_inter_site_distance=mean_distance,
            max_inter_site_distance=max_distance,
            bins_covered=len(covered_bins),
            bins_total=n_bins
        )

    def _analyze_specificity(self) -> SpecificityAnalysis:
        """Analyze specificity and enrichment"""
        # Count total sites
        fg_sites = sum(len(pos['+']) + len(pos['-'])
                      for pos in self.fg_positions.values())
        bg_sites = sum(len(pos['+']) + len(pos['-'])
                      for pos in self.bg_positions.values())

        # Calculate densities (sites per Mbp)
        fg_density = fg_sites / (self.fg_length / 1e6)
        bg_density = bg_sites / (self.bg_length / 1e6)

        return SpecificityAnalysis(
            enrichment=self.result.enrichment,
            target_sites=fg_sites,
            background_sites=bg_sites,
            specificity_score=self.result.specificity_score,
            target_density=fg_density,
            background_density=bg_density
        )

    def _analyze_primer_contributions(self) -> List[PrimerContribution]:
        """Analyze individual primer contributions"""
        contributions = []

        n_bins = (self.fg_length + self.bin_size - 1) // self.bin_size

        # Get global covered bins
        all_covered_bins = set()
        for pos_data in self.fg_positions.values():
            for pos in pos_data['+']:
                bin_idx = int(pos) // self.bin_size
                if 0 <= bin_idx < n_bins:
                    all_covered_bins.add(bin_idx)
            for pos in pos_data['-']:
                bin_idx = int(pos) // self.bin_size
                if 0 <= bin_idx < n_bins:
                    all_covered_bins.add(bin_idx)

        # Analyze each primer
        for primer in self.fg_positions.keys():
            fg_pos = self.fg_positions[primer]
            bg_pos = self.bg_positions.get(primer, {'+': [], '-': []})

            # Count sites
            fg_sites = len(fg_pos['+']) + len(fg_pos['-'])
            bg_sites = len(bg_pos['+']) + len(bg_pos['-'])

            # Calculate specificity
            if bg_sites > 0:
                specificity = fg_sites / bg_sites
            else:
                specificity = float('inf')

            # Bins covered by this primer
            primer_bins = set()
            for pos in fg_pos['+']:
                bin_idx = int(pos) // self.bin_size
                if 0 <= bin_idx < n_bins:
                    primer_bins.add(bin_idx)
            for pos in fg_pos['-']:
                bin_idx = int(pos) // self.bin_size
                if 0 <= bin_idx < n_bins:
                    primer_bins.add(bin_idx)

            # Unique bins (only this primer covers)
            other_bins = set()
            for other_primer, other_pos in self.fg_positions.items():
                if other_primer != primer:
                    for pos in other_pos['+']:
                        bin_idx = int(pos) // self.bin_size
                        if 0 <= bin_idx < n_bins:
                            other_bins.add(bin_idx)
                    for pos in other_pos['-']:
                        bin_idx = int(pos) // self.bin_size
                        if 0 <= bin_idx < n_bins:
                            other_bins.add(bin_idx)

            unique_bins = primer_bins - other_bins

            # Contribution score (combination of coverage and specificity)
            coverage_contribution = len(primer_bins) / n_bins if n_bins > 0 else 0
            spec_contribution = min(specificity / 100, 1.0)  # Normalize to 0-1
            contribution_score = 0.6 * coverage_contribution + 0.4 * spec_contribution

            contributions.append(PrimerContribution(
                primer=primer,
                target_sites=fg_sites,
                background_sites=bg_sites,
                specificity=specificity,
                bins_covered=len(primer_bins),
                unique_bins=len(unique_bins),
                contribution_score=contribution_score
            ))

        # Sort by contribution score
        contributions.sort(key=lambda c: c.contribution_score, reverse=True)

        return contributions

    def _identify_detailed_gaps(self, covered_bins: set, n_bins: int) -> List[Dict]:
        """Identify gaps with detailed characterization"""
        if len(covered_bins) == 0:
            return []

        sorted_bins = sorted(covered_bins)
        gaps = []

        for i in range(len(sorted_bins) - 1):
            gap_size = sorted_bins[i+1] - sorted_bins[i] - 1

            if gap_size > 0:
                start_pos = (sorted_bins[i] + 1) * self.bin_size
                end_pos = sorted_bins[i+1] * self.bin_size
                gap_length = end_pos - start_pos

                # Characterize gap severity
                if gap_length > 100000:  # >100kb
                    severity = 'critical'
                elif gap_length > 50000:  # >50kb
                    severity = 'high'
                elif gap_length > 20000:  # >20kb
                    severity = 'medium'
                else:
                    severity = 'low'

                gaps.append({
                    'start': start_pos,
                    'end': end_pos,
                    'length': gap_length,
                    'bins': gap_size,
                    'severity': severity
                })

        return gaps

    def _calculate_quality_score(self, coverage: CoverageAnalysis,
                                 specificity: SpecificityAnalysis) -> float:
        """Calculate overall quality score (0-1)"""
        # Coverage component
        coverage_score = coverage.overall_coverage

        # Uniformity component
        uniformity_score = coverage.uniformity

        # Specificity component (log scale)
        spec_score = min(np.log10(specificity.enrichment + 1) / 3, 1.0)

        # Gap penalty
        gap_penalty = 0
        if coverage.largest_gap and coverage.largest_gap['length'] > 100000:
            gap_penalty = 0.1
        elif coverage.largest_gap and coverage.largest_gap['length'] > 50000:
            gap_penalty = 0.05

        # Weighted average
        quality = (
            0.35 * coverage_score +
            0.25 * uniformity_score +
            0.40 * spec_score
        ) - gap_penalty

        return max(0.0, min(1.0, quality))

    def _generate_recommendations(self, coverage: CoverageAnalysis,
                                  specificity: SpecificityAnalysis,
                                  primers: List[PrimerContribution]) -> List[str]:
        """Generate actionable recommendations"""
        recommendations = []

        # Coverage recommendations
        if coverage.overall_coverage < 0.4:
            recommendations.append(
                f"Low coverage ({coverage.overall_coverage:.1%}). "
                f"Consider adding {5-len(primers)} more primers to improve coverage."
            )
        elif coverage.overall_coverage > 0.7:
            recommendations.append(
                "Excellent coverage! Primer set is well-distributed."
            )

        # Uniformity recommendations
        if coverage.uniformity < 0.5:
            recommendations.append(
                f"Low uniformity ({coverage.uniformity:.1%}). "
                "Primers are clustering in specific regions. "
                "Add primers in under-covered regions."
            )

        # Gap recommendations
        critical_gaps = [g for g in coverage.gaps if g['severity'] == 'critical']
        if critical_gaps:
            recommendations.append(
                f"Found {len(critical_gaps)} critical gaps (>100kb). "
                f"Largest gap: {critical_gaps[0]['length']:,} bp at position {critical_gaps[0]['start']:,}. "
                "Add primers targeting these regions."
            )

        # Specificity recommendations
        if specificity.enrichment < 10:
            recommendations.append(
                f"Low specificity ({specificity.enrichment:.1f}×). "
                "Primers bind too frequently to background genome. "
                "Consider stricter frequency filters."
            )
        elif specificity.enrichment > 1000:
            recommendations.append(
                f"Excellent specificity ({specificity.enrichment:.0f}×)! "
                "Primers are highly specific to target genome."
            )

        # Primer contribution recommendations
        weak_primers = [p for p in primers if p.contribution_score < 0.3]
        if weak_primers:
            recommendations.append(
                f"{len(weak_primers)} primers have low contribution scores. "
                f"Consider replacing: {', '.join(p.primer for p in weak_primers[:3])}"
            )

        return recommendations

    def _identify_issues(self, coverage: CoverageAnalysis,
                        specificity: SpecificityAnalysis,
                        primers: List[PrimerContribution]) -> List[str]:
        """Identify specific issues with primer set"""
        issues = []

        # Coverage issues
        if coverage.overall_coverage < 0.3:
            issues.append("CRITICAL: Very low coverage (<30%)")

        if coverage.largest_gap and coverage.largest_gap['length'] > 100000:
            issues.append(
                f"CRITICAL: Large gap detected at position {coverage.largest_gap['start']:,} "
                f"({coverage.largest_gap['length']:,} bp)"
            )

        # Specificity issues
        if specificity.enrichment < 5:
            issues.append("CRITICAL: Poor specificity (<5× enrichment)")

        # High background binding
        if specificity.background_density > 10:  # >10 sites per Mbp
            issues.append(
                f"WARNING: High background binding ({specificity.background_density:.1f} sites/Mbp)"
            )

        # Primer redundancy
        redundant = sum(1 for p in primers if p.unique_bins == 0)
        if redundant > 3:
            issues.append(
                f"WARNING: {redundant} primers are redundant (cover no unique bins)"
            )

        return issues


def print_analysis_report(analysis: ComprehensiveAnalysis):
    """Print comprehensive analysis report"""
    print("\n" + "=" * 80)
    print("COMPREHENSIVE SIMULATION ANALYSIS")
    print("=" * 80)

    # Coverage section
    print("\n--- COVERAGE ANALYSIS ---")
    print(f"Overall coverage: {analysis.coverage.overall_coverage:.1%}")
    print(f"Uniformity: {analysis.coverage.uniformity:.1%}")
    print(f"Bins covered: {analysis.coverage.bins_covered}/{analysis.coverage.bins_total}")
    print(f"Mean inter-site distance: {analysis.coverage.mean_inter_site_distance:,.0f} bp")
    print(f"Max inter-site distance: {analysis.coverage.max_inter_site_distance:,} bp")
    print(f"Number of gaps: {len(analysis.coverage.gaps)}")
    if analysis.coverage.largest_gap:
        print(f"Largest gap: {analysis.coverage.largest_gap['length']:,} bp "
              f"({analysis.coverage.largest_gap['severity']} severity)")

    # Specificity section
    print("\n--- SPECIFICITY ANALYSIS ---")
    print(f"Enrichment: {analysis.specificity.enrichment:.1f}×")
    print(f"Target sites: {analysis.specificity.target_sites}")
    print(f"Background sites: {analysis.specificity.background_sites}")
    print(f"Target density: {analysis.specificity.target_density:.1f} sites/Mbp")
    print(f"Background density: {analysis.specificity.background_density:.1f} sites/Mbp")

    # Primer contributions
    print("\n--- TOP 5 PRIMER CONTRIBUTIONS ---")
    for i, primer in enumerate(analysis.primer_contributions[:5], 1):
        print(f"{i}. {primer.primer}")
        print(f"   Sites: {primer.target_sites} target, {primer.background_sites} background")
        print(f"   Specificity: {primer.specificity:.1f}×")
        print(f"   Bins covered: {primer.bins_covered} ({primer.unique_bins} unique)")
        print(f"   Contribution score: {primer.contribution_score:.2f}")

    # Quality and recommendations
    print("\n--- OVERALL ASSESSMENT ---")
    print(f"Quality score: {analysis.quality_score:.2f}/1.00")

    if analysis.issues:
        print("\nISSUES IDENTIFIED:")
        for issue in analysis.issues:
            print(f"  • {issue}")

    print("\nRECOMMENDATIONS:")
    for rec in analysis.recommendations:
        print(f"  • {rec}")

    print("\n" + "=" * 80)
