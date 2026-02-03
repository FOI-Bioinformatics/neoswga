#!/usr/bin/env python3
"""
Results analysis tool for benchmark suite.

Generates comprehensive analysis including:
- Method comparison tables
- Performance vs quality tradeoffs
- Best method recommendations per scenario
- Statistical analysis
- Detailed markdown report
"""

import sys
import csv
import json
import numpy as np
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Tuple
import logging

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)


class ResultsAnalyzer:
    """Analyze benchmark results"""

    def __init__(self, results_file: str = "./benchmarks/benchmark_results.csv"):
        """
        Initialize analyzer.

        Args:
            results_file: Path to benchmark results CSV
        """
        self.results_file = Path(results_file)

        if not self.results_file.exists():
            raise FileNotFoundError(f"Results file not found: {self.results_file}")

        # Load results
        self.results = self._load_results()

        logger.info(f"Loaded {len(self.results)} benchmark results")

    def _load_results(self) -> List[Dict]:
        """Load results from CSV"""
        results = []

        with open(self.results_file) as f:
            reader = csv.DictReader(f)
            for row in reader:
                # Convert numeric fields
                row['runtime'] = float(row['runtime'])
                row['peak_memory'] = float(row['peak_memory'])
                row['n_primers_selected'] = int(row['n_primers_selected'])
                row['coverage_score'] = float(row['coverage_score'])
                row['network_nodes'] = int(row['network_nodes'])
                row['network_edges'] = int(row['network_edges'])
                row['network_components'] = int(row['network_components'])
                row['background_binding'] = int(row['background_binding'])
                row['success'] = row['success'].lower() == 'true'
                row['use_cache'] = row['use_cache'].lower() == 'true'
                row['use_background_filter'] = row['use_background_filter'].lower() == 'true'
                row['gc_tolerance'] = float(row['gc_tolerance'])
                row['num_primers'] = int(row['num_primers'])

                results.append(row)

        return results

    def generate_report(self, output_file: str = "./benchmarks/analysis_report.md"):
        """
        Generate comprehensive analysis report.

        Args:
            output_file: Output markdown file
        """
        logger.info("Generating analysis report...")

        report_lines = []

        # Header
        report_lines.append("# Benchmark Analysis Report")
        report_lines.append("")
        report_lines.append(f"**Generated from**: {self.results_file}")
        report_lines.append(f"**Total tests**: {len(self.results)}")
        report_lines.append(f"**Successful tests**: {sum(1 for r in self.results if r['success'])}")
        report_lines.append("")
        report_lines.append("---")
        report_lines.append("")

        # Overall statistics
        report_lines.extend(self._analyze_overall_statistics())

        # Method comparison
        report_lines.extend(self._analyze_methods())

        # Genome-specific analysis
        report_lines.extend(self._analyze_by_genome())

        # GC tolerance impact
        report_lines.extend(self._analyze_gc_tolerance())

        # Primer count analysis
        report_lines.extend(self._analyze_primer_counts())

        # Cache impact
        report_lines.extend(self._analyze_cache_impact())

        # Recommendations
        report_lines.extend(self._generate_recommendations())

        # Save report
        output_path = Path(output_file)
        with open(output_path, 'w') as f:
            f.write('\n'.join(report_lines))

        logger.info(f"Report saved to: {output_path}")

        return report_lines

    def _analyze_overall_statistics(self) -> List[str]:
        """Generate overall statistics section"""
        lines = []
        lines.append("## Overall Statistics")
        lines.append("")

        successful = [r for r in self.results if r['success']]

        if not successful:
            lines.append("No successful tests to analyze.")
            return lines

        # Runtime stats
        runtimes = [r['runtime'] for r in successful]
        lines.append("### Runtime")
        lines.append(f"- Average: {np.mean(runtimes):.2f}s")
        lines.append(f"- Median: {np.median(runtimes):.2f}s")
        lines.append(f"- Min: {np.min(runtimes):.2f}s")
        lines.append(f"- Max: {np.max(runtimes):.2f}s")
        lines.append(f"- Std dev: {np.std(runtimes):.2f}s")
        lines.append("")

        # Coverage stats
        coverages = [r['coverage_score'] for r in successful]
        lines.append("### Coverage Score")
        lines.append(f"- Average: {np.mean(coverages):.4f}")
        lines.append(f"- Median: {np.median(coverages):.4f}")
        lines.append(f"- Min: {np.min(coverages):.4f}")
        lines.append(f"- Max: {np.max(coverages):.4f}")
        lines.append("")

        lines.append("---")
        lines.append("")

        return lines

    def _analyze_methods(self) -> List[str]:
        """Analyze performance by optimization method"""
        lines = []
        lines.append("## Method Comparison")
        lines.append("")

        successful = [r for r in self.results if r['success']]

        # Group by method
        by_method = defaultdict(list)
        for r in successful:
            by_method[r['optimization_method']].append(r)

        # Create comparison table
        lines.append("| Method | Tests | Avg Runtime (s) | Avg Coverage | Avg Primers | Success Rate |")
        lines.append("|--------|-------|-----------------|--------------|-------------|--------------|")

        for method in sorted(by_method.keys()):
            results = by_method[method]
            total_tests = sum(1 for r in self.results if r['optimization_method'] == method)

            avg_runtime = np.mean([r['runtime'] for r in results])
            avg_coverage = np.mean([r['coverage_score'] for r in results])
            avg_primers = np.mean([r['n_primers_selected'] for r in results])
            success_rate = len(results) / total_tests * 100

            lines.append(f"| {method} | {len(results)} | {avg_runtime:.2f} | "
                        f"{avg_coverage:.4f} | {avg_primers:.1f} | {success_rate:.1f}% |")

        lines.append("")

        # Best method for speed
        fastest_method = min(by_method.items(), key=lambda x: np.mean([r['runtime'] for r in x[1]]))
        lines.append(f"**Fastest method**: {fastest_method[0]} "
                    f"({np.mean([r['runtime'] for r in fastest_method[1]]):.2f}s average)")
        lines.append("")

        # Best method for quality
        best_coverage = max(by_method.items(), key=lambda x: np.mean([r['coverage_score'] for r in x[1]]))
        lines.append(f"**Best coverage**: {best_coverage[0]} "
                    f"({np.mean([r['coverage_score'] for r in best_coverage[1]]):.4f} average)")
        lines.append("")

        lines.append("---")
        lines.append("")

        return lines

    def _analyze_by_genome(self) -> List[str]:
        """Analyze performance by genome"""
        lines = []
        lines.append("## Genome-Specific Analysis")
        lines.append("")

        successful = [r for r in self.results if r['success']]

        for genome in ['francisella', 'burkholderia']:
            genome_results = [r for r in successful if r['genome'] == genome]

            if not genome_results:
                continue

            lines.append(f"### {genome.capitalize()}")
            lines.append("")

            # Best method for this genome
            by_method = defaultdict(list)
            for r in genome_results:
                by_method[r['optimization_method']].append(r)

            best_method = max(by_method.items(),
                            key=lambda x: np.mean([r['coverage_score'] for r in x[1]]))

            lines.append(f"**Best method**: {best_method[0]} "
                        f"({np.mean([r['coverage_score'] for r in best_method[1]]):.4f} coverage)")
            lines.append("")

            # Statistics
            runtimes = [r['runtime'] for r in genome_results]
            coverages = [r['coverage_score'] for r in genome_results]

            lines.append(f"- Tests: {len(genome_results)}")
            lines.append(f"- Avg runtime: {np.mean(runtimes):.2f}s")
            lines.append(f"- Avg coverage: {np.mean(coverages):.4f}")
            lines.append("")

        lines.append("---")
        lines.append("")

        return lines

    def _analyze_gc_tolerance(self) -> List[str]:
        """Analyze impact of GC tolerance"""
        lines = []
        lines.append("## GC Tolerance Impact")
        lines.append("")

        successful = [r for r in self.results if r['success']]

        # Group by GC tolerance
        by_gc = defaultdict(list)
        for r in successful:
            by_gc[r['gc_tolerance']].append(r)

        lines.append("| GC Tolerance | Tests | Avg Coverage | Avg Runtime (s) |")
        lines.append("|--------------|-------|--------------|-----------------|")

        for gc_tol in sorted(by_gc.keys()):
            results = by_gc[gc_tol]
            avg_coverage = np.mean([r['coverage_score'] for r in results])
            avg_runtime = np.mean([r['runtime'] for r in results])

            lines.append(f"| {gc_tol:.2f} | {len(results)} | {avg_coverage:.4f} | {avg_runtime:.2f} |")

        lines.append("")

        # Recommendation
        best_gc = max(by_gc.items(), key=lambda x: np.mean([r['coverage_score'] for r in x[1]]))
        lines.append(f"**Recommended GC tolerance**: {best_gc[0]:.2f} "
                    f"(best average coverage: {np.mean([r['coverage_score'] for r in best_gc[1]]):.4f})")
        lines.append("")

        lines.append("---")
        lines.append("")

        return lines

    def _analyze_primer_counts(self) -> List[str]:
        """Analyze impact of primer count"""
        lines = []
        lines.append("## Primer Count Analysis")
        lines.append("")

        successful = [r for r in self.results if r['success']]

        # Group by primer count
        by_count = defaultdict(list)
        for r in successful:
            by_count[r['num_primers']].append(r)

        lines.append("| Num Primers | Tests | Avg Coverage | Avg Runtime (s) | Coverage/Primer |")
        lines.append("|-------------|-------|--------------|-----------------|-----------------|")

        for count in sorted(by_count.keys()):
            results = by_count[count]
            avg_coverage = np.mean([r['coverage_score'] for r in results])
            avg_runtime = np.mean([r['runtime'] for r in results])
            coverage_per_primer = avg_coverage / count

            lines.append(f"| {count} | {len(results)} | {avg_coverage:.4f} | "
                        f"{avg_runtime:.2f} | {coverage_per_primer:.5f} |")

        lines.append("")

        # Efficiency analysis
        best_efficiency = max(by_count.items(),
                            key=lambda x: np.mean([r['coverage_score']/r['num_primers'] for r in x[1]]))
        lines.append(f"**Most efficient**: {best_efficiency[0]} primers "
                    f"({np.mean([r['coverage_score']/r['num_primers'] for r in best_efficiency[1]]):.5f} coverage/primer)")
        lines.append("")

        lines.append("---")
        lines.append("")

        return lines

    def _analyze_cache_impact(self) -> List[str]:
        """Analyze impact of position cache"""
        lines = []
        lines.append("## Position Cache Impact")
        lines.append("")

        successful = [r for r in self.results if r['success']]

        # Group by cache usage
        with_cache = [r for r in successful if r['use_cache']]
        without_cache = [r for r in successful if not r['use_cache']]

        if with_cache and without_cache:
            avg_runtime_with = np.mean([r['runtime'] for r in with_cache])
            avg_runtime_without = np.mean([r['runtime'] for r in without_cache])

            speedup = avg_runtime_without / avg_runtime_with

            lines.append(f"- With cache: {avg_runtime_with:.2f}s average")
            lines.append(f"- Without cache: {avg_runtime_without:.2f}s average")
            lines.append(f"- **Speedup**: {speedup:.1f}x")
            lines.append("")

            if speedup > 1.5:
                lines.append("✓ **Recommendation**: Use position cache for significant speedup")
            else:
                lines.append("Note: Cache impact minimal for this dataset size")

        lines.append("")
        lines.append("---")
        lines.append("")

        return lines

    def _generate_recommendations(self) -> List[str]:
        """Generate recommendations based on analysis"""
        lines = []
        lines.append("## Recommendations")
        lines.append("")

        successful = [r for r in self.results if r['success']]

        # Best overall method
        by_method = defaultdict(list)
        for r in successful:
            by_method[r['optimization_method']].append(r)

        # Score: balance of speed and quality
        method_scores = {}
        for method, results in by_method.items():
            avg_runtime = np.mean([r['runtime'] for r in results])
            avg_coverage = np.mean([r['coverage_score'] for r in results])

            # Normalize and combine (higher coverage, lower runtime = better)
            score = avg_coverage / (avg_runtime / 10)  # Scale runtime to similar range
            method_scores[method] = (score, avg_runtime, avg_coverage)

        best_method = max(method_scores.items(), key=lambda x: x[1][0])

        lines.append("### Best Overall Method")
        lines.append(f"**{best_method[0]}**")
        lines.append(f"- Average runtime: {best_method[1][1]:.2f}s")
        lines.append(f"- Average coverage: {best_method[1][2]:.4f}")
        lines.append(f"- Balance score: {best_method[1][0]:.2f}")
        lines.append("")

        # For different scenarios
        lines.append("### Scenario-Specific Recommendations")
        lines.append("")

        # Speed priority
        fastest = min(by_method.items(), key=lambda x: np.mean([r['runtime'] for r in x[1]]))
        lines.append(f"**For speed**: Use `{fastest[0]}` "
                    f"({np.mean([r['runtime'] for r in fastest[1]]):.2f}s average)")
        lines.append("")

        # Quality priority
        best_quality = max(by_method.items(), key=lambda x: np.mean([r['coverage_score'] for r in x[1]]))
        lines.append(f"**For best coverage**: Use `{best_quality[0]}` "
                    f"({np.mean([r['coverage_score'] for r in best_quality[1]]):.4f} average)")
        lines.append("")

        # Francisella
        franc_results = [r for r in successful if r['genome'] == 'francisella']
        if franc_results:
            franc_by_method = defaultdict(list)
            for r in franc_results:
                franc_by_method[r['optimization_method']].append(r)
            best_franc = max(franc_by_method.items(),
                           key=lambda x: np.mean([r['coverage_score'] for r in x[1]]))
            lines.append(f"**For Francisella (low GC)**: Use `{best_franc[0]}`")
            lines.append("")

        # Burkholderia
        burkh_results = [r for r in successful if r['genome'] == 'burkholderia']
        if burkh_results:
            burkh_by_method = defaultdict(list)
            for r in burkh_results:
                burkh_by_method[r['optimization_method']].append(r)
            best_burkh = max(burkh_by_method.items(),
                           key=lambda x: np.mean([r['coverage_score'] for r in x[1]]))
            lines.append(f"**For Burkholderia (high GC)**: Use `{best_burkh[0]}`")
            lines.append("")

        lines.append("---")
        lines.append("")

        # Command examples
        lines.append("### Command Examples")
        lines.append("")
        lines.append("```bash")
        lines.append(f"# Best overall")
        lines.append(f"neoswga step4 -j params.json --optimization-method={best_method[0]}")
        lines.append("")
        lines.append(f"# Fastest")
        lines.append(f"neoswga step4 -j params.json --optimization-method={fastest[0]}")
        lines.append("")
        lines.append(f"# Best quality")
        lines.append(f"neoswga step4 -j params.json --optimization-method={best_quality[0]}")
        lines.append("```")
        lines.append("")

        return lines


def main():
    """Main entry point"""
    import argparse

    parser = argparse.ArgumentParser(
        description='Analyze benchmark results and generate report'
    )
    parser.add_argument(
        '--results-file',
        default='./benchmarks/benchmark_results.csv',
        help='Path to benchmark results CSV'
    )
    parser.add_argument(
        '--output-file',
        default='./benchmarks/analysis_report.md',
        help='Output markdown report file'
    )

    args = parser.parse_args()

    # Analyze results
    analyzer = ResultsAnalyzer(results_file=args.results_file)
    report = analyzer.generate_report(output_file=args.output_file)

    # Print summary to console
    logger.info("")
    logger.info("="*80)
    logger.info("ANALYSIS COMPLETE")
    logger.info("="*80)
    logger.info("")
    logger.info(f"Report generated: {args.output_file}")
    logger.info("")
    logger.info("Key findings:")

    # Extract some key lines
    in_recommendations = False
    for line in report:
        if line.startswith("## Recommendations"):
            in_recommendations = True
        elif in_recommendations and line.startswith("**"):
            logger.info(f"  {line}")
        elif in_recommendations and line.startswith("---"):
            break


if __name__ == '__main__':
    main()
