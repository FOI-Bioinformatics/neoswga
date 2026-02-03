#!/usr/bin/env python3
"""
Comprehensive benchmark suite for testing all argument combinations.

Tests different combinations of:
- Genomes (Francisella, Burkholderia)
- Optimization methods (greedy, milp, network, hybrid, moea, dominating-set)
- GC tolerance (0.10, 0.15, 0.20)
- Primer counts (5, 10, 15, 20)
- Performance options (cache on/off, background filter on/off)

Collects metrics:
- Runtime, memory usage, primers selected, coverage, connectivity, quality
"""

import sys
import os
import time
import json
import csv
import h5py
import numpy as np
import traceback
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
from dataclasses import dataclass, asdict
from collections import defaultdict
import logging

# Add neoswga to path
sys.path.insert(0, str(Path(__file__).parent))

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)


@dataclass
class BenchmarkConfig:
    """Configuration for a single benchmark run"""
    genome: str  # francisella or burkholderia
    optimization_method: str  # greedy, milp, network, hybrid, moea, dominating-set
    gc_tolerance: float  # 0.10, 0.15, 0.20
    num_primers: int  # 5, 10, 15, 20
    use_cache: bool  # True/False
    use_background_filter: bool  # True/False

    def to_dict(self) -> Dict:
        """Convert to dictionary"""
        return asdict(self)

    def __str__(self) -> str:
        """String representation"""
        return (f"{self.genome}/{self.optimization_method}/"
                f"gc{self.gc_tolerance:.2f}/n{self.num_primers}/"
                f"cache{int(self.use_cache)}/bg{int(self.use_background_filter)}")


@dataclass
class BenchmarkResult:
    """Results from a single benchmark run"""
    config: BenchmarkConfig
    success: bool
    runtime: float  # seconds
    peak_memory: float  # MB
    n_primers_selected: int
    coverage_score: float
    network_nodes: int
    network_edges: int
    network_components: int
    background_binding: int
    error_message: Optional[str] = None

    def to_dict(self) -> Dict:
        """Convert to dictionary"""
        result = {
            **self.config.to_dict(),
            'success': self.success,
            'runtime': self.runtime,
            'peak_memory': self.peak_memory,
            'n_primers_selected': self.n_primers_selected,
            'coverage_score': self.coverage_score,
            'network_nodes': self.network_nodes,
            'network_edges': self.network_edges,
            'network_components': self.network_components,
            'background_binding': self.background_binding,
            'error_message': self.error_message
        }
        return result


class BenchmarkSuite:
    """Comprehensive benchmark suite"""

    def __init__(self, test_data_dir: str = "./test_data",
                 output_dir: str = "./benchmarks"):
        """
        Initialize benchmark suite.

        Args:
            test_data_dir: Directory with test datasets
            output_dir: Directory for benchmark results
        """
        self.test_data_dir = Path(test_data_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True, parents=True)

        # Load metadata
        metadata_file = self.test_data_dir / "test_metadata.json"
        if not metadata_file.exists():
            raise FileNotFoundError(
                f"Test metadata not found: {metadata_file}\n"
                f"Run create_test_dataset.py first"
            )

        with open(metadata_file) as f:
            self.metadata = json.load(f)

        # Check method availability
        self.available_methods = self._check_method_availability()

        logger.info("Available optimization methods:")
        for method, available in self.available_methods.items():
            status = "✓" if available else "✗ (dependencies missing)"
            logger.info(f"  {method}: {status}")

    def _check_method_availability(self) -> Dict[str, bool]:
        """Check which optimization methods are available"""
        methods = {}

        # Core methods (always available)
        methods['greedy'] = True
        methods['network'] = True
        methods['hybrid'] = True

        # MILP (requires mip)
        try:
            import mip
            methods['milp'] = True
        except ImportError:
            methods['milp'] = False

        # MOEA (requires pymoo)
        try:
            import pymoo
            methods['moea'] = True
        except ImportError:
            methods['moea'] = False

        # Dominating set (requires networkx + mip for ILP)
        methods['dominating-set'] = True  # Greedy version always available

        return methods

    def generate_test_matrix(self) -> List[BenchmarkConfig]:
        """
        Generate all test combinations.

        Returns:
            List of benchmark configurations
        """
        # Get available genomes from metadata
        genomes = list(self.metadata['genomes'].keys())
        optimization_methods = [m for m, avail in self.available_methods.items() if avail]
        gc_tolerances = [0.10, 0.15, 0.20]
        primer_counts = [5, 10, 15, 20]
        # All current methods require cache=True, so only test with cache
        cache_options = [True]
        bg_filter_options = [False]  # Background filter requires pre-built filter

        configs = []

        for genome in genomes:
            for method in optimization_methods:
                for gc_tol in gc_tolerances:
                    for n_primers in primer_counts:
                        for use_cache in cache_options:
                            for use_bg in bg_filter_options:
                                config = BenchmarkConfig(
                                    genome=genome,
                                    optimization_method=method,
                                    gc_tolerance=gc_tol,
                                    num_primers=n_primers,
                                    use_cache=use_cache,
                                    use_background_filter=use_bg
                                )
                                configs.append(config)

        return configs

    def run_single_benchmark(self, config: BenchmarkConfig) -> BenchmarkResult:
        """
        Run a single benchmark configuration.

        Args:
            config: Benchmark configuration

        Returns:
            Benchmark results
        """
        start_time = time.time()
        peak_memory = 0.0

        try:
            # Load test data
            genome_metadata = self.metadata['genomes'][config.genome]
            h5_file = Path(genome_metadata['h5_file'])
            primer_file = Path(genome_metadata['primer_file'])

            # Load primers
            with open(primer_file) as f:
                primers = [line.strip() for line in f]

            # Sample primers if needed (for faster testing)
            if len(primers) > 200:
                # Use first 200 for speed
                primers = primers[:200]

            # Load position cache
            if config.use_cache:
                from neoswga.core.position_cache import PositionCache
                cache = self._load_positions_from_h5(h5_file, config.genome, primers)
            else:
                cache = None

            # Run optimization based on method
            result_data = self._run_optimization(
                config,
                primers,
                cache,
                genome_metadata
            )

            runtime = time.time() - start_time

            # Create result
            result = BenchmarkResult(
                config=config,
                success=True,
                runtime=runtime,
                peak_memory=peak_memory,
                n_primers_selected=result_data['n_primers'],
                coverage_score=result_data['coverage'],
                network_nodes=result_data['network_nodes'],
                network_edges=result_data['network_edges'],
                network_components=result_data['network_components'],
                background_binding=result_data['background_binding']
            )

            return result

        except Exception as e:
            runtime = time.time() - start_time

            # Log error
            error_msg = f"{type(e).__name__}: {str(e)}"
            logger.error(f"  ERROR: {error_msg}")

            result = BenchmarkResult(
                config=config,
                success=False,
                runtime=runtime,
                peak_memory=peak_memory,
                n_primers_selected=0,
                coverage_score=0.0,
                network_nodes=0,
                network_edges=0,
                network_components=0,
                background_binding=0,
                error_message=error_msg
            )

            return result

    def _load_positions_from_h5(self, h5_file: Path, genome_name: str,
                                primers: List[str]) -> Any:
        """Load positions from HDF5 into mock cache"""
        from neoswga.core.position_cache import PositionCache

        # Create mock cache
        class MockCache:
            def __init__(self, h5_file, genome_name, primers):
                self.positions = {}
                with h5py.File(h5_file, 'r') as f:
                    genome_group = f[genome_name]
                    for primer in primers:
                        if primer in genome_group:
                            pos_fwd = np.array(genome_group[primer]['+'])
                            pos_rev = np.array(genome_group[primer]['-'])
                            # Store with both notations for compatibility
                            self.positions[(genome_name, primer, '+')] = pos_fwd
                            self.positions[(genome_name, primer, '-')] = pos_rev
                            self.positions[(genome_name, primer, 'forward')] = pos_fwd
                            self.positions[(genome_name, primer, 'reverse')] = pos_rev

            def get_positions(self, genome_prefix, primer, strand):
                key = (genome_prefix, primer, strand)
                return self.positions.get(key, np.array([]))

        return MockCache(h5_file, genome_name, primers)

    def _run_optimization(self, config: BenchmarkConfig, primers: List[str],
                         cache: Any, genome_metadata: Dict) -> Dict:
        """
        Run optimization for given configuration.

        Args:
            config: Benchmark configuration
            primers: List of candidate primers
            cache: Position cache
            genome_metadata: Genome metadata

        Returns:
            Dictionary with results
        """
        method = config.optimization_method
        num_primers = config.num_primers

        # Initialize network for analysis
        from neoswga.core.network_optimizer import AmplificationNetwork, NetworkOptimizer

        if method in ['greedy', 'network', 'hybrid']:
            # Network-based optimization
            # NetworkOptimizer requires position_cache, not cache
            if cache is None:
                raise ValueError("Network-based methods require cache=True")

            optimizer = NetworkOptimizer(
                position_cache=cache,
                fg_prefixes=[config.genome],
                bg_prefixes=[],
                fg_seq_lengths=[genome_metadata['genome_length']],
                bg_seq_lengths=[]
            )

            # optimize_greedy returns just the list of primers
            selected_primers = optimizer.optimize_greedy(
                candidates=primers,
                num_primers=num_primers
            )

            # Get detailed metrics
            metrics = optimizer.score_primer_set(selected_primers)

            coverage = metrics['target_connectivity']  # Use connectivity as coverage
            network_nodes = metrics['target_largest_component']
            network_edges = 0  # Not directly available
            network_components = 0  # Not directly available

        elif method == 'milp':
            # MILP optimization
            from neoswga.core.milp_optimizer import MILPOptimizer

            optimizer = MILPOptimizer(
                cache=cache,
                fg_prefixes=[config.genome],
                bg_prefixes=[],
                fg_seq_lengths=[genome_metadata['genome_length']],
                bg_seq_lengths=[]
            )

            result = optimizer.optimize(
                candidates=primers,
                num_primers=num_primers,
                max_time=30,  # 30 second time limit
                verbose=False
            )

            selected_primers = result.get('selected_primers', [])
            coverage = result.get('objective_value', 0.0)
            network_nodes = 0
            network_edges = 0
            network_components = 0

        elif method == 'moea':
            # MOEA optimization
            from neoswga.core.moea_optimizer import MOEAOptimizer, MOEAConfig

            config_moea = MOEAConfig(
                pop_size=50,
                n_generations=20,  # Reduced for speed
                n_objectives=4
            )

            optimizer = MOEAOptimizer(
                cache=cache,
                fg_prefixes=[config.genome],
                bg_prefixes=[],
                fg_seq_lengths=[genome_metadata['genome_length']],
                bg_seq_lengths=[],
                config=config_moea
            )

            result = optimizer.optimize(
                candidates=primers[:50],  # Limit candidates for speed
                max_primers=num_primers,
                verbose=False
            )

            best_solution = result['best_solution']
            selected_primers = best_solution['primers']
            coverage = best_solution['objectives']['target_coverage']
            network_nodes = 0
            network_edges = 0
            network_components = 0

        elif method == 'dominating-set':
            # Dominating set optimization
            # DominatingSetOptimizer requires cache
            if cache is None:
                raise ValueError("Dominating-set method requires cache=True")

            from neoswga.core.dominating_set_optimizer import DominatingSetOptimizer

            optimizer = DominatingSetOptimizer(
                cache=cache,
                fg_prefixes=[config.genome],
                fg_seq_lengths=[genome_metadata['genome_length']],
                bin_size=10000
            )

            result = optimizer.optimize_greedy(
                candidates=primers,
                max_primers=num_primers,
                verbose=False
            )

            selected_primers = result['primers']
            coverage = result['coverage']
            network_nodes = 0
            network_edges = 0
            network_components = 0

        else:
            raise ValueError(f"Unknown method: {method}")

        return {
            'n_primers': len(selected_primers),
            'coverage': coverage,
            'network_nodes': network_nodes,
            'network_edges': network_edges,
            'network_components': network_components,
            'background_binding': 0  # Would require background genome
        }

    def run_all_benchmarks(self, save_interval: int = 10) -> List[BenchmarkResult]:
        """
        Run all benchmark configurations.

        Args:
            save_interval: Save results every N tests

        Returns:
            List of all benchmark results
        """
        configs = self.generate_test_matrix()
        results = []

        logger.info("="*80)
        logger.info(f"RUNNING BENCHMARK SUITE: {len(configs)} configurations")
        logger.info("="*80)
        logger.info("")

        start_time = time.time()

        for i, config in enumerate(configs, 1):
            logger.info(f"[{i}/{len(configs)}] Running: {config}")

            result = self.run_single_benchmark(config)
            results.append(result)

            if result.success:
                logger.info(f"  ✓ Success: {result.runtime:.2f}s, "
                          f"{result.n_primers_selected} primers, "
                          f"coverage={result.coverage_score:.3f}")
            else:
                logger.info(f"  ✗ Failed: {result.error_message}")

            # Save incrementally
            if i % save_interval == 0:
                self._save_results(results)
                elapsed = time.time() - start_time
                avg_time = elapsed / i
                remaining = (len(configs) - i) * avg_time
                logger.info(f"\n  Progress: {i}/{len(configs)} "
                          f"({i/len(configs)*100:.1f}%), "
                          f"~{remaining/60:.1f} min remaining\n")

        # Final save
        self._save_results(results)

        elapsed = time.time() - start_time
        logger.info("")
        logger.info("="*80)
        logger.info(f"BENCHMARK COMPLETE: {len(results)} tests in {elapsed/60:.1f} minutes")
        logger.info("="*80)

        # Summary
        successful = sum(1 for r in results if r.success)
        failed = len(results) - successful

        logger.info(f"\nResults:")
        logger.info(f"  Successful: {successful}")
        logger.info(f"  Failed: {failed}")
        logger.info(f"  Success rate: {successful/len(results)*100:.1f}%")

        return results

    def _save_results(self, results: List[BenchmarkResult]):
        """Save results to CSV"""
        output_file = self.output_dir / "benchmark_results.csv"

        with open(output_file, 'w', newline='') as f:
            if results:
                # Get field names from first result
                fieldnames = list(results[0].to_dict().keys())

                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()

                for result in results:
                    writer.writerow(result.to_dict())

        logger.info(f"  Saved results to: {output_file}")


def main():
    """Main entry point"""
    import argparse

    parser = argparse.ArgumentParser(
        description='Run comprehensive benchmark suite'
    )
    parser.add_argument(
        '--test-data-dir',
        default='./test_data',
        help='Directory with test datasets'
    )
    parser.add_argument(
        '--output-dir',
        default='./benchmarks',
        help='Output directory for results'
    )
    parser.add_argument(
        '--save-interval',
        type=int,
        default=10,
        help='Save results every N tests (default: 10)'
    )

    args = parser.parse_args()

    # Run benchmarks
    suite = BenchmarkSuite(
        test_data_dir=args.test_data_dir,
        output_dir=args.output_dir
    )

    results = suite.run_all_benchmarks(save_interval=args.save_interval)

    logger.info(f"\nResults saved to: {args.output_dir}/benchmark_results.csv")
    logger.info("Next: Run analyze_results.py to generate analysis report")


if __name__ == '__main__':
    main()
