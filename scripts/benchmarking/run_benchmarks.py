#!/usr/bin/env python3
"""
Automated test runner with progress tracking and resource monitoring.

Executes the complete benchmark suite with:
- Progress bars
- Time estimation
- Checkpointing (resume from failures)
- Resource monitoring
- Summary statistics
"""

import sys
import os
import time
import json
from pathlib import Path
from datetime import datetime, timedelta
import logging

# Add neoswga to path
sys.path.insert(0, str(Path(__file__).parent))

from benchmark_suite import BenchmarkSuite

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)


class BenchmarkRunner:
    """Automated benchmark runner with progress tracking"""

    def __init__(self, test_data_dir: str = "./test_data",
                 output_dir: str = "./benchmarks",
                 checkpoint_file: str = None):
        """
        Initialize runner.

        Args:
            test_data_dir: Directory with test datasets
            output_dir: Output directory for results
            checkpoint_file: File for saving checkpoint state
        """
        self.suite = BenchmarkSuite(
            test_data_dir=test_data_dir,
            output_dir=output_dir
        )

        self.output_dir = Path(output_dir)
        self.checkpoint_file = checkpoint_file or (self.output_dir / "checkpoint.json")

        self.start_time = None
        self.completed_configs = set()

    def load_checkpoint(self):
        """Load checkpoint to resume from previous run"""
        if self.checkpoint_file.exists():
            logger.info(f"Loading checkpoint from: {self.checkpoint_file}")
            with open(self.checkpoint_file) as f:
                checkpoint = json.load(f)

            self.completed_configs = set(
                tuple(sorted(c.items())) for c in checkpoint['completed_configs']
            )

            logger.info(f"Resuming from checkpoint: {len(self.completed_configs)} tests already complete")
        else:
            logger.info("No checkpoint found, starting from beginning")

    def save_checkpoint(self, completed_config):
        """Save checkpoint after completing a test"""
        self.completed_configs.add(tuple(sorted(completed_config.items())))

        checkpoint = {
            'completed_configs': [dict(c) for c in self.completed_configs],
            'last_update': datetime.now().isoformat()
        }

        with open(self.checkpoint_file, 'w') as f:
            json.dump(checkpoint, f, indent=2)

    def is_completed(self, config):
        """Check if configuration has already been tested"""
        config_tuple = tuple(sorted(config.to_dict().items()))
        return config_tuple in self.completed_configs

    def print_header(self):
        """Print banner"""
        logger.info("")
        logger.info("#" * 80)
        logger.info("# NEOSWGA COMPREHENSIVE BENCHMARK SUITE")
        logger.info("#" * 80)
        logger.info("")
        # Get genome info from metadata
        genome_info = []
        for genome_name, metadata in self.suite.metadata['genomes'].items():
            gc = metadata['gc_content']
            genome_info.append(f"{genome_name.capitalize()} ({gc:.1%} GC)")

        logger.info("Testing all combinations of:")
        logger.info(f"  - Genomes: {', '.join(genome_info)}")
        logger.info("  - Methods: " + ", ".join(
            m for m, a in self.suite.available_methods.items() if a
        ))
        logger.info("  - GC tolerances: 0.10, 0.15, 0.20")
        logger.info("  - Primer counts: 5, 10, 15, 20")
        logger.info("  - Cache options: on, off")
        logger.info("")

    def print_progress(self, current: int, total: int, elapsed: float):
        """
        Print progress bar and time estimates.

        Args:
            current: Current test number
            total: Total number of tests
            elapsed: Elapsed time in seconds
        """
        percent = current / total * 100
        avg_time = elapsed / current if current > 0 else 0
        remaining = (total - current) * avg_time

        # Progress bar
        bar_length = 50
        filled = int(bar_length * current / total)
        bar = "█" * filled + "░" * (bar_length - filled)

        # Time formatting
        elapsed_str = str(timedelta(seconds=int(elapsed)))
        remaining_str = str(timedelta(seconds=int(remaining)))

        logger.info(f"\n[{bar}] {percent:.1f}%")
        logger.info(f"Progress: {current}/{total} tests")
        logger.info(f"Elapsed: {elapsed_str}, Remaining: ~{remaining_str}")
        logger.info(f"Avg time per test: {avg_time:.1f}s")
        logger.info("")

    def run_with_progress(self, save_interval: int = 10,
                         progress_interval: int = 10) -> list:
        """
        Run benchmarks with progress tracking.

        Args:
            save_interval: Save results every N tests
            progress_interval: Print progress every N tests

        Returns:
            List of benchmark results
        """
        self.print_header()

        # Load checkpoint if exists
        self.load_checkpoint()

        # Generate test matrix
        configs = self.suite.generate_test_matrix()

        # Filter out completed tests
        remaining_configs = [c for c in configs if not self.is_completed(c)]

        logger.info(f"Total tests: {len(configs)}")
        logger.info(f"Already completed: {len(self.completed_configs)}")
        logger.info(f"Remaining: {len(remaining_configs)}")
        logger.info("")

        if not remaining_configs:
            logger.info("All tests already completed!")
            logger.info("Delete checkpoint.json to re-run")
            return []

        # Confirm start
        logger.info(f"Starting benchmark suite...")
        logger.info(f"Estimated time: {len(remaining_configs) * 5 / 60:.1f} minutes")
        logger.info("(Ctrl+C to pause, will resume from checkpoint)")
        logger.info("")
        logger.info("="*80)
        logger.info("")

        results = []
        self.start_time = time.time()

        try:
            for i, config in enumerate(remaining_configs, 1):
                # Run test
                logger.info(f"[{i}/{len(remaining_configs)}] Running: {config}")
                result = self.suite.run_single_benchmark(config)
                results.append(result)

                # Log result
                if result.success:
                    logger.info(f"  ✓ Success: {result.runtime:.2f}s, "
                              f"{result.n_primers_selected} primers, "
                              f"coverage={result.coverage_score:.3f}")
                else:
                    logger.info(f"  ✗ Failed: {result.error_message}")

                # Save checkpoint
                self.save_checkpoint(config.to_dict())

                # Save results periodically
                if i % save_interval == 0:
                    self.suite._save_results(results)

                # Print progress
                if i % progress_interval == 0:
                    elapsed = time.time() - self.start_time
                    self.print_progress(i, len(remaining_configs), elapsed)

        except KeyboardInterrupt:
            logger.info("\n\nBenchmark interrupted by user")
            logger.info("Progress saved to checkpoint")
            logger.info("Run again to resume from where you left off")

            # Save partial results
            if results:
                self.suite._save_results(results)

            return results

        # Final save
        self.suite._save_results(results)

        # Print final summary
        elapsed = time.time() - self.start_time
        self.print_summary(results, elapsed)

        return results

    def print_summary(self, results: list, elapsed: float):
        """
        Print final summary.

        Args:
            results: List of benchmark results
            elapsed: Total elapsed time
        """
        logger.info("")
        logger.info("="*80)
        logger.info("BENCHMARK COMPLETE")
        logger.info("="*80)
        logger.info("")

        successful = sum(1 for r in results if r.success)
        failed = len(results) - successful

        logger.info(f"Tests run: {len(results)}")
        logger.info(f"Successful: {successful} ({successful/len(results)*100:.1f}%)")
        logger.info(f"Failed: {failed}")
        logger.info(f"Total time: {timedelta(seconds=int(elapsed))}")
        logger.info("")

        if successful > 0:
            # Statistics
            runtimes = [r.runtime for r in results if r.success]
            avg_runtime = sum(runtimes) / len(runtimes)
            min_runtime = min(runtimes)
            max_runtime = max(runtimes)

            logger.info("Runtime statistics:")
            logger.info(f"  Average: {avg_runtime:.2f}s")
            logger.info(f"  Min: {min_runtime:.2f}s")
            logger.info(f"  Max: {max_runtime:.2f}s")
            logger.info("")

            # Best performers
            logger.info("Top 5 fastest methods:")
            sorted_results = sorted(
                [r for r in results if r.success],
                key=lambda r: r.runtime
            )
            for i, result in enumerate(sorted_results[:5], 1):
                logger.info(f"  {i}. {result.config.optimization_method} "
                          f"({result.config.genome}, n={result.config.num_primers}): "
                          f"{result.runtime:.2f}s")
            logger.info("")

            # Best quality
            logger.info("Top 5 best coverage:")
            sorted_by_coverage = sorted(
                [r for r in results if r.success],
                key=lambda r: r.coverage_score,
                reverse=True
            )
            for i, result in enumerate(sorted_by_coverage[:5], 1):
                logger.info(f"  {i}. {result.config.optimization_method} "
                          f"({result.config.genome}, n={result.config.num_primers}): "
                          f"{result.coverage_score:.4f}")
            logger.info("")

        logger.info(f"Results saved to: {self.output_dir}/benchmark_results.csv")
        logger.info("")
        logger.info("Next steps:")
        logger.info("  1. Run analyze_results.py to generate detailed analysis")
        logger.info("  2. Review benchmark_results.csv for raw data")
        logger.info("")
        logger.info("#"*80)


def main():
    """Main entry point"""
    import argparse

    parser = argparse.ArgumentParser(
        description='Run comprehensive benchmark suite with progress tracking'
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
    parser.add_argument(
        '--progress-interval',
        type=int,
        default=10,
        help='Print progress every N tests (default: 10)'
    )
    parser.add_argument(
        '--reset',
        action='store_true',
        help='Reset checkpoint and start from beginning'
    )

    args = parser.parse_args()

    # Create runner
    runner = BenchmarkRunner(
        test_data_dir=args.test_data_dir,
        output_dir=args.output_dir
    )

    # Reset checkpoint if requested
    if args.reset and runner.checkpoint_file.exists():
        logger.info("Resetting checkpoint...")
        runner.checkpoint_file.unlink()

    # Run benchmarks
    results = runner.run_with_progress(
        save_interval=args.save_interval,
        progress_interval=args.progress_interval
    )


if __name__ == '__main__':
    main()
