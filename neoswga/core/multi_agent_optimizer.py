"""
Multi-Agent Optimization Orchestrator for SWGA primer selection.

Coordinates multiple optimizer agents to find optimal primer sets through:
- Parallel execution of multiple optimization strategies
- Performance profiling and bottleneck identification
- Intelligent result aggregation and ensemble methods
- Cost tracking and resource management

Usage:
    orchestrator = MultiAgentOrchestrator(position_cache, fg_prefixes, fg_seq_lengths)
    result = orchestrator.optimize_parallel(
        candidates,
        methods=['greedy', 'dominating-set', 'hybrid'],
        target_size=10
    )
"""

import time
import logging
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed, TimeoutError as FuturesTimeoutError
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Any, Tuple
from enum import Enum

from .base_optimizer import OptimizationResult, OptimizationStatus, PrimerSetMetrics, OptimizerConfig
from .optimizer_factory import OptimizerFactory, OptimizerRegistry
from .exceptions import NeoSWGAError, OperationTimeoutError

logger = logging.getLogger(__name__)


class AggregationStrategy(Enum):
    """Strategy for combining results from multiple optimizers."""
    BEST_SCORE = "best_score"  # Pick result with highest score
    CONSENSUS = "consensus"  # Pick primers that appear in multiple results
    WEIGHTED_VOTE = "weighted_vote"  # Weight by optimizer performance
    ENSEMBLE = "ensemble"  # Combine metrics from all results


@dataclass
class AgentPerformanceMetrics:
    """Performance metrics for a single optimizer agent."""
    optimizer_name: str
    execution_time: float
    peak_memory_mb: float = 0.0
    iterations_used: int = 0
    success: bool = False
    score: float = float('-inf')
    error_message: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return {
            'optimizer': self.optimizer_name,
            'execution_time_s': round(self.execution_time, 3),
            'peak_memory_mb': round(self.peak_memory_mb, 2),
            'iterations': self.iterations_used,
            'success': self.success,
            'score': round(self.score, 4) if self.score != float('-inf') else None,
            'error': self.error_message or None,
        }


@dataclass
class OrchestratorResult:
    """Result from multi-agent optimization."""
    best_result: OptimizationResult
    all_results: Dict[str, OptimizationResult]
    performance_metrics: Dict[str, AgentPerformanceMetrics]
    aggregation_strategy: AggregationStrategy
    total_time: float
    agents_used: int
    agents_succeeded: int

    @property
    def success_rate(self) -> float:
        """Fraction of agents that succeeded."""
        return self.agents_succeeded / self.agents_used if self.agents_used > 0 else 0.0

    def summary(self) -> str:
        """Generate a summary report."""
        lines = [
            "=" * 60,
            "Multi-Agent Optimization Summary",
            "=" * 60,
            f"Strategy: {self.aggregation_strategy.value}",
            f"Agents used: {self.agents_used}",
            f"Agents succeeded: {self.agents_succeeded} ({self.success_rate:.0%})",
            f"Total time: {self.total_time:.2f}s",
            "",
            "Best Result:",
            f"  Optimizer: {self.best_result.optimizer_name}",
            f"  Score: {self.best_result.score:.4f}",
            f"  Primers: {self.best_result.num_primers}",
            f"  Coverage: {self.best_result.metrics.fg_coverage:.1%}",
            "",
            "Agent Performance:",
        ]

        for name, metrics in sorted(
            self.performance_metrics.items(),
            key=lambda x: x[1].execution_time
        ):
            status = "OK" if metrics.success else "FAIL"
            score_str = f"{metrics.score:.4f}" if metrics.success else "N/A"
            lines.append(
                f"  {name}: {metrics.execution_time:.2f}s, score={score_str} [{status}]"
            )

        lines.append("=" * 60)
        return "\n".join(lines)


class MultiAgentOrchestrator:
    """
    Orchestrator for coordinating multiple optimizer agents.

    Features:
    - Parallel execution: Run multiple optimizers concurrently
    - Performance profiling: Track execution time and resource usage
    - Result aggregation: Combine results using various strategies
    - Fault tolerance: Handle individual optimizer failures gracefully
    - Thread pool reuse: Maintains thread pool across calls for efficiency
    """

    # Class-level default max workers (can be overridden per instance)
    DEFAULT_MAX_WORKERS = 8

    def __init__(
        self,
        position_cache,
        fg_prefixes: List[str],
        fg_seq_lengths: List[int],
        bg_prefixes: Optional[List[str]] = None,
        bg_seq_lengths: Optional[List[int]] = None,
        max_workers: Optional[int] = None,
    ):
        """
        Initialize the orchestrator.

        Args:
            position_cache: PositionCache with primer binding data
            fg_prefixes: Foreground genome prefixes
            fg_seq_lengths: Foreground genome lengths
            bg_prefixes: Background genome prefixes (optional)
            bg_seq_lengths: Background genome lengths (optional)
            max_workers: Maximum parallel workers (defaults to DEFAULT_MAX_WORKERS)
        """
        self.position_cache = position_cache
        self.fg_prefixes = fg_prefixes
        self.fg_seq_lengths = fg_seq_lengths
        self.bg_prefixes = bg_prefixes or []
        self.bg_seq_lengths = bg_seq_lengths or []
        self.max_workers = max_workers or self.DEFAULT_MAX_WORKERS

        # Performance tracking
        self._performance_history: List[Dict[str, AgentPerformanceMetrics]] = []
        self._lock = threading.Lock()

        # Reusable thread pool (lazy initialization)
        self._executor: Optional[ThreadPoolExecutor] = None
        self._executor_lock = threading.Lock()

    def _get_executor(self, workers: int) -> ThreadPoolExecutor:
        """
        Get or create a thread pool executor.

        Reuses existing executor if worker count matches, otherwise creates new one.
        Thread-safe lazy initialization.

        Args:
            workers: Number of workers needed (minimum 1)

        Returns:
            ThreadPoolExecutor instance
        """
        # Ensure at least 1 worker
        workers = max(1, workers)

        with self._executor_lock:
            # Create new executor if none exists or if worker count changed
            if self._executor is None or self._executor._max_workers < workers:
                # Shutdown old executor if exists
                if self._executor is not None:
                    self._executor.shutdown(wait=False)

                self._executor = ThreadPoolExecutor(max_workers=workers)
                logger.debug(f"Created new thread pool with {workers} workers")

            return self._executor

    def shutdown(self, wait: bool = True) -> None:
        """
        Shutdown the thread pool executor.

        Call this when done with the orchestrator to free resources.

        Args:
            wait: If True, wait for pending tasks to complete
        """
        with self._executor_lock:
            if self._executor is not None:
                self._executor.shutdown(wait=wait)
                self._executor = None
                logger.debug("Thread pool executor shutdown")

    def __del__(self):
        """Cleanup executor on garbage collection."""
        try:
            self.shutdown(wait=False)
        except Exception:
            pass  # Ignore errors during cleanup

    def __enter__(self):
        """Context manager entry."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit - shutdown executor."""
        self.shutdown(wait=True)
        return False

    def list_available_agents(self) -> Dict[str, str]:
        """List all available optimizer agents."""
        from .unified_optimizer import _ensure_optimizers_registered
        _ensure_optimizers_registered()
        return OptimizerRegistry.list_all()

    def optimize_parallel(
        self,
        candidates: List[str],
        methods: Optional[List[str]] = None,
        target_size: int = 6,
        aggregation: AggregationStrategy = AggregationStrategy.BEST_SCORE,
        timeout_per_agent: float = 300.0,
        config: Optional[OptimizerConfig] = None,
        verbose: bool = True,
    ) -> OrchestratorResult:
        """
        Run multiple optimizers in parallel and aggregate results.

        Args:
            candidates: List of candidate primer sequences
            methods: Optimizer methods to use (defaults to all available)
            target_size: Target number of primers
            aggregation: Strategy for combining results
            timeout_per_agent: Maximum time per optimizer (seconds)
            config: Base configuration for optimizers
            verbose: Print progress information

        Returns:
            OrchestratorResult with aggregated results and performance metrics
        """
        from .unified_optimizer import _ensure_optimizers_registered
        _ensure_optimizers_registered()

        # Default to commonly used optimizers
        if methods is None:
            methods = ['greedy', 'dominating-set', 'hybrid']

        # Filter to available methods
        available = set(OptimizerRegistry.list_all().keys())
        methods = [m for m in methods if m in available]

        if not methods:
            raise ValueError("No valid optimizer methods specified")

        if verbose:
            logger.info(f"Multi-agent optimization: {len(methods)} agents")
            logger.info(f"  Methods: {', '.join(methods)}")
            logger.info(f"  Candidates: {len(candidates)}")
            logger.info(f"  Target size: {target_size}")

        start_time = time.time()

        # Prepare optimizer config
        base_config = config or OptimizerConfig(
            target_set_size=target_size,
            verbose=False,  # Suppress individual optimizer output
        )

        # Run optimizers in parallel using reusable executor
        results: Dict[str, OptimizationResult] = {}
        metrics: Dict[str, AgentPerformanceMetrics] = {}

        n_workers = min(self.max_workers, len(methods))
        executor = self._get_executor(n_workers)

        futures = {}

        for method in methods:
            future = executor.submit(
                self._run_agent,
                method,
                candidates,
                target_size,
                base_config,
                timeout_per_agent,
            )
            futures[future] = method

        # Use as_completed with total timeout, then individual timeouts
        total_timeout = timeout_per_agent * len(methods) * 1.5  # Allow some margin
        try:
            for future in as_completed(futures, timeout=total_timeout):
                method = futures[future]
                try:
                    # Enforce per-agent timeout
                    result, perf = future.result(timeout=timeout_per_agent)
                    results[method] = result
                    metrics[method] = perf

                    if verbose and result.is_success:
                        logger.info(
                            f"  {method}: score={result.score:.4f}, "
                            f"time={perf.execution_time:.2f}s"
                        )

                except FuturesTimeoutError:
                    logger.warning(f"Agent {method} timed out after {timeout_per_agent}s")
                    metrics[method] = AgentPerformanceMetrics(
                        optimizer_name=method,
                        execution_time=timeout_per_agent,
                        success=False,
                        error_message=f"Timeout after {timeout_per_agent}s",
                    )

                except Exception as e:
                    logger.error(f"Agent {method} failed: {e}")
                    metrics[method] = AgentPerformanceMetrics(
                        optimizer_name=method,
                        execution_time=0.0,
                        success=False,
                        error_message=str(e),
                    )

        except FuturesTimeoutError:
            # Total timeout exceeded - mark remaining futures as timed out
            logger.warning(f"Total timeout exceeded ({total_timeout}s)")
            for future, method in futures.items():
                if method not in metrics:
                    metrics[method] = AgentPerformanceMetrics(
                        optimizer_name=method,
                        execution_time=total_timeout,
                        success=False,
                        error_message=f"Total timeout exceeded",
                    )

        total_time = time.time() - start_time

        # Aggregate results
        best_result = self._aggregate_results(results, aggregation, target_size)

        agents_succeeded = sum(1 for m in metrics.values() if m.success)

        orchestrator_result = OrchestratorResult(
            best_result=best_result,
            all_results=results,
            performance_metrics=metrics,
            aggregation_strategy=aggregation,
            total_time=total_time,
            agents_used=len(methods),
            agents_succeeded=agents_succeeded,
        )

        # Store in history
        with self._lock:
            self._performance_history.append(metrics)

        if verbose:
            logger.info(f"\nCompleted in {total_time:.2f}s")
            logger.info(f"Best: {best_result.optimizer_name} (score={best_result.score:.4f})")

        return orchestrator_result

    def _run_agent(
        self,
        method: str,
        candidates: List[str],
        target_size: int,
        config: OptimizerConfig,
        timeout: float,
    ) -> Tuple[OptimizationResult, AgentPerformanceMetrics]:
        """Run a single optimizer agent with performance tracking."""
        start_time = time.time()

        try:
            optimizer = OptimizerFactory.create(
                name=method,
                position_cache=self.position_cache,
                fg_prefixes=self.fg_prefixes,
                fg_seq_lengths=self.fg_seq_lengths,
                bg_prefixes=self.bg_prefixes,
                bg_seq_lengths=self.bg_seq_lengths,
                config=config,
            )

            result = optimizer.optimize(candidates, target_size)

            execution_time = time.time() - start_time

            metrics = AgentPerformanceMetrics(
                optimizer_name=method,
                execution_time=execution_time,
                iterations_used=result.iterations,
                success=result.is_success,
                score=result.score if result.is_success else float('-inf'),
            )

            return result, metrics

        except Exception as e:
            execution_time = time.time() - start_time
            logger.error(f"Optimizer {method} failed: {e}")

            return (
                OptimizationResult.failure(method, str(e)),
                AgentPerformanceMetrics(
                    optimizer_name=method,
                    execution_time=execution_time,
                    success=False,
                    error_message=str(e),
                ),
            )

    def _aggregate_results(
        self,
        results: Dict[str, OptimizationResult],
        strategy: AggregationStrategy,
        target_size: int,
    ) -> OptimizationResult:
        """Aggregate results from multiple optimizers."""
        successful = {k: v for k, v in results.items() if v.is_success}

        if not successful:
            # All failed - return first failure
            first_result = next(iter(results.values()), None)
            if first_result:
                return first_result
            return OptimizationResult.failure("orchestrator", "All agents failed")

        if strategy == AggregationStrategy.BEST_SCORE:
            # Pick result with highest score
            best_name = max(successful.keys(), key=lambda k: successful[k].score)
            return successful[best_name]

        elif strategy == AggregationStrategy.CONSENSUS:
            # Pick primers that appear in multiple results
            primer_counts: Dict[str, int] = {}
            for result in successful.values():
                for primer in result.primers:
                    primer_counts[primer] = primer_counts.get(primer, 0) + 1

            # Sort by frequency, then select top target_size
            sorted_primers = sorted(
                primer_counts.keys(),
                key=lambda p: primer_counts[p],
                reverse=True
            )
            consensus_primers = tuple(sorted_primers[:target_size])

            # Use best result as template
            best = max(successful.values(), key=lambda r: r.score)

            return OptimizationResult(
                primers=consensus_primers,
                score=best.score,  # Would ideally recalculate
                status=OptimizationStatus.SUCCESS,
                metrics=best.metrics,
                iterations=sum(r.iterations for r in successful.values()),
                optimizer_name="consensus",
                message=f"Consensus from {len(successful)} optimizers",
            )

        elif strategy == AggregationStrategy.WEIGHTED_VOTE:
            # Weight primers by optimizer score
            primer_weights: Dict[str, float] = {}
            total_score = sum(r.score for r in successful.values() if r.score > 0)

            for result in successful.values():
                weight = result.score / total_score if total_score > 0 else 1.0
                for primer in result.primers:
                    primer_weights[primer] = primer_weights.get(primer, 0) + weight

            sorted_primers = sorted(
                primer_weights.keys(),
                key=lambda p: primer_weights[p],
                reverse=True
            )
            weighted_primers = tuple(sorted_primers[:target_size])

            best = max(successful.values(), key=lambda r: r.score)

            return OptimizationResult(
                primers=weighted_primers,
                score=best.score,
                status=OptimizationStatus.SUCCESS,
                metrics=best.metrics,
                iterations=sum(r.iterations for r in successful.values()),
                optimizer_name="weighted-vote",
                message=f"Weighted vote from {len(successful)} optimizers",
            )

        else:
            # Default: ENSEMBLE (same as best for now)
            best_name = max(successful.keys(), key=lambda k: successful[k].score)
            return successful[best_name]

    def get_performance_history(self) -> List[Dict[str, Any]]:
        """Get historical performance data for analysis."""
        with self._lock:
            return [
                {name: m.to_dict() for name, m in run.items()}
                for run in self._performance_history
            ]

    def benchmark_optimizers(
        self,
        candidates: List[str],
        target_size: int = 6,
        methods: Optional[List[str]] = None,
        runs: int = 3,
    ) -> Dict[str, Dict[str, float]]:
        """
        Benchmark optimizer performance with multiple runs.

        Returns statistics (mean, std, min, max) for each optimizer.
        """
        import numpy as np

        if methods is None:
            methods = list(self.list_available_agents().keys())

        results: Dict[str, List[float]] = {m: [] for m in methods}

        for run_idx in range(runs):
            logger.info(f"Benchmark run {run_idx + 1}/{runs}")

            result = self.optimize_parallel(
                candidates,
                methods=methods,
                target_size=target_size,
                verbose=False,
            )

            for method, metrics in result.performance_metrics.items():
                if metrics.success:
                    results[method].append(metrics.execution_time)

        # Calculate statistics
        stats = {}
        for method, times in results.items():
            if times:
                stats[method] = {
                    'mean': float(np.mean(times)),
                    'std': float(np.std(times)),
                    'min': float(np.min(times)),
                    'max': float(np.max(times)),
                    'runs': len(times),
                }
            else:
                stats[method] = {'error': 'No successful runs'}

        return stats


# =============================================================================
# Factory Registration
# =============================================================================

from .base_optimizer import BaseOptimizer


@OptimizerFactory.register('multi-agent', aliases=['ensemble', 'parallel'])
class MultiAgentBaseOptimizer(BaseOptimizer):
    """
    Multi-agent optimizer that runs multiple strategies in parallel.

    This optimizer coordinates multiple underlying optimizers and aggregates
    their results to find the best primer set.
    """

    def __init__(
        self,
        position_cache,
        fg_prefixes: List[str],
        fg_seq_lengths: List[int],
        bg_prefixes: Optional[List[str]] = None,
        bg_seq_lengths: Optional[List[int]] = None,
        config: Optional[OptimizerConfig] = None,
        **kwargs
    ):
        super().__init__(
            position_cache, fg_prefixes, fg_seq_lengths,
            bg_prefixes, bg_seq_lengths, config
        )
        self.orchestrator = MultiAgentOrchestrator(
            position_cache=position_cache,
            fg_prefixes=fg_prefixes,
            fg_seq_lengths=fg_seq_lengths,
            bg_prefixes=bg_prefixes,
            bg_seq_lengths=bg_seq_lengths,
        )
        self.methods = kwargs.get('methods', ['greedy', 'dominating-set', 'hybrid'])
        self.aggregation = kwargs.get(
            'aggregation',
            AggregationStrategy.BEST_SCORE
        )

    @property
    def name(self) -> str:
        return "multi-agent"

    @property
    def description(self) -> str:
        return "Multi-agent parallel optimizer with result aggregation"

    def optimize(
        self,
        candidates: List[str],
        target_size: Optional[int] = None,
        **kwargs
    ) -> OptimizationResult:
        """Run multi-agent optimization."""
        candidates = self._validate_candidates(candidates)
        target = target_size or self.config.target_set_size

        result = self.orchestrator.optimize_parallel(
            candidates=candidates,
            methods=self.methods,
            target_size=target,
            aggregation=self.aggregation,
            config=self.config,
            verbose=self.config.verbose,
        )

        return result.best_result


# =============================================================================
# CLI Entry Point
# =============================================================================

def main():
    """CLI entry point for multi-agent optimization."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Multi-Agent Primer Set Optimization"
    )
    parser.add_argument(
        '-j', '--json-file',
        required=True,
        help='Parameters JSON file'
    )
    parser.add_argument(
        '-m', '--methods',
        nargs='+',
        default=['greedy', 'dominating-set', 'hybrid'],
        help='Optimizer methods to use'
    )
    parser.add_argument(
        '-n', '--num-primers',
        type=int,
        default=6,
        help='Target number of primers'
    )
    parser.add_argument(
        '--aggregation',
        choices=['best_score', 'consensus', 'weighted_vote'],
        default='best_score',
        help='Result aggregation strategy'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Verbose output'
    )

    args = parser.parse_args()

    # Load parameters and run
    from . import parameter
    from . import pipeline as core_pipeline
    from .position_cache import PositionCache
    import pandas as pd

    parameter.json_file = args.json_file
    core_pipeline._initialize()

    # Load candidates
    step3_path = f"{parameter.data_dir}/step3_df.csv"
    step3_df = pd.read_csv(step3_path)
    candidates = step3_df['primer'].tolist()

    # Create position cache
    cache = PositionCache(core_pipeline.fg_prefixes, candidates)

    # Create orchestrator
    orchestrator = MultiAgentOrchestrator(
        position_cache=cache,
        fg_prefixes=core_pipeline.fg_prefixes,
        fg_seq_lengths=core_pipeline.fg_seq_lengths,
        bg_prefixes=getattr(core_pipeline, 'bg_prefixes', []),
        bg_seq_lengths=getattr(core_pipeline, 'bg_seq_lengths', []),
    )

    # Run optimization
    aggregation = AggregationStrategy(args.aggregation)

    result = orchestrator.optimize_parallel(
        candidates=candidates,
        methods=args.methods,
        target_size=args.num_primers,
        aggregation=aggregation,
        verbose=args.verbose,
    )

    # Print summary
    print(result.summary())

    # Print best primers
    print("\nSelected primers:")
    for primer in result.best_result.primers:
        print(f"  {primer}")


if __name__ == '__main__':
    main()
