#!/usr/bin/env python3
"""
Performance benchmarks for NeoSWGA optimizer components.

Measures and compares performance of critical operations:
1. Dimer matrix computation (parallel vs sequential)
2. Network edge building (incremental vs rebuild)
3. BoundedScoreCache hit rates
4. Multi-agent orchestrator thread pool reuse

Usage:
    python benchmarks/performance_benchmarks.py
"""

import time
import random
import string
import numpy as np
from typing import List, Callable
import sys


def generate_random_primers(n: int, length: int = 10) -> List[str]:
    """Generate random primer sequences."""
    bases = 'ACGT'
    return [''.join(random.choices(bases, k=length)) for _ in range(n)]


def generate_random_positions(n_primers: int, genome_length: int = 1_000_000,
                              sites_per_primer: int = 100) -> dict:
    """Generate random binding positions for primers."""
    positions = {}
    for i in range(n_primers):
        primer = f"PRIMER_{i}"
        positions[primer] = {
            'forward': np.random.randint(0, genome_length, size=sites_per_primer),
            'reverse': np.random.randint(0, genome_length, size=sites_per_primer),
        }
    return positions


def benchmark(func: Callable, name: str, iterations: int = 5) -> dict:
    """Run benchmark and return statistics."""
    times = []
    for _ in range(iterations):
        start = time.perf_counter()
        result = func()
        elapsed = time.perf_counter() - start
        times.append(elapsed)

    return {
        'name': name,
        'mean': np.mean(times),
        'std': np.std(times),
        'min': np.min(times),
        'max': np.max(times),
        'iterations': iterations,
    }


def print_benchmark_result(result: dict):
    """Print benchmark result in a formatted way."""
    print(f"  {result['name']}: {result['mean']*1000:.2f}ms "
          f"(+/- {result['std']*1000:.2f}ms, n={result['iterations']})")


def benchmark_dimer_matrix():
    """Benchmark dimer matrix computation methods."""
    print("\n" + "=" * 60)
    print("BENCHMARK: Dimer Matrix Computation")
    print("=" * 60)

    from neoswga.core import dimer, parameter

    # Ensure parameter is initialized with default
    if not hasattr(parameter, 'max_self_dimer_bp'):
        parameter.max_self_dimer_bp = 4

    for n_primers in [50, 100, 200, 500]:
        print(f"\nPrimers: {n_primers}")
        primers = generate_random_primers(n_primers)

        # Original (sequential)
        result_seq = benchmark(
            lambda p=primers: dimer.heterodimer_matrix(p),
            "heterodimer_matrix (original)",
            iterations=3
        )
        print_benchmark_result(result_seq)

        # Fast sequential
        result_fast = benchmark(
            lambda p=primers: dimer.heterodimer_matrix_fast(p),
            "heterodimer_matrix_fast",
            iterations=3
        )
        print_benchmark_result(result_fast)
        speedup = result_seq['mean'] / result_fast['mean'] if result_fast['mean'] > 0 else 0
        print(f"    Speedup: {speedup:.2f}x")

        # Parallel (for larger sets)
        if n_primers >= 100:
            result_parallel = benchmark(
                lambda p=primers: dimer.heterodimer_matrix_parallel(p),
                "heterodimer_matrix_parallel",
                iterations=3
            )
            print_benchmark_result(result_parallel)
            speedup = result_seq['mean'] / result_parallel['mean'] if result_parallel['mean'] > 0 else 0
            print(f"    Speedup: {speedup:.2f}x")


def benchmark_bounded_cache():
    """Benchmark BoundedScoreCache performance."""
    print("\n" + "=" * 60)
    print("BENCHMARK: BoundedScoreCache")
    print("=" * 60)

    from neoswga.core.greedy_optimizer import BoundedScoreCache

    cache = BoundedScoreCache(max_size=10000)
    primers = generate_random_primers(100)

    # Populate cache
    n_entries = 5000
    keys = [frozenset(random.sample(primers, 6)) for _ in range(n_entries)]

    print(f"\nPopulating cache with {n_entries} entries...")
    start = time.perf_counter()
    for key in keys:
        cache.set(key, random.random())
    populate_time = time.perf_counter() - start
    print(f"  Population time: {populate_time*1000:.2f}ms")

    # Cache lookups
    n_lookups = 10000
    lookup_keys = [random.choice(keys) for _ in range(n_lookups)]

    print(f"\nPerforming {n_lookups} lookups...")
    start = time.perf_counter()
    for key in lookup_keys:
        cache.get(key)
    lookup_time = time.perf_counter() - start

    print(f"  Lookup time: {lookup_time*1000:.2f}ms ({lookup_time/n_lookups*1e6:.2f}us per lookup)")
    print(f"  Hit rate: {cache.hit_rate()*100:.1f}%")
    print(f"  Cache size: {len(cache)}")


def benchmark_gini_coefficient():
    """Benchmark Gini coefficient calculation."""
    print("\n" + "=" * 60)
    print("BENCHMARK: Gini Coefficient Calculation")
    print("=" * 60)

    from neoswga.core.network_optimizer import AmplificationNetwork, BindingSite

    for n_sites in [100, 1000, 5000, 10000]:
        print(f"\nSites: {n_sites}")

        # Create network with random sites
        network = AmplificationNetwork()
        for i in range(n_sites):
            pos = random.randint(0, 1_000_000)
            strand = '+' if random.random() > 0.5 else '-'
            site = BindingSite(position=pos, strand=strand, primer=f"P{i}")
            network.binding_sites.append(site)
            network.graph.add_node(site)

        # Benchmark
        result = benchmark(
            lambda n=network: n._calculate_gini_coefficient(),
            f"Gini ({n_sites} sites)",
            iterations=5
        )
        print_benchmark_result(result)


def benchmark_network_edge_building():
    """Benchmark network edge building methods."""
    print("\n" + "=" * 60)
    print("BENCHMARK: Network Edge Building")
    print("=" * 60)

    from neoswga.core.network_optimizer import AmplificationNetwork, BindingSite

    for n_sites in [500, 1000, 2000]:
        print(f"\nSites: {n_sites}")

        # Create sites
        sites = []
        for i in range(n_sites):
            pos = random.randint(0, 1_000_000)
            strand = '+' if random.random() > 0.5 else '-'
            site = BindingSite(position=pos, strand=strand, primer=f"P{i % 10}")
            sites.append(site)

        # Full build
        def full_build():
            network = AmplificationNetwork(max_extension=70000)
            for site in sites:
                network.binding_sites.append(site)
                network.graph.add_node(site)
            network.build_edges()
            return network

        result_full = benchmark(full_build, "build_edges (full)", iterations=3)
        print_benchmark_result(result_full)

        # Incremental build (add 50 sites to existing n-50)
        def incremental_build():
            network = AmplificationNetwork(max_extension=70000)
            # Add base sites
            for site in sites[:-50]:
                network.binding_sites.append(site)
                network.graph.add_node(site)
            network.build_edges()

            # Add new sites incrementally
            new_sites = []
            for site in sites[-50:]:
                network.binding_sites.append(site)
                network.graph.add_node(site)
                new_sites.append(site)
            network.add_edges_for_sites(new_sites)
            return network

        result_incr = benchmark(incremental_build, "add_edges_for_sites (incremental)", iterations=3)
        print_benchmark_result(result_incr)

    # Real-world simulation: adding primers one at a time
    print("\n--- Real-World Simulation: Adding 10 primers incrementally ---")
    n_sites_per_primer = 200

    def rebuild_each_time():
        """Old approach: rebuild all edges after each primer."""
        network = AmplificationNetwork(max_extension=70000)
        for primer_idx in range(10):
            # Add sites for this primer
            for _ in range(n_sites_per_primer):
                pos = random.randint(0, 1_000_000)
                strand = '+' if random.random() > 0.5 else '-'
                site = BindingSite(position=pos, strand=strand, primer=f"P{primer_idx}")
                network.binding_sites.append(site)
                network.graph.add_node(site)
            network._index_dirty = True
            # Clear and rebuild all edges (simulating old behavior)
            network.graph.remove_edges_from(list(network.graph.edges()))
            network.build_edges()
        return network

    def incremental_each_time():
        """New approach: incremental edge building."""
        network = AmplificationNetwork(max_extension=70000)
        for primer_idx in range(10):
            new_sites = []
            for _ in range(n_sites_per_primer):
                pos = random.randint(0, 1_000_000)
                strand = '+' if random.random() > 0.5 else '-'
                site = BindingSite(position=pos, strand=strand, primer=f"P{primer_idx}")
                network.binding_sites.append(site)
                network.graph.add_node(site)
                new_sites.append(site)
            network.add_edges_for_sites(new_sites)
        return network

    result_rebuild = benchmark(rebuild_each_time, "Rebuild all edges each primer", iterations=3)
    print_benchmark_result(result_rebuild)

    result_incremental = benchmark(incremental_each_time, "Incremental edges each primer", iterations=3)
    print_benchmark_result(result_incremental)

    speedup = result_rebuild['mean'] / result_incremental['mean'] if result_incremental['mean'] > 0 else 0
    print(f"    Speedup: {speedup:.2f}x")


def benchmark_orchestrator_pool_reuse():
    """Benchmark thread pool reuse in MultiAgentOrchestrator."""
    print("\n" + "=" * 60)
    print("BENCHMARK: Orchestrator Thread Pool Reuse")
    print("=" * 60)

    from neoswga.core.multi_agent_optimizer import MultiAgentOrchestrator

    # Mock position cache
    class MockCache:
        def get_positions(self, prefix, primer, strand='both'):
            return np.random.randint(0, 1000000, size=50)

    cache = MockCache()

    # Create orchestrator
    orchestrator = MultiAgentOrchestrator(
        position_cache=cache,
        fg_prefixes=['target'],
        fg_seq_lengths=[1000000],
    )

    # Test executor reuse
    n_calls = 10
    print(f"\nGetting executor {n_calls} times (should reuse)...")

    start = time.perf_counter()
    for _ in range(n_calls):
        executor = orchestrator._get_executor(4)
    elapsed = time.perf_counter() - start

    print(f"  Time for {n_calls} calls: {elapsed*1000:.2f}ms")
    print(f"  Per call: {elapsed/n_calls*1000:.3f}ms")
    print(f"  Executor reused: {orchestrator._executor is not None}")

    # Cleanup
    orchestrator.shutdown()


def run_all_benchmarks():
    """Run all performance benchmarks."""
    print("=" * 60)
    print("NeoSWGA Performance Benchmarks")
    print("=" * 60)
    print(f"Python: {sys.version}")
    print(f"NumPy: {np.__version__}")

    benchmark_dimer_matrix()
    benchmark_bounded_cache()
    benchmark_gini_coefficient()
    benchmark_network_edge_building()
    benchmark_orchestrator_pool_reuse()

    print("\n" + "=" * 60)
    print("Benchmarks Complete")
    print("=" * 60)


if __name__ == '__main__':
    run_all_benchmarks()
