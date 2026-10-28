"""The hybrid optimizer's result record.

Extracted from `core/hybrid_optimizer.py` on 2026-10-28, which sat at 1646
lines against a 1650 budget -- the minimum headroom the ceiling ratchet allows,
so the next correct one-line edit there would have failed CI. `hybrid_optimizer`
re-exports this name, so existing imports and `isinstance` are unaffected.

It moved rather than `HybridOptimizer` or `HybridBaseOptimizer` because it is
self-contained: a dataclass that references no module-level name in the module
it came from.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class HybridResult:
    """Result from hybrid optimization"""

    # Final primer set
    primers: list[str]

    # Stage 1 (Coverage) results
    stage1_primers: list[str]
    stage1_coverage: float
    stage1_regions_covered: int

    # Stage 2 (Network) results
    stage2_primers: list[str]
    stage2_connectivity: float
    stage2_predicted_amplification: float
    stage2_largest_component: int

    # Overall metrics
    final_coverage: float
    final_connectivity: float
    final_predicted_amplification: float

    # Stage-1 primers in greedy selection order. `stage1_primers` comes from a
    # set and has no meaningful order; this does. Set-cover picks the largest
    # marginal gain at each step, so its first N picks are the same N whatever
    # `stage1_count` was, and coverage over that prefix is non-decreasing in N.
    # Stage 2 is floored against it, which is what makes final coverage
    # monotone in the requested count.
    stage1_ordered_primers: list[str] = field(default_factory=list)

    # Simulation validation (optional)
    simulation_fitness: object | None = None  # SimulationFitness if validated

    # Metadata
    runtime_stage1: float = 0.0
    runtime_stage2: float = 0.0
    runtime_simulation: float = 0.0
    total_runtime: float = 0.0

    def __str__(self):
        base = f"""Hybrid Optimization Result:
  Final primers: {len(self.primers)}
  Coverage: {self.final_coverage:.1%} (estimated, binned)
  Connectivity: {self.final_connectivity:.2f}
  Amplification: {self.final_predicted_amplification:.1f}×"""

        if self.simulation_fitness:
            base += f"""
  Simulation validation:
    Coverage: {self.simulation_fitness.mean_coverage:.1%} ± {self.simulation_fitness.std_coverage:.1%}
    Uniformity: {self.simulation_fitness.coverage_uniformity:.2f}
    Fitness: {self.simulation_fitness.fitness_score:.3f}"""

        base += f"\n  Runtime: {self.total_runtime:.2f}s"
        return base
