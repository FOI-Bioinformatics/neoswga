"""The resolved configuration `run_optimization_from_config` takes.

Extracted from `core/unified_optimizer.py` on 2026-10-28, which sat at 1646
lines against a 1650 budget. `unified_optimizer` re-exports this name, so
`from neoswga.core.unified_optimizer import OptimizationConfig` still resolves
and `dataclasses.fields()` still sees the same class.

This is NOT `base_optimizer.OptimizerConfig`, which is what an optimizer
instance holds; `test_optimization_config_fields_are_forwarded.py` exists
because the two were confused once already.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class OptimizationConfig:
    """
    Unified configuration for optimization runs.

    Combines all parameters needed for optimization in one place.
    """

    # Optimizer selection
    method: str = "hybrid"

    # Target parameters
    target_set_size: int = 6
    max_iterations: int = 100

    # File paths
    data_dir: str = "."
    step3_file: str = "step3_df.csv"
    output_file: str = "step4_improved_df.csv"

    # Genome info (loaded from pipeline if not provided)
    fg_prefixes: list[str] | None = None
    fg_seq_lengths: list[int] | None = None
    bg_prefixes: list[str] | None = None
    bg_seq_lengths: list[int] | None = None

    # Performance options
    use_position_cache: bool = True
    use_background_filter: bool = True

    # Output options
    verbose: bool = True
    quiet: bool = False

    # Advanced options
    allow_dimer_relaxation: bool = False
    refinement_method: str = "network"
    swap_max_evaluations: int = 10000
    stage1_objective_width: int | None = None
    stage2_tm_weight: float | None = None
    stage2_background_weight: float | None = None
    swap_max_seconds: float = 10.0
    uniformity_weight: float = 0.0
    minimize_primers: bool = False
    target_coverage: float = 0.70
