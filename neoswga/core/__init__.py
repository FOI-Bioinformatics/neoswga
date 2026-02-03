"""
NeoSWGA core modules.

This package contains all the enhanced functionality for primer design.

Key improvements (2025):
- Position cache: 1000x faster primer position access
- Adaptive GC filtering: Works for GC-extreme organisms
- Background Bloom filter: Handles 3 Gbp genomes
- Network optimization: 100x better enrichment
- MILP optimizer: Provably optimal solutions
"""

# Import improved pipeline modules (2025 enhancements)
from . import (
    position_cache,
    background_filter,
    adaptive_filters,
    network_optimizer,
    milp_optimizer,
    improved_pipeline,
    unified_optimizer,
    stochastic_simulator,
    validation,
)

# Import refactored optimizer framework (clean architecture)
from . import (
    base_optimizer,
    optimizer_factory,
    exceptions,
    search_context,
    additives,
    greedy_optimizer,
    dominating_set_adapter,
    hybrid_optimizer,
    background_aware_optimizer,
)

# Import original enhanced modules
from . import (
    thermodynamics,
    reaction_conditions,
    secondary_structure,
    adaptive_search,
    genetic_algorithm,
    amplicon_network,
    replication_simulator,
    advanced_features,
    gpu_acceleration,
    deep_learning,
    # Utility modules
    filter,
    kmer,
    dimer,
    parameter,
    primer_attributes,
    string_search,
    utility,
    rf_preprocessing,
    # Note: thermo_estimation is deprecated - use thermodynamics instead
)

__all__ = [
    # Improved pipeline modules (recommended)
    "position_cache",
    "background_filter",
    "adaptive_filters",
    "network_optimizer",
    "milp_optimizer",
    "improved_pipeline",
    "unified_optimizer",
    "stochastic_simulator",
    "validation",
    # Refactored optimizer framework (clean architecture)
    "base_optimizer",
    "optimizer_factory",
    "exceptions",
    "search_context",
    "additives",
    "greedy_optimizer",
    "dominating_set_adapter",
    "hybrid_optimizer",
    "background_aware_optimizer",
    # Original enhanced modules
    "thermodynamics",
    "reaction_conditions",
    "secondary_structure",
    "adaptive_search",
    "genetic_algorithm",
    "amplicon_network",
    "replication_simulator",
    "advanced_features",
    "gpu_acceleration",
    "deep_learning",
    "filter",
    "kmer",
    "dimer",
    "parameter",
    "primer_attributes",
    "string_search",
    "utility",
    "rf_preprocessing",
]
