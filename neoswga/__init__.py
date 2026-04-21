"""
NeoSWGA - Selective Whole Genome Amplification primer design.

For CLI usage:
    $ neoswga count-kmers -j params.json
    $ neoswga filter -j params.json
    $ neoswga score -j params.json
    $ neoswga optimize -j params.json

For library usage:
    >>> import neoswga
    >>> from neoswga import ReactionConditions, run_optimization
    >>> conditions = ReactionConditions(temp=30.0, polymerase='phi29', betaine_m=1.0)
    >>> result = run_optimization(method='hybrid', candidates=[...],
    ...                            fg_prefixes=['target'], fg_seq_lengths=[1_000_000],
    ...                            conditions=conditions, seed=42)
    >>> print(result.primers, result.metrics.fg_coverage)

Documentation: See README.md
"""

__version__ = "3.6.0"
__author__ = "Andreas Sjodin"
__email__ = "andreas.sjodin@foi.se"
__license__ = "AGPL-3.0-or-later"


# Public API — symbols a library consumer can import directly from
# ``neoswga`` without reaching into submodules. Added in Phase 16 to close
# the library-ergonomics gap surfaced by the post-Phase-15 code review.
# Imports are lazy (via __getattr__) so the module's import-time cost
# stays minimal for CLI-only consumers.
__all__ = [
    "__version__",
    "get_version",
    # Optimization entry points
    "run_optimization",
    "optimize_step4",
    "OptimizationResult",
    "OptimizationStatus",
    "PrimerSetMetrics",
    "OptimizerFactory",
    # Scientific configuration
    "ReactionConditions",
    "get_polymerase_processivity",
    "get_typical_amplicon_length",
    # Coverage computation
    "compute_per_prefix_coverage",
    "polymerase_extension_reach",
    # Position data
    "PositionCache",
]


def get_version():
    """Get version string."""
    return __version__


def print_info():
    """Print package information."""
    print(f"NeoSWGA version {__version__}")
    print(f"Author: {__author__}")
    print(f"License: {__license__}")
    print()
    print("Enhanced primer design for selective whole genome amplification")
    print("Documentation: https://github.com/FOI-Bioinformatics/neoswga")


# Lazy-import table: maps each public symbol to its source module so
# `from neoswga import X` resolves without loading heavy optimizer
# machinery for simple CLI consumers.
_LAZY_EXPORTS = {
    "run_optimization":             "neoswga.core.unified_optimizer",
    "optimize_step4":               "neoswga.core.unified_optimizer",
    "OptimizationResult":           "neoswga.core.base_optimizer",
    "OptimizationStatus":           "neoswga.core.base_optimizer",
    "PrimerSetMetrics":             "neoswga.core.base_optimizer",
    "OptimizerFactory":             "neoswga.core.optimizer_factory",
    "ReactionConditions":           "neoswga.core.reaction_conditions",
    "get_polymerase_processivity":  "neoswga.core.reaction_conditions",
    "get_typical_amplicon_length":  "neoswga.core.reaction_conditions",
    "compute_per_prefix_coverage":  "neoswga.core.coverage",
    "polymerase_extension_reach":   "neoswga.core.coverage",
    "PositionCache":                "neoswga.core.position_cache",
}


def __getattr__(name):
    """PEP 562 lazy import for the public API.

    Keeps `import neoswga` cheap (no heavy optimizer imports) while still
    letting library users write `from neoswga import run_optimization`.
    """
    if name in _LAZY_EXPORTS:
        from importlib import import_module
        module = import_module(_LAZY_EXPORTS[name])
        value = getattr(module, name)
        # Cache on the module so subsequent lookups skip the import.
        globals()[name] = value
        return value
    raise AttributeError(f"module 'neoswga' has no attribute {name!r}")


def __dir__():
    """Surface the lazy exports in `dir(neoswga)` so IDEs and interactive
    Python shells auto-complete the public API."""
    return sorted(set(list(globals().keys()) + list(_LAZY_EXPORTS.keys())))
