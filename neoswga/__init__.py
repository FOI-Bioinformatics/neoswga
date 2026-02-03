"""
NeoSWGA - Enhanced Selective Whole Genome Amplification

Next-generation primer design with advanced thermodynamics and optimization.

Key Features:
- Complete SantaLucia thermodynamics with salt corrections
- Reaction condition optimization (additives, polymerases)
- Adaptive k-mer selection (6-15mers)
- Advanced secondary structure prediction
- Genetic algorithm optimization
- Network-based coverage analysis
- Agent-based replication simulation
- GPU acceleration with CuPy
- Deep learning transformer embeddings

Example usage:
    >>> import neoswga
    >>> from neoswga import thermodynamics as thermo
    >>> tm = thermo.calculate_tm_with_salt('ATCGATCG', na_conc=50)
    >>> print(f"Tm: {tm:.1f}°C")

For CLI usage:
    $ neoswga-enhanced quick-design --fg genome.fasta --output results/

Documentation: See README.md and ENHANCEMENTS_QUICKSTART.md
"""

__version__ = "3.0.0"
__author__ = "Andreas Sjodin"
__email__ = "your-email@domain.com"
__license__ = "MIT"

# Import key modules for convenient access
import sys
import os

# Add parent directory to path for src imports
_parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _parent_dir not in sys.path:
    sys.path.insert(0, _parent_dir)

from neoswga.core import (
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
)

# Convenience imports
from neoswga.core.thermodynamics import (
    calculate_tm_with_salt,
    calculate_free_energy,
    gc_content,
)
from neoswga.core.reaction_conditions import (
    ReactionConditions,
    get_standard_conditions,
    get_enhanced_conditions,
)

__all__ = [
    # Version info
    "__version__",
    "__author__",
    "__email__",
    "__license__",

    # Core modules
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

    # Convenience functions
    "calculate_tm_with_salt",
    "calculate_free_energy",
    "gc_content",
    "ReactionConditions",
    "get_standard_conditions",
    "get_enhanced_conditions",
    "UnifiedPipeline",
    "PipelineConfig",
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
    print("Documentation: https://github.com/yourusername/neoswga")
