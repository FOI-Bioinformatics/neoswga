"""Derived views over the registry, shaped like the tables they replace.

Each function returns exactly what one existing hardcoded table used to contain, so
consumers can swap a literal for a call without changing behaviour. Keeping the old
shapes is deliberate: it lets the registry land with the entire existing test suite
untouched, which is the proof that the seeding was faithful.

New code should read ``PolymeraseSpec`` fields directly rather than going through
these.
"""

from typing import Dict, Tuple

from neoswga.core.registry.polymerases import POLYMERASES

__all__ = [
    "as_characteristics",
    "hard_temp_ranges",
    "warn_temp_ranges",
    "primer_length_ranges",
    "default_primer_lengths",
    "mg_defaults",
    "processivity_map",
    "mechanistic_enzyme_params",
    "additive_optimizer_constraints",
]


def as_characteristics() -> Dict[str, dict]:
    """Shape of ``reaction_conditions.POLYMERASE_CHARACTERISTICS``."""
    return {
        key: {
            "name": spec.display_name,
            "temp_range": spec.temp_optimal_range,
            "optimal_temp": spec.optimal_temp,
            "processivity": spec.processivity_bp,
            "typical_amplicon_length": spec.typical_amplicon_bp,
            "strand_displacement": spec.strand_displacement,
            "exonuclease": spec.exonuclease,
            "error_rate": spec.error_rate,
            "requires_primer": True,
            "description": spec.description,
        }
        for key, spec in POLYMERASES.items()
    }


def hard_temp_ranges() -> Dict[str, Tuple[float, float]]:
    """Shape of the ``polymerase_temp_ranges`` local in ``ReactionConditions._validate``."""
    return {key: spec.temp_hard_range for key, spec in POLYMERASES.items()}


def warn_temp_ranges() -> Dict[str, Tuple[float, float]]:
    """Shape of ``param_validator.POLYMERASE_TEMP_RANGES``."""
    return {key: spec.temp_warn_range for key, spec in POLYMERASES.items()}


def primer_length_ranges() -> Dict[str, Tuple[int, int]]:
    """Shape of ``param_validator.POLYMERASE_PRIMER_LENGTHS`` (and its two copies)."""
    return {key: spec.primer_length_range for key, spec in POLYMERASES.items()}


def default_primer_lengths() -> Dict[str, int]:
    """Shape of the scalar defaults in ``cli/commands.py``.

    Distinct from ``primer_length_ranges``: this is a single suggested length, not
    a design window.
    """
    return {key: spec.default_primer_length for key, spec in POLYMERASES.items()}


def mg_defaults() -> Dict[str, float]:
    """Shape of ``parameter.MG_DEFAULTS_MM``."""
    return {key: spec.mg_default_mm for key, spec in POLYMERASES.items()}


def processivity_map() -> Dict[str, int]:
    """Single-molecule processivity, in bp."""
    return {key: spec.processivity_bp for key, spec in POLYMERASES.items()}


def mechanistic_enzyme_params() -> Dict[str, dict]:
    """Per-polymerase block of ``mechanistic_params.MECHANISTIC_MODEL_PARAMS['enzyme']``."""
    return {
        key: {
            "optimal_temp": spec.optimal_temp,
            "processivity": spec.processivity_bp,
            "processivity_step": spec.processivity_step,
            "extension_rate": spec.extension_rate_nt_s,
            "dmso_threshold": spec.dmso_response["threshold"],
            "dmso_mild_coef": spec.dmso_response["mild_coef"],
            "dmso_steep_coef": spec.dmso_response["steep_coef"],
        }
        for key, spec in POLYMERASES.items()
    }


def additive_optimizer_constraints() -> Dict[str, dict]:
    """Shape of ``additive_optimizer.AdditiveOptimizer.POLYMERASE_CONSTRAINTS``.

    Note ``optimal_temp`` here is ``preset_reaction_temp``, which differs from
    ``optimal_temp`` for bst (60.0 vs 63.0).
    """
    return {
        key: {
            "max_dmso": spec.additive_limits["dmso_percent"],
            "max_betaine": spec.additive_limits["betaine_m"],
            "optimal_temp": spec.preset_reaction_temp,
            "max_primer_length_base": spec.additive_opt_base_length,
        }
        for key, spec in POLYMERASES.items()
    }
