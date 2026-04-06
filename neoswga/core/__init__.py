"""NeoSWGA core modules for primer design and optimization."""

__all__ = [
    "ReactionConditions",
    "PositionCache",
    "MechanisticModel",
    "MechanisticEffects",
    "get_standard_conditions",
    "get_enhanced_conditions",
    "calculate_tm_with_salt",
    "calculate_tm_basic",
]


def __getattr__(name):
    """Lazy imports for public API symbols."""
    if name == "ReactionConditions":
        from neoswga.core.reaction_conditions import ReactionConditions
        return ReactionConditions
    if name == "get_standard_conditions":
        from neoswga.core.reaction_conditions import get_standard_conditions
        return get_standard_conditions
    if name == "get_enhanced_conditions":
        from neoswga.core.reaction_conditions import get_enhanced_conditions
        return get_enhanced_conditions
    if name == "PositionCache":
        from neoswga.core.position_cache import PositionCache
        return PositionCache
    if name == "MechanisticModel":
        from neoswga.core.mechanistic_model import MechanisticModel
        return MechanisticModel
    if name == "MechanisticEffects":
        from neoswga.core.mechanistic_model import MechanisticEffects
        return MechanisticEffects
    if name == "calculate_tm_with_salt":
        from neoswga.core.thermodynamics import calculate_tm_with_salt
        return calculate_tm_with_salt
    if name == "calculate_tm_basic":
        from neoswga.core.thermodynamics import calculate_tm_basic
        return calculate_tm_basic
    raise AttributeError(f"module 'neoswga.core' has no attribute {name!r}")
