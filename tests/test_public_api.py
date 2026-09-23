"""Test that the public API is importable and stable."""


def test_top_level_exports():
    import neoswga

    assert hasattr(neoswga, "__version__")
    assert hasattr(neoswga, "__all__")
    assert "__version__" in neoswga.__all__


def test_core_lazy_imports():
    from neoswga.core import (
        MechanisticEffects,
        MechanisticModel,
        PositionCache,
        ReactionConditions,
        calculate_tm_basic,
        calculate_tm_with_salt,
        get_enhanced_conditions,
        get_standard_conditions,
    )

    assert callable(calculate_tm_with_salt)
    assert callable(calculate_tm_basic)


def test_core_all_is_defined():
    import neoswga.core

    assert hasattr(neoswga.core, "__all__")
    assert "ReactionConditions" in neoswga.core.__all__
    assert len(neoswga.core.__all__) == 8


def test_core_invalid_attr_raises():
    import pytest

    import neoswga.core

    with pytest.raises(AttributeError):
        _ = neoswga.core.NonexistentThing
