"""Test that the public API is importable and stable."""


def test_top_level_exports():
    import neoswga
    assert hasattr(neoswga, '__version__')
    assert hasattr(neoswga, '__all__')
    assert '__version__' in neoswga.__all__


def test_core_lazy_imports():
    from neoswga.core import ReactionConditions
    from neoswga.core import PositionCache
    from neoswga.core import MechanisticModel
    from neoswga.core import MechanisticEffects
    from neoswga.core import calculate_tm_with_salt
    from neoswga.core import calculate_tm_basic
    from neoswga.core import get_standard_conditions
    from neoswga.core import get_enhanced_conditions
    assert callable(calculate_tm_with_salt)
    assert callable(calculate_tm_basic)


def test_core_all_is_defined():
    import neoswga.core
    assert hasattr(neoswga.core, '__all__')
    assert 'ReactionConditions' in neoswga.core.__all__
    assert len(neoswga.core.__all__) == 8


def test_core_invalid_attr_raises():
    import neoswga.core
    import pytest
    with pytest.raises(AttributeError):
        _ = neoswga.core.NonexistentThing
