"""Tests for safe_pickle module -- verifies no unsafe fallbacks exist."""
import pickle
import tempfile
import os
import pytest
from unittest.mock import patch


def test_safe_load_raises_on_disallowed_class():
    """safe_load must raise, not fall back, when encountering disallowed classes."""
    from neoswga.core.safe_pickle import safe_load

    # Create a pickle containing a disallowed class (os.getcwd)
    payload = pickle.dumps({"action": os.getcwd})

    with tempfile.NamedTemporaryFile(suffix='.p', delete=False) as f:
        f.write(payload)
        tmp_path = f.name

    try:
        with pytest.raises(pickle.UnpicklingError):
            safe_load(tmp_path, context='sklearn_model')
    finally:
        os.unlink(tmp_path)


def test_no_unsafe_load_with_warning_function():
    """The unsafe_load_with_warning function should not exist."""
    from neoswga.core import safe_pickle
    assert not hasattr(safe_pickle, 'unsafe_load_with_warning'), \
        "unsafe_load_with_warning must be removed -- it defeats safe_load"


def test_background_bloom_filter_load_no_fallback():
    """BackgroundBloomFilter.load must not fall back to raw pickle.load."""
    import inspect
    from neoswga.core.background_filter import BackgroundBloomFilter
    source = inspect.getsource(BackgroundBloomFilter.load)
    assert 'pickle.load' not in source, \
        "BackgroundBloomFilter.load must not contain pickle.load fallback"


def test_sampled_genome_index_load_no_fallback():
    """SampledGenomeIndex.load must not fall back to raw pickle.load."""
    import inspect
    from neoswga.core.background_filter import SampledGenomeIndex
    source = inspect.getsource(SampledGenomeIndex.load)
    assert 'pickle.load' not in source, \
        "SampledGenomeIndex.load must not contain pickle.load fallback"


def test_rf_preprocessing_load_no_fallback():
    """load_model_safely must not fall back to raw pickle.load."""
    import inspect
    from neoswga.core.rf_preprocessing import load_model_safely
    source = inspect.getsource(load_model_safely)
    assert 'pickle.load' not in source and 'unsafe_load' not in source, \
        "load_model_safely must not contain pickle.load or unsafe_load fallback"
