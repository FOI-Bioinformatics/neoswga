"""Test that GA population evaluation does not crash due to pickling."""
import inspect
import pytest


def test_evaluate_population_does_not_use_process_pool():
    """GA must not use ProcessPoolExecutor with bound methods."""
    from neoswga.core.genetic_algorithm import PrimerSetGA
    source = inspect.getsource(PrimerSetGA._evaluate_population)
    assert 'ProcessPoolExecutor' not in source, \
        "ProcessPoolExecutor with bound methods fails on macOS spawn start method"


def test_dimer_check_excludes_self_pairs():
    """Dimer compatibility check should not compare a primer with itself."""
    from neoswga.core.genetic_algorithm import PrimerSetGA
    source = inspect.getsource(PrimerSetGA._check_dimer_compatibility)
    assert 'range(i + 1,' in source or 'range(i+1,' in source, \
        "Inner loop should start at i+1, not i, to exclude self-pairs"
