"""Coverage for neoswga.core.validation (the installation validation suite).

Exercises the lightweight, self-contained suite checks (position cache and
adaptive GC filter) plus the result dataclass — not the full jellyfish/MILP
end-to-end path, which is covered by the integration suite.
"""

from neoswga.core.validation import ValidationResult, ValidationSuite


def test_validation_result_dataclass():
    r = ValidationResult(test_name="x", passed=True, runtime=0.1, details={"a": 1})
    assert r.passed is True
    assert r.error is None
    assert r.details["a"] == 1


def test_position_cache_check_runs():
    suite = ValidationSuite(verbose=False)
    passed, details = suite.test_position_cache()
    assert isinstance(passed, bool)
    assert isinstance(details, dict)


def test_adaptive_gc_filter_check_runs():
    suite = ValidationSuite(verbose=False)
    passed, details = suite.test_adaptive_gc_filter()
    assert isinstance(passed, bool)
    assert isinstance(details, dict)
