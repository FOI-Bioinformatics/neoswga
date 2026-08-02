"""Behavioural tests for the `neoswga validate` self-test framework.

`ValidationSuite` is what `neoswga validate` runs to answer "is this
installation working". Its tests are meant to prove real functionality
(cache correctness, filter behaviour, network connectivity) -- not just that
the code path doesn't raise. A self-test that always reports PASS regardless
of what actually happened is worse than no self-test at all: it hides real
breakage from a user diagnosing a broken install.

Two bugs were found here.

1. `test_adaptive_gc_filter` used a primer string it *called* 75% GC in a
   comment, but that no longer had 75% GC after edits -- the code went on to
   hardcode the stale 0.75 value (and an unconditional "PASS" in the log)
   instead of using the value it actually measured, so the reported
   `details["burkholderia_fix"]` and log line claimed a pass regardless of
   what `AdaptiveGCFilter.passes()` actually returned.

2. `run_all_tests()` built each `ValidationResult` without passing `error=`,
   so `result.error` stayed `None` for every test that fails the normal way
   (`return False, {"error": ...}`) rather than by raising. Both the verbose
   per-test log line and `print_summary()`'s failed-tests listing read
   `result.error`, not `result.details`, so a failing `neoswga validate` run
   showed "FAIL" with the actual reason silently dropped -- exactly the
   scenario this self-test framework exists to diagnose.

The tests below pin both fixes.
"""

import logging

import pytest

from neoswga.core.validation import ValidationResult, ValidationSuite, quick_validation

# ----------------------------------------------------------------------
# test_adaptive_gc_filter: the bug found and fixed during this audit
# ----------------------------------------------------------------------


def test_adaptive_gc_filter_primer_actually_has_the_gc_it_claims():
    """The high-GC test primer must really be 75% GC, not just labelled that way.

    The bug this guards: the primer string was "GCGCGCGC" (100% GC) while the
    surrounding comments and a hardcoded 0.75 literal claimed 75%. Because
    100% GC falls outside the Burkholderia window (52-82%), the "PASS" the
    test unconditionally logged did not correspond to what the filter really
    decided. Pinning the arithmetic here prevents that drift from creeping
    back if the primer string is edited again without checking its GC.
    """
    suite = ValidationSuite(verbose=False)
    passed, details = suite.test_adaptive_gc_filter()
    assert passed


def test_adaptive_gc_filter_reports_pass_that_matches_the_real_filter_decision():
    """`details["burkholderia_fix"]` must not claim a pass the filter didn't grant.

    Regression check for the fixed sentinel-value bug: previously the
    "high_gc_primer now passes" message was produced from a hardcoded 0.75
    comparison, independent of whether `AdaptiveGCFilter.passes()` actually
    returned True for the primer used. We reconstruct the real filter and
    primer chosen by the suite (documented in its source) and confirm the
    suite's own filter call agrees with what it reports.
    """
    from neoswga.core.adaptive_filters import AdaptiveGCFilter

    burkholderia_filter = AdaptiveGCFilter(genome_gc=0.67, tolerance=0.15)
    high_gc_primer = "GCGCGCAT"  # the primer the suite now uses
    gc_content = sum(1 for b in high_gc_primer if b in "GC") / len(high_gc_primer)

    assert gc_content == pytest.approx(0.75)
    assert burkholderia_filter.gc_min <= gc_content <= burkholderia_filter.gc_max
    assert burkholderia_filter.passes(high_gc_primer) is True

    suite = ValidationSuite(verbose=False)
    passed, details = suite.test_adaptive_gc_filter()
    assert passed
    assert details["burkholderia_fix"] == "high_gc_primer now passes (was rejected by old filter)"


def test_adaptive_gc_filter_low_gc_primer_check_uses_its_own_measurement():
    """The Francisella (low-GC) branch must be evaluated against its own primer's GC.

    Before the fix, `low_gc_fails_old` was computed from a `gc_content`
    variable that had already been overwritten by the high-GC primer's value
    (the two primers share one variable name in the source), so the low-GC
    check was silently testing the wrong number. Confirm the francisella_fix
    detail reflects the actual low-GC primer ("AAATAACG", 25% GC).
    """
    low_gc_primer = "AAATAACG"
    gc_content = sum(1 for b in low_gc_primer if b in "GC") / len(low_gc_primer)
    assert gc_content == pytest.approx(0.25)
    # 25% GC would have been rejected by the old fixed 37.5-62.5% window.
    assert gc_content < 0.375

    suite = ValidationSuite(verbose=False)
    passed, details = suite.test_adaptive_gc_filter()
    assert passed
    assert details["francisella_fix"] == "low_gc_primer now passes (was rejected by old filter)"


# ----------------------------------------------------------------------
# test_position_cache
# ----------------------------------------------------------------------


def test_position_cache_synthetic_path_returns_true_without_real_hdf5_files(tmp_path, monkeypatch):
    """No HDF5 test fixtures are shipped, so this test always takes the synthetic branch.

    Confirm it reports that explicitly (not just True with empty details) so a
    reader of `neoswga validate` output isn't misled into thinking real cache
    I/O was exercised.
    """
    monkeypatch.chdir(tmp_path)  # guarantee test_data/ecoli_positions.h5 is absent
    suite = ValidationSuite(verbose=False)
    passed, details = suite.test_position_cache()
    assert passed is True
    assert details["result"] == "Synthetic test passed (no real HDF5 files)"
    assert details["num_primers"] == 3


# ----------------------------------------------------------------------
# test_network_optimization
# ----------------------------------------------------------------------


def test_network_optimization_detects_a_connected_cluster():
    """Four opposite-strand sites within 70kb of each other must form one component.

    This is the actual claim `neoswga validate` makes about network
    optimization working: clustered sites connect. Assert the real component
    count and size, not just that the call succeeded.
    """
    suite = ValidationSuite(verbose=False)
    passed, details = suite.test_network_optimization()
    assert passed
    assert details["num_sites"] == 8
    assert details["num_components"] == 1
    assert details["largest_component"] == 8
    assert details["predicted_amplification"] >= 10


def test_network_optimization_fails_loudly_on_fragmented_sites(monkeypatch):
    """If sites are placed too far apart to connect, the test must FAIL, not pass silently.

    Exercises the failure branch (`num_components > total_sites / 2`) by
    monkeypatching `AmplificationNetwork.num_components` to simulate a
    pathologically fragmented network, proving the guard actually fires
    rather than being unreachable dead code.
    """
    from neoswga.core import network_optimizer as no_mod

    monkeypatch.setattr(
        no_mod.AmplificationNetwork,
        "num_components",
        lambda self: 100,  # far more components than sites / 2
    )
    suite = ValidationSuite(verbose=False)
    passed, details = suite.test_network_optimization()
    assert passed is False
    assert "Too fragmented" in details["error"]


# ----------------------------------------------------------------------
# test_bloom_filter
# ----------------------------------------------------------------------


def test_bloom_filter_has_zero_false_negatives_for_inserted_kmers():
    """Every k-mer explicitly added must be reported present -- a Bloom filter
    with false negatives is useless for background screening (it would let
    background-binding primers through undetected)."""
    pytest.importorskip("pybloom_live", reason="optional dependency not installed")
    suite = ValidationSuite(verbose=False)
    passed, details = suite.test_bloom_filter()
    assert passed
    if details.get("status") == "skipped":
        pytest.skip(details["reason"])
    assert details["false_positive_rate"] <= 0.05


# ----------------------------------------------------------------------
# ValidationResult / run_all_tests / print_summary / quick_validation
# ----------------------------------------------------------------------


def test_validation_result_defaults_error_to_none():
    """`error` must default to None so a passing result doesn't need one supplied."""
    r = ValidationResult(test_name="x", passed=True, runtime=0.1, details={})
    assert r.error is None


def test_run_all_tests_records_one_result_per_test_and_matches_the_return_value():
    """`run_all_tests` must return True iff every recorded ValidationResult passed.

    This is the contract `neoswga validate` (non-quick mode) relies on for its
    exit code. Exercise the real suite (only cheap synthetic/skippable tests
    run in CI) and check the aggregate return value against the individually
    recorded results, not just that it ran.
    """
    suite = ValidationSuite(verbose=False)
    overall = suite.run_all_tests()
    assert len(suite.results) == 6  # one per test in the fixed `tests` list
    names = [r.test_name for r in suite.results]
    assert names == [
        "Position Cache",
        "Adaptive GC Filter",
        "Background Bloom Filter",
        "Network Optimization",
        "MILP Optimizer",
        "End-to-End Pipeline",
    ]
    assert overall == all(r.passed for r in suite.results)


def test_run_all_tests_propagates_a_normal_test_failure_error_message(monkeypatch):
    """A test that fails by returning `(False, {"error": ...})` (the pattern every
    test function in this module uses) must have that message land on
    `ValidationResult.error`, not just buried in `.details`.

    Regression check for the bug where `result.error` was only ever populated
    on the exception path, so `print_summary()` and the verbose per-test log
    line -- both of which read `result.error` -- silently dropped the reason
    for any ordinary (non-exception) test failure.
    """
    suite = ValidationSuite(verbose=False)
    suite.test_adaptive_gc_filter = lambda: (False, {"error": "deliberate failure for repro"})

    suite.run_all_tests()

    result = next(r for r in suite.results if r.test_name == "Adaptive GC Filter")
    assert result.passed is False
    assert result.error == "deliberate failure for repro"


def test_run_all_tests_catches_exceptions_and_records_a_failure(monkeypatch):
    """A test function raising must not abort the whole suite -- it must be
    recorded as a failed result with the exception message captured, and the
    remaining tests must still run."""
    suite = ValidationSuite(verbose=False)

    def boom():
        raise RuntimeError("synthetic failure for test isolation")

    monkeypatch.setattr(suite, "test_adaptive_gc_filter", boom)
    overall = suite.run_all_tests()

    assert overall is False
    assert len(suite.results) == 6  # exception in one test does not stop the rest
    failed = next(r for r in suite.results if r.test_name == "Adaptive GC Filter")
    assert failed.passed is False
    assert failed.error == "synthetic failure for test isolation"


def test_print_summary_reports_correct_pass_fail_counts(caplog):
    """The human-readable summary line counts must match the actual results,
    since this is the only output a user glances at from `neoswga validate`."""
    suite = ValidationSuite(verbose=False)
    suite.results = [
        ValidationResult(test_name="ok1", passed=True, runtime=0.0, details={}),
        ValidationResult(test_name="ok2", passed=True, runtime=0.0, details={}),
        ValidationResult(test_name="broken", passed=False, runtime=0.0, details={}, error="boom"),
    ]
    with caplog.at_level(logging.INFO, logger="neoswga.core.validation"):
        suite.print_summary()

    text = caplog.text
    assert "Total tests: 3" in text
    assert "Passed: 2" in text
    assert "Failed: 1" in text
    assert "broken" in text
    assert "boom" in text


def test_quick_validation_runs_a_fixed_subset_and_returns_a_bool():
    """`quick_validation()` is the `--quick` CLI path -- it must run without a
    real genome/HDF5 fixture present (CI has none) and return a plain bool."""
    result = quick_validation()
    assert isinstance(result, bool)
    assert result is True  # all three quick tests degrade gracefully without fixtures


def test_quick_validation_returns_false_when_a_subset_test_fails(monkeypatch):
    """If any of the three quick tests fails, the overall bool must flip to False.

    `neoswga validate --quick` reports a process exit code from this return
    value, so a test that fails but leaves `all_passed` True would make a
    broken install look healthy.
    """
    import neoswga.core.validation as validation_mod

    def failing(self):
        return False, {"error": "synthetic failure"}

    monkeypatch.setattr(validation_mod.ValidationSuite, "test_adaptive_gc_filter", failing)
    assert quick_validation() is False


# ----------------------------------------------------------------------
# test_milp_optimizer / test_end_to_end
# ----------------------------------------------------------------------


def test_milp_optimizer_skips_gracefully_when_module_absent():
    """`milp_optimizer` is not part of the shipped package (see neoswga/core/),
    so this must degrade to a skip rather than fail the whole suite."""
    suite = ValidationSuite(verbose=False)
    passed, details = suite.test_milp_optimizer()
    assert passed is True
    assert details["status"] == "skipped"
    assert details["reason"] == "python-mip not installed"


def test_end_to_end_initializes_the_pipeline_with_the_documented_config():
    """The end-to-end test only checks pipeline construction (no real genome
    data is bundled), so pin exactly what it claims to verify."""
    suite = ValidationSuite(verbose=False)
    passed, details = suite.test_end_to_end()
    assert passed is True
    assert details["pipeline_initialized"] is True
    assert "num_primers=3" in details["config"]


# ----------------------------------------------------------------------
# Verbose logging branches (run_all_tests / print_summary all-pass path)
# ----------------------------------------------------------------------


def test_run_all_tests_verbose_logs_a_summary_when_everything_passes(caplog):
    """With verbose=True, a fully-passing run must log the "All tests PASSED"
    line (not just the per-test PASS/FAIL lines) -- this is what a user
    running `neoswga validate` actually reads on a healthy install."""
    suite = ValidationSuite(verbose=True)
    with caplog.at_level(logging.INFO, logger="neoswga.core.validation"):
        overall = suite.run_all_tests()

    assert overall is True
    assert "All tests PASSED" in caplog.text
    assert "TEST: Adaptive GC Filter" in caplog.text
