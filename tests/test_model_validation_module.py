"""Coverage for neoswga.core.model_validation (literature-behaviour checks for
the mechanistic model). Exercised by `neoswga validate-model`."""

from neoswga.core import model_validation as mv


def test_validate_mechanistic_model_returns_result_list():
    results = mv.validate_mechanistic_model()
    assert isinstance(results, list) and results
    for r in results:
        assert "test" in r
        assert "passed" in r
        assert isinstance(r["passed"], bool)


def test_individual_validators_return_expected_shape():
    for fn in (
        mv.validate_dmso_inhibition,
        mv.validate_betaine_biphasic,
        mv.validate_gc_accessibility,
        mv.validate_temperature_optimum,
        mv.validate_binding_kinetics,
        mv.validate_mg_effects,
        mv.validate_formamide_inhibition,
    ):
        result = fn()
        assert set(("test", "passed", "summary")).issubset(result)
        assert isinstance(result["passed"], bool)


def test_format_validation_report_is_text():
    results = mv.validate_mechanistic_model()
    report = mv.format_validation_report(results)
    assert isinstance(report, str)
    assert "PASS" in report or "FAIL" in report
    # every test title appears in the rendered report
    for r in results:
        assert r["test"] in report


def test_quick_validation_returns_bool_and_summary():
    ok, summary = mv.quick_validation()
    assert isinstance(ok, bool)
    assert isinstance(summary, str) and summary
