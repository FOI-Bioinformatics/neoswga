"""Guards for the class of bug where a default never fires.

Three of the four bugs found in this sweep are the same mistake in different
places: a module-level parameter that *exists* holding ``None``, read through
``getattr(parameter, name, default)``. The default is dead code -- getattr only
falls back when the attribute is absent -- so ``None`` propagates until
something tries arithmetic on it, often far from the cause.

The symptom is that anything reaching the filtering or thermodynamic code before
``get_params`` has run raises rather than falling back. That covers library use
of ``neoswga.core``, and any helper invoked outside the CLI pipeline.

These tests exercise the pre-``get_params`` state deliberately, which no other
test did -- every existing test either runs the pipeline or sets the globals by
hand, so the uninitialised path was never touched.
"""

import json
import math
import subprocess
import sys

import pytest

# ----------------------------------------------------------------------
# filter_extra before get_params
# ----------------------------------------------------------------------


def test_filter_extra_runs_before_get_params():
    """`reaction_temp` is None at import, so `getattr(..., 30.0)` returned None.

    ReactionConditions(temp=None) then raised TypeError from its range check.
    """
    code = "from neoswga.core import filter as F;" "print(F.filter_extra('ATCGATCGAT'))"
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=180)
    assert r.returncode == 0, f"filter_extra raised on a fresh import:\n{r.stderr[-1500:]}"
    assert r.stdout.strip() in {"True", "False"}


def test_every_global_filter_extra_reads_has_a_module_default():
    """`max_self_dimer_bp` was read as a bare attribute and did not exist.

    A getattr fallback would not have helped -- this was a plain AttributeError.
    Rather than sprinkle defaults at the call site, the module now defines the
    same values as PipelineParameters.
    """
    from neoswga.core import parameter

    for name in ("gc_min", "gc_max", "max_self_dimer_bp", "max_dimer_bp"):
        assert hasattr(parameter, name), f"parameter.{name} missing at module level"
        assert getattr(parameter, name) is not None, name


def test_module_defaults_match_the_dataclass():
    """The two must not drift; PipelineParameters is the documented contract.

    Checked in a subprocess because these globals are mutable and the pipeline
    (and other tests) legitimately overwrite them. The claim is about the
    *import-time* values, which is exactly the state this file is about.
    """
    names = ("gc_min", "gc_max", "gc_tolerance", "max_dimer_bp", "max_self_dimer_bp")
    code = (
        "import json;"
        "from neoswga.core import parameter as p;"
        "d = p.PipelineParameters();"
        f"names = {names!r};"
        "print(json.dumps({n: [getattr(p, n), getattr(d, n)] for n in names}))"
    )
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=180)
    assert r.returncode == 0, r.stderr[-1500:]

    for name, (module_value, dataclass_value) in json.loads(r.stdout).items():
        assert module_value == dataclass_value, (
            f"parameter.{name} is {module_value} at import but "
            f"PipelineParameters.{name} defaults to {dataclass_value}"
        )


def test_default_reaction_temp_covers_every_polymerase():
    from neoswga.core.parameter import default_reaction_temp
    from neoswga.core.registry import views

    for name in views.as_characteristics():
        temp = default_reaction_temp(name)
        assert isinstance(temp, float) and math.isfinite(temp), name
        assert 20.0 <= temp <= 70.0, f"{name}: implausible default {temp}"

    # Unknown polymerase falls back rather than raising.
    assert default_reaction_temp("not-a-polymerase") == 30.0
    assert default_reaction_temp(None) == 30.0


# ----------------------------------------------------------------------
# Def-time default binding
# ----------------------------------------------------------------------


def test_adaptive_gc_window_reads_tolerance_at_call_time():
    """A module global bound as a signature default freezes at import.

    `def f(x, tol=gc_tolerance)` captures 0.15 once, so later changes to
    `parameter.gc_tolerance` are silently ignored by every default-path caller.
    """
    from neoswga.core import parameter

    original = parameter.gc_tolerance
    try:
        assert parameter.adaptive_gc_window(0.50) == (0.35, 0.65)
        parameter.gc_tolerance = 0.30
        assert parameter.adaptive_gc_window(0.50) == pytest.approx((0.20, 0.80))
    finally:
        parameter.gc_tolerance = original


# ----------------------------------------------------------------------
# Non-finite floats must not reach JSON
# ----------------------------------------------------------------------


def _strict_loads(text):
    """json.loads that rejects Infinity/-Infinity/NaN, as RFC 8259 requires."""

    def reject(token):
        raise ValueError(f"non-standard JSON token {token!r}")

    return json.loads(text, parse_constant=reject)


def test_failed_optimization_serialises_to_valid_json():
    """`step4_improved_df_summary.json` is read by the report and by users.

    A failed run carried score=-inf and inf gaps, which json.dump writes as the
    bare tokens `-Infinity` / `Infinity`. Python reads those back, so the file
    looked fine from inside the project while being invalid for every other
    JSON parser.
    """
    from neoswga.core.base_optimizer import OptimizationResult

    result = OptimizationResult.failure("hybrid", "no feasible set")
    text = json.dumps(result.to_dict(), indent=2, default=str)

    parsed = _strict_loads(text)  # raises if a non-standard token is present
    assert parsed["score"] is None
    assert parsed["metrics"]["mean_gap"] is None
    assert parsed["metrics"]["max_gap"] is None


def test_json_safe_leaves_finite_values_alone():
    from neoswga.core.base_optimizer import json_safe

    for value in (0.0, -1.5, 3000.0, 1, "text", None, [1, 2]):
        assert json_safe(value) == value

    assert json_safe(float("inf")) is None
    assert json_safe(float("-inf")) is None
    assert json_safe(float("nan")) is None


def test_successful_optimization_keeps_its_numbers():
    """The sanitiser must not blank out real measurements."""
    from neoswga.core.base_optimizer import PrimerSetMetrics, json_safe

    metrics = PrimerSetMetrics.empty()
    finite = {
        k: json_safe(v)
        for k, v in metrics.to_dict().items()
        if isinstance(v, float) and math.isfinite(v)
    }
    assert all(v is not None for v in finite.values())


# ----------------------------------------------------------------------
# Downstream: the report must not render "inf kb"
# ----------------------------------------------------------------------


@pytest.mark.parametrize(
    "mean_gap,max_gap,gini,why",
    [
        (3000, float("inf"), 0.5, "inf straight from a failed optimizer"),
        (None, None, None, "round-tripped through JSON null"),
        (float("nan"), 9000, 0.5, "nan from an empty aggregation"),
    ],
)
def test_non_finite_gaps_report_unknown_not_a_number(mean_gap, max_gap, gini, why):
    """It rendered "About inf kb of this gap sits beyond the reach ..." verbatim."""
    from neoswga.core.coverage import interpret_gap_metrics

    rows = interpret_gap_metrics(mean_gap, max_gap, gini)
    assert all(r["verdict"] == "unknown" for r in rows), why
    for row in rows:
        assert "inf" not in row["value"].lower()
        assert "nan" not in row["value"].lower()
        assert "inf" not in row["detail"].lower()


def test_real_measurements_still_get_real_verdicts():
    """Regression guard on the above: do not blanket everything as unknown."""
    from neoswga.core.coverage import interpret_gap_metrics

    rows = interpret_gap_metrics(4100, 31500, 0.533)
    assert [r["verdict"] for r in rows] != ["unknown"] * 3
