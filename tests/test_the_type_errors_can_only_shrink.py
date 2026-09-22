"""The type-error ratchet's rule, tested without measuring anything.

`scripts/check_type_ratchet.py` is the gate and CI's lint job runs it. This
file tests the rule it applies, against synthetic counts, so it needs no mypy
and no particular environment.

That split was forced by two findings, and both are worth keeping.

**The counts depend on what is installed.** mypy sees more when it can resolve
an import, so the same commit measures differently depending on which optional
dependencies are present. CI's six-cell test matrix reported `parameter.py` at
19 errors where the lint environment reports 14, and the first CI run of this
gate failed on that difference rather than on a defect. Adding pytest to an
otherwise mypy-only venv moved the totals from 164 to 200 by itself.

So the gate runs in ONE environment -- the lint job, which installs four pinned
linters and nothing else -- and the baseline is only meaningful there.

**`tests/conftest.py` imports numpy**, so pytest cannot start in that
environment at all. A gate that has to run where pytest cannot run should not
be a test.

What is left here is the part that is genuinely environment-free: given counts
and a baseline, which differences are a regression. That is the rule the gate
turns on, and it is worth pinning because the second of its three cases --
a file NOT in the baseline must be clean -- is the half a repository-wide total
cannot give you. New code is where a type error is cheapest to fix and
likeliest to be a real defect.
"""

import importlib.util
import json
import pathlib

import pytest

ROOT = pathlib.Path(__file__).resolve().parent.parent
SCRIPT = ROOT / "scripts" / "check_type_ratchet.py"
BASELINE = ROOT / "tests" / "type_error_baseline.json"


@pytest.fixture(scope="module")
def ratchet():
    """The gate's own module, loaded by path: `scripts/` is not a package."""
    spec = importlib.util.spec_from_file_location("check_type_ratchet", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# ---------------------------------------------------------------------------
# The rule
# ---------------------------------------------------------------------------


def test_matching_counts_hold(ratchet):
    assert ratchet.compare({"a.py": 3}, {"a.py": 3}) == []


def test_an_improvement_holds(ratchet):
    """A file may shrink freely. Forcing an edit on every fix would make
    lowering the number a chore and the baseline a fiction."""
    assert ratchet.compare({"a.py": 1}, {"a.py": 3}) == []


def test_a_regression_fails(ratchet):
    kinds = [kind for kind, _ in ratchet.compare({"a.py": 4}, {"a.py": 3})]

    assert kinds == ["worse"]


def test_a_file_outside_the_baseline_must_be_clean(ratchet):
    """The half a repository-wide total cannot give you."""
    kinds = [kind for kind, _ in ratchet.compare({"new.py": 1}, {"a.py": 3})]

    assert "unlisted" in kinds


def test_a_clean_file_outside_the_baseline_is_fine(ratchet):
    """Only files WITH errors are measured, so a clean new module is silent."""
    assert ratchet.compare({}, {}) == []


def test_a_file_that_became_clean_must_leave_the_baseline(ratchet):
    """A fix is not finished until it also removes its excuse."""
    kinds = [kind for kind, _ in ratchet.compare({}, {"a.py": 3})]

    assert kinds == ["stale"]


def test_the_message_names_the_file_and_both_numbers(ratchet):
    """A ratchet failure a contributor cannot act on gets suppressed."""
    _kind, detail = ratchet.compare({"core/thing.py": 9}, {"core/thing.py": 4})[0]

    assert "core/thing.py" in detail and "9" in detail and "4" in detail


# ---------------------------------------------------------------------------
# The gate is wired, and to the environment it was measured in
# ---------------------------------------------------------------------------


def test_the_baseline_parses_and_is_not_growing_quietly(ratchet):
    baseline = json.loads(BASELINE.read_text())

    assert baseline, "the baseline is empty, so the gate compares against nothing"
    assert sum(baseline.values()) <= 200, (
        f"the baseline now admits {sum(baseline.values())} type errors. It is a "
        "ratchet; raising the total is a decision, not a fix."
    )


def test_ci_runs_the_gate_in_the_environment_it_was_measured_in():
    """Every part of this gate is outside pytest, so nothing else would notice
    if CI stopped running it: the suite would stay green and the ratchet would
    be gone."""
    workflow = (ROOT / ".github" / "workflows" / "ci.yml").read_text()
    module = SCRIPT.relative_to(ROOT).as_posix()

    assert module in workflow, f"no CI step runs {module}; the ratchet is not running"

    lint_job = workflow.split("  test:")[0]
    assert module in lint_job, (
        "the ratchet must run in the lint job. Its counts depend on which "
        "optional dependencies are importable, and the test matrix installs a "
        "different set -- which is what failed its first CI run."
    )


def test_ci_pins_mypy_to_the_version_the_baseline_was_recorded_under(ratchet):
    """An upgrade moves the counts with nothing in this repository changing --
    the drift that turned every model-loading test red when skops 0.15 narrowed
    its default trust list."""
    workflow = (ROOT / ".github" / "workflows" / "ci.yml").read_text()

    assert (
        f'"mypy=={ratchet.BASELINE_MYPY_VERSION}"' in workflow
    ), f"CI does not pin mypy=={ratchet.BASELINE_MYPY_VERSION}"


def test_the_gate_measures_what_ci_prints(ratchet):
    """The informational mypy step and the gate must see the same thing."""
    workflow = (ROOT / ".github" / "workflows" / "ci.yml").read_text()

    for argument in ratchet.MYPY_ARGS:
        if argument == "--no-color-output":
            continue  # only for parsing; the printed step renders for a human
        assert argument in workflow, f"CI's mypy step does not pass {argument}"
