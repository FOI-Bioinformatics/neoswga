"""Importing the CLI must not drag in the pipeline or scikit-learn.

`cli_unified` imported `neoswga.core.pipeline` at module scope for the single
purpose of making `StepPrerequisiteError` catchable. `core/pipeline.py` imports
`rf_preprocessing`, which imports `sklearn` and `sklearn.ensemble` at module
scope. Measured with `-X importtime` on 2026-09-06: 0.60 s of a 1.14 s warm
import was sklearn, for a model the default path no longer loads.

These run in a subprocess because `sys.modules` in the test process is already
polluted by every other test in the suite.
"""

import json
import subprocess
import sys


def _modules_after(import_stmt):
    """Names in sys.modules after `import_stmt`, from a clean interpreter."""
    code = (
        f"{import_stmt}\n"
        "import sys, json\n"
        "print(json.dumps(sorted(m for m in sys.modules if '.' not in m)))\n"
    )
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True)
    return set(json.loads(out.stdout.strip().splitlines()[-1]))


def test_importing_the_cli_does_not_import_sklearn():
    loaded = _modules_after("import neoswga.cli_unified")
    assert "sklearn" not in loaded, (
        "importing the CLI pulled in scikit-learn; the amplification model was "
        "retired from the default path on 2026-09-05 and costs about 0.6 s here"
    )


def test_importing_the_cli_does_not_import_the_pipeline():
    code = (
        "import neoswga.cli_unified, sys, json\n"
        "print(json.dumps('neoswga.core.pipeline' in sys.modules))\n"
    )
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True)
    assert (
        out.stdout.strip().splitlines()[-1] == "false"
    ), "cli_unified still imports neoswga.core.pipeline at module scope"


def test_the_exception_lives_in_exceptions_and_is_still_reachable_from_pipeline():
    """Two tests import it from `neoswga.core.pipeline`; that must keep working,
    and it must be the same object so `except` clauses still match."""
    from neoswga.core.exceptions import StepPrerequisiteError as FromExceptions
    from neoswga.core.exceptions import StepValidationResult as ResultFromExceptions
    from neoswga.core.pipeline import StepPrerequisiteError as FromPipeline
    from neoswga.core.pipeline import StepValidationResult as ResultFromPipeline

    assert FromExceptions is FromPipeline
    assert ResultFromExceptions is ResultFromPipeline


def test_the_exception_still_formats_its_remediation():
    from neoswga.core.exceptions import StepPrerequisiteError, StepValidationResult

    validation = StepValidationResult(
        valid=False,
        missing_files=["/tmp/a.fasta", "/tmp/b.fasta"],
        error_message="Genome files not found.",
        remediation="Run: neoswga count-kmers -j params.json",
    )
    err = StepPrerequisiteError(step=2, validation=validation)

    assert err.step == 2
    assert err.validation is validation
    text = str(err)
    assert "STEP 2 PREREQUISITE ERROR" in text
    assert "/tmp/a.fasta" in text
    assert "neoswga count-kmers" in text


def test_cli_pipeline_has_no_module_scope_fallback_class():
    """The try/except ImportError fallback defined a *different* class, so a
    raised StepPrerequisiteError would not have been caught by it."""
    from pathlib import Path

    import neoswga.cli.pipeline as cli_pipeline

    source = Path(cli_pipeline.__file__).read_text()
    assert "class StepPrerequisiteError(Exception):" not in source, (
        "cli/pipeline.py still defines a fallback exception class that would "
        "not catch the one core/pipeline.py raises"
    )


def test_the_exception_raised_by_the_pipeline_is_caught_by_the_cli_import():
    """`unified_optimizer` imports the name from `core.pipeline` (the
    re-export) while `cli/pipeline.py` catches the name imported from
    `core.exceptions`. `except` matches on class identity, so the re-export
    must not be a distinct class."""
    from neoswga.core.exceptions import StepPrerequisiteError as Caught
    from neoswga.core.pipeline import StepPrerequisiteError as Raised
    from neoswga.core.pipeline import StepValidationResult

    validation = StepValidationResult(
        valid=False,
        missing_files=[],
        error_message="missing",
        remediation="re-run",
    )
    try:
        raise Raised(step=4, validation=validation)
    except Caught as exc:
        assert exc.step == 4
    else:  # pragma: no cover - the assert above is the point of the test
        raise AssertionError("the re-exported class was not caught")
