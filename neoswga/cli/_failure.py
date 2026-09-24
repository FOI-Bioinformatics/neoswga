"""Structured failure records at the command boundary.

Split out of `cli_unified.py` so that module stays inside its size budget, and
because this is one concern: what a command writes and prints when a design
run does not finish.

The rule it enforces is that an output directory must never look current after
a run that failed. A directory holding last week's `step4_improved_df.csv` and
nothing else is indistinguishable from one holding this morning's, so a failed
run leaves `design_failure.json` beside it saying otherwise.

Kept to the standard library plus `core.design_result`, which is itself
dependency-light: this is imported on every invocation including `--help`, and
Known Issue 12 records what a heavier import chain costs there.
"""

import json
import logging
import os
import sys

from neoswga.core.design_result import RunState, describe_failure

logger = logging.getLogger(__name__)

__all__ = [
    "RunState",
    "exit_on_step_failure",
    "report_design_failure",
    "write_failure_artifact",
]


def _data_dir_from_params_file(json_file):
    """The output directory a params file names, resolved as a run resolves it.

    Relative to the file's own directory, which is how every command reads
    paths out of params.json.
    """
    if not json_file:
        return None
    try:
        with open(json_file) as handle:
            params = json.load(handle)
    except (OSError, ValueError):
        return None
    directory = params.get("data_dir")
    if not directory:
        return None
    if os.path.isabs(directory):
        return directory
    return os.path.join(os.path.dirname(os.path.abspath(json_file)), directory)


def _failure_artifact_path(args):
    """Where a failure record is written, or None when there is no run directory.

    A command with no output directory still fails loudly on stderr; it just
    has nowhere to leave a record. Writing one into the working directory
    instead would litter, and writing none is not a silent success because the
    exit code is nonzero either way.
    """
    directory = getattr(args, "data_dir", None)
    if not directory:
        # The pipeline commands resolve their output directory into the
        # `parameter` module during `get_params`, and several of them carry no
        # `--data-dir` at all. Reading it back here is what puts the record
        # beside the run rather than nowhere.
        try:
            from neoswga.core import parameter

            directory = getattr(parameter, "data_dir", None)
        except ImportError:  # pragma: no cover - the package is always importable
            directory = None
    if not directory:
        # And when the refusal came BEFORE `get_params` ran, that module still
        # holds its default and knows nothing about this run. Every validation
        # added by the design-request contract fires there, which is exactly
        # when a record is most wanted: `optimize -j params.json` with an
        # unusable setting exited nonzero and left nothing behind, so the
        # directory still looked like its previous successful run.
        #
        # The same ordering trap as `warn_on_condition_drift`, reached from
        # the other side.
        directory = _data_dir_from_params_file(getattr(args, "json_file", None))
    if not directory:
        # `--output` last, and only when it is ALREADY a directory. Its meaning
        # varies by command: for `analyze-coverage` it names a directory, for
        # `calibrate-reach`, `predict` and `report` it names a FILE. Creating a
        # directory at a file path is destructive twice over -- the next
        # successful run then cannot write its own output there, and the record
        # lands where `export.export_is_blocked` never looks, so it blocks
        # nothing. Measured: `_failure_artifact_path(Namespace(output="fit.json"))`
        # left a DIRECTORY called `fit.json` behind.
        candidate = getattr(args, "output", None)
        directory = candidate if candidate and os.path.isdir(candidate) else None
    if not directory:
        return None
    try:
        os.makedirs(directory, exist_ok=True)
    except OSError:
        return None
    return os.path.join(directory, "design_failure.json")


def write_failure_artifact(args, error, run_state=RunState.FAILED):
    """Record the failure beside the run so a stale directory cannot mislead.

    An output directory holding last week's `step4_improved_df.csv` and nothing
    else looks exactly like a directory holding a current result. This file is
    what tells a later command, or a person, that the most recent run did not
    finish. Failing to WRITE it must not mask the original error, so every
    error here is swallowed after being logged.
    """
    path = _failure_artifact_path(args)
    if path is None:
        return None
    if error is None:
        record = {
            "run_state": run_state,
            "termination": "interrupted",
            "stage": getattr(args, "command", "unknown"),
            "recommendation_written": False,
            "qualified": False,
        }
    else:
        record = describe_failure(
            error,
            stage=getattr(args, "command", "unknown"),
            request_hash=getattr(args, "request_hash", None),
        )
        record["run_state"] = run_state
    try:
        with open(path, "w") as handle:
            json.dump(record, handle, indent=2, default=str)
    except OSError as exc:
        logger.debug("Could not write the failure record to %s: %s", path, exc)
        return None
    return path


def report_design_failure(args, error, expected=True):
    """Print the named stage, input and model, then leave a machine record."""
    if expected:
        logger.error("Design failed (%s): %s", type(error).__name__, error)
        for label, attribute in (
            ("field", "field"),
            ("artifact", "artifact"),
            ("model", "model"),
            ("quantity", "quantity"),
            ("input", "subject"),
        ):
            value = getattr(error, attribute, None)
            if value is not None:
                logger.error("  %s: %s", label, value)
        remediation = getattr(error, "remediation", None)
        if remediation:
            logger.error("  to fix: %s", remediation)
        logger.error(
            "  No oligo pool was recommended. A failed calculation is not a "
            "screened candidate, so nothing here was substituted for it."
        )
    path = write_failure_artifact(args, error)
    if path:
        logger.error("  failure record: %s", path)


def clear_failure_artifact(directory):
    """Remove a failure record once a run in this directory has succeeded.

    A record that is never cleared is as wrong as one that is never read: the
    user fixes the problem, the next run succeeds, and the export still
    refuses on evidence that no longer describes anything.

    Deliberately quiet about a missing file. Most runs have nothing to clear,
    and the absence is the normal case rather than a condition worth reporting.
    """
    from neoswga.core.export import DESIGN_FAILURE_FILENAME

    if not directory:
        return
    path = os.path.join(str(directory), DESIGN_FAILURE_FILENAME)
    try:
        os.remove(path)
    except FileNotFoundError:
        return
    except OSError as exc:  # pragma: no cover - unwritable directory
        logger.debug("Could not clear the failure record at %s: %s", path, exc)


def exit_on_step_failure(step_name, error, logger, data_dir=None):
    """Report a step's failure the one way, then exit nonzero.

    Steps 2, 3 and 4 each had this block inline and they had drifted: only
    one of them explained an HDF5 lock collision, which is the failure a user
    is most likely to hit and least likely to diagnose. Keeping it in one
    place is what stops the next addition landing in one step and not the
    others.

    `DesignError` does not come here. It carries its own record and is
    re-raised to the command boundary; see `report_design_failure`.
    """
    import logging as _logging

    from neoswga.core.concurrent_runs import locked_file_advice

    logger.error(f"{step_name} failed: {error}")
    advice = locked_file_advice(error, data_dir)
    if advice:
        logger.error(advice)
    if logger.level <= _logging.DEBUG:
        import traceback

        traceback.print_exc()
    else:
        # Step 4 alone printed this hint. Another drift the three inline
        # copies had accumulated, and the reason they are now one.
        logger.error("Run with --verbose for full traceback")
    sys.exit(1)


def report_import_failure(error):
    """Say which import failed before advising a reinstall.

    This handler wraps a whole pipeline step, so it catches an `ImportError`
    raised anywhere inside it -- not only a failure to import the pipeline
    module. The common case is an OPTIONAL dependency: networkx for the clique
    optimizer, pysam for the [bam] extra, pybloom_live, matplotlib. Each of
    those raises with a message naming exactly what to install, and this then
    added "This may indicate a corrupted installation. Try: pip install -e .
    --force-reinstall" underneath it -- correct advice followed by wrong
    advice, with the wrong one last.

    Reinstalling the package does not add an optional extra, so following it
    costs a reinstall and leaves the failure in place.

    `ImportError.name` is the module that could not be imported, so a
    third-party one is reported as what it is and the original message, which
    already names the remedy, is left to stand.
    """
    module = getattr(error, "name", None) or ""
    logger.error("Import failed: %s", error)
    if module and not module.startswith("neoswga"):
        logger.error(
            "%r is an optional dependency rather than part of neoswga. "
            "Install it, or the extra that provides it; reinstalling the "
            "package will not add it.",
            module.split(".")[0],
        )
        return
    logger.error("This may indicate a corrupted installation.")
    logger.error("Try: pip install -e . --force-reinstall")
