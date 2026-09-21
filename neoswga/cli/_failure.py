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

from neoswga.core.design_result import RunState, describe_failure

logger = logging.getLogger(__name__)

__all__ = ["RunState", "write_failure_artifact", "report_design_failure"]


def _failure_artifact_path(args):
    """Where a failure record is written, or None when there is no run directory.

    A command with no output directory still fails loudly on stderr; it just
    has nowhere to leave a record. Writing one into the working directory
    instead would litter, and writing none is not a silent success because the
    exit code is nonzero either way.
    """
    directory = getattr(args, "data_dir", None) or getattr(args, "output", None)
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
