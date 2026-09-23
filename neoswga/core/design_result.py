"""Run status, termination reason, panel qualification and output eligibility.

Three properties of a design run that used to be one. Collapsing them is how a
run that crashed during final validation, and a run that finished having found
nothing acceptable, both ended up writing an oligo pool.

- **Run state** is what happened to the process: it finished, it failed, or it
  was interrupted. Nothing about the panel.
- **Termination reason** is why the search stopped: it qualified, it spent its
  allowance, it ran out of candidates, or it hit an error. A spent allowance is
  a recorded stopping point, not a failure and not a proof of infeasibility.
- **Qualification** is a property of the delivered panel: does it satisfy every
  hard constraint under the request's own chemistry.

Only a finished run with a qualifying panel may write a recommended pool. A
finished run that stopped on its budget may recommend its verified incumbent,
with the limit stated and no completeness claim; that is still `FINISHED` and
still qualified, so the gate here admits it. A failed run writes a failure
artifact and may keep a separately marked diagnostic checkpoint, which is not a
recommendation.

This module holds no dependencies beyond the standard library, for the same
reason `exceptions.py` does not: the CLI boundary imports it on every command.
"""

from __future__ import annotations

import traceback
from typing import Any, Dict, Optional

from .exceptions import DesignError

__all__ = [
    "RunState",
    "TerminationReason",
    "recommendation_allowed",
    "describe_failure",
]


class RunState:
    """What happened to the run. Not what happened to the panel."""

    #: The run completed its declared work, whatever it found.
    FINISHED = "finished"
    #: The run stopped on an error. A panel it held is diagnostic only.
    FAILED = "failed"
    #: The run was stopped from outside (signal, cancellation).
    INTERRUPTED = "interrupted"

    ALL = frozenset({FINISHED, FAILED, INTERRUPTED})


class TerminationReason:
    """Why the search stopped. Recorded on finished and failed runs alike."""

    #: A qualifying panel was found and the size policy was satisfied.
    QUALIFIED = "qualified"
    #: The shared allowance was spent. The incumbent is valid; the search is
    #: incomplete. Not a model failure and not a proof that no panel exists.
    BUDGET_EXHAUSTED = "budget_exhausted"
    #: Every eligible candidate was examined. This exhausts candidates, not
    #: candidate subsets, so it is still not an infeasibility certificate.
    CANDIDATES_EXHAUSTED = "candidates_exhausted"
    #: The frontier could not be widened further within its refill allowance.
    REFILL_EXHAUSTED = "refill_exhausted"
    #: A required calculation or reference answer failed. Always with FAILED.
    ERROR = "error"

    ALL = frozenset({QUALIFIED, BUDGET_EXHAUSTED, CANDIDATES_EXHAUSTED, REFILL_EXHAUSTED, ERROR})


def recommendation_allowed(run_state: str, qualified: bool) -> bool:
    """May this run write a recommended oligo pool?

    Both conditions, because each without the other has produced an invalid
    recommendation: an unqualified panel from a finished run is a panel that
    broke a constraint the user set, and a qualifying panel from a failed run
    was qualified against measurements the failure says we do not have.

    An unknown state raises rather than defaulting either way. A typo that
    returns False silently suppresses a legitimate recommendation, and one that
    returns True is the failure this function exists to prevent.
    """
    if run_state not in RunState.ALL:
        raise ValueError(f"Unknown run state {run_state!r}; expected one of {sorted(RunState.ALL)}")
    return run_state == RunState.FINISHED and bool(qualified)


def describe_failure(
    error: BaseException,
    stage: str,
    request_hash: Optional[str] = None,
    include_traceback: bool = True,
) -> Dict[str, Any]:
    """A JSON-serializable failure record naming stage, input and model.

    Written in place of a result so that a stale output directory cannot leave
    an apparently current recommendation beside a run that failed. Every value
    is a string, a bool or None, so the record round-trips through JSON
    unchanged; a caller that needs the live exception still has it.

    `include_traceback` is on by default because an unexpected exception must
    retain one. A `DesignError` is expected and self-describing, so its
    traceback is informative rather than required.
    """
    record: Dict[str, Any] = {
        "run_state": RunState.FAILED,
        "termination": TerminationReason.ERROR,
        "stage": str(stage),
        "request_hash": str(request_hash) if request_hash is not None else None,
        "error_type": type(error).__name__,
        "message": str(error),
        "expected": isinstance(error, DesignError),
        "recommendation_written": False,
        "qualified": False,
    }

    # Identifiers the four DesignError subclasses carry. Absent on an
    # unexpected exception, in which case the field stays None rather than
    # being filled with a guess.
    for field in ("field", "artifact", "model", "quantity"):
        value = getattr(error, field, None)
        record[field] = str(value) if value is not None else None
    for field in ("subject", "requested", "value"):
        value = getattr(error, field, None)
        if value is not None and record.get("input") is None:
            record["input"] = str(value)
    record.setdefault("input", None)
    remediation = getattr(error, "remediation", None)
    record["remediation"] = str(remediation) if remediation else None

    # A QC rejection and a failed calculation mean opposite things. This family
    # is never the former, and the record says so explicitly rather than
    # leaving a consumer to infer it.
    record["qc_reason"] = None

    cause = error.__cause__ or error.__context__
    record["cause_type"] = type(cause).__name__ if cause is not None else None
    record["cause_message"] = str(cause) if cause is not None else None

    if include_traceback:
        record["traceback"] = "".join(
            traceback.format_exception(type(error), error, error.__traceback__)
        )
    else:
        record["traceback"] = None
    return record


#: Validator codes that make a delivered pool unfit to recommend.
#:
#: Only findings that are defects IN THE POOL belong here. A pool breaking the
#: dimer threshold the user configured is one.
#: `coverage_saturated_on_small_genome` deliberately is NOT: it says a metric
#: cannot be trusted on a small target, which is inherent to designing against a
#: plasmid and not something the user can fix, so blocking on it would refuse
#: every such design and teach people to ignore the line.
#:
#: `panel_limit_not_met` is the second member, added 2026-09-23. A limit the
#: user configured is a defect in the pool by the user's own definition, which
#: is the test this comment already states for membership -- and unlike the
#: saturation code it is something they can act on, by relaxing the limit,
#: widening the pool or passing `--allow-unqualified`. It is raised only when
#: the bounded repair was attempted and did not resolve the violation. Setting
#: no limit cannot produce it: `constraints_from_parameter` returns None and no
#: objective is built.
#:
#: `duplicate_primers` and `blacklist_primer_in_set` joined on 2026-09-23. Both
#: were already recorded at `level="error"` and neither blocked, so a panel
#: holding the same oligo twice, or holding an oligo the user blacklisted,
#: printed "Primers ready for ordering!". `ok` was False in both cases, which is
#: why it survived: that flag says what a reader expects and no command consults
#: it. Both are defects in the pool and both are actionable -- a panel of twelve
#: holding one oligo twice is eleven distinct sequences in twelve tubes, and a
#: blacklisted oligo is the user's own instruction violated.
#:
#: They differ in reach, and the difference is recorded rather than smoothed
#: over. `blacklist_primer_in_set` fires on the real path, through
#: `_collect_forbidden_primers`, whenever `bl_prefixes` is configured.
#: `duplicate_primers` is not reachable from any path checked -- every point
#: that assembles a panel from two sources deduplicates first -- so it is
#: insurance, at the cost of one string here.
#:
#: `set_size_mismatch` stays OUT. A panel shorter than requested is a documented
#: normal outcome, and blocking on it would refuse three of four runs on the
#: shipped plasmid example.
#:
#: This set lived as a bare literal in BOTH `cli/report.py` and
#: `core/results_interpreter.py`, the two commands that tell a user a pool is
#: ready. A new code had to be added twice or they would disagree about the same
#: pool. It lives here because whether a result may be recommended is this
#: module's subject.
BLOCKING_VALIDATOR_CODES = frozenset(
    {
        "delivered_pool_exceeds_max_dimer_bp",
        "panel_limit_not_met",
        "duplicate_primers",
        "blacklist_primer_in_set",
    }
)

#: What step 4 writes about the pool it just delivered.
VALIDATION_FILENAME = "step4_improved_df_validation.json"


def blocking_validator_findings(results_dir) -> list:
    """Details of any recorded finding that makes this pool unfit to order.

    Returns an empty list when the file is ABSENT. That is deliberate and is
    the one place this module fails open: directories written before the
    validator existed are explicitly supported, and refusing them would be a
    refusal on the absence of evidence rather than on evidence.

    A file that exists and cannot be READ blocks, which is the same rule
    `export_is_blocked` applies to a corrupt failure record: there the evidence
    exists and cannot be read, and the reason nobody can tell whether it named
    a blocking finding is exactly why it must not be assumed it did not. This
    docstring stated that rule and the code applied it to the other artifact
    only, so a record truncated by a full disk or a killed process read as a
    clean pool.
    """
    import json
    import os

    path = os.path.join(str(results_dir), VALIDATION_FILENAME)
    if not os.path.exists(path):
        return []
    try:
        with open(path) as handle:
            payload = json.load(handle)
    except (OSError, ValueError) as exc:
        return [
            f"the recorded verdict at {path} exists and could not be read ({exc}). "
            "Whether it named a blocking finding is unknown, so this pool is not "
            "vouched for. Re-run `neoswga optimize` to rewrite it."
        ]
    if not isinstance(payload, dict):
        return [
            f"the recorded verdict at {path} is not a validation record. "
            "Re-run `neoswga optimize` to rewrite it."
        ]

    return [
        str(issue.get("detail") or issue.get("code"))
        for issue in payload.get("issues", []) or []
        if issue.get("code") in BLOCKING_VALIDATOR_CODES
    ]


def panel_validation_is_ok(issues) -> bool:
    """Whether a validator issue list leaves the pool fit to recommend.

    Two ways `ok` came apart from the findings beside it, both of which made a
    saved file say the pool was fine while `export` refused to order it.

    **`ok` was computed before the issues were complete.**
    `base_optimizer.validate` folded `level == "error"` into `ok` and returned,
    and `unified_optimizer` then appended the saturation warnings and the
    delivered-panel dimer finding to the SAME dict. A finding appended after
    the fold could not move the flag whatever its level, so the dimer breach
    was invisible to `ok` by construction rather than by policy.

    **A blocking code could be emitted at warning level.** It was:
    `dimer_validation_issue` recorded `level="warning"` while
    `BLOCKING_VALIDATOR_CODES` held its code, so `export` and `interpret`
    refused a pool the report rendered as having no errors.

    Checking both conditions means the two cannot drift apart again: a code
    that blocks makes `ok` false whatever level it carries, and a level that
    says error makes `ok` false whatever its code. Callers must apply this to
    the FINAL issue list, which is what `unified_optimizer` now does.

    Saturation warnings still leave `ok` true, deliberately. They say a metric
    cannot be trusted on a small target, which is inherent to designing against
    a plasmid rather than something the user can repair.
    """
    return not any(
        issue.get("level") == "error" or issue.get("code") in BLOCKING_VALIDATOR_CODES
        for issue in issues or []
        if isinstance(issue, dict)
    )
