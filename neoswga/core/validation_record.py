"""What goes into the delivered pool's validation record, besides the validator.

`result_validation.validate_result` answers one question -- did the optimizer
misbehave -- and this module holds the other two things the record carries:
the panel assessment, and the configured limits the delivered panel missed.

Extracted from `unified_optimizer.py` on 2026-09-23, which had reached 1672 of
its 1650-line ceiling once both arrived. They are one subject and were two
inline blocks in a function that is not about either.

The three contributors stay separate on purpose, because they answer different
questions about the same panel and none subsumes another:

- `validate_result` -- duplicates, drift from the requested size, a starved
  target, blacklist re-injection. Defects in the SEARCH.
- `panel_assessment` -- every metric with its units and the reach it was
  computed at, an unavailable quantity kept apart from a measured zero, a
  non-finite value refused. What the panel IS.
- `limit_violation_issue` -- the limits the USER configured and the bounded
  repair could not satisfy. What the user ASKED FOR.

Only the third is a blocking finding today. The first is blocking per-issue by
level, and the second records `qualified` and gates nothing until the
disagreement between it and the validator has been measured.
"""

from __future__ import annotations

from typing import Any

__all__ = ["panel_assessment", "limit_violation_issue"]


def panel_assessment(design_request, result, optimizer) -> dict[str, Any]:
    """The one acceptance record for the delivered panel, as a JSON dict.

    Recorded, not consulted. `evaluate_panel` carries each metric WITH its
    units and the reach it was computed at, keeps an unavailable quantity apart
    from a measured zero, and refuses a non-finite required value. Until the
    resolved request reached `run_optimization` it had no production caller at
    all, because it takes a `DesignRequest` and there was none to give it.

    It adds no issue, so it cannot change what `export` refuses: the gate is
    `design_result.BLOCKING_VALIDATOR_CODES`, not the `ok` flag. Letting it
    decide would move verdicts in cases nobody has enumerated, and enumerating
    them is a separate increment.
    """
    from .panel_evaluation import evaluate_panel
    from .panel_refinement import objective_for_optimizer

    return evaluate_panel(
        design_request,
        list(result.primers),
        result.metrics,
        objective=objective_for_optimizer(optimizer) if optimizer is not None else None,
    ).as_dict()


def limit_violation_issue(acceptance) -> dict[str, Any] | None:
    """A blocking issue when the panel misses a limit this run configured.

    None when no limit was configured, when every limit was met, or when the
    bounded repair resolved the violation.

    This is the only channel by which a configured limit reaches `export`. The
    `AcceptanceReport` carrying these violations was previously built inside
    `if verbose:` and then discarded, and `apply_configured_limits` returned a
    bare `None` both for a panel that met every limit and for one whose repair
    failed -- its own docstring states the rule that makes those identical, that
    a repair which does not resolve the violation returns its input. So a panel
    missing a limit somebody set exported under "Primers ready for ordering!".

    Scoped to the configured limits deliberately. The panel assessment's own
    `violations` list is wider: it includes a panel shorter than requested, and
    blocking on that would refuse a normal, documented outcome, since
    `num_primers` is a request rather than a guarantee.
    """
    if acceptance is None or not acceptance.violations:
        return None
    return {
        "level": "error",
        "code": "panel_limit_not_met",
        "detail": (
            "the delivered panel does not meet a configured limit and the bounded "
            "repair did not resolve it: "
            + "; ".join(acceptance.violations)
            + ". Relax the limit, widen the candidate pool, or accept the panel "
            "with --allow-unqualified"
        ),
    }
