"""The post-optimization validator, and what makes a pool fit to recommend.

Extracted from `base_optimizer.py` on 2026-09-22 because that module sat at
exactly its 1600-line ceiling and three PRs each fitting it to exactly 1600
stacked to 1603 once merged. A ceiling with no headroom fails on the merge
rather than on the change, which is a worse place to find out.

The subject is its own: whether a delivered pool may be recommended is a
question about the pool, not about the optimizer that produced it, and
`design_result` already owns the rule that turns the findings into a verdict.

This never mutates the result. It emits `level="warning"` for
degenerate-but-acceptable outcomes and `level="error"` for genuine correctness
bugs, and its caller appends further findings to the SAME issue list before
writing it, which is why `ok` must be recomputed from the final list rather
than folded in here.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

from .design_result import panel_validation_is_ok

__all__ = ["validate_result"]


def validate_result(
    result,
    target_size: Optional[int] = None,
    min_per_target_coverage: float = 0.0,
    forbidden_primers: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """Post-optimization sanity validation.

    Catches the broad class of silent failures an optimizer can produce
    without explicit assertions: duplicates, unexpected set size, a starved
    target, blacklist re-injection. Returns a report dict that
    `unified_optimizer.run_optimization` writes alongside
    ``step4_improved_df.csv``.

    It does NOT check foreground coverage; the comment below says why, and
    "zero coverage" stood in this list for a check that could not fire.

    This never mutates the result. It emits warnings (``level='warning'``)
    for degenerate-but-acceptable outcomes (e.g., target_size underfilled
    when status is already PARTIAL) and errors (``level='error'``) for
    genuine correctness bugs (e.g., duplicate primers).
    """
    issues: List[Dict[str, Any]] = []

    primer_set = list(result.primers)
    n = len(primer_set)

    # Duplicate check — a hard error regardless of status.
    unique = set(primer_set)
    if len(unique) != n:
        dup_counts = {p: primer_set.count(p) for p in unique}
        duplicates = [p for p, c in dup_counts.items() if c > 1]
        issues.append(
            {
                "level": "error",
                "code": "duplicate_primers",
                "detail": f"{len(duplicates)} duplicate primer(s): {duplicates[:5]}",
            }
        )

    # Set size check — warning unless the optimizer already flagged PARTIAL.
    if target_size is not None and n != target_size:
        issues.append(
            {
                "level": ("warning" if result.status.value == "partial" else "error"),
                "code": "set_size_mismatch",
                "detail": f"Requested target_size={target_size}, got {n}",
            }
        )

    # A foreground-coverage floor used to live here and was retired on
    # 2026-09-22. It could not fire: the one production call site passed
    # `min_coverage=0.0` and no caller tightened it, and a coverage is never
    # below zero. The deeper reason it is not simply wired is that this
    # function records optimizer MISBEHAVIOUR -- duplicates, drift from the
    # requested size, blacklist re-injection -- and low coverage is a design
    # outcome rather than a bug in the search. Acceptance criteria belong to
    # `panel_acceptance.LIMIT_KEYS`, and none of those six is a minimum
    # coverage. `tests/test_the_coverage_floor_had_no_reachable_setting.py`
    # holds the decision and the evidence bar for filling that gap.

    # Per-target coverage (multi-genome mode). Optimizers that populate
    # `metrics.per_target_coverage` (Phase 11D) have their minimum
    # checked here. Missing attribute means single-genome mode or the
    # optimizer did not report per-target numbers yet.
    per_target = getattr(result.metrics, "per_target_coverage", None)
    if per_target and min_per_target_coverage > 0.0:
        below = {k: v for k, v in per_target.items() if v < min_per_target_coverage}
        if below:
            issues.append(
                {
                    "level": "warning",
                    "code": "per_target_coverage_below_threshold",
                    "detail": (
                        f"{len(below)} target(s) below {min_per_target_coverage:.2f}: {below}"
                    ),
                }
            )

    # Blacklist re-injection guard. Caller passes forbidden_primers when
    # a blacklist is configured; catches the case where expand-primers
    # or swap-primer fed a non-filtered candidate pool into the optimizer.
    if forbidden_primers:
        hits = set(primer_set) & set(forbidden_primers)
        if hits:
            issues.append(
                {
                    "level": "error",
                    "code": "blacklist_primer_in_set",
                    "detail": f"Primers in blacklist: {sorted(hits)[:5]}",
                }
            )

    # Callers append to `issues` after this returns and must recompute
    # `ok` from the final list; `panel_validation_is_ok` says why.
    from .design_result import panel_validation_is_ok

    return {
        "optimizer": result.optimizer_name,
        "num_primers": n,
        "ok": panel_validation_is_ok(issues),
        "issues": issues,
    }
