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
    min_coverage: float = 0.0,
    min_per_target_coverage: float = 0.0,
    forbidden_primers: Optional[List[str]] = None,
) -> Dict[str, Any]:
    """Post-optimization sanity validation.

    Catches the broad class of silent failures an optimizer can produce
    without explicit assertions: duplicates, unexpected set size, zero
    coverage, blacklist re-injection. Returns a report dict that
    `unified_optimizer.run_optimization` stores on the result and writes
    alongside ``step4_improved_df.csv``.

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

    # Non-empty foreground coverage. ERROR status already carries 0.0
    # coverage implicitly, so skip that case.
    if result.status.value != "error":
        fg_cov = getattr(result.metrics, "fg_coverage", 0.0)
        if fg_cov < min_coverage:
            issues.append(
                {
                    "level": ("warning" if result.status.value == "partial" else "error"),
                    "code": "coverage_below_threshold",
                    "detail": (f"fg_coverage={fg_cov:.3f} < min_coverage={min_coverage:.3f}"),
                }
            )

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
