"""Ranking the members of an ensemble run, and recording the ones that failed.

Split out of `unified_optimizer.py` on 2026-09-21. Two small units that answer
one question -- how the comparison table is built and which row wins -- and
that the valid-design contract keeps separate from running the search itself:
a proposal method reports what it produced, and ranking those proposals is not
the method's decision to make.

Selection is by the shared objective when one is available, because raw
`score` is not comparable across optimizers. The tie-breaks are load-bearing
rather than cosmetic: ties are common here, since `background-aware` wraps the
same `HybridOptimizer` that `hybrid` uses and the two frequently return the
identical set. Selection used to be a bare `max()` over a dict built in
`--ensemble-methods` order, so reordering three tied methods returned three
different winners while the documentation claimed order-independence.
"""

__all__ = ["_select_ensemble_winner", "_ensemble_error_row"]


def _select_ensemble_winner(results, application, objective=None):
    """Rank on shared constraints and objective, with stable size/name ties.

    Result-only library optimizers without a panel evaluator retain the
    application-weighted normalized score as their comparison rule.
    """

    def rank(method):
        result = results[method]
        if objective is not None:
            from .optimization_service import panel_rank

            return (
                not bool(
                    result.stage_history and result.stage_history[-1].get("failed_constraints")
                ),
                *panel_rank(objective, result.primers),
                -len(result.primers),
                tuple(-ord(c) for c in method),
            )
        return (
            result.metrics.normalized_score(application=application),
            -len(result.primers),
            # `max` takes the largest, so invert the name to prefer the earliest.
            tuple(-ord(c) for c in method),
        )

    return max(results, key=rank)


def _ensemble_error_row(method: str) -> dict:
    """The comparison row for a member that could not run.

    `-inf` rather than 0.0 for the raw score: a method that produced nothing
    must never out-rank one that produced a poor panel, and an empty
    `PrimerSetMetrics` does not score zero. `selected` is False so the row is
    visible in `ensemble_comparison` without being eligible to win.
    """
    return {
        "method": method,
        "normalized_score": 0.0,
        "score": float("-inf"),
        "n_primers": 0,
        "fg_coverage": 0.0,
        "bg_coverage": 0.0,
        "status": "error",
        "selected": False,
    }
