"""`--minimize-primers` must not destroy an ensemble run.

`run_optimization` passes `optimizer=optimizer` to the trimming post-process,
and `optimizer` is bound only inside the single-method branch. The ensemble
branch never binds it, so the name was unbound and the step raised:

    UnboundLocalError: cannot access local variable 'optimizer'
    where it is not associated with a value

The exception is raised while the call arguments are being evaluated, so the
`try` inside `_minimize_primer_count` never sees it. It propagated to the CLI
handler, which logged "Step 4 failed" and exited 1.

The cost is the whole run. Every ensemble method has already finished and been
compared by the time this fires, and all of it is discarded: no
step4_improved_df.csv, no summary, no validation report. The `auto` and `all`
aliases reach the same branch.
"""

import pytest

from dataclasses import replace

from neoswga.core.base_optimizer import OptimizationResult, PrimerSetMetrics


@pytest.mark.parametrize("method", ["ensemble", "auto", "all"])
def test_minimize_primers_survives_an_ensemble_run(monkeypatch, method):
    from neoswga.core import unified_optimizer as uo

    primers = ("ACCACAGATAGC", "GTTGTAGATGGA", "ATCAGCAGACCA")
    winner = OptimizationResult(
        primers=primers,
        score=1.0,
        metrics=replace(PrimerSetMetrics.empty(), fg_coverage=0.9),
        status=uo.OptimizationStatus.SUCCESS,
        optimizer_name="hybrid",
        iterations=1,
    )

    class _Cache:
        def __init__(self, *a, **k):
            pass

        def get_positions(self, *a, **k):
            import numpy as np

            return np.array([], dtype=np.int64)

    monkeypatch.setattr(uo, "PositionCache", _Cache)
    monkeypatch.setattr(uo, "StreamingPositionCache", _Cache, raising=False)
    monkeypatch.setattr(uo, "_run_ensemble", lambda *a, **k: winner)
    monkeypatch.setattr(uo, "save_results", lambda *a, **k: None)

    trimmed = []
    monkeypatch.setattr(
        uo,
        "_minimize_primer_count",
        lambda **kw: trimmed.append(kw) or kw["result"],
    )

    result = uo.run_optimization(
        method=method,
        candidates=list(primers),
        fg_prefixes=["fg"],
        fg_seq_lengths=[100_000],
        bg_prefixes=[],
        bg_seq_lengths=[],
        target_size=3,
        verbose=False,
        minimize_primers=True,
    )

    assert result.primers == primers, "the ensemble winner was lost"
    assert trimmed, "the minimisation post-process never ran"
    # It needs SOMETHING to measure coverage with. None is acceptable only if
    # the callee tolerates it; an unbound name never is.
    assert "optimizer" in trimmed[0]
