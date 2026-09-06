"""A method that produced no primers must not be rankable as the winner.

`_run_ensemble` entered every result that did not RAISE into the candidate set,
and ranked on `res.metrics.normalized_score(...)` rather than on
`OptimizationResult.normalized_score`, the property that returns 0.0 for ERROR
and NO_CONVERGENCE.

That distinction matters because an empty `PrimerSetMetrics` does not score
zero. A `(0.0, 0.0)` Tm range reads as perfectly tight, so it collects the whole
Tm term:

    empty metrics, balanced       0.1000
    empty metrics, clinical       0.0500
    empty metrics, metagenomics   0.0500

Hybrid, network and background-aware all return NO_CONVERGENCE with an empty
primer tuple rather than raising, so such a result was a ranked candidate. When
every method fails that way the ensemble returns a "winner" holding no primers,
and `optimize_step4` then logs "Optimization failed: " with an empty message,
instead of the clear all-methods-failed error that already exists.

A real set scores well above 0.10, so an empty result out-ranking a good one is
unlikely rather than impossible. The failure that is not unlikely is the
all-fail case.
"""

from dataclasses import replace

import pytest

from neoswga.core.base_optimizer import OptimizationResult, OptimizationStatus, PrimerSetMetrics


def _result(name, primers, status=OptimizationStatus.SUCCESS, coverage=0.9):
    return OptimizationResult(
        primers=tuple(primers),
        score=1.0,
        metrics=replace(PrimerSetMetrics.empty(), fg_coverage=coverage),
        status=status,
        optimizer_name=name,
        iterations=1,
    )


def _run(results_by_method, monkeypatch):
    from neoswga.core import unified_optimizer as uo

    class _Stub:
        def __init__(self, method):
            # `progress_context` labels the run with this.
            self.name = method

        def optimize(self, candidates, target_size=None, **kw):
            return results_by_method[self.name]

    def fake_create(name, **kwargs):
        return _Stub(name)

    monkeypatch.setattr(uo.OptimizerFactory, "create", staticmethod(fake_create))
    return uo._run_ensemble(
        methods=list(results_by_method),
        cache=None,
        candidates=["AAACCCGGGT"],
        fg_prefixes=["fg"],
        fg_seq_lengths=[10_000],
        bg_prefixes=[],
        bg_seq_lengths=[],
        target_size=2,
        config=None,
        conditions=None,
        verbose=False,
    )


def test_an_empty_result_never_beats_a_real_one(monkeypatch):
    outcome = _run(
        {
            # Listed first, so order alone would hand it the win if it were
            # allowed to tie.
            "hybrid": _result("hybrid", [], status=OptimizationStatus.NO_CONVERGENCE, coverage=0.0),
            "network": _result("network", ["AAACCCGGGT", "ACCCGGGTTT"], coverage=0.8),
        },
        monkeypatch,
    )

    assert outcome.primers, "the ensemble returned a set with no primers"
    assert outcome.optimizer_name == "network"


def test_when_every_method_produces_nothing_the_result_is_a_failure(monkeypatch):
    """The clear error already exists; it was simply unreachable."""
    outcome = _run(
        {
            "hybrid": _result("hybrid", [], status=OptimizationStatus.NO_CONVERGENCE, coverage=0.0),
            "network": _result(
                "network", [], status=OptimizationStatus.NO_CONVERGENCE, coverage=0.0
            ),
        },
        monkeypatch,
    )

    assert not outcome.primers
    assert outcome.status is OptimizationStatus.ERROR, (
        f"status is {outcome.status}; an ensemble in which nothing produced a "
        "set is a failure, not a success with an empty winner"
    )
    assert outcome.message, "the failure carries no message to show the user"
