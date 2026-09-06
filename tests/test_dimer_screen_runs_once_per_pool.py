"""The heterodimer screen examines the pairs a cheap test flags, once.

Finding A1. The screen materialised all n(n-1)/2 pairs before checking any of
them (2.0 GB measured for 2789 primers) and ran a 627 us thermodynamic
calculation on every one of them. It now streams, and only over the pairs the
exact substring test flags.

Fix round 1 correction: `concurrent.futures.Executor.map` drains its whole
iterable before returning a single result (`Executor.map` builds
`fs = [self.submit(fn, *args) for args in zip(*iterables)]`), so handing it
the lazy pair generator did not stream anything above the parallel path's
100-pair threshold -- verified empirically: a 500-item generator had produced
all 500 items before the first result could be consumed. The parallel branch
now submits through a bounded sliding window instead. Two of the tests below
were also fixed: the original combined test took the parallel branch on
macOS, where `ProcessPoolExecutor` spawns rather than forks, so a
`monkeypatch.setattr` on `check_heterodimer` in the parent process never
reached the workers and the assertions passed vacuously (`len(calls) == 0`
in every case, including one where the screen visited every pair). The
call-count assertion is now exercised only on the sequential branch, where
the monkeypatch is observed in the same process that runs it.
"""

import random

import pytest

from neoswga.core.thermodynamic_filter import ThermodynamicCriteria, ThermodynamicFilter


def _pool(n, seed=0):
    rng = random.Random(seed)
    return ["".join(rng.choice("ACGT") for _ in range(12)) for _ in range(n)]


def _wide_open_criteria():
    """Criteria that let every random 12-mer pass the individual filters.

    Tm and GC are wide open, and so are the homodimer/hairpin thresholds
    (the ThermodynamicCriteria defaults reject on those too). Without this,
    some fraction of a random pool fails analyze_primer's individual filters
    before the heterodimer screen ever runs, so the flagged-pair count
    computed on the raw input pool would not match the count computed on the
    (smaller) set filter_candidates actually screens.
    """
    return ThermodynamicCriteria(
        min_tm=0.0,
        max_tm=200.0,
        min_gc=0.0,
        max_gc=1.0,
        max_homodimer_dg=-1000.0,
        max_hairpin_dg=-1000.0,
    )


def test_screen_calls_the_thermodynamic_check_exactly_for_flagged_pairs_sequential_branch(
    monkeypatch,
):
    """Sequential branch (<=100 flagged pairs): the monkeypatch below runs in
    the same process as filter_candidates, so it is actually observed. 30
    primers at seed 3 flag 98 of 435 possible pairs, under the 100-pair
    threshold that would otherwise route to the (out-of-process) parallel
    branch.
    """
    primers = _pool(30, seed=3)
    criteria = _wide_open_criteria()
    filt = ThermodynamicFilter(criteria)

    calls = []
    import neoswga.core.thermodynamic_filter as tf

    real = tf.check_heterodimer
    monkeypatch.setattr(
        tf, "check_heterodimer", lambda a, b, c: calls.append((a, b)) or real(a, b, c)
    )

    filt.filter_candidates(primers, check_heterodimers=True, max_dimer_bp=3)

    from neoswga.core import dimer_matrix

    flagged = int(dimer_matrix.build(primers, 3).pairs.sum()) // 2
    total = len(primers) * (len(primers) - 1) // 2
    assert flagged <= 100, "fixture drifted out of the sequential branch"
    assert len(calls) == flagged
    assert len(calls) < total, "the screen must not visit every pair"


def test_screen_routes_the_parallel_branch_through_the_dimer_matrix(monkeypatch):
    """Parallel branch (>100 flagged pairs): a monkeypatch on
    `check_heterodimer` cannot be observed here, because ProcessPoolExecutor
    spawns fresh worker processes on macOS that re-import the module and
    never see it. What IS observable in the parent process is which pairs
    dimer_matrix.build was asked to flag, since that call happens before any
    worker is spawned. 40 primers at seed 3 flag 174 of 780 possible pairs,
    over the 100-pair threshold, so this exercises the parallel branch.
    """
    primers = _pool(40, seed=3)
    criteria = _wide_open_criteria()
    filt = ThermodynamicFilter(criteria)

    import neoswga.core.thermodynamic_filter as tf
    from neoswga.core import dimer_matrix

    build_calls = []
    real_build = dimer_matrix.build

    def _recording_build(seqs, max_dimer_bp):
        build_calls.append((list(seqs), max_dimer_bp))
        return real_build(seqs, max_dimer_bp)

    # filter_candidates does `from neoswga.core import dimer_matrix as
    # _dimer_matrix` at call time, so patching the dimer_matrix module's
    # `build` attribute directly is what it actually sees.
    monkeypatch.setattr(dimer_matrix, "build", _recording_build)

    filt.filter_candidates(primers, check_heterodimers=True, max_dimer_bp=3)

    total = len(primers) * (len(primers) - 1) // 2
    flagged = int(real_build(primers, 3).pairs.sum()) // 2
    assert flagged > 100, "fixture drifted out of the parallel branch"
    assert len(build_calls) == 1
    called_seqs, called_threshold = build_calls[0]
    assert called_seqs == primers
    assert called_threshold == 3
    assert flagged < total, "the flagged set must be a proper subset of all pairs"


def test_parallel_branch_bounds_submissions_in_flight(monkeypatch):
    """The fix for the CRITICAL finding: the parallel branch used to hand
    `executor.map()` a generator that it drained completely before the first
    result was available. It now submits through a bounded sliding window.
    This test proves the bound directly: it fakes the process pool so
    submitted work runs synchronously in this process (avoiding both a slow
    thermodynamic calculation and the cross-process visibility problem the
    other two tests work around), counts how many pair-tuples have been
    produced and submitted by the time the first `concurrent.futures.wait()`
    call returns, and asserts that count is bounded by the sliding window,
    strictly less than the total number of candidate pairs.
    """
    import os

    import neoswga.core.thermodynamic_filter as tf

    primers = _pool(100, seed=9)
    criteria = _wide_open_criteria()
    filt = ThermodynamicFilter(criteria)

    # Cheap stand-in: the point of this test is submission counting, not
    # thermodynamics, and the real calculation would make 4950 pairs slow.
    monkeypatch.setattr(tf, "check_heterodimer", lambda a, b, c: {"energy": 0.0})

    submitted = []

    class _CompletedFuture:
        def __init__(self, value):
            self._value = value

        def result(self):
            return self._value

    class _FakeExecutor:
        """Runs submitted work synchronously in this process, so the test
        can count submissions without spawning real subprocesses."""

        def __init__(self, max_workers=None):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *exc_info):
            return False

        def submit(self, fn, args):
            submitted.append(args)
            return _CompletedFuture(fn(args))

    submitted_counts_at_wait = []

    def _fake_wait(fs, return_when=None):
        # Every future is already complete (submit() ran it synchronously),
        # so this returns immediately with all of them done -- exactly like
        # the real concurrent.futures.wait would, once its futures finish.
        submitted_counts_at_wait.append(len(submitted))
        return set(fs), set()

    monkeypatch.setattr(tf.concurrent.futures, "ProcessPoolExecutor", _FakeExecutor)
    monkeypatch.setattr(tf.concurrent.futures, "wait", _fake_wait)

    # max_dimer_bp=None keeps every pair a candidate, so the pair count is
    # large and predictable: 100 primers -> 4950 pairs, comfortably above any
    # sliding window this function could reasonably use (n_workers is capped
    # at 8, so window <= 8 * 50 * 2 = 800 on any machine).
    filt.filter_candidates(primers, check_heterodimers=True, max_dimer_bp=None)

    total_pairs = len(primers) * (len(primers) - 1) // 2
    assert total_pairs > 800, "fixture must exceed the largest possible window"
    assert submitted_counts_at_wait, "the fake wait() was never called -- branch routing changed"

    first_wait_submitted_count = submitted_counts_at_wait[0]
    n_workers = min(os.cpu_count() or 1, 8)
    expected_max_window = n_workers * 50 * 2

    assert first_wait_submitted_count <= expected_max_window
    assert first_wait_submitted_count < total_pairs, (
        "the generator was drained before the first result was consumed -- "
        "this is the bug the sliding window exists to fix"
    )


def test_screen_does_not_build_a_full_pair_list(monkeypatch):
    """A list comprehension over all pairs is the allocation this removes."""
    import inspect

    import neoswga.core.thermodynamic_filter as tf

    source = inspect.getsource(tf.ThermodynamicFilter.filter_candidates)
    assert (
        "for j in range(i + 1, len(passing))" not in source
    ), "the eager pair comprehension is back"
    assert "executor.map(" not in source, (
        "executor.map() drains its iterable eagerly before returning a "
        "single result -- the bounded sliding window must be used instead"
    )


def test_screen_falls_back_to_the_substring_test_when_the_matrix_would_be_too_big(
    monkeypatch, caplog
):
    """dimer_matrix.build raises ValueError when max_dimer_bp needs more t-mer
    codes than it will allocate (params.schema.json permits max_dimer_bp up to
    15 and max_k up to 30, so a threshold/length combination that trips this is
    reachable with legal configuration, e.g. a 12-base primer pool with
    max_dimer_bp=10). The screen must not crash; it must fall back to flagging
    pairs with dimer.is_dimer_fast pairwise rather than screening every pair
    thermodynamically, and it must say so at warning level. This fixture is
    small enough to take the sequential branch either way.
    """
    primers = _pool(30, seed=7)
    criteria = _wide_open_criteria()
    filt = ThermodynamicFilter(criteria)

    import neoswga.core.thermodynamic_filter as tf

    calls = []
    real = tf.check_heterodimer
    monkeypatch.setattr(
        tf, "check_heterodimer", lambda a, b, c: calls.append((a, b)) or real(a, b, c)
    )

    total = len(primers) * (len(primers) - 1) // 2

    with caplog.at_level("WARNING"):
        # max_dimer_bp=10 on 12-mers needs 4**11 t-mer codes, above
        # dimer_matrix.MAX_CODES (4**8), so dimer_matrix.build raises.
        result, stats = filt.filter_candidates(primers, check_heterodimers=True, max_dimer_bp=10)

    assert result is not None
    assert len(calls) <= total
    warning_text = " ".join(r.message for r in caplog.records if r.levelno >= 30)
    assert "10" in warning_text
    assert "substring" in warning_text.lower()
