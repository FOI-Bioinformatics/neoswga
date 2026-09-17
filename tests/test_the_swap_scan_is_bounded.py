"""Bound the objective-scored scan with a cheap prescreen in front of it.

Phase 4 increment 4 of the plan for `docs/validation/pipeline_audit_2026-09-16/`,
reframed by measurement.

The plan describes bounding "the per-step scan" of a greedy that scores every
candidate with the objective. That greedy does not exist: `optimize_greedy`
takes `objective=None` and no caller in the package supplies one, so the only
objective-scored searches are this swap loop and the beam. See
`docs/validation/what_actually_bounds_the_search_2026-09-17.md`.

The scan that does need bounding is this one, over `candidates x panel`. At two
thousand candidates and a twelve-primer panel that is twenty-four thousand pairs
per round, and every pair surviving the dimer guard was scored with the full
objective.

The prescreen was already here, as the other branch of the same function: the
binless path computes a bin gain and loss per pair. Bounded mode uses that to
rank the pairs and scores only the leaders. Measured on the real pool, the bin
ranking's top ten pairs were exactly the objective's top ten, in order, so it
filters rather than decides.

What is pinned here is that bounding does not quietly change the answer, that a
width wide enough reproduces the unbounded result, and that a caller asking for
a bound it cannot have is refused rather than silently given an unbounded scan.
"""

import pytest

from neoswga.core.swap_refinement import refine_by_swaps


class _Dimers:
    """No pair dimerises, unless the sequences are named to."""

    def __init__(self, forbidden=()):
        self._forbidden = {frozenset(pair) for pair in forbidden}

    def dimerises(self, incoming, retained):
        return any(frozenset((incoming, held)) in self._forbidden for held in retained)


class _Objective:
    """A counting objective over a table, so evaluations are observable.

    Coverage is a stated number per panel rather than a computation, because
    what is under test is how many times the objective is consulted and which
    panel comes back, not what coverage means.
    """

    def __init__(self, coverage_by_primer, background_by_primer=None, infeasible=()):
        self._coverage = coverage_by_primer
        self._background = background_by_primer or {}
        self._infeasible = {frozenset(p) for p in infeasible}
        self.calls = 0
        self.panels = []

    def _count(self, panel):
        self.calls += 1
        self.panels.append(tuple(sorted(panel)))

    def coverage(self, panel):
        self._count(panel)
        return sum(self._coverage.get(p, 0.0) for p in panel)

    def violations(self, panel):
        self._count(panel)
        return ("infeasible",) if frozenset(panel) in self._infeasible else ()

    def shortfall(self, panel):
        """Ranked on by both searches; see `PoolObjective.shortfall`.

        Feasibility here is a set membership rather than a distance, so the
        count is the honest magnitude: zero when nothing is violated, one when
        something is.
        """
        return float(len(self.violations(panel)))

    def metrics(self, panel):
        self._count(panel)

        class _M:
            total_bg_sites = sum(self._background.get(p, 0) for p in panel)

        return _M()

    @property
    def distinct_panels(self):
        return len(set(self.panels))


def _pool(n):
    return [f"P{i:03d}" for i in range(n)]


def _bins(pool, width=10):
    """Each primer owns one bin, so a swap's bin gain is its own weight."""
    return {p: {i} for i, p in enumerate(pool)}, {i: width for i in range(len(pool))}


# -- the contract ---------------------------------------------------------


def test_a_width_without_bins_is_refused():
    """The prescreen needs bins; running unbounded instead would hide that.

    `plan_pool` passes no bins on the objective path today, so a silent
    fallback here would mean the one production caller that asked to be bounded
    quietly was not.
    """
    pool = _pool(20)
    objective = _Objective({p: 1.0 for p in pool})

    with pytest.raises(ValueError, match="bins"):
        refine_by_swaps(
            pool[:4],
            pool,
            None,
            None,
            _Dimers(),
            objective=objective,
            objective_scan_width=8,
        )


def test_a_width_without_an_objective_is_refused():
    """Nothing to put behind the prescreen."""
    pool = _pool(20)
    bins, weights = _bins(pool)

    with pytest.raises(ValueError, match="objective"):
        refine_by_swaps(pool[:4], pool, bins, weights, _Dimers(), objective_scan_width=8)


@pytest.mark.parametrize("width", [0, -1, 2.5, "8"])
def test_a_width_that_is_not_a_positive_integer_is_refused(width):
    pool = _pool(20)
    bins, weights = _bins(pool)
    objective = _Objective({p: 1.0 for p in pool})

    with pytest.raises(ValueError):
        refine_by_swaps(
            pool[:4],
            pool,
            bins,
            weights,
            _Dimers(),
            objective=objective,
            objective_scan_width=width,
        )


# -- the bound itself -----------------------------------------------------


def test_the_scan_evaluates_at_most_the_width_per_round():
    """The property the whole increment exists for."""
    pool = _pool(60)
    bins, weights = _bins(pool)
    # Coverage rises along the pool, so later primers are always worth taking
    # and the loop keeps finding an improving swap.
    objective = _Objective({p: float(i) for i, p in enumerate(pool)})

    result = refine_by_swaps(
        pool[:4],
        pool,
        bins,
        weights,
        _Dimers(),
        objective=objective,
        objective_scan_width=5,
    )

    assert result.swaps > 0, "nothing was refined, so the bound was not exercised"
    # One round scores at most `width` panels; each panel costs three calls
    # through the lexicographic score. Rounds are bounded by the swaps taken
    # plus the final round that finds nothing.
    assert result.objective_evaluations <= 5 * (result.swaps + 1)


def test_the_bound_is_on_evaluations_not_on_pairs_considered():
    """The prescreen must see the whole frontier, or it is biased.

    The pair budget was calibrated when every pair cost a full evaluation. If
    it also cut the cheap pass short, the prescreen would rank only a prefix of
    the pool and the bound would reintroduce exactly the blindness it removes.
    """
    pool = _pool(40)
    bins, weights = _bins(pool)
    objective = _Objective({p: float(i) for i, p in enumerate(pool)})

    result = refine_by_swaps(
        pool[:4],
        pool,
        bins,
        weights,
        _Dimers(),
        objective=objective,
        objective_scan_width=3,
        max_evaluations=6,
    )

    assert result.pairs_considered > result.objective_evaluations
    assert result.pairs_considered >= 36 * 4 - 16, (
        "the cheap pass stopped early, so the prescreen saw only part of the pool"
    )


def test_a_wide_enough_scan_reproduces_the_unbounded_answer():
    """Bounding is an optimisation, and has to behave like one at the limit."""
    pool = _pool(30)
    bins, weights = _bins(pool)
    coverage = {p: float((i * 7) % 30) for i, p in enumerate(pool)}

    unbounded = refine_by_swaps(
        pool[:5], pool, bins, weights, _Dimers(), objective=_Objective(dict(coverage))
    )
    bounded = refine_by_swaps(
        pool[:5],
        pool,
        bins,
        weights,
        _Dimers(),
        objective=_Objective(dict(coverage)),
        objective_scan_width=len(pool) * 5,
    )

    assert set(bounded.primers) == set(unbounded.primers)


def test_the_prescreen_does_not_hand_over_a_dimerising_pair():
    """The guard runs in the cheap pass, so a rejected pair never ranks."""
    pool = _pool(20)
    bins, weights = _bins(pool)
    objective = _Objective({p: float(i) for i, p in enumerate(pool)})
    # The best candidate by coverage dimerises with a primer that starts in the
    # panel, so an unscreened prescreen would rank it first and waste the width.
    forbidden = [(pool[19], pool[0])]

    result = refine_by_swaps(
        pool[:4],
        pool,
        bins,
        weights,
        _Dimers(forbidden),
        objective=objective,
        objective_scan_width=4,
    )

    if pool[19] in result.primers:
        assert pool[0] not in result.primers


def test_a_feasible_panel_is_not_traded_for_an_infeasible_one():
    """Constraints stay ahead of coverage, bounded or not."""
    pool = _pool(20)
    bins, weights = _bins(pool)
    start = pool[:4]
    tempting = [*start[:3], pool[19]]
    objective = _Objective(
        {p: float(i) for i, p in enumerate(pool)},
        infeasible=[tempting],
    )

    result = refine_by_swaps(
        start,
        pool,
        bins,
        weights,
        _Dimers(),
        objective=objective,
        objective_scan_width=8,
    )

    assert not objective.violations(list(result.primers))


def test_the_result_is_never_worse_than_what_it_started_from():
    """A bounded search may find less; it must not lose ground."""
    pool = _pool(40)
    bins, weights = _bins(pool)
    coverage = {p: float((i * 13) % 40) for i, p in enumerate(pool)}
    start = pool[:6]

    result = refine_by_swaps(
        start,
        pool,
        bins,
        weights,
        _Dimers(),
        objective=_Objective(dict(coverage)),
        objective_scan_width=2,
    )

    before = sum(coverage[p] for p in start)
    after = sum(coverage[p] for p in result.primers)
    assert after >= before


def test_the_stop_reason_names_the_scan_width_when_it_binds():
    """A reader has to be able to tell a narrow scan from a finished one."""
    pool = _pool(40)
    bins, weights = _bins(pool)
    objective = _Objective({p: float(i) for i, p in enumerate(pool)})

    result = refine_by_swaps(
        pool[:4],
        pool,
        bins,
        weights,
        _Dimers(),
        objective=objective,
        objective_scan_width=1,
    )

    assert result.scan_width == 1


def test_unbounded_mode_reports_no_scan_width():
    pool = _pool(20)
    bins, weights = _bins(pool)

    result = refine_by_swaps(
        pool[:4], pool, bins, weights, _Dimers(), objective=_Objective({p: 1.0 for p in pool})
    )

    assert result.scan_width is None
    assert result.objective_evaluations > 0
