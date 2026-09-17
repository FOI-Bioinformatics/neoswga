"""Dimer compatibility computed against the selected panel, not the whole pool.

Task 4 of the condition-aware pool design plan.

`dimer_matrix.build` allocates an n-by-n boolean array over every candidate. At
the shortlist sizes the pipeline used to deliver, that was unremarkable. With
`candidate_retention="all_qc"` the search is handed an order of magnitude more:

===============  ===============
candidates       pairs matrix
===============  ===============
2,000            4 MB
20,670           427 MB
50,000           2,500 MB
===============  ===============

The Wolbachia design has 20,670 candidates clearing the hard gates, so retaining
them all makes that matrix the largest single allocation in the run, and it
grows quadratically from there.

Almost none of it is read. The greedy screens each candidate against the
SELECTED panel, which is tens of primers, so it touches a thin slice of a
quadratic structure. Computing those pairs on demand and caching them bounded
answers the same question for a bounded cost.

This implements the same `dimerises(candidate, selected)` interface
`_would_dimerise` already calls, so it substitutes for the matrix rather than
adding a second screen with its own rules.
"""

from __future__ import annotations

import logging
from typing import Dict, Sequence, Tuple

logger = logging.getLogger(__name__)

# Above this many candidates, compute pairs on demand instead of materialising
# them. 4,000 candidates is a 16-million-entry array; the threshold has lived in
# `dominating_set_optimizer` since that was the only search that respected it.
LAZY_DIMER_POOL_THRESHOLD = 4_000


class LazyDimerCompatibility:
    """Answers `dimerises(candidate, selected)` without materialising all pairs.

    Constructed from a threshold rather than a pool, which is the property that
    makes it independent of pool size.
    """

    def __init__(self, max_dimer_bp: int, cache_size: int = 200_000):
        self.max_dimer_bp = int(max_dimer_bp)
        self.cache_size = int(cache_size)
        self._cache: Dict[Tuple[str, str], bool] = {}
        # How many pairs were actually computed, so a test can show that the
        # cache is used rather than merely present.
        self.computations = 0

    def _pair(self, one: str, two: str) -> bool:
        """Whether these two dimerise, computed once per unordered pair.

        Keyed on the sorted pair: dimerisation is symmetric, and keying on the
        argument order would compute and store each pair twice.
        """
        key = (one, two) if one <= two else (two, one)
        cached = self._cache.get(key)
        if cached is not None:
            return cached

        from neoswga.core.dimer import is_dimer_fast

        result = bool(is_dimer_fast(key[0], key[1], max_dimer_bp=self.max_dimer_bp))
        self.computations += 1
        if len(self._cache) >= self.cache_size:
            # Cleared rather than evicted one at a time. A search's working set
            # is the panel it is building, so the next few hundred lookups
            # repopulate what matters; tracking recency would cost more than it
            # saves here.
            self._cache.clear()
        self._cache[key] = result
        return result

    def dimerises(self, candidate: str, selected: Sequence[str]) -> bool:
        """Whether `candidate` conflicts with anything already chosen."""
        if not selected:
            return False
        upper = str(candidate).upper()
        return any(self._pair(upper, str(other).upper()) for other in selected)


def dimer_screen(pool: Sequence[str], max_dimer_bp: int):
    """The dimer screen suited to this pool's size.

    One decision in one place. `dominating_set_optimizer` made it correctly and
    two other call sites did not: `network_optimizer` and
    `refine_hybrid_stage2` both built the dense matrix unconditionally, and
    `plan-pool` sets ``refinement_method="swap"``, so the second was on the hot
    path. At the 491,836 candidates `candidate_retention="all_qc"` retains, the
    dense array is roughly 242 GB -- a MemoryError, not a slowdown.

    Returns whichever object answers ``dimerises(candidate, selected)``. Both do,
    with the same rule and the same threshold, so a caller never needs to know
    which it got.

    A threshold the dense matrix cannot represent also takes the lazy branch,
    rather than raising. The matrix codes t-mers in a 4**8 space and so cannot
    hold 8 or above; `dimer.is_dimer_fast`, which the lazy screen uses, has no
    such limit. So the pair of them can always enforce what was configured, and
    the screen never has to weaken or refuse. A screen that cannot be enforced
    must not read as "no conflicts": an 11 bp delivered heterodimer against a
    configured 3 is what that looks like.

    params.schema.json refuses 8 and above anyway, so this is a failure mode
    removed rather than a capability added.
    """
    candidates = list(dict.fromkeys(pool))
    if len(candidates) > LAZY_DIMER_POOL_THRESHOLD:
        logger.info(
            "Pool of %d candidates: computing dimer compatibility on demand rather "
            "than materialising %d pairs (%.0f MB).",
            len(candidates),
            len(candidates) ** 2,
            len(candidates) ** 2 / 1e6,
        )
        return LazyDimerCompatibility(max_dimer_bp)
    if not _dense_can_hold(max_dimer_bp, candidates):
        logger.info(
            "max_dimer_bp=%d needs more t-mer codes than the dense matrix allocates; "
            "screening pairwise on demand instead.",
            max_dimer_bp,
        )
        return LazyDimerCompatibility(max_dimer_bp)
    from neoswga.core import dimer_matrix as _dimer_matrix

    return _dimer_matrix.build(candidates, max_dimer_bp)


def _dense_can_hold(max_dimer_bp: int, candidates: Sequence[str]) -> bool:
    """Whether `dimer_matrix.build` can represent this threshold.

    Mirrors the two cases `build` itself distinguishes: a run longer than the
    longest primer cannot occur at all, which it short-circuits to an all-False
    matrix, and anything else needs 4**(max_dimer_bp + 1) codes.
    """
    from neoswga.core import dimer_matrix as _dimer_matrix

    t = int(max_dimer_bp) + 1
    longest = max((len(str(s)) for s in candidates), default=0)
    if t > longest:
        return True
    return 4**t <= _dimer_matrix.MAX_CODES
