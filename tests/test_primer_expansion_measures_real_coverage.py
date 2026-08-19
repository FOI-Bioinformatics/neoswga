"""`expand-primers` measured bin occupancy and called it coverage.

`PrimerExpander._calculate_coverage` built its bipartite graph and then called
`add_primer_coverage(primer, positions, prefix, length)` with **no
`extension_reach`**. With reach 0 that method falls through to "mark the bin the
site sits in", so the number it returned was *bin occupancy* -- "does any site
fall in this bin" -- and its value depended on `bin_size` rather than on the
polymerase. At the 10 kb default a single site claimed 10 kb of coverage; at a
1 kb bin the same design scored a tenth of that.

It also skipped `coverage_bin_size`, the cap that exists because a bin larger
than the reach lets a primer covering part of a bin claim all of it. Both
faults are the ones `hybrid_optimizer._calculate_coverage` was fixed for; this
copy was missed.

Separately, `max_extension` on this class is documented as -- and defaults to --
the *per-primer coverage* reach of 3000 bp. It was then handed to
`HybridOptimizer(max_extension=...)`, which is the amplification-**network**
reach: "could two primers connect via one uninterrupted extension?". That
question is answered by single-molecule processivity, ~70 kb for phi29, so
expansion ran its network stage at a reach 23x too short while leaving
`coverage_reach` to be inferred from the same value. Two different quantities,
one field.
"""

import numpy as np
import pytest

from neoswga.core.primer_expansion import PrimerExpander

GENOME = 400_000
REACH = 3000
# Two sites per primer, far apart, so coverage is exactly 2 * 2 * reach.
PRIMERS = {
    "AAAACCCCGGGG": [50_000, 250_000],
    "TTTTGGGGCCCC": [150_000, 350_000],
}


class _Cache:
    def get_positions(self, prefix, primer, strand="both"):
        pos = PRIMERS.get(primer, [])
        if strand == "reverse":
            return np.array([], dtype=np.int64)
        return np.asarray(pos, dtype=np.int64)


def _expander(**kwargs):
    return PrimerExpander(
        position_cache=_Cache(),
        fg_prefixes=["fg"],
        fg_seq_lengths=[GENOME],
        **kwargs,
    )


def test_coverage_reflects_the_reach_not_the_bin_size():
    """The regression. With no reach threaded, this returned bin occupancy, so
    a 10 kb bin and a 1 kb bin gave answers ten times apart for one design."""
    coarse = _expander(bin_size=10_000)._calculate_coverage(list(PRIMERS))
    fine = _expander(bin_size=1_000)._calculate_coverage(list(PRIMERS))

    assert coarse == pytest.approx(fine, abs=0.05)


@pytest.mark.parametrize("reach", [1000, 3000, 8000])
def test_coverage_tracks_the_configured_reach(reach):
    coverage = _expander(coverage_reach=reach)._calculate_coverage(list(PRIMERS))

    # Four sites, symmetric window, none overlapping at these reaches.
    assert coverage == pytest.approx(4 * 2 * reach / GENOME, abs=0.05)


def test_the_bin_is_capped_against_the_reach():
    """`coverage_bin_size` exists because a bin larger than the reach lets a
    primer covering part of it claim all of it."""
    expander = _expander(bin_size=10_000, coverage_reach=REACH)

    assert expander.coverage_bin_size <= REACH


def test_coverage_never_exceeds_one():
    """What an uncapped bin produces first: a fraction above 1.0."""
    assert _expander(bin_size=10_000)._calculate_coverage(list(PRIMERS)) <= 1.0


def test_the_network_reach_is_not_the_coverage_reach():
    """`HybridOptimizer.max_extension` answers "could two primers connect via
    one extension?", which is single-molecule processivity, not the per-primer
    coverage reach this class defaults to."""
    from neoswga.core.registry import POLYMERASES

    expander = _expander(coverage_reach=REACH)
    optimizer = expander._build_hybrid_optimizer()

    assert optimizer.coverage_reach == REACH
    assert optimizer.max_extension == POLYMERASES["phi29"].processivity_bp


def test_max_extension_still_sets_the_coverage_reach():
    """The old field name meant the coverage reach, and callers pass it. It has
    to keep meaning that rather than silently becoming the network reach."""
    assert _expander(max_extension=5000).coverage_reach == 5000
