"""Forward sites extend downstream, reverse sites extend upstream.

phi29 extends from the 3' end of the primer in one direction along the
template, so a binding site does not cover the genome symmetrically around
itself. A bidirectional model over-estimates coverage by up to 2x in a
strand-imbalanced set.

These four cases were previously written against
`BackgroundAwareOptimizer._calculate_coverage`, a class no dispatch path
reached, which was deleted on 2026-09-10. They are ported here to
`BipartiteGraph.add_primer_coverage`, which is the implementation the four
optimizers that ship actually use, and which had no direct test of this rule at
all -- the only coverage of a real biological rule was against a dead copy of
it.

The bin size matters to the arithmetic and is pinned at 1 bp here, so the
assertions are about the extension rule rather than about binning.
"""

import numpy as np
import pytest

from neoswga.core.dominating_set_optimizer import BipartiteGraph

GENOME = 100_000
REACH = 3_000
SITE = 5_000


def _covered_fraction(forward, reverse):
    """Fraction of a 100 kb genome covered by one primer's sites."""
    graph = BipartiteGraph(bin_size=1)
    graph.add_primer_coverage(
        "PRIMER",
        positions=np.array(sorted(set(forward) | set(reverse)), dtype=np.int64),
        genome_id="fg",
        genome_length=GENOME,
        extension_reach=REACH,
        forward_positions=np.array(forward, dtype=np.int64),
        reverse_positions=np.array(reverse, dtype=np.int64),
    )
    return len(graph.primer_to_regions["PRIMER"]) / GENOME


def test_a_forward_site_extends_downstream_only():
    """[5000, 8000), not [2000, 8000). Bidirectional would read 0.06."""
    assert _covered_fraction(forward=[SITE], reverse=[]) == pytest.approx(0.03, abs=1e-3)


def test_a_reverse_site_extends_upstream_only():
    """[2000, 5000)."""
    assert _covered_fraction(forward=[], reverse=[SITE]) == pytest.approx(0.03, abs=1e-3)


def test_a_site_on_both_strands_covers_both_directions():
    """The union, [2000, 8000). This is the only case where the bidirectional
    answer happens to be right, which is why a test using it alone would not
    have caught the defect."""
    assert _covered_fraction(forward=[SITE], reverse=[SITE]) == pytest.approx(0.06, abs=1e-3)


def test_no_sites_cover_nothing():
    assert _covered_fraction(forward=[], reverse=[]) == 0.0


def test_extension_does_not_run_off_either_end():
    """A site within one reach of a boundary is clipped, not wrapped."""
    near_start = _covered_fraction(forward=[], reverse=[1_000])
    near_end = _covered_fraction(forward=[GENOME - 1_000], reverse=[])

    assert near_start == pytest.approx(0.01, abs=1e-3)
    assert near_end == pytest.approx(0.01, abs=1e-3)
