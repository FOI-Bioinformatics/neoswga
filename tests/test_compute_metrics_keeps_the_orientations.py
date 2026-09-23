"""`compute_metrics` pooled the orientations, so nothing downstream could use them.

`_union_coverage` has taken a geometry since PR #107 and no production path
supplied one, because `compute_metrics` unions the two orientations into a set
per prefix before calling it. This connects the two, and clears the debt three
merged increments have been carrying: `site_spans`, `mark_span` and the
geometry argument were all reachable from nothing.

`coverage_geometry` lives on `OptimizerConfig` beside `extension_reach` and
`fg_circular`, which are the other two settings that decide what a coverage
figure means. It defaults to `symmetric`, and the default path still makes ONE
cache lookup per primer and prefix: the two orientations are fetched separately
only when they are going to be used, because `get_positions` is the hottest
method in the package and a run that does not want the second geometry should
not pay for it.

Site COUNTS must not move under either geometry. `total_fg_sites` is a count of
binding sites, not of covered bases, and it is the denominator of
`selectivity_ratio`. Under the directional geometry the pooled list is derived
locally as `sorted(set(forward) | set(reverse))`, which is what
`get_positions(..., "both")` returns by construction.
"""

import numpy as np
import pytest

from neoswga.core.base_optimizer import OptimizerConfig
from neoswga.core.dominating_set_adapter import DominatingSetAdapter
from neoswga.core.position_cache import PositionCache

GENOME = 40_000
PRIMER = "ACCACAGATAGC"  # not a palindrome; its rc is a distinct dataset


@pytest.fixture
def world(tmp_path):
    """One prefix whose forward and reverse sites are deliberately unbalanced.

    Balanced sites make the two geometries agree by accident -- which is what
    `docs/validation/2026-09-23-directional-coverage.md` measured on real
    panels and why the headline figure barely moves there. A skewed prefix is
    what makes a wiring mistake visible.
    """
    import h5py

    prefix = str(tmp_path / "target")
    from neoswga.core.thermodynamics import reverse_complement

    with h5py.File(prefix + "_12mer_positions.h5", "w") as handle:
        handle.create_dataset(PRIMER, data=np.array([5_000, 9_000, 13_000], dtype=np.int64))
        handle.create_dataset(
            reverse_complement(PRIMER), data=np.array([31_000], dtype=np.int64)
        )
    return prefix


def optimizer(world, geometry):
    return DominatingSetAdapter(
        position_cache=PositionCache([world], [PRIMER]),
        fg_prefixes=[world],
        fg_seq_lengths=[GENOME],
        bg_prefixes=[],
        bg_seq_lengths=[],
        config=OptimizerConfig(
            target_set_size=1, extension_reach=2_000, coverage_geometry=geometry
        ),
        conditions=None,
    )


def test_the_default_geometry_is_symmetric():
    """Every recorded figure was produced under it, and the reach was fitted
    alongside it. A new default would silently recalibrate every saved run."""
    assert OptimizerConfig(target_set_size=1).coverage_geometry == "symmetric"


def test_the_symmetric_coverage_is_what_it_always_was(world):
    """Four sites, `2 * 2000` wide each, none overlapping: 16,000 of 40,000."""
    metrics = optimizer(world, "symmetric").compute_metrics([PRIMER])

    assert metrics.fg_coverage == pytest.approx(16_000 / GENOME)


def test_the_directional_coverage_reaches_one_way(world):
    """Same four sites, `2000` wide each and placed on the reachable side:
    8,000 of 40,000. Half, because none of these spans overlap either."""
    metrics = optimizer(world, "directional").compute_metrics([PRIMER])

    assert metrics.fg_coverage == pytest.approx(8_000 / GENOME)


def test_the_geometry_reaches_the_metric_at_all(world):
    """The point of the increment. Before it, both configs gave one answer."""
    symmetric = optimizer(world, "symmetric").compute_metrics([PRIMER]).fg_coverage
    directional = optimizer(world, "directional").compute_metrics([PRIMER]).fg_coverage

    assert symmetric != directional


def test_the_site_count_does_not_move_with_the_geometry(world):
    """`total_fg_sites` counts binding sites, not covered bases, and it is the
    denominator of `selectivity_ratio`. Deriving the pooled list per orientation
    must reproduce what `get_positions(..., "both")` returned."""
    symmetric = optimizer(world, "symmetric").compute_metrics([PRIMER])
    directional = optimizer(world, "directional").compute_metrics([PRIMER])

    assert symmetric.total_fg_sites == 4
    assert directional.total_fg_sites == symmetric.total_fg_sites


def test_an_unknown_geometry_is_refused(world):
    """A typo must not fall through to a silent default. `fg_coverage` is the
    authoritative figure in every saved summary."""
    with pytest.raises(ValueError):
        optimizer(world, "radial").compute_metrics([PRIMER])
