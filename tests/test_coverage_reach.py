"""Coverage is meaningless without the reach it was computed at.

The published SWGA tools do not share a convention. swga 2.0's
`coverage_ratio` uses phi29's ~70 kb single-molecule processivity as the
per-primer reach; NeoSWGA selects on the ~3 kb realistic per-primer reach,
arguing that one molecule being extendable 70 kb is not the reach at which a
primer reliably covers a genome in a dense multi-primer reaction.

Both positions are defensible and they are not reconcilable by argument,
because the same primer set measures:

    3 kb   0.418
   10 kb   0.836
   70 kb   ~1.0

A bare "coverage" number therefore cannot be compared across tools in either
direction, and a design chosen at one reach is not the design that would be
chosen at another. Two changes follow:

- the reach is explicit and configurable (`coverage_reach` / `--coverage-reach`)
  rather than a constant derived from the polymerase, because the right value is
  an empirical property of the reaction rather than of the enzyme alone;
- coverage is reported at several reaches, so a NeoSWGA figure can be placed
  beside a published one without either convention having to move.

The default is unchanged. Nothing here shifts a result unless a reach is asked
for.
"""

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")

from neoswga.core.base_optimizer import OptimizerConfig
from neoswga.core.coverage import REPORTING_REACHES, resolve_coverage_reach
from neoswga.core.dominating_set_adapter import DominatingSetAdapter
from neoswga.core.position_cache import PositionCache
from neoswga.core.reaction_conditions import ReactionConditions
from neoswga.core.thermodynamics import reverse_complement

GENOME = 200_000
K = 10
PRIMER = "GCGCGGCCGC"


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    """One primer with widely spaced sites, so reach visibly drives coverage."""
    import random

    rng = random.Random(11)
    directory = tmp_path_factory.mktemp("reach")
    seq = [rng.choice("ACGT") for _ in range(GENOME)]
    for j in range(4):
        at = 20_000 + j * 40_000
        seq[at : at + K] = list(PRIMER)

    prefix = str(directory / "fg")
    joined = "".join(seq)
    with h5py.File(f"{prefix}_{K}mer_positions.h5", "w") as handle:
        for key in sorted({PRIMER, reverse_complement(PRIMER)}):
            hits, at = [], joined.find(key)
            while at != -1:
                hits.append(at)
                at = joined.find(key, at + 1)
            if hits:
                handle.create_dataset(key, data=np.array(hits, dtype=np.int32))
    return prefix


def metrics_at(prefix, reach):
    cache = PositionCache([prefix], [PRIMER])
    return DominatingSetAdapter(
        position_cache=cache,
        fg_prefixes=[prefix],
        fg_seq_lengths=[GENOME],
        bg_prefixes=[],
        bg_seq_lengths=[],
        config=OptimizerConfig(target_set_size=1, extension_reach=reach),
        conditions=ReactionConditions(temp=42.0, polymerase="equiphi29"),
    ).compute_metrics([PRIMER])


# ----------------------------------------------------------------------
# Resolution
# ----------------------------------------------------------------------


def test_the_default_is_the_polymerase_reach_unchanged():
    """Nothing moves for a user who does not ask for a reach."""
    from neoswga.core.coverage import polymerase_extension_reach

    for polymerase in ("phi29", "equiphi29"):
        assert resolve_coverage_reach(polymerase) == polymerase_extension_reach(
            polymerase, coverage_metric="realistic"
        )


def test_an_explicit_reach_wins():
    assert resolve_coverage_reach("phi29", override=10000) == 10000


@pytest.mark.parametrize("bad", ["10kb", 0, -1, True, 0.0])
def test_a_nonsense_reach_is_rejected_rather_than_applied(bad):
    """Silently falling back would produce a design at a different reach than
    the user asked for, and coverage numbers that cannot be compared to
    anything. `True` is included because bool is an int in Python."""
    with pytest.raises(ValueError, match="coverage_reach"):
        resolve_coverage_reach("phi29", override=bad)


# ----------------------------------------------------------------------
# Reporting
# ----------------------------------------------------------------------


def test_coverage_rises_with_reach(world):
    """The fact that makes a bare coverage figure uninterpretable."""
    values = [metrics_at(world, r).fg_coverage for r in (3000, 10000, 70000)]

    assert values == sorted(values)
    assert values[0] < values[-1]


def test_every_reporting_reach_is_present(world):
    m = metrics_at(world, 3000)

    for reach in REPORTING_REACHES:
        assert reach in m.coverage_by_reach


def test_the_selection_reach_appears_on_the_curve(world):
    """`fg_coverage` must be locatable on the reported curve. A headline number
    that appears on no reported reach reads as coming from nowhere."""
    m = metrics_at(world, 5000)

    assert 5000 in m.coverage_by_reach
    assert m.coverage_by_reach[5000] == pytest.approx(m.fg_coverage)


def test_the_reported_curve_agrees_with_selecting_at_that_reach(world):
    """The curve and a run actually performed at that reach must agree, or one
    of them is measuring something else."""
    curve = metrics_at(world, 3000).coverage_by_reach

    for reach in (10000, 70000):
        assert curve[reach] == pytest.approx(metrics_at(world, reach).fg_coverage)


def test_it_serializes_for_the_report(world):
    """JSON object keys must be strings; the report reads this dict."""
    payload = metrics_at(world, 3000).to_dict()

    assert "coverage_by_reach" in payload
    assert all(isinstance(k, str) for k in payload["coverage_by_reach"])
    assert payload["extension_reach"] == 3000


def test_one_implementation_of_union_coverage(world):
    """`_compute_coverage` and `_compute_coverage_at` must not drift. Two
    implementations of one quantity is how this codebase has produced
    disagreeing coverage numbers before -- three different semantics were found
    in one audit."""
    cache = PositionCache([world], [PRIMER])
    opt = DominatingSetAdapter(
        position_cache=cache,
        fg_prefixes=[world],
        fg_seq_lengths=[GENOME],
        bg_prefixes=[],
        bg_seq_lengths=[],
        config=OptimizerConfig(target_set_size=1, extension_reach=7000),
        conditions=None,
    )
    positions = sorted({int(x) for x in cache.get_positions(world, PRIMER, "both")})

    assert opt._compute_coverage(positions, GENOME) == pytest.approx(
        opt._compute_coverage_at(positions, GENOME, 7000)
    )
