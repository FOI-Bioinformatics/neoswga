"""`--coverage-reach` changed the selection and not the report.

`run_optimization` resolves the reach once, honouring `coverage_reach` from
params.json or `--coverage-reach`, and threads it into the optimizer config. It
then recomputed a *second* reach from the polymerase alone for
`per_target_coverage`:

    _ext = polymerase_extension_reach(parameter.polymerase, default=extension_reach)

so `--coverage-reach 8000` selected primers at 8 kb and reported their coverage
at phi29's 3 kb default. Coverage is roughly linear in reach over this range, so
the reported figure was about a third of the one the set was chosen for, with
nothing saying the two numbers came from different assumptions.

The reach is the single largest lever on any coverage figure this tool prints --
the same primer set measures 0.418 at 3 kb, 0.836 at 10 kb and ~1.0 at 70 kb --
so two of them disagreeing silently inside one run is worse than either value
being wrong.
"""

import numpy as np
import pytest

from neoswga.core.coverage import compute_per_prefix_coverage


class _Cache:
    """Position cache stand-in with one primer at fixed, sparse positions."""

    def __init__(self, positions):
        self._positions = positions

    def get_positions(self, prefix, primer, strand="both"):
        return np.asarray(self._positions, dtype=np.int64)


# Spaced so that at every reach below, the windows neither overlap each other
# nor clip against the ends -- the expectation is then just 2*reach per site.
GENOME = 400_000
POSITIONS = [50_000, 150_000, 250_000, 350_000]


@pytest.mark.parametrize("reach", [1000, 3000, 8000, 20_000])
def test_coverage_tracks_the_reach_it_is_given(reach):
    """The premise: this is not a number that survives being recomputed."""
    overall, _ = compute_per_prefix_coverage(
        cache=_Cache(POSITIONS),
        primers=["ATGCATGCATGC"],
        prefixes=["fg"],
        seq_lengths=[GENOME],
        extension=reach,
    )

    expected = min(1.0, len(POSITIONS) * 2 * reach / GENOME)
    assert overall == pytest.approx(expected, abs=0.02)


def test_the_reporting_path_uses_the_resolved_reach():
    """Guards the shape: `run_optimization` must not resolve a second reach for
    reporting after it has already resolved one for selection.
    """
    import ast
    from pathlib import Path

    source = (
        Path(__file__).parent.parent / "neoswga" / "core" / "unified_optimizer.py"
    ).read_text()
    tree = ast.parse(source)

    calls = [
        node.lineno
        for node in ast.walk(tree)
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Name)
        and node.func.id == "polymerase_extension_reach"
    ]

    assert not calls, (
        "unified_optimizer resolves the reach once, from `coverage_reach` when "
        f"the user set it; recomputing it from the polymerase at line(s) {calls} "
        "discards the override for whatever is computed there"
    )


def test_an_override_is_carried_into_the_optimizer_config():
    """The override has to reach the config that selection reads."""
    from neoswga.core.coverage import resolve_coverage_reach

    assert resolve_coverage_reach("phi29", override=8000) == 8000
    assert resolve_coverage_reach("bst", override=8000) == 8000
    assert resolve_coverage_reach("phi29", override=None) == 3000
