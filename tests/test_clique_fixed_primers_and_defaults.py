"""Two small clique defects: an off-by-one on fixed primers, and split defaults.

`optimize` reduced its search target by the number of fixed primers with
`max(1, target - len(fixed))`. The floor of 1 meant that fixing as many primers
as were asked for still requested one more from the search, so `--num-primers 4`
with four fixed primers returned five. Zero is the correct floor: fixing the
whole set is a complete answer.

Separately, `build_compatibility_graph`'s own defaults were `max_dimer_bp=3,
max_self_dimer_bp=4` while `OptimizerConfig`'s are 4 and 5. A direct caller of
the helper therefore got a STRICTER graph than the optimizer built from the same
unset configuration -- two answers to "what counts as a dimer here", declared a
few lines apart in the same file.
"""

import dataclasses
import inspect

import pytest

pytest.importorskip("networkx")

from neoswga.core.base_optimizer import OptimizerConfig
from neoswga.core.clique_optimizer import build_compatibility_graph


@pytest.mark.parametrize("field", ["max_dimer_bp", "max_self_dimer_bp"])
def test_the_helper_default_matches_the_config_default(field):
    """One answer to "what counts as a dimer", not two."""
    helper = inspect.signature(build_compatibility_graph).parameters[field].default
    config = {f.name: f.default for f in dataclasses.fields(OptimizerConfig)}[field]

    assert helper == config, (
        f"build_compatibility_graph defaults {field}={helper} while "
        f"OptimizerConfig defaults it to {config}; a direct caller and the "
        "optimizer disagree about the threshold"
    )


def _target_after_fixing(target, n_fixed):
    """The production expression, reached through the optimizer.

    Kept as a direct check of the arithmetic because reaching the branch
    end-to-end needs a genome, a cache and a pool that admits a clique, none of
    which this off-by-one depends on.
    """
    from neoswga.core import clique_optimizer

    source = inspect.getsource(clique_optimizer.CliqueOptimizer.optimize)
    assert "max(0, target - len(fixed))" in source, (
        "the fixed-primer reduction is not floored at 0; `max(1, ...)` asks the "
        "search for one more primer than was requested once the fixed set fills "
        "the target"
    )
    return max(0, target - n_fixed)


@pytest.mark.parametrize(
    "target,n_fixed,expected",
    [(6, 0, 6), (6, 2, 4), (4, 4, 0), (4, 5, 0)],
)
def test_fixing_the_whole_set_asks_the_search_for_nothing_more(target, n_fixed, expected):
    assert _target_after_fixing(target, n_fixed) == expected
