"""`tm_weight` and `uniformity_weight` do nothing on the default optimizer.

Found 2026-09-21 while checking an external report's claim that the Tm term
"collapses a mixed-length panel by construction on every default run". The Tm
term is real and its span across lengths is enormous -- median `tm_score` on
the shipped plasmid pool runs from 1.4e-05 at k=7 to 0.697 at k=9, a
48,488-fold range -- but it never reaches the panel `hybrid` delivers.

`HybridOptimizer.__init__` constructs a `NetworkOptimizer` and passes it both
weights. Nothing in the package or the tests ever reads
`self.network_optimizer` again. Stage 2 is `HybridOptimizer._network_refine`,
a different method on a different object, and it has no Tm or uniformity term.

So this is the Known Issue 8 class once more, in the shape
`attach_search_config` records: both ends exist and the path does not. It is
also the reason the external report's criticism, while pointing at something
real, gets the mechanism backwards -- a mixed-length panel is not penalised by
this term, because this term is not applied.

**These tests pin present behaviour, which is a defect.** Wiring the weights in
would change every delivered panel, so it is a decision rather than a
correction, and it needs the measurement recorded in
`docs/validation/selection_weights_are_inert_2026-09-21.md` first. When that
decision is taken, these tests are what will fail, and their failure is the
signal that the wiring worked.
"""

import json

import pandas as pd
import pytest

from neoswga.core.hybrid_optimizer import HybridOptimizer
from neoswga.core.position_cache import PositionCache
from neoswga.core.reaction_conditions import ReactionConditions

from tests.conftest import plasmid_example_ready

pytestmark = pytest.mark.skipif(
    not plasmid_example_ready(),
    reason="needs the generated plasmid example, which needs jellyfish",
)

EXAMPLE = "examples/plasmid_example"


@pytest.fixture(scope="module")
def pool():
    candidates = pd.read_csv(f"{EXAMPLE}/step3_df.csv")["primer"].astype(str).tolist()
    prefixes = [f"{EXAMPLE}/pcDNA"]
    return candidates, prefixes, PositionCache(prefixes, candidates)


def panel(pool, **weights):
    candidates, prefixes, cache = pool
    optimizer = HybridOptimizer(
        position_cache=cache,
        fg_prefixes=prefixes,
        fg_seq_lengths=[6157],
        bg_prefixes=[],
        bg_seq_lengths=[],
        coverage_reach=500,
        conditions=ReactionConditions(temp=30.0, polymerase="phi29"),
        reaction_temp=30.0,
        **weights,
    )
    return tuple(optimizer.optimize(candidates, final_count=6, verbose=False).primers)


# ---------------------------------------------------------------------------
# The weights are accepted and do nothing
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "weights",
    [
        {"tm_weight": 0.25},
        {"tm_weight": 1.0},
        {"uniformity_weight": 0.5},
        {"tm_weight": 0.3, "uniformity_weight": 0.2},
    ],
)
def test_a_selection_weight_does_not_change_the_hybrid_panel(pool, weights):
    """Measured, not inferred. `--application` routes these from a profile."""
    baseline = panel(pool, tm_weight=0.0, uniformity_weight=0.0)

    assert panel(pool, **weights) == baseline, (
        f"{weights} changed the delivered panel. If that is deliberate, this "
        "file and the validation record behind it are now out of date."
    )


def test_the_object_that_would_read_them_is_never_used():
    """The mechanism, asserted directly so the reason is not lost.

    `HybridOptimizer` builds a `NetworkOptimizer` and hands it both weights.
    If anything ever reads it back, the weights might start mattering, and
    this assertion is what would notice.
    """
    import ast
    import inspect

    source = inspect.getsource(HybridOptimizer)
    tree = ast.parse(source.strip())

    reads = [
        node
        for node in ast.walk(tree)
        if isinstance(node, ast.Attribute)
        and node.attr == "network_optimizer"
        and isinstance(node.ctx, ast.Load)
    ]

    assert not reads, (
        "HybridOptimizer now reads self.network_optimizer. The selection "
        "weights it was handed may have become live; re-measure before "
        "trusting any panel."
    )


def test_the_tm_term_itself_is_real_and_steeply_length_dependent():
    """Not a defence of the term, a statement of what it would do if wired.

    The Gaussian halves 3.7 C from its optimum, and one k step moves an
    effective Tm by 4-7 C. So if these weights were ever connected, a
    mixed-length panel would be penalised hard, and that is a decision to take
    deliberately rather than by reconnecting a wire.
    """
    import math

    from neoswga.core.reaction_conditions import ReactionConditions

    conditions = ReactionConditions(temp=30.0, polymerase="phi29")
    optimal = 30.0 + 5

    def tm_score(primer):
        return math.exp(-0.05 * abs(conditions.calculate_effective_tm(primer) - optimal) ** 2)

    short, mid = tm_score("ACGTACG"), tm_score("ACGTACGTA")

    assert mid > short
    assert mid / max(short, 1e-12) > 100, (mid, short)
