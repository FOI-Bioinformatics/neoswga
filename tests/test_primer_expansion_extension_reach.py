"""Regression tests: PrimerExpander must thread extension_reach to its optimizer.

Without this, the underlying DominatingSetOptimizer falls back to extension_reach=0,
meaning each primer covers only the single 10 kb bin where its binding site lands
and no neighbouring bins via amplicon extension. The reported coverage is then
much smaller than the actual amplification reach.
"""

from unittest.mock import patch, MagicMock

from neoswga.core.primer_expansion import PrimerExpander


def _make_expander(**kwargs):
    cache = MagicMock()
    cache.get_positions.return_value = []
    return PrimerExpander(
        position_cache=cache,
        fg_prefixes=["foo"],
        fg_seq_lengths=[1_000_000],
        **kwargs,
    )


def test_default_max_extension_matches_base_optimizer_realistic_reach():
    """PrimerExpander default reach must align with base_optimizer's measurement
    default (3000 bp, Phase 16 realistic per-primer reach), not theoretical
    processivity (70000 bp).
    """
    expander = _make_expander()
    assert expander.max_extension == 3000


def test_expand_dominating_set_passes_extension_reach_to_optimizer():
    """When PrimerExpander dispatches to DominatingSetOptimizer it must pass
    extension_reach so coverage reflects amplicon reach, not only the single
    bin containing the binding site.
    """
    expander = _make_expander(max_extension=3000)
    captured = {}

    class FakeOptimizer:
        def __init__(self, *args, **kwargs):
            captured.update(kwargs)

        def optimize_greedy(self, *args, **kwargs):
            return {"new_primers": [], "coverage": 0.0}

    with patch(
        "neoswga.core.dominating_set_optimizer.DominatingSetOptimizer",
        FakeOptimizer,
    ):
        # _expand_dominating_set is the path through which PrimerExpander
        # constructs the optimizer; call it directly with empty inputs.
        expander._expand_dominating_set(
            candidates=[],
            fixed_primers=[],
            target_new=0,
            verbose=False,
        )

    assert "extension_reach" in captured, (
        "DominatingSetOptimizer was constructed without extension_reach; "
        "it will default to 0 and report bin-only coverage."
    )
    assert captured["extension_reach"] == 3000
