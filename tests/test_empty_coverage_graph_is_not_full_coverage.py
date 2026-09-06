"""An empty coverage graph is reported as empty, not as full coverage.

Finding A7. The loop breaks when covered_regions equals graph.regions, which is
vacuously true of zero regions, so a graph built from candidates with no cached
positions announced "Full coverage achieved" immediately before reporting 0.0%.
"""

import logging

from neoswga.core.dominating_set_optimizer import BipartiteGraph, DominatingSetOptimizer


class _EmptyCache:
    """A position cache that knows nothing, which is what an absent index is."""

    def get_positions(self, prefix, primer, strand):
        import numpy as np

        return np.array([], dtype=np.int64)


def test_an_empty_graph_reports_no_coverage():
    graph = BipartiteGraph(bin_size=1000)
    assert len(graph.regions) == 0
    assert graph.get_coverage_score(set()) == 0.0


def test_optimizer_warns_instead_of_claiming_full_coverage(caplog):
    # The real signature is DominatingSetOptimizer(cache, fg_prefixes,
    # fg_seq_lengths, bin_size=10000, extension_reach=0).
    optimizer = DominatingSetOptimizer(
        _EmptyCache(), fg_prefixes=["missing_prefix"], fg_seq_lengths=[100000]
    )
    with caplog.at_level(logging.WARNING):
        result = optimizer.optimize_greedy(["AAACCCGGGTTT", "ACCCGGGTTTAA"], max_primers=4)

    messages = " ".join(r.getMessage() for r in caplog.records)
    assert "Full coverage achieved" not in messages
    assert "no coverage regions" in messages
    assert result["n_primers"] == 0
