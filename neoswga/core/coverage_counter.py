"""Per-bin coverage counts that survive a removal without a rebuild.

`HybridOptimizer._calculate_coverage` rebuilds a `BipartiteGraph` over every
primer in the set on every call, and `_prune_background` calls it once per
candidate per removal. Pruning 160 primers to 24 is 136 removals, each
evaluating about 90 candidates, each rebuilding over about 90 primers.

The arithmetic does not need any of that. Hold how many primers cover each bin;
removing a primer decrements only the bins it occupies, and the coverage it was
carrying alone is the number of its bins whose count is one. This is exact, not
an estimate: it computes the same number the rebuild computes.
"""

from typing import Dict, Hashable, Iterable, Set


class CoverageCounter:
    """How many of the current primers cover each bin."""

    def __init__(self, total_bins: int):
        self.total_bins = int(total_bins)
        self._count: Dict[Hashable, int] = {}
        self._bins_of: Dict[str, Set[Hashable]] = {}

    def add(self, primer: str, bins: Iterable[Hashable]) -> None:
        """Record that `primer` covers `bins`. Adding a known primer replaces it."""
        if primer in self._bins_of:
            self.remove(primer)
        owned = set(bins)
        self._bins_of[primer] = owned
        for b in owned:
            self._count[b] = self._count.get(b, 0) + 1

    def remove(self, primer: str) -> None:
        """Drop `primer`. Unknown primers are ignored."""
        owned = self._bins_of.pop(primer, None)
        if owned is None:
            return
        for b in owned:
            remaining = self._count.get(b, 0) - 1
            if remaining > 0:
                self._count[b] = remaining
            else:
                self._count.pop(b, None)

    def covered_count(self) -> int:
        """How many bins at least one current primer covers."""
        return len(self._count)

    def covered_fraction(self) -> float:
        """Covered bins as a fraction of the genome's bins.

        Zero when nothing is covered. An empty counter is not full coverage; see
        finding A7 for what that confusion cost elsewhere.
        """
        if self.total_bins <= 0:
            return 0.0
        return self.covered_count() / self.total_bins

    def loss_if_removed(self, primer: str) -> int:
        """Bins that would become uncovered if `primer` were dropped."""
        owned = self._bins_of.get(primer)
        if not owned:
            return 0
        return sum(1 for b in owned if self._count.get(b, 0) == 1)
