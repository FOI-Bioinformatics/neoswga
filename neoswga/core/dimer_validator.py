"""
Unified dimer validation for primer sets.

Consolidates dimer checking logic that was previously scattered across
greedy_optimizer, network_optimizer, base_optimizer, genetic_algorithm,
and clique_optimizer into a single reusable component.

Usage:
    validator = DimerValidator(max_dimer_bp=4)

    # Check a single pair
    if validator.is_compatible(primer1, primer2):
        ...

    # Build compatibility matrix for a pool
    matrix = validator.build_matrix(candidates)

    # Check if a set is dimer-free
    if validator.is_set_compatible(primer_set):
        ...

    # Get dimer risk score (fraction of pairs with dimers)
    risk = validator.dimer_risk(primer_set)
"""

import logging
import numpy as np
from typing import List, Tuple, Optional

logger = logging.getLogger(__name__)

# Pool size threshold for switching to parallel matrix computation
_PARALLEL_THRESHOLD = 200


class DimerValidator:
    """
    Centralized dimer compatibility checking with caching.

    Wraps the low-level functions in dimer.py and provides a consistent
    interface used by all optimizers.
    """

    def __init__(
        self,
        max_dimer_bp: int = 4,
        max_self_dimer_bp: int = 4,
    ):
        self.max_dimer_bp = max_dimer_bp
        self.max_self_dimer_bp = max_self_dimer_bp
        self._pair_cache = {}

    def is_compatible(self, primer1: str, primer2: str) -> bool:
        """Check whether two primers are dimer-free.

        Results are cached for repeated lookups on the same pair.
        """
        key = (primer1, primer2) if primer1 <= primer2 else (primer2, primer1)
        if key not in self._pair_cache:
            from . import dimer
            self._pair_cache[key] = not dimer.is_dimer_fast(
                primer1, primer2, self.max_dimer_bp
            )
        return self._pair_cache[key]

    def has_self_dimer(self, primer: str) -> bool:
        """Check whether a primer forms a homodimer."""
        from . import dimer
        return dimer.is_dimer_fast(primer, primer, self.max_self_dimer_bp)

    def build_matrix(self, primers: List[str]) -> np.ndarray:
        """Build an n x n binary dimer matrix.

        Entry (i, j) is 1 if primers[i] and primers[j] form a dimer,
        0 otherwise.  Automatically uses parallel computation for pools
        larger than the threshold.

        Args:
            primers: List of primer sequences.

        Returns:
            numpy array of shape (len(primers), len(primers)).
        """
        from . import dimer

        if len(primers) > _PARALLEL_THRESHOLD:
            matrix = dimer.heterodimer_matrix_parallel(
                primers, max_dimer_bp=self.max_dimer_bp
            )
        else:
            matrix = dimer.heterodimer_matrix_fast(
                primers, max_dimer_bp=self.max_dimer_bp
            )

        return matrix.astype(bool)

    def is_set_compatible(self, primers: List[str]) -> bool:
        """Return True if no pair in the set forms a dimer."""
        for i in range(len(primers)):
            for j in range(i + 1, len(primers)):
                if not self.is_compatible(primers[i], primers[j]):
                    return False
        return True

    def dimer_risk(self, primers: List[str]) -> float:
        """Fraction of primer pairs that form dimers (0.0 = dimer-free)."""
        n = len(primers)
        if n < 2:
            return 0.0
        total = 0
        dimers = 0
        for i in range(n):
            for j in range(i + 1, n):
                total += 1
                if not self.is_compatible(primers[i], primers[j]):
                    dimers += 1
        return dimers / total

    def incompatible_pairs(
        self, primers: List[str]
    ) -> List[Tuple[str, str]]:
        """Return list of primer pairs that form dimers."""
        pairs = []
        for i in range(len(primers)):
            for j in range(i + 1, len(primers)):
                if not self.is_compatible(primers[i], primers[j]):
                    pairs.append((primers[i], primers[j]))
        return pairs

    def filter_self_dimers(self, primers: List[str]) -> List[str]:
        """Remove primers that form homodimers."""
        return [p for p in primers if not self.has_self_dimer(p)]

    def clear_cache(self):
        """Clear the pair compatibility cache."""
        self._pair_cache.clear()
