"""Centralized coverage computation for optimizer post-processing and
ad-hoc rescoring CLIs.

The single `compute_per_prefix_coverage` helper consolidates the numpy-
vectorised union-of-extension-windows coverage calculation that previously
lived inline in `cli_unified.run_contract_set` and `run_rescore_set`
(Phase 12D). Extracted here so `unified_optimizer` (Phase 15A) and the
CLIs share the same code path.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np


def compute_per_prefix_coverage(
    cache,
    primers: Sequence[str],
    prefixes: Sequence[str],
    seq_lengths: Sequence[int],
    extension: int = 70000,
    strand: str = "both",
) -> Tuple[float, Dict[str, float]]:
    """Compute union-of-extension-windows coverage across prefixes.

    For each (primer, prefix) pair, fetch positions from the PositionCache
    and mark occupied[start:end] = True on a bool numpy array; coverage
    fraction is the sum of the occupied array divided by the genome length.

    Args:
        cache: PositionCache instance that responds to
            ``get_positions(prefix, primer, strand)``.
        primers: primer sequences to query.
        prefixes: HDF5 file prefixes (usually fg_prefixes or bg_prefixes).
        seq_lengths: genome lengths matching `prefixes` elementwise.
        extension: polymerase extension reach in bp; a position at `p`
            marks `[max(0, p - extension), min(length, p + extension))`
            as occupied.
        strand: 'both' / 'forward' / 'reverse'.

    Returns:
        (aggregate_coverage, per_prefix_coverage_dict). aggregate is the
        ratio of total-bases-covered to total-bases across every prefix;
        per_prefix maps each prefix to its own coverage fraction.

        Empty inputs or an unusable cache yield ``(0.0, {})``.
    """
    if not prefixes or not seq_lengths or cache is None or not primers:
        return 0.0, {}

    per_prefix: Dict[str, float] = {}
    total_cov = 0
    total_len = 0

    for prefix, length in zip(prefixes, seq_lengths):
        if length <= 0:
            per_prefix[prefix] = 0.0
            continue
        occupied = np.zeros(length, dtype=bool)
        for primer in primers:
            try:
                positions = cache.get_positions(prefix, primer, strand)
            except Exception:
                continue
            for pos in positions:
                start = max(0, int(pos) - extension)
                end = min(length, int(pos) + extension)
                if end > start:
                    occupied[start:end] = True
        covered = int(occupied.sum())
        per_prefix[prefix] = covered / length if length else 0.0
        total_cov += covered
        total_len += length

    agg = total_cov / total_len if total_len else 0.0
    return agg, per_prefix


def polymerase_extension_reach(polymerase: str, default: int = 70000) -> int:
    """Resolve the extension reach (in bp) for a polymerase, with a safe
    default when the reaction_conditions helper is unavailable or the
    polymerase name is unknown.
    """
    try:
        from .reaction_conditions import get_polymerase_processivity
        return int(get_polymerase_processivity(polymerase))
    except Exception:
        return default
