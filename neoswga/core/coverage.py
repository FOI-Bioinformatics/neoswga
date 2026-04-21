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
    extension: int = 3000,
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
        extension: extension reach in bp. Default 3000 bp corresponds to
            the realistic mean phi29 amplicon (Phase 16 critical gap #2).
            A position at ``p`` marks
            ``[max(0, p - extension), min(length, p + extension))``
            as occupied. Use :func:`polymerase_extension_reach` to pick
            a per-polymerase value; pass 70000 explicitly if you want
            the theoretical phi29 processivity upper bound rather than
            the practitioner-observed mean amplicon.
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


def polymerase_extension_reach(
    polymerase: str,
    default: int = 3000,
    coverage_metric: str = "realistic",
) -> int:
    """Resolve the extension reach (in bp) for a polymerase.

    Phase 16 critical gap #2: distinguish the two legitimate meanings of
    "extension reach":

    - ``coverage_metric='realistic'`` (default): returns the typical mean
      amplicon length actually observed in MDA / SWGA reactions (phi29
      ~3 kb, equiphi29 ~4 kb, bst ~1 kb, klenow ~1.5 kb). Use this for
      user-facing coverage metrics — "how much genome is amplified to
      meaningful copy number?".
    - ``coverage_metric='processivity'``: returns the theoretical
      single-molecule processivity (phi29 70 kb, equiphi29 80 kb,
      bst 2 kb, klenow 10 kb). Use this when you care about graph-level
      reachability — "can a primer at position X reach position Y in
      principle?".

    Prior to Phase 16 this helper returned processivity unconditionally,
    which inflated user-facing coverage numbers 5-20x over what a
    practitioner observes in the lab. Callers that still need the legacy
    behaviour pass ``coverage_metric='processivity'`` explicitly.

    Args:
        polymerase: Polymerase name (phi29 / equiphi29 / bst / klenow).
        default: Fallback when the polymerase is unknown or the helper
            is unavailable. Default switched to 3000 bp (realistic phi29
            amplicon) in Phase 16; pass 70000 explicitly if the legacy
            processivity value is intended.
        coverage_metric: 'realistic' (default, Phase 16+) or 'processivity'.

    Returns:
        Extension reach in bp.
    """
    try:
        if coverage_metric == "processivity":
            from .reaction_conditions import get_polymerase_processivity
            return int(get_polymerase_processivity(polymerase))
        # realistic (default)
        from .reaction_conditions import get_typical_amplicon_length
        return int(get_typical_amplicon_length(polymerase))
    except Exception:
        return default
