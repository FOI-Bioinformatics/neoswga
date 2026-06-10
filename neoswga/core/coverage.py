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
    circular: bool = False,
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
            the effective per-primer reach in a dense phi29 SWGA design
            (Clarke et al. 2017 post-hoc <5 kb inter-primer-site filter;
            Dwivedi-Yu et al. 2023 1/2-5 kbp successful-set densities).
            A position
            at ``p`` marks
            ``[max(0, p - extension), min(length, p + extension))``
            as occupied. Use :func:`polymerase_extension_reach` to pick
            a per-polymerase value; pass 70000 explicitly if you want
            the theoretical phi29 processivity upper bound (a different
            question — "could two primers connect in principle?").
        strand: 'both' / 'forward' / 'reverse'.
        circular: If True, windows that run off either end of the
            sequence wrap around to the other end. Default False matches
            linear chromosomes, fragmented assemblies, and concatenated
            multi-genome inputs. Set True for single closed-circular
            bacterial chromosomes or plasmids (passes `fg_circular` /
            `bg_circular` from params.json through to here).

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
        # If the window diameter equals or exceeds the genome on a
        # circular target, every position is covered by any single site.
        if circular and 2 * extension >= length:
            occupied[:] = True
        else:
            for primer in primers:
                try:
                    positions = cache.get_positions(prefix, primer, strand)
                except (KeyError, ValueError):
                    # Genuinely-absent primer in this prefix; skip it. Do NOT
                    # swallow I/O / HDF5 errors here — a systemic cache failure
                    # must surface rather than silently undercount coverage.
                    continue
                for pos in positions:
                    _mark_window(occupied, int(pos), extension, length, circular)
        covered = int(occupied.sum())
        per_prefix[prefix] = covered / length if length else 0.0
        total_cov += covered
        total_len += length

    agg = total_cov / total_len if total_len else 0.0
    return agg, per_prefix


def _mark_window(
    occupied: "np.ndarray",
    pos: int,
    extension: int,
    length: int,
    circular: bool,
) -> None:
    """Mark ``occupied[pos-extension:pos+extension]`` as True.

    On circular sequences the window wraps across the origin when it
    runs past either end. On linear sequences the window is clipped.
    Used by :func:`compute_per_prefix_coverage` and
    :meth:`base_optimizer.BaseOptimizer._compute_coverage` to keep the
    two coverage implementations in sync.
    """
    start = pos - extension
    end = pos + extension
    if circular:
        if start < 0 and end > length:
            # Window exceeds genome from both ends; everything is covered.
            occupied[:] = True
        elif start < 0:
            occupied[0:end] = True
            occupied[length + start : length] = True
        elif end > length:
            occupied[start:length] = True
            occupied[0 : end - length] = True
        else:
            occupied[start:end] = True
    else:
        clipped_start = max(0, start)
        clipped_end = min(length, end)
        if clipped_end > clipped_start:
            occupied[clipped_start:clipped_end] = True


def polymerase_extension_reach(
    polymerase: str,
    default: int = 3000,
    coverage_metric: str = "realistic",
) -> int:
    """Resolve the extension reach (in bp) for a polymerase.

    Phase 16 critical gap #2: distinguish two legitimate meanings of
    "extension reach":

    - ``coverage_metric='realistic'`` (default): returns the effective
      per-primer reach in a dense SWGA design (phi29 ~3 kb, equiphi29
      ~4 kb, bst ~1 kb, klenow ~1.5 kb). In dense multi-primer reactions
      extension is truncated by neighbouring primers' strand-displacement
      products after a few kb — Clarke et al. (2017) filter candidate
      sets to <5 kb mean inter-primer-site spacing; Dwivedi-Yu et al.
      (2023) report successful Prevotella sets at 1/2-5 kbp densities.
      Use this for `fg_coverage` / `per_target_coverage`.
    - ``coverage_metric='processivity'``: returns the theoretical
      single-molecule processivity (phi29 70 kb, equiphi29 80 kb,
      bst 2 kb, klenow 10 kb). Use this when the question is graph-level
      reachability — "can primer A and B connect via one extension event
      in principle?" — not "how much genome does the set cover?".

    Prior to Phase 16 this helper returned processivity unconditionally,
    which inflated `fg_coverage` 5-20x over what the selected primer set
    can actually amplify in a dense SWGA reaction. Callers that still
    need the legacy behaviour pass ``coverage_metric='processivity'``
    explicitly.

    Args:
        polymerase: Polymerase name (phi29 / equiphi29 / bst / klenow).
        default: Fallback when the polymerase is unknown or the helper
            is unavailable. Default switched to 3000 bp (phi29 per-primer
            reach) in Phase 16; pass 70000 explicitly if the legacy
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
