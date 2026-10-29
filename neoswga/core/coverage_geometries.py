"""Both coverage geometries for one delivered panel, labelled.

`fg_coverage` is computed under one geometry and reported as a bare number.
Measured on 29 October 2026, the shipped `symmetric` window overstates the
directional figure by 1.27x to 1.85x on three real panels at the shipped 3 kb
reach -- 21.0, 25.8 and 19.5 percentage points
(`docs/validation/2026-10-29-directional-coverage-measured.md`). The ratio
tends to 2.0 as windows stop overlapping and to 1.0 as the genome saturates,
and the worst ABSOLUTE error sits in the middle, which is where a real SWGA
design lives.

The geometry was already settable and the reach already refitted under the
second one (4.4-6.7 kb against 2.9-4.6 kb,
`docs/validation/2026-09-23-reach-refit-directional.md`), so what was missing
was not a capability. It was that a reader of a summary could not tell which
convention produced the number in front of them, or by how much the other one
would differ. This computes both for the DELIVERED panel and labels which one
the run was scored on.

**Why not in the search.** `_compute_metrics` runs per candidate panel, and
`_oriented_positions` fetches the two orientations only when the configured
geometry asks for them, because `get_positions` is the hottest method in the
package. Computing both geometries for every candidate would pay that cost on
every evaluation to produce a figure nothing selects on. This runs once, after
the search, on the panel that is actually delivered.

**The default does not change here.** Every figure recorded in
`docs/validation/` was produced under `symmetric`, the reach was fitted
alongside it, and a directional run now has to name its own reach
(`design_request._resolve_reach`). Switching the default would move every
recorded figure in this repository and is a decision for the maintainer, not a
repair. Reporting both is what the measurement's own record proposed as the
honest intermediate step.
"""

from __future__ import annotations

import logging
from collections.abc import Sequence
from typing import Any

logger = logging.getLogger(__name__)

#: The geometry a figure is labelled with when nothing says otherwise.
DEFAULT_GEOMETRY = "symmetric"


def _positions(cache: Any, prefix: str, primers: Sequence[str]) -> tuple[list[int], list[int]]:
    """Forward and reverse-complement sites for a panel in one sequence.

    `forward` is each oligo's own occurrences and `reverse` its reverse
    complement's, which is `position_cache`'s convention and the split the
    directional model needs. An oligo the cache cannot answer for contributes
    nothing to either list; it is not counted as binding nowhere, because this
    figure is reported beside one the search already produced and a quieter
    number here than there would be worse than no number.
    """
    forward: set[int] = set()
    reverse: set[int] = set()
    for primer in primers:
        for strand, sink in (("forward", forward), ("reverse", reverse)):
            try:
                sink.update(int(p) for p in cache.get_positions(prefix, primer, strand))
            except Exception as exc:
                logger.debug(f"No {strand} positions for {primer} in {prefix}: {exc}")
    return sorted(forward), sorted(reverse)


def _coverage(
    cache: Any,
    prefixes: Sequence[str],
    lengths: Sequence[int],
    primers: Sequence[str],
    reach: int,
    circular: bool,
    geometry: str,
) -> float | None:
    """Length-weighted coverage across the targets, or `None` if unmeasurable.

    Each sequence is measured in its own coordinate space and the covered
    BASES are summed, which is the only aggregate that means anything when the
    parts are separate molecules -- the same rule
    `base_optimizer._coverage_over_prefixes` follows.
    """
    from .coverage import union_coverage

    total = sum(int(n) for n in lengths or ())
    if not prefixes or total <= 0:
        return None
    covered = 0.0
    for prefix, length in zip(prefixes, lengths, strict=True):
        if length <= 0:
            continue
        forward, reverse = _positions(cache, prefix, primers)
        if not forward and not reverse:
            continue
        covered += (
            union_coverage(
                forward,
                int(length),
                int(reach),
                circular,
                reverse=reverse,
                geometry=geometry,
            )
            * length
        )
    return covered / total


def compare_geometries(
    cache: Any,
    prefixes: Sequence[str],
    lengths: Sequence[int],
    primers: Sequence[str],
    *,
    reach: int,
    circular: bool = False,
    reported: str = DEFAULT_GEOMETRY,
) -> dict[str, Any] | None:
    """Both figures for this panel, or `None` when neither can be computed.

    `inflation` is symmetric / directional, the factor by which the wider
    window overstates the narrower one on THIS panel. `None` when the
    directional figure is zero, because a ratio over zero is not a large
    number, it is no number.

    The panel was SELECTED under `reported`, so the other column bounds how
    much that figure overstates this panel rather than how much coverage a run
    under the other geometry would deliver: a set cover that knew each site
    reaches `r` and not `2r` would choose differently.
    `docs/validation/2026-09-23-does-the-panel-move.md` measured that it does.
    """
    if not primers or cache is None:
        return None
    figures: dict[str, float | None] = {}
    for geometry in ("symmetric", "directional"):
        figures[geometry] = _coverage(cache, prefixes, lengths, primers, reach, circular, geometry)
    if figures["symmetric"] is None and figures["directional"] is None:
        return None
    wide, narrow = figures["symmetric"], figures["directional"]
    inflation = None
    if wide is not None and narrow is not None and narrow > 0:
        inflation = wide / narrow
    return {
        "reach": int(reach),
        "reported_geometry": reported,
        "symmetric": wide,
        "directional": narrow,
        "inflation": inflation,
        "note": (
            "The delivered panel was selected and scored under "
            f"'{reported}'. The other figure is the same sites under the other "
            "window, not the coverage a run under that geometry would deliver. "
            "The polymerase reach defaults were fitted under 'symmetric'."
        ),
    }


def log_comparison(comparison: dict[str, Any] | None) -> None:
    """One line, as the hybrid run already prints two coverage figures."""
    if not comparison:
        return
    wide, narrow = comparison.get("symmetric"), comparison.get("directional")
    if wide is None or narrow is None:
        return
    inflation = comparison.get("inflation")
    ratio = "n/a" if inflation is None else f"{inflation:.2f}x"
    logger.info(
        "Coverage geometry at %s bp: %.4f (symmetric) against %.4f (directional), %s. "
        "Scored on '%s'.",
        f"{comparison['reach']:,}",
        wide,
        narrow,
        ratio,
        comparison.get("reported_geometry", DEFAULT_GEOMETRY),
    )


def for_delivered_panel(parameter: Any, result: Any) -> tuple[Any, dict[str, Any] | None]:
    """The cache over the delivered panel, and both figures for it.

    Returns the cache as well because `optimize_step4` hands one back to its
    caller and building it twice would load the same positions twice. Reads
    the references and the geometry off the resolved `parameter` module and the
    reach off the result's own metrics, so the figures are computed at the
    reach the run actually scored on rather than at a default resolved again
    here.
    """
    from .position_cache import PositionCache

    prefixes = getattr(parameter, "fg_prefixes", []) or []
    primers = list(getattr(result, "primers", ()) or ())
    cache = PositionCache(prefixes, primers) if prefixes else None
    comparison = compare_geometries(
        cache,
        prefixes,
        getattr(parameter, "fg_seq_lengths", []) or [],
        primers,
        reach=getattr(result.metrics, "extension_reach", 0) or 0,
        circular=bool(getattr(parameter, "fg_circular", False)),
        reported=str(getattr(parameter, "coverage_geometry", DEFAULT_GEOMETRY) or DEFAULT_GEOMETRY),
    )
    log_comparison(comparison)
    return cache, comparison
