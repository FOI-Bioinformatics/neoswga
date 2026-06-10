"""Strand-convention adapters.

NeoSWGA historically grew three different ways to label primer-binding
strand, and refactoring all call sites at once is high-risk. This module
provides small pure helpers so NEW code (BAM/BED gap work, exports) can
translate between conventions without inventing a fourth one.

Conventions in the codebase:

============================  =========================================
Module / context              Strand encoding
============================  =========================================
position_cache.PositionCache  ``'forward'`` / ``'reverse'`` / ``'both'``
strand_bias_analyzer          ``'+'`` / ``'-'``
amplicon_network              sign-encoded positions (negative = reverse)
gpu_acceleration              sign-encoded positions (negative = reverse)
BED / export.export_to_bed    ``'+'`` / ``'-'``
============================  =========================================

Use these helpers at the boundary between subsystems; do not refactor the
existing internal encodings.
"""

from __future__ import annotations

_CACHE_STRANDS = ("forward", "reverse", "both")


def to_cache_strand(strand: str) -> str:
    """Normalize any supported strand label to the PositionCache vocabulary.

    Accepts ``'+'``/``'-'``, ``'forward'``/``'reverse'``/``'both'``, and the
    single-letter ``'f'``/``'r'`` shorthands. Returns one of
    ``'forward'``/``'reverse'``/``'both'``.

    Raises:
        ValueError: if the label is not recognized.
    """
    s = str(strand).strip().lower()
    if s in _CACHE_STRANDS:
        return s
    if s in ("+", "fwd", "f"):
        return "forward"
    if s in ("-", "rev", "r"):
        return "reverse"
    raise ValueError(f"Unrecognized strand label: {strand!r}")


def to_bed_strand(strand: str) -> str:
    """Translate a strand label to the BED ``'+'``/``'-'`` convention.

    ``'both'`` maps to ``'.'`` (unstranded), matching BED's strand column
    for features without a defined orientation.
    """
    s = to_cache_strand(strand)
    if s == "forward":
        return "+"
    if s == "reverse":
        return "-"
    return "."  # 'both' -> unstranded


def from_sign(position: int) -> str:
    """Decode the sign-encoded position convention.

    Modules that store reverse-strand sites as negative positions
    (amplicon_network, gpu_acceleration) use this. A negative position is
    reverse; zero or positive is forward.

    Returns ``'forward'`` or ``'reverse'`` (PositionCache vocabulary).
    """
    return "reverse" if int(position) < 0 else "forward"
