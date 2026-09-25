"""Which of the alternative primer sets a command is talking about.

`step4_improved_df.csv` holds up to `max_sets` (default 5) primer sets, one per
`set_index`, and they are ALTERNATIVES: each is found by excluding the primers
already chosen and selecting again, so set 0 and set 1 are two answers to the
same question rather than two halves of one answer. Set 0 is the one
`step4_improved_df_summary.json` describes.

Four commands read that file and, until 2026-09-21, none of them looked at
`set_index`. They analysed, reported, simulated and EXPORTED the union of every
set as though it were one panel.

Measured on the bundled plasmid example at 300 bp reach, where the pool is
large enough for alternatives to be found:

- set 0 is 8 oligos with no pair above the configured `max_dimer_bp` of 3;
- the exported FASTA is 18 oligos with five pairs above it, the worst a 9 bp
  complementary run, and every one of the five joins oligos from DIFFERENT
  sets.

Those cross-set pairs have never been screened and by construction never could
be: two alternatives are not meant to share a tube. The delivered panel honours
the dimer limit, the file offered "for ordering" does not, and nothing in it
marks where one set ends and the next begins -- the records are numbered
SWGA_001 upward straight through.

No saved run in this repository exhibits it, because every one holds set 0
alone, which is why it survived. It needs only a pool big enough for
`collect_alternative_sets` to find a second set, which is the default
configuration.

This module holds the selection rule once so the four readers cannot disagree.
"""

from __future__ import annotations

import csv
import logging
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path

from neoswga.core.exceptions import ReferenceDataError

logger = logging.getLogger(__name__)

__all__ = ["DeliveredSet", "SET_INDEX_COLUMN", "read_delivered_set", "select_delivered_set"]

SET_INDEX_COLUMN = "set_index"

#: The set `step4_improved_df_summary.json` describes.
DEFAULT_SET_INDEX = 0


@dataclass(frozen=True)
class DeliveredSet:
    """One primer set, and what else the file held.

    `set_index` is None when the file carries no `set_index` column at all,
    which is a fact about the file rather than a set numbered None. Results
    written before the column existed are that shape, and they hold exactly one
    set, so every row is returned.
    """

    rows: tuple[dict[str, str], ...]
    set_index: int | None
    available: tuple[int, ...]

    @property
    def primers(self) -> tuple[str, ...]:
        """The sequences, in file order, deduplicated.

        A set should not repeat a primer; deduplicating here keeps a
        malformed file from inflating a panel size rather than hiding it,
        and the count is reported by `describe`.
        """
        seen: list[str] = []
        for row in self.rows:
            sequence = str(row.get("primer") or row.get("sequence") or "").strip()
            if sequence and sequence not in seen:
                seen.append(sequence)
        return tuple(seen)

    @property
    def has_alternatives(self) -> bool:
        return len(self.available) > 1

    def describe(self) -> str:
        """One line naming what was selected, for a log or a header."""
        if self.set_index is None:
            return f"{len(self.primers)} primers (file records a single set)"
        others = [index for index in self.available if index != self.set_index]
        suffix = f"; {len(others)} alternative set(s) not included: {others}" if others else ""
        return f"set {self.set_index}, {len(self.primers)} primers{suffix}"


def select_delivered_set(
    rows: Sequence[dict[str, str]],
    set_index: int | None = DEFAULT_SET_INDEX,
    *,
    artifact: str = "step4_improved_df.csv",
) -> DeliveredSet:
    """Pick one set out of already-loaded result rows.

    `set_index=None` asks for every row regardless of set. That is the one way
    to get the old behaviour, and a caller must name it, because pooling
    alternatives is a request nobody should make by accident.

    A requested set that is not in the file raises rather than returning an
    empty list. An empty panel and an absent one read identically downstream,
    which is the distinction `ReferenceDataError` exists to keep.
    """
    rows = list(rows)
    if not rows:
        return DeliveredSet(rows=(), set_index=None, available=())

    if SET_INDEX_COLUMN not in rows[0]:
        return DeliveredSet(rows=tuple(rows), set_index=None, available=())

    available = []
    for row in rows:
        try:
            value = int(float(row[SET_INDEX_COLUMN]))
        except (TypeError, ValueError):
            continue
        if value not in available:
            available.append(value)
    available.sort()

    if set_index is None:
        return DeliveredSet(rows=tuple(rows), set_index=None, available=tuple(available))

    if set_index not in available:
        raise ReferenceDataError(
            artifact=artifact,
            reason=(
                f"no primer set with set_index={set_index}; the file holds "
                f"{available or 'none'}"
            ),
            remediation=(
                "Pass a set index the file holds, or re-run `neoswga optimize` "
                "with a larger `max_sets` if you expected more alternatives."
            ),
        )

    selected = tuple(row for row in rows if _as_int(row.get(SET_INDEX_COLUMN)) == set_index)
    return DeliveredSet(rows=selected, set_index=set_index, available=tuple(available))


def read_delivered_set(
    path,
    set_index: int | None = DEFAULT_SET_INDEX,
) -> DeliveredSet:
    """Read `step4_improved_df.csv` and return one set from it."""
    path = Path(path)
    if not path.exists():
        raise ReferenceDataError(
            artifact=str(path),
            reason="results file not found",
            remediation="Run `neoswga optimize` first.",
        )
    with open(path, newline="") as handle:
        rows = list(csv.DictReader(handle))

    delivered = select_delivered_set(rows, set_index, artifact=str(path))
    if delivered.has_alternatives:
        logger.info(
            "Reading %s from %s. Alternative sets are separate answers to the "
            "same design, not additions to it: pooling them puts oligos in one "
            "tube that no dimer screen has ever compared.",
            delivered.describe(),
            path.name,
        )
    return delivered


def _as_int(value) -> int | None:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None
