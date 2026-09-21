"""How much of the time each oligo length is bound, in one reaction.

A design may mix oligo lengths: `params.schema.json` allows k of 4 to 30, the
scan writes one position index per length, and
`tests/integration/test_variable_oligo_length.py` takes a mixed design through
all four steps and out to an ordering file. That is a mechanical fact and it is
settled.

The chemistry is not so simple, and this module exists to say so rather than
let it be assumed. Occupancy is the fraction of the time a site is bound, and
it rises very steeply with length at a fixed temperature. Measured on 3,000
random k-mers per length, median occupancy:

    k         7       8       9      10      11      12
    phi29 30 C   0.017   0.145   0.661   0.960   0.996   1.000
    equiphi29 42 C 0.001  0.005   0.032   0.205   0.611   0.917

So a phi29 panel spanning k 7 to 11 holds oligos bound under 2% of the time
beside oligos bound essentially always. They are not equal contributors, and a
panel size counts them as though they were.

Length is not the only driver: on the bundled plasmid panel two 8-mers sit at
0.066 and 0.693, which is composition rather than length. What the per-length
view adds is that the pool-wide mean `filter` reports describes no length in a
mixed pool at all -- on that pool it is 2.405 where the lengths run 3.010 down
to 1.015.

**Reported, never enforced**, which is the resolution Known Issue 17 reached
for the same quantity. Gating candidates on occupancy was measured on the
Wolbachia pool and made the delivered panel worse on both axes. The lever that
works is the reaction, not the filter, so this names the temperature and says
what a low-occupancy member means; it removes nothing.

Silent on a single-length design, which is every design in this repository
before 2026-09-21, so no existing run acquires new output.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

logger = logging.getLogger(__name__)

__all__ = [
    "LengthOccupancy",
    "occupancy_by_length",
    "log_occupancy_spread",
]

#: Below this, an oligo spends most of the reaction unbound. Not a threshold
#: anything is rejected for: it is where the message changes from describing a
#: spread to naming a member that contributes little.
WEAK_OCCUPANCY = 0.1


@dataclass(frozen=True)
class LengthOccupancy:
    """Occupancy for the oligos of one length, and how many could not be read.

    `unmeasured` is carried rather than dropped. Silently skipping a primer
    whose thermodynamics would not evaluate leaves a mean over an unknown
    subset reported with the authority of a mean over the whole group, which is
    the failure `occupancy.discrimination_profile` records at its own except
    clause.
    """

    length: int
    n: int
    median: Optional[float]
    lowest: Optional[float]
    highest: Optional[float]
    unmeasured: int


def occupancy_by_length(primers: Sequence[str], conditions) -> Tuple[LengthOccupancy, ...]:
    """Per-length occupancy for `primers` under `conditions`, shortest first.

    Returns an empty tuple when there is no temperature to evaluate at. An
    empty answer says "not measured"; a zero would say "never bound", and the
    two must not share a representation.
    """
    if conditions is None or not primers:
        return ()
    temp = getattr(conditions, "temp", None)
    if temp is None:
        return ()

    from neoswga.core.exceptions import DesignError
    from neoswga.core.occupancy import site_occupancy
    from neoswga.core.thermodynamics import calculate_enthalpy_entropy

    grouped: Dict[int, List[float]] = {}
    unmeasured: Dict[int, int] = {}
    for primer in primers:
        sequence = str(primer)
        length = len(sequence)
        grouped.setdefault(length, [])
        unmeasured.setdefault(length, 0)
        try:
            enthalpy, _ = calculate_enthalpy_entropy(sequence)
            tm = conditions.calculate_effective_tm(sequence)
        except DesignError:
            raise
        except Exception:
            unmeasured[length] += 1
            continue
        grouped[length].append(site_occupancy(enthalpy, tm, temp))

    records = []
    for length in sorted(grouped):
        values = sorted(grouped[length])
        records.append(
            LengthOccupancy(
                length=length,
                n=len(values) + unmeasured[length],
                median=_median(values),
                lowest=values[0] if values else None,
                highest=values[-1] if values else None,
                unmeasured=unmeasured[length],
            )
        )
    return tuple(records)


def log_occupancy_spread(primers: Sequence[str], conditions, label: str = "pool") -> None:
    """Say what a mixed-length design is asking of one reaction.

    Silent unless more than one length is present, so a single-length design
    gains no output. Best-effort: a diagnostic must not fail the step that
    produced what it describes.
    """
    try:
        records = occupancy_by_length(primers, conditions)
    except Exception as exc:  # pragma: no cover - diagnostic only
        logger.debug("Could not profile occupancy by length: %s", exc)
        return

    if len(records) < 2:
        return

    temp = getattr(conditions, "temp", None)
    logger.info(
        "Oligo lengths in this %s span %d to %d. Occupancy -- the fraction of "
        "the time a site is bound -- rises steeply with length at a fixed "
        "temperature, so the members are not equal contributors:",
        label,
        records[0].length,
        records[-1].length,
    )
    for record in records:
        if record.median is None:
            logger.info("    k=%-3d n=%-4d occupancy unmeasurable", record.length, record.n)
            continue
        note = f", {record.unmeasured} unmeasured" if record.unmeasured else ""
        logger.info(
            "    k=%-3d n=%-4d median occupancy %.3f (%.3f to %.3f)%s",
            record.length,
            record.n,
            record.median,
            record.lowest,
            record.highest,
            note,
        )

    measured = [record for record in records if record.median is not None]
    if not measured:
        return
    weakest = min(measured, key=lambda record: record.median)
    strongest = max(measured, key=lambda record: record.median)
    if weakest.median < WEAK_OCCUPANCY <= strongest.median:
        logger.warning(
            "At %s the %d-mers in this %s are bound about %.1f%% of the time "
            "against %.1f%% for the %d-mers. The short members contribute "
            "little at this temperature. Nothing is removed for it: the lever "
            "is the reaction, so lower the temperature or drop the shortest "
            "length if you meant every member to bind.",
            f"{temp} C" if temp is not None else "this temperature",
            weakest.length,
            label,
            weakest.median * 100,
            strongest.median * 100,
            strongest.length,
        )


def _median(values: List[float]) -> Optional[float]:
    if not values:
        return None
    middle = len(values) // 2
    if len(values) % 2:
        return values[middle]
    return (values[middle - 1] + values[middle]) / 2
