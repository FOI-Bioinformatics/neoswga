"""Would a hard gap constraint separate the wet-lab winners? No.

Tests the idea that a coverage hole wider than twice the polymerase's reach
makes a panel defective regardless of which property limits the design. If
winners never carry such a hole and losers sometimes do, a hard constraint is
supported without fitting anything to the outcomes.

Measured against the 18 published sets with wet-lab outcomes already in
`tests/validation/data/`: 12 from Clarke et al. (2017) against
*M. tuberculosis* and 6 from Dwivedi-Yu et al. (2023) against *Prevotella*.

The answer is that every reach-derived threshold rejects all 18, winners
included, so the measurement is evidence AGAINST the constraint. It is more
useful as a description of what a working panel looks like: panels enriching 96
to 120 fold carry 31 to 33 kb holes at a 3 kb selection reach.

`fg_max_distance` is each paper's own published figure under its own
definition, not a quantity recomputed from the genomes.

Usage:
    python scripts/benchmarking/published_gap_thresholds.py
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Callable, Dict, List, Tuple

DATA = Path(__file__).resolve().parents[2] / "tests" / "validation" / "data"

# Reaches to test. The calibrated band is 3.0 to 6.2 kb
# (docs/validation/reach_calibration.md); the rest are included to show where a
# separating threshold would have to sit and how far outside that band it is.
REACHES = (3_000, 4_000, 5_000, 6_200, 10_000, 20_000, 70_000)

Row = Tuple[str, float, int, bool]


def _load(filename: str, is_winner: Callable[[Dict], bool]) -> List[Row]:
    payload = json.loads((DATA / filename).read_text())
    rows = []
    for name, entry in payload["sets"].items():
        published = entry["published"]
        rows.append(
            (
                name,
                float(published["fg_max_distance"]),
                int(published["size"]),
                bool(is_winner(entry["outcome"])),
            )
        )
    return rows


def benchmarks() -> List[Tuple[str, List[Row]]]:
    """The two fixtures, with each paper's own definition of success."""
    return [
        (
            "Clarke M. tuberculosis",
            _load(
                "clarke_2017_mtb.json",
                lambda outcome: outcome.get("most_effective") is True,
            ),
        ),
        (
            "Dwivedi-Yu Prevotella",
            _load(
                "dwivedi_yu_2023_prevotella.json",
                lambda outcome: outcome.get("successful") is True,
            ),
        ),
    ]


def describe(label: str, rows: List[Row]) -> None:
    print(f"\n## {label} (n={len(rows)})")
    print(f"{'set':12s} {'primers':>8s} {'worst hole bp':>14s} {'winner':>7s}")
    for name, gap, size, winner in sorted(rows, key=lambda row: row[1]):
        print(f"{name:12s} {size:8d} {gap:14.0f} {'yes' if winner else '':>7s}")


def test_thresholds(rows: List[Row]) -> None:
    winners = [row for row in rows if row[3]]
    others = [row for row in rows if not row[3]]
    print("\n## Does a reach-derived threshold separate winners from the rest?")
    print(
        f"{'reach':>7s} {'2x reach':>9s} {'winners above':>14s} "
        f"{'others above':>13s} {'separates':>10s}"
    )
    for reach in REACHES:
        threshold = 2 * reach
        above_w = sum(1 for row in winners if row[1] > threshold)
        above_o = sum(1 for row in others if row[1] > threshold)
        separates = "yes" if above_w == 0 and above_o > 0 else "no"
        print(
            f"{reach:7d} {threshold:9d} {f'{above_w} of {len(winners)}':>14s} "
            f"{f'{above_o} of {len(others)}':>13s} {separates:>10s}"
        )


def main() -> None:
    all_rows: List[Row] = []
    for label, rows in benchmarks():
        describe(label, rows)
        all_rows.extend(rows)
    test_thresholds(all_rows)

    tightest = min(all_rows, key=lambda row: row[1])
    print(
        f"\nTightest of all {len(all_rows)} published sets: {tightest[0]} at "
        f"{tightest[1]:.0f} bp, winner={'yes' if tightest[3] else 'no'}."
    )
    print(
        "Every set exceeds twice the calibrated reach, so hole-closing is not "
        "what makes a published SWGA panel work."
    )


if __name__ == "__main__":
    main()
