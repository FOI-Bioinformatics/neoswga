"""Vendor CSV length, Tm and GC describe the oligo, not its order syntax.

`export_vendor_csv` applied the modification syntax first (phosphorothioate
bonds written as `*`) and then computed length, Tm and GC on that string, so
TGTACTGCCAA was listed as 13 nt and 38.5% GC under the default profile
(2026-11-23 sweep, tier 1 item 10).
"""

import csv

import pytest

from neoswga.core.export import (
    PrimerExporter,
    PrimerModifications,
    export_to_vendor_csv,
)

PRIMER = "TGTACTGCCAA"


def _rows(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle))


def test_vendor_csv_figures_use_the_unmodified_sequence(tmp_path):
    exporter = PrimerExporter(
        [PRIMER], modifications=PrimerModifications.from_profile("standard")
    )
    out = tmp_path / "order.csv"
    exporter.export_vendor_csv(str(out), vendor="generic")
    (row,) = _rows(out)
    assert "*" in row["sequence"]
    assert int(row["length"]) == 11
    assert row["gc"] == "45.5%"


def test_vendor_csv_tm_matches_the_unmodified_export(tmp_path):
    exporter = PrimerExporter(
        [PRIMER], modifications=PrimerModifications.from_profile("standard")
    )
    modified = tmp_path / "modified.csv"
    exporter.export_vendor_csv(str(modified), vendor="generic")
    plain = tmp_path / "plain.csv"
    export_to_vendor_csv([PRIMER], str(plain), vendor="generic", conditions=exporter.conditions())
    assert _rows(modified)[0]["tm"] == _rows(plain)[0]["tm"]


def test_ordered_sequences_must_match_the_primers(tmp_path):
    with pytest.raises(ValueError):
        export_to_vendor_csv(
            [PRIMER], str(tmp_path / "x.csv"), ordered_sequences=[PRIMER, PRIMER]
        )
