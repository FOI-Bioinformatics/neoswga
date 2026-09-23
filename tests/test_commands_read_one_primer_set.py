"""The four readers of step4_improved_df.csv must each mean ONE set.

`max_sets` defaults to 5, and `collect_alternative_sets` finds each one by
excluding the primers already chosen and selecting again. So the rows of
`step4_improved_df.csv` are alternatives: set 0 and set 1 are two answers to
the same design, not two halves of one.

Until 2026-09-21, `export`, `interpret`, `report` and `simulate` all read every
row. The costly one is `export`, whose output the tool calls "Primers ready for
ordering!". Measured on the bundled plasmid example at 300 bp reach:

    set 0             :  8 oligos, 0 pairs above max_dimer_bp=3
    exported FASTA    : 18 oligos, 5 pairs above it, worst a 9 bp run

and all five cross a set boundary. Those pairs were never screened, and by
construction never could be, because two alternatives are not meant to share a
tube. The records are numbered SWGA_001 upward with nothing marking where one
set ends.

No saved run in this repository shows it: every one holds set 0 alone. It needs
only a pool large enough for a second set to be found, which is the default
configuration, which is why a unit test rather than a fixture is what holds
this closed.

`tests/test_delivered_panel_honours_the_dimer_limit.py` pins the limit on what
selection DELIVERS. This file pins it on what the commands hand to a user, and
the two were not the same thing.
"""

import csv

import pytest

from neoswga.core.delivered_set import (
    DeliveredSet,
    read_delivered_set,
    select_delivered_set,
)
from neoswga.core.exceptions import ReferenceDataError

# Three sets, deliberately overlapping in nothing, so a pooled read is visible
# in the primer count alone.
ROWS = [
    {"primer": "ACCGCATC", "set_index": "0", "score": "1.0"},
    {"primer": "TCAGCGA", "set_index": "0", "score": "1.0"},
    {"primer": "AGCACTG", "set_index": "1", "score": "0.8"},
    {"primer": "ACGTCAATGG", "set_index": "1", "score": "0.8"},
    {"primer": "AGTGAGC", "set_index": "2", "score": "0.5"},
]


def write_results(tmp_path, rows, columns=("primer", "set_index", "score")):
    path = tmp_path / "step4_improved_df.csv"
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(columns))
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in columns})
    return path


# ---------------------------------------------------------------------------
# The rule itself
# ---------------------------------------------------------------------------


def test_the_default_is_the_set_the_summary_describes():
    """Set 0, because that is the one step4_improved_df_summary.json reports."""
    delivered = select_delivered_set(ROWS)

    assert delivered.set_index == 0
    assert delivered.primers == ("ACCGCATC", "TCAGCGA")


def test_the_other_sets_are_reachable_by_index():
    assert select_delivered_set(ROWS, 1).primers == ("AGCACTG", "ACGTCAATGG")
    assert select_delivered_set(ROWS, 2).primers == ("AGTGAGC",)


def test_every_set_is_reported_even_though_one_is_returned():
    """A caller must be able to say what it left out."""
    delivered = select_delivered_set(ROWS)

    assert delivered.available == (0, 1, 2)
    assert delivered.has_alternatives
    assert "alternative" in delivered.describe()


def test_a_set_that_is_not_there_raises_rather_than_returning_nothing():
    """An empty panel and an absent one read identically downstream.

    This is the distinction `ReferenceDataError` exists to keep, and returning
    `[]` for a set nobody wrote is the silent-zero shape in its cheapest form.
    """
    with pytest.raises(ReferenceDataError) as excinfo:
        select_delivered_set(ROWS, 7)

    assert "7" in str(excinfo.value)
    assert "[0, 1, 2]" in str(excinfo.value)


def test_pooling_every_set_must_be_asked_for_by_name():
    """The old behaviour is still reachable, and it cannot happen by default."""
    pooled = select_delivered_set(ROWS, None)

    assert len(pooled.primers) == 5
    assert pooled.set_index is None


# ---------------------------------------------------------------------------
# Results written before the column existed
# ---------------------------------------------------------------------------


def test_a_file_without_the_column_returns_every_row():
    """Older directories hold exactly one set and no `set_index`."""
    rows = [{"primer": "ACCGCATC"}, {"primer": "TCAGCGA"}]

    delivered = select_delivered_set(rows)

    assert delivered.primers == ("ACCGCATC", "TCAGCGA")
    assert delivered.set_index is None, "absence of the column is not set None"
    assert delivered.available == ()
    assert not delivered.has_alternatives


def test_an_empty_file_is_empty_rather_than_an_error():
    assert select_delivered_set([]) == DeliveredSet(rows=(), set_index=None, available=())


# ---------------------------------------------------------------------------
# Reading from disk
# ---------------------------------------------------------------------------


def test_reading_a_results_file_selects_one_set(tmp_path):
    path = write_results(tmp_path, ROWS)

    assert read_delivered_set(path).primers == ("ACCGCATC", "TCAGCGA")


def test_a_missing_results_file_names_the_remedy(tmp_path):
    with pytest.raises(ReferenceDataError) as excinfo:
        read_delivered_set(tmp_path / "step4_improved_df.csv")

    assert "optimize" in str(excinfo.value)


# ---------------------------------------------------------------------------
# The commands
# ---------------------------------------------------------------------------


def test_the_exporter_ships_one_set(tmp_path):
    """The defect this file exists for, at the site where it costs money."""
    write_results(tmp_path, ROWS)

    from neoswga.core.export import PrimerExporter

    exporter = PrimerExporter.from_results_dir(str(tmp_path))

    assert list(exporter.primers) == ["ACCGCATC", "TCAGCGA"], (
        "the exported panel is not set 0; a user ordering this file gets "
        "oligos from mutually exclusive alternatives in one tube"
    )


def test_the_exporter_can_be_pointed_at_an_alternative(tmp_path):
    write_results(tmp_path, ROWS)

    from neoswga.core.export import PrimerExporter

    exporter = PrimerExporter.from_results_dir(str(tmp_path), set_index=1)

    assert list(exporter.primers) == ["AGCACTG", "ACGTCAATGG"]


def test_the_interpreter_assesses_one_set(tmp_path):
    write_results(tmp_path, ROWS)

    from neoswga.core.results_interpreter import ResultsInterpreter

    loaded = ResultsInterpreter(str(tmp_path))._load_step4_results()

    assert [row["primer"] for row in loaded] == ["ACCGCATC", "TCAGCGA"]


def test_the_fasta_holds_only_the_delivered_set(tmp_path):
    """End to end, because the header numbering is what misleads a reader.

    Records are written SWGA_001 upward with nothing marking a set boundary, so
    a pooled file is not merely wrong, it is unreadable as anything else.
    """
    write_results(tmp_path, ROWS)

    from neoswga.core.export import PrimerExporter

    out = tmp_path / "panel.fasta"
    PrimerExporter.from_results_dir(str(tmp_path)).export_fasta(str(out))

    sequences = [line.strip() for line in out.read_text().splitlines() if not line.startswith(">")]
    assert sequences == ["ACCGCATC", "TCAGCGA"]
