"""Scanning a whole table streams from the counter, not from a text file.

Some consumers genuinely need every k-mer: the candidate loader in
`kmer_counter` and the two Bloom builders. They stream line by line already,
so the change is where the lines come from, not how they are consumed.

Streaming from the database avoids materialising the text table at all. On
Drosophila at k=18 that table is 2.3 GB and takes 9.89 s to write; dumping the
same data to a pipe takes 7.61 s and leaves nothing behind.

The text path is kept and is what every existing data directory uses.
"""

import pytest

from neoswga.core import kmer_tables
from neoswga.core.kmer_backend import select_backend


@pytest.fixture
def counted(tmp_path):
    backend = select_backend()
    if not backend.available():
        pytest.skip(f"the {backend.name} counter is not installed")
    fasta = tmp_path / "g.fna"
    fasta.write_text(">g\n" + "ACGTTGCAAGGCTTAC" * 60 + "\n")
    prefix = str(tmp_path / "g")
    backend.count(str(fasta), 8, prefix)
    return prefix


def test_a_text_table_streams_without_any_counter(tmp_path):
    """The path every existing directory takes. Needs no KMC."""
    (tmp_path / "old_8mer_all.txt").write_text("ACGTTGCA 7\nTGCAACGT 3\n")
    assert dict(kmer_tables.iter_table(str(tmp_path / "old"), 8)) == {
        "ACGTTGCA": 7,
        "TGCAACGT": 3,
    }


def test_an_uncounted_prefix_raises_rather_than_yielding_nothing(tmp_path):
    """An empty iteration reads as "this reference has no k-mers", which is
    the most permissive answer available."""
    with pytest.raises(FileNotFoundError) as excinfo:
        list(kmer_tables.iter_table(str(tmp_path / "never"), 8))
    assert "count-kmers" in str(excinfo.value)


def test_streaming_the_database_yields_what_the_counter_holds(counted):
    """Compared against an independent dump of the same database.

    Not against a text table: the whole point is that counting no longer
    writes one, so a test expecting it beside the database was asserting the
    behaviour this change removes.
    """
    import subprocess

    from neoswga.core.kmer_backend import KmcBackend

    backend = KmcBackend()
    database = kmer_tables._kmc_database(counted, 8)
    if not database:
        pytest.skip("the configured backend is not KMC")

    dumped = subprocess.run(
        [backend.binary("kmc_dump"), database, "/dev/stdout"],
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    expected = {
        line.split()[0]: int(line.split()[1]) for line in dumped.splitlines() if line.strip()
    }
    assert dict(kmer_tables.iter_table(counted, 8)) == expected
    assert expected, "the fixture produced an empty database"


def test_stopping_early_does_not_raise(counted):
    """Consumers stop once they have enough. A closed pipe must not surface
    as BrokenPipeError from the writer."""
    seen = 0
    for _ in kmer_tables.iter_table(counted, 8):
        seen += 1
        if seen > 2:
            break
    assert seen == 3


def test_the_candidate_loader_reads_through_the_table_layer():
    """Pins the PATH. A streaming layer nothing calls is Known Issue 8."""
    import inspect

    from neoswga.core import kmer_counter

    source = inspect.getsource(kmer_counter.get_primer_list_from_kmers)
    assert "iter_table" in source
    assert "mer_all.txt" not in source, "still opening the text table by name"
