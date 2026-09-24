"""`calibrate-reach`'s advice must not send a multi-record user into a wall.

When no prefix matched a BAM contig the command said only "Map one explicitly
with --contig-alias FG=BAMCONTIG". For a MULTI-RECORD reference, following
that advice was worse than the state it was offered from.

A prefix is a FASTA FILE and `fg_seq_lengths` is the concatenated total across
every record in it -- `utility.get_seq_length` is `sum(1 for _ in
read_fasta_file(genome))` over CHARACTERS. So an alias binding the file to one
contig then asked for depth over the whole concatenation, and the tail past
that contig was padded with zeros.

Measured on the shipped `tests/validation/genomes/params.json` (Prevotella,
two chromosomes, 3,168,282 bp total, first chromosome 1,796,408):

    3,168,282 requested - 1,796,408 real = 1,371,874 fabricated, 43.3%

`_require_contig_covers` now refuses that, so the advice leads to a refusal
rather than a wrong answer. These pin the better outcome: it leads to neither,
because the message says up front that this command cannot describe a
multi-record reference.

Note this could not have been caught by the alias fix in
`test_one_contig_alias_means_one_thing.py`. `calibrate-reach` is the one
command where the prefix-keyed alias works exactly as documented.
"""

import pytest

from neoswga.cli.commands import _no_contig_matched_advice
from neoswga.core import parameter


@pytest.fixture
def fasta(tmp_path):
    def build(name, records):
        path = tmp_path / name
        path.write_text("".join(f">{r}\n" + "ACGT" * 100 + "\n" for r in records))
        return str(path)

    return build


def test_a_multi_record_reference_is_told_the_flag_cannot_help(fasta, monkeypatch):
    genome = fasta("prevotella.fna", ["NC_014370.1", "NC_014371.1"])
    monkeypatch.setattr(parameter, "fg_genomes", [genome], raising=False)

    message = _no_contig_matched_advice(["prevotella"])

    assert "2 records" in message
    assert "prevotella.fna" in message
    assert "will not help" in message, "it must say the flag cannot help, not offer it"
    assert "concatenation" in message, "it must say why"


def test_a_single_record_reference_still_gets_the_original_advice(fasta, monkeypatch):
    """Where the alias DOES work, and where it is the right thing to suggest.
    Narrowing the message to the case that needs it must not remove it from
    the case it was written for."""
    genome = fasta("one.fna", ["chr1"])
    monkeypatch.setattr(parameter, "fg_genomes", [genome], raising=False)

    message = _no_contig_matched_advice(["one"])

    assert "--contig-alias" in message
    assert "will not help" not in message


def test_an_unreadable_reference_falls_back_rather_than_raising(monkeypatch):
    """This runs while REPORTING a failure. Raising from inside the advice
    would replace a useful message with a traceback about the advice."""
    monkeypatch.setattr(parameter, "fg_genomes", ["/nonexistent/x.fna"], raising=False)

    assert "--contig-alias" in _no_contig_matched_advice(["x"])
