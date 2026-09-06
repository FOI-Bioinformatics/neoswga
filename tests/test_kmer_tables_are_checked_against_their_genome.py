"""A k-mer table is reused only when it was counted from the same genome.

Finding C1. The existence check keys on the output prefix, which comes from
params.json, not on the genome. Repointing fg_genomes at a new assembly while
leaving fg_prefixes alone silently built the whole design from the previous
organism's counts, with no warning at any step.
"""

import json
import logging
import os

import pytest

from neoswga.core import kmer_counter


@pytest.fixture
def fake_genome(tmp_path):
    path = tmp_path / "genome.fna"
    path.write_text(">chr1\n" + "ACGTACGTAC" * 20 + "\n")
    return str(path)


@pytest.fixture
def other_genome(tmp_path):
    path = tmp_path / "other.fna"
    path.write_text(">chr1\n" + "TTTTGGGGCC" * 20 + "\n")
    return str(path)


def test_fingerprint_differs_between_genomes(fake_genome, other_genome):
    assert kmer_counter.genome_fingerprint(fake_genome) != kmer_counter.genome_fingerprint(
        other_genome
    )


def test_fingerprint_is_stable_for_one_genome(fake_genome):
    assert kmer_counter.genome_fingerprint(fake_genome) == kmer_counter.genome_fingerprint(
        fake_genome
    )


def test_matching_provenance_skips_the_count(tmp_path, fake_genome, monkeypatch):
    prefix = str(tmp_path / "fg")
    open(f"{prefix}_12mer_all.txt", "w").write("ACGTACGTACGT 4\n")
    with open(kmer_counter.table_provenance_path(prefix, 12), "w") as fh:
        json.dump({"genome": os.path.abspath(fake_genome),
                   "fingerprint": kmer_counter.genome_fingerprint(fake_genome)}, fh)

    called = []
    monkeypatch.setattr(kmer_counter.subprocess, "run",
                        lambda *a, **k: called.append(a) or None)

    kmer_counter._run_jellyfish_for_k(prefix, fake_genome, 12, cpus=1, hash_size=1000)
    assert called == [], "an up-to-date table must not be recounted"


def test_stale_provenance_is_detected(tmp_path, fake_genome, other_genome, caplog):
    prefix = str(tmp_path / "fg")
    open(f"{prefix}_12mer_all.txt", "w").write("TTTTGGGGCCTT 4\n")
    with open(kmer_counter.table_provenance_path(prefix, 12), "w") as fh:
        json.dump({"genome": os.path.abspath(other_genome),
                   "fingerprint": kmer_counter.genome_fingerprint(other_genome)}, fh)

    with caplog.at_level(logging.WARNING):
        current = kmer_counter._table_is_current(prefix, fake_genome, 12)

    assert current is False
    assert any(
        "counted from a different genome" in r.getMessage() for r in caplog.records
    )


def test_missing_provenance_is_detected(tmp_path, fake_genome, caplog):
    """A table from before this change has no provenance and cannot be trusted."""
    prefix = str(tmp_path / "fg")
    open(f"{prefix}_12mer_all.txt", "w").write("ACGTACGTACGT 4\n")

    with caplog.at_level(logging.WARNING):
        current = kmer_counter._table_is_current(prefix, fake_genome, 12)

    assert current is False
    assert any("no provenance record" in r.getMessage() for r in caplog.records)
