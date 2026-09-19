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

from neoswga.core import kmer_counter, pipeline


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
        json.dump(
            {
                "genome": os.path.abspath(fake_genome),
                "fingerprint": kmer_counter.genome_fingerprint(fake_genome),
                "digest_algorithm": kmer_counter.DIGEST_ALGORITHM,
            },
            fh,
        )

    called = []
    monkeypatch.setattr(kmer_counter.subprocess, "run", lambda *a, **k: called.append(a) or None)

    kmer_counter._run_jellyfish_for_k(prefix, fake_genome, 12, cpus=1, hash_size=1000)
    assert called == [], "an up-to-date table must not be recounted"


def test_stale_provenance_is_detected(tmp_path, fake_genome, other_genome, caplog):
    prefix = str(tmp_path / "fg")
    open(f"{prefix}_12mer_all.txt", "w").write("TTTTGGGGCCTT 4\n")
    with open(kmer_counter.table_provenance_path(prefix, 12), "w") as fh:
        json.dump(
            {
                "genome": os.path.abspath(other_genome),
                "fingerprint": kmer_counter.genome_fingerprint(other_genome),
                # A CURRENT record naming a different genome, which is the case
                # this test is about. Without the algorithm it would be
                # incomparable instead, and the message would rightly say so.
                "digest_algorithm": kmer_counter.DIGEST_ALGORITHM,
            },
            fh,
        )

    with caplog.at_level(logging.WARNING):
        current = kmer_counter._table_is_current(prefix, fake_genome, 12)

    assert current is False
    assert any("counted from a different genome" in r.getMessage() for r in caplog.records)


def test_missing_provenance_is_detected(tmp_path, fake_genome, caplog):
    """A table from before this change has no provenance and cannot be trusted."""
    prefix = str(tmp_path / "fg")
    open(f"{prefix}_12mer_all.txt", "w").write("ACGTACGTACGT 4\n")

    with caplog.at_level(logging.WARNING):
        current = kmer_counter._table_is_current(prefix, fake_genome, 12)

    assert current is False
    assert any("no provenance record" in r.getMessage() for r in caplog.records)


class TestStep2ChecksTheSameRecord:
    """Step 1 was the only step that read the provenance record.

    A user who repoints fg_genomes at a new assembly and runs `filter` without
    re-running `count-kmers` reaches step 2 with tables from the previous
    organism. validate_step2_prerequisites checked only that the text file
    existed, so the whole design was built from those counts.
    """

    def _table(self, tmp_path, genome, k=12):
        prefix = str(tmp_path / "fg")
        open(f"{prefix}_{k}mer_all.txt", "w").write("ACGTACGTACGT 4\n")
        with open(kmer_counter.table_provenance_path(prefix, k), "w") as fh:
            json.dump(
                {
                    "genome": os.path.abspath(genome),
                    "fingerprint": kmer_counter.genome_fingerprint(genome),
                    # Written as today's code writes it. Without the algorithm
                    # the record is UNKNOWN rather than comparable, and step 2
                    # skips it -- which is right for a record from the partial
                    # hash and wrong for one this helper means to be current.
                    "digest_algorithm": kmer_counter.DIGEST_ALGORITHM,
                    "k": k,
                },
                fh,
            )
        return prefix

    def test_a_table_from_another_genome_is_refused(self, tmp_path, fake_genome, other_genome):
        prefix = self._table(tmp_path, other_genome)
        result = pipeline.validate_step2_prerequisites(
            str(tmp_path), [prefix], [], 12, 12, fg_genomes=[fake_genome]
        )
        assert result.valid is False
        assert "different genome" in result.error_message
        assert "count-kmers" in result.remediation

    def test_a_table_from_this_genome_is_accepted(self, tmp_path, fake_genome):
        prefix = self._table(tmp_path, fake_genome)
        result = pipeline.validate_step2_prerequisites(
            str(tmp_path), [prefix], [], 12, 12, fg_genomes=[fake_genome]
        )
        assert result.valid is True

    def test_a_background_table_from_another_genome_is_refused(
        self, tmp_path, fake_genome, other_genome
    ):
        prefix = self._table(tmp_path, other_genome)
        result = pipeline.validate_step2_prerequisites(
            str(tmp_path), [], [prefix], 12, 12, bg_genomes=[fake_genome]
        )
        assert result.valid is False

    def test_a_table_with_no_provenance_record_is_not_refused(self, tmp_path, fake_genome):
        """Every table written before the record existed lacks one. Absence is
        unknown, not wrong; step 1 recounts and writes the record."""
        prefix = str(tmp_path / "fg")
        open(f"{prefix}_12mer_all.txt", "w").write("ACGTACGTACGT 4\n")
        result = pipeline.validate_step2_prerequisites(
            str(tmp_path), [prefix], [], 12, 12, fg_genomes=[fake_genome]
        )
        assert result.valid is True

    def test_omitting_the_genome_paths_keeps_the_old_behaviour(
        self, tmp_path, fake_genome, other_genome
    ):
        """Callers that pass prefixes only must still validate."""
        prefix = self._table(tmp_path, other_genome)
        result = pipeline.validate_step2_prerequisites(str(tmp_path), [prefix], [], 12, 12)
        assert result.valid is True
