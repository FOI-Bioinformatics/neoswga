"""A fingerprint that misses the middle of a file is not a fingerprint.

Finding F6 of the 2026-09-16 pipeline audit. `genome_fingerprint` hashed the
file size plus the first and last 1 MB, so a substitution anywhere in the
middle of a file over 2 MB left it unchanged. A same-length consensus or a
sample-specific assembly therefore reused the previous genome's k-mer counts,
its position index and its inventory, silently and through the whole design.

That is the same class as repointing `fg_genomes` without re-running
`count-kmers`, which the sidecar was introduced to catch. It caught the case
where the file CHANGED SIZE and missed the one where it did not.

The stated reason for the partial hash was cost, and the measurement does not
support it: SHA-256 runs at about 2.5 GB/s here, so hg38 is around a second,
computed once per input per run rather than once per k. The counting step reads
the whole file anyway.
"""

import hashlib
import os

import pytest

from neoswga.core.kmer_counter import DIGEST_ALGORITHM, genome_fingerprint

WINDOW = 1024 * 1024


def _genome(tmp_path, name, middle_byte):
    """A file over 2 MB whose head and tail are identical across calls."""
    path = tmp_path / name
    with open(path, "wb") as handle:
        handle.write(b">chr1\n")
        handle.write(b"A" * WINDOW)
        handle.write(middle_byte * WINDOW)
        handle.write(b"C" * WINDOW)
    return str(path)


class TestTheMiddleIsCovered:
    def test_two_genomes_differing_only_in_the_middle_differ(self, tmp_path):
        """The defect, pinned. These were identical before: same size, same
        first megabyte, same last megabyte."""
        one = _genome(tmp_path, "a.fna", b"G")
        two = _genome(tmp_path, "b.fna", b"T")

        assert genome_fingerprint(one) != genome_fingerprint(two)

    def test_a_single_byte_in_the_middle_changes_it(self, tmp_path):
        """A consensus differing from its reference at one position is exactly
        the case a sample-specific assembly presents."""
        path = _genome(tmp_path, "c.fna", b"G")
        before = genome_fingerprint(path)

        with open(path, "r+b") as handle:
            handle.seek(WINDOW + 100)
            handle.write(b"T")

        assert genome_fingerprint(path) != before

    def test_identical_files_agree(self, tmp_path):
        one = _genome(tmp_path, "a.fna", b"G")
        two = _genome(tmp_path, "b.fna", b"G")

        assert genome_fingerprint(one) == genome_fingerprint(two)

    def test_it_is_the_digest_of_the_whole_file(self, tmp_path):
        """Not a scheme of its own, so a reader can verify it with `shasum`."""
        path = _genome(tmp_path, "a.fna", b"G")

        with open(path, "rb") as handle:
            expected = hashlib.sha256(handle.read()).hexdigest()

        assert genome_fingerprint(path) == expected


class TestItIsComputedOncePerInput:
    def test_a_second_call_does_not_re_read_the_file(self, tmp_path, monkeypatch):
        """`run_jellyfish` fingerprints once per k. Seven k values meant seven
        passes over the genome, which is where the cost objection came from."""
        path = _genome(tmp_path, "a.fna", b"G")
        genome_fingerprint(path)

        reads = []
        real_open = open

        def counting_open(*args, **kwargs):
            if args and str(args[0]) == path:
                reads.append(args[0])
            return real_open(*args, **kwargs)

        monkeypatch.setattr("builtins.open", counting_open)
        genome_fingerprint(path)

        assert reads == [], "the digest was recomputed rather than reused"

    def test_a_changed_file_is_not_served_from_the_cache(self, tmp_path):
        """Caching on the path alone would make the guard useless: the file it
        exists to notice changing is the one it would stop looking at."""
        path = _genome(tmp_path, "a.fna", b"G")
        before = genome_fingerprint(path)

        with open(path, "r+b") as handle:
            handle.seek(WINDOW + 5)
            handle.write(b"T")
        os.utime(path, (0, 0))

        assert genome_fingerprint(path) != before


class TestTheAlgorithmIsRecorded:
    def test_it_is_named(self):
        """A sidecar written under the old partial hash must be recognisable as
        such, so it can be treated as UNKNOWN and recounted once rather than
        reported as a mismatch that alarms without cause."""
        assert DIGEST_ALGORITHM

    def test_the_provenance_record_carries_it(self, tmp_path):
        import json

        from neoswga.core.kmer_counter import _write_table_provenance, table_provenance_path

        genome = _genome(tmp_path, "a.fna", b"G")
        prefix = str(tmp_path / "out")

        _write_table_provenance(prefix, genome, 12)
        record = json.loads(open(table_provenance_path(prefix, 12)).read())

        assert record["digest_algorithm"] == DIGEST_ALGORITHM

    def test_a_record_without_the_algorithm_is_not_current(self, tmp_path):
        """Written under the partial hash. Its fingerprint cannot be compared,
        so the table is recounted once rather than trusted."""
        import json

        from neoswga.core.kmer_counter import _table_is_current, table_provenance_path

        genome = _genome(tmp_path, "a.fna", b"G")
        prefix = str(tmp_path / "out")
        open(f"{prefix}_12mer_all.txt", "w").write("")
        with open(table_provenance_path(prefix, 12), "w") as handle:
            json.dump(
                {"genome": os.path.abspath(genome), "fingerprint": "old-style", "k": 12},
                handle,
            )

        assert _table_is_current(prefix, genome, 12) is False

    def test_a_record_with_the_algorithm_and_a_matching_digest_is_current(self, tmp_path):
        from neoswga.core.kmer_counter import _table_is_current, _write_table_provenance

        genome = _genome(tmp_path, "a.fna", b"G")
        prefix = str(tmp_path / "out")
        open(f"{prefix}_12mer_all.txt", "w").write("")
        _write_table_provenance(prefix, genome, 12)

        assert _table_is_current(prefix, genome, 12) is True


class TestTheIndexCarriesItToo:
    def test_the_position_index_records_the_same_digest(self):
        """The plan asks for it in the index and the inventory provenance, not
        only the count tables: a stale index is as silent as a stale table."""
        import inspect

        from neoswga.core import string_search

        source = inspect.getsource(string_search)

        assert "genome_fingerprint" in source
        assert "reference_digest" in source


class TestUpgradingDoesNotRefuseAWorkingDirectory:
    """The regression this nearly shipped, pinned.

    Making `_table_is_current` false for an old-algorithm record was right for
    step 1, which recounts. Step 2 uses the same predicate to decide whether a
    table was counted from ANOTHER GENOME, and refuses if so. So every existing
    data directory failed step 2 on upgrade with "counted from a different
    genome", which is both alarming and untrue. 28 tests caught it.

    Unknown is not wrong. An absent record already got that treatment; a record
    whose fingerprint cannot be compared gets it too.
    """

    def test_an_old_style_record_does_not_read_as_another_genome(self, tmp_path):
        import json

        from neoswga.core.kmer_counter import table_provenance_path
        from neoswga.core.pipeline import _tables_counted_from_another_genome

        genome = _genome(tmp_path, "a.fna", b"G")
        prefix = str(tmp_path / "out")
        open(f"{prefix}_12mer_all.txt", "w").write("")
        with open(table_provenance_path(prefix, 12), "w") as handle:
            json.dump(
                {"genome": os.path.abspath(genome), "fingerprint": "old-style", "k": 12},
                handle,
            )

        stale = _tables_counted_from_another_genome([prefix], [genome], 12, 12)

        assert stale == [], (
            "an incomparable record must not be reported as a different "
            "genome; step 2 refuses on that and every existing directory "
            "would fail on upgrade"
        )

    def test_a_genuinely_different_genome_is_still_caught(self, tmp_path):
        """The guard must survive the accommodation."""
        from neoswga.core.kmer_counter import _write_table_provenance
        from neoswga.core.pipeline import _tables_counted_from_another_genome

        one = _genome(tmp_path, "a.fna", b"G")
        two = _genome(tmp_path, "b.fna", b"T")
        prefix = str(tmp_path / "out")
        open(f"{prefix}_12mer_all.txt", "w").write("")
        _write_table_provenance(prefix, one, 12)

        stale = _tables_counted_from_another_genome([prefix], [two], 12, 12)

        assert len(stale) == 1 and stale[0].endswith("_12mer_all.txt"), stale

    def test_comparability_is_its_own_question(self):
        from neoswga.core.kmer_counter import table_provenance_is_comparable

        assert table_provenance_is_comparable("/nonexistent/prefix", 12) is False
