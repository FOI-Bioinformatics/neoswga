"""A reference counted into a binary database HAS been counted.

Nine places ask "has this prefix been counted at this k" by testing for
`{prefix}_{k}mer_all.txt`. Once KMC writes a database and no text dump, that
question gets the wrong answer, and it is wrong in the direction that matters:
a counted reference reads as uncounted, so a run either recounts it or skips
the reference entirely.

The predicate has to be right in both directions, and both are costly:

  too strict  every existing data directory stops working, because they all
              hold text tables and no database
  too loose   a prefix that was never counted reads as counted, and whatever
              reads it next gets an empty answer -- which for a background is
              the silent-zero shape of Known Issues 5, 6, 13 and 15

So the three cases are pinned here, and the call sites are pinned separately:
a predicate nothing asks is Known Issue 8.
"""

import inspect

import pytest

from neoswga.core import kmer_tables


def test_a_text_table_counts_as_a_table(tmp_path):
    """The path every existing data directory takes."""
    (tmp_path / "old_8mer_all.txt").write_text("ACGTTGCA 7\n")
    assert kmer_tables.table_exists(str(tmp_path / "old"), 8)


def test_a_binary_database_counts_as_a_table(tmp_path):
    """Both KMC files must be present; one alone is a half-written count."""
    stem = tmp_path / "g_8mer"
    stem.with_suffix(".kmc_pre").write_bytes(b"x")
    stem.with_suffix(".kmc_suf").write_bytes(b"x")
    assert kmer_tables.table_exists(str(tmp_path / "g"), 8)


def test_half_a_database_does_not_count(tmp_path):
    """KMC writes .kmc_pre and .kmc_suf. One without the other is an
    interrupted run, not a table."""
    (tmp_path / "g_8mer.kmc_pre").write_bytes(b"x")
    assert not kmer_tables.table_exists(str(tmp_path / "g"), 8)


def test_a_jellyfish_database_counts_as_a_table(tmp_path):
    (tmp_path / "g_8mer.jf").write_bytes(b"x")
    assert kmer_tables.table_exists(str(tmp_path / "g"), 8)


def test_nothing_present_is_false(tmp_path):
    assert not kmer_tables.table_exists(str(tmp_path / "never"), 8)


def test_a_different_k_is_not_this_k(tmp_path):
    """A prefix counted at k=8 has not been counted at k=12."""
    (tmp_path / "g_8mer_all.txt").write_text("ACGTTGCA 7\n")
    assert kmer_tables.table_exists(str(tmp_path / "g"), 8)
    assert not kmer_tables.table_exists(str(tmp_path / "g"), 12)


@pytest.mark.parametrize(
    "module_name",
    [
        "neoswga.core.pipeline",
        "neoswga.core.param_validator",
        "neoswga.core.background_registry",
        "neoswga.core.genome_library",
    ],
)
def test_no_module_still_tests_for_the_text_file_by_name(module_name):
    """Pins the PATH. A predicate nothing asks is Known Issue 8's class, and
    these are the modules whose checks decide whether a reference is counted.
    """
    import importlib

    source = inspect.getsource(importlib.import_module(module_name))
    offenders = [
        line.strip()
        for line in source.splitlines()
        if "mer_all.txt" in line and "table_exists" not in line and not line.strip().startswith("#")
    ]
    assert not offenders, (
        f"{module_name} still decides whether a reference is counted by "
        f"looking for a text file: {offenders}"
    )
