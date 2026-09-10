"""The foreground GC fraction was recomputed from the FASTA once per step.

`_apply_gc_adaptive_defaults` needs it, so `parameter.get_params` loaded every
foreground genome to count G and C: 0.30 s on a 4.6 Mb target, four times per
pipeline, scaling with target size. The value cannot change between steps of one
run, so it is cached, keyed on each file's path, size and mtime.

It is cached rather than skipped deliberately. The derivation is what fills
`effective_conditions`, and the optimize step now warns when its conditions
differ from the filter step's, so a `score` step that skipped the derivation
would make the default pipeline warn about itself.
"""

import json
import time

import pytest

from neoswga.core import genome_gc_cache as gc_cache


@pytest.fixture
def genome(tmp_path):
    path = tmp_path / "target.fasta"
    # 40 bases, 20 of them G or C: GC = 0.5 exactly.
    path.write_text(">t\n" + "ACGT" * 10 + "\n")
    return path


def test_a_written_value_reads_back(tmp_path, genome):
    gc_cache.write_cached_gc(str(tmp_path), [str(genome)], 0.5077)
    assert gc_cache.read_cached_gc(str(tmp_path), [str(genome)]) == 0.5077


def test_an_empty_cache_reads_none(tmp_path, genome):
    assert gc_cache.read_cached_gc(str(tmp_path), [str(genome)]) is None


def test_editing_the_genome_invalidates_the_cache(tmp_path, genome):
    gc_cache.write_cached_gc(str(tmp_path), [str(genome)], 0.5077)
    time.sleep(0.01)
    genome.write_text(">t\n" + "AAAA" * 10 + "\n")

    assert (
        gc_cache.read_cached_gc(str(tmp_path), [str(genome)]) is None
    ), "a changed foreground genome must not be scored under the old GC value"


def test_a_same_size_edit_invalidates_the_cache(tmp_path, genome):
    """mtime alone is the load-bearing half of the key here.

    A same-length edit -- swapping bases, which is what a corrected assembly
    often is -- leaves st_size identical. A key built from path and size only
    would hand the new genome the old GC fraction.
    """
    gc_cache.write_cached_gc(str(tmp_path), [str(genome)], 0.5)
    time.sleep(0.01)
    genome.write_text(">t\n" + "GGCC" * 10 + "\n")

    assert gc_cache.read_cached_gc(str(tmp_path), [str(genome)]) is None


def test_a_different_genome_list_invalidates_the_cache(tmp_path, genome):
    other = tmp_path / "other.fasta"
    other.write_text(">o\n" + "GGCC" * 10 + "\n")

    gc_cache.write_cached_gc(str(tmp_path), [str(genome)], 0.5077)
    assert gc_cache.read_cached_gc(str(tmp_path), [str(genome), str(other)]) is None


def test_a_missing_genome_yields_no_key(tmp_path):
    assert gc_cache.cache_key([str(tmp_path / "gone.fasta")]) is None


def test_a_corrupt_cache_file_is_ignored(tmp_path, genome):
    (tmp_path / gc_cache.CACHE_FILENAME).write_text("{not json")
    assert gc_cache.read_cached_gc(str(tmp_path), [str(genome)]) is None


def test_get_params_populates_and_then_reuses_the_cache(tmp_path):
    """The value must be identical whether it was computed or read back."""
    from neoswga.core import parameter

    fasta = tmp_path / "target.fasta"
    fasta.write_text(">t\n" + "ACGTGGCC" * 250 + "\n")
    params_file = tmp_path / "params.json"
    params_file.write_text(
        json.dumps(
            {
                "data_dir": str(tmp_path),
                "fg_genomes": [str(fasta)],
                "fg_prefixes": [str(tmp_path / "target")],
                "bg_genomes": [],
                "bg_prefixes": [],
                "fg_seq_lengths": [2000],
                "bg_seq_lengths": [],
                "cpus": 1,
            }
        )
    )

    class _Args:
        json_file = str(params_file)

        def __getattr__(self, name):
            return None

    first = parameter.get_params(_Args())["genome_gc"]
    assert (tmp_path / gc_cache.CACHE_FILENAME).is_file()

    cached = gc_cache.read_cached_gc(str(tmp_path), [str(fasta)])
    assert cached == first

    second = parameter.get_params(_Args())["genome_gc"]
    assert second == first, (
        "the cached value must match the computed one exactly; the recorded "
        "effective_conditions have to agree across steps"
    )


def test_the_cached_value_matches_a_direct_count(tmp_path):
    """The cache must return the same number the FASTA read produced, not a
    rounded or otherwise mangled one."""
    from neoswga.core import parameter

    fasta = tmp_path / "target.fasta"
    # 3 of every 8 bases are G or C: 0.375, which no rounding preserves by luck.
    fasta.write_text(">t\n" + "ACGTAGTA" * 250 + "\n")
    params_file = tmp_path / "params.json"
    params_file.write_text(
        json.dumps(
            {
                "data_dir": str(tmp_path),
                "fg_genomes": [str(fasta)],
                "fg_prefixes": [str(tmp_path / "target")],
                "bg_genomes": [],
                "bg_prefixes": [],
                "fg_seq_lengths": [2000],
                "bg_seq_lengths": [],
                "cpus": 1,
            }
        )
    )

    class _Args:
        json_file = str(params_file)

        def __getattr__(self, name):
            return None

    assert parameter.get_params(_Args())["genome_gc"] == 0.375
    assert gc_cache.read_cached_gc(str(tmp_path), [str(fasta)]) == 0.375


def test_an_explicit_genome_gc_still_wins(tmp_path):
    """A user who pins genome_gc must not be overridden by the cache."""
    from neoswga.core import parameter

    fasta = tmp_path / "target.fasta"
    fasta.write_text(">t\n" + "ACGT" * 500 + "\n")
    gc_cache.write_cached_gc(str(tmp_path), [str(fasta)], 0.9)

    params_file = tmp_path / "params.json"
    params_file.write_text(
        json.dumps(
            {
                "data_dir": str(tmp_path),
                "fg_genomes": [str(fasta)],
                "fg_prefixes": [str(tmp_path / "target")],
                "bg_genomes": [],
                "bg_prefixes": [],
                "fg_seq_lengths": [2000],
                "bg_seq_lengths": [],
                "genome_gc": 0.42,
                "cpus": 1,
            }
        )
    )

    class _Args:
        json_file = str(params_file)

        def __getattr__(self, name):
            return None

    assert parameter.get_params(_Args())["genome_gc"] == 0.42
