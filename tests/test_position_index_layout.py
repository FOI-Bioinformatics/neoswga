"""The sorted-blocks position index says exactly what the per-dataset one said.

Every question a reader asks -- which positions, is there an entry at all,
where do the records start, which reference was this built from -- must get
the same answer from either layout, because the layout is an encoding and not
a measurement. The distinction most likely to be lost in re-encoding is the
one between an EMPTY entry (scanned, occurs nowhere) and an ABSENT one (never
scanned), so it is checked in every direction.
"""

import subprocess
import sys
import textwrap

import h5py
import numpy as np
import pytest

from neoswga.core import position_index as pi


def _entries():
    return {
        "ACGTACGTACGT": [5, 900, 12],
        "TTTTTTTTTTTT": [],  # scanned, absent from the genome
        "ACGTACGTACGA": [3_000_000_000],  # past the int32 ceiling (Known Issue 7)
        "AAAAAAAAAAAC": [7],
        "GTTTTTTTTTTT": [8, 9],  # the reverse complement of the previous key
    }


def _write_per_dataset(path, entries, record_starts=None, attrs=None):
    """An index as every writer before 2026-09-25 produced it."""
    with h5py.File(path, "w") as f:
        for name, value in (attrs or {}).items():
            f.attrs[name] = value
        if record_starts is not None:
            f.create_dataset(pi.RECORD_STARTS_KEY, data=np.asarray(record_starts))
        for key, positions in entries.items():
            # An empty list becomes float64, which is what the old writer stored.
            f.create_dataset(key, data=positions)


@pytest.fixture(params=["sorted_blocks", "per_dataset"])
def index(request, tmp_path):
    path = str(tmp_path / "g_12mer_positions.h5")
    if request.param == "sorted_blocks":
        pi.write_entries(path, _entries(), record_starts=[0, 100], attrs={"reference_digest": "d"})
    else:
        _write_per_dataset(path, _entries(), [0, 100], {"reference_digest": "d"})
    with pi.open_index(path) as opened:
        assert opened.layout == request.param
        yield opened


def test_positions_are_returned_as_written(index):
    for key, positions in _entries().items():
        got = index.get(key)
        assert got is not None
        assert got.dtype == np.int64
        assert got.tolist() == positions


def test_an_empty_entry_is_not_an_absent_one(index):
    assert "TTTTTTTTTTTT" in index
    assert index.get("TTTTTTTTTTTT").tolist() == []
    assert "CCCCCCCCCCCC" not in index
    assert index.get("CCCCCCCCCCCC") is None


def test_get_many_leaves_out_what_was_never_scanned(index):
    found = index.get_many(["ACGTACGTACGT", "CCCCCCCCCCCC", "TTTTTTTTTTTT"])
    assert set(found) == {"ACGTACGTACGT", "TTTTTTTTTTTT"}
    assert found["TTTTTTTTTTTT"].tolist() == []


def test_metadata_survives(index):
    assert index.record_starts() == [0, 100]
    assert index.attrs["reference_digest"] == "d"
    assert sorted(index.keys()) == sorted(_entries())
    assert len(index) == len(_entries())


def test_keys_are_not_canonicalised(index):
    # A primer and its reverse complement are separate entries: the reverse
    # strand of one IS the forward entry of the other.
    assert index.get("AAAAAAAAAAAC").tolist() == [7]
    assert index.get("GTTTTTTTTTTT").tolist() == [8, 9]


def test_bookkeeping_names_are_never_entries(index):
    assert pi.RECORD_STARTS_KEY not in index
    assert pi.KEYS_KEY not in index
    assert index.get(pi.RECORD_STARTS_KEY) is None


def test_the_old_layout_is_converted_on_its_first_write(tmp_path):
    path = str(tmp_path / "g_12mer_positions.h5")
    _write_per_dataset(path, _entries(), [0, 50], {"provenance_k": 12})

    pi.write_entries(path, {"CCCCCCCCCCCC": [1]})

    with pi.open_index(path) as index:
        assert index.layout == pi.SORTED_BLOCKS
        assert index.attrs[pi.LAYOUT_ATTR] == pi.SORTED_BLOCKS
        assert index.attrs["provenance_k"] == 12
        assert index.record_starts() == [0, 50]
        expected = {**_entries(), "CCCCCCCCCCCC": [1]}
        assert {k: v.tolist() for k, v in index.items()} == expected


def test_a_merge_overwrites_named_entries_and_keeps_the_rest(tmp_path):
    path = str(tmp_path / "g_12mer_positions.h5")
    pi.write_entries(path, _entries(), record_starts=[0])
    pi.write_entries(
        path,
        {"ACGTACGTACGT": [1, 2, 3, 4, 5], "TTTTTTTTTTTT": [42], "AAAAAAAAAAAA": []},
    )
    with pi.open_index(path) as index:
        assert index.get("ACGTACGTACGT").tolist() == [1, 2, 3, 4, 5]
        assert index.get("TTTTTTTTTTTT").tolist() == [42]
        assert index.get("AAAAAAAAAAAA").tolist() == []
        assert index.get("GTTTTTTTTTTT").tolist() == [8, 9]
        assert index.get("ACGTACGTACGA").tolist() == [3_000_000_000]
        assert index.record_starts() == [0]
        assert index.keys() == sorted(index.keys())


def test_replace_discards_entries_attributes_and_geometry(tmp_path):
    path = str(tmp_path / "g_12mer_positions.h5")
    pi.write_entries(path, _entries(), record_starts=[0, 7], attrs={"provenance_k": 12})
    pi.write_entries(path, {"CCCCCCCCCCCC": [1]}, replace=True)
    with pi.open_index(path) as index:
        assert index.keys() == ["CCCCCCCCCCCC"]
        assert "provenance_k" not in index.attrs
        assert index.record_starts() is None


def test_an_empty_file_reads_as_an_index_with_no_entries(tmp_path):
    # `string_search` creates the file with h5py "a" before its first write.
    path = str(tmp_path / "g_12mer_positions.h5")
    with h5py.File(path, "a"):
        pass
    with pi.open_index(path) as index:
        assert len(index) == 0
        assert index.get("ACGTACGTACGT") is None
    pi.write_entries(path, {"ACGTACGTACGT": [1]})
    with pi.open_index(path) as index:
        assert index.get("ACGTACGTACGT").tolist() == [1]


def test_a_failed_write_leaves_the_previous_index_intact(tmp_path, monkeypatch):
    path = str(tmp_path / "g_12mer_positions.h5")
    pi.write_entries(path, _entries())

    def refuse(*_args, **_kwargs):
        raise OSError("disk full")

    monkeypatch.setattr(pi.os, "replace", refuse)
    with pytest.raises(OSError, match="disk full"):
        pi.write_entries(path, {"ACGTACGTACGT": [0]})

    with pi.open_index(path) as index:
        assert index.get("ACGTACGTACGT").tolist() == [5, 900, 12]
    assert list(tmp_path.iterdir()) == [tmp_path / "g_12mer_positions.h5"]


def test_coalesced_reads_agree_with_single_reads(tmp_path, monkeypatch):
    rng = np.random.default_rng(3)
    alphabet = np.array(list("ACGT"))
    entries = {}
    for _ in range(3000):
        key = "".join(rng.choice(alphabet, 10))
        entries[key] = rng.integers(0, 10**9, rng.integers(0, 6)).tolist()
    path = str(tmp_path / "g_10mer_positions.h5")
    pi.write_entries(path, entries)

    # Small limits, so gaps, run breaks and the run cap are all exercised.
    monkeypatch.setattr(pi, "_COALESCE_GAP", 3)
    monkeypatch.setattr(pi, "_MAX_RUN", 20)
    wanted = list(rng.choice(sorted(entries), 900, replace=False)) + ["NNNNNNNNNN"]
    with pi.open_index(path) as index:
        batch = index.get_many(wanted)
        assert "NNNNNNNNNN" not in batch
        for key in wanted[:-1]:
            assert batch[key].tolist() == index.get(key).tolist() == entries[key]


def test_a_second_process_reading_the_index_still_blocks_a_writer(tmp_path):
    """Known Issue 21: two runs in one directory must collide loudly.

    Writing to a temporary file and renaming it would otherwise let both
    writers succeed, and the later rename would drop the earlier run's entries
    without a word.
    """
    path = str(tmp_path / "g_12mer_positions.h5")
    pi.write_entries(path, _entries())
    holder = subprocess.Popen(
        [
            sys.executable,
            "-c",
            textwrap.dedent(f"""
                import sys, h5py
                f = h5py.File({path!r}, "r")
                print("open", flush=True)
                sys.stdin.read()
                """),
        ],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        text=True,
    )
    try:
        assert holder.stdout.readline().strip() == "open"
        with pytest.raises(OSError, match="lock"):
            pi.write_entries(path, {"CCCCCCCCCCCC": [1]})
    finally:
        holder.communicate("")
    with pi.open_index(path) as index:
        assert "CCCCCCCCCCCC" not in index


def test_a_file_an_older_release_wrote_into_is_refused(tmp_path):
    """An older release reads no entry from sorted blocks, rescans everything,
    and writes per-k-mer datasets beside them. Reading only the blocks would
    ignore what it wrote, so the mixed file is refused as unreadable, which
    is what makes the scan path rebuild it."""
    path = str(tmp_path / "g_12mer_positions.h5")
    pi.write_entries(path, _entries())
    with h5py.File(path, "r+") as f:
        f.create_dataset("CCCCCCCCCCCC", data=[1, 2])
    with pytest.raises(pi.MixedLayoutError, match="older NeoSWGA release"):
        pi.open_index(path)
    with pytest.raises(OSError):
        pi.write_entries(path, {"ACGTACGTACGT": [1]})
    pi.write_entries(path, {"ACGTACGTACGT": [1]}, replace=True)
    with pi.open_index(path) as index:
        assert index.keys() == ["ACGTACGTACGT"]
