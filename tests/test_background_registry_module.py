"""Coverage for neoswga.core.background_registry (CRUD over a JSON registry)."""

import pytest

from neoswga.core.background_registry import (
    BackgroundEntry,
    BackgroundRegistry,
    get_common_background,
)


@pytest.fixture
def registry(tmp_path):
    return BackgroundRegistry(registry_path=str(tmp_path / "reg.json"))


def test_add_get_list(registry):
    entry = registry.add("Human", "Homo sapiens", 3_000_000, bloom_path="/tmp/h.bloom")
    assert entry.name == "Human"
    assert registry.get("Human") is entry
    assert registry.get("missing") is None
    names = [e.name for e in registry.list_all()]
    assert "Human" in names


def test_remove(registry):
    registry.add("E. coli", "Escherichia coli", 4_600_000)
    assert registry.remove("E. coli") is True
    assert registry.get("E. coli") is None
    assert registry.remove("E. coli") is False  # already gone


def test_persistence_across_instances(tmp_path):
    path = str(tmp_path / "reg.json")
    r1 = BackgroundRegistry(registry_path=path)
    r1.add("Mouse", "Mus musculus", 2_700_000, kmer_prefix="/data/mouse")
    # A fresh registry reading the same file must see the entry.
    r2 = BackgroundRegistry(registry_path=path)
    got = r2.get("Mouse")
    assert got is not None and got.kmer_prefix == "/data/mouse"


def test_search(registry):
    registry.add("Human", "Homo sapiens", 3_000_000)
    registry.add("Mouse", "Mus musculus", 2_700_000)
    hits = [e.name for e in registry.search("homo")]
    assert "Human" in hits and "Mouse" not in hits


def test_entry_flags_and_roundtrip(tmp_path):
    # has_bloom / has_kmers check that the files actually exist on disk.
    bloom = tmp_path / "x.bloom"
    bloom.write_bytes(b"x")
    kmer_prefix = str(tmp_path / "x")
    (tmp_path / "x_6mer_all.txt").write_text("AAAAAA 1\n")

    e = BackgroundEntry(
        name="X",
        species="Xx",
        genome_size=1000,
        description="desc",
        bloom_path=str(bloom),
        kmer_prefix=kmer_prefix,
        k_range=(6, 12),
    )
    assert e.has_bloom is True
    assert e.has_kmers is True

    blank = BackgroundEntry(name="Y", species="Yy", genome_size=1000, bloom_path="/no/such.bloom")
    assert blank.has_bloom is False and blank.has_kmers is False

    # to_dict / from_dict roundtrip
    restored = BackgroundEntry.from_dict(e.to_dict())
    assert restored.name == "X" and restored.bloom_path == str(bloom)


def test_get_common_background_known_and_unknown():
    unknown = get_common_background("definitely-not-a-real-organism")
    assert unknown is None
    known = get_common_background("human")
    assert known is not None and known["species"] == "Homo sapiens"
