"""Position files that already exist must not be scanned again.

`pipeline.step2` computed `position_files_exist` and fed it to nothing but a log
line reading "Reusing existing position files (incremental update only)". On the
Aho-Corasick branch the automaton was rebuilt from the whole primer list every
time, so re-running `filter` after a threshold tweak re-paid the entire scan
while announcing that it had not.
`check_which_primers_absent_in_h5py` was reached only from the fallback path.

The reuse has to keep the returned position map complete.
`primer_attributes.get_gini_from_txt_for_one_k` reads it with
`position_cache.get((prefix, primer), [])`, so a primer that is skipped and left
out of the map comes back with no positions, which the Gini gate turns into NaN
and drops. On a second run that would empty the entire candidate pool.

Reuse is only correct under the conditions the file was produced under, so it is
gated on a provenance record holding the genome fingerprint and the scan
parameters. Before that gate the change traded a silent rescan for a silent
wrong answer: the default path had rescanned every pattern on every run and so
was immune to a stale file, and reuse made the file authoritative.
"""

import math
import random

import h5py
import pytest

from neoswga.core import parameter, primer_attributes
from neoswga.core import string_search as ss

A12 = "GCATTACGGTAC"
B12 = "TTGACCATGACG"


@pytest.fixture
def genome(tmp_path):
    """Two 12-mers planted at known offsets in a fixed random genome.

    Seed 11 was checked: each primer occurs exactly once and neither reverse
    complement occurs at all, so the expected position lists are exact.
    """
    rng = random.Random(11)
    seq = "".join(rng.choice("ACGT") for _ in range(6_000))
    seq = seq[:1_000] + A12 + seq[1_000 + len(A12) :]
    seq = seq[:5_000] + B12 + seq[5_000 + len(B12) :]

    fasta = tmp_path / "g.fasta"
    fasta.write_text(">g\n" + seq + "\n")
    ss.clear_genome_cache()
    yield {"fasta": str(fasta), "prefix": str(tmp_path / "g")}
    ss.clear_genome_cache()


@pytest.fixture(autouse=True)
def k_window():
    before = (parameter.min_k, parameter.max_k)
    parameter.min_k, parameter.max_k = 6, 12
    yield
    parameter.min_k, parameter.max_k = before


def _record_scans(monkeypatch):
    """Record how many patterns each automaton build was handed."""
    calls = []
    real = ss.get_all_positions_multi_k

    def spy(primer_lists_by_k, seq_fname, circular, chunk_size=None):
        calls.append(sum(len(v) for v in primer_lists_by_k.values()))
        return real(primer_lists_by_k, seq_fname, circular, chunk_size)

    monkeypatch.setattr(ss, "get_all_positions_multi_k", spy)
    return calls


def test_a_repeated_scan_builds_no_automaton(genome, monkeypatch):
    calls = _record_scans(monkeypatch)

    ss.get_positions([A12, B12], [genome["prefix"]], [genome["fasta"]], circular=False)
    ss.get_positions([A12, B12], [genome["prefix"]], [genome["fasta"]], circular=False)

    assert calls == [4], (
        "the second call rebuilt the automaton from the full primer list; "
        f"pattern counts per scan were {calls} (two primers plus two reverse "
        "complements is one scan of 4)"
    )


def test_the_reused_result_is_identical_to_the_scanned_one(genome):
    """The returned map must stay complete, or the Gini gate empties the pool."""
    first = ss.get_positions([A12, B12], [genome["prefix"]], [genome["fasta"]], circular=False)
    second = ss.get_positions([A12, B12], [genome["prefix"]], [genome["fasta"]], circular=False)

    assert second == first, (
        "the reused map differs from the scanned one; a primer missing from it "
        "reads downstream as a primer that binds nowhere"
    )
    assert first[(genome["prefix"], A12)] == [1000]
    assert first[(genome["prefix"], B12)] == [5000]


def test_only_the_primers_not_already_scanned_are_handed_to_the_automaton(genome, monkeypatch):
    """A pool that grows by one primer costs one primer's worth of scanning."""
    ss.get_positions([A12], [genome["prefix"]], [genome["fasta"]], circular=False)

    calls = _record_scans(monkeypatch)
    result = ss.get_positions([A12, B12], [genome["prefix"]], [genome["fasta"]], circular=False)

    assert calls == [2], (
        "expected one scan of B12 and its reverse complement only; pattern "
        f"counts per scan were {calls}"
    )
    assert result[(genome["prefix"], A12)] == [1000], "the cached primer was dropped"
    assert result[(genome["prefix"], B12)] == [5000]


def test_a_second_run_does_not_empty_the_candidate_pool(genome):
    """The failure this reuse could introduce, measured where it would be felt.

    A primer skipped by the reuse and left out of the returned map reaches
    `get_gini_from_txt_for_one_k` as a primer with no binding sites, which
    scores NaN, and `filter.get_gini` drops every NaN. On a second run every
    foreground primer is a skipped primer, so the whole pool would be dropped
    and `filter` would write an empty result while logging success. This test
    fails if the map comes back short, whatever the reason.
    """
    primers = [A12, B12]
    ss.get_positions(primers, [genome["prefix"]], [genome["fasta"]], circular=False)
    reused = ss.get_positions(primers, [genome["prefix"]], [genome["fasta"]], circular=False)

    ginis = primer_attributes.get_gini_from_txt_for_one_k(
        primers,
        genome["prefix"],
        genome["fasta"],
        seq_length=6_000,
        circular=False,
        position_cache=reused,
    )

    survivors = [
        primer
        for primer, (forward, reverse) in ginis.items()
        if not (math.isnan(forward) and math.isnan(reverse))
    ]
    assert survivors == primers, (
        "the second run dropped primers from the Gini gate; the reused position "
        f"map was short and the surviving pool was {survivors}"
    )


def test_overwrite_still_rescans_everything(genome, monkeypatch):
    """`overwrite=True` takes the fallback path and must not be short-circuited."""
    ss.get_positions([A12], [genome["prefix"]], [genome["fasta"]], circular=False)

    calls = _record_scans(monkeypatch)
    result = ss.get_positions(
        [A12], [genome["prefix"]], [genome["fasta"]], circular=False, overwrite=True
    )

    assert calls == [], "overwrite=True must not reach the Aho-Corasick branch"
    assert result is None, "the fallback path returns None by contract"


# ---------------------------------------------------------------------------
# Reuse is only valid under the conditions the file was scanned under.
# ---------------------------------------------------------------------------


def _write_fasta(path, sequence):
    path.write_text(">g\n" + sequence + "\n")
    return str(path)


def _assembly(copies, length=12_000, seed=23):
    """A genome carrying `copies` evenly spaced occurrences of A12."""
    rng = random.Random(seed)
    seq = list("".join(rng.choice("ACGT") for _ in range(length)))
    offsets = [1000 + 3000 * n for n in range(copies)]
    for offset in offsets:
        seq[offset : offset + len(A12)] = list(A12)
    return "".join(seq), offsets


def test_a_position_file_scanned_from_another_genome_is_not_reused(tmp_path):
    """Repointing the genome must not return the previous assembly's sites.

    Before reuse existed this branch rescanned every pattern on every run, so a
    stale file was irrelevant. Reuse made it authoritative, and the position
    files carry no filename of their own: with a single foreground genome the
    prefix comes from params.json and is independent of the genome path, so
    pointing `fg_genomes` at a new assembly and leaving the prefix alone reuses
    the same files by construction.
    """
    ss.clear_genome_cache()
    v1, _ = _assembly(copies=1)
    v2, expected = _assembly(copies=4)
    prefix = str(tmp_path / "fg")
    fasta_v1 = _write_fasta(tmp_path / "v1.fasta", v1)
    fasta_v2 = _write_fasta(tmp_path / "v2.fasta", v2)

    ss.get_positions([A12], [prefix], [fasta_v1], circular=False)
    ss.clear_genome_cache()
    repointed = ss.get_positions([A12], [prefix], [fasta_v2], circular=False)
    ss.clear_genome_cache()

    assert repointed[(prefix, A12)] == expected, (
        "the run reported the previous assembly's binding sites; a position file "
        "may only be reused for the genome it was scanned from"
    )


def test_toggling_circular_is_not_reused(tmp_path):
    """A pool scanned linearly must not be reused for a circular run.

    `get_all_positions_multi_k` appends `sequence[:max_k - 1]` before scanning,
    so a primer straddling the origin is found only when `circular` is true.
    `fg_circular` and `bg_circular` are params.json fields and plasmids are
    exactly this case.
    """
    ss.clear_genome_cache()
    rng = random.Random(5)
    seq = list("".join(rng.choice("ACGT") for _ in range(4_000)))
    split = 6
    seq[-split:] = list(A12[:split])
    seq[:split] = list(A12[split:])
    prefix = str(tmp_path / "fg")
    fasta = _write_fasta(tmp_path / "g.fasta", "".join(seq))

    linear = ss.get_positions([A12], [prefix], [fasta], circular=False)
    ss.clear_genome_cache()
    circular = ss.get_positions([A12], [prefix], [fasta], circular=True)
    ss.clear_genome_cache()

    assert linear[(prefix, A12)] == []
    assert circular[(prefix, A12)] == [4_000 - split], (
        "the circular run reused a linear scan and lost the site spanning the "
        "origin; circular is a scan parameter, not a reporting option"
    )


def test_a_position_file_without_a_provenance_record_is_rescanned(tmp_path):
    """An unrecorded file is unknown, not trustworthy.

    Every position file on disk before this guard existed has no record beside
    it, and the k-mer tables set the policy: absent means rescan once and write
    the record, rather than reuse.
    """
    ss.clear_genome_cache()
    v1, _ = _assembly(copies=1)
    v2, expected = _assembly(copies=4)
    prefix = str(tmp_path / "fg")
    fasta_v1 = _write_fasta(tmp_path / "v1.fasta", v1)
    fasta_v2 = _write_fasta(tmp_path / "v2.fasta", v2)

    ss.get_positions([A12], [prefix], [fasta_v1], circular=False)
    ss.clear_position_provenance(prefix, len(A12))
    ss.clear_genome_cache()

    result = ss.get_positions([A12], [prefix], [fasta_v2], circular=False)
    ss.clear_genome_cache()

    assert result[(prefix, A12)] == expected, (
        "a file with no provenance record was reused; it cannot be matched to "
        "any genome and must be rescanned"
    )


def test_datasets_from_a_replaced_genome_do_not_survive_into_a_later_pool(tmp_path):
    """A rejected file is emptied, not merely bypassed.

    Bypassing rescans only the primers this run asks about, leaving the rest of
    the datasets in place. The file is then recorded as current for the new
    genome while still holding the old genome's positions for those primers, and
    a later run with a different candidate pool reuses them.
    """
    ss.clear_genome_cache()
    v1, _ = _assembly(copies=1)
    v2, expected = _assembly(copies=4)
    prefix = str(tmp_path / "fg")
    fasta_v1 = _write_fasta(tmp_path / "v1.fasta", v1)
    fasta_v2 = _write_fasta(tmp_path / "v2.fasta", v2)

    # Pool 1 against v1 records A12. Pool 2 against v2 does not ask about it.
    ss.get_positions([A12], [prefix], [fasta_v1], circular=False)
    ss.clear_genome_cache()
    ss.get_positions([B12], [prefix], [fasta_v2], circular=False)
    ss.clear_genome_cache()

    # Pool 3 against v2 asks about it again.
    result = ss.get_positions([A12], [prefix], [fasta_v2], circular=False)
    ss.clear_genome_cache()

    assert result[(prefix, A12)] == expected, (
        "a dataset written from the replaced genome survived and was reused "
        "once the candidate pool asked about it again"
    )


def test_the_fallback_path_leaves_no_record_for_the_reuse_path_to_trust(tmp_path):
    """`overwrite=True` writes without checking provenance, so it clears it.

    That path decides what to rescan with `check_which_primers_absent_in_h5py`,
    which consults no record, so the datasets it leaves behind are of unverified
    origin and must not be inherited by the guarded path.
    """
    ss.clear_genome_cache()
    v1, _ = _assembly(copies=1)
    prefix = str(tmp_path / "fg")
    fasta = _write_fasta(tmp_path / "v1.fasta", v1)

    ss.get_positions([A12], [prefix], [fasta], circular=False)
    ss.get_positions([A12], [prefix], [fasta], circular=False, overwrite=True)
    ss.clear_genome_cache()

    with h5py.File(ss.position_file_path(prefix, len(A12)), "r") as handle:
        assert "provenance_fingerprint" not in handle.attrs

    calls = []
    real = ss.get_all_positions_multi_k

    def spy(primer_lists_by_k, seq_fname, circular, chunk_size=None):
        calls.append(sum(len(v) for v in primer_lists_by_k.values()))
        return real(primer_lists_by_k, seq_fname, circular, chunk_size)

    ss.get_all_positions_multi_k = spy
    try:
        result = ss.get_positions([A12], [prefix], [fasta], circular=False)
    finally:
        ss.get_all_positions_multi_k = real
    ss.clear_genome_cache()

    assert calls == [2], f"expected a rescan after the fallback path; scans were {calls}"
    assert result[(prefix, A12)] == [1000]


def test_a_primer_that_binds_nowhere_is_still_recorded_as_scanned(genome):
    """Key presence must mean "scanned", not "found".

    The reuse gate decides a primer needs no rescan by asking whether the HDF5
    file holds a dataset for it. That is only sound because the scan writes an
    entry for every pattern it was handed, including the ones that occur
    nowhere: `write_to_h5py` is given `all_positions.get(p, [])` for each
    pattern, so a primer with no binding sites gets a key with an empty
    dataset.

    Nothing else enforces that. A future scan path that recorded only the
    primers it found would leave the absent ones without keys, the reuse gate
    would classify them as needing a scan forever, and the far worse case is
    the mirror: a caller that trusts key ABSENCE to mean "not scanned" while
    some other path has written only hits. Both the implementer and the
    reviewer of this change identified the invariant independently and neither
    found a test holding it, which is why this one exists.

    `ABSENT12` was chosen so that neither it nor its reverse complement occurs
    in the fixture genome; the assertions below check that rather than assume
    it.
    """
    absent = "CGCGATATCGCG"
    contents = open(genome["fasta"]).read()
    assert absent not in contents
    assert ss.reverse_complement(absent) not in contents

    result = ss.get_positions([A12, absent], [genome["prefix"]], [genome["fasta"]], circular=False)

    # The returned map is complete: the absent primer is present with no sites.
    assert result[(genome["prefix"], absent)] == []
    assert result[(genome["prefix"], A12)] == [1_000]

    # And the file records it, which is what the reuse gate reads.
    with h5py.File(ss.position_file_path(genome["prefix"], 12), "r") as handle:
        assert absent in handle, "a scanned primer with no sites must still be a key"
        assert len(handle[absent][:]) == 0


def test_the_same_genome_by_another_path_spelling_is_still_reused(tmp_path):
    """The fingerprint identifies the genome, so the path spelling must not.

    The record also carries the genome path, and comparing whole records made a
    path difference a mismatch: the file was rescanned in full and the warning
    blamed `circular`, printing the same value twice, because a path difference
    fell through to that arm. `os.path.abspath` does not resolve symlinks, so
    the spellings differ after a symlinked parent directory, an automount, or
    `/tmp` against `/private/tmp` on macOS.
    """
    ss.clear_genome_cache()
    v1, expected = _assembly(copies=1)
    real_dir = tmp_path / "real"
    real_dir.mkdir()
    link_dir = tmp_path / "link"
    link_dir.symlink_to(real_dir)

    fasta = _write_fasta(real_dir / "g.fasta", v1)
    prefix = str(real_dir / "fg")
    linked_fasta = str(link_dir / "g.fasta")

    ss.get_positions([A12], [prefix], [fasta], circular=False)
    ss.clear_genome_cache()

    calls = []
    real = ss.get_all_positions_multi_k

    def spy(primer_lists_by_k, seq_fname, circular, chunk_size=None):
        calls.append(sum(len(v) for v in primer_lists_by_k.values()))
        return real(primer_lists_by_k, seq_fname, circular, chunk_size)

    ss.get_all_positions_multi_k = spy
    try:
        result = ss.get_positions([A12], [prefix], [linked_fasta], circular=False)
    finally:
        ss.get_all_positions_multi_k = real
    ss.clear_genome_cache()

    assert calls == [], (
        "the same genome reached by a second path spelling was rescanned; the "
        f"path is not the identity, the fingerprint is. Scans were {calls}"
    )
    assert result[(prefix, A12)] == expected


def test_a_run_that_does_not_finish_leaves_the_previous_position_file_intact(tmp_path):
    """A rejected file is discarded when the replacement is in hand, not before.

    Discarding at detection time left an empty position file if the rescan did
    not complete. On the foreground the step-4 prerequisite check catches that
    and refuses; on the background nothing does, and a background file holding
    nothing reports `total_bg_sites: 0` and a perfect selectivity ratio, which
    is Known Issue 5's symptom by another route.
    """
    ss.clear_genome_cache()
    v1, before = _assembly(copies=1)
    v2, _ = _assembly(copies=4)
    prefix = str(tmp_path / "fg")
    fasta_v1 = _write_fasta(tmp_path / "v1.fasta", v1)
    fasta_v2 = _write_fasta(tmp_path / "v2.fasta", v2)

    ss.get_positions([A12], [prefix], [fasta_v1], circular=False)
    ss.clear_genome_cache()

    real = ss.get_all_positions_multi_k

    def die(primer_lists_by_k, seq_fname, circular, chunk_size=None):
        raise KeyboardInterrupt("interrupted mid-scan")

    ss.get_all_positions_multi_k = die
    try:
        with pytest.raises(KeyboardInterrupt):
            ss.get_positions([A12], [prefix], [fasta_v2], circular=False)
    finally:
        ss.get_all_positions_multi_k = real
    ss.clear_genome_cache()

    with h5py.File(ss.position_file_path(prefix, len(A12)), "r") as handle:
        assert A12 in handle, (
            "the interrupted run left an empty position file; the previous one "
            "is stale but complete, and an empty background file is read "
            "downstream as a panel that binds the host nowhere"
        )
        assert handle[A12][:].tolist() == before


def test_a_position_file_that_cannot_be_opened_is_replaced_rather_than_fatal(tmp_path):
    """An unreadable file is treated as nothing scanned, and then written over.

    The reuse gate always returned everything for scanning here, but the write
    that followed opened the same file "r+" and died on it. Discarding at write
    time covers this case for free.
    """
    ss.clear_genome_cache()
    v1, expected = _assembly(copies=1)
    prefix = str(tmp_path / "fg")
    fasta = _write_fasta(tmp_path / "v1.fasta", v1)

    with open(ss.position_file_path(prefix, len(A12)), "wb") as fh:
        fh.write(b"not an HDF5 file")

    result = ss.get_positions([A12], [prefix], [fasta], circular=False)
    ss.clear_genome_cache()

    assert result[(prefix, A12)] == expected
