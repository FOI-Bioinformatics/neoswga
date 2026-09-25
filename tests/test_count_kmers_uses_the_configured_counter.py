"""count-kmers runs the counter params.json names.

`kmer_counter` was declared in the schema, validated, documented and given a
default, and `count-kmers` ignored it: every call site ran jellyfish directly.
That is Known Issue 8's class, on the branch that introduced the key, and the
tests covering the backend could not see it because they called the backend
themselves rather than walking from the step that counts.

These tests walk from `run_jellyfish` and `MultiGenomeKmerCounter`, the two
counting entry points, and replace only the tools. They need neither KMC nor
jellyfish installed, so CI runs them.
"""

import json

import pytest

from neoswga.core import kmer_counter, kmer_tables, parameter
from neoswga.core.kmer_backend import JellyfishBackend, KmcBackend, select_backend


@pytest.fixture
def genome(tmp_path):
    fasta = tmp_path / "g.fna"
    fasta.write_text(">g\n" + "ACGTTGCAAGGCTTAC" * 20 + "\n")
    return str(fasta)


@pytest.fixture
def fake_kmc(monkeypatch):
    calls = []

    def count(self, genome_fname, k, out_prefix, threads=4):
        calls.append(k)
        stem = self.database(out_prefix, k)
        for suffix in (".kmc_pre", ".kmc_suf"):
            with open(stem + suffix, "wb") as fh:
                fh.write(b"x")

    monkeypatch.setattr(KmcBackend, "count", count)
    monkeypatch.setattr(KmcBackend, "available", lambda self: True)
    return calls


@pytest.fixture
def fake_jellyfish(monkeypatch):
    calls = []

    def one_k(output_prefix, genome_fname, k, cpus, hash_size):
        calls.append(k)
        with open(kmer_tables.text_table_path(output_prefix, k), "w") as fh:
            fh.write("ACGTTGCA 3\n")

    monkeypatch.setattr(kmer_counter, "_count_one_k_with_jellyfish", one_k)
    monkeypatch.setattr(kmer_counter, "require_jellyfish", lambda: None)
    monkeypatch.setattr(JellyfishBackend, "available", lambda self: True)
    return calls


def _configure(monkeypatch, value):
    monkeypatch.setattr(parameter, "kmer_counter", value, raising=False)


def test_an_explicit_kmc_runs_kmc(monkeypatch, tmp_path, genome, fake_kmc, fake_jellyfish):
    _configure(monkeypatch, "kmc")
    kmer_counter.run_jellyfish(genome, str(tmp_path / "out"), min_k=6, max_k=8, cpus=1)
    assert sorted(fake_kmc) == [6, 7, 8]
    assert fake_jellyfish == []


def test_an_explicit_jellyfish_runs_jellyfish(
    monkeypatch, tmp_path, genome, fake_kmc, fake_jellyfish
):
    _configure(monkeypatch, "jellyfish")
    kmer_counter.run_jellyfish(genome, str(tmp_path / "out"), min_k=6, max_k=8, cpus=1)
    assert sorted(fake_jellyfish) == [6, 7, 8]
    assert fake_kmc == []


def test_unset_prefers_kmc_when_it_is_installed(
    monkeypatch, tmp_path, genome, fake_kmc, fake_jellyfish
):
    _configure(monkeypatch, None)
    kmer_counter.run_jellyfish(genome, str(tmp_path / "out"), min_k=6, max_k=6, cpus=1)
    assert fake_kmc == [6] and fake_jellyfish == []


def test_unset_falls_back_to_jellyfish_without_kmc(
    monkeypatch, tmp_path, genome, fake_kmc, fake_jellyfish
):
    _configure(monkeypatch, None)
    monkeypatch.setattr(KmcBackend, "available", lambda self: False)
    kmer_counter.run_jellyfish(genome, str(tmp_path / "out"), min_k=6, max_k=6, cpus=1)
    assert fake_jellyfish == [6] and fake_kmc == []


def test_an_explicit_kmc_without_kmc_refuses(monkeypatch, tmp_path, genome, fake_jellyfish):
    """Asked for by name, it is required. Quietly running jellyfish instead
    would make the setting say one thing and do another."""
    _configure(monkeypatch, "kmc")
    monkeypatch.setattr(KmcBackend, "available", lambda self: False)
    with pytest.raises(RuntimeError) as excinfo:
        kmer_counter.run_jellyfish(genome, str(tmp_path / "out"), min_k=6, max_k=6, cpus=1)
    assert "conda install -c bioconda kmc" in str(excinfo.value)
    assert fake_jellyfish == []


def test_a_database_only_run_passes_the_output_check(monkeypatch, tmp_path, genome, fake_kmc):
    """The check after counting demanded text files, which a correct KMC run
    does not write."""
    _configure(monkeypatch, "kmc")
    prefix = str(tmp_path / "out")
    kmer_counter.run_jellyfish(genome, prefix, min_k=6, max_k=7, cpus=1)
    assert kmer_tables.table_exists(prefix, 6) and kmer_tables.table_exists(prefix, 7)
    assert not (tmp_path / "out_6mer_all.txt").exists()


def test_a_current_table_is_not_recounted(monkeypatch, tmp_path, genome, fake_kmc):
    _configure(monkeypatch, "kmc")
    prefix = str(tmp_path / "out")
    kmer_counter.run_jellyfish(genome, prefix, min_k=6, max_k=7, cpus=1)
    kmer_counter.run_jellyfish(genome, prefix, min_k=6, max_k=7, cpus=1)
    assert sorted(fake_kmc) == [6, 7]


def test_changing_the_counter_does_not_recount_an_existing_directory(
    monkeypatch, tmp_path, genome, fake_kmc, fake_jellyfish
):
    """Every existing data directory holds jellyfish text tables. Switching the
    counter must not throw them away: the counts are the same quantity."""
    prefix = str(tmp_path / "out")
    _configure(monkeypatch, "jellyfish")
    kmer_counter.run_jellyfish(genome, prefix, min_k=6, max_k=6, cpus=1)
    _configure(monkeypatch, "kmc")
    kmer_counter.run_jellyfish(genome, prefix, min_k=6, max_k=6, cpus=1)
    assert fake_jellyfish == [6] and fake_kmc == []


def test_provenance_records_which_counter_ran(monkeypatch, tmp_path, genome, fake_kmc):
    _configure(monkeypatch, "kmc")
    prefix = str(tmp_path / "out")
    kmer_counter.run_jellyfish(genome, prefix, min_k=6, max_k=6, cpus=1)
    with open(kmer_counter.table_provenance_path(prefix, 6)) as fh:
        assert json.load(fh)["counter"] == "kmc"


def test_kmc_runs_one_k_at_a_time(monkeypatch, tmp_path, genome, fake_kmc):
    """KMC is multithreaded and needs at least 2 GB per run. Running the k
    values concurrently, as the jellyfish path does, would multiply that by
    the number of k -- 14 GB for the default six-to-twelve range."""
    _configure(monkeypatch, "kmc")
    import concurrent.futures

    seen = []
    real = concurrent.futures.ThreadPoolExecutor

    def recording(max_workers=None, **kwargs):
        seen.append(max_workers)
        return real(max_workers=max_workers, **kwargs)

    monkeypatch.setattr(concurrent.futures, "ThreadPoolExecutor", recording)
    kmer_counter.run_jellyfish(genome, str(tmp_path / "out"), min_k=6, max_k=12, cpus=1)
    assert seen == [1]


def test_the_multi_genome_counter_uses_it_too(monkeypatch, tmp_path, genome, fake_kmc):
    """The second counting entry point, used by `multi-genome`, had its own
    jellyfish command line."""
    _configure(monkeypatch, "kmc")
    monkeypatch.setattr(kmer_tables, "iter_table", lambda prefix, k: iter([("ACGTTGCA", 3)]))
    counter = kmer_counter.MultiGenomeKmerCounter(cpus=1, output_dir=str(tmp_path))
    counter.add_genome("g", genome)
    assert counter.count_kmers("g", 8) == {"ACGTTGCA": 3}
    assert fake_kmc == [8]


def test_select_backend_is_what_the_counting_step_asks():
    """Pins the PATH: the counting module must consult the selection."""
    import inspect

    source = inspect.getsource(kmer_counter)
    assert "select_backend" in source
    assert select_backend("jellyfish").name == "jellyfish"
