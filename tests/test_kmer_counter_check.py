"""`count-kmers` and `design` accept whichever k-mer counter will run.

The pre-flight check asked only for `jellyfish`, so a machine with KMC3 alone
was refused before the backend that prefers KMC was ever consulted
(2026-11-23 sweep, tier 1 item 8).
"""

import json
import shutil

import pytest

from neoswga.cli._common import check_kmer_counter_available


def _installed(monkeypatch, *names):
    real = set(names)
    monkeypatch.setattr(shutil, "which", lambda name: f"/bin/{name}" if name in real else None)
    for var in ("KMC_BIN", "JELLYFISH_BIN"):
        monkeypatch.delenv(var, raising=False)


def _params(tmp_path, **keys):
    path = tmp_path / "params.json"
    path.write_text(json.dumps(keys))
    return str(path)


def test_kmc_only_passes(tmp_path, monkeypatch):
    _installed(monkeypatch, "kmc", "kmc_dump")
    check_kmer_counter_available(_params(tmp_path))


def test_jellyfish_only_passes(tmp_path, monkeypatch):
    _installed(monkeypatch, "jellyfish")
    monkeypatch.setattr("neoswga.core.kmer_counter.get_jellyfish_version", lambda: "2.3.1")
    check_kmer_counter_available(_params(tmp_path))


def test_named_counter_missing_exits(tmp_path, monkeypatch):
    _installed(monkeypatch, "kmc", "kmc_dump")
    with pytest.raises(SystemExit):
        check_kmer_counter_available(_params(tmp_path, kmer_counter="jellyfish"))


def test_blacklist_needs_jellyfish(tmp_path, monkeypatch):
    _installed(monkeypatch, "kmc", "kmc_dump")
    with pytest.raises(SystemExit):
        check_kmer_counter_available(_params(tmp_path), need_jellyfish=True)


def test_neither_installed_exits(tmp_path, monkeypatch):
    _installed(monkeypatch)
    with pytest.raises(SystemExit):
        check_kmer_counter_available(_params(tmp_path))


def test_jellyfish_one_is_still_refused(tmp_path, monkeypatch):
    _installed(monkeypatch, "jellyfish")
    monkeypatch.setattr("neoswga.core.kmer_counter.get_jellyfish_version", lambda: "1.1.11")
    with pytest.raises(SystemExit):
        check_kmer_counter_available(_params(tmp_path))
