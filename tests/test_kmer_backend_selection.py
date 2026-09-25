"""The k-mer counter is chosen, and the choice is named.

KMC3 becomes the default and jellyfish stays selectable. The parity flags are
the substance of this module rather than an implementation detail: KMC's
defaults compute a different quantity from jellyfish's, and both differences
are failures this repository already carries.

  `-ci2`   excludes k-mers occurring once. Measured on the shipped genomes
           that is 96.9% of wMel's 18-mers and 95.6% of Drosophila's. A
           background counted that way reports a host as almost k-mer-free and
           every candidate as specific, which is the silent-zero shape of
           Known Issues 5, 6, 13 and 15.
  `-cs255` saturates the counter at 255, so a host k-mer occurring 10,000
           times reads 255. That is Known Issue 7 exactly, where an integer
           ceiling hid a primer's host load and changed the delivered panel
           once the true figure was visible.

Measurements: docs/validation/kmer_counter_comparison_2026-09-25.md
"""

import inspect

import pytest

from neoswga.core import kmer_backend


def test_kmc_is_preferred_when_unset_and_installed(monkeypatch):
    from neoswga.core import parameter

    monkeypatch.setattr(parameter, "kmer_counter", None, raising=False)
    monkeypatch.setattr(kmer_backend.KmcBackend, "available", lambda self: True)
    assert kmer_backend.select_backend().name == "kmc"


def test_unset_falls_back_to_jellyfish_when_kmc_is_absent(monkeypatch):
    """Both produce identical tables, so the fallback changes time and memory,
    never a result -- which is what makes it acceptable to do unasked."""
    from neoswga.core import parameter

    monkeypatch.setattr(parameter, "kmer_counter", None, raising=False)
    monkeypatch.setattr(kmer_backend.KmcBackend, "available", lambda self: False)
    monkeypatch.setattr(kmer_backend.JellyfishBackend, "available", lambda self: True)
    assert kmer_backend.select_backend().name == "jellyfish"


def test_an_explicit_setting_is_not_second_guessed(monkeypatch):
    """Named in params.json, KMC is returned even when absent, so the caller's
    availability check refuses with the install command rather than quietly
    running jellyfish."""
    from neoswga.core import parameter

    monkeypatch.setattr(parameter, "kmer_counter", "kmc", raising=False)
    monkeypatch.setattr(kmer_backend.KmcBackend, "available", lambda self: False)
    assert kmer_backend.select_backend().name == "kmc"


def test_jellyfish_is_selectable():
    assert kmer_backend.select_backend("jellyfish").name == "jellyfish"


def test_an_unknown_backend_is_refused_and_the_message_names_the_real_ones():
    with pytest.raises(ValueError) as excinfo:
        kmer_backend.select_backend("dsk")
    message = str(excinfo.value)
    assert "dsk" in message
    assert "kmc" in message and "jellyfish" in message


def test_a_missing_backend_names_how_to_install_it(monkeypatch):
    """The default changes, so a machine carrying only jellyfish must be told
    what to install rather than getting a bare FileNotFoundError."""
    backend = kmer_backend.select_backend("kmc")
    monkeypatch.setattr(backend, "available", lambda: False)
    with pytest.raises(Exception) as excinfo:
        backend.require_available()
    message = str(excinfo.value)
    assert "kmc" in message
    assert "conda" in message or "install" in message


def test_kmc_counts_singletons_and_does_not_saturate():
    """Neither KMC default may come back. See the module docstring."""
    source = inspect.getsource(kmer_backend)
    assert "-ci1" in source, "without -ci1 KMC excludes k-mers occurring once"
    assert "-cs255" not in source, "the default counter ceiling saturates at 255"


def test_jellyfish_still_counts_canonically():
    """Jellyfish needs -C; KMC is canonical unless -b is passed."""
    source = inspect.getsource(kmer_backend)
    assert '"-C"' in source
    assert '"-b"' not in source, "-b would turn KMC's canonical form off"
