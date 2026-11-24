"""`validate` reports both k-mer counters, and which one counting will use.

It printed `Jellyfish: NOT FOUND` on a machine with KMC3 alone, which counts
fine and is the preferred counter (2026-11-23 sweep, tier 3).
"""

import pytest

from neoswga.cli.setup import _counter_capability_lines
from neoswga.core import kmer_backend


@pytest.mark.parametrize(
    "kmc,jellyfish,used",
    [(True, False, "kmc"), (False, True, "jellyfish"), (True, True, "kmc")],
)
def test_the_counter_lines_name_both_and_the_one_used(monkeypatch, kmc, jellyfish, used):
    monkeypatch.setattr(kmer_backend.KmcBackend, "available", lambda self: kmc)
    monkeypatch.setattr(kmer_backend.JellyfishBackend, "available", lambda self: jellyfish)
    text = "\n".join(_counter_capability_lines())
    assert "KMC" in text and "Jellyfish" in text
    assert f"will use: {used}" in text
    assert "NOT FOUND" not in text


def test_no_counter_says_counting_cannot_run(monkeypatch):
    monkeypatch.setattr(kmer_backend.KmcBackend, "available", lambda self: False)
    monkeypatch.setattr(kmer_backend.JellyfishBackend, "available", lambda self: False)
    assert "cannot run" in "\n".join(_counter_capability_lines())
