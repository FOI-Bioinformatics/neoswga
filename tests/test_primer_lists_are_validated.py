"""Primer lists are validated where they are read.

`analyze-set --primers ACGTXXQQ` printed a Tm, a free energy and a GC content
for a sequence that is not DNA, and `swap-primer`, `contract-set` and
`rescore-set` read their lists through `_load_primer_list`, which upper-cased
and checked nothing, so a duplicate was kept twice in the set it analysed.
`analyze-dimers` and `analyze-stability` already refused the same input
(2026-11-23 sweep, tier 2).
"""

import sys

import pytest

from neoswga import cli_unified
from neoswga.cli.iterate import _load_primer_list


def test_analyze_set_refuses_a_sequence_that_is_not_dna(monkeypatch, capsys):
    monkeypatch.setattr(sys, "argv", ["neoswga", "analyze-set", "--primers", "ACGTXXQQ"])
    with pytest.raises(SystemExit) as info:
        cli_unified.main()
    assert info.value.code == 1
    assert "Tm:" not in capsys.readouterr().out


def test_the_set_editing_loader_refuses_a_sequence_that_is_not_dna():
    with pytest.raises(ValueError, match="ACGTXXQQ"):
        _load_primer_list(["ACGTACGTAC", "acgtxxqq"], None)


def test_the_set_editing_loader_drops_a_duplicate():
    assert _load_primer_list(["ACGTACGTAC", "acgtacgtac", "TTGACCATGA"], None) == [
        "ACGTACGTAC",
        "TTGACCATGA",
    ]
