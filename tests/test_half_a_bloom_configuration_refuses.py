"""`use_bloom_filter` without a path screened nothing, and said nothing.

`get_all_rates` warns when a path is set without the flag. The opposite
pairing was silent, and it is the damaging one:

    use_bloom_filter: true, bloom_filter_path: null
      -> the `use_bloom and bloom_path` branch is not taken
      -> falls through to get_rates_for_one_species(primer_list, bg_prefixes)
      -> parameter.py sets bg_seq_lengths = [] when use_bloom_filter is set
         and bg_prefixes is empty, which is the configuration this flag is
         documented to produce
      -> no tasks, an empty dict, every bg_count is None
      -> bg_bool = (bg_count is None) = True

Every candidate clears the background gate. The None-passes rule is deliberate
-- a k-mer missing from a jellyfish table may still have sites the string
search finds -- but here it is not one k-mer missing, it is the whole gate off.
"""

import pandas as pd
import pytest

from neoswga.core import filter as filter_mod
from neoswga.core import parameter
from neoswga.core.exceptions import InvalidDesignRequest


@pytest.fixture(autouse=True)
def gates(monkeypatch):
    """`get_all_rates` reads the frequency gates off the parameter module,
    which `get_params` populates. Seed them so these tests exercise the Bloom
    branch rather than an unrelated AttributeError."""
    monkeypatch.setattr(parameter, "min_fg_freq", 1e-5, raising=False)
    monkeypatch.setattr(parameter, "max_bg_freq", 5e-6, raising=False)


def test_the_flag_without_a_path_refuses(monkeypatch):
    monkeypatch.setattr(parameter, "use_bloom_filter", True, raising=False)
    monkeypatch.setattr(parameter, "bloom_filter_path", None, raising=False)

    with pytest.raises(InvalidDesignRequest) as excinfo:
        filter_mod.get_all_rates(
            ["ACGTACGTAC"],
            fg_prefixes=[],
            bg_prefixes=[],
            fg_total_length=1000,
            bg_total_length=1000,
        )
    message = str(excinfo.value)
    assert "bloom_filter_path" in message
    assert "build-filter" in message


def test_the_refusal_comes_before_any_counting(monkeypatch):
    """It must not depend on reaching the empty-prefix path to notice."""
    called = []
    monkeypatch.setattr(parameter, "use_bloom_filter", True, raising=False)
    monkeypatch.setattr(parameter, "bloom_filter_path", None, raising=False)
    monkeypatch.setattr(
        filter_mod,
        "get_rates_for_one_species",
        lambda *a, **kw: called.append(a) or {},
    )

    with pytest.raises(InvalidDesignRequest):
        filter_mod.get_all_rates(
            ["ACGTACGTAC"],
            fg_prefixes=["fg"],
            bg_prefixes=["bg"],
            fg_total_length=1000,
            bg_total_length=1000,
        )
    assert called == [], "counting started before the configuration was checked"


def test_neither_set_is_untouched(monkeypatch):
    """The default path must be unaffected by this check."""
    monkeypatch.setattr(parameter, "use_bloom_filter", False, raising=False)
    monkeypatch.setattr(parameter, "bloom_filter_path", None, raising=False)
    monkeypatch.setattr(filter_mod, "get_rates_for_one_species", lambda *a, **kw: {})

    df = filter_mod.get_all_rates(
        ["ACGTACGTAC"],
        fg_prefixes=[],
        bg_prefixes=[],
        fg_total_length=1000,
        bg_total_length=1000,
    )
    assert isinstance(df, pd.DataFrame)
    assert list(df["primer"]) == ["ACGTACGTAC"]


def test_a_path_without_the_flag_still_only_warns(monkeypatch, caplog):
    """That pairing is recoverable -- exact counting is a real answer -- so it
    keeps its warning rather than acquiring a refusal."""
    monkeypatch.setattr(parameter, "use_bloom_filter", False, raising=False)
    monkeypatch.setattr(parameter, "bloom_filter_path", "/nonexistent.pkl", raising=False)
    monkeypatch.setattr(filter_mod, "get_rates_for_one_species", lambda *a, **kw: {})

    with caplog.at_level("WARNING"):
        filter_mod.get_all_rates(
            ["ACGTACGTAC"],
            fg_prefixes=[],
            bg_prefixes=[],
            fg_total_length=1000,
            bg_total_length=1000,
        )
    assert "use_bloom_filter" in caplog.text
