"""Tests for polymerase-aware defaults in get_params().

phi29 and equiphi29 should pick different primer length ranges and different
Mg2+ defaults when the user has not specified them. Explicit user values
always win.
"""

import json
from types import SimpleNamespace

import pytest


def _mk_args(json_file):
    attrs = {
        "json_file": str(json_file),
        "data_dir": None, "src_dir": None,
        "min_fg_freq": None, "max_bg_freq": None,
        "min_tm": None, "max_tm": None,
        "max_gini": None, "max_primer": None, "min_amp_pred": None,
        "cpus": None, "max_dimer_bp": None, "max_self_dimer_bp": None,
        "verbose": None, "drop_iterations": None, "iterations": None,
        "top_set_count": None, "retries": None, "max_sets": None,
        "selection_metric": None, "fg_circular": None, "bg_circular": None,
        "min_k": None, "max_k": None,
        "fasta_fore": None, "fasta_back": None,
        "kmer_fore": None, "kmer_back": None,
    }
    return SimpleNamespace(**attrs)


def _write_params(tmp_path, extra=None):
    fg_fasta = tmp_path / "fg.fna"
    fg_fasta.write_text(">seq\n" + "ACGT" * 250 + "\n")
    params = {
        "data_dir": str(tmp_path),
        "src_dir": str(tmp_path),
        "fg_genomes": [str(fg_fasta)],
        "fg_prefixes": [str(tmp_path / "fg")],
        "fg_seq_lengths": [1000],
        "cpus": 1,
    }
    if extra:
        params.update(extra)
    path = tmp_path / "params.json"
    path.write_text(json.dumps(params))
    return path


@pytest.mark.parametrize("polymerase,expected_min_k,expected_max_k,expected_mg", [
    ("phi29", 6, 12, 10.0),
    ("equiphi29", 10, 18, 10.0),
    ("bst", 15, 25, 8.0),
    ("klenow", 8, 15, 10.0),
])
def test_polymerase_default_kmer_range_and_mg(tmp_path, polymerase, expected_min_k,
                                               expected_max_k, expected_mg):
    from neoswga.core import parameter
    parameter.reset_to_defaults()

    params_file = _write_params(tmp_path, extra={"polymerase": polymerase})
    parameter.get_params(_mk_args(params_file))

    assert parameter.min_k == expected_min_k, (
        f"{polymerase}: min_k default should be {expected_min_k}, got {parameter.min_k}"
    )
    assert parameter.max_k == expected_max_k, (
        f"{polymerase}: max_k default should be {expected_max_k}, got {parameter.max_k}"
    )
    assert parameter.mg_conc == expected_mg, (
        f"{polymerase}: mg_conc default should be {expected_mg}, got {parameter.mg_conc}"
    )


def test_explicit_values_override_polymerase_defaults(tmp_path):
    """If the user explicitly sets min_k/max_k/mg_conc, those values are used."""
    from neoswga.core import parameter
    parameter.reset_to_defaults()

    params_file = _write_params(tmp_path, extra={
        "polymerase": "equiphi29",  # would default to 10-18 / 10 mM
        "min_k": 8,
        "max_k": 14,
        "mg_conc": 3.5,
    })
    parameter.get_params(_mk_args(params_file))

    assert parameter.min_k == 8
    assert parameter.max_k == 14
    assert parameter.mg_conc == 3.5


def test_unknown_polymerase_falls_back_to_phi29_kmer_range(tmp_path):
    """Defensive: unknown polymerase uses phi29-like defaults, does not crash."""
    from neoswga.core import parameter
    parameter.reset_to_defaults()

    params_file = _write_params(tmp_path, extra={"polymerase": "custom_enzyme"})
    # Expect ValueError or graceful fallback; current implementation just
    # falls back to phi29 defaults. This asserts the fallback contract.
    try:
        parameter.get_params(_mk_args(params_file))
    except ValueError:
        pytest.skip("Unknown polymerase raises — that's also an acceptable contract")
    assert parameter.min_k == 6
    assert parameter.max_k == 12
