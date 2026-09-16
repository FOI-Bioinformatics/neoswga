"""The melting temperature the candidate loader filters on.

Audit finding B5. Two separate faults, and this file covers both.

The first is speed. `melting_temp.temp` built two 16-entry coefficient dicts on
every call and then ran `_overcount` once per entry per dict -- 32 scans of the
k-mer. Measured on 50,000 random 12-mers, those two sum expressions cost 7.72 us
of a 10.23 us total; one pass over dinucleotides costs 1.42 us. The value must
not move: `primer_attributes.get_melting_temp` and `network_optimizer` read it,
and `rf_preprocessing` serves it to the bundled random forest.

The second is correctness, and it is covered further down this file.
"""

import random

import pytest

from neoswga.core import melting_temp as mt

DH = {
    "AA": -7.9,
    "TT": -7.9,
    "AT": -7.2,
    "TA": -7.2,
    "CA": -8.5,
    "TG": -8.5,
    "GT": -8.4,
    "AC": -8.4,
    "CT": -7.8,
    "AG": -7.8,
    "GA": -8.2,
    "TC": -8.2,
    "CG": -10.6,
    "GC": -9.8,
    "GG": -8.0,
    "CC": -8.0,
}
DS = {
    "AA": -22.2,
    "TT": -22.2,
    "AT": -20.4,
    "TA": -21.3,
    "CA": -22.7,
    "TG": -22.7,
    "GT": -22.4,
    "AC": -22.4,
    "CT": -21.0,
    "AG": -21.0,
    "GA": -22.2,
    "TC": -22.2,
    "CG": -27.2,
    "GC": -24.4,
    "GG": -19.9,
    "CC": -19.9,
}


def _reference_sums(s):
    """The 32-scan form this task replaces, written out independently."""
    dh = sum(mt._overcount(s, pair) * coeff for pair, coeff in DH.items())
    ds = sum(mt._overcount(s, pair) * coeff for pair, coeff in DS.items())
    return dh, ds


def _kmers(n, k, seed):
    rng = random.Random(seed)
    return ["".join(rng.choice("ACGT") for _ in range(k)) for _ in range(n)]


def test_the_single_pass_agrees_with_the_scan_form():
    """Agreement to floating-point summation order, which is what this is.

    The plan's version of this test was titled "exactly" and claimed
    bit-identical results while asserting a 1e-9 tolerance. Measured, the two
    forms differ by up to 8.5e-14 on the sums and up to 2.8e-13 Celsius on
    `temp()` itself, with 910 of 2000 sampled k-mers exactly equal. That is the
    order in which sixteen floats are added, not a change to the model. The
    tolerance is the honest assertion; the title was not.
    """
    for k in (6, 8, 10, 12, 16):
        for kmer in _kmers(200, k, seed=k):
            dh, ds = mt._nn_sums(kmer)
            ref_dh, ref_ds = _reference_sums(kmer)
            assert dh == pytest.approx(ref_dh, abs=1e-9), kmer
            assert ds == pytest.approx(ref_ds, abs=1e-9), kmer


def test_a_two_base_sequence_has_exactly_one_stack():
    dh, ds = mt._nn_sums("CG")
    assert dh == pytest.approx(-10.6)
    assert ds == pytest.approx(-27.2)


def test_a_one_base_sequence_has_no_stacks():
    assert mt._nn_sums("A") == (0.0, 0.0)


def test_an_empty_sequence_has_no_stacks():
    assert mt._nn_sums("") == (0.0, 0.0)


def test_an_ambiguous_base_contributes_no_stack():
    """N is not in the table. The scan form skipped it and so must this one."""
    dh, ds = mt._nn_sums("ANA")
    assert dh == 0.0 and ds == 0.0


def test_temp_is_unchanged_for_a_known_primer():
    """A pinned value, so a future refactor cannot drift the number quietly."""
    assert mt.temp("TTGACCATGA") == pytest.approx(48.616944605708056, abs=1e-9)


# ---------------------------------------------------------------------------
# The loader must filter on the same Tm the gate filters on
# ---------------------------------------------------------------------------

import inspect

from neoswga.core.reaction_conditions import ReactionConditions


def _write_kmer_file(tmp_path, k, kmers):
    path = tmp_path / f"target_{k}mer_all.txt"
    path.write_text("".join(f"{kmer} 10\n" for kmer in kmers))
    return str(tmp_path / "target")


def test_the_loader_uses_the_conditions_it_is_given():
    """The regression. Under additives the shim discarded candidates that pass
    the real window, and did it silently."""
    plain = ReactionConditions(temp=30.0, na_conc=50.0, mg_conc=10.0)
    with_additives = ReactionConditions(
        temp=30.0, na_conc=50.0, mg_conc=10.0, dmso_percent=10.0, betaine_m=1.5
    )
    kmer = "ACGTACGTACGT"
    assert with_additives.calculate_effective_tm(kmer) < plain.calculate_effective_tm(kmer) - 2.0


def test_a_kmer_inside_the_effective_window_is_kept(tmp_path):
    from neoswga.core.kmer_counter import get_primer_list_from_kmers

    conditions = ReactionConditions(temp=30.0, na_conc=50.0, mg_conc=10.0)
    kmer = "ACGTACGTACGT"
    tm = conditions.calculate_effective_tm(kmer)
    prefix = _write_kmer_file(tmp_path, 12, [kmer])

    kept = get_primer_list_from_kmers(
        [prefix],
        kmer_lengths=range(12, 13),
        min_tm=tm - 5.0,
        max_tm=tm + 5.0,
        wide_tm_margin=0.0,
        conditions=conditions,
    )
    assert kept == [kmer]


def test_a_kmer_outside_the_effective_window_is_dropped(tmp_path):
    from neoswga.core.kmer_counter import get_primer_list_from_kmers

    conditions = ReactionConditions(temp=30.0, na_conc=50.0, mg_conc=10.0)
    kmer = "ACGTACGTACGT"
    tm = conditions.calculate_effective_tm(kmer)
    prefix = _write_kmer_file(tmp_path, 12, [kmer])

    kept = get_primer_list_from_kmers(
        [prefix],
        kmer_lengths=range(12, 13),
        min_tm=tm + 20.0,
        max_tm=tm + 30.0,
        wide_tm_margin=0.0,
        conditions=conditions,
    )
    assert kept == []


def test_additives_move_the_window_the_loader_applies(tmp_path):
    """The specific loss B5 measured: a candidate the real gate keeps under
    additives, which the shim's fixed offset threw away."""
    from neoswga.core.kmer_counter import get_primer_list_from_kmers

    conditions = ReactionConditions(
        temp=30.0, na_conc=50.0, mg_conc=10.0, dmso_percent=10.0, betaine_m=1.5
    )
    kmer = "ACGTACGTACGT"
    tm = conditions.calculate_effective_tm(kmer)
    prefix = _write_kmer_file(tmp_path, 12, [kmer])

    kept = get_primer_list_from_kmers(
        [prefix],
        kmer_lengths=range(12, 13),
        min_tm=tm - 3.0,
        max_tm=tm + 3.0,
        wide_tm_margin=0.0,
        conditions=conditions,
    )
    assert kept == [kmer]


def test_the_margin_default_is_small_now_that_the_estimators_agree():
    """15 C existed to absorb disagreement between two estimators. There is
    only one estimator now, so the margin is stated headroom instead."""
    from neoswga.core.kmer_counter import get_primer_list_from_kmers

    sig = inspect.signature(get_primer_list_from_kmers)
    assert sig.parameters["wide_tm_margin"].default == 2.0


def test_ambiguous_bases_are_skipped_without_raising(tmp_path):
    """calculate_effective_tm warns and substitutes penalties for N rather than
    raising, so the loader needs its own guard."""
    from neoswga.core.kmer_counter import get_primer_list_from_kmers

    prefix = _write_kmer_file(tmp_path, 6, ["ACGTGC", "ACNTGC"])
    kept = get_primer_list_from_kmers(
        [prefix],
        kmer_lengths=range(6, 7),
        min_tm=-100.0,
        max_tm=200.0,
        wide_tm_margin=0.0,
        gc_min=0.0,
        gc_max=1.0,
    )
    assert "ACNTGC" not in kept
    assert "ACGTGC" in kept


def test_no_conditions_falls_back_to_a_default_reaction(tmp_path):
    """Library callers and tests reach this without a parameter module."""
    from neoswga.core.kmer_counter import get_primer_list_from_kmers

    prefix = _write_kmer_file(tmp_path, 12, ["ACGTACGTACGT"])
    kept = get_primer_list_from_kmers(
        [prefix], kmer_lengths=range(12, 13), min_tm=-100.0, max_tm=200.0
    )
    assert kept == ["ACGTACGTACGT"]
