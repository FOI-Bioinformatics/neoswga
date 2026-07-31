"""Per-primer validation against measured plasmid amplification.

Data: Dwivedi-Yu et al. (2023) PLOS Comput Biol 19(4):e1010137, S1 Table -- 310
qPCR measurements of fold amplification for 284 primers on the pcDNA and pLTR
plasmids. Those are the same two plasmids shipped in `examples/plasmid_example`,
so the binding sites can be recomputed with this codebase rather than trusted
from the paper.

**What this dataset can and cannot test.** Every primer in it has exactly one
on-target binding site and zero off-target sites -- they were selected that way.
Site geometry is therefore *constant*, and no amount of analysis can make this
dataset speak to whether foreground/background site ratio predicts selectivity:
the predictor has no variance.

That is worth stating because the naive analyses are actively misleading here. A
Spearman correlation of site ratio against measured selectivity returns about
-0.01, which reads as "no relationship" but actually means "no variance". And
splitting into deciles by site ratio appears to show a four-fold effect, which is
purely an artefact of how ties get ordered.

What the dataset *does* test is whether sequence-level properties predict
measured amplification when site geometry is held fixed. They do: primer length
and melting temperature both track measured fold amplification. That is a direct
check on the premise behind the pipeline's thermodynamic filtering.
"""

import json
from pathlib import Path

import pytest

DATA = Path(__file__).parent / "data" / "dwivedi_yu_2023_plasmid_amplification.json"
PLASMIDS = Path(__file__).resolve().parents[2] / "examples" / "plasmid_example"


@pytest.fixture(scope="module")
def measurements():
    return json.loads(DATA.read_text())["measurements"]


def _spearman(x, y):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0] * len(v)
        for pos, i in enumerate(order):
            r[i] = pos + 1
        return r

    rx, ry = rank(x), rank(y)
    n = len(x)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = (sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry)) ** 0.5
    return num / den if den else 0.0


# ----------------------------------------------------------------------
# Fixture integrity
# ----------------------------------------------------------------------


def test_fixture_is_intact(measurements):
    assert len(measurements) == 310
    assert len({m["primer"] for m in measurements}) == 284  # some tested at 2 molarities
    for m in measurements:
        assert set(m["primer"]) <= set("ACGT")
        assert 6 <= len(m["primer"]) <= 12


def test_site_geometry_is_constant_by_construction(measurements):
    """The premise of every other test in this file.

    If this ever fails, the dataset has gained variance in site counts and can
    then legitimately be used to test selectivity prediction.
    """
    geometry = {(m["sites_pcDNA"], m["sites_pLTR"]) for m in measurements}
    on_off = {
        (
            m["sites_pcDNA"] if m["designed_for"] == "pcDNA" else m["sites_pLTR"],
            m["sites_pLTR"] if m["designed_for"] == "pcDNA" else m["sites_pcDNA"],
        )
        for m in measurements
    }
    assert all(off == 0 for _, off in on_off), "some primer binds its off-target plasmid"
    assert all(on >= 1 for on, _ in on_off), "some primer does not bind its own target"
    assert len(geometry) <= 4, f"site geometry gained variance: {sorted(geometry)}"


# ----------------------------------------------------------------------
# Our scanner against the shipped plasmids
# ----------------------------------------------------------------------


def _naive_count(primer, sequence, circular=True):
    rc = primer.translate(str.maketrans("ACGT", "TGCA"))[::-1]
    text = sequence + sequence[: len(primer) - 1] if circular else sequence
    total = 0
    for pattern in {primer, rc}:
        i = text.find(pattern)
        while i != -1:
            total += 1
            i = text.find(pattern, i + 1)
    return total


def _load(name):
    path = PLASMIDS / f"{name}.fasta"
    return "".join(
        line.strip() for line in path.read_text().splitlines() if not line.startswith(">")
    ).upper()


@pytest.mark.skipif(not PLASMIDS.exists(), reason="example plasmids not present")
def test_scanner_agrees_with_a_naive_search(measurements):
    """Guards the reverse-strand bug that made this benchmark necessary.

    `string_search.get_all_positions_per_k` scans the forward strand only, so
    the reverse complement has to be queried as a separate k-mer. The scan
    fallback in PositionCache originally omitted that, which silently zeroed
    every reverse-strand binding site -- 147 of 568 counts on this dataset,
    each reading as "this primer does not bind" rather than "we only looked at
    one strand".
    """
    from neoswga.core.position_cache import PositionCache

    sequences = {"pcDNA": _load("pcDNA"), "pLTR": _load("pLTR")}
    primers = sorted({m["primer"] for m in measurements})

    for name, sequence in sequences.items():
        cache = PositionCache(
            [f"/tmp/_validation_{name}"],
            primers,
            genome_paths=[str(PLASMIDS / f"{name}.fasta")],
            circular=True,
            on_missing="scan",
        )
        for primer in primers:
            ours = len(cache.get_positions(f"/tmp/_validation_{name}", primer))
            assert ours == _naive_count(primer, sequence), f"{name}/{primer}"


@pytest.mark.skipif(not PLASMIDS.exists(), reason="example plasmids not present")
def test_every_primer_binds_its_designed_target(measurements):
    """A primer designed against a plasmid must be found in that plasmid."""
    sequences = {"pcDNA": _load("pcDNA"), "pLTR": _load("pLTR")}
    for m in measurements:
        assert _naive_count(m["primer"], sequences[m["designed_for"]]) >= 1, m["name"]


# ----------------------------------------------------------------------
# Do sequence-level properties predict measured amplification?
# ----------------------------------------------------------------------


@pytest.fixture(scope="module")
def features(measurements):
    from neoswga.core import thermodynamics as thermo
    from neoswga.core.reaction_conditions import ReactionConditions

    conditions = ReactionConditions(temp=30.0, polymerase="phi29")
    return {
        "tm": [conditions.calculate_effective_tm(m["primer"]) for m in measurements],
        "gc": [thermo.gc_content(m["primer"]) for m in measurements],
        "length": [len(m["primer"]) for m in measurements],
        "amp": [m["on_target_amp"] for m in measurements],
    }


def test_primer_length_predicts_amplification(features):
    rho = _spearman(features["length"], features["amp"])
    assert rho > 0.4, f"expected a clear positive relationship, got {rho:+.3f}"


def test_melting_temperature_predicts_amplification(features):
    """Direct check on the premise behind thermodynamic filtering."""
    rho = _spearman(features["tm"], features["amp"])
    assert rho > 0.3, f"expected Tm to track amplification, got {rho:+.3f}"


def test_site_ratio_cannot_be_evaluated_on_this_dataset(measurements):
    """Documents why the obvious analysis is not run here.

    With site geometry constant, a rank correlation of site ratio against
    measured selectivity is dominated by ties and returns ~0. That reads as
    "site ratio does not predict selectivity", which would be a false and
    consequential conclusion drawn from a predictor with no variance.
    """
    ratios = {(m["sites_pcDNA"] + 0.5) / (m["sites_pLTR"] + 0.5) for m in measurements}
    assert len(ratios) <= 4, (
        "site ratio has gained variance; this dataset may now be usable for "
        "selectivity validation and this test should be revisited"
    )
