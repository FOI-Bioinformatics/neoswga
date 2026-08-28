"""Behavioural tests for `SwgaSimulator`.

This class validates a primer set by simulating amplification against target and
background, and `simulation_plots` and `simulation_report` are both built around
its `SimulationResult`. It sat at 20% coverage, and the two tests that existed
constructed it with `__new__` to bypass `__init__` -- which turned out to be
hiding the reason: its position loader could not read the position files this
pipeline writes, and said nothing about it.

Everything downstream of loading (coverage, uniformity, gaps, enrichment, the
composite score and the recommendation) is computed from those positions, so a
silent zero there is a confident-looking verdict derived from nothing.
"""

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")

from neoswga.core.swga_simulator import SimulationResult, SwgaSimulator
from neoswga.core.thermodynamics import reverse_complement

PRIMERS = ["TTGACCATGA", "GCATTACGGT", "AACCGGTTAC"]


@pytest.fixture
def genome(tmp_path):
    """A 20 kb genome with the primers planted at known positions."""
    import random

    rng = random.Random(31337)
    seq = list("".join(rng.choice("ACGT") for _ in range(20_000)))
    for offset, primer in enumerate(PRIMERS):
        for pos in (1_000, 7_000, 13_000):
            start = pos + offset * 100
            seq[start : start + len(primer)] = list(primer)
    sequence = "".join(seq)

    path = tmp_path / "target.fasta"
    path.write_text(
        ">target\n" + "\n".join(sequence[i : i + 70] for i in range(0, len(sequence), 70)) + "\n"
    )
    return {"path": str(path), "seq": sequence}


def _positions_of(sequence, pattern):
    out, i = [], sequence.find(pattern)
    while i != -1:
        out.append(i)
        i = sequence.find(pattern, i + 1)
    return out


@pytest.fixture
def flat_h5(tmp_path, genome):
    """The layout this pipeline actually writes: primer -> dataset of positions.

    Reverse-strand hits live under the reverse complement, exactly as
    `PositionCache` reads them.
    """
    path = tmp_path / "flat_positions.h5"
    with h5py.File(path, "w") as f:
        for primer in PRIMERS:
            for key in {primer, reverse_complement(primer)}:
                pos = _positions_of(genome["seq"], key)
                if pos:
                    f.create_dataset(key, data=np.array(pos, dtype=np.int32))
    return str(path)


@pytest.fixture
def nested_h5(tmp_path, genome):
    """The legacy layout the loader used to assume: genome -> primer -> '+'/'-'."""
    path = tmp_path / "nested_positions.h5"
    with h5py.File(path, "w") as f:
        group = f.create_group("target")
        for primer in PRIMERS:
            pg = group.create_group(primer)
            pg.create_dataset("+", data=np.array(_positions_of(genome["seq"], primer)))
            pg.create_dataset(
                "-", data=np.array(_positions_of(genome["seq"], reverse_complement(primer)))
            )
    return str(path)


def _simulator(genome, h5_path):
    return SwgaSimulator(
        primers=PRIMERS,
        fg_genome=genome["path"],
        bg_genome=genome["path"],
        fg_positions_h5=h5_path,
        bg_positions_h5=h5_path,
        bin_size=1_000,
    )


# ----------------------------------------------------------------------
# Position loading -- the bug everything else was hiding behind
# ----------------------------------------------------------------------


def test_reads_the_layout_this_pipeline_writes(genome, flat_h5):
    """The loader assumed a genome -> primer -> strand hierarchy that nothing
    produces. Against a real file it picked a PRIMER as the genome group, every
    lookup missed, and it returned zero positions without raising."""
    sim = _simulator(genome, flat_h5)

    total = sum(len(v["+"]) + len(v["-"]) for v in sim.fg_positions.values())
    assert total > 0, "no positions loaded from the pipeline's own HDF5 layout"

    for primer in PRIMERS:
        assert len(sim.fg_positions[primer]["+"]) == 3, primer


def test_reverse_strand_positions_come_from_the_reverse_complement(genome, flat_h5):
    """Forward hits are stored under the primer, reverse hits under its revcomp.

    Reading only the primer key would drop every reverse-strand binding site,
    the same blind spot found twice already in this audit.
    """
    sim = _simulator(genome, flat_h5)
    for primer in PRIMERS:
        expected = _positions_of(genome["seq"], reverse_complement(primer))
        assert sorted(sim.fg_positions[primer]["-"].tolist()) == sorted(expected), primer


def test_legacy_nested_layout_still_works(genome, nested_h5):
    """Accepting the old shape as well keeps any external file readable."""
    sim = _simulator(genome, nested_h5)
    for primer in PRIMERS:
        assert len(sim.fg_positions[primer]["+"]) == 3, primer


def test_both_layouts_agree(genome, flat_h5, nested_h5):
    flat = _simulator(genome, flat_h5)
    nested = _simulator(genome, nested_h5)
    for primer in PRIMERS:
        assert flat.fg_positions[primer]["+"].tolist() == nested.fg_positions[primer]["+"].tolist()


def test_a_primer_absent_from_the_file_gets_empty_arrays(genome, flat_h5):
    """Absent must be empty, not missing -- downstream code indexes these."""
    sim = SwgaSimulator(
        primers=PRIMERS + ["CGCGCGCGCG"],
        fg_genome=genome["path"],
        bg_genome=genome["path"],
        fg_positions_h5=flat_h5,
        bg_positions_h5=flat_h5,
        bin_size=1_000,
    )
    entry = sim.fg_positions["CGCGCGCGCG"]
    assert len(entry["+"]) == 0 and len(entry["-"]) == 0
    assert set(entry) == {"+", "-", "forward", "reverse"}


def test_strand_aliases_are_consistent(genome, flat_h5):
    """Both naming conventions are exposed; they must not disagree.

    CLAUDE.md warns that PositionCache uses 'forward'/'reverse' while this
    module uses '+'/'-', so the aliases exist to bridge them.
    """
    sim = _simulator(genome, flat_h5)
    for primer in PRIMERS:
        entry = sim.fg_positions[primer]
        assert entry["forward"].tolist() == entry["+"].tolist()
        assert entry["reverse"].tolist() == entry["-"].tolist()


# ----------------------------------------------------------------------
# The fast simulation
# ----------------------------------------------------------------------


def test_fast_simulation_produces_bounded_metrics(genome, flat_h5):
    sim = _simulator(genome, flat_h5)
    result = sim.simulate_fast()

    assert isinstance(result, SimulationResult)
    assert 0.0 <= result.target_coverage <= 1.0
    assert 0.0 <= result.target_uniformity <= 1.0
    assert 0.0 <= result.composite_score <= 1.0
    assert result.recommendation in {"EXCELLENT", "GOOD", "FAIR", "POOR"}
    assert result.runtime >= 0


def test_coverage_is_not_zero_for_a_binding_set(genome, flat_h5):
    """The failure mode the loader bug produced: a confident zero."""
    result = _simulator(genome, flat_h5).simulate_fast()
    assert result.target_coverage > 0


def test_result_serialises(genome, flat_h5):
    """`to_dict` feeds the report and the plots."""
    import json
    import math

    result = _simulator(genome, flat_h5).simulate_fast()
    payload = result.to_dict()

    assert payload["mode"] == "fast"
    non_finite = [k for k, v in payload.items() if isinstance(v, float) and not math.isfinite(v)]
    assert not non_finite, f"non-finite values would not survive JSON: {non_finite}"
    json.dumps(payload, default=str)


# ----------------------------------------------------------------------
# Scoring helpers
# ----------------------------------------------------------------------


def test_gini_of_uniform_coverage_is_zero(genome, flat_h5):
    sim = _simulator(genome, flat_h5)
    uniform = np.ones(10_000)
    assert sim._calculate_gini(uniform) == pytest.approx(0.0, abs=1e-9)


def test_gini_rises_with_unevenness(genome, flat_h5):
    sim = _simulator(genome, flat_h5)
    uniform = np.ones(10_000)

    lumpy = np.zeros(10_000)
    lumpy[:1_000] = 10.0

    assert sim._calculate_gini(lumpy) > sim._calculate_gini(uniform)


def test_gini_handles_empty_input(genome, flat_h5):
    assert _simulator(genome, flat_h5)._calculate_gini(np.array([])) == 0.0


@pytest.mark.parametrize(
    "coverage,uniformity,enrichment,expected",
    [
        (1.0, 1.0, 1000.0, "EXCELLENT"),
        (0.1, 0.1, 1.0, "POOR"),
    ],
)
def test_composite_score_ranks_the_obvious_cases(
    genome, flat_h5, coverage, uniformity, enrichment, expected
):
    sim = _simulator(genome, flat_h5)
    scored = sim._calculate_composite_score(coverage, uniformity, enrichment, "fast")
    assert scored["recommendation"] == expected


def test_composite_score_is_monotonic_in_coverage(genome, flat_h5):
    sim = _simulator(genome, flat_h5)
    low = sim._calculate_composite_score(0.2, 0.8, 10.0, "fast")["composite_score"]
    high = sim._calculate_composite_score(0.9, 0.8, 10.0, "fast")["composite_score"]
    assert high > low


def test_enrichment_below_one_does_not_score_negative(genome, flat_h5):
    """Specificity is log-scaled; enrichment under 1x must floor, not go negative."""
    sim = _simulator(genome, flat_h5)
    scored = sim._calculate_composite_score(0.5, 0.5, 0.01, "fast")
    assert scored["composite_score"] >= 0.0
