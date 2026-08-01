"""Behavioural tests for `SimulationBasedEvaluator`.

This is documented as providing "ground truth for validation" and is meant for
comparing optimization methods. It sat at 21% coverage, and two independent
faults meant it returned zero for every primer set -- so it would have ranked
every method identically while looking like it had measured something.

Both faults are the same shape as others found in this audit: a lookup that
misses returns empty rather than raising, and nothing downstream can tell an
empty result from a real one. The tests below therefore assert on non-zero
outcomes for a set that demonstrably binds, not merely that the code runs.
"""

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")

from neoswga.core.position_cache import PositionCache
from neoswga.core.simulation_fitness import SimulationBasedEvaluator, SimulationFitness
from neoswga.core.thermodynamics import reverse_complement

PRIMERS = ["TTGACCATGA", "GCATTACGGT", "AACCGGTTAC"]


def _occurrences(sequence, pattern):
    out, i = [], sequence.find(pattern)
    while i != -1:
        out.append(i)
        i = sequence.find(pattern, i + 1)
    return out


@pytest.fixture(scope="module")
def genome_sequence():
    import random

    rng = random.Random(31337)
    seq = list("".join(rng.choice("ACGT") for _ in range(20_000)))
    for offset, primer in enumerate(PRIMERS):
        for pos in (1_000, 7_000, 13_000):
            start = pos + offset * 100
            seq[start : start + len(primer)] = list(primer)
    return "".join(seq)


@pytest.fixture
def cache(tmp_path, genome_sequence):
    """A PositionCache over a real HDF5, in the layout the pipeline writes."""
    prefix = str(tmp_path / "target")
    with h5py.File(prefix + "_10mer_positions.h5", "w") as f:
        for primer in PRIMERS:
            for key in {primer, reverse_complement(primer)}:
                pos = _occurrences(genome_sequence, key)
                if pos:
                    f.create_dataset(key, data=np.array(pos, dtype=np.int32))
    return PositionCache([prefix], PRIMERS)


@pytest.fixture
def evaluator(genome_sequence, cache):
    return SimulationBasedEvaluator(
        genome_sequence=genome_sequence,
        genome_length=len(genome_sequence),
        position_cache=cache,
        n_replicates=2,
        simulation_duration=3600.0,
    )


# ----------------------------------------------------------------------
# Fault 1: the position lookup used the wrong key
# ----------------------------------------------------------------------


def test_positions_are_read_from_the_cache(evaluator):
    """It queried `get_positions("", primer, ...)`.

    PositionCache is keyed by the real path prefix, so the empty string matched
    nothing. `get_positions` returns an empty array rather than raising, so
    every primer silently came back with no binding sites.
    """
    positions = evaluator._build_position_dict(PRIMERS)

    total = sum(len(v["forward"]) + len(v["reverse"]) for v in positions.values())
    assert total > 0, "no positions read from a cache that demonstrably has them"

    for primer in PRIMERS:
        assert len(positions[primer]["forward"]) == 3, primer


def test_missing_primer_yields_empty_lists_not_an_error(evaluator):
    positions = evaluator._build_position_dict(["CGCGCGCGCG"])
    assert positions["CGCGCGCGCG"] == {"forward": [], "reverse": []}


def test_a_cache_with_no_prefixes_is_warned_about(genome_sequence, caplog):
    """Silence here is what made the original fault invisible."""
    import logging

    class _Empty:
        fname_prefixes = []

        def get_positions(self, *_args, **_kwargs):
            return np.array([])

    evaluator = SimulationBasedEvaluator(
        genome_sequence=genome_sequence,
        genome_length=len(genome_sequence),
        position_cache=_Empty(),
        n_replicates=1,
        simulation_duration=10.0,
    )
    with caplog.at_level(logging.WARNING):
        evaluator._build_position_dict(PRIMERS)

    assert "no genome prefixes" in caplog.text or "No binding sites" in caplog.text


# ----------------------------------------------------------------------
# Fault 2: the simulated enzyme and buffer disagreed
# ----------------------------------------------------------------------


def test_default_conditions_and_default_polymerase_agree(evaluator):
    """The config named phi29 while the defaults were equiphi29 at 42 C.

    phi29-length primers have Tm around 26-31 C and cannot prime at 42 C, so
    every set produced zero forks -- for a reason unrelated to the primers.
    """
    conditions = evaluator._get_default_conditions()
    polymerase = evaluator._get_polymerase_type(conditions)

    assert polymerase == conditions.polymerase, (
        f"simulating {polymerase} in a {conditions.polymerase} buffer at " f"{conditions.temp} C"
    )


def test_default_conditions_match_the_project_default_polymerase(evaluator):
    """phi29 at 30 C, the same default the pipeline designs for."""
    conditions = evaluator._get_default_conditions()
    assert conditions.polymerase == "phi29"
    assert conditions.temp == pytest.approx(30.0)


def test_polymerase_is_taken_from_supplied_conditions(evaluator):
    from neoswga.core.reaction_conditions import ReactionConditions

    hot = ReactionConditions(temp=42.0, polymerase="equiphi29")
    assert evaluator._get_polymerase_type(hot) == "equiphi29"


# ----------------------------------------------------------------------
# The evaluation as a whole
# ----------------------------------------------------------------------


def test_a_binding_set_gets_non_zero_fitness(evaluator):
    """The headline: both faults independently produced zero here."""
    fitness = evaluator.evaluate(PRIMERS, verbose=False)

    assert isinstance(fitness, SimulationFitness)
    assert fitness.mean_coverage > 0, "a set with 9 binding sites covered nothing"
    assert fitness.mean_forks_created > 0, "no replication forks were created"
    assert fitness.fitness_score > 0


def test_fitness_metrics_are_bounded(evaluator):
    fitness = evaluator.evaluate(PRIMERS, verbose=False)

    assert 0.0 <= fitness.mean_coverage <= 1.0
    assert 0.0 <= fitness.coverage_uniformity <= 1.0
    assert 0.0 <= fitness.fitness_score <= 1.0
    assert fitness.n_replicates == 2


def test_a_set_that_binds_nowhere_scores_zero_and_says_so(evaluator, caplog):
    """Zero is the right answer here; it must be distinguishable from the bug."""
    import logging

    with caplog.at_level(logging.WARNING):
        fitness = evaluator.evaluate(["CGCGCGCGCG"], verbose=False)

    assert fitness.fitness_score == 0.0
    assert fitness.n_replicates == 0
    assert "No primer positions" in caplog.text or "No binding sites" in caplog.text


def test_comparing_sets_can_distinguish_them(evaluator):
    """The stated use case. With everything scoring zero it could not."""
    results = evaluator.compare_sets(
        [("binds", PRIMERS), ("does-not-bind", ["CGCGCGCGCG"])], verbose=False
    )

    ranked = [name for name, _ in results]
    assert ranked[0] == "binds", "the binding set did not rank first"
    assert results[0][1].fitness_score > results[1][1].fitness_score
