"""Two more coverage paths: one live and slow, one dead and wrong.

`AdvancedFeatureEngineer._estimate_coverage` is a per-primer random-forest
feature -- "what fraction of the genome lies within reach of a site of this
primer". The semantics match the authoritative metric; the implementation did
not. It built a Python `set` of every covered base index, so a primer with 500
sites on a 3.2 Mb genome inserted three million integers one at a time. Measured
at 61 ms per primer against 1.1 ms for the same answer via a numpy bool array --
56x, and it runs once per candidate, so a 2,000-primer pool spent about two
minutes there for nothing.

The reach was also the literal 3000 rather than the polymerase's, which is the
same fault fixed in `background_aware_optimizer` and `primer_expansion`.

`thermodynamics.get_reaction_conditions_from_params` is the dead one. It has no
callers anywhere, hand-lists ten of `ReactionConditions`' twenty fields --
dropping propanediol, glycerol, PEG, BSA, SSB and every buffer species -- and
defaults magnesium to 2.0 mM, a value retired when the model was recalibrated
for isothermal amplification. It is also outside the `PRODUCTION_SITES` guard
in tests/test_reaction_conditions_completeness.py, so nothing would have
stopped it being revived. `build_reaction_conditions` is the sanctioned
replacement and reads the field list off the constructor.
"""

import random

import numpy as np
import pytest


def _union_fraction(positions, genome_length, reach):
    """The authoritative shape, computed independently here."""
    occupied = np.zeros(genome_length, dtype=bool)
    for pos in positions:
        occupied[max(0, pos - reach) : min(genome_length, pos + reach)] = True
    return occupied.sum() / genome_length


@pytest.fixture
def engineer():
    from neoswga.core.advanced_features import AdvancedFeatureEngineer

    instance = AdvancedFeatureEngineer.__new__(AdvancedFeatureEngineer)
    instance.genome_length = 500_000
    return instance


@pytest.mark.parametrize("n_sites", [1, 5, 200])
def test_the_estimate_matches_a_union_of_windows(engineer, n_sites):
    rng = random.Random(n_sites)
    positions = sorted(rng.sample(range(engineer.genome_length), n_sites))

    assert engineer._estimate_coverage(positions, max_gap=3000) == pytest.approx(
        _union_fraction(positions, engineer.genome_length, 3000)
    )


def test_overlapping_sites_are_counted_once(engineer):
    """The property a union has and a sum does not."""
    assert engineer._estimate_coverage([1000, 1001, 1002], max_gap=3000) == pytest.approx(
        _union_fraction([1000, 1001, 1002], engineer.genome_length, 3000)
    )


def test_windows_are_clipped_at_the_genome_ends(engineer):
    assert engineer._estimate_coverage([0], max_gap=3000) == pytest.approx(
        3000 / engineer.genome_length
    )


def test_no_positions_is_no_coverage(engineer):
    assert engineer._estimate_coverage([], max_gap=3000) == 0.0


def test_it_does_not_build_a_set_of_every_covered_base(engineer):
    """The performance defect, pinned by the mechanism rather than by a clock.

    A wall-time assertion would be flaky on a loaded machine; this fails only if
    someone reintroduces per-base Python iteration.
    """
    import ast
    import inspect
    import textwrap

    from neoswga.core.advanced_features import AdvancedFeatureEngineer

    source = textwrap.dedent(inspect.getsource(AdvancedFeatureEngineer._estimate_coverage))
    calls = {
        node.func.id
        for node in ast.walk(ast.parse(source))
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
    }

    assert "range" not in calls, "per-base iteration is back; use a numpy mask"


def test_the_dead_conditions_builder_is_gone():
    """It dropped ten of twenty fields and defaulted magnesium to the retired
    2.0 mM. Nothing called it, and nothing guarded it against being revived."""
    import neoswga.core.thermodynamics as thermodynamics

    assert not hasattr(thermodynamics, "get_reaction_conditions_from_params")
