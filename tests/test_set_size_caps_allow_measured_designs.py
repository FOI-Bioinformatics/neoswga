"""The set-size ceilings sat below the sizes coverage actually needs.

`num_primers` and `target_set_size` were schema-capped at 50, and
`estimate_optimal_set_size` clamped its own answer to 20. Neither number is
sourced; both predate the measurements below.

Calibrating the per-primer reach against Clarke et al. (2017) *Wolbachia* --
the one published set with a quantitative breadth outcome -- puts it between
3.0 and 6.2 kb (see docs/validation/reach_calibration.md). Re-running the
Prevotella sweep inside that band:

    reach    >=90%   >=95%   >=99%
    3.0 kb      76      96     126
    4.5 kb      40      51      70
    6.2 kb      24      31      43

So a 95% design needs 31-96 primers depending where in the band the true reach
sits, and the cap of 50 forbids expressing most of that range. A user asking
for 96 primers got a schema validation error; one relying on `--auto-size` got
at most 20 whatever the genome.

The caps are kept -- they are a guard against a runaway config, and synthesis
cost is a real constraint the `oligo_cost` model exists to price -- but raised
past the sizes the measurement calls for.
"""

import json
from pathlib import Path

import pytest

SCHEMA = json.loads(
    (
        Path(__file__).parent.parent / "neoswga" / "core" / "schema" / "params.schema.json"
    ).read_text()
)

# The largest set size the reach calibration implies, at the low end of the
# calibrated band (3.0 kb) for 99% coverage of a 3.2 Mb target.
MEASURED_MAX = 126


@pytest.mark.parametrize("field", ["num_primers", "target_set_size"])
def test_the_schema_admits_a_measured_design(field):
    assert SCHEMA["properties"][field]["maximum"] >= MEASURED_MAX


@pytest.mark.parametrize("field", ["num_primers", "target_set_size"])
def test_the_schema_still_bounds_a_runaway_config(field):
    """A ceiling that is gone is not a ceiling. Every primer is an oligo
    someone pays for."""
    assert SCHEMA["properties"][field]["maximum"] <= 500


@pytest.mark.parametrize("size", [64, 96, 126])
def test_a_large_set_size_validates(size):
    from neoswga.core.param_validator import ParamValidator, ValidationLevel

    messages = ParamValidator().validate_params(
        {"num_primers": size, "fg_genomes": ["x.fna"], "data_dir": "."}
    )
    errors = [
        m for m in messages if m.level == ValidationLevel.ERROR and "num_primers" in m.message
    ]

    assert not errors, [m.message for m in errors]


def test_auto_size_is_genome_independent_and_pinned_as_such():
    """A known defect, recorded rather than hidden behind a raised ceiling.

    `genome_length` cancels out of `coverage_per_primer` -- it appears in
    `expected_sites` and again in the denominator -- so the estimate depends
    only on primer length. It is the same for a 6 kb plasmid and a 50 Mb genome,
    and at every SWGA primer length the clamp is the whole answer.

    The cause is the random-sequence site model: `expected_sites` assumes a
    k-mer occurs 4^-k times per base, but SWGA candidates are selected for being
    abundant. The Prevotella 12-mer pool has a median of 8 sites per primer
    against the 0.38 this predicts.

    The ceiling stays at 20 until the model is recalibrated against
    docs/validation/reach_calibration.md, because `--auto-size` output is spent
    on oligos and a wrong number that is small does less damage than a wrong
    number that is large. `num_primers` in params.json is the way to ask for
    more, and the schema now admits it.
    """
    from neoswga.core.set_size_optimizer import (
        create_baseline_effects,
        estimate_optimal_set_size,
    )

    effects = create_baseline_effects()
    sizes = {
        length: estimate_optimal_set_size(length, 12, 0.95, 3000, effects)
        for length in (6_000, 100_000, 3_168_282, 50_000_000)
    }

    assert len(set(sizes.values())) == 1, f"the defect has changed shape: {sizes}"
    assert set(sizes.values()) == {20}, sizes
