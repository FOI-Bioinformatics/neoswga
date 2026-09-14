"""A delivered panel must either honour `max_dimer_bp` or say that it does not.

Both halves of this were tested apart and never together: the validator against
a hand-built primer list, and the optimizers without the validator. So nothing
related what an optimizer actually delivered to the limit it was configured
with, and the shipped panels carried an 11 bp heterodimer against a configured
3 for as long as they did partly because no test could see it.

This is deliberately written to hold under either answer to the open question
of whether `--num-primers` is a guarantee or a request. What it forbids is the
third thing: quietly returning the requested count with a violating pair in it.

- Strict (the current default): every delivered pair is within the limit. The
  panel may be shorter than requested, and the run says so.
- Relaxed (`allow_dimer_relaxation=True`): the panel may carry a violating
  pair, but the run must warn and name the primer it admitted.
"""

import logging
import random

import pytest

h5py = pytest.importorskip("h5py")
pytest.importorskip("networkx")

import numpy as np

from neoswga.core.dimer_validator import DimerValidator
from neoswga.core.dominating_set_optimizer import DominatingSetOptimizer
from neoswga.core.position_cache import PositionCache
from neoswga.core.thermodynamics import reverse_complement

GENOME_LENGTH = 100_000
PRIMER_LENGTH = 10
COVERAGE_REACH = 3_000
N_CANDIDATES = 40
MAX_DIMER_BP = 3
REQUESTED = 20


@pytest.fixture(scope="module")
def pool(tmp_path_factory):
    """A target big enough that coverage does not saturate before the request.

    The bundled 6 kb plasmid cannot exercise this: a 3 kb reach covers it with
    two primers, so the requested count never binds and the greedy never
    stalls on a dimer.
    """
    rng = random.Random(4242)
    seq = list("".join(rng.choice("ACGT") for _ in range(GENOME_LENGTH)))

    primers = []
    for i in range(N_CANDIDATES):
        primer = "".join(rng.choice("ACGT") for _ in range(PRIMER_LENGTH))
        if primer in primers:
            continue
        primers.append(primer)
        for j in range(3):
            pos = (i * (GENOME_LENGTH // N_CANDIDATES) + j * 900) % (GENOME_LENGTH - 20)
            seq[pos : pos + PRIMER_LENGTH] = list(primer)

    prefix = str(tmp_path_factory.mktemp("dimer") / "target")
    keys = set()
    for primer in primers:
        keys.update({primer, reverse_complement(primer)})
    text = "".join(seq)
    with h5py.File(f"{prefix}_{PRIMER_LENGTH}mer_positions.h5", "w") as handle:
        for key in sorted(keys):
            positions, i = [], text.find(key)
            while i != -1:
                positions.append(i)
                i = text.find(key, i + 1)
            if positions:
                handle.create_dataset(key, data=np.array(positions, dtype=np.int64))
    return {"primers": primers, "prefix": prefix}


def _optimize(pool, *, allow_relaxation):
    optimizer = DominatingSetOptimizer(
        PositionCache([pool["prefix"]], pool["primers"]),
        fg_prefixes=[pool["prefix"]],
        fg_seq_lengths=[GENOME_LENGTH],
        extension_reach=COVERAGE_REACH,
        max_dimer_bp=MAX_DIMER_BP,
        allow_dimer_relaxation=allow_relaxation,
    )
    # `optimize_greedy` returns a result dict; the panel is under "primers".
    return optimizer.optimize_greedy(pool["primers"], max_primers=REQUESTED, verbose=False)[
        "primers"
    ]


def _violations(primers):
    return DimerValidator(MAX_DIMER_BP, MAX_DIMER_BP + 1).incompatible_pairs(list(primers))


def test_the_fixture_actually_stalls_on_the_dimer_limit(pool):
    """Guard the guard: a pool that never stalls would prove nothing below."""
    strict = _optimize(pool, allow_relaxation=False)
    assert len(strict) < REQUESTED, (
        "the candidate pool supports the whole request at "
        f"max_dimer_bp={MAX_DIMER_BP}, so neither contract below is exercised"
    )


def test_a_strict_panel_carries_no_pair_above_the_limit(pool, caplog):
    with caplog.at_level(logging.WARNING):
        primers = _optimize(pool, allow_relaxation=False)

    assert _violations(primers) == [], (
        "the delivered panel carries a pair above the configured " f"max_dimer_bp={MAX_DIMER_BP}"
    )
    # Short of what was asked for, and not silently so.
    assert len(primers) < REQUESTED
    assert any(
        "max_dimer_bp" in record.message for record in caplog.records
    ), "a panel short of the request must say why"


def test_a_relaxed_panel_may_violate_the_limit_but_must_report_it(pool, caplog):
    with caplog.at_level(logging.WARNING):
        primers = _optimize(pool, allow_relaxation=True)

    assert _violations(primers), (
        "the relaxed branch is not exercised on this fixture, so the contract " "below is vacuous"
    )
    assert any(
        "unscreened" in record.message for record in caplog.records
    ), "a primer admitted under the relaxed limit must be reported, not delivered quietly"


def test_relaxation_is_what_separates_the_two_panels(pool):
    """The flag must actually change the outcome on this fixture.

    If both sides returned the same panel the two contracts above would be
    vacuous, and the open `--num-primers` question would not be a real choice.
    """
    strict = _optimize(pool, allow_relaxation=False)
    relaxed = _optimize(pool, allow_relaxation=True)
    assert len(relaxed) > len(strict)
