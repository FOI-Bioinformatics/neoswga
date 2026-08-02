"""Design options must change the design.

A primer designer exposes a lot of knobs. The failure mode this file exists to
catch is not a knob that computes the wrong answer -- it is a knob that is
accepted, logged, stored, and then never read. Everything still runs and still
emits a primer set, so nothing distinguishes "this option is working" from
"this option is decoration". The only assertion that separates them is that
changing the option changes the output.

Two distinct faults are pinned here.

The first is the coverage-granularity fault, which makes `--num-primers`
inoperative on the default path. Stage-1 set cover selects over `bin_size` bins
(default 10 kb) while extension is modelled at the realistic per-primer reach
(~3 kb for phi29). `BipartiteGraph.add_primer_coverage` marks a bin covered when
ANY part of it is reached, so at the default settings a primer that reaches 30%
of a bin claims the whole bin. Two primers then appear to cover the genome,
greedy selection stops, and every request from 2 to 12 primers returns the same
two. The reported coverage is 1.000 -- the same set measured at a granularity
finer than the reach covers 0.433. That is the pattern this audit has found
repeatedly: the sentinel for "no information" is indistinguishable from a real
result, and here it is a perfect score.

The second is the parsed-then-ignored flag. `use_cooperative_binding` and
`primer_strategy` are accepted by the optimize parser and written to the global
parameter object, and no module under `neoswga/core/` reads either one.
"""

import random
import subprocess
import sys
from pathlib import Path

import pytest

h5py = pytest.importorskip("h5py")

import numpy as np

from neoswga.core.hybrid_optimizer import HybridOptimizer
from neoswga.core.position_cache import PositionCache
from neoswga.core.thermodynamics import reverse_complement

GENOME_LENGTH = 30_000
COVERAGE_REACH = 3_000
PRIMER_LENGTH = 10

# The production default (hybrid_optimizer.py:1047). It is deliberately used
# here rather than a test-friendly smaller value: a bug that only appears at the
# shipped default is exactly the bug worth pinning.
DEFAULT_BIN_SIZE = 10_000


@pytest.fixture(scope="module")
def genome():
    """A genome with 12 planted primers spread across it.

    Spread matters. A set-cover objective needs candidates whose coverage
    genuinely differs, or every selection is equally good and a test that
    asserts "more primers, more coverage" proves nothing.
    """
    rng = random.Random(9182)
    seq = list("".join(rng.choice("ACGT") for _ in range(GENOME_LENGTH)))

    primers = []
    for i in range(12):
        primer = "".join(rng.choice("ACGT") for _ in range(PRIMER_LENGTH))
        if primer in primers:
            continue
        primers.append(primer)
        for j in range(2 + (i % 3)):
            pos = (i * 2_400 + j * 700) % (GENOME_LENGTH - 20)
            seq[pos : pos + PRIMER_LENGTH] = list(primer)

    return {"seq": "".join(seq), "primers": primers}


@pytest.fixture(scope="module")
def cache(tmp_path_factory, genome):
    prefix = str(tmp_path_factory.mktemp("positions") / "target")
    with h5py.File(f"{prefix}_{PRIMER_LENGTH}mer_positions.h5", "w") as f:
        for primer in genome["primers"]:
            for key in {primer, reverse_complement(primer)}:
                positions, i = [], genome["seq"].find(key)
                while i != -1:
                    positions.append(i)
                    i = genome["seq"].find(key, i + 1)
                if positions:
                    f.create_dataset(key, data=np.array(positions, dtype=np.int32))
    return PositionCache([prefix], genome["primers"]), prefix


def build_optimizer(cache, bin_size):
    position_cache, prefix = cache
    return HybridOptimizer(
        position_cache=position_cache,
        fg_prefixes=[prefix],
        fg_seq_lengths=[GENOME_LENGTH],
        bin_size=bin_size,
        coverage_reach=COVERAGE_REACH,
    )


# ----------------------------------------------------------------------
# --num-primers must select the number of primers it was asked for
# ----------------------------------------------------------------------


@pytest.mark.parametrize("requested", [2, 4, 6, 8])
def test_num_primers_is_honoured_at_the_default_bin_size(cache, genome, requested):
    """The headline defect.

    All 12 candidates together cover the genome, so there is real coverage left
    to buy at every one of these sizes -- nothing here asks the optimizer to
    select primers that add nothing. At `bin_size=1000` the optimizer returns
    exactly what it is asked for; at the shipped default of 10 kb it returns two
    primers for every request.
    """
    optimizer = build_optimizer(cache, DEFAULT_BIN_SIZE)
    result = optimizer.optimize(genome["primers"], final_count=requested, verbose=False)

    assert len(result.primers) == requested


def test_asking_for_more_primers_yields_more_coverage(cache, genome):
    """`--minimize-primers` is a separate opt-in defaulting to False.

    So the default path is not supposed to be minimising, and a larger request
    that buys real coverage must actually buy it.
    """
    optimizer = build_optimizer(cache, DEFAULT_BIN_SIZE)

    small = optimizer.optimize(genome["primers"], final_count=2, verbose=False)
    large = optimizer.optimize(genome["primers"], final_count=8, verbose=False)

    assert len(large.primers) > len(small.primers)
    assert large.final_coverage > small.final_coverage


# ----------------------------------------------------------------------
# Coverage must mean the same thing regardless of bin size
# ----------------------------------------------------------------------


def test_reported_coverage_does_not_depend_on_the_bin_size(cache, genome):
    """Bin size is an implementation detail of the search, not a property of
    the genome. If the number reported to the user moves when it changes, the
    number is measuring the bins rather than the genome.

    The same two primers report 1.000 at a 10 kb bin and 0.433 at a 1 kb bin.
    The 1 kb figure is the honest one: with a 3 kb reach a primer cannot cover
    10 kb, and a coarse bin simply rounds its own error up to a perfect score.
    """
    subset = genome["primers"][:2]

    coarse = build_optimizer(cache, DEFAULT_BIN_SIZE)._calculate_coverage(subset)
    fine = build_optimizer(cache, 1_000)._calculate_coverage(subset)

    assert coarse == pytest.approx(fine, abs=0.05)


# `coverage_bin_size` bins at a quarter of the reach. A site marks every bin it
# touches, so it claims up to `reach + bin_size` bases -- a bounded 25%
# overstatement. Pinning the bound here means a change to that divisor has to be
# made deliberately rather than drifting.
BIN_SLOP = 1.25


def test_coverage_of_a_single_primer_stays_within_its_reach(cache, genome):
    """A physical bound, and the cleanest statement of the original fault.

    One primer with `n` binding sites reaches at most `n * COVERAGE_REACH`
    bases. The binned estimate may exceed that by the documented slop and no
    more. Before the fix a single primer could report 1.000 against a ceiling of
    0.20 -- five times what the chemistry allows, which is how two primers came
    to look like full genome coverage.
    """
    position_cache, prefix = cache
    optimizer = build_optimizer(cache, DEFAULT_BIN_SIZE)

    for primer in genome["primers"]:
        n_sites = len(position_cache.get_positions(prefix, primer, "both"))
        if n_sites == 0:
            continue
        ceiling = min(1.0, n_sites * COVERAGE_REACH * BIN_SLOP / GENOME_LENGTH)

        assert optimizer._calculate_coverage([primer]) <= ceiling + 1e-9, (
            f"{primer} with {n_sites} sites reports more coverage than "
            f"{COVERAGE_REACH} bp of reach per site can produce"
        )


# ----------------------------------------------------------------------
# Options that are parsed and then ignored
# ----------------------------------------------------------------------

CORE = Path(__file__).resolve().parent.parent / "neoswga" / "core"

# Flags the optimize parser accepts that are meant to steer the design.
# Each must be read by something under neoswga/core/, or it is decoration.
DESIGN_FLAGS = [
    "uniformity_weight",
    "minimize_primers",
    "target_coverage",
    "min_per_target_coverage",
]


@pytest.mark.parametrize("flag", DESIGN_FLAGS)
def test_design_flag_is_consumed_by_a_core_module(flag):
    """A flag stored on the parameter object and read by nobody changes nothing.

    These four were each traced end-to-end through `optimize_step4` into the
    optimizers, so they earn their place on the command line.
    """
    hits = subprocess.run(
        ["grep", "-rl", flag, str(CORE), "--include=*.py"],
        capture_output=True,
        text=True,
    ).stdout.split()

    assert hits, f"no module under neoswga/core/ reads {flag!r}"


@pytest.mark.parametrize("flag", ["--use-cooperative-binding", "--primer-strategy"])
def test_removed_flags_are_not_offered(flag):
    """Both were accepted, logged, stored on the parameter object -- and read by
    nothing under `neoswga/core/`, so neither could ever change a design.

    `--use-cooperative-binding` said so in its own help text ("not yet fully
    integrated"). They were removed rather than wired, because a flag that
    silently does nothing is worse than an absent one: it reads as a control
    that was tried and found not to matter. This test keeps them from
    reappearing without a consumer.
    """
    from neoswga.cli_unified import create_parser

    help_text = create_parser().format_help()
    assert flag not in help_text
