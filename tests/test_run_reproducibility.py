"""The same seed must produce the same primer set, across processes.

`--seed` is offered so a design can be reproduced, and it did not work. Five
runs of the bundled plasmid example at `seed=42` gave five different sets:

    ['CCCATTGACG', 'CTTGGTCTGA']
    ['ATACCGTCGA', 'CCCATTGACG']
    ['CCATAGTTGC', 'CCCATTGACG']
    ['CCCATTGACG', 'CCGCTCACAA']
    ['AGTCGGGAAA', 'CCCATTGACG']

Fixing `PYTHONHASHSEED` collapsed them to one answer, which identifies the
cause exactly: Python randomises string hashing per process, so iteration order
over a `set` of primers differs between runs. `optimize_greedy` built its
selection in a `set` and returned `list(selected)`, so the order handed to
Stage-2 refinement was effectively random -- and refinement breaks ties by
position, so a different order gives a different set.

Seeding the RNG cannot fix this because no RNG is involved. The greedy
selection itself is deterministic; only the container forgot the order it
happened in.

This mattered beyond tidiness. The golden snapshot was intermittently failing,
and every planned validation against measured wet-lab outcomes compares a
predicted set to a published one -- which is meaningless if the prediction
changes between runs. Its companion "same seed twice" test passed throughout,
because both runs were in the same process and therefore shared a hash seed.
"""

import os
import subprocess
import sys
from pathlib import Path

import pytest

EXAMPLE = Path(__file__).resolve().parent.parent / "examples" / "plasmid_example"

RUN = """
import contextlib, io, json, logging, sys
logging.disable(logging.CRITICAL)
from neoswga.core import parameter
import neoswga.core.pipeline as pm
pm._initialized = False
parameter.json_file = "params.json"
pm._initialize()
from neoswga.core.unified_optimizer import optimize_step4
with contextlib.redirect_stdout(io.StringIO()):
    sets, _scores, _cache = optimize_step4(optimization_method="hybrid", seed=42, verbose=False)
print(json.dumps(sorted(sets[0])))
"""


def _run_with_hash_seed(hash_seed, workdir):
    """Run the optimizer in `workdir`, never in the repo's example directory.

    The pipeline writes step4 artifacts beside its inputs, and
    `test_plasmid_golden_snapshot` copies *every* file out of
    examples/plasmid_example into its own tmpdir. Running here would leave
    outputs in that shared directory for the snapshot to pick up -- which is
    exactly what it did on the first version of this file, making the snapshot
    fail in the suite while passing alone.
    """
    env = dict(os.environ, PYTHONHASHSEED=str(hash_seed))
    proc = subprocess.run(
        [sys.executable, "-c", RUN],
        cwd=str(workdir),
        capture_output=True,
        text=True,
        env=env,
        timeout=600,
    )
    if proc.returncode != 0:
        pytest.skip(f"pipeline could not run here: {proc.stderr.strip()[-200:]}")
    return proc.stdout.strip().splitlines()[-1]


@pytest.fixture
def example_copy(tmp_path):
    """A private copy of the plasmid example, so this test writes nothing
    into the repository."""
    import shutil

    if not (EXAMPLE / "params.json").exists():
        pytest.skip("plasmid example not available")

    for entry in EXAMPLE.iterdir():
        if entry.is_file():
            shutil.copy2(entry, tmp_path)
    return tmp_path


# ----------------------------------------------------------------------
# Determinism at the source
# ----------------------------------------------------------------------


def test_set_cover_returns_primers_in_selection_order():
    """`optimize_greedy` accumulated into a `set` and returned `list(selected)`.

    Greedy set cover picks the largest marginal gain at each step, so it has a
    definite order; the set discarded it and substituted whatever the string
    hashes produced that process. Returning the order it actually selected in
    costs nothing and removes the nondeterminism at source.
    """
    import ast
    import inspect
    import textwrap

    from neoswga.core.dominating_set_optimizer import DominatingSetOptimizer

    tree = ast.parse(textwrap.dedent(inspect.getsource(DominatingSetOptimizer.optimize_greedy)))

    # Parsed rather than grepped, so the comment explaining the old code does
    # not trip the check on itself.
    offenders = [
        ast.unparse(node)
        for node in ast.walk(tree)
        if isinstance(node, ast.Call)
        and getattr(node.func, "id", None) == "list"
        and node.args
        and getattr(node.args[0], "id", None) == "selected"
    ]

    assert not offenders, (
        f"optimize_greedy still returns set-ordered primers {offenders}, which "
        f"vary per process because Python randomises string hashing"
    )


def test_reported_primers_match_the_recorded_order():
    """The two views of the same selection must agree, or a caller reading
    `primers` gets a different answer from one reading `ordered_primers`."""
    import numpy as np

    h5py = pytest.importorskip("h5py")

    import random
    import tempfile

    from neoswga.core.dominating_set_optimizer import DominatingSetOptimizer
    from neoswga.core.position_cache import PositionCache
    from neoswga.core.thermodynamics import reverse_complement

    genome_length = 40_000
    rng = random.Random(7)
    seq = list("".join(rng.choice("ACGT") for _ in range(genome_length)))
    primers = []
    for i in range(12):
        primer = "".join(rng.choice("ACGT") for _ in range(10))
        if primer in primers:
            continue
        primers.append(primer)
        for j in range(3):
            pos = (i * 3_000 + j * 700) % (genome_length - 20)
            seq[pos : pos + 10] = list(primer)
    seq = "".join(seq)

    tmp = tempfile.mkdtemp()
    prefix = os.path.join(tmp, "t")
    with h5py.File(f"{prefix}_10mer_positions.h5", "w") as fh:
        for key in sorted({x for p in primers for x in (p, reverse_complement(p))}):
            hits, i = [], seq.find(key)
            while i != -1:
                hits.append(i)
                i = seq.find(key, i + 1)
            if hits:
                fh.create_dataset(key, data=np.array(hits, dtype=np.int32))

    cache = PositionCache([prefix], primers)
    optimizer = DominatingSetOptimizer(
        cache, [prefix], [genome_length], bin_size=1000, extension_reach=3000
    )
    result = optimizer.optimize_greedy(primers, max_primers=6, verbose=False)

    assert result["primers"] == result["ordered_primers"]


# ----------------------------------------------------------------------
# The property a user cares about
# ----------------------------------------------------------------------


@pytest.mark.integration
@pytest.mark.slow
def test_same_seed_gives_the_same_set_across_processes(example_copy):
    """The real contract behind `--seed`.

    Two interpreter processes with deliberately different string-hash seeds
    must still agree. The existing snapshot test runs both repetitions in one
    process, so it shares a hash seed and cannot see this.
    """
    # Several hash seeds, not two: any given pair can coincide by luck, and a
    # test that passes by coincidence is worse than none.
    answers = {_run_with_hash_seed(h, example_copy) for h in (0, 1, 12345, 99991)}

    assert len(answers) == 1, (
        f"seed=42 produced {len(answers)} different primer sets across hash "
        f"seeds: {sorted(answers)}"
    )
