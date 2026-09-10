"""`expand-primers` scored every candidate as if the host genome were empty.

The four call sites below build a `PositionCache` over `fg_prefixes` alone and
then hand `bg_prefixes` to a `HybridOptimizer`. `PositionCache.get_positions`
returns an empty array for a prefix it was not built over, silently, so every
background lookup came back empty.

The manifestation is in `NetworkOptimizer._evaluate_primer_addition`, whose
score is `fg_improvement / (1.0 + bg_added)`. With an fg-only cache `bg_added`
is always 0.0, so every candidate scores as perfectly selective and the
selectivity denominator does nothing.

That is worse here than elsewhere: `expand-primers` exists to add primers to an
existing panel, so specificity is exactly the property the user is asking it to
preserve.

The existing expansion tests cannot catch this. Their fake caches ignore the
prefix argument entirely, so an fg-only cache and a two-prefix one behave
identically under them.
"""

import random

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")

from neoswga.core.position_cache import PositionCache
from neoswga.core.utility import reverse_complement

GENOME_LENGTH = 60_000


@pytest.fixture(scope="module")
def planted():
    """A target and a background sharing some primers, written as HDF5.

    No jellyfish and no FASTA: `tests/test_hybrid_optimizer_run.py` establishes
    this pattern, and a test about background counting must not itself depend
    on an external tool to produce the background.
    """
    rng = random.Random(31415)
    primers = []
    while len(primers) < 8:
        p = "".join(rng.choice("ACGT") for _ in range(10))
        if p not in primers:
            primers.append(p)

    fg = list("".join(rng.choice("ACGT") for _ in range(GENOME_LENGTH)))
    bg = list("".join(rng.choice("ACGT") for _ in range(GENOME_LENGTH)))

    for i, primer in enumerate(primers):
        for j in range(3):
            pos = (i * 3_100 + j * 900) % (GENOME_LENGTH - 20)
            fg[pos : pos + 10] = list(primer)

    # The first four also bind the background, densely. Sharing matters: with a
    # background that binds nothing, the count is zero either way and the test
    # would pass without exercising anything.
    shared = primers[:4]
    for i, primer in enumerate(shared):
        for j in range(5):
            pos = (i * 1_700 + j * 2_300) % (GENOME_LENGTH - 20)
            bg[pos : pos + 10] = list(primer)

    return {"primers": primers, "shared": shared, "fg": "".join(fg), "bg": "".join(bg)}


@pytest.fixture
def prefixes(tmp_path, planted):
    def write(prefix, sequence):
        with h5py.File(f"{prefix}_10mer_positions.h5", "w") as fh:
            for primer in planted["primers"]:
                for key in {primer, reverse_complement(primer)}:
                    hits, i = [], sequence.find(key)
                    while i != -1:
                        hits.append(i)
                        i = sequence.find(key, i + 1)
                    if hits:
                        fh.create_dataset(key, data=np.array(hits, dtype=np.int64))

    fg_prefix = str(tmp_path / "fg")
    bg_prefix = str(tmp_path / "bg")
    write(fg_prefix, planted["fg"])
    write(bg_prefix, planted["bg"])
    return fg_prefix, bg_prefix


def test_the_fixture_really_plants_background_hits(prefixes, planted):
    """Guard the guard. If the background genome held none of these primers,
    every assertion below would pass for the wrong reason."""
    _fg, bg = prefixes
    cache = PositionCache([bg], planted["primers"])
    total = sum(len(cache.get_positions(bg, p, "both")) for p in planted["shared"])
    assert total > 0


def test_an_fg_only_cache_now_refuses_the_background_query(prefixes, planted):
    """The defect itself, stated as a property of the cache.

    This is what every one of the five call sites was handing the optimizer. It
    used to return an empty array, silently. It now raises, which is what stops
    a sixth site being written: fixing the call sites removed the bug, and this
    makes the same mistake impossible to make quietly.
    """
    from neoswga.core.position_cache import MissingPositionsError

    fg, bg = prefixes
    fg_only = PositionCache([fg], planted["primers"])

    with pytest.raises(MissingPositionsError):
        fg_only.get_positions(bg, planted["shared"][0], "both")


def test_a_two_prefix_cache_counts_the_background(prefixes, planted):
    fg, bg = prefixes
    cache = PositionCache([fg, bg], planted["primers"])

    assert sum(len(cache.get_positions(bg, p, "both")) for p in planted["shared"]) > 0


def _background_seen_by_the_optimizer(cache, fg, bg, primers, lengths):
    from neoswga.core.hybrid_optimizer import HybridOptimizer

    opt = HybridOptimizer(
        position_cache=cache,
        fg_prefixes=[fg],
        fg_seq_lengths=[lengths],
        bg_prefixes=[bg],
        bg_seq_lengths=[lengths],
    )
    return opt._count_background_sites(primers)


def test_the_optimizer_sees_the_background_when_the_cache_holds_it(prefixes, planted):
    """`_count_background_sites` is what pruning and the selectivity
    denominator both read."""
    fg, bg = prefixes
    cache = PositionCache([fg, bg], planted["primers"])

    assert _background_seen_by_the_optimizer(cache, fg, bg, planted["shared"], GENOME_LENGTH) > 0


def test_the_optimizer_is_stopped_rather_than_told_zero(prefixes, planted):
    """The same call, the same primers, the same background prefix.

    Before both changes this returned 0, indistinguishable downstream from a
    perfectly specific panel. An optimizer handed a cache that cannot answer for
    its background now fails loudly instead of scoring every candidate as
    perfectly selective.
    """
    from neoswga.core.position_cache import MissingPositionsError

    fg, bg = prefixes
    fg_only = PositionCache([fg], planted["primers"])

    with pytest.raises(MissingPositionsError):
        _background_seen_by_the_optimizer(fg_only, fg, bg, planted["shared"], GENOME_LENGTH)


def test_no_site_forwards_background_prefixes_it_did_not_index():
    """The defect as a property of the package, not of one file.

    A function that hands `bg_prefixes` to something which then queries the
    cache by prefix, while building that cache without them, gets zero for
    every background lookup. Five functions did.

    This walks the AST rather than grepping for a string: it asks whether each
    function both forwards `bg_prefixes` and constructs a `PositionCache` whose
    prefix argument does not mention them. Written this way because a grep for
    the fixed spelling would pass on a sixth site written differently, and this
    check is what found the fifth site after a manual review had settled on
    four.

    A function that builds a SEPARATE background cache is not a hit: its
    `PositionCache(bg_prefixes, ...)` call names them. `cli/evaluate.py` and
    `run_contract_set` both do that and are correctly excluded.
    """
    import ast
    import pathlib

    offenders = []
    for path in sorted(pathlib.Path("neoswga").rglob("*.py")):
        src = path.read_text()
        if "PositionCache(" not in src:
            continue
        for fn in ast.walk(ast.parse(src)):
            if not isinstance(fn, (ast.FunctionDef, ast.AsyncFunctionDef)):
                continue
            body = ast.get_source_segment(src, fn) or ""
            if "PositionCache(" not in body:
                continue
            if "bg_prefixes=" not in body and "bg_prefixes," not in body:
                continue
            # Resolve one level of aliasing: `all_prefixes = fg + bg` then
            # `PositionCache(all_prefixes, ...)` is correct, and reading the
            # argument text alone would call it a defect.
            aliases = {}
            for node in ast.walk(fn):
                if isinstance(node, ast.Assign) and len(node.targets) == 1:
                    target = node.targets[0]
                    if isinstance(target, ast.Name):
                        aliases[target.id] = ast.unparse(node.value)

            names = []
            for node in ast.walk(fn):
                if (
                    isinstance(node, ast.Call)
                    and isinstance(node.func, ast.Name)
                    and node.func.id == "PositionCache"
                    and node.args
                ):
                    arg = ast.unparse(node.args[0])
                    names.append(arg + " " + aliases.get(arg, ""))

            if names and not any("bg_prefixes" in n for n in names):
                offenders.append(
                    f"{path}:{fn.lineno} {fn.name}() builds PositionCache({names[0].strip()})"
                )

    assert not offenders, "\n".join(offenders)
