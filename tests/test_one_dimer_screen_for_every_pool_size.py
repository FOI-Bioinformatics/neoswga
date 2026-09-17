"""One decision about how to screen dimers, in one place.

Phase 4 increment 5 of the plan for `docs/validation/pipeline_audit_2026-09-16/`,
audit finding F11.

`dimer_matrix.build` allocates an n-by-n boolean array over the whole pool.
`dominating_set_optimizer` learned to switch to the lazy screen above a
threshold; two other call sites did not. `refine_hybrid_stage2` built the matrix
unconditionally and `plan-pool` sets `refinement_method="swap"`, so that was on
the hot path: at 491,836 retained candidates the array is roughly 242 GB, a
MemoryError rather than a slowdown.

Three copies of one decision drift. `lazy_dimer.dimer_screen` makes it once,
and the tests below are what stops a fourth copy appearing.
"""

import ast
import pathlib

import pytest

from neoswga.core.lazy_dimer import (
    LAZY_DIMER_POOL_THRESHOLD,
    LazyDimerCompatibility,
    dimer_screen,
)

ROOT = pathlib.Path(__file__).resolve().parent.parent
PACKAGE = ROOT / "neoswga"

# The one module allowed to reach for the dense matrix, because deciding
# between it and the lazy screen is its job.
SCREEN_MODULE = "core/lazy_dimer.py"

# Sites that build a dense matrix for something other than screening a search:
# a report figure, an exhaustive clique enumeration that needs the whole
# relation, or a filter with no production caller. Each needs a reason.
ALLOWED_DENSE = {
    "core/dimer.py": "the module that defines these helpers; the calls are its own "
    "internal composition, not a search screening a pool",
    "core/dimer_validator.py": "works on a delivered panel, not a candidate pool: "
    "`incompatible_pairs` is called from `pool_planner._assess` with the panel, and "
    "`build_matrix` has no caller outside its own docstring",
    "core/thermodynamic_filter.py": "bounded sliding window; no production caller "
    "passes check_heterodimers=True (see the parallelism audit)",
    "core/clique_optimizer.py": "max-clique needs the whole compatibility relation "
    "by construction, and is documented for pools of about 200",
}


def _dense_matrix_callers():
    """Modules that call `dimer_matrix.build` or a heterodimer matrix helper."""
    offenders = {}
    for path in sorted(PACKAGE.rglob("*.py")):
        rel = path.relative_to(PACKAGE).as_posix()
        if rel == SCREEN_MODULE:
            continue
        tree = ast.parse(path.read_text())
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            func = node.func
            name = None
            if isinstance(func, ast.Attribute):
                name = func.attr
            elif isinstance(func, ast.Name):
                name = func.id
            if name in {"build", "heterodimer_matrix", "heterodimer_matrix_fast"}:
                # `build` is a common name; only count it when the module
                # actually imports the dimer matrix.
                source = path.read_text()
                if name != "build" or "dimer_matrix" in source:
                    offenders.setdefault(rel, set()).add(name)
    return offenders


def test_no_search_builds_a_dense_dimer_matrix_of_its_own():
    """The ratchet. A fourth copy of the decision fails here."""
    offenders = {
        rel: names for rel, names in _dense_matrix_callers().items() if rel not in ALLOWED_DENSE
    }

    assert not offenders, (
        "These modules build a dense dimer matrix directly instead of asking "
        "lazy_dimer.dimer_screen, so they will allocate n-by-n over whatever pool "
        "they are handed:\n"
        + "\n".join(f"  {rel}: {sorted(names)}" for rel, names in sorted(offenders.items()))
        + "\n\nRoute through dimer_screen, or add an entry to ALLOWED_DENSE with the reason."
    )


def test_the_allowed_list_has_no_stale_entries():
    """It can only shrink, like the other allowlists in this suite."""
    actual = set(_dense_matrix_callers())
    stale = set(ALLOWED_DENSE) - actual

    assert not stale, (
        f"ALLOWED_DENSE names {sorted(stale)}, which no longer build a dense "
        "matrix. Remove the entry; a fix is not finished until it drops its excuse."
    )


# -- what the screen chooses ----------------------------------------------


def _pool(n):
    """Distinct 12-mers, deterministic, without relying on a generator."""
    letters = "ACGT"
    out = []
    for i in range(n):
        digits = []
        value = i
        for _ in range(12):
            digits.append(letters[value % 4])
            value //= 4
        out.append("".join(digits))
    return out


def test_a_small_pool_gets_the_dense_matrix():
    """Below the threshold the array is small and answers without recomputation."""
    screen = dimer_screen(_pool(50), max_dimer_bp=3)

    assert not isinstance(screen, LazyDimerCompatibility)


def test_a_pool_past_the_threshold_gets_the_lazy_screen():
    """The case that was a MemoryError on three of four call sites."""
    screen = dimer_screen(_pool(LAZY_DIMER_POOL_THRESHOLD + 1), max_dimer_bp=3)

    assert isinstance(screen, LazyDimerCompatibility)


def test_the_two_screens_agree_on_the_same_pool():
    """Switching on pool size must not switch the answer."""
    pool = _pool(80)
    dense = dimer_screen(pool, max_dimer_bp=3)
    lazy = LazyDimerCompatibility(3)

    panel = pool[:6]
    for candidate in pool[6:40]:
        assert dense.dimerises(candidate, panel) == lazy.dimerises(candidate, panel), candidate


def test_the_lazy_screen_is_not_quadratic_in_the_pool():
    """It is constructed from a threshold, not from a pool."""
    screen = dimer_screen(_pool(LAZY_DIMER_POOL_THRESHOLD + 1), max_dimer_bp=3)
    panel = _pool(8)

    screen.dimerises("ACGTACGTACGT", panel)

    assert screen.computations <= len(panel)


def test_a_threshold_the_dense_matrix_cannot_hold_still_gets_screened():
    """Not refused, and above all not silently weakened.

    The dense matrix codes t-mers in a 4**8 space and cannot hold 8 or above.
    `dimer.is_dimer_fast`, which the lazy screen uses, has no such limit, so the
    pair can always enforce what was configured. A screen that cannot be
    enforced must not read as "no conflicts": an 11 bp delivered heterodimer
    against a configured 3 is what that looks like.
    """
    screen = dimer_screen(_pool(50), max_dimer_bp=8)

    assert isinstance(screen, LazyDimerCompatibility)
    assert screen.max_dimer_bp == 8


def test_the_loose_threshold_is_actually_applied():
    """Choosing the lazy branch is no use if it then screens on something else."""
    # Eight complementary bases and no more, so a threshold of 8 flags it and a
    # threshold of 3 flags it too, while a pair sharing nothing flags at neither.
    screen = dimer_screen(["AAAAAAAAAAAA", "TTTTTTTTTTTT", "ACACACACACAC"], max_dimer_bp=8)

    assert screen.dimerises("AAAAAAAAAAAA", ["TTTTTTTTTTTT"])


def test_the_screen_is_never_silently_weaker_than_asked():
    """Whatever branch it took, the threshold it reports is the one requested."""
    for pool_size in (50, LAZY_DIMER_POOL_THRESHOLD + 1):
        for threshold in (3, 4, 7):
            screen = dimer_screen(_pool(pool_size), max_dimer_bp=threshold)
            assert screen.max_dimer_bp == threshold, (pool_size, threshold)


# -- the unread mirror ----------------------------------------------------


def test_the_bipartite_graph_keeps_no_networkx_mirror():
    """It was written on every edge and read by nothing.

    `BipartiteGraph` maintained a `networkx.Graph` alongside its own
    `primer_to_regions` and `region_to_primers` dicts, adding an edge per
    primer-bin pair. Nothing in the package or the suite read it. The
    `full_network.graph` reads in `hybrid_optimizer` are an
    `AmplificationNetwork`, a different class whose graph is load-bearing.
    """
    from neoswga.core.dominating_set_optimizer import BipartiteGraph

    graph = BipartiteGraph(bin_size=1000)

    assert not hasattr(graph, "graph"), (
        "BipartiteGraph carries a networkx mirror again. It costs an edge per "
        "primer-bin pair and nothing reads it."
    )


def test_coverage_bookkeeping_survives_without_the_mirror():
    """The structures that are actually read still answer."""
    import numpy as np

    from neoswga.core.dominating_set_optimizer import BipartiteGraph

    graph = BipartiteGraph(bin_size=1000)
    graph.add_primer_coverage(
        "ACGTACGTACGT",
        np.array([500, 2500], dtype=np.int64),
        "target",
        10_000,
        extension_reach=100,
    )

    assert "ACGTACGTACGT" in graph.primers
    assert graph.primer_to_regions["ACGTACGTACGT"], "no regions recorded"
    for region in graph.primer_to_regions["ACGTACGTACGT"]:
        assert "ACGTACGTACGT" in graph.region_to_primers[region]
