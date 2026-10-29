"""The run log must not report a weight the selecting engine does not hold.

`--application` resolved three selection weights and logged "selection weights
applied: tm=0.30, uniformity=0.05, dimer_penalty=0.45" whatever the method,
while on `hybrid` -- the DEFAULT -- the object doing the selecting used 0.0 for
two of the three, and on `dominating-set` and `clique` for all three.

A log line that asserts the opposite of what the code does is worse than
silence. It is part of why the inert-weights defect survived several audits:
anyone checking whether `--application` worked read the line and stopped.

These tests pin the message against the engines rather than against a copy of
the expected string, so the two cannot drift: `test_the_inert_map_matches_the_engines`
derives the answer from each engine's own `__init__` and fails if the map says
a weight is applied where the engine cannot accept it.
"""

import ast
import pathlib

import pytest

from neoswga.core.unified_optimizer import (
    _ALL_SELECTION_WEIGHTS,
    _WEIGHTS_INERT_ON,
    _selection_weight_note,
)

_CORE = pathlib.Path(__file__).resolve().parent.parent / "neoswga" / "core"

#: The engine that makes the selection decision for each method, and the module
#: it lives in. `hybrid` and `background-aware` share one.
_SELECTING_ENGINE = {
    "hybrid": ("hybrid_optimizer", "HybridOptimizer"),
    "background-aware": ("hybrid_optimizer", "HybridOptimizer"),
    "network": ("network_optimizer", "NetworkOptimizer"),
    "dominating-set": ("dominating_set_optimizer", "DominatingSetOptimizer"),
    "clique": ("clique_optimizer", "CliqueOptimizer"),
}

_FIGURES = {"tm_weight": 0.30, "uniformity_weight": 0.05, "dimer_penalty": 0.45}


def _accepted_by(module, cls):
    """The weight names `cls.__init__` accepts at all."""
    tree = ast.parse((_CORE / f"{module}.py").read_text())
    for node in ast.walk(tree):
        if isinstance(node, ast.ClassDef) and node.name == cls:
            for item in node.body:
                if isinstance(item, ast.FunctionDef) and item.name == "__init__":
                    args = item.args
                    names = {a.arg for a in args.args} | {a.arg for a in args.kwonlyargs}
                    return {w for w in _ALL_SELECTION_WEIGHTS if w in names}
    raise AssertionError(f"{module}.{cls} has no __init__")


@pytest.mark.parametrize("method", sorted(_SELECTING_ENGINE))
def test_the_inert_map_matches_the_engines(method):
    """A weight the engine cannot even accept must be listed as inert.

    This is the half that keeps the message honest as the code moves. It does
    not assert the converse -- that an accepted weight is actually READ -- which
    is why `hybrid` carries two entries despite accepting all three: that half
    rests on the measurement in
    `docs/validation/selection_weights_are_inert_2026-09-21.md`.
    """
    module, cls = _SELECTING_ENGINE[method]
    accepted = _accepted_by(module, cls)
    claimed_applied = set(_ALL_SELECTION_WEIGHTS) - set(_WEIGHTS_INERT_ON.get(method, ()))

    overclaimed = sorted(claimed_applied - accepted)
    assert not overclaimed, (
        f"the log would report {overclaimed} as applied on '{method}', but "
        f"{cls} does not accept {overclaimed} at all. Add them to "
        f"_WEIGHTS_INERT_ON."
    )


def test_hybrid_does_not_claim_the_inert_weights_are_applied():
    note = _selection_weight_note("hybrid", _FIGURES)
    assert "applied:" not in note, f"still claims the weights were applied: {note}"
    assert "tm_weight" in note and "uniformity_weight" in note, note
    assert "no selection decision" in note, note
    # dimer_penalty IS read by hybrid's Stage 2, so it must not be disclaimed.
    assert "dimer_penalty reach no" not in note, note


def test_a_method_that_reads_none_of_them_says_so():
    for method in ("dominating-set", "clique"):
        note = _selection_weight_note(method, _FIGURES)
        assert "none of them reach any selection decision" in note, note


def test_network_still_reports_them_as_applied():
    """`network` scores with all three, so the plain message is correct there.

    Without this the fix could be "disclaim everywhere", which trades one
    inaccurate message for another.
    """
    note = _selection_weight_note("network", _FIGURES)
    assert note.startswith("selection weights applied:"), note
    assert "no selection decision" not in note, note


def test_the_resolved_figures_are_always_reported():
    """Whatever the method, the numbers the profile resolved are still shown.

    The point is to stop overclaiming, not to hide what was configured: a user
    debugging `--application` needs to see the figures either way.
    """
    for method in sorted(_SELECTING_ENGINE):
        note = _selection_weight_note(method, _FIGURES)
        assert "tm=0.30" in note and "uniformity=0.05" in note, (method, note)
        assert "dimer_penalty=0.45" in note, (method, note)
