"""An adapter must not drop a parameter on the way to the engine that decides.

Three of the five optimization methods are an outer adapter over an inner engine
that actually selects the panel. A value can reach the adapter's
`OptimizerConfig`, be asserted there by a test, and still be dropped when the
adapter constructs the inner -- and every one of the six reachability ratchets
named in CLAUDE.md passes, because they check the two ends and not the seam.

Four instances were found by hand on 2026-10-26 and a fifth on 2026-10-29:

- `stage2_tm_weight` reached `HybridBaseOptimizer`'s config and was never set on
  the inner, so the Stage-2 Tm term shipped days earlier was inert on every CLI
  `hybrid` run.
- `dimer_temp` was never forwarded by `DominatingSetAdapter` while
  `max_dimer_dg` beside it was, so the free-energy dimer screen ran at 37 C
  under every reaction.
- `polymerase` and `reaction_temp` reached no optimizer at all, so `hybrid` and
  `background-aware` ran every reaction as phi29 at 30 C.
- `max_dimer_dg`, `stage1_objective_width`, `stage2_tm_weight` and
  `stage2_background_weight` were all dropped by
  `BackgroundAwareBaseOptimizer`, which builds the same inner class as the
  `hybrid` wrapper from a second, separately maintained argument list.

The check this file encodes is mechanical and cheap: diff the inner's
`__init__` parameters against what the adapter passes or assigns. That is what
found the second and the fifth, so it is worth having a machine do it.

**What it does NOT catch**, stated so nobody reads a green run as the class
being closed (Known Issue 8 is never recorded as closed):

- A parameter passed with a REAL default instead of the configured value. That
  is the `polymerase=kwargs.get("polymerase", "phi29")` shape, which this rule
  sees as "passed" and which caused the widest of the five. It is pinned
  separately by `test_the_configured_reaction_reaches_the_hybrid_inner`.
- Anything about whether the inner then USES the value.
- Resolution logic inside the inner, such as `resolve_reaction`: it is not part
  of the adapter's argument list, so removing it leaves this test green. That
  case is pinned by the reaction test named above.
"""

import ast
import pathlib

import pytest

_CORE = pathlib.Path(__file__).resolve().parent.parent / "neoswga" / "core"

#: (adapter module, adapter class, inner class, module the inner is defined in).
#: `test_every_adapter_is_covered` fails if a class constructs an inner
#: optimizer and is not listed here, so a new adapter cannot escape the check.
ADAPTERS = [
    ("hybrid_optimizer", "HybridBaseOptimizer", "HybridOptimizer", "hybrid_optimizer"),
    (
        "background_aware_optimizer",
        "BackgroundAwareBaseOptimizer",
        "HybridOptimizer",
        "hybrid_optimizer",
    ),
    (
        "dominating_set_adapter",
        "DominatingSetAdapter",
        "DominatingSetOptimizer",
        "dominating_set_optimizer",
    ),
    ("network_optimizer", "NetworkBaseOptimizer", "NetworkOptimizer", "network_optimizer"),
    # Not a method adapter: `hybrid`'s own Stage 1 wraps the dominating-set
    # engine. Listed because it is the same kind of seam, and it is clean.
    (
        "hybrid_optimizer",
        "HybridOptimizer",
        "DominatingSetOptimizer",
        "dominating_set_optimizer",
    ),
    # Two further seams the completeness test below found on 2026-10-29. Both
    # carry entries in ALLOWED_UNREACHED; see the reasons there.
    ("primer_expansion", "PrimerExpander", "HybridOptimizer", "hybrid_optimizer"),
    (
        "primer_expansion",
        "PrimerExpander",
        "DominatingSetOptimizer",
        "dominating_set_optimizer",
    ),
    ("improved_pipeline", "ImprovedPipeline", "NetworkOptimizer", "network_optimizer"),
]

#: Inner parameters an adapter is currently known to leave at the inner's
#: default, keyed by `(adapter, inner, parameter)`. Each is a documented HOLE,
#: not a sanction: `test_no_allowed_entry_is_stale` fails on an entry that has
#: since been wired, so like the other reachability allowlists this one may only
#: shrink. The five method adapters have no entries and must keep none.
_EXPANSION = (
    "expand-primers builds a lean optimizer and does not forward the reaction, "
    "the Tm window or the dimer floor to it. Found 2026-10-29. Known Issue 14 "
    "records this path's own history of reading the background without acting "
    "on it; what a full forwarding changes to a delivered expansion is not "
    "measured, so it is recorded here rather than silently fixed."
)
_IMPROVED = (
    "improved_pipeline is reached only from core/validation.py, not from any "
    "design command, and constructs a NetworkOptimizer for a utility call. Not "
    "on the design path, so a dropped parameter here moves no delivered panel."
)

ALLOWED_UNREACHED: dict[tuple[str, str, str], str] = {
    **{
        ("PrimerExpander", "HybridOptimizer", name): _EXPANSION
        for name in (
            "absolute_coverage_floor",
            "allow_dimer_relaxation",
            "background_weight",
            "dimer_penalty",
            "genome_gc_content",
            "max_dimer_dg",
            "max_tm",
            "mechanistic_weight",
            "min_coverage_drop",
            "min_coverage_threshold",
            "min_tm",
            "polymerase",
            "reaction_temp",
            "stage1_objective_width",
            "stage2_tm_weight",
            "swap_max_evaluations",
            "swap_max_seconds",
            "template_gc",
            "tm_weight",
            "uniformity_weight",
        )
    },
    **{
        ("PrimerExpander", "DominatingSetOptimizer", name): _EXPANSION
        for name in (
            "allow_dimer_relaxation",
            "dimer_temp",
            "max_dimer_dg",
            "stage1_objective_width",
        )
    },
    **{
        ("ImprovedPipeline", "NetworkOptimizer", name): _IMPROVED
        for name in (
            "allow_dimer_relaxation",
            "bg_prefixes",
            "bg_seq_lengths",
            "conditions",
            "fg_prefixes",
            "fg_seq_lengths",
            "max_dimer_dg",
            "mechanistic_weight",
            "position_cache",
            "template_gc",
        )
    },
}

#: The adapters the factory returns for a `--optimization-method`, plus
#: hybrid's own Stage 1. These decide a delivered panel and may have no
#: allowlist entry at all.
STRICT_ADAPTERS = {
    "HybridBaseOptimizer",
    "BackgroundAwareBaseOptimizer",
    "DominatingSetAdapter",
    "NetworkBaseOptimizer",
    "HybridOptimizer",
}


def _tree(module):
    return ast.parse((_CORE / f"{module}.py").read_text())


def _init_parameters(module, cls):
    """Every named parameter of `cls.__init__`, positional and keyword-only."""
    for node in ast.walk(_tree(module)):
        if isinstance(node, ast.ClassDef) and node.name == cls:
            for item in node.body:
                if isinstance(item, ast.FunctionDef) and item.name == "__init__":
                    args = item.args
                    return [a.arg for a in args.args if a.arg != "self"] + [
                        a.arg for a in args.kwonlyargs
                    ]
    raise AssertionError(f"{module}.{cls} has no __init__; update ADAPTERS")


def _reached(adapter_module, adapter_cls, inner_cls):
    """Inner parameter names the adapter supplies, however it supplies them.

    Three routes count, because all three are in use: a keyword in the
    construction call, a `setattr` with a literal name (including the names in a
    `for` loop's literal tuple, which is how `HybridBaseOptimizer` sets the two
    Stage-2 weights), and a plain attribute assignment on the instance.
    """
    passed, assigned, splat = set(), set(), False
    for node in ast.walk(_tree(adapter_module)):
        if not (isinstance(node, ast.ClassDef) and node.name == adapter_cls):
            continue
        for call in ast.walk(node):
            if isinstance(call, ast.Call) and getattr(call.func, "id", None) == inner_cls:
                passed |= {kw.arg for kw in call.keywords if kw.arg}
                splat = splat or any(kw.arg is None for kw in call.keywords)
            if isinstance(call, ast.Call) and getattr(call.func, "id", None) == "setattr":
                if len(call.args) >= 2 and isinstance(call.args[1], ast.Constant):
                    assigned.add(call.args[1].value)
            if isinstance(call, ast.For) and isinstance(call.iter, ast.Tuple):
                assigned |= {e.value for e in call.iter.elts if isinstance(e, ast.Constant)}
            if isinstance(call, ast.Assign):
                assigned |= {t.attr for t in call.targets if isinstance(t, ast.Attribute)}
    return passed, assigned, splat


@pytest.mark.parametrize(
    "adapter_module,adapter_cls,inner_cls,inner_module",
    ADAPTERS,
    ids=[a[1] for a in ADAPTERS],
)
def test_no_inner_parameter_is_dropped(adapter_module, adapter_cls, inner_cls, inner_module):
    parameters = _init_parameters(inner_module, inner_cls)
    passed, assigned, splat = _reached(adapter_module, adapter_cls, inner_cls)
    if splat:
        pytest.skip(
            f"{adapter_cls} forwards **kwargs into {inner_cls}, so this rule "
            "cannot tell which parameters arrive"
        )

    unreached = [
        p
        for p in parameters
        if p not in passed
        and p not in assigned
        and (adapter_cls, inner_cls, p) not in ALLOWED_UNREACHED
    ]

    assert not unreached, (
        f"{adapter_cls} constructs {inner_cls} without supplying "
        f"{unreached}, so each of those silently takes the inner's default and "
        f"no configuration can reach it. Pass it in the constructor call, assign "
        f"it on the instance, or -- only if it genuinely cannot reach a "
        f"delivered panel -- add it to ALLOWED_UNREACHED with a reason."
    )


def test_a_method_adapter_has_no_allowlist_entry():
    """The adapters that decide a delivered panel carry no documented holes.

    `ALLOWED_UNREACHED` exists for seams off the design path. An entry for one
    of the five adapters the factory returns would mean a configured value
    cannot reach the engine that picks the panel, which is the whole defect
    class this file guards.
    """
    offenders = sorted({(a, i, p) for (a, i, p) in ALLOWED_UNREACHED if a in STRICT_ADAPTERS})
    assert (
        not offenders
    ), f"these are method adapters and may not have allowlist entries: {offenders}"


def test_no_allowed_entry_is_stale():
    """An entry that has since been wired must be removed, so the list shrinks.

    Without this the allowlist silently sanctions parameters that are now
    forwarded, and the next genuine drop hides among them -- the same way the
    six reachability ratchets' allowlists are required to only shrink.
    """
    stale = []
    for adapter_module, adapter_cls, inner_cls, inner_module in ADAPTERS:
        parameters = set(_init_parameters(inner_module, inner_cls))
        passed, assigned, _ = _reached(adapter_module, adapter_cls, inner_cls)
        for key in ALLOWED_UNREACHED:
            if key[:2] != (adapter_cls, inner_cls):
                continue
            name = key[2]
            if name not in parameters or name in passed or name in assigned:
                stale.append(key)
    assert not stale, (
        f"remove these ALLOWED_UNREACHED entries; they are forwarded now or no "
        f"longer parameters: {sorted(stale)}"
    )


def test_every_adapter_is_covered():
    """A new adapter must not escape the check by not being listed.

    Scans `core/` for a class that constructs one of the known inner engines and
    asserts it appears in ADAPTERS. Without this the ratchet silently covers
    only what someone remembered to add.
    """
    inner_names = {inner for _, _, inner, _ in ADAPTERS}
    listed = {(mod, cls) for mod, cls, _, _ in ADAPTERS}

    found = set()
    for path in sorted(_CORE.glob("*.py")):
        try:
            tree = ast.parse(path.read_text())
        except SyntaxError:  # pragma: no cover
            continue
        for node in ast.walk(tree):
            if not isinstance(node, ast.ClassDef):
                continue
            for call in ast.walk(node):
                if (
                    isinstance(call, ast.Call)
                    and getattr(call.func, "id", None) in inner_names
                    # The engine's own module may construct it in a factory or a
                    # test helper; only a DIFFERENT class wrapping it counts.
                    and node.name != getattr(call.func, "id", None)
                ):
                    found.add((path.stem, node.name))

    missing = sorted(found - listed)
    assert not missing, (
        f"these classes construct an inner optimizer and are not in ADAPTERS: "
        f"{missing}. Add each one, so the seam is checked."
    )
