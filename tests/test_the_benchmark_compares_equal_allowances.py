"""The comparison harness must give every arm the same allowance.

`sequential_panel_search.py` conceded in its own docstring that "budgets are per
stage, not equal total compute", and ran one unseeded configuration. Both flaws
are the same kind: a difference in the delivered panel could not be attributed
to the search rather than to the compute each arm was handed, or to the one seed
that happened to run.

A benchmark whose method drifts is worse than no benchmark, because its numbers
get quoted. These are source and behaviour checks on the replacement.
"""

import ast
import pathlib

import pytest

SCRIPT = pathlib.Path("scripts/benchmarking/equal_allowance_comparison.py")


@pytest.fixture(scope="module")
def tree():
    return ast.parse(SCRIPT.read_text())


@pytest.fixture(scope="module")
def source():
    return SCRIPT.read_text()


def test_every_arm_receives_one_shared_total_allowance(source):
    """A `SearchBudget` per arm, from the same `--evaluations` integer.

    Per-STAGE budgets are what the replaced script used. The give-away there was
    `swap_max_evaluations` and `beam_max_evaluations` in the optimizer config,
    which bound one loop each and reset on entry.
    """
    assert "SearchBudget(max_evaluations=evaluations)" in source
    assert "budget=budget" in source, "the allowance must reach the request"
    assert (
        "swap_max_evaluations" not in source
    ), "a per-stage budget in the harness reintroduces the flaw it replaces"
    assert "beam_max_evaluations" not in source


def test_the_allowance_is_required_rather_than_defaulted(tree):
    """A default would let two runs of the harness disagree silently."""
    calls = [
        node
        for node in ast.walk(tree)
        if isinstance(node, ast.Call)
        and getattr(node.func, "attr", None) == "add_argument"
        and node.args
        and isinstance(node.args[0], ast.Constant)
        and node.args[0].value == "--evaluations"
    ]
    assert calls, "the harness must take an explicit allowance"
    required = any(
        keyword.arg == "required" and keyword.value.value is True
        for call in calls
        for keyword in call.keywords
    )
    assert required, "--evaluations must be required, not defaulted"


def test_several_seeds_run_by_default(source):
    assert "default=[1, 2, 3]" in source, "one seed cannot show a tendency"
    assert "_seed_everything(seed)" in source


def test_proposal_generation_is_inside_the_timed_region(source):
    """The replaced script built one proposal outside the timing and shared it.

    The optimizer must be constructed and the search started after the clock,
    so an arm that spends its allowance generating proposals is charged for it.
    """
    body = source[source.index("def run_arm(") : source.index("def summarise(")]
    started = body.index("started = time.monotonic()")
    created = body.index("OptimizerFactory.create(")
    assert started < created, "the optimizer is built before the clock starts"


def test_memory_is_reported_once_and_says_why(source):
    """Per-arm RSS was measured, found meaningless, and removed.

    `ru_maxrss` never falls, so the first arm absorbs the shared cache build:
    1,112 MB against 54 MB for the same work on a smoke run.
    """
    assert "rss_delta_mb" not in source, "a per-arm memory delta is not a measurement here"
    assert "process_peak_rss_mb" in source
    assert "never falls" in source


def test_the_shared_warm_cache_is_disclosed(source):
    """Task 8 requires a fresh cache per arm or the sharing stated explicitly."""
    assert "shared_across_arms" in source
    # Whitespace-normalised, because prose wraps and a test that pins a line
    # break fails on reformatting rather than on a change of meaning.
    flat = " ".join(source.lower().split())
    assert "warm cache" in flat


def test_what_the_allowance_covers_is_recorded(source):
    """A count is uninterpretable without its scope."""
    assert "what_the_allowance_covers" in source
    assert "uncharged" in source


def test_the_output_states_that_coverage_is_a_proxy(source):
    assert "geometric proxy" in source
    assert "design_release_gates.md" in source


def test_an_arm_that_finishes_inside_its_allowance_is_distinguishable(source):
    """`stop_reason` separates "spent the budget" from "did not need it"."""
    assert '"stop_reason"' in source
    assert "stop_reasons" in source, "the summary must carry it too"


def test_distinct_panels_are_counted_so_a_median_is_not_oversold(source):
    """Three identical answers are determinism, not an uncertainty estimate."""
    assert "distinct_panels" in source


def test_a_failing_arm_is_reported_rather_than_aborting_the_sweep(source):
    assert '"failed"' in source
    assert '"failures"' in source, "the summary must count them"


def test_the_coverage_target_reaches_the_search(source):
    """`--target` was accepted, recorded in the output, and read by nothing.

    The first sweep therefore ran every arm at `OptimizationRequest`'s default
    of 0.7 while its JSON reported whatever was asked for. That is the
    inert-option defect this repository keeps a ratchet for, in the harness
    written to measure the search.
    """
    assert (
        "target_coverage=target_coverage" in source
    ), "the requested coverage target must reach OptimizationRequest"
    assert "target_coverage=args.target" in source, "and the CLI value must reach run_arm"
