"""One implementation of "does the delivered panel break its own threshold".

The question -- is there a pair in the pool we are about to order whose
complementary run exceeds `max_dimer_bp` -- had two answers in the package:

- `dimer.dimer_validation_issue`, reached from `unified_optimizer`, which finds
  the worst pair through `worst_heterodimer` and names it.
- a nested loop written out again inside `panel_evaluation._panel_violations`.

Both compared `max_complementary_run(a, b) > max_dimer_bp` over every unordered
pair, so they agreed on every panel a resolved request can produce. That is what
makes removing one a refactor and not a change, and it is also the whole risk:
they agreed because they were written the same way, not because anything held
them to each other. This repository has already paid for that shape once, in
`string_search`, where two scanners computed one quantity two ways and disagreed
across a record join for as long as both existed.

`optimization_service.panel_violations` is deliberately NOT folded in, and the
last test here is the evidence rather than the assertion. It is the search-time
screen and asks a wider question -- self-dimers against `max_self_dimer_bp`, and
the optional `max_dimer_dg` free-energy floor -- and it answers with the opaque
string "dimer constraint" without naming a pair. A delivered-panel record that
adopted it would fault a pool on a threshold that does not govern the pool,
which is the reason `worst_heterodimer` already gives for excluding self-dimers.
"""

import ast
import itertools
import random
from pathlib import Path

from neoswga.core.design_request import resolve_design_request
from neoswga.core.dimer import dimer_validation_issue, max_complementary_run
from neoswga.core.dimer_validator import DimerValidator
from neoswga.core.panel_evaluation import evaluate_panel

MAPPING = {
    "fg_genomes": ["a.fna"],
    "fg_prefixes": ["a"],
    "polymerase": "phi29",
    "reaction_temp": 30.0,
    "min_k": 12,
    "max_k": 12,
    "num_primers": 12,
    "max_dimer_bp": 3,
}

# 10 of 12 bases complementary, the pair this project's dimer tests use.
DIMERISING = ["TGTACTGCCAAG", "CTTGGCAGTACA"]
CLEAN = ["ACCACAGATAGC", "GTTGTAGATGGA", "ATCAGCAGACCA"]


class Metrics:
    """The evaluator's own measurements, as the assessment receives them."""

    fg_coverage = 0.5
    bg_coverage = 0.01
    total_fg_sites = 100
    total_bg_sites = 5
    mean_gap = 1000.0
    max_gap = 5000.0
    gap_gini = 0.3
    mean_tm = 35.0
    selectivity_ratio = 20.0
    per_target_coverage: dict = {}
    effective_fg_coverage = 0.45


def _pools():
    """Panels a request can actually deliver, plus the two extremes.

    Random 12-mers under a fixed seed rather than a handful of chosen pools:
    the chosen ones are the cases whose answer the author already knows, and
    this check exists for the pair nobody thought of.
    """
    rng = random.Random(20260922)
    pools = [CLEAN, DIMERISING, CLEAN + DIMERISING, CLEAN[:1], []]
    for size in (2, 3, 6, 12):
        for _ in range(12):
            pools.append(["".join(rng.choice("ACGT") for _ in range(12)) for _ in range(size)])
    return pools


def _assessment_faults_a_dimer(pool, threshold):
    request = resolve_design_request(
        {**MAPPING, "num_primers": len(pool), "max_dimer_bp": threshold}
    )
    assessment = evaluate_panel(request, pool, Metrics())
    return any("heterodimer" in violation for violation in assessment.violations)


def test_the_two_answers_agree_on_every_pool():
    """The refactor is only safe while these cannot disagree, so pin it."""
    disagreements = []
    for pool in _pools():
        if not pool:
            continue
        for threshold in (1, 2, 3, 4, 5, 6, 7):
            shared = dimer_validation_issue(list(pool), threshold) is not None
            recorded = _assessment_faults_a_dimer(pool, threshold)
            if shared != recorded:
                disagreements.append((pool, threshold, shared, recorded))

    assert not disagreements, disagreements


def test_the_pool_with_no_pair_faults_nothing():
    """Guard the guard. An answer of "no violation" everywhere would pass the
    agreement check above without measuring anything."""
    assert not _assessment_faults_a_dimer(CLEAN[:1], 3)
    assert _assessment_faults_a_dimer(DIMERISING, 3)


def test_the_rule_has_one_implementation():
    """`panel_evaluation` must ask the shared rule rather than rewrite it.

    A source check, not a behavioural one, because the behaviour is identical
    today: the defect this prevents is the two drifting apart tomorrow, and by
    then a behavioural check has already gone red for a reason nobody planned.
    """
    source = (
        Path(__file__).resolve().parent.parent / "neoswga" / "core" / "panel_evaluation.py"
    ).read_text()
    tree = ast.parse(source)
    called = {
        node.func.id
        for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
    }

    assert "max_complementary_run" not in called, (
        "panel_evaluation is measuring complementary runs itself; the delivered-panel "
        "dimer rule belongs to dimer.dimer_validation_issue"
    )


def test_the_search_time_screen_is_a_wider_rule_and_stays_separate():
    """Evidence for why `panel_violations` is not the shared implementation.

    A lone self-dimering primer breaks the search-time screen and is not a
    delivered-panel heterodimer breach at all -- there is no pair to measure.
    Folding the two together would fault the pool on `max_self_dimer_bp`, an
    admission threshold applied during `filter`, under a finding that names
    `max_dimer_bp`.
    """
    lone = ["GCGCGCGCATAT"]

    assert DimerValidator(3, 4).has_self_dimer(lone[0])
    assert dimer_validation_issue(lone, 3) is None
    assert not _assessment_faults_a_dimer(lone, 3)


def test_the_shared_rule_still_measures_what_it_claims_to():
    """An independent count, so agreement is not two wrappers over one bug."""
    worst = max(
        (max_complementary_run(a, b) for a, b in itertools.combinations(DIMERISING, 2)),
        default=0,
    )

    assert worst >= 8
    assert dimer_validation_issue(DIMERISING, worst) is None
    assert dimer_validation_issue(DIMERISING, worst - 1) is not None
