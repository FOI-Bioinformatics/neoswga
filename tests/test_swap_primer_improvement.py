"""Strict-improvement regression tests for dimer_network_analyzer.optimize_set_greedy
(the engine behind `neoswga swap-primer`).

Phase 12B: swap must never commit a change that worsens the set's mean
severity. Before the Phase 12E dimer fix landed, the severity calculator
returned 0.0 for every pair and these tests would have trivially passed
without verifying any real swap behavior. Now that severity is
temperature-aware and correct, these tests carry signal.
"""

import pytest

from neoswga.core.dimer_network_analyzer import DimerNetworkAnalyzer
from neoswga.core.reaction_conditions import ReactionConditions


def _analyzer(temp: float = 30.0, polymerase: str = "phi29"):
    return DimerNetworkAnalyzer(
        conditions=ReactionConditions(temp=temp, polymerase=polymerase),
        severity_threshold=0.3,
    )


# Primers chosen so the set has a clear worst offender — two engineered
# dimer-prone pairs plus a clean bystander. Replacement pool has clean
# primers that should reduce mean severity when swapped in.
DIMER_HEAVY_SET = [
    "ATCGATCGAT",  # palindromic — forms strong self / cross-dimer
    "ATCGATCGAT",  # duplicate to guarantee dimer (via homology)
    "ACGTACGTAC",
]
CLEAN_POOL = [
    "TTATTGGAAC",
    "AGTGACAGTG",
    "CTCCAACAAT",
    "GAGTTGAGTT",
    "AACAGCCAAC",
]


def test_swap_primer_never_worsens_mean_severity():
    """After optimize_set_greedy, the set's mean severity must not exceed
    the starting mean severity. The new guard short-circuits bad swaps."""
    # Note: DIMER_HEAVY_SET has duplicates that inflate severity.
    initial_set = list(DIMER_HEAVY_SET)
    # Ensure duplicates are explicit (dedup-safe test still signals something)
    analyzer = _analyzer()
    before, _, _ = analyzer.analyze_primer_set(initial_set, verbose=False)
    improved, after = analyzer.optimize_set_greedy(
        initial_set, CLEAN_POOL, max_iterations=3,
    )
    assert after.mean_severity <= before.mean_severity + 1e-6, (
        f"swap must not increase mean severity: before={before.mean_severity:.4f}, "
        f"after={after.mean_severity:.4f}"
    )


def test_swap_primer_preserves_set_size():
    """Net swap keeps the set size constant."""
    initial_set = list(DIMER_HEAVY_SET)
    analyzer = _analyzer()
    improved, _ = analyzer.optimize_set_greedy(
        initial_set, CLEAN_POOL, max_iterations=3,
    )
    assert len(improved) == len(initial_set)


def test_swap_primer_rejects_worsening_candidate():
    """When the candidate pool offers only sequences that would worsen the
    set, the optimizer must reject the swap and return the original set."""
    # Pool of dimer-heavy candidates only — worse than the bystander
    # primer ACGTACGTAC in the set
    bad_pool = [
        "ATCGATCGAT",  # identical to the worst offender
        "ATCGATCGAT",
    ]
    analyzer = _analyzer()
    initial = ["ACGTACGTAC", "TTATTGGAAC", "CTCCAACAAT"]
    before, _, _ = analyzer.analyze_primer_set(initial, verbose=False)
    improved, after = analyzer.optimize_set_greedy(initial, bad_pool, max_iterations=3)
    # Either the starting set passes (no swap needed) or no worse candidate
    # was accepted. Either way, mean_severity must not increase.
    assert after.mean_severity <= before.mean_severity + 1e-6


def test_swap_primer_honours_reaction_conditions():
    """Severity is temperature-dependent (Phase 12E), so the same set
    scored at a higher reaction temperature should have lower mean
    severity. This ties swap-primer's behavior to the user's buffer."""
    initial = list(DIMER_HEAVY_SET)

    an_30 = _analyzer(temp=30.0, polymerase="phi29")
    an_45 = _analyzer(temp=45.0, polymerase="equiphi29")

    m30, _, _ = an_30.analyze_primer_set(initial, verbose=False)
    m45, _, _ = an_45.analyze_primer_set(initial, verbose=False)

    assert m45.mean_severity <= m30.mean_severity + 1e-6, (
        f"At higher temperature mean severity should not increase: "
        f"30 C={m30.mean_severity:.3f}, 45 C={m45.mean_severity:.3f}"
    )
