"""Regression test: --auto-size dispatch must use realistic per-primer reach.

Prior to Phase 16, the CLI dispatch passed theoretical polymerase processivity
(phi29 = 70 kb) to recommend_set_size while base_optimizer measures coverage
with realistic per-primer reach (3 kb). The mismatch caused the recommender
to under-count required primers by an order of magnitude.

This test exercises the underlying helper and confirms the recommender is
sensitive to the chosen reach value, so a future regression that swaps the
helper out for a fixed 70 kb default is caught.
"""

from neoswga.core.coverage import polymerase_extension_reach
from neoswga.core.set_size_optimizer import recommend_set_size, create_baseline_effects


def test_realistic_reach_values_match_phase16_design():
    assert polymerase_extension_reach("phi29", coverage_metric="realistic") == 3000
    assert polymerase_extension_reach("equiphi29", coverage_metric="realistic") == 4000
    assert polymerase_extension_reach("phi29", coverage_metric="processivity") == 70000


def test_recommend_set_size_responds_to_reach_choice():
    """A smaller per-primer reach must not yield a smaller recommendation."""
    effects = create_baseline_effects()
    common = dict(
        application="discovery",
        genome_length=5_000_000,
        primer_length=15,
        mech_effects=effects,
    )

    rec_realistic = recommend_set_size(processivity=3000, **common)
    rec_processivity = recommend_set_size(processivity=70000, **common)

    # The recommended size is bounded to [4, 20] in the current implementation,
    # so we cannot assert a strict ratio. We only require monotonicity:
    # smaller reach must require at least as many primers.
    assert rec_realistic["recommended_size"] >= rec_processivity["recommended_size"]


def test_cli_dispatch_does_not_hardcode_processivity_map():
    """Source-level guard: the auto-size dispatch must derive per-primer reach
    from polymerase_extension_reach(), not a hardcoded processivity_map.
    """
    from pathlib import Path

    # The optimize/auto-size dispatch was extracted from cli_unified.py into the
    # neoswga/cli/ package; scan the whole CLI surface so this guard survives
    # the ongoing god-file split.
    cli_dir = Path(__file__).resolve().parents[1] / "neoswga"
    cli_sources = [cli_dir / "cli_unified.py", *(cli_dir / "cli").glob("*.py")]
    cli_src = "\n".join(p.read_text() for p in cli_sources if p.exists())

    assert "processivity_map = {'phi29': 70000" not in cli_src, (
        "Hardcoded processivity_map at 70000 bp under-counts required primers "
        "by ~10-20x (it is theoretical phi29 processivity, not the realistic "
        "per-primer reach base_optimizer measures). Use polymerase_extension_reach()."
    )
    assert "polymerase_extension_reach" in cli_src, (
        "Auto-size dispatch must import polymerase_extension_reach to derive "
        "the right per-primer reach for the recommender."
    )
