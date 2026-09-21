"""Standalone contraction measures external oligos and respects panel limits."""

from dataclasses import replace
from types import SimpleNamespace

import pytest

from neoswga.core.base_optimizer import OptimizerConfig, PrimerSetMetrics
from neoswga.core.panel_contraction import contract_panel, scan_panel_positions
from neoswga.core.pool_objective import PoolConstraints

A, B = "ACGGACGGACGG", "AGGAGGAGGAGG"


def test_scans_do_not_match_oligos_across_record_joins(tmp_path):
    fasta = tmp_path / "two.fna"
    fasta.write_text(">first\nAAAAAA\n>second\nCCCCCC\n")
    cache = scan_panel_positions(["AAAACCCC", "AAAA"], [("fg", str(fasta), 12, False)])
    assert cache.get_positions("fg", "AAAACCCC", "both").size == 0
    assert cache.get_positions("fg", "AAAA", "forward").tolist() == [0, 1, 2]
    assert cache.get_record_starts("fg") == [0, 6]
    assert not list(tmp_path.glob("*.h5"))


def test_each_reference_keeps_its_own_circularity(tmp_path):
    fasta = tmp_path / "one.fna"
    fasta.write_text(">one\nAAAACCCC\n")
    cache = scan_panel_positions(
        ["CCAAAA"], [("target", str(fasta), 8, True), ("host", str(fasta), 8, False)]
    )
    assert cache.get_positions("target", "CCAAAA", "forward").tolist() == [6]
    assert cache.get_positions("host", "CCAAAA", "both").size == 0


@pytest.mark.parametrize("single_coverage,single_density", [(0.5, 20), (0.95, 2)])
def test_contraction_preserves_effective_coverage_and_specificity(single_coverage, single_density):
    def metrics(panel):
        return replace(
            PrimerSetMetrics.empty(),
            fg_coverage=0.99,
            effective_fg_coverage=single_coverage if len(panel) == 1 else 0.9,
            selectivity_density=single_density if len(panel) == 1 else 20,
            total_bg_sites=2,
        )

    optimizer = SimpleNamespace(
        compute_metrics=metrics,
        config=OptimizerConfig(),
        conditions=SimpleNamespace(fingerprint=lambda: "test"),
        bg_prefixes=["bg"],
        bg_seq_lengths=[1000],
    )
    result = contract_panel(optimizer, [A, B], 0.8, PoolConstraints(min_selectivity_density=10))
    assert result["contracted_set"] == [A, B]
    assert result["meets_target"]
    assert result["final_coverage"] == 0.9
    assert result["final_raw_coverage"] == 0.99
