"""Small reproductions for the follow-up audit; no production files are changed.

Run from the repository root:
    PYTHONPATH=. python docs/validation/pipeline_audit_2026-09-16/probes.py
"""

import json
import tempfile
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

from neoswga.core.bam_coverage import compute_bam_depth, match_contigs
from neoswga.core.base_optimizer import OptimizationStatus, OptimizerConfig
from neoswga.core.candidate_inventory import CandidateInventory, record_stage2_inventory
from neoswga.core.candidate_provider import CandidateProvider
from neoswga.core.filter import check_gini_stage_kept_something
from neoswga.core.kmer_counter import genome_fingerprint
from neoswga.core.pool_planner import plan_pool
from neoswga.core.primer_expansion import CoverageGap, PrimerExpander
from neoswga.core.reach_calibration import fit_reach

A, B = "ACGGACGGACGG", "AGGAGGAGGAGG"


def main():
    results = {}
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        both = pd.DataFrame({"primer": [A, B]})
        only_a = pd.DataFrame({"primer": [A]})
        path = record_stage2_inventory(root, "same-condition", both, both, only_a)
        # Simulate a stricter QC rerun with unchanged reaction conditions.
        record_stage2_inventory(root, "same-condition", only_a, only_a, only_a)
        with CandidateInventory(path) as inventory:
            results["stale_eligible_after_stricter_rerun"] = list(
                inventory.iter_eligible("same-condition", [12])
            )
        assert B in results["stale_eligible_after_stricter_rerun"]

        # An all-QC inventory under post_gini retains unindexed candidates too.
        post = root / "post"
        path = record_stage2_inventory(post, "condition", both, only_a, only_a, indexed=[A])
        with CandidateInventory(path) as inventory:
            provider = CandidateProvider(inventory, "condition", [12])
            results["provider_expands_to_unindexed_candidate"] = provider.expand([A], 10)
            assert B in results["provider_expands_to_unindexed_candidate"]

        reference = root / "ref.fa"
        reference.write_bytes(b">ref\n" + b"A" * (3 * 1024 * 1024))
        before = genome_fingerprint(str(reference))
        with reference.open("r+b") as handle:
            handle.seek(1500000)
            handle.write(b"C")
        results["same_fingerprint_after_interior_edit"] = before == genome_fingerprint(str(reference))
        assert results["same_fingerprint_after_interior_edit"]

        import pysam

        bam_path = root / "reads.bam"
        header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "target", "LN": 100}]}
        with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
            for i, (flag, mapq) in enumerate([(0, 60), (0, 0), (2048, 60)]):
                read = pysam.AlignedSegment()
                read.query_name = f"read{i}"
                read.query_sequence = "A" * 10
                read.query_qualities = pysam.qualitystring_to_array("I" * 10)
                read.flag = flag
                read.reference_id = 0
                read.reference_start = 10
                read.mapping_quality = mapq
                read.cigar = [(0, 10)]
                bam.write(read)
        pysam.index(str(bam_path))
        results["bam_depth_includes_mapq_zero_and_supplementary"] = int(
            compute_bam_depth(str(bam_path), "target", 100)[12]
        )
        assert results["bam_depth_includes_mapq_zero_and_supplementary"] == 3

    try:
        check_gini_stage_kept_something(pd.DataFrame({"primer": [A]}), pd.DataFrame())
    except ValueError as error:
        results["gini_empty_blocks_before_inventory"] = str(error)

    class Optimizer:
        name = "audit-fixed-return"
        bg_prefixes = ["bg"]
        bg_seq_lengths = [1000]
        conditions = object()
        config = OptimizerConfig(max_dimer_bp=3, max_self_dimer_bp=4)

        def optimize(self, candidates, target_size):
            return SimpleNamespace(primers=[A], status=OptimizationStatus.SUCCESS, message="")

        def compute_metrics(self, primers):
            coverage = 0.95 if B in primers else 0.5
            return SimpleNamespace(
                fg_coverage=coverage, effective_fg_coverage=coverage,
                selectivity_density=20.0, total_bg_sites=1, max_gap=100,
            )

    plan = plan_pool(Optimizer(), [A, B], [1], [0.9], min_selectivity_density=10)
    results["coverage_miss_no_repair"] = {
        "recommendation": plan["recommendations"][0],
        "repair": plan["rows"][0]["repair"],
        "available_alternative_coverage": 0.95,
    }
    assert not plan["rows"][0]["repair"]["attempted"]
    assert plan["recommendations"][0]["status"] == "not_found"

    class Cache:
        def get_positions(self, prefix, primer, strand="both"):
            return np.array([900] if primer == A else [1500])

    expander = PrimerExpander(Cache(), ["fg"], [10000], coverage_reach=3000)
    gap = CoverageGap(chromosome="fg", start=1000, end=2000, size=1000)
    kept = expander._filter_candidates_to_gaps([A, B], [gap])
    results["gap_screen_excludes_nearby_reaching_primer"] = kept
    assert A not in kept  # Its reach covers the entire gap, but its site is outside.

    results["multi_record_bam_mapping"] = match_contigs(
        ["contig1", "contig2"], [1000, 2000], ["target"], [3000]
    )
    assert results["multi_record_bam_mapping"] == {}

    try:
        fit_reach(np.ones(500), [100], bin_size=1000)
    except ValueError as error:
        results["short_contig_reach_fit_error"] = str(error)
    assert "short_contig_reach_fit_error" in results
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
