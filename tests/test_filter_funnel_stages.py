"""The filtering funnel must name the stages that actually filter.

`after_thermodynamic` and `after_background` were set from the same dataframe
with only configuration-gated blocks between them, so they were equal on every
run without a blacklist -- 458/458, 1231/1231, 337/337, 20329/20329 on four
independent runs. The report labelled the inert one "After background/blacklist",
which read as a background result.

The real background gate is the bg_bool term folded into after_frequency, and
the max_primer cut -- the largest single reduction, 85% of survivors on one run
-- appeared only as final_candidates, under no stage name.
"""

import json

import pytest

from neoswga.core.report.metrics import FilteringStats


@pytest.fixture(scope="module")
def counted_workspace(tmp_path_factory):
    """Run count-kmers and filter once against a small synthetic target."""
    import os
    import random
    import shutil

    if not shutil.which("jellyfish"):
        pytest.skip("jellyfish not available (required for count-kmers)")

    root = tmp_path_factory.mktemp("funnel")
    rng = random.Random(20260906)
    seq = "".join(rng.choice("ACGT") for _ in range(60000))
    fasta = root / "target.fasta"
    fasta.write_text(
        ">target\n" + "\n".join(seq[i : i + 70] for i in range(0, len(seq), 70)) + "\n"
    )

    params = {
        "fg_genomes": [str(fasta)],
        "bg_genomes": [],
        "fg_prefixes": [str(root / "target")],
        "bg_prefixes": [],
        "data_dir": str(root / "results"),
        "min_k": 10,
        "max_k": 10,
        "polymerase": "phi29",
        "reaction_temp": 30.0,
        "min_fg_freq": 1e-6,
        "max_bg_freq": 1.0,
        "max_gini": 1.0,
        "max_primer": 40,
        "min_tm": 0,
        "max_tm": 100,
        "gc_min": 0.0,
        "gc_max": 1.0,
        "num_primers": 4,
        "target_set_size": 4,
        "max_sets": 2,
        "iterations": 2,
        "cpus": 1,
        "fg_circular": True,
        "schema_version": 2,
    }
    params_file = root / "params.json"
    params_file.write_text(json.dumps(params, indent=2))

    from neoswga import cli_unified

    previous = os.getcwd()
    os.chdir(root)
    try:
        for command, handler in (
            ("count-kmers", cli_unified.run_step1),
            ("filter", cli_unified.run_step2),
        ):
            args = cli_unified.create_parser().parse_args([command, "-j", str(params_file)])
            handler(args)
    finally:
        os.chdir(previous)

    return {"root": root, "data_dir": root / "results"}


def test_the_frequency_row_is_split_into_foreground_and_background():
    stats = FilteringStats(
        total_kmers=2542354,
        after_fg_frequency=180000,
        after_bg_frequency=1639,
        after_thermodynamic=458,
        after_gini=449,
        after_max_primer_cut=449,
        final_candidates=449,
    )
    labels = [label for label, _ in stats.as_funnel()]

    assert "After foreground frequency" in labels
    assert "After background frequency" in labels
    assert (
        "After background/blacklist" not in labels
    ), "the stage that could not filter is still labelled as a background result"


def test_the_max_primer_cut_has_a_stage_name():
    stats = FilteringStats(
        total_kmers=20329000,
        after_fg_frequency=100000,
        after_bg_frequency=20329,
        after_thermodynamic=20329,
        after_gini=20301,
        after_max_primer_cut=3000,
        final_candidates=3000,
    )
    funnel = dict(stats.as_funnel())

    assert funnel["After max_primer cut"] == 3000
    assert funnel["After Gini filter"] == 20301


def test_the_funnel_is_monotonically_non_increasing():
    stats = FilteringStats(
        total_kmers=2542354,
        after_fg_frequency=180000,
        after_bg_frequency=1639,
        after_thermodynamic=458,
        after_exclusion_blacklist=458,
        after_gini=449,
        after_max_primer_cut=449,
        final_candidates=449,
    )
    counts = [count for _, count in stats.as_funnel()]
    assert counts == sorted(counts, reverse=True), counts


def test_an_unconfigured_exclusion_stage_is_omitted():
    """Reporting a stage that is not configured is what made the old funnel
    read as if background filtering did nothing."""
    stats = FilteringStats(
        total_kmers=1000,
        after_fg_frequency=800,
        after_bg_frequency=500,
        after_thermodynamic=200,
        after_gini=190,
        after_max_primer_cut=100,
        final_candidates=100,
    )
    labels = [label for label, _ in stats.as_funnel()]
    assert "After exclusion/blacklist" not in labels


def test_a_legacy_filter_stats_file_still_renders(tmp_path):
    """Directories produced before this change carry the old key names."""
    from neoswga.core.report.metrics import collect_pipeline_metrics

    (tmp_path / "step4_improved_df.csv").write_text("sequence\nACGTACGTACGT\n")
    (tmp_path / "filter_stats.json").write_text(
        json.dumps(
            {
                "total_kmers": 2542354,
                "after_frequency": 1639,
                "after_thermodynamic": 458,
                "after_background": 458,
                "after_gini": 449,
                "final_candidates": 449,
            }
        )
    )

    metrics = collect_pipeline_metrics(str(tmp_path))

    assert metrics.filtering is not None
    assert metrics.filtering.total_kmers == 2542354
    assert metrics.filtering.final_candidates == 449
    funnel = dict(metrics.filtering.as_funnel())
    assert funnel["Total k-mers"] == 2542354
    assert funnel["Final candidates"] == 449


def test_a_new_filter_stats_file_renders_the_split_stages(tmp_path):
    """The mirror of the legacy case: a file written by the current filter step
    must not fall back to the legacy labels."""
    from neoswga.core.report.metrics import collect_pipeline_metrics

    (tmp_path / "step4_improved_df.csv").write_text("sequence\nACGTACGTACGT\n")
    (tmp_path / "filter_stats.json").write_text(
        json.dumps(
            {
                "total_kmers": 2542354,
                "after_fg_frequency": 180000,
                "after_bg_frequency": 1639,
                "after_thermodynamic": 458,
                "after_gini": 449,
                "after_max_primer_cut": 449,
                "final_candidates": 449,
            }
        )
    )

    funnel = dict(collect_pipeline_metrics(str(tmp_path)).filtering.as_funnel())

    assert funnel["After foreground frequency"] == 180000
    assert funnel["After background frequency"] == 1639
    assert funnel["After max_primer cut"] == 449
    assert "After frequency filter" not in funnel
    assert "After background/blacklist" not in funnel


@pytest.mark.parametrize(
    "key",
    [
        "total_kmers",
        "after_fg_frequency",
        "after_bg_frequency",
        "after_thermodynamic",
        "after_gini",
        "after_max_primer_cut",
        "final_candidates",
    ],
)
def test_the_filter_step_records_every_stage(counted_workspace, key):
    """An end-to-end filter run must write each stage the report can render."""
    data_dir = counted_workspace["data_dir"]
    stats = json.loads((data_dir / "filter_stats.json").read_text())
    assert key in stats, f"filter_stats.json is missing {key}: {sorted(stats)}"


def test_the_end_to_end_funnel_never_increases(counted_workspace):
    """The counts must be consistent with each other, not merely present.

    A stage written from the wrong dataframe passes the key-presence test above
    and still reports more survivors than the stage before it.
    """
    from neoswga.core.report.metrics import FilteringStats as _Stats

    stats = json.loads((counted_workspace["data_dir"] / "filter_stats.json").read_text())
    counts = [c for _, c in _Stats(**stats).as_funnel()]
    assert counts == sorted(counts, reverse=True), stats
