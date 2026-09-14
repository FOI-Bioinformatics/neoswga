"""Pool recommendations require coverage, specificity and compatibility together."""

from types import SimpleNamespace

import pytest

from neoswga.core.base_optimizer import OptimizationStatus, OptimizerConfig
from neoswga.core.pool_planner import plan_pool

A, C, AC = "AAAAAAAAAAAA", "CCCCCCCCCCCC", "ACACACACACAC"


class Optimizer:
    name = "test"
    bg_prefixes = ["bg"]
    bg_seq_lengths = [1000]
    conditions = object()
    config = OptimizerConfig()

    def __init__(self, panels=None):
        self.panels = panels or {1: [A], 2: [A, C], 3: [A, C, AC]}
        self.metrics = {1: (0.6, 50, 5), 2: (0.92, 4, 30), 3: (0.96, 20, 10)}
        self.called = []

    def optimize(self, candidates, target_size):
        self.called.append((target_size, candidates))
        return SimpleNamespace(
            primers=self.panels[target_size], status=OptimizationStatus.PARTIAL, message=""
        )

    def compute_metrics(self, primers):
        coverage, density, bg = self.metrics[len(primers)]
        return SimpleNamespace(
            effective_fg_coverage=coverage,
            fg_coverage=coverage,
            selectivity_density=density,
            total_bg_sites=bg,
            max_gap=100,
        )


def test_smallest_pool_must_meet_both_coverage_and_specificity():
    opt = Optimizer()
    result = plan_pool(opt, [A, C, AC], [3, 1, 2], [0.5, 0.9, 0.99], min_selectivity_density=10)
    assert [r["size"] for r in result["recommendations"]] == [1, 3, None]
    assert result["rows"][1]["failed_constraints"] == ["selectivity below minimum"]
    assert [n for n, _ in opt.called] == [1, 2, 3]


def test_exact_background_limit_can_be_used_separately():
    result = plan_pool(Optimizer(), [A, C, AC], [1, 2, 3], [0.9], max_background_sites=5)
    assert result["recommendations"][0]["status"] == "not_found"


def test_partial_panel_is_measured_at_delivered_size():
    opt = Optimizer({4: [A]})
    result = plan_pool(opt, [A, C, AC], [4], [0.5], min_selectivity_density=10)
    assert result["rows"][0]["requested_size"] == 4
    assert result["recommendations"][0]["size"] == 1


def test_no_specificity_claim_without_background():
    opt = Optimizer()
    opt.bg_prefixes = []
    opt.bg_seq_lengths = []
    with pytest.raises(ValueError, match="background genome"):
        plan_pool(opt, [A], [1], [0.5], min_selectivity_density=10)
    result = plan_pool(opt, [A], [1], [0.5])
    assert result["background_assessed"] is False
    assert result["rows"][0]["selectivity_density"] is None


def test_background_requires_explicit_specificity_target():
    with pytest.raises(ValueError, match="Specify"):
        plan_pool(Optimizer(), [A], [1], [0.5])


def test_uses_only_requested_length():
    opt = Optimizer()
    plan_pool(opt, [A, "AAA"], [1], [0.5], min_selectivity_density=10)
    assert opt.called[0][1] == [A]


def test_non_finite_target_is_rejected():
    with pytest.raises(ValueError, match="Coverage targets"):
        plan_pool(Optimizer(), [A], [1], [float("nan")], min_selectivity_density=10)


@pytest.mark.parametrize("coverage,density", [(float("nan"), 20), (0.9, float("nan"))])
def test_non_finite_metrics_cannot_produce_a_recommendation(coverage, density):
    opt = Optimizer()
    opt.metrics[1] = (coverage, density, 0)
    with pytest.raises(ValueError, match="Optimizer returned"):
        plan_pool(opt, [A], [1], [0.5], min_selectivity_density=10)


def test_dimer_pair_cannot_qualify():
    t = "TTTTTTTTTTTT"
    opt = Optimizer({2: [A, t]})
    result = plan_pool(opt, [A, t], [2], [0.5], min_selectivity_density=1)
    assert not result["rows"][0]["eligible"]
    assert result["recommendations"][0]["size"] is None


def test_exports_only_qualifying_panels(tmp_path):
    from neoswga.core.pool_plan_report import write_pool_plan

    result = plan_pool(Optimizer(), [A, C, AC], [1, 2, 3], [0.9, 0.99], min_selectivity_density=10)
    path = write_pool_plan(result, tmp_path)
    assert path.exists()
    # The figure needs matplotlib, which lives in the `viz` extra rather than
    # the base install, so it is drawn when available and skipped when not.
    # See test_the_report_is_written_without_matplotlib.
    try:
        import matplotlib  # noqa: F401

        drawable = True
    except ImportError:
        drawable = False
    if drawable:
        assert (tmp_path / "pool_sizes.png").stat().st_size > 0
    assert (tmp_path / "target_90pct_oligos.fasta").read_text().count(">") == 3
    assert not (tmp_path / "target_99pct_oligos.fasta").exists()


def test_cli_registers_pool_planner():
    from neoswga.cli_unified import create_parser

    args = create_parser().parse_args(
        ["plan-pool", "-j", "params.json", "--min-selectivity-density", "10"]
    )
    assert args.primer_length == 12
    assert args.coverage_targets == [0.9, 0.95]


def test_report_pool_renders_saved_results_without_genome_files(tmp_path):
    import json

    from neoswga.cli.plan_pool import run_report_pool
    from neoswga.cli_unified import create_parser

    plan = plan_pool(Optimizer(), [A, C, AC], [1, 2, 3], [0.9], min_selectivity_density=10)
    plan["design_parameters"] = {
        "fg_genomes": ["/missing/arbitrary_target.fna"],
        "bg_genomes": ["/missing/arbitrary_background.fna"],
    }
    source = tmp_path / "pool_plan.json"
    source.write_text(json.dumps(plan))
    original = source.read_bytes()
    output = tmp_path / "report"
    args = create_parser().parse_args(
        [
            "report-pool",
            "--input",
            str(tmp_path),
            "--title",
            "Target <comparison>",
            "-o",
            str(output),
        ]
    )
    run_report_pool(args)
    page = (output / "pool_plan.html").read_text()
    assert "arbitrary_target.fna" in page
    assert "arbitrary_background.fna" in page
    assert "Target &lt;comparison&gt;" in page
    assert (output / "target_90pct_oligos.fasta").exists()
    assert source.read_bytes() == original


def test_report_preserves_existing_outputs(tmp_path):
    from neoswga.core.pool_plan_report import write_pool_plan

    sentinel = tmp_path / "previous.fasta"
    sentinel.write_text("previous design")
    plan = plan_pool(Optimizer(), [A], [1], [0.5], min_selectivity_density=10)
    with pytest.raises(ValueError, match="not empty"):
        write_pool_plan(plan, tmp_path)
    assert sentinel.read_text() == "previous design"


def test_report_rejects_non_plan_json_before_writing(tmp_path):
    from neoswga.core.pool_plan_report import write_pool_plan

    destination = tmp_path / "report"
    with pytest.raises(ValueError, match="saved plan-pool"):
        write_pool_plan({"fg_genomes": []}, destination)
    assert not destination.exists()


def _saved_plan():
    """A well-formed plan, as `plan-pool` would have written it."""
    return plan_pool(Optimizer(), [A, C, AC], [1, 2, 3], [0.9], min_selectivity_density=10)


def test_report_escapes_every_value_taken_from_the_saved_plan(tmp_path):
    """`report-pool` renders a JSON file the tool did not necessarily write.

    The title was escaped and the scalar fields around it were not, so a plan
    carrying markup in `coverage_metric`, `primer_length` or a specificity limit
    put that markup straight into the page. Nothing here is a privilege
    boundary, but the sibling report already escapes and a report that alters
    its own structure on hostile input cannot be trusted to show what it was
    given.
    """
    from neoswga.core.pool_plan_report import write_pool_plan

    plan = _saved_plan()
    plan["coverage_metric"] = "<script>alert(1)</script>"
    plan["primer_length"] = "<img src=x onerror=alert(2)>"
    plan["min_selectivity_density"] = "<b>one</b>"
    plan["max_background_sites"] = "<b>two</b>"
    plan["rows"][0]["requested_size"] = "<b>three</b>"
    plan["rows"][0]["background_sites"] = "<b>four</b>"

    page = write_pool_plan(plan, tmp_path / "report").read_text()

    assert "<script>" not in page
    assert "<img src=x" not in page
    assert "<b>" not in page
    for escaped in ("&lt;script&gt;", "&lt;img src=x", "&lt;b&gt;one&lt;/b&gt;"):
        assert escaped in page


def test_report_rejects_a_plan_whose_rows_are_malformed(tmp_path):
    """A row missing a field must fail like any other unusable plan.

    It used to surface as a bare KeyError, which the CLI reports as a missing
    *parameter* and answers with advice about `neoswga init` and params.json --
    neither of which has anything to do with the file actually at fault.
    """
    from neoswga.core.pool_plan_report import write_pool_plan

    plan = _saved_plan()
    plan["rows"][0].pop("eligible")
    destination = tmp_path / "report"

    with pytest.raises(ValueError, match="saved plan-pool"):
        write_pool_plan(plan, destination)
    assert not destination.exists()


def test_report_rejects_a_recommendation_pointing_outside_the_rows(tmp_path):
    """`row_index` is an index into `rows`; an unusable one must say so."""
    from neoswga.core.pool_plan_report import write_pool_plan

    plan = _saved_plan()
    plan["recommendations"][0]["row_index"] = 9999
    destination = tmp_path / "report"

    with pytest.raises(ValueError, match="saved plan-pool"):
        write_pool_plan(plan, destination)
    assert not destination.exists()


def test_a_row_carrying_failed_constraints_does_not_also_need_a_status(tmp_path):
    """The fallback was an eager default, so `status` was read either way."""
    from neoswga.core.pool_plan_report import write_pool_plan

    plan = _saved_plan()
    row = plan["rows"][0]
    row["eligible"], row["failed_constraints"] = False, ["selectivity below minimum"]
    row.pop("status")

    page = write_pool_plan(plan, tmp_path / "report").read_text()
    assert "selectivity below minimum" in page


def test_report_says_what_is_wrong_when_the_output_path_is_a_file(tmp_path):
    """ "Directory is not empty" describes the wrong problem for a file."""
    from neoswga.core.pool_plan_report import write_pool_plan

    occupied = tmp_path / "design_report"
    occupied.write_text("not a directory")

    with pytest.raises(ValueError, match="not a directory"):
        write_pool_plan(_saved_plan(), occupied)
    assert occupied.read_text() == "not a directory"


def test_report_pool_names_the_file_it_could_not_find(tmp_path):
    """The shared handler answers FileNotFoundError with genome-file advice.

    This command reads no genome and no params.json, so the message has to
    carry the whole diagnosis itself.
    """
    from neoswga.cli.plan_pool import run_report_pool
    from neoswga.cli_unified import create_parser

    args = create_parser().parse_args(
        ["report-pool", "--input", str(tmp_path), "-o", str(tmp_path / "report")]
    )
    with pytest.raises(FileNotFoundError, match="plan-pool"):
        run_report_pool(args)
    assert not (tmp_path / "report").exists()


def test_report_refuses_to_export_something_that_is_not_an_oligo(tmp_path):
    """The FASTA is what goes to a synthesis vendor.

    `plan_pool` already rejects a candidate that is not A/C/G/T, but
    `report-pool` renders a saved file and never goes through it. A primer
    string carrying a newline split into a second, invented record, and a
    string that is not a sequence at all was exported as though it were.
    """
    from neoswga.core.pool_plan_report import write_pool_plan

    plan = _saved_plan()
    index = next(r["row_index"] for r in plan["recommendations"] if r["row_index"] is not None)
    plan["rows"][index]["primers"] = [A, "AAAA\n>injected_contig\nTTTT"]
    destination = tmp_path / "report"

    with pytest.raises(ValueError, match="A, C, G and T"):
        write_pool_plan(plan, destination)
    assert not destination.exists()


def test_the_report_is_written_without_matplotlib(tmp_path, monkeypatch):
    """matplotlib is in the `viz` extra, not the base install.

    `pool_plan_report` imported it unconditionally, so `plan-pool` and
    `report-pool` raised ModuleNotFoundError for anyone who had not installed
    that extra. Continuous integration installs only the dev extra and caught
    it; a local environment with matplotlib present could not.

    The figure is the only output that needs plotting. Everything a reader acts
    on -- the recommended panel sizes, the oligo FASTAs, the CSV and the JSON --
    is computed without it.
    """
    import builtins

    from neoswga.core.pool_plan_report import write_pool_plan

    real_import = builtins.__import__

    def without_matplotlib(name, *args, **kwargs):
        if name == "matplotlib" or name.startswith("matplotlib."):
            raise ImportError("No module named 'matplotlib'")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", without_matplotlib)

    output = tmp_path / "report"
    page = write_pool_plan(_saved_plan(), output).read_text()

    assert not (output / "pool_sizes.png").exists()
    assert "<img" not in page, "the page must not point at a figure it has not written"
    assert "matplotlib is not installed" in page
    # The outputs that carry the actual result are all present.
    assert (output / "pool_plan.json").is_file()
    assert (output / "pool_sizes.csv").is_file()
    assert (output / "target_90pct_oligos.fasta").is_file()
