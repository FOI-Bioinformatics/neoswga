"""`background_aggregate` decides the step-2 background gate (Phase 7b).

Phase 7a plumbed the key and asserted that nothing read it. This file is the
replacement for the step-2 half of that claim: the gate reads it, `worst-case`
removes candidates the pooled gate keeps, and `sum` is the pooled gate
unchanged.

**What makes the `sum` comparison real.** "Unchanged" is asserted against
figures recorded from the PRE-CHANGE tree -- commit 8e0ff18, the tip of
`phase7a-aggregate-key` -- not against a second run of the new code. The
synthetic world below was built there first and its `step2_df.csv` and
`filter_stats.json` written down; `BASELINE_*` are those numbers. The gate's
`sum` arithmetic is additionally checked against the pooled formula transcribed
from the code it replaced, so two independent oracles cover the default.

**The world.** One 14 kb target and two hosts, 40 kb and 2.0 Mb, a factor of
50 apart. A 60 bp stretch of the target is repeated three times in the SMALL
host. At `max_bg_freq` 5e-6 and k = 10 (where the 4**(10-k) scaling is 1) the
pooled ceiling is 10.2 sites over 2,040,000 bp and the small host's own ceiling
is 0.2 sites over 40,000 bp. A candidate with three small-host sites and a few
incidental large-host ones is therefore comfortably under the pooled ceiling
and far over the small host's. That is the measured shape of the thing: with
the three real hosts pooled, hg38 was 95.7 percent of the pooled length, so the
pooled ceiling sat within 5 percent of hg38's own while the numerator still
charged a candidate for its sites on the small hosts
(`docs/validation/2026-10-02-diversity-baseline.md`, Q4 and Q5).
"""

import hashlib
import json
import os
import random
import shutil

import pytest

from neoswga.core import background_gate
from neoswga.core.background_aggregate import SUM, WORST_CASE
from neoswga.core.exceptions import InvalidDesignRequest

# ---------------------------------------------------------------------------
# The synthetic world, and the baseline recorded from the pre-change tree
# ---------------------------------------------------------------------------

#: The columns `step2_df.csv` carried before 7b. The digest below is over these
#: and the row index only, so the per-host columns 7b adds cannot mask a change
#: to any figure that existed already.
BASELINE_COLUMNS = [
    "primer",
    "fg_count",
    "bg_count",
    "gini",
    "gini_bool",
    "ratio",
    "occupancy_ratio",
]

#: sha256 of `step2_df.csv` restricted to BASELINE_COLUMNS, as the default path
#: produced it on 8e0ff18 (jellyfish; the counts are exact, so the counter does
#: not enter).
BASELINE_TABLE_SHA = "000d2f457cfd5ea56fcc595000c1359fcf576f3369dbc075dd2c2d72f98de156"

#: `filter_stats.json` on 8e0ff18, stage counts only.
BASELINE_FUNNEL = {
    "total_kmers": 13741,
    "after_fg_frequency": 13741,
    "after_bg_frequency": 13708,
    "after_thermodynamic": 9636,
    "after_gini": 9636,
    "after_max_primer_cut": 9636,
    "final_candidates": 9636,
}

#: How many copies of the shared stretch the small host carries.
SHARED_COPIES = 3

pytestmark = pytest.mark.filterwarnings("ignore::DeprecationWarning")


def _write_fasta(path, name, sequence):
    path.write_text(
        f">{name}\n" + "\n".join(sequence[i : i + 70] for i in range(0, len(sequence), 70)) + "\n"
    )


def _build_world(root):
    rng = random.Random(70020)

    def randseq(length):
        return "".join(rng.choice("ACGT") for _ in range(length))

    target = randseq(14_000)
    shared = target[6_000:6_060]
    sequences = {
        "target": target,
        "small": randseq(20_000) + (shared * SHARED_COPIES) + randseq(19_820),
        "large": randseq(2_000_000),
    }
    paths = {}
    for name, sequence in sequences.items():
        path = root / f"{name}.fasta"
        _write_fasta(path, name, sequence)
        paths[name] = str(path)
    return paths, {shared[i : i + 10] for i in range(51)}


def _params(root, paths, **overrides):
    config = {
        "fg_genomes": [paths["target"]],
        "bg_genomes": [paths["small"], paths["large"]],
        "fg_prefixes": [str(root / "target")],
        "bg_prefixes": [str(root / "small"), str(root / "large")],
        "data_dir": str(root / "results"),
        "min_k": 10,
        "max_k": 10,
        "polymerase": "phi29",
        "reaction_temp": 30.0,
        "min_fg_freq": 1e-6,
        "max_bg_freq": 5e-6,
        "max_gini": 1.0,
        "min_gini_sites": 1,
        "max_primer": 10000,
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
    config.update(overrides)
    return config


def _run_filter(root, config):
    """One real `count-kmers` plus `filter`, with the globals put back.

    `get_params` assigns `parameter` module globals that outlive the call, and
    `background_aggregate` is one of them. Leaving `worst-case` behind made
    `resolve_design_request` read it in another test file, where a request with
    the key absent then resolved to `worst-case` and its hash moved. The
    product path does not have that problem -- there the global and the file
    agree, because the global came from the file -- but a test that leaves it
    set is a test that decides another file's result.
    """
    from neoswga import cli_unified
    from neoswga.core import parameter, pipeline

    params_file = root / "params.json"
    params_file.write_text(json.dumps(config, indent=2))
    previous = os.getcwd()
    previous_aggregate = getattr(parameter, "background_aggregate", None)
    os.chdir(root)
    try:
        pipeline.reset_pipeline_state()
        for command, handler in (
            ("count-kmers", cli_unified.run_step1),
            ("filter", cli_unified.run_step2),
        ):
            args = cli_unified.create_parser().parse_args([command, "-j", str(params_file)])
            handler(args)
    finally:
        os.chdir(previous)
        parameter.background_aggregate = previous_aggregate
        pipeline.reset_pipeline_state()
    return root / "results"


def _read_outputs(data_dir):
    import pandas as pd

    table = pd.read_csv(data_dir / "step2_df.csv", index_col=0)
    stats = json.loads((data_dir / "filter_stats.json").read_text())
    return table, stats


def _baseline_digest(table):
    payload = table[BASELINE_COLUMNS].to_csv(index=True, float_format="%.10g")
    return hashlib.sha256(payload.encode()).hexdigest()


@pytest.fixture(scope="module")
def filtered_worlds(tmp_path_factory):
    """Step 2 run three ways on the same world: key absent, `sum`, `worst-case`.

    One fixture, because `count-kmers` on the 2.0 Mb host is the cost and the
    three runs share nothing else. Each run gets its own directory, so no run
    reads another's tables or indexes.
    """
    if not shutil.which("jellyfish") and not shutil.which("kmc"):
        pytest.skip("no k-mer counter available (required for count-kmers)")

    root = tmp_path_factory.mktemp("bg_aggregate")
    worlds = {}
    for label, overrides in (
        ("absent", {}),
        ("sum", {"background_aggregate": "sum"}),
        ("worst-case", {"background_aggregate": "worst-case"}),
    ):
        directory = root / label.replace("-", "_")
        directory.mkdir()
        paths, shared = _build_world(directory)
        data_dir = _run_filter(directory, _params(directory, paths, **overrides))
        table, stats = _read_outputs(data_dir)
        worlds[label] = {
            "table": table,
            "stats": stats,
            "shared": shared,
            "data_dir": data_dir,
        }
    return worlds


# ---------------------------------------------------------------------------
# The default, against the pre-change tree
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("label", ["absent", "sum"])
def test_the_default_reproduces_the_pre_change_step2_table(filtered_worlds, label):
    """Every figure step 2 produced before 7b, unchanged, row order included.

    Compared against BASELINE_TABLE_SHA, recorded from 8e0ff18, rather than
    against the other new run: two new paths agreeing says nothing about
    whether either matches what shipped.
    """
    table = filtered_worlds[label]["table"]

    assert len(table) == BASELINE_FUNNEL["final_candidates"]
    assert _baseline_digest(table) == BASELINE_TABLE_SHA


@pytest.mark.parametrize("label", ["absent", "sum"])
def test_the_default_reproduces_the_pre_change_funnel(filtered_worlds, label):
    stats = filtered_worlds[label]["stats"]

    assert {key: stats[key] for key in BASELINE_FUNNEL} == BASELINE_FUNNEL


def test_an_absent_key_and_sum_are_one_design_on_a_real_run():
    """The request-level claim, which `test_background_aggregate_is_plumbed`
    makes on the hash, restated where the gate can disagree with it."""
    assert background_gate.pooled_verdict(5, 1_000_000, 1e-5) is True
    assert background_gate.host_panel(["a"], [100], None) == background_gate.host_panel(
        ["a"], [100], SUM
    )


def test_the_stage_is_recorded_in_filter_stats(filtered_worlds):
    """A funnel that does not say what the background was describes numbers it
    cannot account for."""
    assert filtered_worlds["absent"]["stats"]["background_aggregate"] == SUM
    assert filtered_worlds["sum"]["stats"]["background_aggregate"] == SUM
    assert filtered_worlds["worst-case"]["stats"]["background_aggregate"] == WORST_CASE


def test_every_key_the_filter_writes_is_read_by_the_report(filtered_worlds):
    """`FilteringStats` is built from `filter_stats.json`, so a key the filter
    writes must be a declared field AND pass the loader's filter, or it is a
    TypeError in `report` on one side and inert on the other. Both happened:
    adding `occupancy_ranking_mode` broke a funnel test, and the fix left it
    declared but dropped by a hand-maintained whitelist."""
    from neoswga.core.report.metrics import FILTERING_STATS_FIELDS, FilteringStats

    written = set(filtered_worlds["worst-case"]["stats"])

    assert written <= FILTERING_STATS_FIELDS, sorted(written - FILTERING_STATS_FIELDS)
    # Built, not merely accepted: the whitelist is derived from the fields, so a
    # declared field can no longer be dropped on the way in.
    stats = FilteringStats(**filtered_worlds["worst-case"]["stats"])
    assert stats.background_aggregate == WORST_CASE
    assert "background_aggregate" not in {label for label, _ in stats.as_funnel()}


def test_a_filter_stats_file_with_the_qa_stage_still_renders(tmp_path):
    """`after_qa` is the one key written to this file that is deliberately not
    a `FilteringStats` field: `pipeline_qa_integration` adds it and the
    rendered funnel does not show it. It must be dropped, not raise."""
    from neoswga.core.report.metrics import collect_pipeline_metrics

    (tmp_path / "step4_improved_df.csv").write_text("sequence\nACGTACGTACGT\n")
    (tmp_path / "filter_stats.json").write_text(
        json.dumps({"total_kmers": 100, "after_qa": 7, "final_candidates": 7})
    )

    metrics = collect_pipeline_metrics(str(tmp_path))

    assert metrics.filtering is not None
    assert metrics.filtering.final_candidates == 7
    assert metrics.filtering.background_aggregate == SUM


# ---------------------------------------------------------------------------
# `worst-case` acts, and on the case the measurement found
# ---------------------------------------------------------------------------


def test_worst_case_removes_more_at_the_background_gate(filtered_worlds):
    pooled = filtered_worlds["sum"]["stats"]
    worst = filtered_worlds["worst-case"]["stats"]

    assert worst["total_kmers"] == pooled["total_kmers"]
    assert worst["after_fg_frequency"] == pooled["after_fg_frequency"]
    assert worst["after_bg_frequency"] < pooled["after_bg_frequency"]
    assert len(filtered_worlds["worst-case"]["table"]) < len(filtered_worlds["sum"]["table"])


def test_a_candidate_dense_on_the_small_host_is_kept_by_sum_and_dropped_by_worst_case(
    filtered_worlds,
):
    """The case `docs/validation/2026-10-02-diversity-baseline.md` found.

    The stretch planted in the small host carries three sites there and a
    handful incidentally in the 2.0 Mb host. Its pooled frequency is below the
    pooled ceiling and its small-host frequency is fifteen times the small
    host's own ceiling, so pooling is what admits it.
    """
    pooled = filtered_worlds["sum"]
    shared = pooled["shared"]
    kept = set(pooled["table"]["primer"]) & shared
    dropped = set(filtered_worlds["worst-case"]["table"]["primer"]) & shared

    assert kept, "the planted stretch must survive the pooled gate, or there is nothing to show"
    assert dropped == set()

    # And the reason, from the columns the same run wrote: under the pooled
    # ceiling, over the small host's own.
    rows = pooled["table"][pooled["table"]["primer"].isin(kept)]
    pooled_ceiling = 5e-6 * (40_000 + 2_000_000)
    small_ceiling = 5e-6 * 40_000
    assert (rows["bg_count"] < pooled_ceiling).all()
    assert (rows["bg_count_small"] >= SHARED_COPIES).all()
    assert (rows["bg_count_small"] > small_ceiling).all()


def test_the_per_host_columns_appear_and_sum_to_the_pooled_count(filtered_worlds):
    table = filtered_worlds["sum"]["table"]

    assert "bg_count_small" in table.columns
    assert "bg_count_large" in table.columns
    assert (table["bg_count_small"] + table["bg_count_large"] == table["bg_count"]).all()


def test_the_sort_still_ends_on_the_primer_sequence(filtered_worlds):
    """So the counter's emission order cannot change a result.

    The digest test above is the strong form of this, since it pins the row
    order. This adds the structural statement: the per-host columns are not
    ranking keys, and `primer` is still the last one.
    """
    from neoswga.core.pipeline import _RANK_BY_OCCUPANCY, _RANK_BY_RATIO

    for spec in (_RANK_BY_RATIO, _RANK_BY_OCCUPANCY):
        assert spec["by"][-1] == "primer"
        assert spec["ascending"][-1] is True
        assert not [key for key in spec["by"] if key.startswith(background_gate.COLUMN_PREFIX)]


# ---------------------------------------------------------------------------
# The gate arithmetic, with the counts controlled
# ---------------------------------------------------------------------------

SMALL = "/tmp/small"
LARGE = "/tmp/large"


def _fake_counts(tables):
    """Stand in for `filter._get_rate_for_one_file`, reading a prefix table.

    `tables` maps a prefix to {primer: count}. A primer the table does not
    mention is reported as 0, which is the contract the real function keeps
    (every requested primer appears, because an absent key reads as unknown). A
    prefix mapped to None answers nothing at all, which is how an absent count
    is produced below.
    """

    def counter(task):
        primer_list, prefix, _k, _genome = task
        table = tables.get(prefix, {})
        if table is None:
            return {}
        return {primer: int(table.get(primer, 0)) for primer in primer_list}

    return counter


def _run_gate(monkeypatch, primers, fg, tables, *, aggregate, lengths, max_bg_freq=5e-6):
    from neoswga.core import filter as filter_module
    from neoswga.core import parameter

    monkeypatch.setattr(parameter, "min_fg_freq", 1e-9, raising=False)
    monkeypatch.setattr(parameter, "max_bg_freq", max_bg_freq, raising=False)
    monkeypatch.setattr(parameter, "use_bloom_filter", False, raising=False)
    monkeypatch.setattr(parameter, "bloom_filter_path", None, raising=False)
    monkeypatch.setattr(parameter, "background_aggregate", aggregate, raising=False)
    monkeypatch.setattr(
        filter_module, "_get_rate_for_one_file", _fake_counts({"/tmp/fg": fg, **tables})
    )
    return filter_module.get_all_rates(
        primers,
        fg_prefixes=["/tmp/fg"],
        bg_prefixes=list(tables),
        fg_total_length=100_000,
        bg_total_length=sum(lengths),
        bg_seq_lengths=list(lengths),
    )


def test_sum_reproduces_the_pooled_formula_it_replaced(monkeypatch):
    """The oracle for the default: the rule as the code read before 7b.

        if bg_count is None or bg_total_length <= 0:
            bg_bool = bg_count is None
        else:
            bg_bool = bg_count / bg_total_length < scaled_max_bg

    At k = 10 the 4**(10-k) scaling is 1, so `scaled_max_bg` is `max_bg_freq`.
    """
    primers = ["AAAACCCCGG", "AAAACCCCGT", "AAAACCCCTT", "AAAAGGGGCC"]
    fg = dict.fromkeys(primers, 10)
    tables = {SMALL: {primers[0]: 3, primers[1]: 0}, LARGE: {primers[0]: 2, primers[2]: 40}}
    lengths = (40_000, 2_000_000)

    frame = _run_gate(monkeypatch, primers, fg, tables, aggregate=SUM, lengths=lengths).set_index(
        "primer"
    )

    for primer in primers:
        pooled = int(frame.loc[primer, "bg_count"])
        expected = pooled / sum(lengths) < 5e-6
        assert bool(frame.loc[primer, "bg_bool"]) is expected, primer
    # And the pooled count is still the sum over the hosts.
    assert int(frame.loc[primers[0], "bg_count"]) == 5
    assert bool(frame.loc[primers[2], "bg_bool"]) is False  # 40 / 2.04 Mb = 1.96e-5


def test_worst_case_judges_each_host_on_its_own_length(monkeypatch):
    primers = ["AAAACCCCGG", "AAAACCCCGT"]
    fg = dict.fromkeys(primers, 10)
    # One site on the small host is 2.5e-5, five times its own 5e-6 ceiling,
    # while 1/2,040,000 pooled is 4.9e-7 and passes.
    tables = {SMALL: {primers[0]: 1}, LARGE: {}}
    lengths = (40_000, 2_000_000)

    pooled = _run_gate(monkeypatch, primers, fg, tables, aggregate=SUM, lengths=lengths).set_index(
        "primer"
    )
    worst = _run_gate(
        monkeypatch, primers, fg, tables, aggregate=WORST_CASE, lengths=lengths
    ).set_index("primer")

    assert bool(pooled.loc[primers[0], "bg_bool"]) is True
    assert bool(worst.loc[primers[0], "bg_bool"]) is False
    # A candidate on neither host passes under both.
    assert bool(pooled.loc[primers[1], "bg_bool"]) is True
    assert bool(worst.loc[primers[1], "bg_bool"]) is True


def test_worst_case_with_one_host_equals_sum_with_one_host(monkeypatch):
    """Exactly, frame for frame. With one host the pooled length IS that host's
    own length, so the two readings are the same arithmetic; anything that
    differed would be a difference the key invented."""
    import pandas.testing as pdt

    primers = ["AAAACCCCGG", "AAAACCCCGT", "AAAACCCCTT"]
    fg = dict.fromkeys(primers, 10)
    tables = {SMALL: {primers[0]: 1, primers[1]: 0, primers[2]: 7}}
    lengths = (400_000,)

    pooled = _run_gate(monkeypatch, primers, fg, tables, aggregate=SUM, lengths=lengths)
    worst = _run_gate(monkeypatch, primers, fg, tables, aggregate=WORST_CASE, lengths=lengths)

    pdt.assert_frame_equal(pooled, worst)
    assert list(pooled.columns) == ["primer", "fg_count", "bg_count", "fg_bool", "bg_bool"]
    # Both verdicts are exercised, or the equality is vacuous.
    assert set(pooled["bg_bool"]) == {True, False}


def test_an_absent_per_host_count_blocks_under_worst_case_and_passes_under_sum(monkeypatch):
    """The deliberate difference. Today an absent POOLED count passes, because
    a k-mer missing from a counter's table may still have sites the string
    search finds. The worst case over a set with an unmeasured member is not a
    measured worst, so it cannot pass on the same argument."""
    primers = ["AAAACCCCGG"]
    fg = dict.fromkeys(primers, 10)
    # The large host answers nothing at all, so this primer has no count there.
    tables = {SMALL: {primers[0]: 0}, LARGE: None}
    lengths = (40_000, 2_000_000)

    pooled = _run_gate(monkeypatch, primers, fg, tables, aggregate=SUM, lengths=lengths)
    worst = _run_gate(monkeypatch, primers, fg, tables, aggregate=WORST_CASE, lengths=lengths)

    assert bool(pooled.loc[0, "bg_bool"]) is True
    assert bool(worst.loc[0, "bg_bool"]) is False
    # And the column says which host had no answer, as None rather than 0.
    assert pooled.loc[0, "bg_count_small"] == 0
    assert pooled.loc[0, "bg_count_large"] is None


def test_an_absent_pooled_count_still_passes_under_sum(monkeypatch):
    """Today's rule, unchanged, on the path that produced it: with no
    background prefix at all every count is absent and the gate passes."""
    primers = ["AAAACCCCGG"]
    frame = _run_gate(
        monkeypatch, primers, dict.fromkeys(primers, 10), {}, aggregate=SUM, lengths=()
    )

    assert frame.loc[0, "bg_count"] is None
    assert bool(frame.loc[0, "bg_bool"]) is True


def test_no_host_configured_is_a_configuration_not_an_unmeasured_host(monkeypatch):
    """`worst-case` with no background listed passes, like `sum`. There is no
    host to judge, which is different from a host that could not be measured."""
    primers = ["AAAACCCCGG"]
    frame = _run_gate(
        monkeypatch, primers, dict.fromkeys(primers, 10), {}, aggregate=WORST_CASE, lengths=()
    )

    assert bool(frame.loc[0, "bg_bool"]) is True
    assert background_gate.host_panel([], [], WORST_CASE) == ()


# ---------------------------------------------------------------------------
# What refuses rather than guesses
# ---------------------------------------------------------------------------


def test_unpairable_prefixes_and_lengths_refuse_under_worst_case():
    """A count divided by another host's length is a figure about nothing. The
    pairing is made where both lists are in hand and refused when it cannot
    be."""
    with pytest.raises(InvalidDesignRequest) as caught:
        background_gate.host_panel(["a", "b"], [100], WORST_CASE)

    message = str(caught.value)
    assert "background_aggregate" in message
    assert "bg_seq_lengths" in message


def test_a_host_with_no_positive_length_refuses_under_worst_case():
    for lengths in ([100, 0], [100, None]):
        with pytest.raises(InvalidDesignRequest):
            background_gate.host_panel(["a", "b"], lengths, WORST_CASE)


def test_the_same_input_leaves_sum_alone():
    """`sum` needs only the summed length, so an unpairable list is not its
    problem and must not become one. It gets no per-host columns, which is the
    only consequence."""
    assert background_gate.host_panel(["a", "b"], [100], SUM) == ()
    assert background_gate.host_panel(["a", "b"], [], SUM) == ()


def test_a_bloom_filter_background_refuses_under_worst_case(monkeypatch):
    """The filter answers for the pooled background, and `parameter`
    deliberately leaves `bg_prefixes` empty for it, so there is no per-host
    count to take a worst case over. Passing everything there would be the
    'unknown reported as success' shape the Bloom guard in
    `_resolve_background_source` already exists to stop."""
    from neoswga.core import filter as filter_module
    from neoswga.core import parameter

    monkeypatch.setattr(parameter, "use_bloom_filter", True, raising=False)
    monkeypatch.setattr(parameter, "bloom_filter_path", "/tmp/does-not-exist.pkl", raising=False)
    monkeypatch.setattr(parameter, "background_aggregate", WORST_CASE, raising=False)
    monkeypatch.setattr(parameter, "min_fg_freq", 1e-9, raising=False)
    monkeypatch.setattr(parameter, "max_bg_freq", 5e-6, raising=False)
    monkeypatch.setattr(
        filter_module, "_get_rate_for_one_file", _fake_counts({"/tmp/fg": {"AAAACCCCGG": 10}})
    )

    with pytest.raises(InvalidDesignRequest) as caught:
        filter_module.get_all_rates(
            ["AAAACCCCGG"],
            fg_prefixes=["/tmp/fg"],
            bg_prefixes=[],
            fg_total_length=100_000,
            bg_total_length=0,
        )

    assert "Bloom" in str(caught.value)


def test_an_explicit_argument_beats_the_parameter_global(monkeypatch):
    """The 7a seam, exercised through the gate: a caller that resolved the
    run's configuration once must not have a mutable global re-read under it."""
    from neoswga.core import filter as filter_module
    from neoswga.core import parameter

    primers = ["AAAACCCCGG"]
    monkeypatch.setattr(parameter, "min_fg_freq", 1e-9, raising=False)
    monkeypatch.setattr(parameter, "max_bg_freq", 5e-6, raising=False)
    monkeypatch.setattr(parameter, "use_bloom_filter", False, raising=False)
    monkeypatch.setattr(parameter, "bloom_filter_path", None, raising=False)
    monkeypatch.setattr(parameter, "background_aggregate", WORST_CASE, raising=False)
    monkeypatch.setattr(
        filter_module,
        "_get_rate_for_one_file",
        _fake_counts({"/tmp/fg": {primers[0]: 10}, SMALL: {primers[0]: 1}}),
    )

    frame = filter_module.get_all_rates(
        primers,
        fg_prefixes=["/tmp/fg"],
        bg_prefixes=[SMALL],
        fg_total_length=100_000,
        bg_total_length=40_000,
        bg_seq_lengths=[40_000],
        background_aggregate=SUM,
    )

    assert frame.attrs["background_aggregate"] == SUM
    # 1 / 40,000 = 2.5e-5, over the 5e-6 ceiling, so this is the arithmetic of
    # one host either way; what is asserted is which value was applied.
    assert bool(frame.loc[0, "bg_bool"]) is False


def test_an_unrecognised_value_refuses_at_the_gate(monkeypatch):
    from neoswga.core import filter as filter_module
    from neoswga.core import parameter

    monkeypatch.setattr(parameter, "background_aggregate", "mean", raising=False)

    with pytest.raises(InvalidDesignRequest):
        filter_module.get_all_rates(
            ["AAAACCCCGG"],
            fg_prefixes=["/tmp/fg"],
            bg_prefixes=[],
            fg_total_length=100_000,
            bg_total_length=0,
        )


# ---------------------------------------------------------------------------
# Bookkeeping that is easy to get subtly wrong
# ---------------------------------------------------------------------------


def test_the_per_prefix_counts_sum_to_what_the_pooled_call_returns(monkeypatch):
    """`get_rates_for_one_species` is `rates_by_prefix` followed by
    `pool_counts`, so the two cannot report different counts for one primer.
    Pinned because the pooled function is what every existing caller uses."""
    from neoswga.core import filter as filter_module

    primers = ["AAAACCCCGG", "AAAACCCCGT"]
    tables = {SMALL: {primers[0]: 3}, LARGE: {primers[0]: 2, primers[1]: 9}}
    monkeypatch.setattr(filter_module, "_get_rate_for_one_file", _fake_counts(tables))

    by_prefix = filter_module.rates_by_prefix(primers, [SMALL, LARGE])
    pooled = filter_module.get_rates_for_one_species(primers, [SMALL, LARGE])

    assert set(by_prefix) == {SMALL, LARGE}
    assert pooled == background_gate.pool_counts(by_prefix)
    assert pooled[primers[0]] == 5
    assert pooled[primers[1]] == 9


def test_a_prefix_that_answered_nothing_is_still_a_key():
    """ "This host answered nothing" and "this host was not asked" have to stay
    apart, because the second is what blocks under `worst-case`."""
    panel = background_gate.host_panel([SMALL, LARGE], [40_000, 2_000_000], WORST_CASE)

    asked = {SMALL: {}, LARGE: {"AAAACCCCGG": 0}}
    not_asked = {LARGE: {"AAAACCCCGG": 0}}

    assert background_gate.worst_case_verdict(asked, "AAAACCCCGG", panel, 5e-6) is False
    assert background_gate.worst_case_verdict(not_asked, "AAAACCCCGG", panel, 5e-6) is False
    both = {SMALL: {"AAAACCCCGG": 0}, LARGE: {"AAAACCCCGG": 0}}
    assert background_gate.worst_case_verdict(both, "AAAACCCCGG", panel, 5e-6) is True


def test_column_names_are_deterministic_and_collision_free():
    names = background_gate.host_column_names(["/a/human", "/b/human", "/c/fly"])

    assert names["/a/human"] == "bg_count_human"
    assert names["/b/human"] == "bg_count_human__1"
    assert names["/c/fly"] == "bg_count_fly"
    assert len(set(names.values())) == 3


def test_one_host_gets_no_per_host_column():
    """It would repeat `bg_count` value for value, and withholding it is what
    makes a single-host table byte-identical under either reading."""
    one = background_gate.host_panel([SMALL], [40_000], WORST_CASE)
    two = background_gate.host_panel([SMALL, LARGE], [40_000, 2_000_000], WORST_CASE)

    assert background_gate.columns_are_written(one) is False
    assert background_gate.columns_are_written(two) is True
