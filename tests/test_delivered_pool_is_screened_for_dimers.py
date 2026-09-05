"""A pool that breaks the configured dimer threshold must say so.

Every stage between the optimizer and the order form used to call a set clean
that was not. On the plasmid example with `max_dimer_bp: 3`, the network
optimizer returned six oligos whose worst pair shared a 10 bp duplex, and:

  - `neoswga interpret` printed "Overall rating: EXCELLENT", "Ready for
    synthesis", and "1. Order primers for experimental validation";
  - `neoswga export --format fasta` printed "Primers ready for ordering!";
  - `neoswga report` graded Dimer Risk "Low / Excellent".

The capability to catch it already existed and worked: `neoswga analyze-dimers`
scored that exact set FAIL at severity 1.00. It is a separate command that takes
sequences as arguments, and nothing in the pipeline called it.

That a soft multiplicative penalty cannot exclude a strongly scoring pair is a
documented design limitation, not a bug -- only the clique optimizer constrains
the set structurally. This file is about the other half: whatever the optimizer
returns, the run has to measure the delivered pool against the threshold the user
set, record the answer, and not describe a set as ready to order when it is not.
"""

import pytest

from neoswga.core.dimer import max_complementary_run
from neoswga.core.thermodynamics import reverse_complement

# A pair that is close to a perfect duplex: 10 of 12 bases complementary.
DIMERISING = ("TGTACTGCCAAG", reverse_complement("TGTACTGCCAAG"))
CLEAN = ["ACCACAGATAGC", "GTTGTAGATGGA", "ATCAGCAGACCA"]


def test_the_fixture_pair_really_dimerises():
    """Guard the guard: a pair that does not dimerise would make this file
    pass without testing anything."""
    assert max_complementary_run(*DIMERISING) >= 10


# ---------------------------------------------------------------------------
# The measurement
# ---------------------------------------------------------------------------


def test_the_worst_pair_is_found_and_named():
    from neoswga.core.unified_optimizer import worst_heterodimer

    primers = CLEAN + list(DIMERISING)
    length, pair = worst_heterodimer(primers)

    assert length >= 10, f"the worst pair measured {length} bp, expected at least 10"
    assert set(pair) == set(
        DIMERISING
    ), f"the worst pair was reported as {pair}, but the planted duplex is {DIMERISING}"


def test_a_clean_pool_reports_its_worst_pair_too():
    """The measurement is unconditional; only the verdict depends on the threshold."""
    from neoswga.core.unified_optimizer import worst_heterodimer

    length, pair = worst_heterodimer(CLEAN)
    assert length >= 0
    assert pair is None or len(pair) == 2


@pytest.mark.parametrize("n_primers", [0, 1])
def test_a_pool_too_small_to_have_a_pair_is_not_an_error(n_primers):
    from neoswga.core.unified_optimizer import worst_heterodimer

    length, pair = worst_heterodimer(CLEAN[:n_primers])
    assert length == 0
    assert pair is None


# ---------------------------------------------------------------------------
# The verdict
# ---------------------------------------------------------------------------


def test_a_pool_over_the_threshold_raises_a_validation_issue():
    from neoswga.core.unified_optimizer import dimer_validation_issue

    issue = dimer_validation_issue(CLEAN + list(DIMERISING), max_dimer_bp=3)

    assert issue is not None, "a 10 bp duplex against a threshold of 3 raised nothing"
    assert issue["code"] == "delivered_pool_exceeds_max_dimer_bp"
    assert issue["level"] == "warning"
    # The detail has to name the pair, because "your pool has a dimer" is not
    # actionable and "break this pair" is.
    for primer in DIMERISING:
        assert primer in issue["detail"], f"{primer} is not named in: {issue['detail']}"
    assert "3" in issue["detail"], "the configured threshold is not stated"


def test_a_pool_within_the_threshold_raises_nothing():
    from neoswga.core.unified_optimizer import dimer_validation_issue

    assert dimer_validation_issue(CLEAN, max_dimer_bp=3) is None


def test_the_threshold_is_inclusive_allowed():
    """Matches `is_dimer_fast`, which triggers at `max_dimer_bp + 1`.

    A pool whose worst pair is exactly the configured maximum is compliant, so
    the clique optimizer returning a set at exactly `max_dimer_bp` is the
    intended output rather than a near miss.
    """
    from neoswga.core.unified_optimizer import dimer_validation_issue, worst_heterodimer

    length, _pair = worst_heterodimer(CLEAN)
    assert dimer_validation_issue(CLEAN, max_dimer_bp=length) is None
    if length > 0:
        assert dimer_validation_issue(CLEAN, max_dimer_bp=length - 1) is not None


# ---------------------------------------------------------------------------
# What the commands that hand over oligos do with the finding
# ---------------------------------------------------------------------------


def _write_validation(tmp_path, codes):
    import json

    (tmp_path / "step4_improved_df_validation.json").write_text(
        json.dumps(
            {
                "optimizer": "network",
                "ok": True,
                "issues": [
                    {"level": "warning", "code": c, "detail": f"{c} happened"} for c in codes
                ],
            }
        )
    )
    return tmp_path


def test_a_dimer_violation_blocks_the_order_readiness_line(tmp_path):
    from neoswga.cli.report import _blocking_validator_findings

    _write_validation(tmp_path, ["delivered_pool_exceeds_max_dimer_bp"])
    assert _blocking_validator_findings(str(tmp_path))


def test_a_saturation_warning_does_not_block(tmp_path):
    """It says a metric is untrustworthy on a small target, which is inherent
    to plasmid-scale design. Blocking on it would refuse every such design and
    teach people to ignore the line."""
    from neoswga.cli.report import _blocking_validator_findings

    _write_validation(tmp_path, ["coverage_saturated_on_small_genome"])
    assert not _blocking_validator_findings(str(tmp_path))


def test_a_missing_validation_file_never_blocks(tmp_path):
    """An export must not fail closed because a validator did not run."""
    from neoswga.cli.report import _blocking_validator_findings

    assert not _blocking_validator_findings(str(tmp_path))


@pytest.mark.parametrize(
    "code,should_block",
    [
        ("delivered_pool_exceeds_max_dimer_bp", True),
        ("coverage_saturated_on_small_genome", False),
        ("reaction_conditions_init_failed", False),
    ],
)
def test_the_interpreter_and_the_exporter_agree_on_what_blocks(tmp_path, code, should_block):
    """Two code paths, one rule, checked by behaviour rather than by reading
    each other's source.

    `export` and `interpret` are the two commands that tell a user a pool is
    ready. They already drifted apart once over which column holds the dimer
    score, and each of them passing its own test while disagreeing with the
    other is exactly how that survived.
    """
    from neoswga.cli.report import _blocking_validator_findings
    from neoswga.core.results_interpreter import ResultsInterpreter

    _write_validation(tmp_path, [code])

    exporter_blocks = bool(_blocking_validator_findings(str(tmp_path)))
    interpreter = ResultsInterpreter.__new__(ResultsInterpreter)
    interpreter.validation_file = tmp_path / "step4_improved_df_validation.json"
    interpreter_blocks = interpreter._blocks_synthesis()

    assert exporter_blocks == interpreter_blocks == should_block, (
        f"{code}: export blocks={exporter_blocks}, interpret blocks="
        f"{interpreter_blocks}, expected {should_block}"
    )
