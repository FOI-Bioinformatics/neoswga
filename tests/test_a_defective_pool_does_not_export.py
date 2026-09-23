"""Two defects were recorded at error level and exported anyway.

`BLOCKING_VALIDATOR_CODES` held two members. `validate_result` emits four codes
and two of them -- `duplicate_primers` and `blacklist_primer_in_set` -- were not
among them, so a panel holding the same oligo twice, or holding an oligo the
user blacklisted, printed "Primers ready for ordering!".

Measured before the change:

    code                                   ok       export blocked
    delivered_pool_exceeds_max_dimer_bp    False    True
    panel_limit_not_met                    False    True
    duplicate_primers                      False    False
    blacklist_primer_in_set                False    False
    set_size_mismatch                      False    False

`ok` was False in every row, which is why this survived: the flag says what a
reader expects and no command consults it. The gate is the code set.

Both additions pass the membership test that constant's own comment states --
findings that are defects IN THE POOL -- and both are things the user can act
on. A panel of twelve holding one oligo twice is eleven oligos: twelve tubes
ordered, eleven distinct sequences delivered. A blacklisted oligo in the panel
is the user's own instruction violated.

Reachability differs between them and the difference is worth recording rather
than smoothing over.

`blacklist_primer_in_set` fires on the real path.
`_collect_forbidden_primers` runs inside `run_optimization` whenever
`bl_prefixes` and `bl_seq_lengths` are configured, `examples/multi_genome_
blacklist/params.json` configures them, and
`tests/test_core_optimize_blacklist_guard.py` drives `optimize_step4` through
it. That is a real hole being closed.

`duplicate_primers` is not reachable from any path checked. Every point that
assembles a panel from two sources deduplicates first:
`hybrid_optimizer.py:720` filters the fixed oligos out of the new ones,
`dominating_set_optimizer.py:482` and `:1149` and `panel_refinement.py:70` use
`dict.fromkeys`, and the ensemble union pool at `unified_optimizer.py:340` does
too. So this one is insurance rather than a repair, and it costs a string in a
frozenset. Not reachable is not the same as impossible -- checking the assembly
points is checking candidates, not candidate subsets.

`set_size_mismatch` stays OUT, and that is the decision this file also pins. A
panel shorter than requested is a documented normal outcome: `num_primers` is a
request rather than a guarantee, and selection stops rather than admitting a
pair above `max_dimer_bp`. Blocking on it would refuse three of four runs on
the shipped plasmid example, which is measured in
docs/validation/2026-09-23-assessment-vs-validator.md.
"""

import json

import pytest

from neoswga.core.design_result import (
    BLOCKING_VALIDATOR_CODES,
    VALIDATION_FILENAME,
    blocking_validator_findings,
    panel_validation_is_ok,
)
from neoswga.core.export import export_is_blocked

BLOCKS = ("delivered_pool_exceeds_max_dimer_bp", "panel_limit_not_met")
NOW_BLOCKS = ("duplicate_primers", "blacklist_primer_in_set")
MUST_NOT_BLOCK = ("set_size_mismatch", "coverage_saturated_on_small_genome")


def _record(directory, code, level="error"):
    issues = [{"level": level, "code": code, "detail": f"{code} was recorded"}]
    payload = {
        "optimizer": "test",
        "num_primers": 2,
        "ok": panel_validation_is_ok(issues),
        "issues": issues,
    }
    (directory / VALIDATION_FILENAME).write_text(json.dumps(payload))
    return directory


@pytest.mark.parametrize("code", NOW_BLOCKS)
def test_a_pool_defect_recorded_at_error_level_blocks_the_export(tmp_path, code):
    """The defect. Both were recorded, both exported."""
    _record(tmp_path, code)

    assert code in BLOCKING_VALIDATOR_CODES
    assert export_is_blocked(str(tmp_path)), f"{code} did not block the export"


@pytest.mark.parametrize("code", BLOCKS)
def test_the_codes_that_already_blocked_still_do(tmp_path, code):
    _record(tmp_path, code)

    assert export_is_blocked(str(tmp_path)), f"{code} stopped blocking"


@pytest.mark.parametrize("code", MUST_NOT_BLOCK)
def test_a_normal_outcome_is_not_a_defect(tmp_path, code):
    """Guard the guard, and pin the decision.

    A short panel and a saturation warning are both things the user cannot fix
    and should not be refused for. Widening the set until it catches these
    would teach people to pass --allow-unqualified by reflex, which is the
    reason `coverage_saturated_on_small_genome` was kept out in the first
    place.
    """
    _record(tmp_path, code)

    assert code not in BLOCKING_VALIDATOR_CODES
    assert export_is_blocked(str(tmp_path)) is None, f"{code} now blocks the export"


def test_a_clean_pool_still_exports(tmp_path):
    payload = {"optimizer": "test", "num_primers": 2, "ok": True, "issues": []}
    (tmp_path / VALIDATION_FILENAME).write_text(json.dumps(payload))

    assert export_is_blocked(str(tmp_path)) is None


def test_the_refusal_names_the_cause(tmp_path):
    """A refusal a user cannot act on is barely better than none, and this is
    what the single-bucket alternative would have cost: routing the panel
    assessment's prose violations through one code would have said only that
    something was wrong."""
    _record(tmp_path, "duplicate_primers")

    findings = blocking_validator_findings(tmp_path)

    assert findings
    assert "duplicate_primers" in findings[0]


def test_interpret_and_export_agree_about_the_same_pool(tmp_path):
    """Two commands tell a user a pool is ready. They read one code set, so
    adding a member must move both at once."""
    from neoswga.core.results_interpreter import ResultsInterpreter

    _record(tmp_path, "blacklist_primer_in_set")

    interpreter = ResultsInterpreter.__new__(ResultsInterpreter)
    interpreter.validation_file = tmp_path / VALIDATION_FILENAME

    assert interpreter._blocks_synthesis() is True
    assert export_is_blocked(str(tmp_path))
