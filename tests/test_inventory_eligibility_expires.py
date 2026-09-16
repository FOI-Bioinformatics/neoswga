"""A verdict must not outlive the rules that produced it.

Phase 3 of the plan for `docs/validation/pipeline_audit_2026-09-16/`, audit
finding F3.

`record_stage2_inventory` upserts the current survivors as passing and never
touches a row from an earlier run. There is no delete anywhere in the package
and no run generation. The assessment key is
`(sequence, condition_id, policy_version)`, where `condition_id` is a chemistry
fingerprint and `policy_version` was a hand-bumped constant. The resolved QC
thresholds -- GC window, Tm window, frequency limits, Gini settings, exclusion
configuration -- appeared in none of it.

So tightening a rule and re-running into the same directory left every
previously admitted candidate eligible, under a verdict reached by rules that no
longer exist. That was merely untidy while nothing read the inventory. It stops
being untidy as soon as the inventory becomes the search universe, which is what
Phase 4 does, so this lands first.

Two mechanisms, answering two different questions.

The **policy digest** answers "under which rules was this decided". It goes in
`policy_version`, so a verdict reached under one set of thresholds is never
returned to a caller asking about another.

The **generation** answers "was this decided by the current run". Each recording
increments it, and eligibility reads only the newest, so a candidate absent from
the latest run stops being eligible even when the rules did not move at all.
The digest alone cannot do this: two runs with identical thresholds and
different inputs share a digest.
"""

import pathlib

import pandas as pd
import pytest

from neoswga.core.candidate_inventory import (
    STAGE2_INVENTORY_NAME,
    CandidateInventory,
    qc_policy_fingerprint,
    record_stage2_inventory,
)

A, C, G, T = "AAAAAAAAAAAA", "CCCCCCCCCCCC", "GGGGGGGGGGGG", "TTTTTTTTTTTT"
COND = "tm-2026-09-14:abc123"

STRICT = {"max_gini": 0.6, "min_tm": 25.0, "max_tm": 55.0, "min_fg_freq": 1e-5}
LOOSE = {"max_gini": 0.7, "min_tm": 25.0, "max_tm": 55.0, "min_fg_freq": 1e-5}


def _frame(sequences):
    return pd.DataFrame({"primer": list(sequences)})


def _record(path, survivors, policy):
    frame = _frame(survivors)
    return record_stage2_inventory(
        path,
        condition_id=COND,
        cleared_hard_gates=frame,
        after_gini=frame,
        shortlisted=frame,
        indexed=list(survivors),
        qc_policy=policy,
    )


def _eligible(data_dir, policy):
    """`record_stage2_inventory` returns the file; the tests hold the directory."""
    with CandidateInventory(pathlib.Path(data_dir) / STAGE2_INVENTORY_NAME) as inventory:
        return sorted(inventory.iter_eligible(COND, [12], policy_version=policy))


# -- the digest ------------------------------------------------------------


def test_the_same_thresholds_give_the_same_digest():
    assert qc_policy_fingerprint(STRICT) == qc_policy_fingerprint(dict(STRICT))


def test_a_changed_threshold_changes_the_digest():
    assert qc_policy_fingerprint(STRICT) != qc_policy_fingerprint(LOOSE)


def test_key_order_does_not_change_the_digest():
    """Otherwise a dict built in a different order would look like a new policy."""
    reordered = {key: STRICT[key] for key in reversed(list(STRICT))}

    assert qc_policy_fingerprint(reordered) == qc_policy_fingerprint(STRICT)


def test_the_digest_names_the_policy_version_it_extends():
    """A reader should be able to tell which generation of rules this is."""
    assert qc_policy_fingerprint(STRICT).startswith("qc-2026-09-15:")


def test_an_absent_threshold_is_not_the_same_as_a_set_one():
    """A missing key must not silently share a digest with a present one."""
    without = {k: v for k, v in STRICT.items() if k != "max_gini"}

    assert qc_policy_fingerprint(without) != qc_policy_fingerprint(STRICT)


# -- the generation --------------------------------------------------------


def test_a_stricter_rerun_drops_the_candidates_it_no_longer_admits(tmp_path):
    """The audit's reproduction.

    Two candidates admitted, then a stricter run admitting one. Both used to
    come back eligible.
    """
    _record(tmp_path, [A, C], LOOSE)
    strict_policy = qc_policy_fingerprint(STRICT)
    _record(tmp_path, [A], STRICT)

    assert _eligible(tmp_path, strict_policy) == [A]


def test_a_rerun_under_the_same_rules_also_supersedes(tmp_path):
    """The case a policy digest alone cannot catch.

    Same thresholds, different inputs -- a corrected reference, a different
    k-mer table. The digest is identical, so only the generation separates them.
    """
    policy = qc_policy_fingerprint(LOOSE)
    _record(tmp_path, [A, C, G], LOOSE)
    _record(tmp_path, [A, G], LOOSE)

    assert _eligible(tmp_path, policy) == sorted([A, G])


def test_the_earlier_verdict_is_kept_rather_than_deleted(tmp_path):
    """History is preserved; only eligibility moves.

    A dropped candidate's earlier assessment is still on record, because "this
    was admitted once, under these rules" is worth knowing when a design is
    being explained.
    """
    _record(tmp_path, [A, C], LOOSE)
    _record(tmp_path, [A], LOOSE)

    with CandidateInventory(tmp_path / STAGE2_INVENTORY_NAME) as inventory:
        assert inventory.has_length(12)
        rows = inventory.assessment_history(C, COND)

    assert rows, "the superseded verdict was deleted rather than retired"


def test_a_different_policy_does_not_see_the_other_policy_s_verdicts(tmp_path):
    """Two policies, two answers, from one database."""
    loose_policy = qc_policy_fingerprint(LOOSE)
    _record(tmp_path, [A, C, G], LOOSE)
    _record(tmp_path, [A], STRICT)

    assert _eligible(tmp_path, loose_policy) == sorted([A, C, G])
    assert _eligible(tmp_path, qc_policy_fingerprint(STRICT)) == [A]


def test_recording_nothing_does_not_wipe_the_previous_run(tmp_path):
    """A failed or empty run must not silently empty the inventory."""
    policy = qc_policy_fingerprint(LOOSE)
    _record(tmp_path, [A, C], LOOSE)

    with pytest.raises(ValueError, match="no candidates"):
        _record(tmp_path, [], LOOSE)

    assert _eligible(tmp_path, policy) == sorted([A, C])


def test_the_default_policy_version_still_works_for_a_caller_without_one(tmp_path):
    """Existing callers that name no policy keep a single coherent view."""
    frame = _frame([A, C])
    record_stage2_inventory(
        tmp_path,
        condition_id=COND,
        cleared_hard_gates=frame,
        after_gini=frame,
        shortlisted=frame,
        indexed=[A, C],
    )

    with CandidateInventory(tmp_path / STAGE2_INVENTORY_NAME) as inventory:
        assert sorted(inventory.iter_eligible(COND, [12])) == sorted([A, C])


# -- the policy has to be discoverable -------------------------------------


def test_the_current_policy_can_be_read_back(tmp_path):
    """A caller cannot reconstruct the digest, so it must be findable.

    It hashes thresholds resolved at run time, after the GC-adaptive strategy
    has had its say, so a reader asking under the bare constant would get an
    empty eligible set. That reads as "no candidates qualify" rather than "you
    asked the wrong question", which is the silent-zero shape of Known Issues
    5, 6 and 13.
    """
    _record(tmp_path, [A, C], STRICT)

    with CandidateInventory(tmp_path / STAGE2_INVENTORY_NAME) as inventory:
        assert inventory.current_policy(COND) == qc_policy_fingerprint(STRICT)


def test_the_current_policy_follows_the_newest_recording(tmp_path):
    _record(tmp_path, [A, C], LOOSE)
    _record(tmp_path, [A], STRICT)

    with CandidateInventory(tmp_path / STAGE2_INVENTORY_NAME) as inventory:
        assert inventory.current_policy(COND) == qc_policy_fingerprint(STRICT)


def test_a_stage_three_policy_is_not_mistaken_for_an_admission_rule(tmp_path):
    """Stage 3 records what an efficacy filter carried forward, not who qualifies."""
    from neoswga.core.candidate_inventory import record_stage3_policy

    path = _record(tmp_path, [A, C], STRICT)
    record_stage3_policy(path, COND, carried=[A], rejected={}, policy="carry_forward")

    with CandidateInventory(path) as inventory:
        assert inventory.current_policy(COND) == qc_policy_fingerprint(STRICT)


def test_the_provider_finds_the_candidates_the_pipeline_recorded(tmp_path):
    """The seam this closes.

    The provider used the bare constant. Once the filter started recording
    under a digest, that combination returned nothing at all -- and the
    provider's whole job is to say which candidates a design may reach.
    """
    from neoswga.core.candidate_provider import CandidateProvider

    path = _record(tmp_path, [A, C, G], STRICT)

    with CandidateInventory(path) as inventory:
        provider = CandidateProvider(inventory, COND, [12])
        assert sorted(provider.initial(10)) == sorted([A, C, G])
