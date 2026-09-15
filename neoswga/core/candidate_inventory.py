"""A durable record of every candidate, and why each one is or is not eligible.

Task 3 of the condition-aware pool design plan.

The pipeline used to hand the optimizer a CSV that had already been truncated by
a ranking cap. On the measured Wolbachia design 20,670 candidates survived the
Gini filter and `max_primer` kept 2,000, so 90.3% were deleted by a quality
ORDER rather than by any stated requirement. An optimizer cannot select what it
was never given, and that CSV cannot be described as the candidate universe
because it is a truncation of one.

This module separates two things the pipeline had merged:

hard gate
    A declared requirement: sequence composition, self-dimer limits, exclusion
    genomes, and the configured frequency and background limits. Failing one
    makes a candidate ineligible, and the reason is recorded rather than
    implied by its absence.

measurement
    Gini, abundance, rank, occupancy. These order a search. They are not
    grounds for permanent deletion: a primer that only covers a difficult
    region ranks poorly and may still be the one a pool needs.

Eligibility is stored per reaction, because it is a property of a candidate
UNDER A CHEMISTRY. The same oligo can fail a Tm window under one condition and
pass under another, and both verdicts are worth keeping. Assessments therefore
key on the condition fingerprint and the QC-policy version together, so a
revised policy does not overwrite the verdict the old one reached.

Storage is SQLite from the standard library. The point is streamed inserts and
transactions rather than holding the whole inventory in a dataframe, which is
what made truncation look necessary in the first place.
"""

from __future__ import annotations

import json
import sqlite3
from pathlib import Path
from typing import Any, Dict, Iterator, Sequence

# Bumped when the meaning of a hard gate changes, so a verdict reached under the
# old rules is not silently reused under the new ones.
DEFAULT_POLICY_VERSION = "qc-2026-09-15"

_SCHEMA = """
CREATE TABLE IF NOT EXISTS candidates (
    sequence TEXT PRIMARY KEY,
    length INTEGER NOT NULL,
    metrics_json TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS assessments (
    sequence TEXT NOT NULL REFERENCES candidates(sequence),
    condition_id TEXT NOT NULL,
    policy_version TEXT NOT NULL,
    passed INTEGER NOT NULL,
    reasons_json TEXT NOT NULL,
    metrics_json TEXT NOT NULL,
    PRIMARY KEY (sequence, condition_id, policy_version)
);
CREATE INDEX IF NOT EXISTS assessments_lookup
    ON assessments (condition_id, policy_version, passed);
"""


class CandidateInventory:
    """Durable candidate storage for one design.

    Use as a context manager; it owns the connection and commits on exit.

        with CandidateInventory(path) as inventory:
            inventory.record_candidate(sequence, metrics)
    """

    def __init__(self, path: Path | str):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._connection = sqlite3.connect(str(self.path))
        self._connection.executescript(_SCHEMA)
        self._connection.commit()

    # -- lifecycle ---------------------------------------------------------

    def __enter__(self) -> "CandidateInventory":
        return self

    def __exit__(self, *exc_info) -> None:
        self.close()

    def close(self) -> None:
        try:
            self._connection.commit()
        finally:
            self._connection.close()

    # -- writing -----------------------------------------------------------

    def record_candidate(self, sequence: str, metrics: Dict[str, Any]) -> None:
        """Record a candidate and its condition-independent measurements.

        Called before any condition-specific admission, so that what was
        enumerated is known even for candidates that no chemistry admits.
        Re-recording replaces the measurements: a later pass has better numbers
        than an earlier one, and keeping the first would make the inventory
        depend on scan order.
        """
        sequence = sequence.upper()
        self._connection.execute(
            "INSERT INTO candidates (sequence, length, metrics_json) VALUES (?, ?, ?) "
            "ON CONFLICT(sequence) DO UPDATE SET metrics_json = excluded.metrics_json",
            (sequence, len(sequence), json.dumps(metrics, sort_keys=True)),
        )

    def record_assessment(
        self,
        sequence: str,
        condition_id: str,
        passed: bool,
        reasons: Sequence[str],
        metrics: Dict[str, Any],
        policy_version: str = DEFAULT_POLICY_VERSION,
    ) -> None:
        """Record whether a candidate clears the hard gates under one reaction.

        `reasons` is why it did not, and is kept even when it did, so a later
        reader can tell "passed" from "never assessed" without inferring it from
        a missing row.

        `policy_version` is not in the plan's stated signature but is in its
        schema, where it forms part of the primary key. It is a keyword with a
        default so existing callers need not pass it, and so that a verdict
        reached under revised QC rules does not overwrite the old one.
        """
        self._connection.execute(
            "INSERT INTO assessments "
            "(sequence, condition_id, policy_version, passed, reasons_json, metrics_json) "
            "VALUES (?, ?, ?, ?, ?, ?) "
            "ON CONFLICT(sequence, condition_id, policy_version) DO UPDATE SET "
            "passed = excluded.passed, reasons_json = excluded.reasons_json, "
            "metrics_json = excluded.metrics_json",
            (
                sequence.upper(),
                condition_id,
                policy_version,
                1 if passed else 0,
                json.dumps(list(reasons)),
                json.dumps(metrics, sort_keys=True),
            ),
        )

    def commit(self) -> None:
        """End the current batch. A committed batch survives an interrupted run."""
        self._connection.commit()

    # -- reading -----------------------------------------------------------

    def iter_eligible(
        self,
        condition_id: str,
        lengths: Sequence[int],
        policy_version: str = DEFAULT_POLICY_VERSION,
    ) -> Iterator[str]:
        """Every candidate of the requested lengths passing the hard gates.

        Ordered by sequence. An unordered scan would make an otherwise
        deterministic run depend on SQLite's row order, and the optimizers this
        feeds are order-sensitive.
        """
        if not lengths:
            return
        placeholders = ",".join("?" for _ in lengths)
        rows = self._connection.execute(
            "SELECT c.sequence FROM candidates c "
            "JOIN assessments a ON a.sequence = c.sequence "
            f"WHERE a.condition_id = ? AND a.policy_version = ? AND a.passed = 1 "
            f"AND c.length IN ({placeholders}) ORDER BY c.sequence",
            (condition_id, policy_version, *[int(n) for n in lengths]),
        )
        for (sequence,) in rows:
            yield sequence

    def metrics(self, sequence: str) -> Dict[str, Any]:
        """The condition-independent measurements recorded for one candidate."""
        row = self._connection.execute(
            "SELECT metrics_json FROM candidates WHERE sequence = ?", (sequence.upper(),)
        ).fetchone()
        return json.loads(row[0]) if row else {}

    def reasons(
        self,
        sequence: str,
        condition_id: str,
        policy_version: str = DEFAULT_POLICY_VERSION,
    ) -> list:
        """Why a candidate is ineligible under one reaction."""
        row = self._connection.execute(
            "SELECT reasons_json FROM assessments WHERE sequence = ? AND condition_id = ? "
            "AND policy_version = ?",
            (sequence.upper(), condition_id, policy_version),
        ).fetchone()
        return json.loads(row[0]) if row else []

    def counts(self) -> Dict[str, int]:
        """Counted and assessed are different numbers, and both are reported.

        The funnel this replaces conflated them: its first row was labelled
        `total_kmers` but began after several loading prefilters, so the number
        of candidates enumerated was never stated anywhere.
        """
        counted = self._connection.execute("SELECT COUNT(*) FROM candidates").fetchone()[0]
        assessed = self._connection.execute(
            "SELECT COUNT(DISTINCT sequence) FROM assessments"
        ).fetchone()[0]
        passed = self._connection.execute(
            "SELECT COUNT(DISTINCT sequence) FROM assessments WHERE passed = 1"
        ).fetchone()[0]
        return {"candidates": counted, "assessed": assessed, "hard_qc_passed": passed}


STAGE2_INVENTORY_NAME = "candidate_inventory.sqlite"


def record_stage2_inventory(
    data_dir,
    condition_id: str,
    cleared_hard_gates,
    after_gini,
    shortlisted,
    indexed=None,
    policy_version: str = DEFAULT_POLICY_VERSION,
) -> Path:
    """Write what stage 2 enumerated, and what merely ordered the search.

    Called with three frames from the filter step, narrowing in that order.
    Every row of `cleared_hard_gates` has already passed the declared
    requirements -- frequency limits, the Tm and sequence rules, exclusion
    genomes -- so every one of them is recorded as ELIGIBLE.

    What happens after that is ranking, and this is the distinction the CSV
    could not express. The Gini filter and the `max_primer` cut both remove
    rows, but neither is a stated requirement: Gini measures how evenly a
    primer binds and the cap keeps the top of an order. On the measured
    Wolbachia design the two together removed 90.3% of the survivors, and the
    optimizer could not select any of them. They are recorded here as
    measurements, `passed_gini` and `shortlisted`, against candidates that
    remain addressable.

    Returns the inventory path. Failure to write is not swallowed: a design
    that cannot record what it enumerated should say so rather than leave a
    partial inventory that later reads as complete.
    """
    path = Path(data_dir) / STAGE2_INVENTORY_NAME
    survived_gini = {str(s).upper() for s in after_gini["primer"]}
    kept = {str(s).upper() for s in shortlisted["primer"]}
    # Which candidates got a background position index. Without one, a
    # specificity number computed for a candidate is a zero rather than a
    # measurement, so this is worth recording per candidate and not only as a
    # mode setting.
    have_index = (
        {str(s).upper() for s in indexed}
        if indexed is not None
        else {str(s).upper() for s in cleared_hard_gates["primer"]}
    )

    with CandidateInventory(path) as inventory:
        for row in cleared_hard_gates.to_dict("records"):
            sequence = str(row["primer"]).upper()
            metrics = {
                key: value
                for key, value in row.items()
                if key != "primer" and isinstance(value, (int, float, str, bool))
            }
            metrics["passed_gini"] = sequence in survived_gini
            metrics["shortlisted"] = sequence in kept
            metrics["indexed"] = sequence in have_index
            inventory.record_candidate(sequence, metrics)
            inventory.record_assessment(
                sequence,
                condition_id,
                passed=True,
                reasons=[],
                metrics={},
                policy_version=policy_version,
            )
        inventory.commit()

    import logging

    logging.getLogger(__name__).info(
        "Candidate inventory: %d cleared the hard gates, %d survived the Gini "
        "filter, %d were shortlisted by max_primer. All %d remain addressable (%s).",
        len(cleared_hard_gates),
        len(after_gini),
        len(shortlisted),
        len(cleared_hard_gates),
        path,
    )
    report_design_counts(path, condition_id)
    return path


def background_scan_pool(cleared_hard_gates, shortlisted, retention):
    """Which candidates get a background position index.

    Foreground positions are already written for every hard-QC survivor, because
    the Gini gate reads them. Background positions were written only for the
    shortlist, so a candidate the `max_primer` ranking cut had no background
    index at all -- and a design reaching for it later would score it against an
    empty background, which reads as perfect specificity rather than as a
    missing measurement.

    `all_qc` closes that. It changes what is INDEXED, not what the optimizer
    searches: `step2_df.csv` remains the shortlist, so runtime moves with a
    search setting rather than with a retention setting.
    """
    if retention == "all_qc":
        return list(cleared_hard_gates)
    if retention == "legacy":
        return list(shortlisted)
    raise ValueError(f"candidate_retention must be 'all_qc' or 'legacy', not {retention!r}")


STAGE3_POLICY_PREFIX = "stage3:"


def record_stage3_policy(path, condition_id: str, carried, rejected, policy: str) -> None:
    """Record what stage 3 carried forward, and why it dropped the rest.

    Stage 3's efficacy filter is opt-in and is a POLICY, not a gate. A candidate
    it drops still cleared every declared requirement, so the drop is recorded
    against that candidate under its own policy version rather than by removing
    the hard-QC verdict. Reading a low model score as a QC failure is how a
    ranking came to look like a requirement in the first place.

    `carried` is what reached the optimizer: the `examined` count. `rejected`
    maps a sequence to the reason it did not.
    """
    version = f"{STAGE3_POLICY_PREFIX}{policy}"
    with CandidateInventory(path) as inventory:
        for sequence in carried:
            inventory.record_assessment(
                sequence,
                condition_id,
                passed=True,
                reasons=[],
                metrics={},
                policy_version=version,
            )
        for sequence, reason in dict(rejected).items():
            inventory.record_assessment(
                sequence,
                condition_id,
                passed=False,
                reasons=[reason],
                metrics={},
                policy_version=version,
            )
        inventory.commit()


def design_counts(path, condition_id: str) -> Dict[str, int]:
    """The six distinctions the funnel could not express.

    `counted` and `assessed` differ when enumeration outruns judgement.
    `hard_qc_passed` is eligibility. `shortlisted` and `indexed` are what the
    ranking and the retention policy did. `examined` is what actually reached
    the optimizer. Reporting one of these and labelling it as the pool is what
    made a 90% ranking loss invisible.
    """
    with CandidateInventory(path) as inventory:
        base = inventory.counts()
        shortlisted = indexed = 0
        for (metrics_json,) in inventory._connection.execute("SELECT metrics_json FROM candidates"):
            metrics = json.loads(metrics_json)
            shortlisted += 1 if metrics.get("shortlisted") else 0
            indexed += 1 if metrics.get("indexed") else 0
        examined = inventory._connection.execute(
            "SELECT COUNT(DISTINCT sequence) FROM assessments "
            "WHERE condition_id = ? AND passed = 1 AND policy_version LIKE ?",
            (condition_id, f"{STAGE3_POLICY_PREFIX}%"),
        ).fetchone()[0]
    return {
        "counted": base["candidates"],
        "assessed": base["assessed"],
        "hard_qc_passed": base["hard_qc_passed"],
        "shortlisted": shortlisted,
        "indexed": indexed,
        "examined": examined,
    }


def report_design_counts(path, condition_id: str) -> Dict[str, int]:
    """The six distinctions, by name, in one line.

    The funnel reported a handful of counts and none of them answered how many
    candidates the optimizer could actually have chosen from.
    """
    import logging

    counts = design_counts(path, condition_id)
    logging.getLogger(__name__).info(
        "Design counts: counted=%d assessed=%d hard_qc_passed=%d shortlisted=%d "
        "indexed=%d examined=%d",
        counts["counted"],
        counts["assessed"],
        counts["hard_qc_passed"],
        counts["shortlisted"],
        counts["indexed"],
        counts["examined"],
    )
    return counts


def record_stage3_from_frames(data_dir, use_amp_model, step2_df, carried_df):
    """Record what stage 3 carried, and why it dropped anything.

    The efficacy filter is opt-in and is a POLICY, not a gate: a candidate it
    drops still cleared every declared requirement. Recording the drop against
    the candidate, under its own policy version, keeps that distinction. Reading
    a low model score as a QC failure is how a ranking came to look like a
    requirement.

    A data directory with no inventory -- an older one, or a caller supplying
    its own candidate list -- is not a reason to fail the step, so that case
    returns. A malformed inventory still raises.
    """
    import os

    path = os.path.join(str(data_dir), STAGE2_INVENTORY_NAME)
    if not os.path.exists(path):
        return

    # Resolved only once an inventory is known to exist. Stage 3 never needed
    # reaction conditions before, and resolving them unconditionally coupled it
    # to a chemistry that a caller running stage 3 alone may not have set.
    from neoswga.core.filter import _get_reaction_conditions

    condition_id = _get_reaction_conditions().fingerprint()

    def _primers(frame):
        # `order_step3_rows` leaves the sequence as the index rather than a
        # column, so a frame here may carry it either way.
        if "primer" in frame.columns:
            return [str(s).upper() for s in frame["primer"]]
        return [str(s).upper() for s in frame.index]

    carried = _primers(carried_df)
    carried_set = set(carried)
    dropped = {
        sequence: "removed by the opt-in efficacy filter (--amp-model)"
        for sequence in _primers(step2_df)
        if sequence not in carried_set
    }
    policy = "amp_model" if use_amp_model else "carry_forward"
    record_stage3_policy(path, condition_id, carried=carried, rejected=dropped, policy=policy)
    report_design_counts(path, condition_id)
