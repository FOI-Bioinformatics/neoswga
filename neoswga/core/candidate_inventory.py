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

import hashlib
import json
import sqlite3
from pathlib import Path
from typing import Any, Dict, Iterator, Sequence

# Bumped when the meaning of a hard gate changes, so a verdict reached under the
# old rules is not silently reused under the new ones.
DEFAULT_POLICY_VERSION = "qc-2026-09-15"


def qc_policy_fingerprint(thresholds: Dict[str, Any]) -> str:
    """A stable name for one set of resolved hard-QC thresholds.

    `DEFAULT_POLICY_VERSION` alone said which generation of RULES was in force,
    not which SETTINGS. Two runs with different GC windows, Tm windows or
    frequency limits shared it, so a verdict reached under one was returned to a
    caller asking about the other.

    Sorted keys, so a mapping built in a different order is the same policy. An
    absent threshold is encoded distinctly from a present one, because "not
    configured" and "configured to the default" are different statements about
    what was enforced.
    """
    material = json.dumps(
        {key: thresholds[key] for key in sorted(thresholds)},
        sort_keys=True,
        default=str,
    )
    digest = hashlib.sha256(material.encode()).hexdigest()[:16]
    return f"{DEFAULT_POLICY_VERSION}:{digest}"


STAGE3_POLICY_PREFIX = "stage3:"

_SCHEMA = """
CREATE TABLE IF NOT EXISTS run_facts (
    condition_id TEXT NOT NULL,
    name TEXT NOT NULL,
    value INTEGER,
    PRIMARY KEY (condition_id, name)
);
CREATE TABLE IF NOT EXISTS candidates (
    sequence TEXT PRIMARY KEY,
    length INTEGER NOT NULL,
    metrics_json TEXT NOT NULL,
    search_rank REAL
);
CREATE TABLE IF NOT EXISTS assessments (
    sequence TEXT NOT NULL REFERENCES candidates(sequence),
    condition_id TEXT NOT NULL,
    policy_version TEXT NOT NULL,
    passed INTEGER NOT NULL,
    reasons_json TEXT NOT NULL,
    metrics_json TEXT NOT NULL,
    generation INTEGER NOT NULL DEFAULT 1,
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
        # An inventory written before generations exists in the wild. Adding the
        # column with a default of 1 makes every old row one coherent
        # generation, which is what it was.
        self._connection.executescript(_SCHEMA)
        candidate_columns = {
            row[1] for row in self._connection.execute("PRAGMA table_info(candidates)")
        }
        if "search_rank" not in candidate_columns:
            self._connection.execute("ALTER TABLE candidates ADD COLUMN search_rank REAL")
        columns = {row[1] for row in self._connection.execute("PRAGMA table_info(assessments)")}
        if "generation" not in columns:
            self._connection.execute(
                "ALTER TABLE assessments ADD COLUMN generation INTEGER NOT NULL DEFAULT 1"
            )
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

    def record_candidate(
        self, sequence: str, metrics: Dict[str, Any], search_rank: float | None = None
    ) -> None:
        """Record a candidate and its condition-independent measurements.

        Called before any condition-specific admission, so that what was
        enumerated is known even for candidates that no chemistry admits.
        Re-recording replaces the measurements: a later pass has better numbers
        than an earlier one, and keeping the first would make the inventory
        depend on scan order.

        `search_rank` is the order a search should walk these in, lowest first.
        It is a COLUMN rather than a metric because the traversal sorts by it
        over the whole inventory: reading it out of `metrics_json` meant one
        query and one JSON parse per candidate, which is 491,836 round trips on
        the Wolbachia design before the search has looked at anything.
        """
        sequence = sequence.upper()
        self._connection.execute(
            "INSERT INTO candidates (sequence, length, metrics_json, search_rank) "
            "VALUES (?, ?, ?, ?) "
            "ON CONFLICT(sequence) DO UPDATE SET metrics_json = excluded.metrics_json, "
            "search_rank = excluded.search_rank",
            (
                sequence,
                len(sequence),
                json.dumps(metrics, sort_keys=True),
                None if search_rank is None else float(search_rank),
            ),
        )

    def record_assessment(
        self,
        sequence: str,
        condition_id: str,
        passed: bool,
        reasons: Sequence[str],
        metrics: Dict[str, Any],
        policy_version: str = DEFAULT_POLICY_VERSION,
        generation: int = 1,
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
            "(sequence, condition_id, policy_version, passed, reasons_json, "
            "metrics_json, generation) "
            "VALUES (?, ?, ?, ?, ?, ?, ?) "
            "ON CONFLICT(sequence, condition_id, policy_version) DO UPDATE SET "
            "passed = excluded.passed, reasons_json = excluded.reasons_json, "
            "metrics_json = excluded.metrics_json, generation = excluded.generation",
            (
                sequence.upper(),
                condition_id,
                policy_version,
                1 if passed else 0,
                json.dumps(list(reasons)),
                json.dumps(metrics, sort_keys=True),
                int(generation),
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

        Ordered by `search_rank`, then by sequence. An unordered scan would make
        an otherwise deterministic run depend on SQLite's row order, and the
        optimizers this feeds are order-sensitive; the sequence tie-break is
        what makes the order total when two candidates rank equally.

        An unranked candidate sorts last rather than disappearing. It cleared
        the hard gates, so it is eligible, and eligibility is not the ranking's
        to revoke.

        Only the NEWEST generation for this condition and policy counts. A
        candidate the latest run did not admit is not eligible, even though its
        earlier verdict is still on record: re-running the filter with a
        stricter rule used to leave everything the looser rule had admitted
        selectable, under rules that no longer existed.
        """
        if not lengths:
            return
        placeholders = ",".join("?" for _ in lengths)
        rows = self._connection.execute(
            "SELECT c.sequence FROM candidates c "
            "JOIN assessments a ON a.sequence = c.sequence "
            "WHERE a.condition_id = ? AND a.policy_version = ? AND a.passed = 1 "
            "AND a.generation = (SELECT MAX(generation) FROM assessments "
            "                    WHERE condition_id = ? AND policy_version = ?) "
            f"AND c.length IN ({placeholders}) "
            "ORDER BY c.search_rank IS NULL, c.search_rank, c.sequence",
            (
                condition_id,
                policy_version,
                condition_id,
                policy_version,
                *[int(n) for n in lengths],
            ),
        )
        for (sequence,) in rows:
            yield sequence

    def record_run_fact(self, condition_id: str, name: str, value) -> None:
        """A whole-run number that is not a property of any one candidate.

        How many k-mers were enumerated is the motivating case: candidates are
        only written once they clear hard QC, so the size of the universe they
        were drawn from cannot be recovered by counting rows.
        """
        self._connection.execute(
            "INSERT INTO run_facts (condition_id, name, value) VALUES (?, ?, ?) "
            "ON CONFLICT(condition_id, name) DO UPDATE SET value = excluded.value",
            (condition_id, name, None if value is None else int(value)),
        )

    def run_fact(self, condition_id: str, name: str):
        """One recorded run number, or `None` if nobody recorded it.

        `None` means unmeasured. It is deliberately not filled in with the
        nearest available number: a plausible substitute cannot be told apart
        from a measurement, which is how the filter funnel hid a 90% loss.
        """
        row = self._connection.execute(
            "SELECT value FROM run_facts WHERE condition_id = ? AND name = ?",
            (condition_id, name),
        ).fetchone()
        return None if row is None else row[0]

    def current_policy(self, condition_id: str) -> str | None:
        """The hard-QC policy the newest recording for this condition used.

        A reader cannot reconstruct the digest: it is a hash of thresholds
        resolved at run time, after the GC-adaptive strategy has had its say. So
        it has to be discoverable, or a caller asking under the old constant
        would get an empty eligible set -- a silent zero, from the same family
        as Known Issues 5, 6 and 13.

        Stage-3 policies are excluded: those record what an efficacy filter
        carried forward, under their own prefix, and are not admission rules.
        """
        row = self._connection.execute(
            "SELECT policy_version FROM assessments "
            "WHERE condition_id = ? AND policy_version NOT LIKE ? "
            "ORDER BY generation DESC, policy_version DESC LIMIT 1",
            (condition_id, f"{STAGE3_POLICY_PREFIX}%"),
        ).fetchone()
        return row[0] if row else None

    def next_generation(self, condition_id: str, policy_version: str) -> int:
        """The generation number a new recording should write.

        One more than the highest already present for this condition and policy,
        so the newest run supersedes without deleting what came before.
        """
        row = self._connection.execute(
            "SELECT MAX(generation) FROM assessments "
            "WHERE condition_id = ? AND policy_version = ?",
            (condition_id, policy_version),
        ).fetchone()
        return int(row[0] or 0) + 1

    def assessment_history(self, sequence: str, condition_id: str):
        """Every verdict recorded for one candidate under one reaction.

        Superseded rows are retired rather than deleted, because "this was
        admitted once, under these rules" is worth knowing when a design has to
        be explained. Returns newest first.
        """
        rows = self._connection.execute(
            "SELECT policy_version, passed, reasons_json, generation FROM assessments "
            "WHERE sequence = ? AND condition_id = ? ORDER BY generation DESC",
            (sequence.upper(), condition_id),
        )
        return [
            {
                "policy_version": policy,
                "passed": bool(passed),
                "reasons": json.loads(reasons),
                "generation": generation,
            }
            for policy, passed, reasons, generation in rows
        ]

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

    def has_length(self, length: int) -> bool:
        """Whether any candidate of this length was ever enumerated.

        Distinct from "any candidate of this length is eligible". A length
        nothing was counted at is a missing input; a length where every
        candidate failed a reaction's gates is a result.
        """
        row = self._connection.execute(
            "SELECT 1 FROM candidates WHERE length = ? LIMIT 1", (int(length),)
        ).fetchone()
        return row is not None

    def counts(self, condition_id: str | None = None) -> Dict[str, int]:
        """How many candidates were written, and how many currently qualify.

        `hard_qc_passed` is scoped to one condition, its current policy and the
        newest generation, because that is what eligibility means. Counting
        every passing row in the database instead made it a running total over
        the run's own history: a re-run that admitted fewer candidates reported
        more, since the superseded verdicts were still there.

        Note what `candidates` is not. Rows are written only once a candidate
        has cleared hard QC, so this is the survivors, and the size of the
        universe they came from is a run fact rather than a row count.
        """
        counted = self._connection.execute("SELECT COUNT(*) FROM candidates").fetchone()[0]
        assessed = self._connection.execute(
            "SELECT COUNT(DISTINCT sequence) FROM assessments"
        ).fetchone()[0]
        if condition_id is None:
            passed = self._connection.execute(
                "SELECT COUNT(DISTINCT sequence) FROM assessments WHERE passed = 1"
            ).fetchone()[0]
        else:
            policy = self.current_policy(condition_id)
            passed = self._connection.execute(
                "SELECT COUNT(DISTINCT sequence) FROM assessments "
                "WHERE condition_id = ? AND policy_version = ? AND passed = 1 "
                "AND generation = (SELECT MAX(generation) FROM assessments "
                "                  WHERE condition_id = ? AND policy_version = ?)",
                (condition_id, policy, condition_id, policy),
            ).fetchone()[0]
        return {"candidates": counted, "assessed": assessed, "hard_qc_passed": passed}


STAGE2_INVENTORY_NAME = "candidate_inventory.sqlite"


def _search_ranks(cleared_hard_gates, shortlisted) -> Dict[str, float]:
    """The order a search should walk the whole hard-QC set in, lowest first.

    Two bands. The shortlist leads, in the order stage 2 put it in, because that
    ranking already weighed occupancy and background load over the candidates it
    considered and a search should start from the working set rather than
    rediscover it. Everything else follows, ordered by background sites per
    target site -- the same quantity `_rank_and_cut_candidates` ranks on, and
    the only one available for candidates the cut never scored.

    Written for EVERY hard-QC candidate. A rank that exists only for the
    shortlist cannot order the candidates outside it, which is precisely the set
    an expansion draws from.

    A candidate with no usable counts gets no rank and sorts last, rather than
    being given a number that would place it somewhere it has not earned.
    """
    ranks: Dict[str, float] = {}
    for position, sequence in enumerate(shortlisted["primer"]):
        ranks[str(sequence).upper()] = float(position)

    offset = float(len(ranks))
    rest = []
    for row in cleared_hard_gates.to_dict("records"):
        sequence = str(row["primer"]).upper()
        if sequence in ranks:
            continue
        foreground = row.get("fg_count")
        background = row.get("bg_count")
        try:
            if not foreground:
                continue
            rest.append((float(background or 0) / float(foreground), sequence))
        except (TypeError, ValueError):
            continue
    for position, (_, sequence) in enumerate(sorted(rest)):
        ranks[sequence] = offset + position
    return ranks


def record_stage2_inventory(
    data_dir,
    condition_id: str,
    cleared_hard_gates,
    after_gini,
    shortlisted,
    indexed=None,
    policy_version: str = DEFAULT_POLICY_VERSION,
    qc_policy: Dict[str, Any] | None = None,
    enumerated: int | None = None,
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
    if len(cleared_hard_gates) == 0:
        raise ValueError(
            "Refusing to record an inventory with no candidates. An empty run "
            "would open a new generation holding nothing, retiring every verdict "
            "the previous run reached and leaving the design with no eligible "
            "candidates at all."
        )
    # Which rules this verdict was reached under. Without it a stricter re-run
    # writes into the same key space as the looser one it replaced.
    if qc_policy is not None:
        policy_version = qc_policy_fingerprint(qc_policy)
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

    ranks = _search_ranks(cleared_hard_gates, shortlisted)

    with CandidateInventory(path) as inventory:
        generation = inventory.next_generation(condition_id, policy_version)
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
            metrics["search_rank"] = ranks.get(sequence)
            inventory.record_candidate(sequence, metrics, search_rank=ranks.get(sequence))
            inventory.record_assessment(
                sequence,
                condition_id,
                passed=True,
                reasons=[],
                metrics={},
                policy_version=policy_version,
                generation=generation,
            )
        if enumerated is not None:
            inventory.record_run_fact(condition_id, "enumerated", enumerated)
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


def background_scan_pool(cleared_hard_gates, retention, *, after_gini=None):
    """Which candidates get a background position index.

    Foreground positions are already written for every hard-QC survivor, because
    the Gini gate reads them. Background positions used to be written only for
    the `max_primer` shortlist, so a candidate the ranking cut had no background
    index at all -- and a design reaching for it later would score it against an
    empty background, which reads as perfect specificity rather than as a
    missing measurement.

    Both surviving modes close that, and both change what is INDEXED rather than
    what the optimizer searches: `step2_df.csv` remains the shortlist, so
    runtime moves with a search setting and not with a retention setting.

    `all_qc` indexes every candidate that cleared the declared hard gates.

    `post_gini` also requires the evenness gate, and the distinction it draws is
    between a GATE and a RANKING. The Gini gate is a declared requirement that a
    candidate either meets or does not. `max_primer` is a cut through a ranking,
    chosen for the size of the working set rather than for any property of the
    candidates below the line. Keeping the first and dropping the second retains
    what was set aside arbitrarily and not what was set aside on a stated rule.
    It does not reopen the silent zero, because the provider does not expand
    past a hard gate.

    Measured on the Wolbachia design, 20,670 candidates are indexed under
    `post_gini` and 491,836 under `all_qc`, costing about 19 MB against 443 MB.

    The shortlist-only `legacy` mode was removed on 2026-09-16. It existed to
    reproduce the historical truncation for comparison, that comparison has been
    made and published, and keeping a mode whose only property is a measurement
    fault is how the fault comes back.
    """
    if retention == "all_qc":
        return list(cleared_hard_gates)
    if retention == "post_gini":
        if after_gini is None:
            raise ValueError(
                "candidate_retention='post_gini' needs the post-Gini candidates; "
                "without them the mode would quietly become a different one"
            )
        return list(after_gini)
    if retention == "legacy":
        raise ValueError(
            "candidate_retention='legacy' was removed on 2026-09-16. It indexed "
            "only the max_primer shortlist, so any candidate the ranking cut "
            "scored against an empty background and read as perfectly specific. "
            "Use 'post_gini' for the nearest smaller index, or 'all_qc'."
        )
    raise ValueError(f"candidate_retention must be 'all_qc' or 'post_gini', not {retention!r}")


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

    `counted` is the enumerated universe and `assessed` the survivors written
    here, so they differ when enumeration outruns judgement. `hard_qc_passed` is
    eligibility under this condition, its current policy and the newest
    generation. `shortlisted` and `indexed` are what the ranking and the
    retention policy did. `carried_to_stage3` is what an efficacy filter kept.
    `examined` is what a search actually evaluated. Reporting one of these and
    labelling it as the pool is what made a 90% ranking loss invisible.

    `counted` and `examined` are `None` when nobody recorded them, which is the
    current state of `examined`: no search reports its evaluations yet. An
    absent count says nobody measured it; a substituted one cannot be told from
    a measurement.
    """
    with CandidateInventory(path) as inventory:
        base = inventory.counts(condition_id)
        shortlisted = indexed = 0
        for (metrics_json,) in inventory._connection.execute("SELECT metrics_json FROM candidates"):
            metrics = json.loads(metrics_json)
            shortlisted += 1 if metrics.get("shortlisted") else 0
            indexed += 1 if metrics.get("indexed") else 0
        carried = inventory._connection.execute(
            "SELECT COUNT(DISTINCT sequence) FROM assessments "
            "WHERE condition_id = ? AND passed = 1 AND policy_version LIKE ?",
            (condition_id, f"{STAGE3_POLICY_PREFIX}%"),
        ).fetchone()[0]
        counted = inventory.run_fact(condition_id, "enumerated")
        examined = inventory.run_fact(condition_id, "examined")
    return {
        "counted": counted,
        "assessed": base["assessed"],
        "hard_qc_passed": base["hard_qc_passed"],
        "shortlisted": shortlisted,
        "indexed": indexed,
        "carried_to_stage3": carried,
        "examined": examined,
    }


def report_design_counts(path, condition_id: str) -> Dict[str, int]:
    """The six distinctions, by name, in one line.

    The funnel reported a handful of counts and none of them answered how many
    candidates the optimizer could actually have chosen from.
    """
    import logging

    counts = design_counts(path, condition_id)

    def shown(name):
        """`unknown` rather than a zero, which would read as a measurement."""
        value = counts[name]
        return "unknown" if value is None else str(value)

    logging.getLogger(__name__).info(
        "Design counts: counted=%s assessed=%s hard_qc_passed=%s shortlisted=%s "
        "indexed=%s carried_to_stage3=%s examined=%s",
        shown("counted"),
        shown("assessed"),
        shown("hard_qc_passed"),
        shown("shortlisted"),
        shown("indexed"),
        shown("carried_to_stage3"),
        shown("examined"),
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
