"""A working shortlist that can grow, over an inventory that already holds more.

Task 4 of the condition-aware pool design plan.

Task 3 made every candidate clearing the declared hard gates durable and
background-indexed. This makes them reachable. `max_primer` seeds a shortlist,
and when that shortlist cannot meet the constraints the search asks for another
batch rather than concluding the design is infeasible.

That distinction is the point. From inside a greedy that has run out of
candidates, "nothing left that would help" and "nothing left in the batch I was
handed" produce the same stall, and the pipeline reported both as the second.
`unexamined` and `exhausted` let a caller say which one happened.

Ordering is deterministic throughout. The optimizers this feeds are
order-sensitive, so a traversal that depended on SQLite's row order would make
an otherwise reproducible run depend on storage layout.
"""

from __future__ import annotations

from typing import Iterable, List, Sequence, Set


class CandidateProvider:
    """Hands out eligible candidates in batches, seeded by the existing ranking.

    The inventory is the source of truth for what exists. This class decides
    only what the search looks at next, so a budget bounds examination and never
    eligibility.
    """

    def __init__(self, inventory, condition_id: str, lengths: Sequence[int]):
        self.inventory = inventory
        self.condition_id = condition_id
        self.lengths = list(lengths)
        self._ordered: List[str] | None = None

    # -- ordering ----------------------------------------------------------

    def _eligible_in_search_order(self) -> List[str]:
        """Every eligible candidate, best-ranked first.

        Seeded by `step2_rank`, the order stage 2 established, so the shortlist
        this yields first is the one the pipeline would previously have written
        and nothing about a default run changes. Rank still orders the search;
        it no longer bounds it.

        Candidates with no recorded rank sort after those with one, and ties
        break on the sequence, so the traversal is total and reproducible rather
        than dependent on storage order.
        """
        if self._ordered is None:
            sequences = list(self.inventory.iter_eligible(self.condition_id, self.lengths))
            ranked = []
            for sequence in sequences:
                rank = self.inventory.metrics(sequence).get("step2_rank")
                ranked.append(
                    (0, float(rank), sequence)
                    if isinstance(rank, (int, float))
                    else (1, 0.0, sequence)
                )
            self._ordered = [sequence for _, _, sequence in sorted(ranked)]
        return self._ordered

    # -- batches -----------------------------------------------------------

    def initial(self, limit: int) -> List[str]:
        """The starting shortlist: the best-ranked `limit` eligible candidates."""
        return self._eligible_in_search_order()[: max(0, int(limit))]

    def expand(self, excluded: Iterable[str], limit: int) -> List[str]:
        """The next batch of eligible candidates not yet examined.

        Returns an empty list only when the eligible inventory is exhausted, so
        an empty batch is a real answer rather than a budget artefact.
        """
        seen: Set[str] = {str(s).upper() for s in excluded}
        batch = []
        for sequence in self._eligible_in_search_order():
            if sequence.upper() in seen:
                continue
            batch.append(sequence)
            if len(batch) >= max(0, int(limit)):
                break
        return batch

    # -- reporting ---------------------------------------------------------

    def unexamined(self, excluded: Iterable[str]) -> int:
        """How many eligible candidates the search has not looked at.

        A stop with this above zero is a budget stop. Reporting it is what
        separates "the search ran out of time" from "no panel exists", which the
        pipeline previously conflated.
        """
        seen = {str(s).upper() for s in excluded}
        return sum(1 for s in self._eligible_in_search_order() if s.upper() not in seen)

    def exhausted(self, excluded: Iterable[str]) -> bool:
        """True when every eligible candidate has been examined."""
        return self.unexamined(excluded) == 0

    def ensure_positions(self, sequences: Sequence[str]) -> None:
        """Make sure a batch has position data before it is scored.

        Under `candidate_retention="all_qc"` every eligible candidate is indexed
        at stage 2, so this is a no-op for the common path and exists so a
        caller need not know that. For an index built under a narrower policy,
        or for a candidate supplied from outside the pipeline, it raises rather
        than letting a missing index read as a primer that binds nowhere -- the
        silent-zero shape recorded in Known Issues 5, 6 and 13.
        """
        cache = getattr(self, "position_cache", None)
        if cache is None:
            return
        missing = [
            sequence
            for sequence in sequences
            if not any(
                len(cache.get_positions(prefix, sequence, "both"))
                for prefix in getattr(cache, "fname_prefixes", [])
            )
        ]
        if missing:
            raise ValueError(
                f"{len(missing)} candidate(s) in this batch have no position data "
                f"(first: {missing[0]}). Under candidate_retention='all_qc' stage 2 "
                f"indexes every eligible candidate; 'post_gini' covers only those "
                f"that also cleared the evenness gate, and an index written before "
                f"2026-09-16 may cover only the max_primer shortlist. Re-run "
                f"'neoswga filter' with all_qc retention, or the coverage computed "
                f"for these is a zero rather than a measurement."
            )
