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

from neoswga.core.position_cache import MissingPositionsError


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
        # The provider owns the position data for the candidates it hands out,
        # because it is the only thing that knows a batch was admitted. Nothing
        # may be scored before `attach_positions`.
        self.position_cache = None

    # -- ordering ----------------------------------------------------------

    def _eligible_in_search_order(self):
        """Every eligible candidate, in the order the inventory recorded.

        The order is the inventory's to decide and this walks it. It used to
        sort here on a `step2_rank` metric that stage 2 never wrote -- it is
        created in stage 3, on the shortlist, after the inventory is recorded --
        so every candidate fell into one bucket and the traversal was
        alphabetical while the docstring claimed otherwise.
        """
        if self._ordered is None:
            # Ask the inventory which policy it last recorded under. The digest
            # is over thresholds resolved at run time, so a caller cannot
            # reconstruct it, and asking under the bare constant would return an
            # empty set that reads as "no candidates qualify".
            policy = self.inventory.current_policy(self.condition_id)
            arguments = [policy] if policy else []
            self._ordered = list(
                self.inventory.iter_eligible(self.condition_id, self.lengths, *arguments)
            )
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

    # -- positions ---------------------------------------------------------

    def attach_positions(self, cache) -> None:
        """Hand the provider the cache the candidates it serves are scored on.

        Mandatory before any batch is scored. Attaching is also what makes
        loading on demand possible, since the provider is the only party that
        knows a candidate has just been admitted.
        """
        self.position_cache = cache

    def ensure_positions(self, sequences: Sequence[str]) -> None:
        """Make sure a batch has position data before it is scored.

        Two things had to change here.

        It used to return quietly when no cache was attached, which is the one
        configuration that cannot be checked at all. Nothing in production
        attached one, so the guard was inert on every run it was supposed to
        protect. It now raises.

        And it asked whether a candidate had a hit on ANY prefix. That is the
        wrong question twice over. A candidate with fifty foreground sites and
        no background entry passed, which is unknown specificity reported as
        perfect specificity -- the silent zero of Known Issues 5, 6 and 13. A
        candidate indexed against a host it binds nowhere failed, though its
        zero is a measurement and a good one. The question is whether there is
        an ENTRY on EVERY prefix the design will score against.

        A candidate the cache does not hold yet is loaded rather than refused:
        that is the normal case once the frontier moves, and the cache applies
        its own `on_missing` policy to the batch.
        """
        cache = self.position_cache
        if cache is None:
            raise MissingPositionsError(
                "This provider has no position cache, so it cannot tell a candidate "
                "that binds nowhere from one that was never indexed. Call "
                "attach_positions(cache) before scoring a batch."
            )
        cache.require_entries([str(sequence) for sequence in sequences])
