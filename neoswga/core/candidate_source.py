"""Where a design's candidates come from.

Phase 4 of the plan for `docs/validation/pipeline_audit_2026-09-16/`, audit
finding F1.

`plan-pool`, `optimize` and `expand-primers` each read `step3_df.csv` and build
a list. The SQLite inventory beside it holds every candidate that cleared hard
QC -- 491,836 on the Wolbachia design against the CSV's 2,000 -- and nothing in
production reads it. So the retained candidates are stored, indexed, and unable
to affect any delivered panel.

This is the seam that closes it: one thing a design asks for candidates,
whichever they come from. Two implementations today.

`ListCandidateSource` wraps a list. It backs `--candidates`, a run directory
with no inventory, and every existing caller and test, so nothing has to change
to keep working.

`InventoryCandidateSource` reads the inventory in the order stage 2 recorded.
The search works over a FRONTIER, a bounded ranked window, rather than the whole
universe: scoring every candidate at every greedy step is not affordable at this
scale, and `advance()` is how the window moves. The universe is what a design
MAY reach; the frontier is what it has looked at so far. Those are different
numbers and a report that quotes one as the other is the failure this whole
audit is about.

Expansion arrives in a later increment. This one starts the frontier at exactly
the shortlist the CSV held, so no delivered panel moves, and adds the counts
that make the difference visible.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, List, Optional, Sequence

from neoswga.core.candidate_inventory import STAGE2_INVENTORY_NAME
from neoswga.core.position_cache import MissingPositionsError

# Why a search ran out of candidates, and they call for opposite responses.
#
# `inventory_exhausted` used to be reported for both, which is true of the batch
# a search was handed and false of the inventory behind it. A caller told the
# inventory is empty stops; a caller told the frontier is empty refills and
# continues. On the Wolbachia design those differ by 18,670 candidates.
FRONTIER_EXHAUSTED = "frontier_exhausted"
INVENTORY_EXHAUSTED = "inventory_exhausted"

# How many candidates a refill adds when the frontier was opened at a size.
# A refill doubles rather than adding a fixed batch, so reaching a distant
# candidate costs a logarithmic number of refills rather than a linear one.
FRONTIER_GROWTH = 2


class ListCandidateSource:
    """A fixed list of candidates, exhausted as soon as it is handed over.

    The universe and the frontier are the same set here, which is exactly what
    a `--candidates` file means: the user named the pool, and there is nothing
    else to reach for.
    """

    kind = "list"

    def __init__(self, sequences: Sequence[str]):
        self._sequences = list(dict.fromkeys(str(s).upper() for s in sequences))
        self._examined: List[str] = []

    def initial(self, limit: Optional[int] = None) -> List[str]:
        """The starting frontier."""
        self._examined = self._sequences if limit is None else self._sequences[: max(0, limit)]
        return list(self._examined)

    def frontier(self) -> List[str]:
        return list(self._examined)

    def advance(self, keep: Iterable[str] = ()) -> bool:
        """There is nothing beyond what was handed over."""
        return False

    def universe_size(self) -> int:
        return len(self._sequences)

    def examined(self) -> int:
        return len(self._examined)

    def exhausted(self) -> bool:
        return True

    def exhaustion(self) -> str:
        """Always the inventory: a named pool has nothing behind it."""
        return INVENTORY_EXHAUSTED

    def attach_positions(self, cache) -> None:
        self.position_cache = cache

    def ensure_positions(self, sequences: Sequence[str]) -> None:
        """Refuse a candidate this design cannot score on every reference.

        A `--candidates` list is the case where this matters most: those
        sequences came from outside the pipeline, so nothing has established
        that either reference was ever scanned for them.
        """
        cache = getattr(self, "position_cache", None)
        if cache is None:
            raise MissingPositionsError(
                "No position cache is attached to this candidate source, so a "
                "candidate that binds nowhere cannot be told from one that was "
                "never indexed. Call attach_positions(cache) first."
            )
        cache.require_entries([str(sequence) for sequence in sequences])

    def describe(self) -> dict:
        return {
            "kind": self.kind,
            "universe": self.universe_size(),
            "frontier": len(self._examined),
            "examined": self.examined(),
            "exhausted": True,
            "rank_key": "as supplied",
        }


class InventoryCandidateSource:
    """Every candidate the inventory says a design may reach, in recorded order.

    The order is the inventory's `search_rank`: the shortlist first, in the
    order stage 2 put it in, then the rest by background load per target site.
    The source does not re-rank, because two rankings for one traversal drift
    while each stays self-consistent.
    """

    kind = "inventory"

    def __init__(self, provider, frontier: Optional[int] = None):
        self._provider = provider
        self._frontier_size = frontier
        self._examined: List[str] = []

    def initial(self, limit: Optional[int] = None) -> List[str]:
        size = limit if limit is not None else self._frontier_size
        ordered = self._provider._eligible_in_search_order()
        self._examined = list(ordered if size is None else ordered[: max(0, size)])
        return list(self._examined)

    def frontier(self) -> List[str]:
        return list(self._examined)

    def advance(self, keep: Iterable[str] = ()) -> bool:
        """Widen the frontier over the next candidates in recorded order.

        Returns False only when the frontier already covers the universe, so a
        False is a real answer about the inventory rather than about a batch.

        The order is the inventory's and this walks it, so a refill extends the
        window rather than reshuffling it. That matters because the optimizers
        are order-sensitive: a traversal that changed under refill would make an
        otherwise reproducible run depend on how many refills it happened to
        take.

        `keep` names candidates that must remain reachable whatever the window
        does, which is how a caller protects the panel it has already chosen.
        They are normally inside the frontier already, having been selected from
        it, and are added back when they are not.
        """
        ordered = self._provider._eligible_in_search_order()
        if len(self._examined) >= len(ordered):
            return False

        target = min(len(ordered), max(len(self._examined) * FRONTIER_GROWTH, 1))
        if target <= len(self._examined):
            target = min(len(ordered), len(self._examined) + 1)
        widened = list(ordered[:target])

        seen = set(widened)
        for sequence in keep:
            if sequence not in seen:
                widened.append(sequence)
                seen.add(sequence)
        self._examined = widened
        return True

    def exhaustion(self) -> str:
        """Which kind of "nothing left" this source is reporting."""
        return INVENTORY_EXHAUSTED if self.exhausted() else FRONTIER_EXHAUSTED

    def universe_size(self) -> int:
        return len(self._provider._eligible_in_search_order())

    def examined(self) -> int:
        return len(self._examined)

    def exhausted(self) -> bool:
        return self.examined() >= self.universe_size()

    def attach_positions(self, cache) -> None:
        """Give the provider the cache its position checks need.

        Mandatory before a batch is scored: a provider without one cannot tell a
        candidate with no binding sites from one that was never indexed, and
        that is the silent zero of Known Issues 5, 6 and 13.
        """
        self.position_cache = cache
        self._provider.attach_positions(cache)

    def ensure_positions(self, sequences: Sequence[str]) -> None:
        """Delegate to the provider, which owns the positions it hands out."""
        self._provider.ensure_positions(sequences)

    def describe(self) -> dict:
        universe = self.universe_size()
        return {
            "kind": self.kind,
            "universe": universe,
            "frontier": self._frontier_size,
            "examined": self.examined(),
            "unexamined": max(0, universe - self.examined()),
            "exhausted": self.exhausted(),
            "rank_key": "search_rank",
        }


def as_candidate_source(candidates) -> object:
    """Wrap a bare list, so every existing caller keeps working."""
    if hasattr(candidates, "initial") and hasattr(candidates, "describe"):
        return candidates
    return ListCandidateSource(candidates)


def open_candidate_source(
    data_dir,
    condition_id: str,
    lengths: Sequence[int],
    candidates: Optional[Sequence[str]] = None,
    frontier: Optional[int] = None,
):
    """The one way a command asks for candidates.

    An explicit list wins: a user who named a `--candidates` file has said which
    pool to design over, and quietly searching a different one would be worse
    than useless.

    Otherwise the inventory, when the directory has one. A directory written
    before the inventory existed has none, and a caller that supplied no list
    either has asked for something impossible rather than something empty.
    """
    if candidates is not None:
        return ListCandidateSource(candidates)

    path = Path(data_dir) / STAGE2_INVENTORY_NAME
    if not path.is_file():
        raise ValueError(
            f"No candidate inventory at {path} and no candidate list supplied. "
            "Run `neoswga filter` to build one, or pass an explicit list."
        )

    from neoswga.core.candidate_inventory import CandidateInventory
    from neoswga.core.candidate_provider import CandidateProvider

    inventory = CandidateInventory(path)
    provider = CandidateProvider(inventory, condition_id, lengths)
    if not provider._eligible_in_search_order():
        inventory.close()
        raise ValueError(
            f"The inventory at {path} has no eligible candidates of length(s) "
            f"{sorted(lengths)} under this reaction. Re-run `neoswga filter` "
            "under the chemistry you are designing for."
        )
    return InventoryCandidateSource(provider, frontier=frontier)
