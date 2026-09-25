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

import logging
import math
from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import Any

from neoswga.core.candidate_inventory import STAGE2_INVENTORY_NAME
from neoswga.core.exceptions import ReferenceDataError
from neoswga.core.position_cache import MissingPositionsError

logger = logging.getLogger(__name__)

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


class CandidateFrontier(list):
    """Initial candidates with their source retained for later refills."""

    def __init__(self, source):
        super().__init__(source.initial())
        self.source = source


def describe_reach(candidates) -> dict[str, Any] | None:
    """How much of the available pool a run could reach, or None if unknown.

    A search that qualifies on its opening frontier stops there, so a design
    over a 2,000-candidate shortlist may never look at the other 363,073 the
    inventory holds. Both numbers are known at the moment the pool is opened
    and neither reached any artifact, so a reader had no way to tell a search
    over everything from a search over half a per cent of it.

    `examined` is what the run actually looked at, `universe` is what it could
    have. `complete` is the honest headline: it is true when they agree, which
    is also the ordinary case for a `--candidates` list, where the user named
    the pool and there is nothing behind it.

    Returns None rather than a guess when the caller passed a plain list. That
    is absence, not a reach of zero, and the two must not render alike.
    """
    source = getattr(candidates, "source", None)
    if source is None:
        return None
    try:
        described = source.describe()
    except Exception:  # pragma: no cover - a description must not fail a run
        return None

    universe = int(described.get("universe") or 0)
    examined = int(described.get("examined") or described.get("frontier") or 0)
    return {
        "kind": described.get("kind") or getattr(source, "kind", "unknown"),
        "universe": universe,
        "examined": examined,
        "complete": bool(universe and examined >= universe),
        "fraction": round(examined / universe, 6) if universe else None,
    }


class ListCandidateSource:
    """A fixed list of candidates, exhausted as soon as it is handed over.

    The universe and the frontier are the same set here, which is exactly what
    a `--candidates` file means: the user named the pool, and there is nothing
    else to reach for.
    """

    kind = "list"

    def __init__(self, sequences: Sequence[str]):
        self._sequences = list(dict.fromkeys(str(s).upper() for s in sequences))
        self._examined: list[str] = []

    def initial(self, limit: int | None = None) -> list[str]:
        """The starting frontier."""
        self._examined = self._sequences if limit is None else self._sequences[: max(0, limit)]
        return list(self._examined)

    def frontier(self) -> list[str]:
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

    def __init__(self, provider, frontier: int | None = None):
        self._provider = provider
        self._frontier_size = frontier
        self._examined: list[str] = []

    def initial(self, limit: int | None = None) -> list[str]:
        size = limit if limit is not None else self._frontier_size
        ordered = self._provider._eligible_in_search_order()
        self._examined = list(ordered if size is None else ordered[: max(0, size)])
        return list(self._examined)

    def frontier(self) -> list[str]:
        return list(self._examined)

    def universe(self) -> list[str]:
        """Every eligible candidate, including the ones past the frontier.

        `universe_size` already reports the length. A caller asking whether a
        delivered primer was a legitimate choice needs the sequences instead:
        `step3_df.csv` holds only the `max_primer` shortlist, so it is a strict
        subset of the correct answers whenever an inventory is present.
        """
        return list(self._provider._eligible_in_search_order())

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
    candidates: Sequence[str] | None = None,
    frontier: int | None = None,
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


def open_explicit_source(candidates):
    """The pool the user named, searched exactly as given.

    A `--candidates` file is a statement about which pool to design over.
    Searching a different one is worse than useless, and it has already cost
    this project a whole experiment: the occupancy-gate measurement recorded in
    Known Issue 17 appeared to improve selectivity density to 44.28, and did
    not. `open_source_or_list` had preferred the inventory over the capped list
    that was handed to it, so the run searched the inventory as usual and the
    apparent gain was the smaller frontier. None of the delivered primers were
    in the capped pool.
    """
    pool = list(candidates)
    logger.info("Searching the %d candidates named on the command line.", len(pool))
    return ListCandidateSource(pool)


def open_inventory_source(data_dir, condition_id, lengths, frontier=None):
    """The inventory for this directory, or None when it does not have one.

    The return value and the exception mean different things, and conflating
    them is what this function was split out to stop.

    `None` is a fact about the directory: it holds no inventory, because it was
    written before inventories existed or because `filter` has not run here. A
    caller with a candidate list of its own may use that list, and the choice
    is the caller's to declare.

    A `ReferenceDataError` means the directory HAS an inventory that cannot
    answer this request, which is almost always a reaction fingerprint that no
    candidate was assessed under. Falling back to a list there would design
    over a pool chosen by an error rather than by anyone, and the candidates
    the inventory holds would be silently unreachable.
    """
    path = Path(data_dir) / STAGE2_INVENTORY_NAME
    if not path.is_file():
        return None

    from neoswga.core.candidate_inventory import CandidateInventory
    from neoswga.core.candidate_provider import CandidateProvider

    inventory = CandidateInventory(path)
    provider = CandidateProvider(inventory, condition_id, lengths)
    if not provider._eligible_in_search_order():
        inventory.close()
        raise ReferenceDataError(
            f"candidate inventory {path}",
            f"holds no eligible candidate of length(s) {sorted(lengths)} under "
            f"reaction fingerprint {condition_id}",
            "Re-run `neoswga filter` under the chemistry you are designing for.",
        )
    return InventoryCandidateSource(provider, frontier=frontier)


def open_design_source(data_dir, condition_id, lengths, fallback, explicit=False):
    """Which pool this command searches, decided once and said out loud.

    Replaces `open_source_or_list` on 2026-09-21. The rule is the same in the
    ordinary case -- the inventory when the directory has one, because the
    inventory holds every candidate that cleared hard QC and a `max_primer`
    shortlist hides the rest (audit finding F1) -- but the two ways it used to
    be reached are now distinguishable.

    What changed is the error path. Opening the inventory used to be wrapped in
    `except ValueError`, so a directory with an unusable inventory silently
    became a run over the caller's CSV, at the caller's smaller frontier, with
    one `logger.info` line to say so. That is a pool chosen by an exception.

    The frontier still opens at the size of `fallback`, so no delivered panel
    moves; `plan-pool` is the only caller that reaches past it.
    """
    pool = list(fallback)
    if explicit:
        return open_explicit_source(pool)

    source = open_inventory_source(data_dir, condition_id, lengths, frontier=len(pool) or None)
    if source is None:
        logger.info(
            "No candidate inventory in %s; searching the candidate list supplied "
            "(%d candidates).",
            data_dir,
            len(pool),
        )
        return ListCandidateSource(pool)

    described = source.describe()
    frontier = described["frontier"] or described["universe"]
    logger.info(
        "Candidate universe: %d eligible in the inventory; this run will search a "
        "frontier of %d, leaving %d unexamined.",
        described["universe"],
        frontier,
        max(0, described["universe"] - frontier),
    )
    return source


def order_candidates_by_background(
    cache,
    candidates,
    fg_prefixes,
    bg_prefixes,
    min_ratio=1.0,
    verbose=False,
):
    """Put the candidates that bind the host least in front, delete none.

    Replaces `_prefilter_by_background`, which deleted, on 2026-09-17. Two
    measurements retired it.

    `--min-fg-bg-ratio` was inert. The old function kept every candidate at or
    above the ratio and then, if that removed more than `max_removal_fraction`
    of them, discarded the threshold and kept the top 80% by ratio instead. On
    the real Wolbachia shortlist the threshold removes 64.8% at its default of
    1.0, so the clause fired and exactly 400 of 2,000 went; it fired identically
    at 2.0 and at 5.0. The flag changed nothing above about 1.0, and the rule
    actually in force was "always drop the worst 20%".

    And which candidates survived depended on how many others happened to be in
    the same batch, because a bound on the fraction of a batch is not a rule
    about candidates.

    Ordering has neither problem. After Phase 4 increment 5 a candidate at the
    back of the scan is still reachable, because the frontier refills; deleting
    one is final, which is the shape of Known Issue 9, where an evenness gate
    removed primers the optimizer had already selected.

    The partition is STABLE, so the order within each group is the order the
    candidates arrived in. That order is the inventory's `search_rank`, which
    increment 5 established as the traversal, and the optimizers are
    order-sensitive; a full re-sort by ratio would discard it.

    Returns `(ordered, rejected)`, where `rejected` maps each deprioritised
    candidate to the reason. Nothing is removed from `ordered`.
    """
    if not math.isfinite(min_ratio) or min_ratio < 0:
        raise ValueError(f"min_ratio must be a finite non-negative number, not {min_ratio!r}")

    pool = list(candidates)
    if not pool or not bg_prefixes:
        # With no background there is no ratio to order on, and inventing one
        # would be worse than leaving the caller's order alone.
        return pool, {}

    ahead, behind, rejected = [], [], {}
    for primer in pool:
        fg_count = sum(len(cache.get_positions(p, primer, "both")) for p in fg_prefixes)
        bg_count = sum(len(cache.get_positions(p, primer, "both")) for p in bg_prefixes)
        ratio = fg_count / (bg_count + 1)
        if ratio >= min_ratio:
            ahead.append(primer)
        else:
            behind.append(primer)
            rejected[primer] = (
                f"fg/bg binding ratio {ratio:.3f} is below the {min_ratio} threshold; "
                f"searched last rather than removed"
            )

    if rejected and verbose:
        logger.info(
            "  Background ordering: %d of %d candidates fall below the %s fg/bg "
            "ratio and are searched last. None are removed.",
            len(rejected),
            len(pool),
            min_ratio,
        )
    return ahead + behind, rejected
