"""Exact additive accounting for the specificity constraints.

A diagnostic, not part of the search. It lives here rather than in `neoswga.core`
because no command uses it, and a capability in the package that nothing reaches
is the defect the 2026-09-16 audit is named after.

What it is for. `occupancy.weighted_site_load` is a SUM of per-primer terms that
depend only on the primer, the references and the reaction, so a panel's
foreground and background loads are additive and the selectivity density floor
rearranges into a linear condition:

    density(panel) >= D   <=>   SUM_i (f_i * bg_len - D * b_i * fg_len) >= 0

Each candidate therefore carries a fixed `slack`. That makes two things exactly
computable rather than searchable: the highest density floor any N-primer panel
from a pool can satisfy, and whether a particular panel can still reach one.

Those are what showed the search leaves real headroom on the table. See
docs/validation/stage_one_constraint_awareness_2026-09-18.md.

`tier` and `completable` are kept because they are what three attempted Stage 1
rules were built from, and the note records why none of them shipped.
"""

from __future__ import annotations

import math


class SelectivityBudget:
    """Exact, additive feasibility for the specificity constraints.

    Stage 1's greedy set cover chose purely on marginal coverage, so on a row
    with a selectivity floor it spent candidates on coverage and the delivered
    panel missed the floor. Phase 5 measured the cost: the search reached a
    density of 60.11 where the best 12-primer panel from the same pool reaches
    79.807.

    No objective evaluation is needed to fix that, and no projection from a
    partial panel either. `occupancy.weighted_site_load` is a SUM of per-primer
    terms that depend only on the primer, the references and the reaction, so a
    panel's foreground and background loads are additive. The density floor then
    rearranges into a linear condition:

        density >= D  <=>  SUM_i (f_i * bg_len - D * b_i * fg_len) >= 0

    Each candidate therefore carries a fixed `slack`, positive when it helps the
    panel meet the floor and negative when it costs. Feasibility is the sign of
    a running sum, which is O(1) per candidate.

    `max_background_sites` is additive in the same way, except that
    `total_bg_sites` is the size of a UNION over primers and two primers can
    share a host position. Summing therefore bounds it from above, so a panel
    this admits certainly satisfies the real ceiling; the error is in the safe
    direction.

    `build` returns None when no specificity constraint is configured, so a
    design without one is steered exactly as it was before this existed.
    """

    def __init__(self, constraints, fg_length, bg_length, loads, background_sites):
        self.constraints = constraints
        self._fg_length = float(fg_length)
        self._bg_length = float(bg_length)
        self._loads = loads
        self._sites = background_sites
        # The best slack any k unselected candidates could still contribute,
        # as a prefix sum over the slacks sorted best first. Used as an
        # OPTIMISTIC bound, so it never forbids a choice that could still be
        # completed; see `completable`.
        slacks = sorted((self.slack(p) for p in loads), reverse=True)
        self._best_future = [0.0]
        for value in slacks:
            self._best_future.append(self._best_future[-1] + max(value, 0.0))

    def _reachable(self, slots):
        """The most slack `slots` further candidates could add."""
        if slots <= 0:
            return 0.0
        return self._best_future[min(slots, len(self._best_future) - 1)]

    def completable(self, panel, primer, slots_left):
        """Whether the floor is still REACHABLE after taking `primer`.

        The difference from `admits` is lookahead, and it is the difference
        between a useful rule and a harmful one. Requiring the panel to satisfy
        the floor at every intermediate step forbids ever taking a candidate
        with negative slack, however much coverage it adds and however much
        headroom later candidates would restore. Measured on the real pool that
        made the delivered density WORSE at strict floors, 24.97 against 42.62,
        because the greedy spent its headroom early on coverage and then could
        afford nothing good.

        So the test is not "is the panel feasible now" but "could it still end
        feasible": the running slack, plus this candidate's, plus the best the
        remaining slots could contribute. The bound is optimistic, computed over
        all candidates rather than only unselected ones, which means it never
        wrongly forbids and may admit a choice that turns out infeasible. That
        residue is what Stage 2 and the repair exist for, and an optimistic
        bound is the safe direction: a pessimistic one would silently shrink
        the panel.
        """
        if self.constraints.min_selectivity_density is None:
            return True
        running = sum(self.slack(p) for p in panel)
        return running + self.slack(primer) + self._reachable(slots_left - 1) >= 0

    @classmethod
    def build(cls, constraints, optimizer, candidates):
        """Per-candidate loads for this pool, or None when nothing constrains it.

        The loads cost one `mismatch_class_counts` per candidate, measured at
        about 50 microseconds for a 12-mer, so this is seconds for a shortlist
        and under a minute for the whole retained inventory. It is paid once per
        run rather than once per greedy step, which is what makes the approach
        affordable where scoring the full objective per candidate is not.
        """
        if constraints is None or not constraints.needs_background:
            return None
        # Everything the loads need, or nothing. A caller passing its own
        # optimizer-shaped object with a metrics table has no references and no
        # index, and a fabricated load would steer selection by a number
        # nobody measured -- the same reason `_effective_site_load` falls back
        # rather than inventing one.
        needed = ("bg_prefixes", "fg_prefixes", "fg_seq_lengths", "bg_seq_lengths")
        if optimizer is None or not all(getattr(optimizer, name, None) for name in needed):
            return None
        # `cache`, which is what `BaseOptimizer.__init__` assigns. Asking for
        # `position_cache` -- the constructor's PARAMETER name -- returned None
        # from every real optimizer and made this whole budget silently inert,
        # which the measurement caught and a unit test would not have.
        cache = getattr(optimizer, "cache", None)
        if cache is None:
            return None

        from neoswga.core.mismatch_counts import mismatch_class_counts
        from neoswga.core.occupancy import (
            default_mismatch_penalty,
            mismatch_tm,
            site_occupancy,
        )
        from neoswga.core.thermodynamics import calculate_enthalpy_entropy

        conditions = getattr(optimizer, "conditions", None)
        if conditions is None:
            # Without a temperature there is no occupancy, and a fabricated
            # load would steer the search by a number nobody measured.
            return None

        penalty = default_mismatch_penalty()
        depth = int(getattr(optimizer.config, "max_mismatches", 1) or 0)

        def load(primer, prefixes):
            dh, _ = calculate_enthalpy_entropy(primer)
            tm = conditions.calculate_effective_tm(primer)
            total = 0.0
            for distance, count in mismatch_class_counts(primer, prefixes, depth).items():
                if count:
                    total += count * site_occupancy(
                        dh, mismatch_tm(tm, distance, penalty), conditions.temp
                    )
            return total

        loads, sites = {}, {}
        try:
            for primer in dict.fromkeys(candidates):
                loads[primer] = (
                    load(primer, optimizer.fg_prefixes),
                    load(primer, optimizer.bg_prefixes),
                )
                sites[primer] = sum(
                    len(cache.get_positions(prefix, primer, "both"))
                    for prefix in optimizer.bg_prefixes
                )
        except (FileNotFoundError, OSError):
            # The same fallback `_effective_site_load` takes: without the count
            # tables the occupancy loads are not available, and guessing them
            # would steer selection by a different definition than acceptance.
            return None

        return cls(
            constraints=constraints,
            fg_length=float(sum(optimizer.fg_seq_lengths or [])) or 1.0,
            bg_length=float(sum(optimizer.bg_seq_lengths or [])) or 1.0,
            loads=loads,
            background_sites=sites,
        )

    def slack(self, primer):
        """This candidate's signed contribution to meeting the density floor."""
        floor = self.constraints.min_selectivity_density
        if floor is None:
            return 0.0
        fg_load, bg_load = self._loads.get(primer, (0.0, 0.0))
        return fg_load * self._bg_length - floor * bg_load * self._fg_length

    # How strongly a candidate serves the floor. Lower is preferred, and the
    # ordering is what the measurements support rather than a guess.
    #
    # AFFORDABLE means the candidate meets the floor on its own, so a panel
    # built from such candidates satisfies it however they combine. On the
    # measured pool, taking the best-coverage twelve from just those gave
    # density 75.6 at floor 60, 78.2 at 65 and 79.5 at 70, every one with
    # coverage above the 0.5 target -- against a search that delivered 60.11
    # and failed at 65.
    #
    # COMPLETABLE means it does not, but enough headroom remains that the panel
    # could still end feasible. That tier is needed because at strict floors
    # too few candidates are affordable to fill a panel: 8 at floor 75 and 6 at
    # 79, against a requested 12. Without it the greedy would run out.
    AFFORDABLE, COMPLETABLE, UNAFFORDABLE = 0, 1, 2

    def tier(self, panel, primer, slots_left):
        """Which of the three a candidate falls into. Lower is better."""
        if self.constraints.max_background_sites is not None:
            used = sum(self._sites.get(p, 0) for p in panel)
            if used + self._sites.get(primer, 0) > self.constraints.max_background_sites:
                return self.UNAFFORDABLE
        if self.constraints.min_selectivity_density is None:
            return self.AFFORDABLE
        if self.slack(primer) >= 0:
            return self.AFFORDABLE
        if self.completable(panel, primer, slots_left):
            return self.COMPLETABLE
        return self.UNAFFORDABLE

    def admits(self, panel, primer, slots_left=None):
        """Whether adding `primer` keeps every limit satisfiable.

        The density floor is checked with lookahead through `completable`, since
        it is not monotone: a later candidate can restore slack an earlier one
        spent. The site ceiling IS monotone -- sites only accumulate -- so no
        lookahead could help and it is checked directly.

        `slots_left=None` means "this is the whole panel", which checks the
        floor as it stands.
        """
        if self.constraints.min_selectivity_density is not None:
            if slots_left is None:
                running = sum(self.slack(p) for p in panel)
                if running + self.slack(primer) < 0:
                    return False
            elif not self.completable(panel, primer, slots_left):
                return False

        ceiling = self.constraints.max_background_sites
        if ceiling is not None:
            used = sum(self._sites.get(p, 0) for p in panel)
            if used + self._sites.get(primer, 0) > ceiling:
                return False
        return True
