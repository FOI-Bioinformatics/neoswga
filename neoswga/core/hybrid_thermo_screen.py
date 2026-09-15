"""Stage 0 of the hybrid optimizer: the thermodynamic screen and its cache.

Extracted from `hybrid_optimizer.py`, which was at its module-size budget. The
three methods here form one unit -- the criteria the screen runs under, the
screen itself, and the per-pool cache in front of it -- and they touch nothing
else in the optimizer beyond the configuration attributes named below.

`ThermoScreenMixin` expects its host to provide `polymerase`, `poly_config`,
`min_tm`, `max_tm`, `max_dimer_bp`, `conditions` and `_thermo_filter_cache`.
`HybridOptimizer.__init__` sets all of them.
"""

import logging
from typing import List

logger = logging.getLogger(__name__)


class ThermoScreenMixin:
    """The Stage-0 thermodynamic screen, mixed into `HybridOptimizer`."""

    def _thermo_criteria(self):
        """Criteria for the Stage-0 thermodynamic screen.

        Two things used to be wrong here and both were invisible.

        The Tm bounds came from the polymerase preset alone, so a window
        widened in params.json was honoured by the filter step and then
        narrowed back by the optimizer. An explicitly configured bound now
        wins; the preset remains the default when nothing is configured.

        `na_conc=50.0, mg_conc=0.0` were hardcoded. Every hairpin, homodimer
        and heterodimer energy this screen rejects primers on was therefore
        computed for a reaction containing no magnesium -- the same
        zero-magnesium fault schema v2 corrected in the presets, and Mg2+ is
        the dominant term in the salt correction. The configured buffer is
        used instead.
        """
        from neoswga.core.thermodynamic_filter import ThermodynamicCriteria

        conditions = getattr(self, "conditions", None)
        if conditions is not None:
            na_conc = conditions.na_conc
            mg_conc = conditions.mg_conc
        else:
            from neoswga.core.parameter import default_mg_conc

            na_conc = 50.0
            mg_conc = default_mg_conc(self.polymerase)

        return ThermodynamicCriteria(
            min_tm=self.min_tm if self.min_tm is not None else self.poly_config.min_primer_tm,
            max_tm=self.max_tm if self.max_tm is not None else self.poly_config.max_primer_tm,
            target_tm=self.poly_config.reaction_temp + 5,
            na_conc=na_conc,
            mg_conc=mg_conc,
            max_homodimer_dg=-10.0,
            max_heterodimer_dg=-10.0,
            max_hairpin_dg=-3.0,
            min_gc=self.poly_config.min_gc,
            max_gc=self.poly_config.max_gc,
            reaction_temp=self.poly_config.reaction_temp,
            polymerase=self.polymerase,
        )

    def _thermo_filter_candidates(self, candidates: List[str], verbose: bool = True) -> List[str]:
        """
        Apply thermodynamic filtering based on polymerase requirements.

        Filters candidates by Tm range and GC content appropriate for
        the configured polymerase.
        """
        if verbose:
            logger.info("\n" + "-" * 80)
            logger.info(
                f"PRE-STAGE: Thermodynamic Filtering ({self.polymerase}, "
                f"{self.poly_config.reaction_temp}C)"
            )
            logger.info("-" * 80)

        try:
            from neoswga.core.thermodynamic_filter import ThermodynamicCriteria, ThermodynamicFilter

            criteria = self._thermo_criteria()

            # The resolved conditions, so this screen judges a primer on the
            # same Tm the gate that admitted it used. Without them it applied
            # salt but no additive correction, and rejected primers that the
            # configured chemistry had brought into the window.
            thermo_filter = ThermodynamicFilter(
                criteria, conditions=getattr(self, "conditions", None)
            )
            # `check_heterodimers=False`: the pool-wide hub count is not a
            # property of a panel. It removed a primer conflicting with more
            # than a fraction of the WHOLE pool, while what matters is whether
            # it conflicts with the primers actually selected -- which the
            # greedy enforces pairwise against the panel it is building, and
            # which tests/test_delivered_panel_honours_the_dimer_limit.py pins.
            #
            # Counting across the pool also made the verdict depend on pool
            # size. Under candidate_retention='all_qc' the pool is an order of
            # magnitude larger, so a primer that survived in a 2,000-candidate
            # pool could be removed from a 20,000-candidate one with nothing
            # about the primer having changed.
            filtered, stats = thermo_filter.filter_candidates(
                candidates,
                check_heterodimers=False,
                max_dimer_bp=self.max_dimer_bp,
            )

            if verbose:
                logger.info(f"Filtered: {len(filtered)}/{len(candidates)} passed")
                if stats.get("mean_tm") is not None:
                    logger.info(f"  Mean Tm: {stats['mean_tm']:.1f}C")
                if stats.get("mean_gc") is not None:
                    logger.info(f"  Mean GC: {stats['mean_gc']:.1%}")

            if len(filtered) == 0:
                # Returning the input here turned a QC failure into a pass: a
                # configuration under which nothing survives thermodynamic
                # screening proceeded with every candidate just rejected, which
                # is not a screen. An empty result is now reported as one.
                logger.warning(
                    "No candidate passed thermodynamic screening under the "
                    "configured reaction (Tm window %.1f-%.1f C at %s). The "
                    "empty result stands; widen min_tm/max_tm, change the "
                    "additives, or supply a different candidate pool.",
                    criteria.min_tm,
                    criteria.max_tm,
                    self.polymerase,
                )

            return filtered

        except ImportError:
            logger.warning("Thermodynamic filter not available, skipping")
            return candidates

    def _thermo_filter_with_cache(self, candidates: List[str], verbose: bool = True) -> List[str]:
        """The Stage-0 screen, computed once per pool and reused for subsets.

        Returns the members of `candidates` that passed. A pool that is not a
        subset of the cached one recomputes, which is what a caller supplying a
        genuinely different candidate list should get.

        The reuse is an approximation, not an identity. The screen is not a
        pure function of the candidate sequences and the reaction conditions:
        `ThermodynamicFilter.filter_candidates` removes heterodimer hubs at
        `max_issues = int(len(passing) * max_heterodimer_fraction)`, counted
        over the pool it was handed, so both the threshold and the per-primer
        counts depend on pool composition. A subset has a smaller `max_issues`
        and can therefore remove a primer the full pool kept, which means
        reusing the full-pool verdict on a subset under-removes.

        What that costs: set 0 is always screened from the first pool and is
        unaffected, so the metrics and the summary describe a set the cache did
        not change. Alternative sets 1 to 4 are drawn from subsets and can
        differ from what a per-pool screen would give. On the pools measured
        here the hub removal is dormant either way -- `int(972 * 0.3)` is 291
        and no primer in these pools has 291 dimer partners -- but that is a
        property of these pools, not a guarantee from the code.
        """
        # Identity is the pool AND the reaction AND the QC limits. Keyed by
        # the pool alone, a verdict computed under one chemistry was reused
        # after the conditions changed.
        criteria = self._thermo_criteria()
        conditions = getattr(self, "conditions", None)
        identity = (
            frozenset(c.upper() for c in candidates),
            conditions.fingerprint() if conditions is not None else "no-conditions",
            (criteria.min_tm, criteria.max_tm, criteria.na_conc, criteria.mg_conc),
            getattr(self, "max_dimer_bp", None),
        )
        wanted = identity

        if self._thermo_filter_cache is not None:
            screened_pool, passed = self._thermo_filter_cache
            if screened_pool[1:] == identity[1:] and identity[0] <= screened_pool[0]:
                kept = [c for c in candidates if c.upper() in passed]
                if verbose:
                    logger.info(
                        "PRE-STAGE: reusing the thermodynamic screen computed over "
                        "%d candidates; %d of %d in this pool passed it",
                        len(screened_pool[0]),
                        len(kept),
                        len(candidates),
                    )
                return kept if kept else list(candidates)

        result = self._thermo_filter_candidates(candidates, verbose=verbose)
        self._thermo_filter_cache = (wanted, frozenset(c.upper() for c in result))
        return result
