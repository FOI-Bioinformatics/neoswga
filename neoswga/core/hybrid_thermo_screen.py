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

            thermo_filter = ThermodynamicFilter(criteria)
            filtered, stats = thermo_filter.filter_candidates(
                candidates,
                check_heterodimers=True,
                max_heterodimer_fraction=0.3,
                max_dimer_bp=self.max_dimer_bp,
            )

            if verbose:
                logger.info(f"Filtered: {len(filtered)}/{len(candidates)} passed")
                if stats.get("mean_tm") is not None:
                    logger.info(f"  Mean Tm: {stats['mean_tm']:.1f}C")
                if stats.get("mean_gc") is not None:
                    logger.info(f"  Mean GC: {stats['mean_gc']:.1%}")

            if len(filtered) == 0:
                logger.warning(
                    "No primers passed thermodynamic filtering, " "using unfiltered candidates"
                )
                return candidates

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
        wanted = frozenset(c.upper() for c in candidates)

        if self._thermo_filter_cache is not None:
            screened_pool, passed = self._thermo_filter_cache
            if wanted <= screened_pool:
                kept = [c for c in candidates if c.upper() in passed]
                if verbose:
                    logger.info(
                        "PRE-STAGE: reusing the thermodynamic screen computed over "
                        "%d candidates; %d of %d in this pool passed it",
                        len(screened_pool),
                        len(kept),
                        len(candidates),
                    )
                return kept if kept else list(candidates)

        result = self._thermo_filter_candidates(candidates, verbose=verbose)
        self._thermo_filter_cache = (wanted, frozenset(c.upper() for c in result))
        return result
