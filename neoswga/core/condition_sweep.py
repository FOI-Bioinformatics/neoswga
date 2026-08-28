"""Choosing reaction conditions for specificity, from a measurement.

Two things in this codebase claimed to optimise conditions for specificity and
neither looked at a genome.

``additive_optimizer._evaluate_combination`` scored its ``"specificity"``
objective as ``effective_binding_rate * (1 / koff_factor)`` of a *synthetic*
primer assembled from a length and a GC fraction. No primer set, no background,
nothing that could distinguish a specific design from an indiscriminate one --
it measured binding stability, and a primer that binds everything scores well on
that. The ``optimize-conditions`` CLI did not use that optimiser at all; it was
a GC-content lookup table taking ``--fg`` with no ``--bg``, so it could not have
measured specificity even in principle.

With occupancy-weighted selectivity available the question is answerable. Sweep
a temperature x additive grid, compute selectivity against the actual background
for the actual primer set, and report it beside the amplification it costs.

Both halves are necessary. Stringency buys discrimination by pushing duplexes
off their melting transition, and that is the same motion as pushing them toward
not being bound at all. A search on selectivity alone walks to the hottest,
most-additive corner of any grid and recommends a reaction that amplifies
nothing. The optimum is interior, and it can only be found by weighing the two.

The useful window is near Tm. Far below it both perfect and mismatched sites are
saturated and there is nothing to discriminate; far above, neither is bound.
A set whose primers melt far apart has no single temperature where all of them
sit in that window, which is why Tm tightness is a precondition for conditions
being able to help at all.
"""

import logging
from dataclasses import dataclass
from typing import List, Optional, Sequence

logger = logging.getLogger(__name__)

# Default grids. Deliberately coarse: the point is to show the shape of the
# trade-off and locate the region, not to pretend to a precision the underlying
# additive coefficients do not have (see docs/SCIENCE_CITATIONS.md, where the
# enzyme-pathway magnitudes are marked EMPIRICAL).
DEFAULT_TEMPERATURES = (28.0, 30.0, 32.0, 35.0, 38.0, 42.0)
DEFAULT_DMSO = (0.0, 2.5, 5.0, 7.5)
DEFAULT_BETAINE = (0.0, 0.5, 1.0)

# How much the amplification cost weighs against the specificity gain when
# ranking points. Equal footing: both enter as ratios to the grid's best, so
# neither unit dominates, and a point has to be good at both to win. Exposed
# rather than buried because it encodes a preference, not a fact.
DEFAULT_AMPLIFICATION_WEIGHT = 0.5


@dataclass(frozen=True)
class ConditionPoint:
    """One point on the condition grid, with both halves of the trade."""

    temp: float
    dmso_percent: float
    betaine_m: float
    selectivity: float
    amplification: float
    effective_fg_sites: float
    effective_bg_sites: float

    @property
    def label(self) -> str:
        parts = [f"{self.temp:.0f}C"]
        if self.dmso_percent:
            parts.append(f"{self.dmso_percent:g}% DMSO")
        if self.betaine_m:
            parts.append(f"{self.betaine_m:g}M betaine")
        return ", ".join(parts)


def sweep_conditions(
    primers: Sequence[str],
    fg_prefixes: Sequence[str],
    bg_prefixes: Sequence[str],
    polymerase: str = "phi29",
    temperatures: Optional[Sequence[float]] = None,
    dmso_percentages: Optional[Sequence[float]] = None,
    betaine_concentrations: Optional[Sequence[float]] = None,
    max_mismatches: int = 1,
) -> List[ConditionPoint]:
    """Selectivity and amplification across a temperature x additive grid.

    Selectivity uses the same `occupancy.weighted_site_load` the optimizer's
    metrics and `evaluate-set` use, so a sweep and a scored design cannot
    disagree about what a given condition is worth.
    """
    from neoswga.core.base_optimizer import _selectivity_from_loads
    from neoswga.core.mechanistic_model import MechanisticModel
    from neoswga.core.occupancy import weighted_site_load
    from neoswga.core.reaction_conditions import ReactionConditions

    temperatures = list(temperatures if temperatures is not None else DEFAULT_TEMPERATURES)
    dmso_grid = list(dmso_percentages if dmso_percentages is not None else DEFAULT_DMSO)
    betaine_grid = list(betaine_concentrations if betaine_concentrations is not None else (0.0,))

    points: List[ConditionPoint] = []
    skipped: List[tuple] = []
    for temp in temperatures:
        for dmso in dmso_grid:
            for betaine in betaine_grid:
                try:
                    conditions = ReactionConditions(
                        temp=temp,
                        polymerase=polymerase,
                        dmso_percent=dmso,
                        betaine_m=betaine,
                    )
                except ValueError as exc:
                    # An enzyme has an operating range and the grid may reach
                    # past it -- phi29 stops at 40 C. Skipping is right, but a
                    # silently shorter grid would let a user read an optimum at
                    # the edge as a real one when the search simply stopped
                    # there.
                    skipped.append((temp, str(exc)))
                    continue
                fg_load = weighted_site_load(primers, fg_prefixes, conditions, max_mismatches)
                bg_load = weighted_site_load(primers, bg_prefixes, conditions, max_mismatches)

                # Amplification from the mechanistic model, averaged over the
                # set. It carries the enzyme's own response to temperature and
                # additives -- processivity, speed, stability -- which is the
                # cost that stringency incurs and the reason an optimum exists.
                model = MechanisticModel(conditions)
                amplification = sum(
                    model.calculate_effects(p, template_gc=0.5).predicted_amplification_factor
                    for p in primers
                ) / max(len(primers), 1)

                points.append(
                    ConditionPoint(
                        temp=temp,
                        dmso_percent=dmso,
                        betaine_m=betaine,
                        selectivity=_selectivity_from_loads(fg_load, bg_load),
                        amplification=amplification,
                        effective_fg_sites=fg_load,
                        effective_bg_sites=bg_load,
                    )
                )

    if skipped:
        outside = sorted({temp for temp, _ in skipped})
        logger.warning(
            f"Skipped {len(skipped)} grid point(s) outside {polymerase}'s operating "
            f"range at {outside} C. The grid searched is narrower than the one "
            f"requested, so an optimum at its edge may be the edge rather than an "
            f"optimum."
        )
    if not points:
        raise ValueError(
            f"No grid point was valid for {polymerase}. Every requested "
            f"temperature lies outside its operating range."
        )
    return points


def score_point(point: ConditionPoint, points: Sequence[ConditionPoint], weight: float) -> float:
    """Rank a point against the grid it came from.

    Both terms enter as a fraction of the grid's best, so a selectivity of 12
    and an amplification of 0.01 can be compared without one unit swamping the
    other. `weight` is how much the amplification cost counts; at 0 this
    degenerates to ranking on selectivity alone, which walks to the hottest
    corner of any grid and recommends a reaction that amplifies nothing.
    """
    best_selectivity = max((p.selectivity for p in points), default=0.0) or 1.0
    best_amplification = max((p.amplification for p in points), default=0.0) or 1.0

    return (1.0 - weight) * (point.selectivity / best_selectivity) + weight * (
        point.amplification / best_amplification
    )


def recommend_conditions_for_set(
    primers: Sequence[str],
    fg_prefixes: Sequence[str],
    bg_prefixes: Sequence[str],
    polymerase: str = "phi29",
    amplification_weight: float = DEFAULT_AMPLIFICATION_WEIGHT,
    **grid,
) -> ConditionPoint:
    """The best point on the grid, weighing specificity against amplification."""
    points = sweep_conditions(
        primers=primers,
        fg_prefixes=fg_prefixes,
        bg_prefixes=bg_prefixes,
        polymerase=polymerase,
        **grid,
    )
    if not points:
        raise ValueError("Condition grid was empty; nothing to recommend.")

    return max(points, key=lambda p: score_point(p, points, amplification_weight))


def format_sweep_table(
    points: Sequence[ConditionPoint], amplification_weight: float = DEFAULT_AMPLIFICATION_WEIGHT
) -> str:
    """Render the grid so the trade-off is visible rather than summarised.

    A single recommended condition hides whether the optimum was sharp or flat,
    and that matters: a flat optimum means the chemistry is not the binding
    constraint on this design and effort belongs elsewhere.
    """
    from neoswga.core.base_optimizer import MAX_SELECTIVITY

    lines = [
        f"{'conditions':<28}{'selectivity':>12}{'amplification':>15}{'score':>9}",
        "-" * 64,
    ]
    ranked = sorted(
        points, key=lambda p: score_point(p, points, amplification_weight), reverse=True
    )
    for point in ranked:
        # MAX_SELECTIVITY stands in for "unbounded", not a measurement, and
        # printing 1000000.00 in a column of real ratios reads as one.
        shown = (
            "no bg binding" if point.selectivity >= MAX_SELECTIVITY else f"{point.selectivity:.2f}"
        )
        lines.append(
            f"{point.label:<28}{shown:>12}"
            f"{point.amplification:>15.4f}"
            f"{score_point(point, points, amplification_weight):>9.3f}"
        )

    if all(p.selectivity >= MAX_SELECTIVITY for p in points):
        # Every point saturating means the ranking above is decided entirely by
        # amplification. That is a legitimate answer -- there is no background
        # binding to discriminate against -- but a reader comparing conditions
        # would otherwise take the ordering as a specificity result.
        lines.append("")
        lines.append(
            "No background binding was detected at any condition, so selectivity "
            "cannot separate them and this ranking reflects amplification alone. "
            "A larger or closer background is needed for conditions to be a "
            "specificity lever here."
        )
    return "\n".join(lines)
