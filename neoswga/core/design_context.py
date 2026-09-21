"""One resolution of a params file, shared by every command that designs.

Finding F7 of the 2026-09-16 pipeline audit, the configuration half.
`plan-pool` resolved the reaction, the coverage reach and the dimer limits from
params.json. `expand-primers` resolved none of them: it ran at a hard-coded
3 kb reach with `conditions=None`, `tm_weight=0.0`, `dimer_penalty=0.0` and
`max_dimer_bp=None`.

So one params file produced designs under different chemistry depending on
which command was run, and the command that exists to ADD to an existing panel
was the one running without that panel's chemistry. A primer chosen at 3 kb
against a design made at 8 kb is not chosen for the same genome geometry, and
one chosen with no additive correction is not chosen for the same reaction.

The resolution lives here once. Two resolutions that merely look similar are
how they drift, so the test that matters is the one asserting both commands
agree, not the one asserting each is individually reasonable.

This carries CONFIGURATION, not state: no cache, no candidates, no genome
lengths beyond what the reach needs. A context is safe to build early and to
pass anywhere.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from types import SimpleNamespace
from typing import Any, Mapping, Optional

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class DesignContext:
    """The chemistry, reach and dimer limits one params file resolves to."""

    conditions: Any
    coverage_reach: int
    max_dimer_bp: int
    max_self_dimer_bp: int
    min_tm: float
    max_tm: float
    fg_circular: bool
    polymerase: str
    max_dimer_dg: Optional[float] = None
    constraints: Any = None
    optimizer_settings: tuple[tuple[str, Any], ...] = ()

    def optimizer_config(self, **overrides):
        """An `OptimizerConfig` carrying exactly these values.

        Built here rather than at each call site, because the fields a config
        needs and the fields a params file resolves to are the same fields, and
        listing them twice is how `max_dimer_bp` came to be 3 in one place and
        4 in another.
        """
        from neoswga.core.base_optimizer import OptimizerConfig

        settings = dict(self.optimizer_settings)
        settings.update(
            max_dimer_bp=self.max_dimer_bp,
            max_self_dimer_bp=self.max_self_dimer_bp,
            max_dimer_dg=self.max_dimer_dg,
            min_tm=self.min_tm,
            max_tm=self.max_tm,
            extension_reach=self.coverage_reach,
            fg_circular=self.fg_circular,
        )
        settings.update(overrides)
        return OptimizerConfig(**settings)


def design_context_from_params(
    params: Mapping[str, Any],
    coverage_reach_override: Optional[int] = None,
) -> DesignContext:
    """Resolve a params mapping into the settings every design command needs.

    Precedence for the reach is the explicit override, then `coverage_reach` in
    the file, then the polymerase's realistic per-primer reach. That is the
    order `plan-pool` already used; lifting it here is what lets
    `expand-primers` use the same one instead of 3 kb.
    """
    from dataclasses import fields

    from neoswga.core.base_optimizer import OptimizerConfig
    from neoswga.core.coverage import resolve_coverage_reach
    from neoswga.core.reaction_conditions import build_reaction_conditions

    polymerase = params.get("polymerase", "phi29") or "phi29"
    conditions = build_reaction_conditions(SimpleNamespace(**dict(params)))
    # `override if override is not None else ...` rather than `or`: an explicit
    # 0 is falsy, so the old form silently replaced it with the polymerase
    # default and reported coverage at 3 kb for a request that said otherwise.
    # `resolve_coverage_reach` refuses 0, which is the answer the user needs.
    override = (
        coverage_reach_override
        if coverage_reach_override is not None
        else params.get("coverage_reach")
    )
    reach = resolve_coverage_reach(polymerase, override=override)
    from .panel_acceptance import constraints_from_parameter

    return DesignContext(
        constraints=constraints_from_parameter(SimpleNamespace(**dict(params))),
        conditions=conditions,
        coverage_reach=int(reach),
        max_dimer_bp=int(params.get("max_dimer_bp", 3)),
        max_self_dimer_bp=int(params.get("max_self_dimer_bp", 4)),
        min_tm=float(params.get("min_tm", 20)),
        max_tm=float(params.get("max_tm", 50)),
        fg_circular=bool(params.get("fg_circular", False)),
        polymerase=str(polymerase),
        max_dimer_dg=params.get("max_dimer_dg"),
        optimizer_settings=tuple(
            (field.name, params[field.name])
            for field in fields(OptimizerConfig)
            if field.name in params and params[field.name] is not None
        ),
    )
