"""Single source of truth for polymerase characteristics.

Before this module, polymerase facts were spread across ~20 hardcoded tables in
``reaction_conditions``, ``mechanistic_params``, ``hybrid_optimizer``,
``param_validator``, ``parameter``, ``condition_suggester``,
``additive_optimizer``, ``results_interpreter``, the JSON schema and five argparse
``choices=`` lists. They had already drifted.

Two design rules make this a faithful refactor rather than a silent science change:

1. **Distinct concepts get distinct fields, seeded verbatim.** Several tables
   disagreed *because they meant different things*, not because one was stale.
   Three temperature bands (vendor-optimal / hard-reject / warn), two reach
   quantities (per-primer coverage reach vs single-molecule processivity), and
   three primer-length notions (design range / scalar default / additive-optimizer
   base) are all kept apart. Collapsing them would re-rank primer sets.

2. **Genuine contradictions are recorded, not quietly fixed.** Where the old tables
   are mutually impossible, the values are still seeded as-is and the conflict is
   asserted in ``tests/test_registry_invariants.py`` as a strict xfail. See
   ``INCONSISTENCIES.md``. Fixing them changes optimizer output and belongs in a
   separate change with integration baselines regenerated deliberately.

This module is deliberately stdlib-only: ``param_validator``, the argparse
``choices=`` lists and the schema renderer all import it, and none of them should
have to pull in numpy/scipy just to learn the polymerase names.
"""

from dataclasses import dataclass, field
from typing import Dict, Mapping, Optional, Tuple

__all__ = [
    "PolymeraseSpec",
    "POLYMERASES",
    "get_polymerase",
    "polymerase_names",
    "resolve_polymerase_name",
]


@dataclass(frozen=True)
class PolymeraseSpec:
    """Everything the codebase knows about one polymerase.

    Field groups are documented where two similar-looking numbers are genuinely
    different quantities -- that distinction is the whole point of this dataclass.
    """

    key: str
    display_name: str
    description: str

    # -- Thermal -----------------------------------------------------------
    # THREE different bands. They are not redundant:
    #   optimal_temp        single best operating temperature
    #   temp_optimal_range  vendor-recommended band (POLYMERASE_CHARACTERISTICS)
    #   temp_hard_range     outside this, ReactionConditions raises (_validate)
    #   temp_warn_range     outside this, param_validator warns
    optimal_temp: float
    temp_optimal_range: Tuple[float, float]
    temp_hard_range: Tuple[float, float]
    temp_warn_range: Tuple[float, float]

    # -- Reach -------------------------------------------------------------
    # processivity_bp: single-molecule maximum, used for amplification-network
    #   reachability.
    # typical_amplicon_bp: realistic per-primer reach in a dense SWGA design,
    #   used for fg_coverage scoring. This is the smaller, coverage-relevant one.
    # legacy_hybrid_max_extension: what hybrid_optimizer's set-cover currently
    #   uses. Seeded verbatim; for bst it exceeds processivity_bp, which is
    #   physically impossible. See INCONSISTENCIES.md.
    processivity_bp: int
    typical_amplicon_bp: int
    legacy_hybrid_max_extension: int
    extension_rate_nt_s: float
    processivity_step: float

    # -- Chemistry ---------------------------------------------------------
    strand_displacement: bool
    exonuclease: str  # "3to5" | "none"
    error_rate: float
    mg_default_mm: float

    # -- Primer design guidance -------------------------------------------
    # THREE different notions, all previously separate tables:
    #   primer_length_range      design k-mer window (min_k, max_k)
    #   default_primer_length    single scalar default (cli/commands.py)
    #   additive_opt_base_length additive-optimizer base length
    primer_length_range: Tuple[int, int]
    default_primer_length: int
    additive_opt_base_length: int
    primer_tm_range: Tuple[float, float]
    gc_range: Tuple[float, float]

    # -- Optimizer / preset knobs -----------------------------------------
    # preset_reaction_temp is the operating temperature hybrid_optimizer and
    # additive_optimizer use. It equals optimal_temp for every polymerase except
    # bst (60.0 vs 63.0) -- another seeded-verbatim divergence.
    preset_reaction_temp: float
    thermo_filter: bool
    primer_multiplier: float

    # A distributive enzyme dissociates after a short run and re-binds, so its
    # effective reach over a long incubation exceeds its single-event
    # processivity. For those, typical_amplicon_bp > processivity_bp is correct
    # rather than a contradiction, and the invariant tests treat them separately.
    distributive: bool = False

    additive_limits: Mapping[str, float] = field(default_factory=dict)
    dmso_response: Mapping[str, float] = field(default_factory=dict)

    aliases: Tuple[str, ...] = ()
    references: Tuple[str, ...] = ()

    @property
    def temp_range(self) -> Tuple[float, float]:
        """Back-compat alias for the vendor-optimal band."""
        return self.temp_optimal_range


POLYMERASES: Dict[str, PolymeraseSpec] = {
    "phi29": PolymeraseSpec(
        key="phi29",
        display_name="Phi29 DNA Polymerase",
        description="High-fidelity, high-processivity polymerase for MDA/SWGA",
        optimal_temp=30.0,
        temp_optimal_range=(30.0, 40.0),
        temp_hard_range=(20.0, 40.0),
        temp_warn_range=(20.0, 40.0),
        processivity_bp=70000,
        typical_amplicon_bp=3000,
        legacy_hybrid_max_extension=70000,
        # Blanco et al. (1989) -- the same paper cited for the 70 kb
        # processivity -- reports ~53 nt/s; single-molecule work agrees at
        # ~50 nt/s. The previous 150 nt/s compressed every simulated timing ~3x.
        extension_rate_nt_s=53,
        processivity_step=0.99857,
        strand_displacement=True,
        exonuclease="3to5",
        error_rate=1e-6,
        mg_default_mm=10.0,
        primer_length_range=(6, 12),
        default_primer_length=9,
        additive_opt_base_length=12,
        primer_tm_range=(20.0, 50.0),
        gc_range=(0.25, 0.75),
        preset_reaction_temp=30.0,
        thermo_filter=False,
        primer_multiplier=1.0,
        additive_limits={"dmso_percent": 8.0, "betaine_m": 2.5},
        dmso_response={"threshold": 5.0, "mild_coef": 0.02, "steep_coef": 0.12},
        references=("Blanco et al. (1989) JBC 264:8935",),
    ),
    "equiphi29": PolymeraseSpec(
        key="equiphi29",
        display_name="EquiPhi29 DNA Polymerase",
        description="Thermostable phi29 variant for higher specificity",
        optimal_temp=42.0,
        temp_optimal_range=(42.0, 45.0),
        temp_hard_range=(30.0, 50.0),
        temp_warn_range=(40.0, 47.0),
        processivity_bp=80000,
        typical_amplicon_bp=4000,
        legacy_hybrid_max_extension=80000,
        extension_rate_nt_s=200,
        processivity_step=0.99875,
        strand_displacement=True,
        exonuclease="3to5",
        error_rate=1e-6,
        mg_default_mm=10.0,
        primer_length_range=(10, 18),
        default_primer_length=15,
        additive_opt_base_length=15,
        primer_tm_range=(37.0, 62.0),
        gc_range=(0.25, 0.75),
        preset_reaction_temp=42.0,
        thermo_filter=True,
        primer_multiplier=0.85,
        additive_limits={"dmso_percent": 6.0, "betaine_m": 2.0},
        dmso_response={"threshold": 4.0, "mild_coef": 0.025, "steep_coef": 0.15},
        references=("Thermo Scientific EquiPhi29 manual MAN0030288",),
    ),
    "bst": PolymeraseSpec(
        key="bst",
        display_name="Bst 2.0/3.0 DNA Polymerase",
        description="LAMP-compatible, high-temperature isothermal amplification",
        optimal_temp=63.0,
        temp_optimal_range=(60.0, 65.0),
        temp_hard_range=(50.0, 72.0),
        temp_warn_range=(55.0, 70.0),
        processivity_bp=2000,
        typical_amplicon_bp=1000,
        # NOTE: exceeds processivity_bp. Seeded verbatim; see INCONSISTENCIES.md.
        legacy_hybrid_max_extension=10000,
        extension_rate_nt_s=100,
        processivity_step=0.9512,
        strand_displacement=True,
        exonuclease="none",
        error_rate=1e-4,
        mg_default_mm=8.0,
        primer_length_range=(15, 25),
        default_primer_length=20,
        additive_opt_base_length=18,
        primer_tm_range=(50.0, 75.0),
        gc_range=(0.25, 0.75),
        # NOTE: 60.0, not optimal_temp 63.0. Seeded verbatim.
        preset_reaction_temp=60.0,
        thermo_filter=True,
        primer_multiplier=1.0,
        additive_limits={"dmso_percent": 5.0, "betaine_m": 1.5},
        dmso_response={"threshold": 3.0, "mild_coef": 0.03, "steep_coef": 0.18},
        references=("Notomi et al. (2000) NAR 28:e63",),
    ),
    "klenow": PolymeraseSpec(
        key="klenow",
        display_name="Klenow Fragment (exo-)",
        description="Budget-friendly alternative with moderate processivity",
        optimal_temp=37.0,
        temp_optimal_range=(25.0, 40.0),
        temp_hard_range=(25.0, 40.0),
        # NOTE: wider than temp_hard_range, so 41 C passes validation and then
        # raises. Seeded verbatim; see INCONSISTENCIES.md.
        temp_warn_range=(20.0, 42.0),
        # Klenow is DISTRIBUTIVE: it dissociates after ~5-40 nt. 42 nt was
        # measured directly by single-molecule nanocircuit recording (JACS 2013,
        # PMC3738269). The previous 10000 bp value cited Bambara et al. (1978),
        # which does not support it, and overstated reach by ~250x.
        processivity_bp=40,
        # Effective per-primer reach is NOT bounded by single-event processivity
        # for a distributive enzyme -- the polymerase re-binds and continues, so
        # over a long isothermal incubation synthesis extends well past 40 nt.
        # 500 bp is a conservative estimate for a multi-primer reaction, not a
        # measurement; it is deliberately far below phi29's 3000 bp because low
        # processivity plus moderate strand displacement is exactly why Klenow
        # cannot carry a sparse primer set.
        typical_amplicon_bp=500,
        distributive=True,
        legacy_hybrid_max_extension=5000,
        extension_rate_nt_s=50,
        processivity_step=0.99005,
        strand_displacement=True,
        exonuclease="none",
        error_rate=1e-4,
        mg_default_mm=10.0,
        primer_length_range=(8, 15),
        default_primer_length=12,
        additive_opt_base_length=12,
        primer_tm_range=(20.0, 55.0),
        gc_range=(0.25, 0.75),
        preset_reaction_temp=37.0,
        thermo_filter=False,
        primer_multiplier=1.0,
        additive_limits={"dmso_percent": 10.0, "betaine_m": 2.0},
        dmso_response={"threshold": 5.0, "mild_coef": 0.02, "steep_coef": 0.10},
        references=("Bambara et al. (1978) JBC 253:413",),
    ),
}


_ALIAS_INDEX: Dict[str, str] = {
    alias.lower(): spec.key for spec in POLYMERASES.values() for alias in spec.aliases
}


def polymerase_names() -> list:
    """Canonical polymerase keys, in registry order (schema/CLI choices order)."""
    return list(POLYMERASES)


def resolve_polymerase_name(name: Optional[str]) -> Optional[str]:
    """Map a name or alias to its canonical key, or None if unknown."""
    if not name:
        return None
    lowered = name.lower()
    if lowered in POLYMERASES:
        return lowered
    return _ALIAS_INDEX.get(lowered)


def get_polymerase(name: Optional[str], default: Optional[str] = None) -> PolymeraseSpec:
    """Look up a polymerase spec by name or alias.

    Args:
        name: Polymerase name or alias.
        default: Fallback key used when ``name`` is unknown. When None, an
            unknown name raises -- matching ``ReactionConditions._validate``.

    Raises:
        ValueError: If ``name`` is unknown and no ``default`` was given.
    """
    resolved = resolve_polymerase_name(name)
    if resolved is not None:
        return POLYMERASES[resolved]
    if default is not None:
        return POLYMERASES[default]
    raise ValueError(f"Unknown polymerase '{name}'. Use: {', '.join(POLYMERASES)}")
