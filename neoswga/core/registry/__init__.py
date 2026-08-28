"""Single source of truth for polymerase and buffer-component definitions.

Stdlib-only by design: ``param_validator``, the argparse ``choices=`` lists and the
schema renderer all import this, and none of them should have to pull in the
numpy/scipy scientific stack just to learn the polymerase names.

Typical use::

    from neoswga.core.registry import POLYMERASES, get_polymerase

    spec = get_polymerase("phi29")
    reach = spec.typical_amplicon_bp      # coverage-scoring reach
    limit = spec.processivity_bp          # single-molecule maximum

See ``INCONSISTENCIES.md`` for values that are seeded verbatim despite being
mutually contradictory, and ``tests/test_registry_invariants.py`` for the assertions
that track them.
"""

from neoswga.core.registry.polymerases import (
    POLYMERASES,
    PolymeraseSpec,
    get_polymerase,
    polymerase_names,
    resolve_polymerase_name,
)
from neoswga.core.registry.views import (
    additive_optimizer_constraints,
    as_characteristics,
    default_primer_lengths,
    hard_temp_ranges,
    mechanistic_enzyme_params,
    mg_defaults,
    primer_length_ranges,
    processivity_map,
    warn_temp_ranges,
)

__all__ = [
    "POLYMERASES",
    "PolymeraseSpec",
    "get_polymerase",
    "polymerase_names",
    "resolve_polymerase_name",
    "as_characteristics",
    "hard_temp_ranges",
    "warn_temp_ranges",
    "primer_length_ranges",
    "default_primer_lengths",
    "mg_defaults",
    "processivity_map",
    "mechanistic_enzyme_params",
    "additive_optimizer_constraints",
]
