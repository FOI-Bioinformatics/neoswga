"""The per-run audit record written beside the step 4 results.

Split out of `unified_optimizer` because it is a serialisation concern
rather than an optimisation one, and because that module sits against the
size ratchet in `tests/test_module_size_ratchet.py`.
"""

import hashlib
import json
import logging
import os
from typing import List, Optional

from . import parameter
from .base_optimizer import OptimizationResult

logger = logging.getLogger(__name__)


def write_audit_trail(output_path: str, result: OptimizationResult) -> None:
    """Write run metadata alongside pipeline results for reproducibility.

    Creates a JSON file with version, timestamp, parameters hash, and
    optimizer configuration so that results can be traced back to the
    exact conditions that produced them.
    """
    import hashlib
    from datetime import datetime, timezone

    audit_path = os.path.splitext(output_path)[0] + "_audit.json"
    try:
        # Collect parameters hash
        params_dict = {}
        for attr in [
            "fg_genomes",
            "bg_genomes",
            "min_k",
            "max_k",
            "min_fg_freq",
            "max_bg_freq",
            "max_gini",
            "num_primers",
            "target_set_size",
            "polymerase",
            "reaction_temp",
            "na_conc",
            "mg_conc",
            # Below this line: settings that change the delivered set and were
            # absent from the hash, so two runs differing only in, say, the
            # dimer threshold or the coverage reach produced the same
            # `params_hash` and the audit record could not tell them apart.
            "coverage_reach",
            "max_dimer_bp",
            "max_self_dimer_bp",
            "max_mismatches",
            "min_tm",
            "max_tm",
        ]:
            val = getattr(parameter, attr, None)
            if val is not None:
                params_dict[attr] = str(val)

        # The optimizer is the single largest determinant of the result and has
        # no `parameter` global to read -- it arrives as a CLI argument -- so it
        # is taken from the result itself. `seed` and `application` are the same
        # shape and are recorded in the run manifest rather than here.
        if getattr(result, "optimizer_name", None):
            params_dict["optimizer"] = str(result.optimizer_name)

        params_str = json.dumps(params_dict, sort_keys=True)
        params_hash = hashlib.sha256(params_str.encode()).hexdigest()[:16]

        try:
            from neoswga import __version__

            version = __version__
        except ImportError:
            version = "unknown"

        audit = {
            "neoswga_version": version,
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "params_hash": params_hash,
            "optimizer": result.optimizer_name,
            "num_primers": result.num_primers,
            "score": float(result.score) if result.score else None,
            "status": result.status.value if result.status else None,
            "parameters": params_dict,
        }

        with open(audit_path, "w") as f:
            json.dump(audit, f, indent=2)
        logger.info(f"Audit trail saved to {audit_path}")
    except Exception as e:
        logger.debug(f"Could not write audit trail: {e}")
