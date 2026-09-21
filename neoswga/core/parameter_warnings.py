"""Warnings issued while loading params.json.

Split out of `parameter.py` on 2026-09-21 to keep that module inside its size
budget. Neither function changes a value; both exist so that a configuration
which will be applied differently than the user expects says so once, at load
time, rather than showing up as an unexplained difference in a delivered panel.
"""

import logging

logger = logging.getLogger(__name__)

__all__ = ["_warn_about_schema_version", "_warn_about_unknown_keys"]


def _warn_about_schema_version(data):
    """Warn when params.json declares no schema version, or a different one.

    Extracted from ``get_params`` unchanged, to keep that function inside its
    length budget. The messages and the three branches are as they were.
    """
    from .parameter import CURRENT_SCHEMA_VERSION, SCHEMA_V2_MIGRATION_NOTE

    schema_version = data.get("schema_version", None) if isinstance(data, dict) else None
    if schema_version is None:
        logger.warning(
            "params.json has no 'schema_version' field. "
            "Defaults may differ between NeoSWGA versions. "
            f"Add '\"schema_version\": {CURRENT_SCHEMA_VERSION}' to your "
            "params.json for reproducibility."
        )
    elif schema_version < CURRENT_SCHEMA_VERSION:
        logger.warning(
            f"params.json declares schema_version {schema_version}; this "
            f"NeoSWGA uses version {CURRENT_SCHEMA_VERSION}. Several "
            f"scientific constants were corrected in v2 and results will "
            f"differ from a v1 run:\n"
            f"{SCHEMA_V2_MIGRATION_NOTE}"
        )
    elif schema_version > CURRENT_SCHEMA_VERSION:
        logger.warning(
            f"params.json schema_version {schema_version} is newer than "
            f"this NeoSWGA version supports (max: {CURRENT_SCHEMA_VERSION}). "
            f"Some parameters may not be recognized."
        )


def _warn_about_unknown_keys(data):
    """Warn about params.json keys the schema does not declare.

    `validate params` runs the same check, but a user who never runs it still
    gets one line here rather than a silently applied default: the schema sets
    `additionalProperties: true`, so `max_bg_freqency` was accepted in silence
    and the default for `max_bg_freq` applied, changing the design.

    Never raises. A failure to check is not a reason to fail the run.
    """
    try:
        from neoswga.core.param_validator import unknown_param_keys

        for key, suggestion in unknown_param_keys(data):
            if suggestion:
                logger.warning(
                    "Unknown parameter '%s' in params.json; it will be ignored. "
                    "Did you mean '%s'?",
                    key,
                    suggestion,
                )
            else:
                logger.warning("Unknown parameter '%s' in params.json; it will be ignored.", key)
    except Exception as e:  # pragma: no cover - defensive
        logger.debug(f"Unknown-key check skipped: {e}")
