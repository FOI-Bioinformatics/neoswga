"""Canonical params.json schema packaged with neoswga.

The JSON Schema file here is the authoritative description of a params.json
configuration. It is loaded by :class:`neoswga.core.param_validator.ParamValidator`
and exposed to users via ``neoswga schema --dump``.
"""

from pathlib import Path

SCHEMA_PATH = Path(__file__).parent / "params.schema.json"


def load_schema() -> dict:
    """Return the params.json schema as a dict.

    Raises FileNotFoundError if the schema has not been shipped.
    """
    import json
    with SCHEMA_PATH.open("r", encoding="utf-8") as fh:
        return json.load(fh)
