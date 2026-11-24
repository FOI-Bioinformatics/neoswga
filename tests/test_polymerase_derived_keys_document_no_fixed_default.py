"""Keys resolved from the polymerase document no fixed default.

`min_tm`/`max_tm` were documented as 15/45 and `dtt_mm` as 0.0, while all three
are filled from the polymerase when absent (`POLYMERASE_DERIVED_FIELDS`). A
user copying the documented 0.0 into params.json pins DTT as omitted, which
the mechanistic model penalises (2026-11-23 sweep, tier 3).
"""

import json
import pathlib

import pytest

from neoswga.core import parameter

ROOT = pathlib.Path(__file__).resolve().parent.parent
SCHEMA = json.loads((ROOT / "neoswga/core/schema/params.schema.json").read_text())["properties"]


@pytest.mark.parametrize(
    "key", sorted(k for k in parameter.POLYMERASE_DERIVED_FIELDS if k in SCHEMA)
)
def test_the_schema_states_no_fixed_default(key):
    entry = SCHEMA[key]
    assert "default" not in entry
    assert "Polymerase-aware default" in entry.get("description", "")


@pytest.mark.parametrize("key", ["min_tm", "max_tm", "dtt_mm"])
def test_the_reference_shows_no_fixed_default(key):
    rows = (ROOT / "docs/params-reference.md").read_text().splitlines()
    row = next(r for r in rows if r.startswith(f"| `{key}` "))
    assert row.split("|")[4].strip() == "-"
