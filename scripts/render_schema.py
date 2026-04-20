"""Render the canonical params.json JSON Schema into a human-readable
Markdown reference. The output file (`docs/params-reference.md`) is
regenerated from `neoswga/core/schema/params.schema.json`, so the schema
stays the single source of truth.

Run manually:

    python scripts/render_schema.py

CI or pre-commit can wire the script so the generated doc never drifts.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
SCHEMA_PATH = ROOT / "neoswga" / "core" / "schema" / "params.schema.json"
OUTPUT_PATH = ROOT / "docs" / "params-reference.md"

HEADER = """# params.json Reference

This document is generated from `neoswga/core/schema/params.schema.json`.
Do not edit by hand — run `python scripts/render_schema.py` after
changing the schema.

To dump the schema JSON directly:

```bash
neoswga schema --dump > params.schema.json
```

## Required parameters

"""

FOOTER = """
---

## Additional guidance

- See `docs/production-scenarios.md` for scenario-specific recipes
  (phi29 plasmid, equiphi29 bacterium, extreme GC, multi-genome, etc.).
- See `docs/SWGA_SCIENCE.md` for theoretical background.
- Use `neoswga validate-params -j params.json` to check a configuration
  against all validation layers (schema + ranges + interdependencies).
"""


def _type_label(prop: dict) -> str:
    t = prop.get("type")
    if isinstance(t, list):
        return " or ".join(t)
    if t == "array":
        item_type = prop.get("items", {}).get("type", "any")
        return f"array of {item_type}"
    return t or "mixed"


def _range_label(prop: dict) -> str:
    parts = []
    if "minimum" in prop:
        parts.append(f"min: {prop['minimum']}")
    if "maximum" in prop:
        parts.append(f"max: {prop['maximum']}")
    if "enum" in prop:
        parts.append(f"one of: {', '.join(str(x) for x in prop['enum'])}")
    if "minItems" in prop:
        parts.append(f"minItems: {prop['minItems']}")
    return "; ".join(parts)


def render_param_table(schema: dict, required: set[str], only_required: bool) -> list[str]:
    rows = [
        "| Parameter | Type | Range / allowed | Default | Description |",
        "|---|---|---|---|---|",
    ]
    props = schema.get("properties", {})
    for name, prop in sorted(props.items()):
        is_required = name in required
        if only_required and not is_required:
            continue
        if not only_required and is_required:
            continue
        type_label = _type_label(prop)
        range_label = _range_label(prop) or "-"
        default = prop.get("default", "-")
        if default == "-":
            default_cell = "-"
        else:
            default_cell = f"`{default}`"
        description = prop.get("description", "").replace("\n", " ").strip() or "-"
        rows.append(
            f"| `{name}` | {type_label} | {range_label} | {default_cell} | {description} |"
        )
    return rows


def render(schema: dict) -> str:
    required = set(schema.get("required", []))
    lines = [HEADER.rstrip()]

    lines.append("")
    lines.extend(render_param_table(schema, required, only_required=True))
    lines.append("")
    lines.append("## Optional parameters")
    lines.append("")
    lines.extend(render_param_table(schema, required, only_required=False))
    lines.append(FOOTER)
    return "\n".join(lines) + "\n"


def main() -> int:
    if not SCHEMA_PATH.is_file():
        print(f"ERROR: schema not found at {SCHEMA_PATH}", file=sys.stderr)
        return 1
    with SCHEMA_PATH.open("r", encoding="utf-8") as fh:
        schema = json.load(fh)
    rendered = render(schema)
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_PATH.write_text(rendered, encoding="utf-8")
    print(f"Wrote {OUTPUT_PATH} ({len(rendered)} bytes)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
