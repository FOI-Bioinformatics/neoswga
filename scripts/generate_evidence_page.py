#!/usr/bin/env python3
"""Render docs/EVIDENCE.md from docs/capability_evidence.json.

The page is generated so that a claim and the evidence for it cannot drift
apart in prose. `tests/test_every_claim_names_its_evidence.py` regenerates it
and compares, so an edit to the page alone fails; an edit to the metadata that
names a file which does not exist fails too.

Usage:

    python scripts/generate_evidence_page.py           # write the page
    python scripts/generate_evidence_page.py --check   # exit 1 if stale
"""

import argparse
import json
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
METADATA = ROOT / "docs" / "capability_evidence.json"
PAGE = ROOT / "docs" / "EVIDENCE.md"

LABELS = {
    "implemented": "implemented",
    "connected": "connected",
    "oracle_tested": "oracle-tested",
    "retrospectively_evaluated": "retrospectively evaluated",
    "prospectively_validated": "prospectively validated",
}

PREAMBLE = """# What each capability's evidence actually is

Generated from `capability_evidence.json`. Do not edit this page by hand; run
`scripts/generate_evidence_page.py`.

Every capability this tool advertises is placed on one of five tiers. The tiers
are ordered, and the line a reader should look for first is the one at the top:
**nothing here is prospectively validated**, because no pool this tool designed
has been synthesised, run and sequenced. [What NeoSWGA does not
establish](LIMITATIONS.md) says that first and in more detail.

The point of a tier is to stop "we built it" and "we measured it" from reading
alike. This repository has shipped a capability that was built, unit-tested,
documented and reachable by no command at all -- six of them at once -- which is
why `connected` is a tier of its own rather than an assumption.
"""

TIER_ORDER = list(LABELS)


def render(data) -> str:
    out = [PREAMBLE, "\n## The tiers\n"]
    for tier in data["tiers"]:
        out.append(f"**{LABELS[tier['name']]}** -- {tier['means']}\n")

    out.append("\n## Where each capability sits\n")
    out.append("| Capability | Tier | Reached by |")
    out.append("|---|---|---|")
    for item in sorted(
        data["capabilities"], key=lambda c: (-TIER_ORDER.index(c["tier"]), c["capability"])
    ):
        out.append(f"| {item['capability']} | {LABELS[item['tier']]} | `{item['reached_by']}` |")

    out.append("\n## The evidence, per capability\n")
    for item in sorted(
        data["capabilities"], key=lambda c: (-TIER_ORDER.index(c["tier"]), c["capability"])
    ):
        out.append(f"### {item['capability']}\n")
        out.append(f"*{item['claim']}*\n")
        out.append(f"**Tier:** {LABELS[item['tier']]}. **Reached by:** `{item['reached_by']}`\n")
        if item["tests"]:
            out.append("Tests: " + ", ".join(f"`{t}`" for t in item["tests"]) + "\n")
        if item["measurements"]:
            out.append(
                "Measurements: "
                + ", ".join(
                    f"[{pathlib.Path(m).name}]({pathlib.Path(m).relative_to('docs')})"
                    for m in item["measurements"]
                )
                + "\n"
            )
        else:
            out.append("Measurements: none.\n")
        out.append(f"{item['note']}\n")

    return "\n".join(out).rstrip() + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="exit 1 if the page is stale")
    args = parser.parse_args()

    rendered = render(json.loads(METADATA.read_text()))

    if args.check:
        current = PAGE.read_text() if PAGE.exists() else ""
        if current != rendered:
            print(f"{PAGE} is stale; run {pathlib.Path(__file__).name}", file=sys.stderr)
            return 1
        return 0

    PAGE.write_text(rendered)
    print(f"wrote {PAGE}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
