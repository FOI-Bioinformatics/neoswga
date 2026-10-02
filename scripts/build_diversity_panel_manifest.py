#!/usr/bin/env python3
"""Measure the diversity panel genomes and write the manifest every later script reads.

The panel is two sets over the same files: a target panel of five *Wolbachia*
strains spanning supergroups A, B and D, and a host panel spanning three orders
of magnitude of genome size. The manifest is the one place that records, per
genome, where it is and what it measures, so no later script re-derives a length
or guesses a role.

    python scripts/build_diversity_panel_manifest.py
    python scripts/build_diversity_panel_manifest.py --report-only

The genomes are gitignored and fetched by `scripts/fetch_reference_genomes.py`.
A file that is absent is reported as absent; nothing is estimated for it. A
measurement that could not be made is recorded as null with a note, never as 0.

hg38 is 3.3 Gb, so its SHA-256 costs a full read. The time taken is recorded
alongside the digest, because a reader deciding whether to re-verify needs it.
"""

import argparse
import hashlib
import json
import os
import sys
import time
from collections import Counter
from dataclasses import dataclass

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MANIFEST_PATH = os.path.join("tests", "validation", "genomes", "diversity_panel.json")

# Residues that need no per-character accounting. Everything else in a sequence
# line is counted individually and reported.
EXPECTED = frozenset(b"ACGTN")

# A SHA-256 above this cost is reported with the time it took, so the figure can
# be read as a measurement rather than a free check.
SLOW_CHECKSUM_SECONDS = 180.0


@dataclass(frozen=True)
class PanelEntry:
    key: str
    path: str  # relative to the repository root
    role: str  # "target" or "host"
    group: str  # the panel this entry belongs to
    reference: bool = False  # the one target the design is made against


PANEL: tuple[PanelEntry, ...] = (
    PanelEntry(
        key="wolbachia",
        path="tests/validation/genomes/wolbachia.fna",
        role="target",
        group="wolbachia",
        reference=True,
    ),
    PanelEntry(
        key="wolbachia_wri",
        path="tests/validation/genomes/wolbachia_wri.fna",
        role="target",
        group="wolbachia",
    ),
    PanelEntry(
        key="wolbachia_wpip",
        path="tests/validation/genomes/wolbachia_wpip.fna",
        role="target",
        group="wolbachia",
    ),
    PanelEntry(
        key="wolbachia_walbb",
        path="tests/validation/genomes/wolbachia_walbb.fna",
        role="target",
        group="wolbachia",
    ),
    PanelEntry(
        key="wolbachia_wbm",
        path="tests/validation/genomes/wolbachia_wbm.fna",
        role="target",
        group="wolbachia",
    ),
    PanelEntry(
        key="lactobacillus",
        path="tests/validation/genomes/lactobacillus.fna",
        role="host",
        group="host",
    ),
    PanelEntry(
        key="drosophila",
        path="examples/wolbachia_pool_design/input/drosophila.fna",
        role="host",
        group="host",
    ),
    PanelEntry(
        key="human_full",
        path="tests/validation/genomes/human_full.fna",
        role="host",
        group="host",
    ),
)


def measure_fasta(path: str) -> dict:
    """Count records, sequence length and every non-ACGTN residue in one pass.

    Reads the file as bytes a line at a time, which keeps a 3.3 Gb genome to a
    single pass without loading it. Sequence length is the count of residue
    characters after newlines are stripped, which is what
    `utility.get_seq_length` returns for the same file.
    """
    records = 0
    length = 0
    lowercase = 0
    other: Counter = Counter()

    with open(path, "rb") as handle:
        for line in handle:
            if line.startswith(b">"):
                records += 1
                continue
            seq = line.rstrip(b"\r\n")
            if not seq:
                continue
            length += len(seq)
            upper = seq.upper()
            if upper != seq:
                lowercase += sum(1 for byte in seq if 97 <= byte <= 122)
            accounted = sum(upper.count(residue) for residue in EXPECTED)
            if accounted != len(upper):
                other.update(chr(byte) for byte in upper if byte not in EXPECTED)

    return {
        "records": records,
        "length": length,
        "lowercase_residues": lowercase,
        "non_acgtn": dict(sorted(other.items())),
    }


def sha256_of(path: str) -> tuple[str, float]:
    digest = hashlib.sha256()
    started = time.monotonic()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest(), time.monotonic() - started


def describe(entry: PanelEntry) -> dict:
    absolute = os.path.join(REPO_ROOT, entry.path)
    record: dict = {
        "key": entry.key,
        "path": entry.path,
        "role": entry.role,
        "group": entry.group,
    }
    if entry.reference:
        record["reference"] = True

    if not os.path.exists(absolute):
        # Absent is not empty. Nothing downstream may read a length from this.
        record.update(
            {
                "length": None,
                "records": None,
                "sha256": None,
                "bytes": None,
                "missing": True,
                "sha256_note": f"file not present at {entry.path}; nothing measured",
            }
        )
        return record

    measured = measure_fasta(absolute)
    digest, elapsed = sha256_of(absolute)
    record.update(
        {
            "length": measured["length"],
            "records": measured["records"],
            "sha256": digest,
            "bytes": os.path.getsize(absolute),
            "lowercase_residues": measured["lowercase_residues"],
            "non_acgtn": measured["non_acgtn"],
        }
    )
    if elapsed >= SLOW_CHECKSUM_SECONDS:
        record["sha256_note"] = (
            f"computed over the whole file in {elapsed:.0f} s; re-verification is not free"
        )
    return record


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--report-only",
        action="store_true",
        help="Print the measurements without writing the manifest",
    )
    parser.add_argument("-o", "--output", default=MANIFEST_PATH)
    args = parser.parse_args()

    entries = []
    print(f"{'key':<17}{'bytes':>14}{'records':>9}{'length':>14}  non-ACGTN")
    for entry in PANEL:
        record = describe(entry)
        entries.append(record)
        if record.get("missing"):
            print(f"{entry.key:<17}{'absent':>14}")
            continue
        odd = record["non_acgtn"] or "-"
        print(
            f"{entry.key:<17}{record['bytes']:>14,}{record['records']:>9,}"
            f"{record['length']:>14,}  {odd}"
        )

    manifest = {
        "panel": entries,
        "note": (
            "Lengths and record counts are measured from the files named here. "
            "A figure that could not be measured is null, never 0."
        ),
    }

    if args.report_only:
        return 0

    output = args.output if os.path.isabs(args.output) else os.path.join(REPO_ROOT, args.output)
    os.makedirs(os.path.dirname(output), exist_ok=True)
    with open(output, "w") as handle:
        json.dump(manifest, handle, indent=2)
        handle.write("\n")
    print(f"\nWrote {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
