#!/usr/bin/env python
"""Does a variant table substitute for a strain genome? (Phase 5, gate B)

Phase 5 of `docs/superpowers/plans/2026-10-02-genomic-diversity-and-multi-host.md`.

The question. `evaluate-set --variants` (Phase 4) reports which binding sites of
a panel survive in a strain, from a variant table stated against one reference.
`evaluate-set --genome` (Phase 2) reports what a panel finds by scanning the
strain's own genome. This script measures how far the two answers differ, on the
five *Wolbachia* genomes of `tests/validation/genomes/diversity_panel.json` with
wMel as the reference, for the panels Phase 3 delivered.

Three routes are computed, and the record must not mix them up:

1. **variant route** -- `variant_table.open_variants` plus
   `variant_sites.evaluate_strain_panel`, exactly as `evaluate-set --variants`
   drives them. Denominator: the REFERENCE. A reference site is intact when no
   variant in the table falls inside it. A site in reference sequence the
   aligner never aligned carries no variant, so this route calls it intact.
   That is the route's central weakness and this script counts it.
2. **homologous route** -- derived here from the aligner's base-level `cs` tag,
   independently of the variant caller. A reference site survives when every one
   of its k reference bases is aligned to an identical base and no insertion
   falls inside the site. Denominator: the REFERENCE. This is the like-for-like
   comparison against route 1, and it is the main table: same sites, same
   denominator, same coverage function.
3. **genome route** -- `reference_panel_evaluation.evaluate_reference_panel`
   over the strain's own FASTA. Denominator: the STRAIN. It counts every site
   found in the strain, including sites gained there and sites in sequence the
   reference does not have. Its coverage is NOT comparable term by term with
   routes 1 and 2; it is reported beside them with its own denominator stated.

What this script does not establish is listed in the record,
`docs/validation/2026-10-02-variant-route-agreement.md`.

External tools. `minimap2` and `paftools.js` are not dependencies of this
package and are not on PATH here. Their paths are arguments. Nothing is
installed, and the script stops with the paths it tried if they are absent.

Usage:
    python scripts/benchmarking/variant_route_agreement.py \\
        --minimap2 /path/to/minimap2 --paftools /path/to/paftools.js \\
        --results tests/validation/genomes/diversity_baseline/results.json \\
        --workdir tests/validation/genomes/variant_route_agreement \\
        --output tests/validation/genomes/variant_route_agreement/results.json

    # Verify the cs-tag reader on a hand-checked case, and stop.
    python scripts/benchmarking/variant_route_agreement.py \\
        --minimap2 ... --paftools ... --self-check

    # Regenerate the record's tables from a finished results file.
    python scripts/benchmarking/variant_route_agreement.py \\
        --tables tests/validation/genomes/variant_route_agreement/results.json
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import re
import shutil
import subprocess
import sys
import time

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(os.path.dirname(__file__))))
if REPO not in sys.path:
    sys.path.insert(0, REPO)

logger = logging.getLogger("variant_route_agreement")

#: The reference every variant table here is stated against, and the four
#: strains called against it. Keys are the manifest keys Phase 1 and Phase 3
#: use, so a figure can be read against the divergence record.
REFERENCE_KEY = "wolbachia"
STRAIN_KEYS = ("wolbachia_wri", "wolbachia_wpip", "wolbachia_walbb", "wolbachia_wbm")

#: Short names for the tables; the keys above are what the JSON carries.
SHORT_NAME = {
    "wolbachia": "wMel",
    "wolbachia_wri": "wRi",
    "wolbachia_wpip": "wPip",
    "wolbachia_walbb": "wAlbB",
    "wolbachia_wbm": "wBm",
}

#: Pairwise Jaccard over canonical 12-mers against wMel, from
#: `docs/validation/2026-10-02-wolbachia-panel-divergence.md`. Quoted, not
#: recomputed here; it is the divergence axis gate B is read against.
JACCARD_VS_REFERENCE = {
    "wolbachia_wri": 0.6652,
    "wolbachia_wpip": 0.2571,
    "wolbachia_walbb": 0.2530,
    "wolbachia_wbm": 0.2114,
}

#: Assembly-to-assembly presets to measure. minimap2's own documentation ties
#: them to divergence: asm5 to about 0.1 percent, asm10 to about 1 percent,
#: asm20 to about 5 percent. This panel spans within-supergroup to
#: cross-supergroup, so no single preset is assumed to suit every pair: all
#: three are run and the aligned fraction of each is reported.
PRESETS = ("asm5", "asm10", "asm20")

#: paftools.js call filters. Its defaults are -l 10000 -L 50000, written for
#: large eukaryotic assemblies. A cross-supergroup *Wolbachia* pair aligns in
#: blocks far shorter than 50 kb, so at the default -L almost no block is
#: eligible and almost no variant is called -- which would read as "these
#: genomes are identical". Both are lowered and the value is recorded; the
#: callable region is then computed under the SAME filters, so the region the
#: variant route can see is the region the caller actually used.
CALL_MIN_COV_LEN = 1000
CALL_MIN_VAR_LEN = 1000
CALL_MIN_MAPQ = 5

#: Panel groups taken from Phase 3. D1 is designed on wMel alone, D2 on all
#: five strains pooled. Both are evaluated against every strain, so they give a
#: single-reference panel and a pooled panel at each size. C1 and C2 are Phase
#: 3 controls on the design side and add no new route comparison; H1 and H3
#: involve hg38 and are deliberately untouched here.
DEFAULT_GROUPS = ("D1", "D2")

#: Per-primer reach and geometry, as Phase 3 ran them. Both must match or the
#: coverage figures are not comparable with the ones in results.json.
EXTENSION_BP = 3000
CIRCULAR = True

#: The coverage difference gate B is read against: an absolute difference in
#: coverage fraction on the reference denominator, route 1 minus route 2.
#:
#: It is a STATED choice, not a derived one, and nothing in the package reads
#: it. The reasoning: Phase 3 measured per-strain coverage falling by 0.25 to
#: 0.40 between the reference and a cross-supergroup strain, so a route
#: disagreement of 0.05 is about an eighth of the effect the two routes are
#: being used to measure. Below it the two routes would rank the strains the
#: same way and a reader would draw the same conclusion from either; above it
#: they would not. A different tolerance moves where the boundary falls, and the
#: record reports the raw differences so another one can be applied.
TOLERANCE = 0.05
TOLERANCE_BASIS = (
    "a stated absolute difference in coverage fraction on the reference "
    "denominator, about an eighth of the 0.25-0.40 per-strain coverage loss "
    "Phase 3 measured across this panel; not derived and read by no code"
)


# ----------------------------------------------------------------------
# External tools
# ----------------------------------------------------------------------


class ToolMissing(RuntimeError):
    """A required external tool is absent. Nothing is installed in response."""


def resolve_tool(path: str | None, name: str, hint: str) -> str:
    """The absolute path of an external tool, or a refusal naming what failed."""
    if not path:
        raise ToolMissing(
            f"{name} is needed and no --{name.replace('.', '-')} path was given. "
            f"It is not a dependency of this package and is not on PATH. {hint}"
        )
    resolved = shutil.which(path) or (path if os.path.isfile(path) else None)
    if not resolved:
        raise ToolMissing(
            f"{name} was given as {path!r}, which is neither an executable on "
            f"PATH nor a file that exists. {hint}"
        )
    return os.path.abspath(resolved)


def tool_version(command: list[str]) -> str:
    """One line of version text, or a reason it could not be read."""
    try:
        done = subprocess.run(command, capture_output=True, text=True, timeout=60)
    except Exception as exc:  # noqa: BLE001 - reported, not repaired
        return f"unavailable: {exc}"
    text = (done.stdout or done.stderr or "").strip().splitlines()
    return text[0] if text else "unavailable: the tool printed nothing"


def run(command: list[str], *, stdout_path: str | None = None, env: dict | None = None) -> None:
    """Run a command, raising with its stderr if it fails."""
    logger.info("running: %s", " ".join(command))
    if stdout_path:
        with open(stdout_path, "w") as handle:
            done = subprocess.run(command, stdout=handle, stderr=subprocess.PIPE, text=True)
    else:
        done = subprocess.run(command, capture_output=True, text=True)
    if done.returncode != 0:
        raise RuntimeError(
            f"{command[0]} exited {done.returncode}: "
            f"{(done.stderr or '').strip()[:2000] or 'no stderr'}"
        )


# ----------------------------------------------------------------------
# FASTA and PAF
# ----------------------------------------------------------------------


def read_fasta(path: str) -> tuple[list[tuple[str, int]], str]:
    """Record names with lengths, and the concatenated upper-case sequence.

    Concatenated with no separator, which is the convention every position in
    this package uses (`docs` and `core/position_index.py`).
    """
    names: list[tuple[str, int]] = []
    chunks: list[str] = []
    current: list[str] = []
    name = None
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    seq = "".join(current)
                    names.append((name, len(seq)))
                    chunks.append(seq)
                name = line[1:].split()[0]
                current = []
            else:
                current.append(line.strip().upper())
    if name is not None:
        seq = "".join(current)
        names.append((name, len(seq)))
        chunks.append(seq)
    return names, "".join(chunks)


_CS_TOKEN = re.compile(r"([=*+\-:~])([A-Za-z0-9\[\]]+)")


def parse_paf(path: str) -> list[dict]:
    """Every PAF line as a dict, with the tags it carries."""
    rows: list[dict] = []
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            tags = {}
            for field in fields[12:]:
                parts = field.split(":", 2)
                if len(parts) == 3:
                    tags[parts[0]] = parts[2]
            rows.append(
                {
                    "qname": fields[0],
                    "qlen": int(fields[1]),
                    "qstart": int(fields[2]),
                    "qend": int(fields[3]),
                    "strand": fields[4],
                    "tname": fields[5],
                    "tlen": int(fields[6]),
                    "tstart": int(fields[7]),
                    "tend": int(fields[8]),
                    "nmatch": int(fields[9]),
                    "blocklen": int(fields[10]),
                    "mapq": int(fields[11]),
                    "tags": tags,
                }
            )
    return rows


def apply_cs(
    cs: str,
    tstart: int,
    identical: np.ndarray,
    insertion_after: np.ndarray,
) -> int:
    """Mark reference bases this alignment aligns to an identical base.

    Walks the `cs` tag along the TARGET (reference), which is what the tag
    does for either strand: minimap2 writes it from `tstart` to `tend` in
    reference orientation, with query bases reverse-complemented when the
    alignment is on the minus strand. So the arrays below are in reference
    forward coordinates whichever strand the query aligned on.

    Operations, in `--cs=long` form:
      `=SEQ`  identical run; marks every reference base it spans.
      `:N`    identical run as a count (short `cs`); same effect.
      `*xy`   one substituted reference base; consumes it, marks nothing.
      `-SEQ`  reference bases absent from the query; consumes, marks nothing.
      `+SEQ`  query bases absent from the reference; consumes no reference
              base and records a break between the previous reference base and
              the next, because a site spanning it is not contiguous in the
              strain.
      `~ab##cd` an intron. Not expected between two bacterial assemblies;
              treated as consuming the stated reference length and marking
              nothing, and counted by the caller so a run that meets one says so.

    Args:
        cs: the tag's value, without the `cs:Z:` prefix.
        tstart: this alignment's 0-based reference start.
        identical: bool array over the reference, updated in place with OR.
        insertion_after: bool array over the reference, updated in place; index
            `i` true means a query insertion sits between reference bases `i`
            and `i + 1`.

    Returns:
        The reference position the walk ended on. The caller asserts it equals
        the alignment's `tend`; a disagreement means the tag was misread and
        every figure downstream would be wrong, so it is an error and not a
        warning.
    """
    position = int(tstart)
    length = int(identical.size)
    for operation, value in _CS_TOKEN.findall(cs):
        if operation == "=":
            span = len(value)
            identical[position : position + span] = True
            position += span
        elif operation == ":":
            span = int(value)
            identical[position : position + span] = True
            position += span
        elif operation == "*":
            position += 1
        elif operation == "-":
            position += len(value)
        elif operation == "+":
            if 0 <= position - 1 < length:
                insertion_after[position - 1] = True
        elif operation == "~":
            match = re.match(r"^[a-z]{2}(\d+)[a-z]{2}$", value)
            if not match:
                raise ValueError(f"unparsable cs intron operation {value!r}")
            position += int(match.group(1))
        else:  # pragma: no cover - the regex admits nothing else
            raise ValueError(f"unknown cs operation {operation!r}")
    return position


def alignment_arrays(rows: list[dict], reference_length: int, *, primary_only: bool = True):
    """Reference-wide alignment state from a whole PAF.

    Returns `(aligned, identical, insertion_after, used, introns)`:
      aligned          bool, reference base lies in some alignment block
      identical        bool, reference base is aligned to an identical base
      insertion_after  bool, a query insertion sits after this reference base
      used             how many alignment rows were walked
      introns          how many `~` operations were met (expected: 0)
    """
    aligned = np.zeros(reference_length, dtype=bool)
    identical = np.zeros(reference_length, dtype=bool)
    insertion_after = np.zeros(reference_length, dtype=bool)
    used = 0
    introns = 0
    for row in rows:
        if primary_only and row["tags"].get("tp") not in (None, "P"):
            continue
        cs = row["tags"].get("cs")
        if cs is None:
            raise ValueError("a PAF line carries no cs tag; minimap2 must be run with --cs=long")
        introns += cs.count("~")
        aligned[row["tstart"] : row["tend"]] = True
        end = apply_cs(cs, row["tstart"], identical, insertion_after)
        if end != row["tend"]:
            raise ValueError(
                f"the cs tag of {row['qname']}:{row['qstart']}-{row['qend']} walked to "
                f"reference {end} but the PAF says the alignment ends at {row['tend']}; "
                f"the tag was misread and no figure derived from it can be trusted"
            )
        used += 1
    return aligned, identical, insertion_after, used, introns


def callable_mask(rows: list[dict], reference_length: int) -> np.ndarray:
    """Reference bases inside a block paftools.js call was allowed to use.

    The same two filters the caller is given: a block at least
    `CALL_MIN_VAR_LEN` long and mapping quality at least `CALL_MIN_MAPQ`.
    Outside this mask the variant table carries no variant because no variant
    was looked for, which is not the same as a region with none.
    """
    mask = np.zeros(reference_length, dtype=bool)
    for row in rows:
        if row["tags"].get("tp") not in (None, "P"):
            continue
        if row["blocklen"] < CALL_MIN_VAR_LEN or row["mapq"] < CALL_MIN_MAPQ:
            continue
        mask[row["tstart"] : row["tend"]] = True
    return mask


def interval_union_length(rows: list[dict], *, primary_only: bool = True) -> int:
    """Reference bases covered by at least one block, counted once."""
    spans = sorted(
        (row["tstart"], row["tend"])
        for row in rows
        if not primary_only or row["tags"].get("tp") in (None, "P")
    )
    total = 0
    current_start = current_end = -1
    for start, end in spans:
        if start > current_end:
            total += max(0, current_end - current_start)
            current_start, current_end = start, end
        else:
            current_end = max(current_end, end)
    total += max(0, current_end - current_start)
    return total


# ----------------------------------------------------------------------
# Site-level survival
# ----------------------------------------------------------------------


def site_survives(
    positions: np.ndarray,
    k: int,
    identical: np.ndarray,
    insertion_after: np.ndarray,
    *,
    circular: bool,
) -> np.ndarray:
    """Which reference sites survive at their homologous locus in the strain.

    A site at `p` occupies `[p, p + k)`. It survives when all k reference bases
    are aligned to an identical base and no insertion falls strictly inside the
    site, so the strain holds the same k bases contiguously at that locus.

    Circular geometry wraps the index, matching the convention
    `variant_sites.SiteGeometry` and `string_search` use for a circular
    reference. Under `circular=False` a site that would run past the end cannot
    be assessed and is returned as not surviving; no such site exists in this
    measurement, since the panel's sites come from a circular scan of a
    single-record genome.
    """
    length = int(identical.size)
    positions = np.asarray(positions, dtype=np.int64)
    if positions.size == 0:
        return np.zeros(0, dtype=bool)
    offsets = np.arange(k, dtype=np.int64)
    index = positions[:, None] + offsets[None, :]
    if circular:
        index = index % length
    else:
        if np.any(index >= length):
            inside = ~np.any(index >= length, axis=1)
            index = np.clip(index, 0, length - 1)
            ok = np.all(identical[index], axis=1) & ~np.any(insertion_after[index[:, :-1]], axis=1)
            return ok & inside
    all_identical = np.all(identical[index], axis=1)
    no_insertion = ~np.any(insertion_after[index[:, :-1]], axis=1)
    return all_identical & no_insertion


def uniform_window_fraction(
    identical: np.ndarray,
    insertion_after: np.ndarray | None,
    k: int,
    *,
    circular: bool,
) -> float:
    """The same verdict over EVERY reference position, as a control.

    The panel's sites are not placed uniformly: a SWGA candidate is chosen for
    foreground frequency, so its sites cluster where the composition it likes
    does. Comparing the panel's survival rate against this figure is what
    separates "the aligner aligns half the reference" from "the panel's sites
    sit in the half it does not align".

    With `insertion_after` None this is simply the fraction of k-windows whose
    every base satisfies `identical`; with it, the window must also carry no
    insertion inside, which is the homologous route's own rule.
    """
    length = int(identical.size)
    extended = np.concatenate([identical, identical[: k - 1]]) if circular else identical
    bad = np.concatenate([[0], np.cumsum(~extended)])
    window_all = (bad[k : k + length] - bad[0:length]) == 0
    if insertion_after is not None:
        extended_insert = (
            np.concatenate([insertion_after, insertion_after[: k - 1]])
            if circular
            else insertion_after
        )
        inserts = np.concatenate([[0], np.cumsum(extended_insert)])
        # An insertion strictly inside the window: after base p .. p + k - 2.
        window_all &= (inserts[k - 1 : k - 1 + length] - inserts[0:length]) == 0
    return float(window_all[:length].mean())


def sites_fully_inside(
    positions: np.ndarray, k: int, mask: np.ndarray, *, circular: bool
) -> np.ndarray:
    """Which sites lie entirely inside a reference mask."""
    length = int(mask.size)
    positions = np.asarray(positions, dtype=np.int64)
    if positions.size == 0:
        return np.zeros(0, dtype=bool)
    index = positions[:, None] + np.arange(k, dtype=np.int64)[None, :]
    index = index % length if circular else np.clip(index, 0, length - 1)
    return np.all(mask[index], axis=1)


# ----------------------------------------------------------------------
# The variant table
# ----------------------------------------------------------------------


def call_variants(
    minimap2: str,
    paftools: str,
    reference: str,
    strain_fasta: str,
    preset: str,
    workdir: str,
    sample: str,
) -> tuple[str, str, list[dict]]:
    """Align one strain to the reference and call variants against it.

    Returns `(paf_path, vcf_path, paf_rows)`. The PAF is kept because every
    alignment-derived figure in this script is read from it, so the record's
    aligned fractions and its variant counts come from one alignment and not
    from two runs of the aligner.
    """
    paf = os.path.join(workdir, f"{sample}.{preset}.paf")
    run(
        [
            minimap2,
            "-c",
            "-x",
            preset,
            "--cs=long",
            "-t",
            "2",
            reference,
            strain_fasta,
        ],
        stdout_path=paf,
    )
    rows = parse_paf(paf)

    # paftools.js call reads a PAF sorted by target name then target start,
    # which its own usage line states.
    sorted_paf = os.path.join(workdir, f"{sample}.{preset}.sorted.paf")
    with open(paf) as handle:
        lines = handle.read().splitlines()
    with open(sorted_paf, "w") as handle:
        for line in sorted(lines, key=lambda text: (text.split("\t")[5], int(text.split("\t")[7]))):
            handle.write(line + "\n")

    vcf = os.path.join(workdir, f"{sample}.{preset}.vcf")
    run(
        [
            paftools,
            "call",
            "-l",
            str(CALL_MIN_COV_LEN),
            "-L",
            str(CALL_MIN_VAR_LEN),
            "-q",
            str(CALL_MIN_MAPQ),
            "-f",
            reference,
            "-s",
            sample,
            sorted_paf,
        ],
        stdout_path=vcf,
    )
    return paf, vcf, rows


def count_vcf_variants(path: str) -> dict:
    """SNPs, insertions and deletions in a VCF, by REF and ALT length."""
    snps = insertions = deletions = other = records = 0
    with open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            records += 1
            ref, alt = fields[3], fields[4].split(",")[0]
            if len(ref) == 1 and len(alt) == 1:
                snps += 1
            elif len(ref) > len(alt):
                deletions += 1
            elif len(alt) > len(ref):
                insertions += 1
            else:
                other += 1
    return {
        "records": records,
        "snps": snps,
        "insertions": insertions,
        "deletions": deletions,
        "other": other,
    }


# ----------------------------------------------------------------------
# Panels
# ----------------------------------------------------------------------


def read_panels(results_path: str, groups: tuple[str, ...], sizes: tuple[int, ...]) -> list[dict]:
    """The delivered sets Phase 3 recorded, for the groups and sizes asked for.

    Each entry carries the design id, the group, the requested size, the
    delivered primers, and the genome-route per-strain figures Phase 3 already
    measured, so this script can compare its own recomputation against them.
    """
    with open(results_path) as handle:
        results = json.load(handle)
    panels: list[dict] = []
    for design in results["designs"]:
        if design["group"] not in groups:
            continue
        for entry in design["sizes"]:
            if int(entry["requested_size"]) not in sizes:
                continue
            evaluation = entry.get("evaluation") or {}
            panels.append(
                {
                    "design": design["id"],
                    "group": design["group"],
                    "requested_size": int(entry["requested_size"]),
                    "delivered_size": entry.get("delivered_size"),
                    "primers": list(entry["delivered_primers"]),
                    "phase3_per_target": {
                        key: {
                            "sites": (record.get("sites") or {}).get("value"),
                            "coverage": (record.get("coverage") or {}).get("value"),
                            "length_bp": record.get("length_bp"),
                            "status": record.get("status"),
                            "per_primer_sites": dict(record.get("per_primer_sites") or {}),
                        }
                        for key, record in (evaluation.get("per_target") or {}).items()
                    },
                }
            )
    return panels


# ----------------------------------------------------------------------
# The measurement
# ----------------------------------------------------------------------


def measure(args: argparse.Namespace) -> dict:
    from neoswga.core import coverage as coverage_module
    from neoswga.core.exceptions import DesignError
    from neoswga.core.position_cache import PositionCache
    from neoswga.core.reference_panel_evaluation import ReferenceSpec, evaluate_reference_panel
    from neoswga.core.variant_sites import (
        IntactPositions,
        evaluate_strain_panel,
        split_sites_by_strand,
    )
    from neoswga.core.variant_table import open_variants

    minimap2 = resolve_tool(
        args.minimap2,
        "minimap2",
        "Pass --minimap2 with a path; conda environments on this machine carry one.",
    )
    paftools = resolve_tool(
        args.paftools,
        "paftools.js",
        "Pass --paftools with a path. It needs the k8 interpreter beside it; "
        "pass --k8 if it is not found automatically.",
    )
    k8 = shutil.which("k8") or (args.k8 if args.k8 and os.path.isfile(args.k8) else None)
    if args.k8:
        k8 = resolve_tool(args.k8, "k8", "paftools.js is a k8 script and cannot run without it.")
    env_path = os.environ.get("PATH", "")
    if k8:
        os.environ["PATH"] = os.path.dirname(os.path.abspath(k8)) + os.pathsep + env_path
    else:
        os.environ["PATH"] = os.path.dirname(paftools) + os.pathsep + env_path

    tools = {
        "minimap2": {"path": minimap2, "version": tool_version([minimap2, "--version"])},
        "paftools.js": {"path": paftools, "version": tool_version([paftools, "version"])},
        "k8": {"path": k8, "version": "bundled with paftools.js" if k8 else "not resolved"},
        "presets_measured": list(args.presets),
        "call_preset": args.call_preset,
        "paftools_call_options": {
            "-l": CALL_MIN_COV_LEN,
            "-L": CALL_MIN_VAR_LEN,
            "-q": CALL_MIN_MAPQ,
            "basis": (
                "the defaults, -l 10000 -L 50000, leave almost no eligible block "
                "in a cross-supergroup Wolbachia pair, which would read as "
                "'no variants'; the callable region is computed under the same "
                "two filters"
            ),
        },
    }

    os.makedirs(args.workdir, exist_ok=True)
    genomes = args.genomes
    reference_fasta = os.path.join(genomes, "wolbachia.fna")
    strain_fasta = {
        "wolbachia_wri": os.path.join(genomes, "wolbachia_wri.fna"),
        "wolbachia_wpip": os.path.join(genomes, "wolbachia_wpip.fna"),
        "wolbachia_walbb": os.path.join(genomes, "wolbachia_walbb.fna"),
        "wolbachia_wbm": os.path.join(genomes, "wolbachia_wbm.fna"),
    }
    for path in [reference_fasta, *strain_fasta.values()]:
        if not os.path.isfile(path):
            raise FileNotFoundError(f"{path} is missing; this measurement needs all five genomes")

    records, reference_sequence = read_fasta(reference_fasta)
    reference_length = len(reference_sequence)
    logger.info(
        "reference %s: %d record(s), %d bp", reference_fasta, len(records), reference_length
    )

    panels = read_panels(args.results, tuple(args.groups), tuple(args.sizes))
    if not panels:
        raise SystemExit("no panel matched --groups / --sizes")
    every_primer = sorted({primer for panel in panels for primer in panel["primers"]})
    lengths = sorted({len(primer) for primer in every_primer})
    logger.info(
        "%d panel(s), %d distinct primer(s), length(s) %s", len(panels), len(every_primer), lengths
    )

    # One scan of the reference for every primer any panel uses. The prefix is
    # inside the work directory and carries no index, so PositionCache scans
    # the FASTA; Phase 3's indexes are not opened.
    reference_prefix = os.path.join(args.workdir, "reference_scan", REFERENCE_KEY)
    os.makedirs(os.path.dirname(reference_prefix), exist_ok=True)
    cache = PositionCache(
        [reference_prefix],
        every_primer,
        genome_paths=[reference_fasta],
        circular=CIRCULAR,
        on_missing="scan",
    )

    # ---- per strain: alignment, variant table, derived arrays
    strains: dict[str, dict] = {}
    for key in STRAIN_KEYS:
        name = SHORT_NAME[key]
        started = time.time()
        per_preset = {}
        for preset in args.presets:
            paf = os.path.join(args.workdir, f"{name}.{preset}.paf")
            run(
                [
                    minimap2,
                    "-c",
                    "-x",
                    preset,
                    "--cs=long",
                    "-t",
                    "2",
                    reference_fasta,
                    strain_fasta[key],
                ],
                stdout_path=paf,
            )
            rows = parse_paf(paf)
            primary = [row for row in rows if row["tags"].get("tp") in (None, "P")]
            union = interval_union_length(rows)
            per_preset[preset] = {
                "paf": os.path.relpath(paf, REPO),
                "alignment_rows": len(rows),
                "primary_rows": len(primary),
                "secondary_rows": len(rows) - len(primary),
                "reference_bases_aligned": union,
                "reference_fraction_aligned": union / reference_length,
                "longest_block_bp": max((row["blocklen"] for row in primary), default=0),
                "median_block_bp": (
                    float(np.median([row["blocklen"] for row in primary])) if primary else None
                ),
                "callable_reference_bases": int(callable_mask(rows, reference_length).sum()),
            }
            per_preset[preset]["callable_reference_fraction"] = (
                per_preset[preset]["callable_reference_bases"] / reference_length
            )

        preset = args.call_preset
        paf, vcf, rows = call_variants(
            minimap2, paftools, reference_fasta, strain_fasta[key], preset, args.workdir, name
        )
        aligned, identical, insertion_after, used, introns = alignment_arrays(
            rows, reference_length
        )
        visible = callable_mask(rows, reference_length)
        counts = count_vcf_variants(vcf)

        refusal = None
        table = None
        try:
            table = open_variants(vcf, reference_fasta, prefix=reference_prefix)
        except DesignError as exc:
            refusal = f"{type(exc).__name__}: {exc}"
            logger.error("open_variants refused %s: %s", vcf, refusal)

        strains[key] = {
            "strain": key,
            "short_name": name,
            "fasta": os.path.relpath(strain_fasta[key], REPO),
            "jaccard_12mer_vs_reference": JACCARD_VS_REFERENCE[key],
            "presets": per_preset,
            "call_preset": preset,
            "vcf": os.path.relpath(vcf, REPO),
            "variant_counts": counts,
            "alignment_rows_walked": used,
            "cs_intron_operations": introns,
            "reference_bases_aligned": int(aligned.sum()),
            "reference_fraction_aligned": float(aligned.sum()) / reference_length,
            "reference_bases_identical": int(identical.sum()),
            "reference_fraction_identical": float(identical.sum()) / reference_length,
            "reference_bases_callable": int(visible.sum()),
            "reference_fraction_callable": float(visible.sum()) / reference_length,
            "uniform_window_controls": {
                str(k): {
                    "fraction_of_all_windows_inside_the_callable_region": uniform_window_fraction(
                        visible, None, k, circular=CIRCULAR
                    ),
                    "fraction_of_all_windows_surviving_by_alignment": uniform_window_fraction(
                        identical, insertion_after, k, circular=CIRCULAR
                    ),
                    "basis": (
                        "every reference position, not the panel's sites; the "
                        "control the panel's own rates are read against"
                    ),
                }
                for k in lengths
            },
            "open_variants_refused": refusal,
            "variant_table": table.as_dict() if table is not None else None,
            "seconds": round(time.time() - started, 2),
            "_arrays": {
                "identical": identical,
                "insertion_after": insertion_after,
                "visible": visible,
                "table": table,
            },
        }

    # ---- genome-route recomputation on a sample, to confirm like with like
    genome_route_checks = []
    for panel in panels[: args.genome_route_sample]:
        specs = [
            ReferenceSpec(
                prefix=os.path.join(args.workdir, "genome_route", key),
                genome=path,
                length=len(read_fasta(path)[1]),
                role="target",
                circular=CIRCULAR,
                scan=True,
            )
            for key, path in [(REFERENCE_KEY, reference_fasta), *strain_fasta.items()]
        ]
        os.makedirs(os.path.join(args.workdir, "genome_route"), exist_ok=True)
        assessment = evaluate_reference_panel(panel["primers"], specs, None, extension=EXTENSION_BP)
        recomputed = {}
        for record in assessment.references:
            key = os.path.basename(record.prefix)
            recorded = panel["phase3_per_target"].get(key) or {}
            recomputed[key] = {
                "sites_recomputed": record.sites.value,
                "sites_phase3": recorded.get("sites"),
                "coverage_recomputed": record.coverage.value,
                "coverage_phase3": recorded.get("coverage"),
                "coverage_difference": (
                    None
                    if record.coverage.value is None or recorded.get("coverage") is None
                    else record.coverage.value - recorded["coverage"]
                ),
            }
        genome_route_checks.append(
            {
                "design": panel["design"],
                "requested_size": panel["requested_size"],
                "per_target": recomputed,
            }
        )

    # ---- per panel, per strain: the three routes
    panel_records = []
    for panel in panels:
        primers = panel["primers"]
        reference_sites = {
            primer: split_sites_by_strand(cache, reference_prefix, primer) for primer in primers
        }
        union_sites = {
            primer: np.unique(np.concatenate([per_strand["forward"], per_strand["reverse"]]))
            for primer, per_strand in reference_sites.items()
        }
        reference_total = int(sum(array.size for array in union_sites.values()))

        reference_coverage, _ = coverage_module.compute_per_prefix_coverage(
            cache,
            primers,
            [reference_prefix],
            [reference_length],
            extension=EXTENSION_BP,
            circular=CIRCULAR,
        )

        per_strain = {}
        for key in STRAIN_KEYS:
            state = strains[key]
            arrays = state["_arrays"]
            table = arrays["table"]

            # route 1, exactly as evaluate-set --variants drives it
            if table is None:
                variant_block = None
            else:
                variant_block = evaluate_strain_panel(
                    cache,
                    reference_prefix,
                    primers,
                    reference_length,
                    table,
                    extension=EXTENSION_BP,
                    circular=CIRCULAR,
                )

            # route 2, from the cs tag, independently of the caller
            surviving = {}
            homologous_total = 0
            blind_total = 0
            blind_surviving = 0
            variant_intact_in_blind = 0
            agree = disagree_variant_optimistic = disagree_variant_pessimistic = 0
            for primer in primers:
                k = len(primer)
                per_strand = reference_sites[primer]
                kept = {}
                for strand in ("forward", "reverse"):
                    positions = per_strand[strand]
                    mask = site_survives(
                        positions,
                        k,
                        arrays["identical"],
                        arrays["insertion_after"],
                        circular=CIRCULAR,
                    )
                    kept[strand] = np.asarray(positions, dtype=np.int64)[mask]
                surviving[primer] = kept

                union = union_sites[primer]
                survive_union = site_survives(
                    union, k, arrays["identical"], arrays["insertion_after"], circular=CIRCULAR
                )
                homologous_total += int(survive_union.sum())

                inside = sites_fully_inside(union, k, arrays["visible"], circular=CIRCULAR)
                blind_total += int((~inside).sum())
                blind_surviving += int((survive_union & ~inside).sum())

                if table is not None:
                    intact = variant_route_intact_mask(union, k, table, CIRCULAR, reference_length)
                    variant_intact_in_blind += int((intact & ~inside).sum())
                    agree += int((intact == survive_union).sum())
                    disagree_variant_optimistic += int((intact & ~survive_union).sum())
                    disagree_variant_pessimistic += int((~intact & survive_union).sum())

            # Consistency: a reference site that survives at its homologous
            # locus means the strain holds that exact k-mer, so the genome route
            # must find at least as many sites of that primer in the strain as
            # the alignment says survive. A failure here would mean the cs walk
            # placed a site in the wrong coordinate or the wrong orientation,
            # and every figure in this row would be wrong. It is reported as a
            # count rather than assumed.
            strain_per_primer = (panel["phase3_per_target"].get(key) or {}).get(
                "per_primer_sites"
            ) or {}
            inconsistent = []
            for primer in primers:
                survive_here = int(
                    np.unique(
                        np.concatenate([surviving[primer]["forward"], surviving[primer]["reverse"]])
                    ).size
                )
                found = strain_per_primer.get(primer)
                if found is not None and survive_here > int(found):
                    inconsistent.append(
                        {"primer": primer, "surviving": survive_here, "sites_in_strain": int(found)}
                    )

            homologous_view = IntactPositions(cache, reference_prefix, surviving)
            homologous_coverage, _ = coverage_module.compute_per_prefix_coverage(
                homologous_view,
                primers,
                [reference_prefix],
                [reference_length],
                extension=EXTENSION_BP,
                circular=CIRCULAR,
            )

            variant_intact = None
            variant_coverage = None
            if variant_block is not None:
                strain_name = next(iter(variant_block["per_strain"]))
                record = variant_block["per_strain"][strain_name]
                variant_intact = record["intact_sites"]["value"]
                variant_coverage = record["coverage_on_intact_sites"]["value"]

            genome = panel["phase3_per_target"].get(key) or {}
            per_strain[key] = {
                "short_name": SHORT_NAME[key],
                "jaccard_12mer_vs_reference": JACCARD_VS_REFERENCE[key],
                "reference_sites": reference_total,
                "variant_route": {
                    "intact_sites": variant_intact,
                    "lost_sites": (
                        None if variant_intact is None else reference_total - int(variant_intact)
                    ),
                    "coverage_of_reference": variant_coverage,
                    "denominator": "the reference genome, wMel, 1267782 bp",
                },
                "homologous_route": {
                    "surviving_sites": homologous_total,
                    "primers_surviving_more_than_the_strain_holds": inconsistent,
                    "lost_sites": reference_total - homologous_total,
                    "coverage_of_reference": homologous_coverage,
                    "denominator": "the reference genome, wMel, 1267782 bp",
                },
                "genome_route": {
                    "sites_in_strain": genome.get("sites"),
                    "coverage_of_strain": genome.get("coverage"),
                    "strain_length_bp": genome.get("length_bp"),
                    "denominator": "the strain's own genome",
                    "sites_not_explained_by_a_surviving_reference_site": (
                        None
                        if genome.get("sites") is None
                        else int(genome["sites"]) - homologous_total
                    ),
                },
                "like_for_like": {
                    "coverage_difference_variant_minus_homologous": (
                        None if variant_coverage is None else variant_coverage - homologous_coverage
                    ),
                    "site_difference_variant_minus_homologous": (
                        None if variant_intact is None else int(variant_intact) - homologous_total
                    ),
                    "sites_agreeing": agree if table is not None else None,
                    "sites_variant_route_calls_intact_and_alignment_does_not": (
                        disagree_variant_optimistic if table is not None else None
                    ),
                    "sites_alignment_keeps_and_variant_route_calls_affected": (
                        disagree_variant_pessimistic if table is not None else None
                    ),
                },
                "unaligned": {
                    "reference_sites_outside_the_callable_region": blind_total,
                    "of_those_surviving_by_alignment": blind_surviving,
                    "of_those_called_intact_by_the_variant_route": (
                        variant_intact_in_blind if table is not None else None
                    ),
                    "basis": (
                        "a site is outside the callable region when any of its k "
                        "bases lies outside every alignment block that passed "
                        "paftools.js call's own -L and -q filters; the table "
                        "carries no variant there because none was looked for"
                    ),
                },
            }

        panel_records.append(
            {
                "design": panel["design"],
                "group": panel["group"],
                "requested_size": panel["requested_size"],
                "delivered_size": panel["delivered_size"],
                "primers": primers,
                "reference_sites": reference_total,
                "reference_coverage": reference_coverage,
                "per_strain": per_strain,
            }
        )
        logger.info(
            "%s n=%s: %d reference site(s), coverage %.4f",
            panel["design"],
            panel["requested_size"],
            reference_total,
            reference_coverage,
        )

    for state in strains.values():
        state.pop("_arrays", None)

    return {
        "script": os.path.relpath(os.path.abspath(__file__), REPO),
        "written": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "python": sys.version.split()[0],
        "tools": tools,
        "reference": {
            "key": REFERENCE_KEY,
            "short_name": SHORT_NAME[REFERENCE_KEY],
            "fasta": os.path.relpath(reference_fasta, REPO),
            "records": [{"name": name, "length": length} for name, length in records],
            "length_bp": reference_length,
            "circular": CIRCULAR,
        },
        "extension_reach_bp": EXTENSION_BP,
        "tolerance": TOLERANCE,
        "tolerance_basis": TOLERANCE_BASIS,
        "groups": list(args.groups),
        "sizes": list(args.sizes),
        "phase3_results": os.path.relpath(os.path.abspath(args.results), REPO),
        "strains": strains,
        "genome_route_checks": genome_route_checks,
        "panels": panel_records,
        "routes": {
            "variant_route": "variant_table.open_variants + variant_sites.evaluate_strain_panel, reference denominator",
            "homologous_route": "derived here from the cs tag, reference denominator, independent of the variant caller",
            "genome_route": "reference_panel_evaluation.evaluate_reference_panel over the strain FASTA, strain denominator",
        },
        "limits": [
            "Five genomes of one genus, one reference, one aligner, one variant caller.",
            "The homologous route is derived from the SAME alignment the variant "
            "table is called from, so it is independent of the caller but not of "
            "the aligner. A different aligner or preset would move both.",
            "Coverage by the genome route has the strain's own length as its "
            "denominator and is not comparable term by term with the two "
            "reference-denominator routes.",
            "No clonal pair is on disk; the closest pair here is within one "
            "supergroup at Jaccard 0.6652.",
        ],
    }


def variant_route_intact_mask(
    positions: np.ndarray, k: int, table, circular: bool, length: int
) -> np.ndarray:
    """The variant route's intact verdict per site, for the cross-tabulation.

    Calls `variant_sites.intact_mask` with the same geometry
    `evaluate_strain_panel` uses, on the one strain the table describes, so the
    per-site cross-tabulation agrees with the aggregate the product function
    reports.
    """
    from neoswga.core.variant_sites import SiteGeometry, intact_mask

    strain = table.measured_strains[0]
    geometry = SiteGeometry(length=int(length), circular=bool(circular))
    return np.asarray(
        intact_mask(np.asarray(positions, dtype=np.int64), k, strain.starts, strain.ends, geometry),
        dtype=bool,
    )


# ----------------------------------------------------------------------
# Self-check on the cs reader
# ----------------------------------------------------------------------


def self_check(args: argparse.Namespace) -> int:
    """Verify the cs reader on a hand-written tag and on a real alignment.

    The conversion from the aligner's output to the arrays every figure rests
    on lives in this script, so it is checked here rather than assumed.
    """
    failures = []

    # 1. A hand-written cs tag with a hand-computed expectation.
    #
    #    reference: positions 0..19
    #    cs:  =ACGTA *ag =CCCC +tt =GGGG -acg =TT
    #    walk: 0-4 identical; 5 substituted; 6-9 identical; insertion between
    #          9 and 10; 10-13 identical; 14-16 deleted from the query;
    #          17-18 identical. End at 19.
    identical = np.zeros(20, dtype=bool)
    insertion_after = np.zeros(20, dtype=bool)
    end = apply_cs("=ACGTA*ag=CCCC+tt=GGGG-acg=TT", 0, identical, insertion_after)
    expected_identical = np.zeros(20, dtype=bool)
    expected_identical[[0, 1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 17, 18]] = True
    expected_insertion = np.zeros(20, dtype=bool)
    expected_insertion[9] = True
    if end != 19:
        failures.append(f"hand-written cs walked to {end}, expected 19")
    if not np.array_equal(identical, expected_identical):
        failures.append(
            "hand-written cs identical mask disagrees: got "
            f"{np.flatnonzero(identical).tolist()}, expected "
            f"{np.flatnonzero(expected_identical).tolist()}"
        )
    if not np.array_equal(insertion_after, expected_insertion):
        failures.append(
            "hand-written cs insertion mask disagrees: got "
            f"{np.flatnonzero(insertion_after).tolist()}"
        )

    # 2. site_survives over that hand-computed state, k = 4.
    #    position 0: bases 0-3 all identical, no insertion inside -> survives
    #    position 3: base 5 substituted -> does not survive
    #    position 6: bases 6-9 identical but an insertion sits after 9, which is
    #                NOT inside [6, 10): the break is at the site's right edge,
    #                so the site survives
    #    position 7: bases 7-10 identical, insertion after 9 is inside -> no
    #    position 13: base 14 deleted -> no
    verdicts = site_survives(
        np.array([0, 3, 6, 7, 13]), 4, identical, insertion_after, circular=False
    )
    expected = np.array([True, False, True, False, False])
    if not np.array_equal(verdicts, expected):
        failures.append(
            f"site_survives disagrees: got {verdicts.tolist()}, expected {expected.tolist()}"
        )

    # 3. End to end on a real alignment of a 5 kb slice of the reference with
    #    one substitution and one 3 bp deletion at known offsets.
    minimap2 = resolve_tool(
        args.minimap2, "minimap2", "Pass --minimap2 to run the end-to-end check."
    )
    workdir = os.path.join(args.workdir, "self_check")
    os.makedirs(workdir, exist_ok=True)
    reference_fasta = os.path.join(args.genomes, "wolbachia.fna")
    _, sequence = read_fasta(reference_fasta)
    slice_start, slice_length = 200_000, 5_000
    piece = list(sequence[slice_start : slice_start + slice_length])
    substitution_at = 1_234
    deletion_at, deletion_length = 3_000, 3
    original = piece[substitution_at]
    piece[substitution_at] = {"A": "C", "C": "A", "G": "T", "T": "G"}.get(original, "A")
    mutated = "".join(piece[:deletion_at] + piece[deletion_at + deletion_length :])

    reference_slice = os.path.join(workdir, "slice.fna")
    with open(reference_slice, "w") as handle:
        handle.write(">slice\n" + sequence[slice_start : slice_start + slice_length] + "\n")
    query = os.path.join(workdir, "mutated.fna")
    with open(query, "w") as handle:
        handle.write(">mutated\n" + mutated + "\n")

    paf = os.path.join(workdir, "self_check.paf")
    run([minimap2, "-c", "-x", "asm5", "--cs=long", reference_slice, query], stdout_path=paf)
    rows = parse_paf(paf)
    aligned, ident, ins, used, introns = alignment_arrays(rows, slice_length)
    not_identical = sorted(np.flatnonzero(~ident).tolist())
    expected_not_identical = [substitution_at] + list(
        range(deletion_at, deletion_at + deletion_length)
    )
    if not_identical != expected_not_identical:
        failures.append(
            "end-to-end cs reader: reference bases not aligned identically are "
            f"{not_identical[:20]} (n={len(not_identical)}), expected exactly "
            f"{expected_not_identical}"
        )
    if ins.any():
        failures.append(
            f"end-to-end cs reader: insertions reported at {np.flatnonzero(ins).tolist()}, "
            "expected none for a substitution-plus-deletion query"
        )

    for line in (
        f"hand-written cs tag: {'ok' if len(failures) == 0 or 'hand-written' not in ' '.join(failures) else 'FAILED'}",
        f"alignment rows walked end to end: {used}, introns met: {introns}",
        f"reference bases not aligned identically: {len(not_identical)}",
    ):
        print(line)
    if failures:
        print("\nSELF-CHECK FAILED")
        for failure in failures:
            print(f"  - {failure}")
        return 1
    print(
        "\nself-check passed: the cs reader reproduces a hand-computed case and "
        "recovers exactly the two mutations planted in a 5 kb slice"
    )
    return 0


# ----------------------------------------------------------------------
# Tables
# ----------------------------------------------------------------------


def fmt(value, digits=4):
    if value is None:
        return "-"
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return f"{value:,}"


def tables(path: str) -> None:
    with open(path) as handle:
        results = json.load(handle)

    print("## Tools\n")
    tools = results["tools"]
    print(f"minimap2 {tools['minimap2']['version']}, paftools.js {tools['paftools.js']['version']}")
    print(
        f"presets measured: {', '.join(tools['presets_measured'])}; variants called from "
        f"{tools['call_preset']}"
    )
    print(f"paftools.js call options: {json.dumps(tools['paftools_call_options'])}\n")

    print("## Aligned fraction of the reference, per preset\n")
    print(
        "| strain | Jaccard vs wMel | preset | primary blocks | median block bp | "
        "longest block bp | reference aligned | callable |"
    )
    print("|---|---|---|---|---|---|---|---|")
    for _key, strain in results["strains"].items():
        for preset, record in strain["presets"].items():
            print(
                f"| {strain['short_name']} | {strain['jaccard_12mer_vs_reference']} | {preset} "
                f"| {record['primary_rows']} | {fmt(record['median_block_bp'], 0)} "
                f"| {fmt(record['longest_block_bp'])} "
                f"| {record['reference_fraction_aligned']:.4f} "
                f"| {record['callable_reference_fraction']:.4f} |"
            )
    print()

    print(f"## Variants called, {results['tools']['call_preset']}\n")
    print(
        "| strain | Jaccard vs wMel | reference aligned | reference identical | callable | "
        "VCF records | SNPs | insertions | deletions |"
    )
    print("|---|---|---|---|---|---|---|---|---|")
    for _key, strain in results["strains"].items():
        counts = strain["variant_counts"]
        print(
            f"| {strain['short_name']} | {strain['jaccard_12mer_vs_reference']} "
            f"| {strain['reference_fraction_aligned']:.4f} "
            f"| {strain['reference_fraction_identical']:.4f} "
            f"| {strain['reference_fraction_callable']:.4f} "
            f"| {fmt(counts['records'])} | {fmt(counts['snps'])} "
            f"| {fmt(counts['insertions'])} | {fmt(counts['deletions'])} |"
        )
    print()

    print("## Control: the same verdict over every reference 12-mer position\n")
    print("| strain | windows inside the callable region | windows surviving by alignment |")
    print("|---|---|---|")
    for _key, strain in results["strains"].items():
        for k, control in (strain.get("uniform_window_controls") or {}).items():
            print(
                f"| {strain['short_name']} (k={k}) "
                f"| {control['fraction_of_all_windows_inside_the_callable_region']:.4f} "
                f"| {control['fraction_of_all_windows_surviving_by_alignment']:.4f} |"
            )
    print()

    for panel in results["panels"]:
        print(
            f"## {panel['design']}, n={panel['requested_size']} "
            f"(delivered {panel['delivered_size']}), {panel['reference_sites']} reference sites, "
            f"reference coverage {panel['reference_coverage']:.4f}\n"
        )
        print(
            "| strain | Jaccard | intact (variant) | surviving (alignment) | "
            "cov ref (variant) | cov ref (alignment) | difference | "
            "sites in strain (genome) | cov strain (genome) |"
        )
        print("|---|---|---|---|---|---|---|---|---|")
        for _key, record in panel["per_strain"].items():
            variant = record["variant_route"]
            homologous = record["homologous_route"]
            genome = record["genome_route"]
            print(
                f"| {record['short_name']} | {record['jaccard_12mer_vs_reference']} "
                f"| {fmt(variant['intact_sites'], 0)} | {fmt(homologous['surviving_sites'])} "
                f"| {fmt(variant['coverage_of_reference'])} "
                f"| {fmt(homologous['coverage_of_reference'])} "
                f"| {fmt(record['like_for_like']['coverage_difference_variant_minus_homologous'])} "
                f"| {fmt(genome['sites_in_strain'], 0)} | {fmt(genome['coverage_of_strain'])} |"
            )
        print()
        print(
            "| strain | sites outside the callable region | of those, surviving by alignment | "
            "of those, called intact by the variant route |"
        )
        print("|---|---|---|---|")
        for _key, record in panel["per_strain"].items():
            unaligned = record["unaligned"]
            print(
                f"| {record['short_name']} "
                f"| {fmt(unaligned['reference_sites_outside_the_callable_region'])} "
                f"| {fmt(unaligned['of_those_surviving_by_alignment'])} "
                f"| {fmt(unaligned['of_those_called_intact_by_the_variant_route'])} |"
            )
        print()

    print("## Coverage difference against divergence, every panel\n")
    print(
        "| strain | Jaccard | panels | difference min | difference max | "
        f"panels above tolerance {results['tolerance']} |"
    )
    print("|---|---|---|---|---|---|")
    for key in results["strains"]:
        differences = [
            panel["per_strain"][key]["like_for_like"][
                "coverage_difference_variant_minus_homologous"
            ]
            for panel in results["panels"]
            if panel["per_strain"].get(key)
        ]
        differences = [value for value in differences if value is not None]
        if not differences:
            continue
        above = sum(1 for value in differences if abs(value) > results["tolerance"])
        print(
            f"| {SHORT_NAME[key]} | {JACCARD_VS_REFERENCE[key]} | {len(differences)} "
            f"| {min(differences):.4f} | {max(differences):.4f} | {above} |"
        )
    print()

    print("## Which strain each route calls the worst, and where it puts wBm\n")
    keys = list(results["strains"])
    print(
        "| design | n | worst by variant route | worst by alignment | worst by genome route | "
        "wBm rank, variant | wBm rank, alignment | wBm rank, genome |"
    )
    print("|---|---|---|---|---|---|---|---|")
    agreements = 0
    for panel in results["panels"]:
        by_route = {
            "variant": {
                key: panel["per_strain"][key]["variant_route"]["coverage_of_reference"]
                for key in keys
            },
            "alignment": {
                key: panel["per_strain"][key]["homologous_route"]["coverage_of_reference"]
                for key in keys
            },
            "genome": {
                key: panel["per_strain"][key]["genome_route"]["coverage_of_strain"] for key in keys
            },
        }
        worst = {}
        rank = {}
        for route, figures in by_route.items():
            usable = {key: value for key, value in figures.items() if value is not None}
            worst[route] = SHORT_NAME[min(usable, key=usable.get)] if usable else "-"
            order = sorted(usable, key=usable.get, reverse=True)
            rank[route] = order.index("wolbachia_wbm") + 1 if "wolbachia_wbm" in order else None
        if worst["variant"] == worst["genome"]:
            agreements += 1
        print(
            f"| {panel['design']} | {panel['requested_size']} | {worst['variant']} "
            f"| {worst['alignment']} | {worst['genome']} | {rank['variant']} "
            f"| {rank['alignment']} | {rank['genome']} |"
        )
    print(
        f"\nThe variant route and the genome route name the same worst strain in "
        f"{agreements} of {len(results['panels'])} panels.\n"
    )

    print("## Genome-route recomputation against Phase 3\n")
    print(
        "| design | n | strain | sites recomputed | sites Phase 3 | coverage recomputed | "
        "coverage Phase 3 | difference |"
    )
    print("|---|---|---|---|---|---|---|---|")
    for check in results["genome_route_checks"]:
        for _key, record in check["per_target"].items():
            print(
                f"| {check['design']} | {check['requested_size']} | {SHORT_NAME.get(key, key)} "
                f"| {fmt(record['sites_recomputed'], 0)} | {fmt(record['sites_phase3'], 0)} "
                f"| {fmt(record['coverage_recomputed'])} | {fmt(record['coverage_phase3'])} "
                f"| {fmt(record['coverage_difference'], 6)} |"
            )
    print()


# ----------------------------------------------------------------------


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--minimap2", help="path to the minimap2 executable")
    parser.add_argument("--paftools", help="path to paftools.js")
    parser.add_argument("--k8", help="path to the k8 interpreter paftools.js needs")
    parser.add_argument(
        "--genomes",
        default=os.path.join(REPO, "tests", "validation", "genomes"),
        help="directory holding the five Wolbachia FASTAs",
    )
    parser.add_argument(
        "--results",
        default=os.path.join(
            REPO, "tests", "validation", "genomes", "diversity_baseline", "results.json"
        ),
        help="Phase 3's results.json, read only, for the delivered panels",
    )
    parser.add_argument(
        "--workdir",
        default=os.path.join(REPO, "tests", "validation", "genomes", "variant_route_agreement"),
        help="where alignments, VCFs and the results file are written (gitignored)",
    )
    parser.add_argument("--output", help="results JSON (default: <workdir>/results.json)")
    parser.add_argument("--presets", default=",".join(PRESETS), help="minimap2 presets to measure")
    parser.add_argument(
        "--call-preset", default="asm20", help="the preset the variant table is called from"
    )
    parser.add_argument("--groups", default=",".join(DEFAULT_GROUPS), help="Phase 3 design groups")
    parser.add_argument("--sizes", default="6,12,24", help="requested panel sizes")
    parser.add_argument(
        "--genome-route-sample",
        type=int,
        default=2,
        help="how many panels to recompute through evaluate_reference_panel, to "
        "confirm the genome-route figures quoted from Phase 3",
    )
    parser.add_argument("--self-check", action="store_true", help="verify the cs reader and stop")
    parser.add_argument("--tables", help="print the record's tables from a finished results JSON")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="%(asctime)s %(levelname)s %(message)s",
    )

    if args.tables:
        tables(args.tables)
        return 0

    args.presets = tuple(value.strip() for value in args.presets.split(",") if value.strip())
    args.groups = tuple(value.strip() for value in args.groups.split(",") if value.strip())
    args.sizes = tuple(int(value) for value in args.sizes.split(",") if value.strip())

    if args.self_check:
        return self_check(args)

    try:
        results = measure(args)
    except ToolMissing as exc:
        print(f"stopping: {exc}", file=sys.stderr)
        return 2

    output = args.output or os.path.join(args.workdir, "results.json")
    os.makedirs(os.path.dirname(os.path.abspath(output)), exist_ok=True)
    with open(output, "w") as handle:
        json.dump(results, handle, indent=1, default=str)
    print(f"wrote {output}")
    tables(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
