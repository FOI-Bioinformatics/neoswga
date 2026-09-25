"""Cache for the foreground GC fraction.

`_apply_gc_adaptive_defaults` needs the target's GC content, so
`parameter.get_params` loaded every foreground FASTA and counted G and C. That
read costs about 0.30 s on a 4.6 Mb genome and scales with target size, and it
happens once per pipeline step: four times for a value that cannot change
between steps of one run.

The obvious alternative -- skipping the derivation on the steps that do not
apply it -- was rejected. The derivation is what fills the manifest's
`effective_conditions`, and `optimize` now warns when its conditions differ from
the filter step's, so a `score` step that skipped it would make the default
pipeline warn about itself. Caching keeps every step's recorded conditions
identical.

The key is each foreground file's absolute path, byte size and modification
time in nanoseconds. Any mismatch, a missing file or an unreadable cache
recomputes.
"""

from __future__ import annotations

import hashlib
import json
import logging
import os

logger = logging.getLogger(__name__)

CACHE_FILENAME = "genome_gc.json"


def cache_key(fg_genomes: list[str]) -> str | None:
    """A key over the foreground file list, or None if any file is missing."""
    parts = []
    for path in fg_genomes or []:
        try:
            st = os.stat(path)
        except OSError:
            return None
        parts.append(f"{os.path.abspath(path)}:{st.st_size}:{st.st_mtime_ns}")
    if not parts:
        return None
    return hashlib.sha256("|".join(parts).encode("utf-8")).hexdigest()


def read_cached_gc(data_dir: str | None, fg_genomes: list[str]) -> float | None:
    """The cached GC fraction for this exact set of files, or None."""
    if not data_dir:
        return None
    key = cache_key(fg_genomes)
    if key is None:
        return None
    path = os.path.join(data_dir, CACHE_FILENAME)
    if not os.path.exists(path):
        return None
    try:
        with open(path) as fh:
            data = json.load(fh)
    except (OSError, json.JSONDecodeError) as e:
        logger.debug(f"Ignored unreadable genome GC cache: {e}")
        return None
    if not isinstance(data, dict) or data.get("key") != key:
        return None
    value = data.get("genome_gc")
    return value if isinstance(value, (int, float)) else None


def write_cached_gc(data_dir: str | None, fg_genomes: list[str], genome_gc: float) -> None:
    """Record the GC fraction. Best-effort: never fails the caller."""
    if not data_dir:
        return
    key = cache_key(fg_genomes)
    if key is None:
        return
    try:
        os.makedirs(data_dir, exist_ok=True)
        with open(os.path.join(data_dir, CACHE_FILENAME), "w") as fh:
            json.dump({"key": key, "genome_gc": genome_gc}, fh, indent=2)
    except OSError as e:
        logger.debug(f"Could not write genome GC cache: {e}")
