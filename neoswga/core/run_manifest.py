"""Run manifest writer for reproducibility.

Captures version, git SHA, jellyfish version, seed, input checksums, resolved
params, CLI invocation, and UTC timestamp into ``<data_dir>/run_manifest.json``
after each pipeline step.

Multiple step runs append to the same manifest file as a list under the
``steps`` key so reruns can be traced. Failures to write the manifest are
logged but do not raise: reproducibility metadata should never abort a
pipeline that otherwise succeeded.
"""

from __future__ import annotations

import hashlib
import json
import logging
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)

MANIFEST_FILENAME = "run_manifest.json"


def _sha256(path: str) -> Optional[str]:
    try:
        h = hashlib.sha256()
        with open(path, "rb") as f:
            for chunk in iter(lambda: f.read(65536), b""):
                h.update(chunk)
        return h.hexdigest()
    except OSError as e:
        logger.debug(f"sha256 failed for {path}: {e}")
        return None


def _git_sha() -> Optional[str]:
    try:
        out = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=Path(__file__).resolve().parent,
            check=True,
            capture_output=True,
            text=True,
            timeout=5,
        )
        return out.stdout.strip() or None
    except (subprocess.SubprocessError, OSError, FileNotFoundError):
        return None


def _jellyfish_version() -> Optional[str]:
    try:
        out = subprocess.run(
            ["jellyfish", "--version"],
            check=True,
            capture_output=True,
            text=True,
            timeout=5,
        )
        return out.stdout.strip() or None
    except (subprocess.SubprocessError, OSError, FileNotFoundError):
        return None


def _neoswga_version() -> str:
    try:
        from neoswga import __version__
        return __version__
    except Exception:
        return "unknown"


def write_manifest(
    step: str,
    data_dir: Optional[str],
    params_path: Optional[str] = None,
    resolved_params: Optional[Dict[str, Any]] = None,
    input_files: Optional[List[str]] = None,
    seed: Optional[int] = None,
    extra: Optional[Dict[str, Any]] = None,
) -> Optional[str]:
    """Append a step entry to ``<data_dir>/run_manifest.json``.

    Args:
        step: Step identifier (e.g. ``"count-kmers"``, ``"filter"``).
        data_dir: Directory to write the manifest into. ``None`` skips writing.
        params_path: Path to the user's params.json.
        resolved_params: Resolved parameter dict; loaded from ``params_path``
            if not provided.
        input_files: Paths to checksum. Missing paths are skipped silently.
        seed: Resolved RNG seed used for this step.
        extra: Step-specific fields to merge into the entry.

    Returns:
        The manifest path on success, or ``None`` if writing was skipped or
        failed.
    """
    if not data_dir:
        return None
    try:
        os.makedirs(data_dir, exist_ok=True)
    except OSError as e:
        logger.warning(f"Could not create data_dir for manifest: {e}")
        return None

    manifest_path = os.path.join(data_dir, MANIFEST_FILENAME)

    if params_path and resolved_params is None and os.path.exists(params_path):
        try:
            with open(params_path) as f:
                resolved_params = json.load(f)
        except (OSError, json.JSONDecodeError) as e:
            logger.debug(f"Could not load params for manifest: {e}")

    input_checksums: Dict[str, Optional[str]] = {}
    for path in input_files or []:
        if path and os.path.exists(path):
            input_checksums[path] = _sha256(path)

    entry: Dict[str, Any] = {
        "step": step,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "neoswga_version": _neoswga_version(),
        "git_sha": _git_sha(),
        "jellyfish_version": _jellyfish_version(),
        "python_version": sys.version.split()[0],
        "platform": sys.platform,
        "cli_invocation": list(sys.argv),
        "seed": seed,
        "params_path": params_path,
        "resolved_params": resolved_params,
        "input_checksums": input_checksums,
    }
    if extra:
        entry["extra"] = extra

    existing: Dict[str, Any] = {"steps": []}
    if os.path.exists(manifest_path):
        try:
            with open(manifest_path) as f:
                loaded = json.load(f)
            if isinstance(loaded, dict) and isinstance(loaded.get("steps"), list):
                existing = loaded
        except (OSError, json.JSONDecodeError):
            pass

    existing["steps"].append(entry)

    try:
        with open(manifest_path, "w") as f:
            json.dump(existing, f, indent=2, sort_keys=True)
    except OSError as e:
        logger.warning(f"Could not write run_manifest.json: {e}")
        return None

    logger.info(f"Run manifest updated: {manifest_path}")
    return manifest_path
