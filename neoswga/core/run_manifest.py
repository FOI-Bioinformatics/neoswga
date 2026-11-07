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
from datetime import UTC, datetime, timezone
from pathlib import Path
from typing import Any

try:
    import fcntl

    HAS_FCNTL = True
except ImportError:  # pragma: no cover - Windows
    HAS_FCNTL = False

logger = logging.getLogger(__name__)

MANIFEST_FILENAME = "run_manifest.json"


def _sha256(path: str) -> str | None:
    try:
        h = hashlib.sha256()
        with open(path, "rb") as f:
            for chunk in iter(lambda: f.read(65536), b""):
                h.update(chunk)
        return h.hexdigest()
    except OSError as e:
        logger.debug(f"sha256 failed for {path}: {e}")
        return None


def _git_sha() -> str | None:
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


def _jellyfish_version() -> str | None:
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
    data_dir: str | None,
    params_path: str | None = None,
    resolved_params: dict[str, Any] | None = None,
    input_files: list[str] | None = None,
    output_files: list[str] | None = None,
    seed: int | None = None,
    extra: dict[str, Any] | None = None,
    effective_conditions: dict[str, Any] | None = None,
) -> str | None:
    """Append a step entry to ``<data_dir>/run_manifest.json``.

    Args:
        step: Step identifier (e.g. ``"count-kmers"``, ``"filter"``).
        data_dir: Directory to write the manifest into. ``None`` skips writing.
        params_path: Path to the user's params.json.
        resolved_params: Resolved parameter dict; loaded from ``params_path``
            if not provided.
        input_files: Paths to checksum. Missing paths are skipped silently.
        output_files: Paths the step wrote, checksummed into
            ``output_checksums``. Separate from ``input_files`` because the
            step handlers used to pass their output there, which is what made
            three ``optimize`` entries under two git SHAs all hash-match the
            single result file beside them.
        seed: Resolved RNG seed used for this step.
        extra: Step-specific fields to merge into the entry.
        effective_conditions: The reaction conditions this step actually ran
            under, after any run-time adjustment. ``resolved_params`` is a copy
            of params.json, so on its own it does not describe the reaction:
            the GC-adaptive strategy sets betaine and DMSO from the target's GC
            and neither appears in the file. Recording them is what lets
            ``export`` and ``report`` state a Tm corrected for the buffer the
            design was optimized under.

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

    input_checksums: dict[str, str | None] = {}
    for path in input_files or []:
        if path and os.path.exists(path):
            input_checksums[path] = _sha256(path)

    output_checksums: dict[str, str | None] = {}
    for path in output_files or []:
        if path and os.path.exists(path):
            output_checksums[path] = _sha256(path)

    entry: dict[str, Any] = {
        "step": step,
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "neoswga_version": _neoswga_version(),
        "git_sha": _git_sha(),
        "jellyfish_version": _jellyfish_version(),
        "python_version": sys.version.split()[0],
        "platform": sys.platform,
        "cli_invocation": list(sys.argv),
        "seed": seed,
        "params_path": params_path,
        "resolved_params": resolved_params,
        "effective_conditions": effective_conditions,
        "input_checksums": input_checksums,
        "output_checksums": output_checksums,
    }
    if extra:
        entry["extra"] = extra

    # Read, append and rewrite under an exclusive lock. Without it two runs
    # sharing a data_dir both read the same list and both write their own
    # append, losing one entry. The lock is held on the manifest file itself
    # for the whole read-modify-write; ``experimental_tracker.py`` uses the
    # same fcntl.flock shape. On a platform without fcntl this degrades to the
    # previous unlocked behaviour rather than failing.
    try:
        with open(manifest_path, "a+") as f:
            if HAS_FCNTL:
                fcntl.flock(f.fileno(), fcntl.LOCK_EX)
            try:
                # "a+" positions at end of file, so rewind before reading.
                f.seek(0)
                raw = f.read()
                existing: dict[str, Any] = {"steps": []}
                if raw.strip():
                    try:
                        loaded = json.loads(raw)
                        if isinstance(loaded, dict) and isinstance(loaded.get("steps"), list):
                            existing = loaded
                    except json.JSONDecodeError:
                        pass

                existing["steps"].append(entry)

                f.seek(0)
                f.truncate()
                json.dump(existing, f, indent=2, sort_keys=True)
                f.flush()
                os.fsync(f.fileno())
            finally:
                if HAS_FCNTL:
                    fcntl.flock(f.fileno(), fcntl.LOCK_UN)
    except OSError as e:
        logger.warning(f"Could not write run_manifest.json: {e}")
        return None

    logger.info(f"Run manifest updated: {manifest_path}")
    return manifest_path


def read_effective_conditions(data_dir: str, step: str | None = None) -> dict[str, Any] | None:
    """The reaction conditions a recorded step ran under, if any.

    `export` and `report` otherwise reconstruct conditions from params.json,
    which omits anything the run decided for itself -- so they reported a Tm
    corrected for a different buffer than the one the design was optimized
    under.

    With no ``step``, the latest entry carrying conditions wins, which is the
    historical behaviour and what `export` and `report` use. Pass a step name
    to read that step's own reaction: the manifest is append-only and a rerun
    interleaves steps, so in one real run two ``score`` entries sat after an
    ``optimize`` entry and the unfiltered read described the optimize result
    under the score step's reaction.
    """
    manifest_path = os.path.join(data_dir, MANIFEST_FILENAME)
    if not os.path.exists(manifest_path):
        return None
    try:
        with open(manifest_path) as f:
            data = json.load(f)
    except (OSError, json.JSONDecodeError) as e:
        logger.debug(f"Could not read manifest for conditions: {e}")
        return None

    steps = data.get("steps") if isinstance(data, dict) else None
    if not isinstance(steps, list):
        return None
    for entry in reversed(steps):
        if not isinstance(entry, dict):
            continue
        if step is not None and entry.get("step") != step:
            continue
        if isinstance(entry.get("effective_conditions"), dict):
            return entry["effective_conditions"]
    return None


#: The date the configured polymerase and reaction temperature began reaching
#: the optimizer. Before it, `unified_optimizer` passed neither to any
#: optimizer -- it mentions neither name -- and both hybrid wrappers defaulted
#: `polymerase` to the literal string "phi29", so every `hybrid` and
#: `background-aware` run built a phi29 preset at 30 C whatever params.json
#: declared. The fix landed on 2026-10-26; this boundary is the day after, so a
#: run stamped on the 26th itself is reported rather than cleared. Over-
#: reporting is the safe direction for a provenance warning.
#: `docs/validation/2026-10-26-hybrid-ran-every-reaction-as-phi29.md`.
REACTION_REACHED_OPTIMIZER = datetime(2026, 10, 27, tzinfo=UTC)

#: The reaction every pre-fix optimizer run actually used, whatever was asked
#: for: phi29's preset temperature.
_PRE_FIX_TEMP = 30.0

#: Steps whose delivered panel the optimizer chose. `count-kmers`, `filter` and
#: `prepare-candidates` resolved their own Tm window from `conditions` and are
#: not affected by this defect.
_OPTIMIZER_STEPS = ("optimize", "expand-primers")


def steps_predating_the_reaction_fix(data_dir: str) -> list[dict[str, Any]]:
    """Manifest steps whose panel was selected before the reaction reached it.

    One entry per affected step, newest last, each naming what the step asked
    for and what it got. An empty list means this directory holds nothing
    affected -- including the common case of a phi29 design at 30 C, which the
    defect left correct.

    Two distinct consequences are reported separately, because they differ in
    severity:

    `screen_skipped` -- the polymerase's preset sets ``thermo_filter=True``
    (equiphi29 and bst) while phi29's does not, and `optimize` gates its whole
    Stage-0 thermodynamic pre-screen on that flag. So the screen that rejects
    candidates on hairpin, homodimer and heterodimer energies and on the Tm
    window did not run at all. Measured under the shipped
    `examples/equiphi29_scenario` reaction, that screen keeps 299 of 500
    candidates, so 500/500 reached selection unscreened.

    `temp_used` against `temp_intended` -- Stage 2's dimer screen and its Tm
    term both read the reaction temperature, so both ran at 30 C. This bites
    any design whose reaction was not 30 C, phi29 included: an explicitly
    configured `reaction_temp` did not reach the optimizer either.

    The boundary is a timestamp, not an ancestry check: the manifest records
    `git_sha`, but resolving whether a SHA predates the fix needs the
    repository, and a date needs only the manifest.
    """
    from neoswga.core.registry import POLYMERASES

    manifest_path = os.path.join(data_dir or ".", MANIFEST_FILENAME)
    if not os.path.exists(manifest_path):
        return []
    try:
        with open(manifest_path) as f:
            data = json.load(f)
    except (OSError, json.JSONDecodeError) as e:
        logger.debug(f"Could not read manifest for reaction provenance: {e}")
        return []

    steps = data.get("steps") if isinstance(data, dict) else None
    if not isinstance(steps, list):
        return []

    affected: list[dict[str, Any]] = []
    for entry in steps:
        if not isinstance(entry, dict) or entry.get("step") not in _OPTIMIZER_STEPS:
            continue
        stamp = _parse_timestamp(entry.get("timestamp_utc"))
        if stamp is None or stamp >= REACTION_REACHED_OPTIMIZER:
            continue

        conditions = entry.get("effective_conditions")
        conditions = conditions if isinstance(conditions, dict) else {}
        # An absent polymerase means the default, which is the one case the
        # defect left correct -- so absence must read as phi29, not as unknown.
        polymerase = conditions.get("polymerase") or "phi29"
        spec = POLYMERASES.get(polymerase)
        intended = conditions.get("reaction_temp")
        if intended is None:
            intended = getattr(spec, "preset_reaction_temp", None)

        screen_skipped = bool(getattr(spec, "thermo_filter", False))
        temp_wrong = intended is not None and float(intended) != _PRE_FIX_TEMP
        if not (screen_skipped or temp_wrong):
            continue

        affected.append(
            {
                "step": entry.get("step"),
                "timestamp_utc": entry.get("timestamp_utc"),
                "polymerase": polymerase,
                "temp_intended": intended,
                "temp_used": _PRE_FIX_TEMP,
                "screen_skipped": screen_skipped,
            }
        )
    return affected


def _parse_timestamp(value: Any) -> datetime | None:
    """A manifest timestamp as an aware datetime, or None if unreadable.

    An entry whose stamp cannot be parsed is left out of the provenance report
    rather than assumed old: a fabricated date would be a fabricated finding.
    """
    if not isinstance(value, str):
        return None
    try:
        parsed = datetime.fromisoformat(value)
    except ValueError:
        return None
    return parsed if parsed.tzinfo else parsed.replace(tzinfo=UTC)
