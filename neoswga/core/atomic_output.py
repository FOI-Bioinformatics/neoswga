"""Publish a result file so a partial write cannot look like a complete one.

An output directory is the only thing a later command sees, and nothing in it
carries a timestamp anybody compares. `export` and the report both read these
files back and neither can tell a truncated one from a short one: a summary
JSON cut off mid-object fails to parse, which the report treats as "no
optimizer summary" and falls back to its own estimate, and a results CSV cut
off mid-row is a shorter panel that parses cleanly.

So the rule is the one `position_index` already follows for the same reason:
write beside the target and rename into place. `os.replace` is atomic on POSIX
and on Windows, so a reader sees either the previous file or the new one, never
half of either. A serialisation error then leaves the previous file intact
rather than destroying it, which matters most for exactly the run that failed.

This is deliberately NOT a general file utility. It is for the artifacts a
later command consumes as authoritative.
"""

from __future__ import annotations

import contextlib
import json
import os
from collections.abc import Callable
from typing import TextIO


def atomic_write_text(path: str, render: Callable[[TextIO], None]) -> None:
    """Write `path` by rendering into a temporary file and renaming it.

    `render` receives an open text handle. It may raise: the target is left as
    it was, and the temporary file is removed. The process id is in the
    temporary name so two runs sharing a directory cannot scribble on each
    other's partial file -- they have a provenance problem either way, and this
    at least keeps it from corrupting the survivor.
    """
    # The directory is NOT created. A write into a path whose directory does
    # not exist is a caller error and must say so: inventing the location is
    # how a failure record once left a DIRECTORY where the user had asked for a
    # file, and `tests/test_a_missing_verdict_is_not_a_passing_one.py` pins
    # that the validation record's own write is not swallowed.
    directory = os.path.dirname(os.path.abspath(path)) or "."
    tmp = os.path.join(directory, f".{os.path.basename(path)}.{os.getpid()}.tmp")
    try:
        with open(tmp, "w", encoding="utf-8") as handle:
            render(handle)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(tmp, path)
    finally:
        with contextlib.suppress(FileNotFoundError):
            os.remove(tmp)


def atomic_write_json(path: str, payload, **dump_kwargs) -> None:
    """Serialise `payload` to `path` atomically."""
    atomic_write_text(path, lambda handle: json.dump(payload, handle, **dump_kwargs))


def atomic_write_dataframe(path: str, frame, **to_csv_kwargs) -> None:
    """Write a pandas frame to `path` atomically."""
    atomic_write_text(path, lambda handle: frame.to_csv(handle, **to_csv_kwargs))
