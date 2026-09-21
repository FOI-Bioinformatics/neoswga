"""Saying so when two runs are sharing one data directory.

HDF5 takes a file lock. When a second neoswga process writes position files
into a directory a first one is already writing, the loser raises

    BlockingIOError: [Errno 35] Unable to synchronously open file
    (unable to lock file, errno = 35, ...)

and nothing in the package translated it. What a user saw was an h5py
traceback naming `h5py/h5f.pyx`, with no hint that another process was the
cause or that retrying alone would fix it.

**Measured 2026-09-21, and it settled a misattribution.**
`tests/validation/genomes/f_mixed.log` records this failure from a
mixed-length run, and mixed length was the suspect because that is what the
config changed. It is not the cause.

Four things, each measured:

- **This errno is cross-process by construction.** A foreign process holding
  the file, even READ-ONLY, produces exactly the recorded message. The
  same-process case -- one handle open for reading, then an `r+` open --
  produces a different one, `OSError: ... file is already open for read-only`,
  with no errno 35. So a leaked handle inside one process cannot be the cause
  of THIS message.
- **A reader is enough to stop a writer**, which is the shape that actually
  bites: `optimize` holds read handles open for a whole run through
  `StreamingPositionCache`, so an `optimize` and a `filter` on one directory
  collide even though only one of them writes.
- **Two runs reproduce it and one does not.** Concurrent pairs failed 6 of 7
  and 4 of 4 in two separate trials, at the same frame as the recorded log.
  Single-process runs failed 0 of 6 and 0 of 5, including multi-k ones over
  the same references at realistic scale.
- **The recorded directory shows the collision directly.** Its
  `run_manifest.json` records `score` on that same config completing 2.2 s
  before the failing filter's last log write and `optimize` completing 8.7 s
  after, with no `filter` entry at all, because the manifest is written on
  completion. A `prevotella_13mer_positions.h5` of 800 bytes holding zero
  datasets sits there with the same mtime, for a k that run never requested,
  which is a third process still running under a config edited minutes
  earlier.

The recorded traceback differs from the reproduction only in `h5f.open`
against `h5f.create`, which is whether the target file already existed.

No lock is taken here and no retry is attempted. Serialising runs is the
user's call: a retry loop would hide a genuine second process rather than
report it, and two designs writing one directory have a provenance problem
that outlives the lock.
"""

from __future__ import annotations

import errno

__all__ = ["locked_file_advice"]

_LOCK_MARKERS = ("unable to lock file", "resource temporarily unavailable")


def locked_file_advice(exc: BaseException, data_dir=None) -> str:
    """Advice for an HDF5 lock failure, or an empty string for anything else.

    Matches on the errno AND on the message, because `BlockingIOError` with
    EAGAIN is raised by several unrelated things and a message match alone
    would claim a concurrent run for any error mentioning a lock.
    """
    if not isinstance(exc, OSError):
        return ""
    if exc.errno not in (errno.EAGAIN, errno.EWOULDBLOCK):
        return ""
    if not any(marker in str(exc).lower() for marker in _LOCK_MARKERS):
        return ""

    where = f" in {data_dir}" if data_dir else ""
    return (
        f"A position file{where} is locked by another process. A data "
        "directory takes one neoswga run at a time. The other run does not "
        "have to be writing: `optimize` keeps these files open for READING "
        "for its whole duration, and that alone stops a `filter` writing "
        "them, so an optimize and a filter on one directory collide. "
        "Measured across processes only -- a second handle inside one process "
        "reports something else entirely -- and reproduced in 10 of 11 "
        "concurrent attempts against 0 of 11 single-process runs. "
        "Wait for the other run to finish, or give this one its own data_dir. "
        "Nothing was written for the file that failed, so re-running this step "
        "afterwards is safe. If no other neoswga process is running, a stale "
        "lock can persist on a network filesystem, where HDF5 locking is "
        "unreliable; copy the directory to local storage rather than setting "
        "HDF5_USE_FILE_LOCKING=FALSE, which risks a corrupt index."
    )
