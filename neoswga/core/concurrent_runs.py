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
config changed. It is not the cause. Two neoswga processes against one
directory reproduced the identical error in 6 of 7 attempts, while 6 of 6
single-process runs -- including a multi-k one over the same references --
completed cleanly. The recorded traceback differs from the reproduction only
in `h5f.open` against `h5f.create`, which is whether the target file already
existed.

That directory holds about 60 run logs, which is what made a concurrent run
the likelier explanation once it was tested rather than assumed.

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
        f"A position file{where} is locked by another process. This is what a "
        "second neoswga run against the same data_dir looks like: HDF5 takes a "
        "file lock and the second writer is refused. Measured in 6 of 7 "
        "concurrent attempts, and in none of 6 single-process runs. "
        "Wait for the other run to finish, or give this one its own data_dir. "
        "Nothing was written for the file that failed, so re-running this step "
        "afterwards is safe. If no other neoswga process is running, a stale "
        "lock can persist on a network filesystem, where HDF5 locking is "
        "unreliable; copy the directory to local storage rather than setting "
        "HDF5_USE_FILE_LOCKING=FALSE, which risks a corrupt index."
    )
