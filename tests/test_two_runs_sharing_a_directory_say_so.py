"""An HDF5 lock failure names the cause instead of showing an h5py traceback.

`tests/validation/genomes/f_mixed.log` records a `BlockingIOError: [Errno 35]
unable to lock file` from a mixed-length run, and mixed length was blamed
because that is what the config changed.

**It is not the cause.** Measured 2026-09-21: two neoswga processes against one
data directory reproduced the identical error in 6 of 7 attempts, while 6 of 6
single-process runs -- including a multi-k run over the same references --
completed cleanly. The recorded traceback differs from the reproduction only in
`h5f.open` against `h5f.create`, which is whether the target file existed yet.

So the remedy is a message, not a code change to the scan. Nothing here takes a
lock or retries: a retry loop would hide a genuine second process, and two
designs writing one directory have a provenance problem that outlasts the lock.

These tests pin the predicate rather than drive a real lock. Two processes
racing for an HDF5 file is not something a unit test should stage, and the
predicate is where the mistakes live -- matching on the message alone would
claim a concurrent run for any error that mentions a lock.
"""

import errno

import pytest

from neoswga.core.concurrent_runs import locked_file_advice

LOCK_MESSAGE = (
    "[Errno 35] Unable to synchronously open file (unable to lock file, "
    "errno = 35, error message = 'Resource temporarily unavailable')"
)


def lock_error(message=LOCK_MESSAGE, code=errno.EAGAIN):
    return BlockingIOError(code, message)


# ---------------------------------------------------------------------------
# It fires on the real thing
# ---------------------------------------------------------------------------


def test_the_recorded_failure_is_recognised():
    advice = locked_file_advice(lock_error())

    assert advice
    assert "another process" in advice


def test_the_create_variant_is_recognised_too():
    """The reproduction hit `h5f.create`; the recorded log hit `h5f.open`."""
    advice = locked_file_advice(
        lock_error("[Errno 35] Unable to synchronously create file (unable to lock file)")
    )

    assert advice


def test_the_advice_names_the_directory_when_it_is_known():
    advice = locked_file_advice(lock_error(), data_dir="/tmp/design")

    assert "/tmp/design" in advice


def test_the_advice_names_the_reader_blocks_writer_shape():
    """The surprising half, and the one a user will actually hit.

    Measured: a foreign process holding the file READ-ONLY produces exactly
    the recorded error, so a command that only reads the index can stop a
    `filter`. A message saying only "another writer" would send someone
    looking for a second filter.
    """
    advice = locked_file_advice(lock_error())

    assert "reading" in advice.lower()
    assert "WRITING" in advice, "the point is that the other run need not write"


def test_the_advice_does_not_overstate_how_long_a_reader_holds_the_file():
    """A first version claimed `optimize` holds these open for its whole run.

    That is true only of `StreamingPositionCache`, which is selected when the
    in-memory cache is DISABLED. The default `PositionCache` opens each file
    inside a `with` block and closes it promptly, so on the default path the
    collision is a race rather than a certainty. Naming one command flatly
    told a user something untrue about their own run.
    """
    advice = locked_file_advice(lock_error())

    assert "brief" in advice and "race" in advice


def test_the_two_caches_really_do_differ_in_how_long_they_hold_a_handle():
    """The distinction the message now rests on, asserted against the source.

    Driving this would mean opening real HDF5 files and inspecting handle
    lifetimes, which is what the caches exist to manage; what matters here is
    only that one keeps a handle dictionary and the other does not.
    """
    import inspect

    from neoswga.core.position_cache import PositionCache, StreamingPositionCache

    streaming = inspect.getsource(StreamingPositionCache)
    assert "self.file_handles" in streaming, (
        "StreamingPositionCache no longer retains handles; if it now closes "
        "promptly, the lock message's 'whole run' case is stale"
    )
    assert "def close" in streaming

    assert "file_handles" not in inspect.getsource(PositionCache), (
        "PositionCache now retains handles too, so the default path holds "
        "these files for a whole run and the lock message understates it"
    )


def test_the_advice_says_rerunning_is_safe():
    """The one thing a user needs to know before touching anything."""
    advice = locked_file_advice(lock_error())

    assert "re-running" in advice


def test_the_advice_does_not_recommend_disabling_locking():
    """HDF5_USE_FILE_LOCKING=FALSE is the first hit for this error online.

    It risks a corrupt index, which is worse than the failure it silences, so
    the message must mention it only to steer away from it.
    """
    advice = locked_file_advice(lock_error())

    assert "HDF5_USE_FILE_LOCKING" in advice
    index = advice.index("HDF5_USE_FILE_LOCKING")
    assert "rather than" in advice[:index], "it must be named as the thing NOT to do"


# ---------------------------------------------------------------------------
# It does not fire on anything else
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "exc",
    [
        ValueError("candidates list cannot be empty"),
        FileNotFoundError(errno.ENOENT, "no such file"),
        OSError(errno.ENOSPC, "No space left on device"),
        MemoryError(),
    ],
)
def test_an_unrelated_failure_gets_no_advice(exc):
    assert locked_file_advice(exc) == ""


def test_the_right_errno_with_the_wrong_message_gets_no_advice():
    """EAGAIN is raised by plenty of things that are not a file lock."""
    assert locked_file_advice(BlockingIOError(errno.EAGAIN, "socket would block")) == ""


def test_the_right_message_with_the_wrong_errno_gets_no_advice():
    """Matching on text alone would claim a concurrent run for anything."""
    assert locked_file_advice(OSError(errno.EIO, "unable to lock file")) == ""


# ---------------------------------------------------------------------------
# It reaches the user
# ---------------------------------------------------------------------------


def test_the_shared_step_handler_prints_the_advice_and_exits():
    """Driven, not read off a source location.

    An earlier version of this test asserted that `locked_file_advice` appears
    by name inside `cli/pipeline.py`. It broke the moment the three duplicated
    except blocks were unified into one handler, while the advice still
    reached the user -- the same failure `attach_search_config` records, and
    the same one that broke the dimer-floor test on 2026-09-21.
    """
    import logging as stdlib_logging

    from neoswga.cli._failure import exit_on_step_failure

    logger = stdlib_logging.getLogger("test_lock_advice")
    records = []
    logger.error = lambda message, *a: records.append(str(message))

    with pytest.raises(SystemExit) as excinfo:
        exit_on_step_failure("Step 2", lock_error(), logger, data_dir="/tmp/design")

    assert excinfo.value.code != 0, "a failed step must exit nonzero"
    assert any("locked by another process" in line for line in records), records
    assert any("/tmp/design" in line for line in records), records


def test_an_unrelated_step_failure_gets_no_lock_advice():
    """The handler must not attach this explanation to every failure."""
    import logging as stdlib_logging

    from neoswga.cli._failure import exit_on_step_failure

    logger = stdlib_logging.getLogger("test_lock_advice")
    records = []
    logger.error = lambda message, *a: records.append(str(message))

    with pytest.raises(SystemExit):
        exit_on_step_failure("Step 2", ValueError("pool is empty"), logger)

    assert not any("locked by another process" in line for line in records), records


def test_every_pipeline_step_routes_through_that_handler():
    """Steps 2, 3 and 4 all write HDF5, so all three must use it.

    They each carried this block inline and had drifted: before the handler
    existed, only one of them explained a lock collision.
    """
    import ast
    import pathlib

    source = pathlib.Path(__file__).resolve().parent.parent / "neoswga" / "cli" / "pipeline.py"
    tree = ast.parse(source.read_text(encoding="utf-8"))

    calls = [
        node
        for node in ast.walk(tree)
        if isinstance(node, ast.Call)
        and getattr(node.func, "id", None) == "exit_on_step_failure"
    ]

    assert len(calls) >= 3, f"found {len(calls)} step(s) using the shared handler"
