"""Which MILP backend python-mip is handed.

python-mip's default is CBC. On Python 3.13 (measured on macOS arm64 with mip
2.0.0 and cbcbox 2.935) constructing a CBC model terminates the interpreter
with SIGKILL: no exception, no traceback, no stderr. The same versions on
Python 3.11 build and solve the same models in a third of a second.

mip 2.0 ships a second backend. HiGHS solves those models on 3.13 and on 3.11,
so it is preferred wherever its runtime is importable.

**Without HiGHS this refuses rather than substituting CBC.** It used to warn and
fall back, and a warning that precedes a SIGKILL is never read in context: the
process dies with no traceback, so the user sees a killed command and has no
reason to connect it to a log line. Refusing costs one `pip install highsbox`
and says so; substituting costs a dead process and an unexplained one.

The measurement is ONE platform -- macOS arm64, mip 2.0.0, cbcbox 2.935 -- so
an installation whose CBC works can ask for it with `NEOSWGA_ALLOW_CBC=1`. That
is an explicit request rather than a substitution nobody asked for, which is the
distinction Task 6 of the valid-design plan draws.
"""

import logging
import os

logger = logging.getLogger(__name__)

#: Set to opt in to CBC where it is known to work. Named rather than inferred,
#: because no in-process probe can tell a working CBC from a fatal one.
_ALLOW_CBC_ENV = "NEOSWGA_ALLOW_CBC"

#: The runtime package mip's HiGHS backend loads. Importable means usable;
#: mip raises a plain ImportError naming this module when it is absent, which
#: is why availability can be decided by import rather than by construction.
_HIGHS_RUNTIME = "highsbox"


def highs_is_available() -> bool:
    """Whether mip's HiGHS backend can load its runtime here."""
    import importlib.util

    return importlib.util.find_spec(_HIGHS_RUNTIME) is not None


def cbc_was_asked_for() -> bool:
    """Whether the operator explicitly accepted CBC."""
    return os.environ.get(_ALLOW_CBC_ENV, "").strip().lower() in {"1", "true", "yes"}


def select_solver_name():
    """The solver name to pass to `mip.Model`, preferring HiGHS.

    CBC is NOT probed, and that is deliberate. A probe would have to construct a
    CBC model, which is the operation that kills the process, so no in-process
    test can tell a working CBC from a fatal one.

    Raises `UnsupportedModelError` when HiGHS is unavailable and CBC has not been
    asked for. Substituting a backend this project has measured to terminate the
    interpreter is worse than refusing: the SIGKILL carries no traceback, so the
    fallback's own warning arrives attached to nothing the user can act on.
    """
    import mip

    if highs_is_available():
        return mip.HIGHS

    if not cbc_was_asked_for():
        from .exceptions import UnsupportedModelError

        raise UnsupportedModelError(
            "mip MILP backend",
            "CBC",
            f"HiGHS. Install its runtime with 'pip install {_HIGHS_RUNTIME}' or "
            f"'pip install neoswga[improved]'. CBC has been measured to "
            f"terminate the interpreter (SIGKILL, no traceback, nothing "
            f"catchable) on Python 3.13 with mip 2.0.0 on macOS arm64, which is "
            f"why it is not substituted silently. If CBC works on your platform, "
            f"set {_ALLOW_CBC_ENV}=1 to use it",
        )

    logger.warning(
        "Using the CBC solver because %s is set. CBC has been measured to "
        "terminate the interpreter (SIGKILL, no traceback, nothing catchable) "
        "on Python 3.13/macOS arm64 with mip 2.0.0.",
        _ALLOW_CBC_ENV,
    )
    return mip.CBC
