"""Which MILP backend python-mip is handed.

python-mip's default is CBC. On Python 3.13 (measured on macOS arm64 with mip
2.0.0 and cbcbox 2.935) constructing a CBC model terminates the interpreter
with SIGKILL: no exception, no traceback, no stderr. The same versions on
Python 3.11 build and solve the same models in a third of a second.

mip 2.0 ships a second backend. HiGHS solves those models on 3.13 and on 3.11,
so it is preferred wherever its runtime is importable, and CBC remains the
fallback so nothing changes for an installation that does not have it.
"""

import logging

logger = logging.getLogger(__name__)

#: The runtime package mip's HiGHS backend loads. Importable means usable;
#: mip raises a plain ImportError naming this module when it is absent, which
#: is why availability can be decided by import rather than by construction.
_HIGHS_RUNTIME = "highsbox"


def highs_is_available() -> bool:
    """Whether mip's HiGHS backend can load its runtime here."""
    import importlib.util

    return importlib.util.find_spec(_HIGHS_RUNTIME) is not None


def select_solver_name():
    """The solver name to pass to `mip.Model`, preferring HiGHS.

    The CBC fallback is NOT probed, and that is deliberate. A probe would have
    to construct a CBC model, which is the operation that kills the process,
    so no in-process test can tell a working CBC from a fatal one. The only
    honest options are to prefer the backend that is known to work and to say
    so when falling back on an interpreter where CBC has been measured to die.
    """
    import mip

    if highs_is_available():
        return mip.HIGHS

    # Warned unconditionally. This used to be guarded by a 3.13 version check,
    # which ruff correctly called dead once `requires-python` became ">=3.13":
    # every supported interpreter is one where CBC has been measured to die.
    logger.warning(
        "Falling back to the CBC solver. CBC has been measured to terminate "
        "the interpreter (SIGKILL, no traceback, nothing catchable) on Python "
        "3.13/macOS arm64 with mip 2.0.0. Install the HiGHS runtime to avoid "
        "it: pip install %s, or pip install 'neoswga[improved]'.",
        _HIGHS_RUNTIME,
    )
    return mip.CBC
