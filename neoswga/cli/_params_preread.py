"""Reading params.json before ``get_params()`` has run.

`neoswga.core.parameter` is a module of globals, populated lazily: `get_params()`
runs deep inside the step it belongs to, not when the command starts. Every
global therefore carries a module-level default until then -- `num_primers` is
6, `polymerase` is "phi29" -- and a command handler that reads one *before* that
point silently gets the default rather than what the user wrote.

That is not a hypothetical. It has now happened twice in `run_step4`, twenty
lines apart:

- `num_primers`: a params.json asking for 8 primers got 6, and the completion
  check then compared the 6 produced against the 8 it read later and advised
  relaxing the filters for a shortfall that did not exist.
- `polymerase`: a params.json asking for bst was optimized with phi29's presets,
  Tm window and coverage reach, changing which primers were selected.

Both are the same defect, and the fix is the same: read params.json here, and
prefer the module global only once `get_params()` has populated it -- at which
point it is authoritative, because it carries any adaptation applied on top of
the file (the GC-adaptive path can swap the enzyme, for instance).

These helpers deliberately never raise. `get_params()` and `param_validator`
report a malformed or invalid file far better than a pre-read can, and
pre-empting them with a partial failure names the wrong problem.

New entries belong here rather than in the command module, so the pattern stays
visible; `tests/test_cli_optimize_target_size.py` and
`tests/test_cli_optimize_polymerase.py` pin the two that exist.
"""

import logging
import os

logger = logging.getLogger(__name__)


def _preread(parameter, field):
    """Return ``(already_loaded, json_data)`` for a params.json lookup.

    `already_loaded` says whether `get_params()` has run. It is the load-bearing
    half: only then does a module global carry the user's value rather than its
    own default, and only then may a caller prefer the global over the file.
    """
    loaded = dict(getattr(parameter, "_json_data", None) or {})
    if loaded:
        return True, loaded

    json_file = getattr(parameter, "json_file", None)
    if json_file and os.path.isfile(json_file):
        from neoswga.core.parameter import read_args_from_json

        try:
            return False, read_args_from_json(json_file)
        except (OSError, ValueError) as exc:
            logger.debug(f"Could not pre-read {json_file} for {field}: {exc}")

    return False, {}


def target_size_from_params(parameter, default=6):
    """Resolve the requested set size from params.json, at optimize time.

    Reading the file here rather than waiting for `_json_data` keeps params.json
    authoritative for the set size, which is what a user editing it expects.
    `--num-primers` still wins; this is the no-flag path.
    """
    _, json_data = _preread(parameter, "num_primers")

    size = json_data.get("num_primers", json_data.get("target_set_size", default))
    if not isinstance(size, int) or isinstance(size, bool) or size < 1:
        logger.warning(
            f"Ignoring num_primers={size!r} from params.json: expected a "
            f"positive integer. Using {default}."
        )
        return default
    return size


def apply_polymerase_choice(parameter, args):
    """Honour `--polymerase` / `--preset` on a command that has already inited.

    `filter` runs `pipeline._initialize()` first, so `get_params` has already
    derived the Tm window, primer-length range and buffer from whatever
    params.json named. Setting the enzyme afterwards without re-deriving them
    filtered a bst design through phi29's 20-50 C window and 6-12 bp lengths.

    Anything the flags merged is the user's and must survive the retune;
    `get_params` could not see those, having run with its own args object.
    """
    from neoswga.cli._common import load_preset_conditions

    requested = getattr(args, "polymerase", None)
    if not requested and getattr(args, "preset", None):
        requested = load_preset_conditions(args.preset).get("polymerase")
    if not requested:
        return

    parameter.mark_user_set(
        parameter,
        [
            name
            for name in parameter.POLYMERASE_DERIVED_FIELDS
            if getattr(args, name, None) is not None
        ],
    )
    parameter.retune_for_polymerase(parameter, requested)


def polymerase_from_params(parameter, default="phi29", args=None):
    """Resolve the polymerase for a command, at command time.

    Precedence: an explicit `--polymerase` flag, then -- once `get_params()` has
    run -- `parameter.polymerase`, because it carries any GC-adaptive adaptation
    applied on top of params.json, and before that params.json itself. An
    unknown name falls back to `default`.

    The flag has to be handled here rather than by the usual merge because
    `run_step4` is decorated `@params_command(merge=None)`: the subparser offered
    all six enzymes and nothing ever read the value.
    """
    from neoswga.core.registry.polymerases import resolve_polymerase_name

    already_loaded, json_data = _preread(parameter, "polymerase")

    flag = getattr(args, "polymerase", None) if args is not None else None
    if flag is not None:
        resolved = resolve_polymerase_name(flag) if isinstance(flag, str) else None
        if resolved is not None:
            return resolved
        logger.warning(
            f"Ignoring --polymerase {flag!r}: not a known enzyme. "
            f"Using the configured polymerase instead."
        )

    if already_loaded:
        adapted = resolve_polymerase_name(getattr(parameter, "polymerase", None))
        if adapted:
            return adapted

    raw = json_data.get("polymerase")
    if raw is None:
        return default

    resolved = resolve_polymerase_name(raw) if isinstance(raw, str) else None
    if resolved is None:
        if raw:
            logger.warning(
                f"Ignoring polymerase={raw!r} from params.json: not a known "
                f"enzyme. Using {default}."
            )
        return default
    return resolved
