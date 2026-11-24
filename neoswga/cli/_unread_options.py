"""Options a command declares and does not read refuse when given.

Each entry was accepted and ignored. Several read as choosing something the run
then did not use: `optimize --polymerase bst` optimized the params.json enzyme,
and `count-kmers -x FASTA` counted the params.json genomes. An option a command cannot honour now
stops the command and names what does what the option says (2026-11-24).

`tests/test_every_cli_option_has_an_effect.py` holds the list of options no
code path reads; an option there either says so in its help or appears here.
"""

import argparse


def _params_key(key):
    return f"set `{key}` in params.json"


UNREAD_OPTIONS = {
    "count-kmers": {
        "data_dir": _params_key("data_dir"),
        "polymerase": _params_key("polymerase"),
        "fasta_fore": _params_key("fg_genomes"),
        "fasta_back": _params_key("bg_genomes"),
        "kmer_fore": _params_key("fg_prefixes"),
        "kmer_back": _params_key("bg_prefixes"),
    },
    "filter": {
        "data_dir": _params_key("data_dir"),
        "min_gini_sites": _params_key("min_gini_sites"),
    },
    "prepare-candidates": {
        "data_dir": _params_key("data_dir"),
        "polymerase": _params_key("polymerase"),
    },
    "optimize": {
        "data_dir": _params_key("data_dir"),
        "polymerase": _params_key("polymerase"),
    },
    "design": {"data_dir": _params_key("data_dir")},
    "evaluate-set": {
        "data_dir": _params_key("data_dir"),
        "enable_qa": "this command runs no QA step",
    },
    "analyze-dimers": {"visualize": "no network plot is drawn; the report goes to --output"},
    "analyze-set": {"simulate": "run `neoswga simulate` on the set"},
    "analyze-stability": {
        "temp": "the analysis evaluates terminal stability without a reaction temperature"
    },
}


class _Unread(argparse.Action):
    """Refuses the option, naming the route that does what it says."""

    def __init__(self, option_strings, dest, route, **kwargs):
        super().__init__(option_strings, dest, **kwargs)
        self.route = route

    def __call__(self, parser, namespace, values, option_string=None):
        parser.error(f"{option_string} has no effect on this command: {self.route}")


def _replace(parser, old, new):
    parser._actions[parser._actions.index(old)] = new
    for option in old.option_strings:
        parser._option_string_actions[option] = new
    for group in parser._action_groups + parser._mutually_exclusive_groups:
        if old in group._group_actions:
            group._group_actions[group._group_actions.index(old)] = new


def refuse_unread_options(subparsers):
    """Swap each unread option's action for one that refuses it."""
    for command, routes in UNREAD_OPTIONS.items():
        parser = subparsers.choices[command]
        for action in list(parser._actions):
            route = routes.get(action.dest)
            if route is None:
                continue
            help_text = f"{action.help or ''} Not read by this command: {route}.".strip()
            replacement = _Unread(
                action.option_strings,
                action.dest,
                route,
                nargs=action.nargs,
                default=action.default,
                metavar=action.metavar,
                help=help_text,
            )
            _replace(parser, action, replacement)
