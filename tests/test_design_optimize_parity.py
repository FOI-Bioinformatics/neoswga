"""`design` must resolve step 4 the way `optimize` does from the same params.json.

`run_design` back-fills optimize-only attributes onto the `design` namespace.
It back-filled `mechanistic_weight=0.3`, which `run_step4` reads as an explicit
request for the mechanistic model, so `design` scored primers with a term that
`optimize` on the same file leaves off (2026-11-23 sweep, tier 1 item 2).
"""

from neoswga.cli.commands import _apply_design_defaults
from neoswga.cli_unified import create_parser


def _design_args():
    args = create_parser().parse_args(["design", "-j", "params.json"])
    _apply_design_defaults(args)
    return args


def test_design_defaults_do_not_enable_mechanistic_model():
    args = _design_args()
    assert args.mechanistic_weight is None
    assert args.use_mechanistic_model is False


def test_design_defaults_use_the_attribute_step4_reads_for_simulation():
    args = _design_args()
    assert "validate_with_simulation" not in vars(args)
    assert args.validate_simulation is False
