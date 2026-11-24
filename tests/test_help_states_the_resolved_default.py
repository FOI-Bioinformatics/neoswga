"""Help texts state the defaults the code resolves.

`filter --max-gini` said 0.6 where every other statement says 0.7, and
`--min-tm`/`--max-tm` said 15/45 where the window comes from the polymerase
(2026-11-23 sweep, tier 3).
"""

import re

import pytest

from neoswga.cli_unified import create_parser
from neoswga.core.parameter_config import PipelineParameters


def _filter_help():
    parser = create_parser()
    sub = next(a for a in parser._actions if a.__class__.__name__ == "_SubParsersAction")
    return {a.dest: a.help or "" for a in sub.choices["filter"]._actions}


@pytest.mark.parametrize("key", ["max_gini", "min_fg_freq", "max_bg_freq", "max_primer"])
def test_a_stated_default_is_the_resolved_one(key):
    stated = re.search(r"default:\s*([0-9.eE+-]+)", _filter_help()[key])
    assert stated, _filter_help()[key]
    assert float(stated.group(1)) == float(PipelineParameters.__dataclass_fields__[key].default)


@pytest.mark.parametrize("key", ["min_tm", "max_tm"])
def test_the_tm_window_help_names_the_polymerase(key):
    text = _filter_help()[key]
    assert "polymerase" in text and not re.search(r"default:\s*[0-9]", text), text
