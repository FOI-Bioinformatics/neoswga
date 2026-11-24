"""`--bl-penalty` is refused.

The `bl_penalty` key was retired on 2026-09-14 because no scoring code read
it. The flag survived on `count-kmers` and `filter`, wrote the global, and
changed nothing (2026-11-23 sweep, tier 3).
"""

import pytest

from neoswga.cli_unified import create_parser


@pytest.mark.parametrize("command", ["count-kmers", "filter"])
def test_the_flag_is_refused(command):
    with pytest.raises(SystemExit) as info:
        create_parser().parse_args([command, "-j", "p.json", "--bl-penalty", "3"])
    assert info.value.code == 2
