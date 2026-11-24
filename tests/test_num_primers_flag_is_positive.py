"""`--num-primers` must be a positive integer.

`-n 0` was read as absent (the handler tests the flag for truthiness), so the
params.json size was used silently, and `-n -3` was written into
`parameter.num_primers`. `_resolve_panel_size` refuses a size below 1 from the
file, so the flag was checked less strictly than the file
(2026-11-23 sweep, tier 2).
"""

import pytest

from neoswga.cli_unified import create_parser


@pytest.mark.parametrize("command", ["optimize", "multi-genome"])
@pytest.mark.parametrize("value", ["0", "-3"])
def test_a_non_positive_size_is_refused(command, value, capsys):
    argv = [command, "-n", value] if command == "optimize" else [command, "--num-primers", value]
    if command == "optimize":
        argv += ["-j", "params.json"]
    else:
        argv += ["--genomes", "a.fasta", "--output", "out"]
    with pytest.raises(SystemExit) as info:
        create_parser().parse_args(argv)
    assert info.value.code == 2
    assert "positive" in capsys.readouterr().err


def test_a_positive_size_parses():
    assert create_parser().parse_args(["optimize", "-j", "p.json", "-n", "4"]).num_primers == 4
