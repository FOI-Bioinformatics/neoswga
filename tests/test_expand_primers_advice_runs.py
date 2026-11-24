"""The commands `expand-primers` suggests next work on what it wrote.

It suggested `interpret`, `report` and `export -d OUT`, which read
`step4_improved_df.csv`; it writes `expanded_primers.csv`, whose
`primer,type` rows the `--primers-file` readers skip line by line, so the
expanded set had no working next step (2026-11-23 sweep, tier 3).
"""

import argparse
import logging
import shlex
from types import SimpleNamespace

from neoswga.cli._common import load_primers_from_file
from neoswga.cli.iterate import _report_expansion_summary, _save_expanded_primers
from neoswga.cli_unified import create_parser

RESULT = SimpleNamespace(
    n_fixed=2,
    n_new=1,
    n_total=3,
    coverage_before=0.5,
    coverage_after=0.6,
    predicted_improvement=1.2,
    gaps_remaining=4,
    new_primers=["TTGACCATGA"],
    combined_set=["ACGTACGTAC", "GGCCAATTGC", "TTGACCATGA"],
)


def test_every_suggested_command_parses_and_reads_a_primer_list(tmp_path, caplog):
    args = argparse.Namespace(output=str(tmp_path), json_file="params.json", bam=None)
    with caplog.at_level(logging.INFO, logger="neoswga.cli.iterate"):
        _report_expansion_summary(RESULT, args, quiet=False)
    suggested = [
        line.split("#")[0].strip()
        for line in (record.getMessage() for record in caplog.records)
        if line.strip().startswith("neoswga ")
    ]
    assert suggested
    for command in suggested:
        argv = shlex.split(command)[1:]
        assert argv[0] not in ("interpret", "report", "export"), command
        create_parser().parse_args(argv)


def test_the_suggested_primer_file_loads_the_expanded_set(tmp_path):
    path = _save_expanded_primers(RESULT, str(tmp_path))
    assert load_primers_from_file(path) == RESULT.combined_set
