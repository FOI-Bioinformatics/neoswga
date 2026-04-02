"""CLI smoke tests for neoswga.

Verifies that:
- The main parser creates successfully
- --version works
- Every registered subcommand accepts --help without error
- Basic argument parsing works for key commands
"""

import subprocess
import sys

import pytest

from neoswga.cli_unified import create_parser


class TestParserCreation:
    """The parser should create without errors."""

    def test_create_parser(self):
        parser = create_parser()
        assert parser is not None
        assert parser.prog == "neoswga"


class TestVersionFlag:

    def test_version_output(self):
        result = subprocess.run(
            [sys.executable, "-m", "neoswga.cli_unified", "--version"],
            capture_output=True, text=True
        )
        assert result.returncode == 0
        assert "neoswga" in result.stdout
        assert "3." in result.stdout  # version 3.x


class TestSubcommandHelp:
    """Every subcommand should accept --help and exit 0."""

    SUBCOMMANDS = [
        "count-kmers",
        "filter",
        "score",
        "optimize",
        "design",
        "build-filter",
        "validate",
        "show-presets",
        "init",
        "validate-params",
        "validate-model",
        "interpret",
        "report",
        "export",
        "start",
        "suggest",
        "optimize-conditions",
        "analyze-set",
        "analyze-genome",
        "analyze-dimers",
        "analyze-stability",
        "auto-pipeline",
        "multi-genome",
        "simulate",
        "active-learn",
        "expand-primers",
        "predict-efficiency",
        "background-list",
        "background-add",
        "genome-add",
        "genome-list",
        "genome-remove",
    ]

    @pytest.mark.parametrize("cmd", SUBCOMMANDS)
    def test_help_exits_zero(self, cmd):
        """'neoswga <cmd> --help' should exit 0 and produce output."""
        parser = create_parser()
        with pytest.raises(SystemExit) as exc_info:
            parser.parse_args([cmd, "--help"])
        assert exc_info.value.code == 0


class TestArgumentParsing:
    """Basic argument parsing for key commands."""

    def test_count_kmers_requires_json(self):
        parser = create_parser()
        args = parser.parse_args(["count-kmers", "-j", "params.json"])
        assert args.command == "count-kmers"
        assert args.json_file == "params.json"

    def test_filter_requires_json(self):
        parser = create_parser()
        args = parser.parse_args(["filter", "-j", "params.json"])
        assert args.command == "filter"
        assert args.json_file == "params.json"

    def test_optimize_method_choices(self):
        parser = create_parser()
        args = parser.parse_args([
            "optimize", "-j", "params.json",
            "--optimization-method", "hybrid"
        ])
        assert args.optimization_method == "hybrid"

    def test_optimize_accepts_all_methods(self):
        parser = create_parser()
        methods = [
            "hybrid", "greedy", "milp", "network",
            "genetic", "dominating-set", "background-aware", "moea",
        ]
        for method in methods:
            args = parser.parse_args([
                "optimize", "-j", "params.json",
                "--optimization-method", method
            ])
            assert args.optimization_method == method

    def test_report_level_choices(self):
        parser = create_parser()
        args = parser.parse_args(["report", "-d", "results/", "--level", "full"])
        assert args.command == "report"

    def test_validate_quick_flag(self):
        parser = create_parser()
        args = parser.parse_args(["validate", "--quick"])
        assert args.command == "validate"
        assert args.quick is True

    def test_verbose_and_quiet_flags(self):
        parser = create_parser()
        args = parser.parse_args(["count-kmers", "-j", "p.json", "-v"])
        assert args.verbose is True

        args = parser.parse_args(["count-kmers", "-j", "p.json", "-q"])
        assert args.quiet is True

    def test_no_command_returns_none(self):
        parser = create_parser()
        args = parser.parse_args([])
        assert args.command is None
