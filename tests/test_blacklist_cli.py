"""Tests for CLI blacklist argument parsing."""

import pytest
import sys


class TestBlacklistCLIArguments:
    """Test that CLI parsers accept blacklist arguments."""

    def test_count_kmers_blacklist_args(self):
        """count-kmers parser accepts --blacklist."""
        from neoswga.cli_unified import create_parser

        parser = create_parser()
        args = parser.parse_args(
            [
                "count-kmers",
                "-j",
                "params.json",
                "--blacklist",
                "genome1.fna",
                "genome2.fna",
                "--bl-penalty",
                "3.0",
                "--max-bl-freq",
                "0.001",
            ]
        )
        assert args.blacklist == ["genome1.fna", "genome2.fna"]
        assert args.bl_penalty == 3.0
        assert args.max_bl_freq == 0.001

    def test_filter_blacklist_args(self):
        """filter parser accepts --blacklist."""
        from neoswga.cli_unified import create_parser

        parser = create_parser()
        args = parser.parse_args(
            [
                "filter",
                "-j",
                "params.json",
                "--blacklist",
                "contaminant.fna",
            ]
        )
        assert args.blacklist == ["contaminant.fna"]

    def test_init_blacklist_args(self):
        """init parser accepts --blacklist."""
        from neoswga.cli_unified import create_parser

        parser = create_parser()
        args = parser.parse_args(
            [
                "init",
                "--genome",
                "target.fna",
                "--blacklist",
                "contaminant1.fna",
                "contaminant2.fna",
            ]
        )
        assert args.blacklist == ["contaminant1.fna", "contaminant2.fna"]

    def test_genome_add_args(self):
        """genome-add parser accepts required arguments."""
        from neoswga.cli_unified import create_parser

        parser = create_parser()
        args = parser.parse_args(
            [
                "genome-add",
                "human-grch38",
                "human.fna",
                "--role",
                "bg",
                "--k-ranges",
                "6-12",
            ]
        )
        assert args.name == "human-grch38"
        assert args.fasta == "human.fna"
        assert args.role == "bg"
        assert args.k_ranges == "6-12"

    def test_genome_list_args(self):
        """genome-list parser works."""
        from neoswga.cli_unified import create_parser

        parser = create_parser()
        args = parser.parse_args(["genome-list", "--verbose"])
        assert args.command == "genome-list"

    def test_genome_remove_args(self):
        """genome-remove parser accepts name argument."""
        from neoswga.cli_unified import create_parser

        parser = create_parser()
        args = parser.parse_args(["genome-remove", "human-grch38"])
        assert args.name == "human-grch38"

    def test_count_kmers_no_blacklist(self):
        """count-kmers works without --blacklist (backward compat)."""
        from neoswga.cli_unified import create_parser

        parser = create_parser()
        args = parser.parse_args(["count-kmers", "-j", "params.json"])
        assert args.blacklist is None
        assert args.bl_penalty is None
        assert args.max_bl_freq is None
