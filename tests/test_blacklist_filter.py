"""Tests for blacklist genome filtering."""

import os
import tempfile
import pytest
import pandas as pd
from unittest.mock import patch


class TestFilterBlacklistPenalty:
    """Test _filter_blacklist_penalty function."""

    def test_no_blacklist_files(self, tmp_path):
        """Primers pass when no blacklist k-mer files exist."""
        from neoswga.core.pipeline import _filter_blacklist_penalty

        primers = ['ATCGATCG', 'GCTAGCTA']
        mask, freqs = _filter_blacklist_penalty(
            primers, [str(tmp_path / 'nonexistent')], [1000], max_bl_freq=0.0
        )
        assert all(mask)
        assert all(f == 0.0 for f in freqs)

    def test_zero_tolerance_rejects_hits(self, tmp_path):
        """With max_bl_freq=0, any hit in blacklist rejects primer."""
        from neoswga.core.pipeline import _filter_blacklist_penalty

        # Create a k-mer file with a hit
        kmer_file = tmp_path / 'bl_test_8mer_all.txt'
        kmer_file.write_text('ATCGATCG 5\nGCTAGCTA 0\n')

        primers = ['ATCGATCG', 'GCTAGCTA']
        mask, freqs = _filter_blacklist_penalty(
            primers, [str(tmp_path / 'bl_test')], [10000], max_bl_freq=0.0
        )
        assert mask[0] is False  # ATCGATCG has hits
        assert mask[1] is True   # GCTAGCTA has 0 hits

    def test_frequency_threshold(self, tmp_path):
        """Primers below max_bl_freq pass."""
        from neoswga.core.pipeline import _filter_blacklist_penalty

        kmer_file = tmp_path / 'bl_test_8mer_all.txt'
        kmer_file.write_text('ATCGATCG 1\nGCTAGCTA 100\n')

        primers = ['ATCGATCG', 'GCTAGCTA']
        mask, freqs = _filter_blacklist_penalty(
            primers, [str(tmp_path / 'bl_test')], [10000], max_bl_freq=0.001
        )
        # ATCGATCG: freq = 1/10000 = 0.0001 <= 0.001 -> pass
        assert mask[0] is True
        # GCTAGCTA: freq = 100/10000 = 0.01 > 0.001 -> reject
        assert mask[1] is False

    def test_bl_freq_values(self, tmp_path):
        """Blacklist frequencies are calculated correctly."""
        from neoswga.core.pipeline import _filter_blacklist_penalty

        kmer_file = tmp_path / 'bl_test_6mer_all.txt'
        kmer_file.write_text('ATCGAT 50\n')

        primers = ['ATCGAT']
        _, freqs = _filter_blacklist_penalty(
            primers, [str(tmp_path / 'bl_test')], [5000], max_bl_freq=1.0
        )
        assert abs(freqs[0] - 50 / 5000) < 1e-10

    def test_multiple_blacklist_genomes(self, tmp_path):
        """Counts are summed across multiple blacklist genomes."""
        from neoswga.core.pipeline import _filter_blacklist_penalty

        # Two blacklist genomes
        bl1 = tmp_path / 'bl1_8mer_all.txt'
        bl1.write_text('ATCGATCG 10\n')
        bl2 = tmp_path / 'bl2_8mer_all.txt'
        bl2.write_text('ATCGATCG 20\n')

        primers = ['ATCGATCG']
        _, freqs = _filter_blacklist_penalty(
            primers,
            [str(tmp_path / 'bl1'), str(tmp_path / 'bl2')],
            [5000, 5000],
            max_bl_freq=1.0
        )
        # freq = (10 + 20) / (5000 + 5000) = 0.003
        assert abs(freqs[0] - 30 / 10000) < 1e-10


class TestBlacklistParameterFields:
    """Test that blacklist fields exist in parameter system."""

    def test_pipeline_parameters_has_bl_fields(self):
        """PipelineParameters has all blacklist fields."""
        from neoswga.core.parameter import PipelineParameters

        p = PipelineParameters()
        assert p.bl_genomes == []
        assert p.bl_prefixes == []
        assert p.bl_seq_lengths == []
        assert p.bl_penalty == 5.0
        assert p.max_bl_freq == 0.0

    def test_get_current_config_includes_bl(self):
        """get_current_config returns bl_ fields."""
        from neoswga.core.parameter import get_current_config

        config = get_current_config()
        assert hasattr(config, 'bl_genomes')
        assert hasattr(config, 'bl_prefixes')
        assert hasattr(config, 'bl_penalty')
        assert hasattr(config, 'max_bl_freq')

    def test_set_from_config_roundtrip(self):
        """set_from_config preserves bl_ fields."""
        from neoswga.core.parameter import (
            PipelineParameters, set_from_config, get_current_config
        )

        config = PipelineParameters(
            bl_genomes=['/tmp/test.fna'],
            bl_prefixes=['/tmp/bl_test'],
            bl_penalty=3.0,
            max_bl_freq=0.001,
        )
        set_from_config(config)
        restored = get_current_config()
        assert restored.bl_genomes == ['/tmp/test.fna']
        assert restored.bl_prefixes == ['/tmp/bl_test']
        assert restored.bl_penalty == 3.0
        assert restored.max_bl_freq == 0.001

        # Reset to defaults
        from neoswga.core.parameter import reset_to_defaults
        reset_to_defaults()

    def test_backward_compat_no_bl_fields(self):
        """params.json without bl_ fields works (backward compatibility)."""
        from neoswga.core.parameter import PipelineParameters

        p = PipelineParameters()
        assert p.bl_genomes == []
        assert p.bl_penalty == 5.0
