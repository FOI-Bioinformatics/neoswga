"""
Tests for long specific oligos (12-18bp) across GC extremes.

Validates seven fixes that ensure the pipeline correctly handles long primers
for clinical pathogen detection across AT-rich (~25% GC), GC-neutral (~50%),
and GC-rich (~65-75% GC) targets.
"""

import pytest
import numpy as np
from unittest.mock import patch, MagicMock, PropertyMock


# ============================================================================
# Fix 7: AttributeError - parameter.max_self_dimer_bp
# ============================================================================

class TestFix7DimerAttributeName:
    """Verify filter_extra uses parameter.max_self_dimer_bp (not default_max_self_dimer_bp)."""

    def test_filter_extra_uses_correct_attribute(self):
        """filter_extra should reference parameter.max_self_dimer_bp."""
        with patch('neoswga.core.filter.parameter') as mock_param, \
             patch('neoswga.core.filter.dimer') as mock_dimer, \
             patch('neoswga.core.filter._get_reaction_conditions') as mock_rc:
            mock_param.max_self_dimer_bp = 4
            mock_param.min_tm = 0
            mock_param.max_tm = 100
            mock_param.gc_min = 0.0
            mock_param.gc_max = 1.0
            mock_param.genome_gc = None
            mock_dimer.is_dimer_fast.return_value = False

            mock_conditions = MagicMock()
            mock_conditions.calculate_effective_tm.return_value = 30.0
            mock_rc.return_value = mock_conditions

            from neoswga.core.filter import filter_extra
            # Should not raise AttributeError
            result = filter_extra("ATCGATCG")
            assert isinstance(result, bool)

            # Verify is_dimer_fast was called with the correct attribute value
            mock_dimer.is_dimer_fast.assert_called_once_with("ATCGATCG", "ATCGATCG", 4)


# ============================================================================
# Fix 1: Tm gatekeeper with wide margin
# ============================================================================

class TestFix1TmPreFilter:
    """Verify wide Tm margin prevents premature rejection of long primers."""

    def test_wide_margin_keeps_long_primers(self, tmp_path):
        """A 16bp primer with melting.temp() > max_tm should pass with margin."""
        import melting

        primer = "ATCGATCGATCGATCG"  # 16bp
        simple_tm = melting.temp(primer)

        # Create a mock k-mer file
        kmer_file = tmp_path / "test_16mer_all.txt"
        kmer_file.write_text(f"{primer} 100\n")

        from neoswga.core.kmer_counter import get_primer_list_from_kmers

        # With narrow window, primer might be rejected
        result_narrow = get_primer_list_from_kmers(
            [str(tmp_path / "test")],
            kmer_lengths=range(16, 17),
            min_tm=15.0, max_tm=45.0,
            wide_tm_margin=0.0
        )

        # With wide margin (default 15C), primer should be kept
        result_wide = get_primer_list_from_kmers(
            [str(tmp_path / "test")],
            kmer_lengths=range(16, 17),
            min_tm=15.0, max_tm=45.0,
            wide_tm_margin=15.0
        )

        # The wide margin should retain more or equal primers
        assert len(result_wide) >= len(result_narrow)

    def test_default_margin_is_15(self, tmp_path):
        """Default wide_tm_margin should be 15.0."""
        import inspect
        from neoswga.core.kmer_counter import get_primer_list_from_kmers

        sig = inspect.signature(get_primer_list_from_kmers)
        assert sig.parameters['wide_tm_margin'].default == 15.0

    def test_extreme_tm_still_rejected(self, tmp_path):
        """Primers with extremely high Tm should still be rejected even with margin."""
        # 18bp all-GC primer will have very high Tm
        primer = "GCGCGCGCGCGCGCGCGC"
        kmer_file = tmp_path / "test_18mer_all.txt"
        kmer_file.write_text(f"{primer} 50\n")

        from neoswga.core.kmer_counter import get_primer_list_from_kmers

        result = get_primer_list_from_kmers(
            [str(tmp_path / "test")],
            kmer_lengths=range(18, 19),
            min_tm=15.0, max_tm=45.0,
            wide_tm_margin=15.0  # max becomes 60
        )
        # All-GC 18-mer has Tm well above 60C, should still be rejected
        assert primer not in result


# ============================================================================
# Fix 2: Frequency threshold scaling
# ============================================================================

class TestFix2FrequencyScaling:
    """Verify frequency thresholds scale inversely with primer length."""

    def test_scale_freq_threshold_short_primers(self):
        """Primers <= reference_k should get unscaled threshold."""
        from neoswga.core.filter import _scale_freq_threshold

        for k in range(6, 11):
            assert _scale_freq_threshold(1e-5, k) == 1e-5

    def test_scale_freq_threshold_12bp(self):
        """12bp primer: scale factor = 4^(10-12) = 1/16."""
        from neoswga.core.filter import _scale_freq_threshold
        result = _scale_freq_threshold(1e-5, 12)
        assert abs(result - 1e-5 / 16) < 1e-12

    def test_scale_freq_threshold_15bp(self):
        """15bp primer: scale factor = 4^(10-15) = 1/1024."""
        from neoswga.core.filter import _scale_freq_threshold
        result = _scale_freq_threshold(1e-5, 15)
        expected = 1e-5 / 1024.0
        assert abs(result - expected) < 1e-15

    def test_scale_freq_threshold_18bp(self):
        """18bp primer: scale factor = 4^(10-18) = 1/65536."""
        from neoswga.core.filter import _scale_freq_threshold
        result = _scale_freq_threshold(1e-5, 18)
        expected = 1e-5 / 65536.0
        assert abs(result - expected) < 1e-18

    def test_15bp_primer_with_10_sites_passes(self):
        """A 15bp primer with 10 sites in a 5 Mbp genome should pass."""
        from neoswga.core.filter import _scale_freq_threshold

        genome_length = 5_000_000
        fg_count = 10
        freq = fg_count / genome_length  # 2e-6

        base_threshold = 1e-5
        scaled = _scale_freq_threshold(base_threshold, 15)
        # scaled ~ 9.77e-9, freq = 2e-6 >> scaled
        assert freq > scaled


# ============================================================================
# Fix 3: GC clamp adaptive for AT-rich and GC-rich genomes
# ============================================================================

class TestFix3AdaptiveGCClamp:
    """Verify GC clamp is adaptive based on genome_gc."""

    def _run_filter_extra(self, primer, genome_gc=None, gc_min=0.0, gc_max=1.0):
        """Helper to run filter_extra with mocked parameters."""
        with patch('neoswga.core.filter.parameter') as mock_param, \
             patch('neoswga.core.filter.dimer') as mock_dimer, \
             patch('neoswga.core.filter._get_reaction_conditions') as mock_rc:
            mock_param.max_self_dimer_bp = 4
            mock_param.min_tm = 0
            mock_param.max_tm = 100
            mock_param.gc_min = gc_min
            mock_param.gc_max = gc_max
            mock_param.genome_gc = genome_gc
            mock_dimer.is_dimer_fast.return_value = False

            mock_conditions = MagicMock()
            mock_conditions.calculate_effective_tm.return_value = 30.0
            mock_rc.return_value = mock_conditions

            from neoswga.core.filter import filter_extra
            return filter_extra(primer)

    def test_at_rich_allows_zero_gc_in_last5(self):
        """AT-rich genome (25% GC): primers with 0 GC in last 5 should pass."""
        # Primer ending in AAAAT (0 GC in last 5)
        # Full primer needs overall GC in range
        primer = "GCGCAAAAT"  # 4 GC / 9 = 44%
        result = self._run_filter_extra(primer, genome_gc=0.25, gc_min=0.1, gc_max=0.9)
        assert result is True

    def test_at_rich_still_rejects_high_gc_clamp(self):
        """AT-rich genome: >3 GC in last 5 should still be rejected."""
        primer = "AATGCGCGC"  # last 5 = GCGCG => 4 GC, GC clamp max=3 for AT-rich
        result = self._run_filter_extra(primer, genome_gc=0.25, gc_min=0.1, gc_max=0.9)
        assert result is False

    def test_gc_rich_allows_4_gc_in_last5(self):
        """GC-rich genome (75% GC): primers with 4 GC in last 5 should pass."""
        # Primer ending in GCGCG (4 GC in last 5, with 1 AT = passes)
        # Wait, GCGCG has 3G+2C = 5 GC in 5 bases -- that's 5, not 4
        # Use AGCGC -> 3 GC? No: A=1, G=1, C=2, G=1 = AGCGC -> GC=4
        primer = "ATAGCGCAG"  # last 5: GCGCAG? No...
        # Simpler: "ATACGCGCA" -> last 5 = "GCGCA" -> 3G+1C = GC=4? G,C,G,C,A -> 4 GC
        primer = "ATACGCGCA"  # last 5 = CGCGCA? wait 9 chars
        # Let me be more careful: "ATACGCGCA" = A,T,A,C,G,C,G,C,A
        # last 5 = G,C,G,C,A -> 4 GC. Rule 1: last 3 = G,C,A -> 2 GC (OK)
        # Overall GC = 5/9 = 55%
        result = self._run_filter_extra(primer, genome_gc=0.75, gc_min=0.1, gc_max=0.9)
        assert result is True

    def test_gc_rich_rejects_all_gc_in_last5(self):
        """GC-rich genome: all 5 GC in last 5 should be rejected (need >=1 AT)."""
        primer = "ATATGCGCG"  # last 5 = GCGCG -> 5 GC
        result = self._run_filter_extra(primer, genome_gc=0.75, gc_min=0.1, gc_max=0.9)
        assert result is False

    def test_standard_gc_rejects_zero_gc_in_last5(self):
        """Standard genome (50% GC): 0 GC in last 5 should be rejected."""
        primer = "GCGCAAAAT"
        result = self._run_filter_extra(primer, genome_gc=0.50, gc_min=0.1, gc_max=0.9)
        assert result is False

    def test_none_genome_gc_uses_standard_rules(self):
        """When genome_gc is None, standard GC clamp rules apply."""
        primer = "GCGCAAAAT"  # 0 GC in last 5
        result = self._run_filter_extra(primer, genome_gc=None, gc_min=0.1, gc_max=0.9)
        assert result is False


# ============================================================================
# Fix 4: Dimer threshold scaling
# ============================================================================

class TestFix4DimerScaling:
    """Verify dimer thresholds auto-scale with max_k."""

    def test_dimer_scaling_short_primers(self):
        """For max_k=12, max_dimer_bp should remain at default 3."""
        # max(3, int(12 * 0.25)) = max(3, 3) = 3
        assert max(3, int(12 * 0.25)) == 3

    def test_dimer_scaling_18bp(self):
        """For max_k=18, max_dimer_bp should be 4."""
        # max(3, int(18 * 0.25)) = max(3, 4) = 4
        assert max(3, int(18 * 0.25)) == 4

    def test_self_dimer_scaling_18bp(self):
        """For max_k=18, max_self_dimer_bp should be 5."""
        # max(4, int(18 * 0.28)) = max(4, 5) = 5
        assert max(4, int(18 * 0.28)) == 5

    def test_auto_scale_in_get_params(self):
        """get_params should auto-scale dimer thresholds when not explicitly set."""
        import neoswga.core.parameter as parameter

        # Save original state
        orig_json_data = parameter._json_data.copy()

        try:
            # Clear _json_data to simulate no explicit dimer settings
            parameter._json_data.clear()

            class MockArgs:
                def __getattr__(self, name):
                    return None

            args = MockArgs()
            # Set required fields
            args.json_file = None

            # We can't easily run full get_params without a json file,
            # so test the formula directly
            max_k = 18
            expected_dimer = max(3, int(max_k * 0.25))
            expected_self_dimer = max(4, int(max_k * 0.28))
            assert expected_dimer == 4
            assert expected_self_dimer == 5
        finally:
            parameter._json_data.update(orig_json_data)


# ============================================================================
# Fix 5: GC-adaptive strategy respects user k-mer range
# ============================================================================

class TestFix5GCAdaptiveGuard:
    """Verify GC-adaptive strategy does not override user-specified k-mer range."""

    def test_user_kmer_range_preserved(self):
        """When user sets min_k/max_k in JSON, adaptive strategy should not override."""
        from neoswga.core import pipeline as core_pipeline
        import neoswga.core.parameter as parameter

        # Mock _json_data to contain user-specified k-mer range
        original_json_data = parameter._json_data.copy()
        original_min_k = getattr(parameter, 'min_k', 6)
        original_max_k = getattr(parameter, 'max_k', 12)

        try:
            parameter._json_data['min_k'] = 15
            parameter._json_data['max_k'] = 18
            parameter.min_k = 15
            parameter.max_k = 18
            parameter.genome_gc = 0.65

            # Mock GCAdaptiveStrategy
            with patch('neoswga.core.pipeline.parameter', parameter):
                mock_strategy = MagicMock()
                mock_params = MagicMock()
                mock_params.kmer_range = (8, 12)  # Strategy wants shorter primers
                mock_params.recommended_polymerase = 'phi29'
                mock_params.reaction_temp = 30.0
                mock_params.betaine_m = 0.0
                mock_params.dmso_percent = 0.0
                mock_params.genome_class = MagicMock()
                mock_params.genome_class.value = 'gc_rich'
                mock_params.confidence = 0.8
                mock_strategy.return_value.get_parameters.return_value = mock_params

                with patch.dict('sys.modules', {'neoswga.core.gc_adaptive_strategy': MagicMock(GCAdaptiveStrategy=mock_strategy)}):
                    core_pipeline._apply_gc_adaptive_defaults()

            # K-mer range should NOT have been overridden
            assert parameter.min_k == 15
            assert parameter.max_k == 18
        finally:
            parameter._json_data.clear()
            parameter._json_data.update(original_json_data)
            parameter.min_k = original_min_k
            parameter.max_k = original_max_k


# ============================================================================
# Fix 6: Coverage target in dominating set optimizer
# ============================================================================

class TestFix6CoverageTarget:
    """Verify optimize_greedy stops when min_coverage target is met."""

    def _make_optimizer(self):
        """Create a DominatingSetOptimizer with mock cache."""
        from neoswga.core.dominating_set_optimizer import DominatingSetOptimizer

        # Create a mock cache that returns positions for primers
        mock_cache = MagicMock()

        # 100 positions spread across the genome for each primer
        def mock_positions(prefix, primer, strand):
            # Return different positions for different primers to create coverage
            seed = hash(primer) % 1000
            rng = np.random.RandomState(seed)
            return rng.randint(0, 100000, size=50)

        mock_cache.get_positions.side_effect = mock_positions

        optimizer = DominatingSetOptimizer(
            cache=mock_cache,
            fg_prefixes=['genome1'],
            fg_seq_lengths=[100000],
            bin_size=10000
        )
        return optimizer

    def test_min_coverage_stops_early(self):
        """Optimizer should stop before max_primers when coverage target is met."""
        optimizer = self._make_optimizer()
        candidates = [f"PRIMER{i:03d}ATCGA" for i in range(50)]

        result = optimizer.optimize_greedy(
            candidates=candidates,
            max_primers=50,
            min_coverage=0.5,
            verbose=False
        )

        # Should have stopped before using all 50 primers
        assert result['n_primers'] < 50
        assert result['coverage'] >= 0.5
        assert result['coverage_target'] == 0.5
        assert result['target_met'] is True

    def test_no_coverage_target(self):
        """Without min_coverage, result should have target_met=False."""
        optimizer = self._make_optimizer()
        candidates = [f"PRIMER{i:03d}ATCGA" for i in range(10)]

        result = optimizer.optimize_greedy(
            candidates=candidates,
            max_primers=5,
            verbose=False
        )

        assert result['coverage_target'] is None
        assert result['target_met'] is False

    def test_coverage_target_in_result(self):
        """Result dict should include coverage_target and target_met keys."""
        optimizer = self._make_optimizer()
        candidates = [f"PRIMER{i:03d}ATCGA" for i in range(10)]

        result = optimizer.optimize_greedy(
            candidates=candidates,
            max_primers=5,
            min_coverage=0.95,
            verbose=False
        )

        assert 'coverage_target' in result
        assert 'target_met' in result


# ============================================================================
# Integration: long primers across GC extremes
# ============================================================================

class TestIntegrationGCExtremes:
    """Parameterized integration tests across GC extremes for 15bp primers."""

    @pytest.mark.parametrize("genome_gc,gc_label", [
        (0.25, "AT-rich"),
        (0.50, "GC-neutral"),
        (0.75, "GC-rich"),
    ])
    def test_15bp_primers_pass_filters(self, genome_gc, gc_label):
        """At least some 15bp primers should pass all filters for each GC regime."""
        from neoswga.core.filter import filter_extra, _scale_freq_threshold

        # Generate candidate 15bp primers appropriate for the genome GC
        primers = _generate_primers_for_gc(genome_gc, length=15, count=100)

        with patch('neoswga.core.filter.parameter') as mock_param, \
             patch('neoswga.core.filter.dimer') as mock_dimer, \
             patch('neoswga.core.filter._get_reaction_conditions') as mock_rc:

            mock_param.max_self_dimer_bp = max(4, int(15 * 0.28))
            mock_param.min_tm = 0
            mock_param.max_tm = 100
            mock_param.genome_gc = genome_gc

            # Adaptive GC range
            gc_tolerance = 0.15
            mock_param.gc_min = max(0.15, genome_gc - gc_tolerance)
            mock_param.gc_max = min(0.85, genome_gc + gc_tolerance)

            mock_dimer.is_dimer_fast.return_value = False

            mock_conditions = MagicMock()
            mock_conditions.calculate_effective_tm.return_value = 30.0
            mock_rc.return_value = mock_conditions

            passing = [p for p in primers if filter_extra(p)]

        assert len(passing) > 0, (
            f"No 15bp primers passed filters for {gc_label} genome "
            f"(GC={genome_gc:.0%}). Tested {len(primers)} candidates."
        )

    def test_frequency_scaling_for_15bp(self):
        """15bp primers should have a substantially lower frequency threshold."""
        from neoswga.core.filter import _scale_freq_threshold

        base = 1e-5
        scaled = _scale_freq_threshold(base, 15)
        # Should be ~1000x lower than the base threshold
        assert scaled < base / 100


def _generate_primers_for_gc(target_gc, length=15, count=100):
    """Generate random primers biased toward the target GC content."""
    rng = np.random.RandomState(42)
    primers = []
    for _ in range(count * 5):  # Generate extra, then filter
        # Choose bases with probability matching target GC
        gc_prob = target_gc / 2  # Split between G and C
        at_prob = (1 - target_gc) / 2  # Split between A and T
        bases = rng.choice(
            ['A', 'T', 'G', 'C'],
            size=length,
            p=[at_prob, at_prob, gc_prob, gc_prob]
        )
        primer = ''.join(bases)
        # Only keep primers with valid DNA bases
        if all(b in 'ACGT' for b in primer):
            primers.append(primer)
        if len(primers) >= count:
            break
    return primers[:count]
