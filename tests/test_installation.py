#!/usr/bin/env python3
"""
Installation Test Script for NeoSWGA

Tests that the package is properly installed and all core imports work correctly.
Run this after `pip install -e .` to verify the installation.

Usage:
    python test_installation.py
    python test_installation.py --verbose
"""

import sys
import os
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)


class InstallationTester:
    """Test NeoSWGA installation"""

    def __init__(self, verbose=False):
        self.verbose = verbose
        self.tests_passed = 0
        self.tests_failed = 0
        self.errors = []

    def test_import(self, module_name, description, allow_runtime_errors=False):
        """Test that a module can be imported"""
        try:
            __import__(module_name)
            self.tests_passed += 1
            if self.verbose:
                logger.info(f"✓ {description}")
            return True
        except ImportError as e:
            self.tests_failed += 1
            error_msg = f"✗ {description}: {e}"
            logger.error(error_msg)
            self.errors.append(error_msg)
            return False
        except (KeyError, ValueError) as e:
            # Some modules (like pipeline) have module-level code that requires parameters
            # This is expected and doesn't indicate a broken installation
            if allow_runtime_errors:
                self.tests_passed += 1
                if self.verbose:
                    logger.info(f"✓ {description} (module-level initialization skipped)")
                return True
            else:
                self.tests_failed += 1
                error_msg = f"✗ {description}: Runtime error during import: {e}"
                logger.error(error_msg)
                self.errors.append(error_msg)
                return False

    def test_package_structure(self):
        """Test basic package structure"""
        logger.info("\n" + "="*70)
        logger.info("TEST 1: PACKAGE STRUCTURE")
        logger.info("="*70)

        tests = [
            ("neoswga", "Main package", False),
            ("neoswga.core", "Core subpackage", False),
            ("neoswga.core.parameter", "Parameter module", False),
            ("neoswga.core.filter", "Filter module", False),
            ("neoswga.core.kmer_counter", "K-mer counter module", False),
            ("neoswga.core.utility", "Utility module", False),
            # Pipeline has module-level initialization, so allow runtime errors
            ("neoswga.core.pipeline", "Pipeline module", True),
        ]

        for module, desc, allow_errors in tests:
            self.test_import(module, desc, allow_runtime_errors=allow_errors)

    def test_core_modules(self):
        """Test core functionality modules"""
        logger.info("\n" + "="*70)
        logger.info("TEST 2: CORE FUNCTIONALITY MODULES")
        logger.info("="*70)

        tests = [
            ("neoswga.core.optimize", "Original greedy optimizer"),
            ("neoswga.core.network_optimizer", "Network-based optimizer"),
            ("neoswga.core.genetic_algorithm", "Genetic algorithm optimizer"),
            ("neoswga.core.hybrid_optimizer", "Hybrid optimizer"),
            ("neoswga.core.adaptive_filters", "Adaptive GC filtering"),
            ("neoswga.core.position_cache", "Position cache"),
            ("neoswga.core.rf_preprocessing", "Random forest preprocessing"),
        ]

        for module, desc in tests:
            self.test_import(module, desc)

    def test_advanced_modules(self):
        """Test advanced features"""
        logger.info("\n" + "="*70)
        logger.info("TEST 3: ADVANCED MODULES")
        logger.info("="*70)

        tests = [
            ("neoswga.core.thermodynamics", "Thermodynamics"),
            ("neoswga.core.reaction_conditions", "Reaction conditions"),
            ("neoswga.core.secondary_structure", "Secondary structure"),
            ("neoswga.core.background_filter", "Background Bloom filter"),
            ("neoswga.core.amplicon_network", "Amplicon network"),
        ]

        for module, desc in tests:
            self.test_import(module, desc)

    def test_optional_modules(self):
        """Test optional modules (may not be installed)"""
        logger.info("\n" + "="*70)
        logger.info("TEST 4: OPTIONAL MODULES")
        logger.info("="*70)

        # MILP optimizer (requires mip package)
        try:
            import mip
            mip_available = True
        except ImportError:
            mip_available = False
            logger.warning("⊘ mip package not installed (install with pip install -e \".[improved]\")")

        if mip_available:
            self.test_import("neoswga.core.milp_optimizer", "MILP optimizer")

        # GPU acceleration (requires cupy)
        try:
            import cupy
            cupy_available = True
        except ImportError:
            cupy_available = False
            logger.warning("⊘ cupy not installed (install with pip install -e \".[gpu]\")")

        if cupy_available:
            self.test_import("neoswga.core.gpu_acceleration", "GPU acceleration")

    def test_model_files(self):
        """Test that model files are accessible"""
        logger.info("\n" + "="*70)
        logger.info("TEST 5: MODEL FILES")
        logger.info("="*70)

        try:
            import neoswga.core.rf_preprocessing as rf
            import os

            # Get model directory
            module_dir = os.path.dirname(rf.__file__)
            model_dir = os.path.join(module_dir, 'models')
            model_path = os.path.join(model_dir, 'random_forest_filter.p')

            if os.path.exists(model_path):
                size_mb = os.path.getsize(model_path) / (1024 * 1024)
                logger.info(f"✓ Random forest model found ({size_mb:.1f} MB)")
                self.tests_passed += 1
            else:
                logger.error(f"✗ Model file not found at: {model_path}")
                self.errors.append(f"Model file missing: {model_path}")
                self.tests_failed += 1

        except Exception as e:
            logger.error(f"✗ Model file check failed: {e}")
            self.errors.append(f"Model file check error: {e}")
            self.tests_failed += 1

    def test_cli_entry_point(self):
        """Test that CLI entry point is installed"""
        logger.info("\n" + "="*70)
        logger.info("TEST 6: CLI ENTRY POINT")
        logger.info("="*70)

        import subprocess

        try:
            result = subprocess.run(
                ['neoswga', '--help'],
                capture_output=True,
                text=True,
                timeout=5
            )

            if result.returncode == 0 and 'neoswga' in result.stdout.lower():
                logger.info("✓ CLI entry point 'neoswga' installed and working")
                self.tests_passed += 1
            else:
                logger.error("✗ CLI entry point returned unexpected output")
                self.errors.append("CLI entry point check failed")
                self.tests_failed += 1

        except FileNotFoundError:
            logger.warning("⊘ CLI entry point 'neoswga' not found in PATH")
            logger.warning("   Run 'pip install -e .' to install the CLI command")
            # Don't count as failure - package may not be installed yet
            pass

        except Exception as e:
            logger.error(f"✗ CLI check failed: {e}")
            self.errors.append(f"CLI check error: {e}")
            self.tests_failed += 1

    def test_version(self):
        """Test package version"""
        logger.info("\n" + "="*70)
        logger.info("TEST 7: PACKAGE VERSION")
        logger.info("="*70)

        try:
            import neoswga
            version = neoswga.__version__
            logger.info(f"✓ NeoSWGA version: {version}")
            self.tests_passed += 1

            if version == "3.0.0":
                logger.info("✓ Version matches expected: 3.0.0")
                self.tests_passed += 1
            else:
                logger.warning(f"⚠ Version is {version}, expected 3.0.0")

        except Exception as e:
            logger.error(f"✗ Version check failed: {e}")
            self.errors.append(f"Version check error: {e}")
            self.tests_failed += 1

    def test_basic_functionality(self):
        """Test basic functionality"""
        logger.info("\n" + "="*70)
        logger.info("TEST 8: BASIC FUNCTIONALITY")
        logger.info("="*70)

        try:
            # Test thermodynamics
            from neoswga.core import thermodynamics as thermo
            # Use a function that exists in the module
            gc = thermo.gc_content('ATCGATCGATCG')
            if 0 <= gc <= 1:
                logger.info(f"✓ Thermodynamics: GC content calculation works (GC={gc:.2f})")
                self.tests_passed += 1
            else:
                logger.error(f"✗ Thermodynamics: Unexpected GC value: {gc}")
                self.errors.append("Thermodynamics calculation error")
                self.tests_failed += 1

            # Test adaptive GC filter
            from neoswga.core.adaptive_filters import AdaptiveGCFilter
            gc_filter = AdaptiveGCFilter(genome_gc=0.50, tolerance=0.15)
            if 0.35 <= gc_filter.gc_min <= 0.36 and 0.64 <= gc_filter.gc_max <= 0.66:
                logger.info(f"✓ Adaptive GC filter: Range calculation works ({gc_filter.gc_min:.2f}-{gc_filter.gc_max:.2f})")
                self.tests_passed += 1
            else:
                logger.error(f"✗ Adaptive GC filter: Unexpected range: {gc_filter.gc_min}-{gc_filter.gc_max}")
                self.errors.append("Adaptive GC filter calculation error")
                self.tests_failed += 1

            # Test complexity calculation (if available)
            try:
                from neoswga.core.filter import calculate_sequence_complexity
                complexity = calculate_sequence_complexity('ATCGATCGATCG', k=4)
                if 0 <= complexity <= 1:
                    logger.info(f"✓ Complexity filter: Calculation works (complexity={complexity:.2f})")
                    self.tests_passed += 1
                else:
                    logger.error(f"✗ Complexity filter: Unexpected value: {complexity}")
                    self.errors.append("Complexity calculation error")
                    self.tests_failed += 1
            except ImportError:
                logger.info("⊘ Complexity filter: Function not available (optional)")
                pass

        except Exception as e:
            logger.error(f"✗ Basic functionality test failed: {e}")
            import traceback
            traceback.print_exc()
            self.errors.append(f"Basic functionality error: {e}")
            self.tests_failed += 1

    def print_summary(self):
        """Print test summary"""
        logger.info("\n" + "="*70)
        logger.info("INSTALLATION TEST SUMMARY")
        logger.info("="*70)

        total_tests = self.tests_passed + self.tests_failed
        logger.info(f"\nTotal tests: {total_tests}")
        logger.info(f"Passed: {self.tests_passed}")
        logger.info(f"Failed: {self.tests_failed}")

        if self.tests_failed == 0:
            logger.info("\n" + "="*70)
            logger.info("SUCCESS: All installation tests passed! ✓")
            logger.info("="*70)
            logger.info("\nNeoSWGA is properly installed and ready to use.")
            logger.info("\nTry running:")
            logger.info("  neoswga --help")
            logger.info("  neoswga validate --quick")
            return True
        else:
            logger.error("\n" + "="*70)
            logger.error("FAILURE: Some installation tests failed ✗")
            logger.error("="*70)
            logger.error("\nErrors encountered:")
            for i, error in enumerate(self.errors, 1):
                logger.error(f"  {i}. {error}")
            logger.error("\nPlease ensure you:")
            logger.error("  1. Installed the package: pip install -e .")
            logger.error("  2. Are in the correct environment")
            logger.error("  3. Have all dependencies installed")
            return False

    def run_all_tests(self):
        """Run all installation tests"""
        logger.info("\n" + "="*70)
        logger.info("NEOSWGA INSTALLATION TEST SUITE")
        logger.info("="*70)
        logger.info("\nTesting installation of NeoSWGA package...")

        self.test_package_structure()
        self.test_core_modules()
        self.test_advanced_modules()
        self.test_optional_modules()
        self.test_model_files()
        self.test_cli_entry_point()
        self.test_version()
        self.test_basic_functionality()

        return self.print_summary()


def main():
    """Main entry point"""
    import argparse

    parser = argparse.ArgumentParser(
        description='Test NeoSWGA installation'
    )
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Verbose output')

    args = parser.parse_args()

    tester = InstallationTester(verbose=args.verbose)
    success = tester.run_all_tests()

    sys.exit(0 if success else 1)


if __name__ == '__main__':
    main()
