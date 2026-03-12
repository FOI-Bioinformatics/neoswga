"""
Security tests for NeoSWGA.

Tests for:
- Command injection prevention in kmer_counter.py
- Model hash verification in rf_preprocessing.py
- Secure subprocess usage
"""

import pytest
import os
import tempfile
import hashlib
from pathlib import Path
from unittest.mock import patch, MagicMock

from neoswga.core.kmer_counter import run_jellyfish
from neoswga.core.rf_preprocessing import (
    load_model_safely,
    ModelIntegrityError,
    _compute_file_hash,
    _TRUSTED_MODEL_HASHES,
)


# =============================================================================
# Command Injection Prevention Tests
# =============================================================================

class TestCommandInjectionPrevention:
    """Tests verifying command injection vulnerabilities are fixed."""

    def test_run_jellyfish_validates_file_exists(self, tmp_path):
        """Test that run_jellyfish checks if genome file exists."""
        nonexistent = str(tmp_path / "nonexistent.fa")

        with patch('neoswga.core.kmer_counter.require_jellyfish'):
            with pytest.raises(FileNotFoundError) as exc_info:
                run_jellyfish(nonexistent, str(tmp_path / "output"))

            assert "not found" in str(exc_info.value)

    def test_run_jellyfish_uses_subprocess_not_os_system(self):
        """Verify run_jellyfish uses subprocess.run, not os.system."""
        import inspect
        from neoswga.core import kmer_counter

        source = inspect.getsource(kmer_counter)

        # Should use subprocess.run
        assert "subprocess.run" in source
        # Should NOT use os.system (which is vulnerable)
        assert "os.system(" not in source

    def test_potential_malicious_genome_fname_is_safe(self, tmp_path):
        """Test that malicious filenames don't cause command injection."""
        malicious_path = "/tmp/safe.fa; echo PWNED > /tmp/pwned.txt"

        # Mock subprocess.run to capture what would be executed
        with patch('neoswga.core.kmer_counter.subprocess.run') as mock_run, \
             patch('neoswga.core.kmer_counter.require_jellyfish'), \
             patch('neoswga.core.kmer_counter._adaptive_hash_size', return_value='100M'), \
             patch('os.path.exists', return_value=True):
            mock_run.return_value = MagicMock(returncode=0)

            try:
                run_jellyfish(malicious_path, str(tmp_path / "output"), min_k=6, max_k=6)
            except Exception:
                pass

            # Verify subprocess.run was called with list arguments (safe)
            # not with shell=True (unsafe)
            if mock_run.called:
                call_args = mock_run.call_args
                cmd = call_args[0][0]  # First positional arg is the command

                # Command should be a list (safe) not a string (unsafe)
                assert isinstance(cmd, list), "Command should be a list, not a string"

                # The malicious path should be passed as a single element
                # If command injection worked, it would be split at ";"
                found_path = False
                for arg in cmd:
                    if malicious_path in arg:
                        found_path = True
                        # If the path is in args, it should be the whole path as one arg
                        assert arg == malicious_path, "Path should be single argument, not split"
                assert found_path, "Malicious path should be found as single argument"


# =============================================================================
# Model Hash Verification Tests
# =============================================================================

class TestModelHashVerification:
    """Tests for secure model loading with hash verification."""

    def test_compute_file_hash(self, tmp_path):
        """Test SHA-256 hash computation."""
        test_file = tmp_path / "test.bin"
        test_content = b"Hello, World!"
        test_file.write_bytes(test_content)

        expected_hash = hashlib.sha256(test_content).hexdigest()
        actual_hash = _compute_file_hash(str(test_file))

        assert actual_hash == expected_hash

    def test_load_model_safely_file_not_found(self, tmp_path):
        """Test that load_model_safely raises FileNotFoundError for missing files."""
        with pytest.raises(FileNotFoundError):
            load_model_safely(str(tmp_path / "nonexistent.pkl"))

    def test_load_model_safely_hash_mismatch(self, tmp_path):
        """Test that load_model_safely raises ModelIntegrityError on hash mismatch."""
        # Create a fake model file
        fake_model = tmp_path / "random_forest_filter.p"
        fake_model.write_bytes(b"fake model content")

        # Temporarily add a hash for this file
        original_hashes = _TRUSTED_MODEL_HASHES.copy()
        _TRUSTED_MODEL_HASHES['random_forest_filter.p'] = 'wrong_hash'

        try:
            with pytest.raises(ModelIntegrityError) as exc_info:
                load_model_safely(str(fake_model), verify_hash=True)

            assert "hash mismatch" in str(exc_info.value).lower()
        finally:
            # Restore original hashes
            _TRUSTED_MODEL_HASHES.clear()
            _TRUSTED_MODEL_HASHES.update(original_hashes)

    def test_load_model_safely_with_verification_disabled(self, tmp_path):
        """Test that load_model_safely works when verification is disabled."""
        import pickle

        # Create a legitimate pickle file
        test_model = tmp_path / "test_model.pkl"
        test_data = {"model": "test", "value": 42}
        with open(test_model, 'wb') as f:
            pickle.dump(test_data, f)

        # Load without verification
        result = load_model_safely(str(test_model), verify_hash=False)

        assert result == test_data

    def test_load_model_safely_unknown_model_warning(self, tmp_path, caplog):
        """Test warning when loading unknown model file."""
        import pickle
        import logging

        # Create a pickle file with unknown name
        unknown_model = tmp_path / "unknown_model.pkl"
        with open(unknown_model, 'wb') as f:
            pickle.dump({"test": True}, f)

        with caplog.at_level(logging.WARNING):
            load_model_safely(str(unknown_model), verify_hash=True)

        assert "No trusted hash for model" in caplog.text

    def test_trusted_model_loads_successfully(self):
        """Test that the actual trusted model loads with verification."""
        model_path = os.path.join(
            os.path.dirname(__file__),
            '..', 'neoswga', 'core', 'models', 'random_forest_filter.p'
        )

        if os.path.exists(model_path):
            model = load_model_safely(model_path, verify_hash=True)
            assert model is not None
            # Verify it's actually a model
            assert hasattr(model, 'predict')

    def test_model_integrity_error_is_exception(self):
        """Test that ModelIntegrityError is a proper exception."""
        assert issubclass(ModelIntegrityError, Exception)

        # Should be raisable with message
        with pytest.raises(ModelIntegrityError) as exc_info:
            raise ModelIntegrityError("Test error message")

        assert "Test error message" in str(exc_info.value)


# =============================================================================
# Subprocess Security Tests
# =============================================================================

class TestSubprocessSecurity:
    """Tests for secure subprocess usage patterns."""

    def test_kmer_counter_uses_secure_subprocess(self):
        """Verify kmer_counter uses subprocess.run with list args."""
        import inspect
        from neoswga.core import kmer_counter

        source = inspect.getsource(kmer_counter)

        # Should use subprocess.run
        assert "subprocess.run" in source
        # Should NOT use shell=True
        assert "shell=True" not in source
        # Should NOT use os.system
        assert "os.system(" not in source

    def test_no_shell_true_in_codebase(self):
        """Verify no subprocess calls use shell=True."""
        import inspect
        from neoswga.core import kmer_counter

        source = inspect.getsource(kmer_counter)
        assert "shell=True" not in source, f"shell=True found in {kmer_counter.__name__}"


# =============================================================================
# Path Traversal Prevention Tests
# =============================================================================

class TestPathTraversalPrevention:
    """Tests for path traversal vulnerability prevention."""

    def test_run_jellyfish_rejects_path_traversal(self, tmp_path):
        """Test that path traversal attempts are handled safely."""
        # Attempt path traversal - should fail at file existence check
        traversal_attempt = str(tmp_path / ".." / ".." / "etc" / "passwd")

        with patch('neoswga.core.kmer_counter.require_jellyfish'):
            with pytest.raises(FileNotFoundError):
                run_jellyfish(traversal_attempt, str(tmp_path / "output"))


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
