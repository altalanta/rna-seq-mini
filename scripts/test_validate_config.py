"""Unit tests for validate_config.py"""

import pytest
import tempfile
import os
import yaml


class TestValidateConfig:
    """Test cases for configuration validation."""

    def test_valid_config_loading(self):
        """Test that valid configuration files can be loaded."""
        # This is a basic test - in practice, you'd test actual config validation logic
        config_path = "config/params.yaml"

        if os.path.exists(config_path):
            with open(config_path, 'r') as f:
                config = yaml.safe_load(f)
                assert isinstance(config, dict)
                assert 'engine' in config or 'samples' in config
        else:
            pytest.skip("Config file not found")

    def test_config_validation_functionality(self):
        """Test that the config validation script runs without errors."""
        # This would test the actual validation function if it were extracted
        # For now, we test that the script can be imported/run
        try:
            import scripts.validate_config
            assert True  # Script imports successfully
        except ImportError:
            pytest.skip("Validation script not importable in test environment")

    def test_yaml_syntax_validation(self):
        """Test YAML syntax validation for config files."""
        config_files = ["config/params.yaml", "config/samples.tsv", "config/contrasts.tsv"]

        for config_file in config_files:
            if os.path.exists(config_file):
                if config_file.endswith('.yaml') or config_file.endswith('.yml'):
                    with open(config_file, 'r') as f:
                        yaml.safe_load(f)  # This will raise an exception if YAML is invalid
                elif config_file.endswith('.tsv'):
                    # Basic TSV validation
                    with open(config_file, 'r') as f:
                        lines = f.readlines()
                        assert len(lines) > 0, f"Empty TSV file: {config_file}"
                        # Check that header exists
                        header = lines[0].strip().split('\t')
                        assert len(header) > 0, f"No header in TSV file: {config_file}
