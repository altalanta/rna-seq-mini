import pytest
import json
import yaml
from pathlib import Path

# Adjust the path to import from the src directory
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'src'))

from rnaseq_mini.validation import validate_config_against_schema, check_paths

@pytest.fixture
def temp_config_files(tmp_path):
    """Create temporary config, schema, and data files for testing."""
    schema = {
        "$schema": "http://json-schema.org/draft-07/schema#",
        "type": "object",
        "properties": {
            "project": {"type": "string"},
            "threads": {"type": "integer", "minimum": 1},
            "reference": {
                "type": "object",
                "properties": { "transcripts_fa": {"type": "string"} },
                "required": ["transcripts_fa"]
            }
        },
        "required": ["project", "threads", "reference"]
    }
    
    valid_config = {
        "project": "test-project",
        "threads": 4,
        "reference": { "transcripts_fa": str(tmp_path / "transcripts.fa") }
    }
    
    invalid_type_config = {
        "project": "test-project",
        "threads": "four", # Incorrect type
        "reference": { "transcripts_fa": str(tmp_path / "transcripts.fa") }
    }
    
    missing_key_config = {
        "project": "test-project",
        "reference": { "transcripts_fa": str(tmp_path / "transcripts.fa") }
    }

    # Write files
    schema_path = tmp_path / "schema.json"
    valid_config_path = tmp_path / "valid_config.yaml"
    invalid_type_path = tmp_path / "invalid_type.yaml"
    missing_key_path = tmp_path / "missing_key.yaml"
    
    with open(schema_path, 'w') as f:
        json.dump(schema, f)
    with open(valid_config_path, 'w') as f:
        yaml.dump(valid_config, f)
    with open(invalid_type_path, 'w') as f:
        yaml.dump(invalid_type_config, f)
    with open(missing_key_path, 'w') as f:
        yaml.dump(missing_key_config, f)
        
    # Create a dummy data file for the path check
    (tmp_path / "transcripts.fa").touch()

    return {
        "schema": schema_path,
        "valid": valid_config_path,
        "invalid_type": invalid_type_path,
        "missing_key": missing_key_path,
        "data_file": tmp_path / "transcripts.fa"
    }

def test_valid_config(temp_config_files):
    """Test that a valid config passes schema validation."""
    is_valid, config = validate_config_against_schema(
        temp_config_files["valid"], temp_config_files["schema"]
    )
    assert is_valid is True
    assert config is not None

def test_invalid_type_config(temp_config_files):
    """Test that a config with an incorrect data type fails validation."""
    is_valid, _ = validate_config_against_schema(
        temp_config_files["invalid_type"], temp_config_files["schema"]
    )
    assert is_valid is False

def test_missing_key_config(temp_config_files):
    """Test that a config with a missing required key fails validation."""
    is_valid, _ = validate_config_against_schema(
        temp_config_files["missing_key"], temp_config_files["schema"]
    )
    assert is_valid is False

def test_check_paths_success(temp_config_files):
    """Test that check_paths returns True when all files exist."""
    # We need a mock config object similar to what the validator would return
    with open(temp_config_files["valid"]) as f:
        config = yaml.safe_load(f)
    
    # Manually create mock paths for other checks in the function
    config['paths'] = {'samples': temp_config_files['data_file']}
    config['r'] = {'contrasts_file': temp_config_files['data_file']}

    assert check_paths(config) is True

def test_check_paths_failure():
    """Test that check_paths returns False if a file is missing."""
    # Create a config that points to a non-existent file
    config = {
        "reference": {"transcripts_fa": "nonexistent.fa"},
        "paths": {"samples": "nonexistent.tsv"},
        "r": {"contrasts_file": "nonexistent.tsv"}
    }
    assert check_paths(config) is False


