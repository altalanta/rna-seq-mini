#!/usr/bin/env python3
"""
Config validation script for RNASEQ-MINI.
Validates config/params.yaml against a formal JSON schema.
"""

import yaml
import sys
import json
import argparse
from pathlib import Path
from jsonschema import validate, exceptions

def validate_config_against_schema(config_path, schema_path):
    """Validate a YAML config file against a JSON schema."""
    try:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
    except FileNotFoundError:
        print(f"❌ Error: Config file '{config_path}' not found.")
        return False, None
    except yaml.YAMLError as e:
        print(f"❌ Error parsing YAML in '{config_path}': {e}")
        return False, None

    try:
        with open(schema_path, 'r') as f:
            schema = json.load(f)
    except FileNotFoundError:
        print(f"❌ Error: Schema file '{schema_path}' not found.")
        return False, None
    except json.JSONDecodeError as e:
        print(f"❌ Error parsing JSON in '{schema_path}': {e}")
        return False, None
        
    try:
        validate(instance=config, schema=schema)
        print(f"✅ Schema validation passed for '{config_path}'.")
        return True, config
    except exceptions.ValidationError as e:
        print(f"❌ Schema validation failed for '{config_path}':")
        print(f"  - Message: {e.message}")
        print(f"  - Path: {list(e.path)}")
        print(f"  - Validator: {e.validator} = {e.validator_value}")
        return False, None

def check_paths(config):
    """Check that all file paths in the config exist."""
    print("🔍 Checking file paths...")
    all_paths_exist = True
    
    paths_to_check = [
        config['reference']['transcripts_fa'],
        config['reference']['annotation_gtf'],
        config['paths']['samples'],
        config['r']['contrasts_file']
    ]
    
    # Add decoy fasta if it's specified
    if config.get('reference', {}).get('decoy_fasta'):
        paths_to_check.append(config['reference']['decoy_fasta'])

    for file_path in paths_to_check:
        if not Path(file_path).exists():
            print(f"  - ❌ Missing file: {file_path}")
            all_paths_exist = False
        else:
            print(f"  - ✅ Found file: {file_path}")
            
    return all_paths_exist

def main():
    parser = argparse.ArgumentParser(description="Configuration validation for RNASEQ-MINI")
    parser.add_argument('config_file', nargs='?', default='config/params.yaml', help='Configuration file to validate')
    parser.add_argument('--schema', default='config/schema.json', help='JSON schema to validate against')
    args = parser.parse_args()

    schema_valid, config = validate_config_against_schema(args.config_file, args.schema)
    
    if not schema_valid:
        return 1

    paths_valid = check_paths(config)

    if not paths_valid:
        print("\n❌ Some file paths in the configuration do not exist.")
        return 1
        
    print("\n🎉 Configuration is valid and all required files are present.")
    return 0

if __name__ == "__main__":
    sys.exit(main())

