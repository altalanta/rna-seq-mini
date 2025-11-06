#!/usr/bin/env python3
"""
Validates the rnaseq-mini configuration files using the centralized Pydantic models.
"""

import argparse
import sys
from pathlib import Path
from pydantic import ValidationError

# Add src to the Python path to allow importing the validation module
sys.path.append(str(Path(__file__).resolve().parent.parent / "src"))

from rnaseq_mini.validation import load_and_validate_config

def main():
    """Main entry point for the validation script."""
    parser = argparse.ArgumentParser(description="Configuration validation for RNASEQ-MINI")
    parser.add_argument('config_file', nargs='?', default='config/params.yaml', help='Main configuration file to validate (params.yaml)')
    args = parser.parse_args()

    try:
        print(f"🔍 Validating configuration using {args.config_file}...")
        load_and_validate_config(Path(args.config_file))
        print("✅ Configuration is valid!")
        return 0
    except FileNotFoundError as e:
        print(f"❌ Validation Error: A required file was not found.")
        print(f"   File: {e.filename}")
        return 1
    except ValidationError as e:
        print(f"❌ Configuration validation failed with {len(e.errors())} error(s):")
        for error in e.errors():
            loc = " -> ".join(map(str, error['loc']))
            print(f"  - Location: {loc}")
            print(f"    Message: {error['msg']}")
        return 1
    except ValueError as e:
        print(f"❌ Cross-validation failed:")
        print(f"   {e}")
        return 1
    except Exception as e:
        print(f"An unexpected error occurred during validation: {e}")
        return 1

if __name__ == "__main__":
    sys.exit(main())
