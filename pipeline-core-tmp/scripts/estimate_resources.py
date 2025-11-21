#!/usr/bin/env python3
"""
Estimates the appropriate resource profile for a given task based on input file size.
"""

import argparse
from pathlib import Path
import yaml
import json
import sys

def get_yaml_loader():
    """Returns a YAML loader that supports the !include directive."""
    from yaml.loader import SafeLoader

    class YamlIncludeLoader(SafeLoader):
        def __init__(self, stream):
            self._root = Path(stream.name).parent
            super().__init__(stream)

        def include(self, node):
            filename = self._root / self.construct_scalar(node)
            with open(filename, 'r') as f:
                return yaml.load(f, YamlIncludeLoader)

    YamlIncludeLoader.add_constructor("!include", YamlIncludeLoader.include)
    return YamlIncludeLoader

def get_resource_profile(fastq_path: Path, thresholds: dict) -> str:
    """
    Determines the resource profile (small, medium, large) based on FASTQ file size.
    """
    if not fastq_path.exists():
        return "small"

    file_size = fastq_path.stat().st_size
    
    if file_size < thresholds["small"]:
        return "small"
    elif file_size < thresholds["medium"]:
        return "medium"
    else:
        return "large"

def main():
    """Main entry point for the script."""
    parser = argparse.ArgumentParser(description="Estimate resource profiles for a list of FASTQ files.")
    parser.add_argument("sample_sheet", type=Path, help="Path to the sample sheet TSV file.")
    parser.add_argument("--config", default="../config/params.yaml", type=Path, help="Path to the params.yaml config file.")
    args = parser.parse_args()
    
    with open(args.config, 'r') as f:
        config = yaml.load(f, Loader=YamlIncludeLoader)
    
    thresholds = config["execution"]["resource_profiles"]["size_thresholds"]
    
    profiles = {}
    with open(args.sample_sheet, 'r') as f:
        # Skip header
        next(f)
        for line in f:
            fields = line.strip().split('\t')
            sample_id = fields[0]
            fastq_path = Path(fields[1])
            profiles[sample_id] = get_resource_profile(fastq_path, thresholds)
            
    print(json.dumps(profiles, indent=4))

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
