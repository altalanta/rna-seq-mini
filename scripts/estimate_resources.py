#!/usr/bin/env python3
"""
Estimates the appropriate resource profile for a given task based on input file size.
"""

import argparse
from pathlib import Path
import yaml

def get_resource_profile(fastq_path: Path, config_path: Path = Path("config/params.yaml")) -> str:
    """
    Determines the resource profile (small, medium, large) based on FASTQ file size.
    """
    if not fastq_path.exists():
        # Return a default if the file doesn't exist to avoid errors in dry-runs
        return "small"

    file_size = fastq_path.stat().st_size
    
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    thresholds = config["resource_profiles"]["size_thresholds"]
    
    if file_size < thresholds["small"]:
        return "small"
    elif file_size < thresholds["medium"]:
        return "medium"
    else:
        return "large"

def main():
    """Main entry point for the script."""
    parser = argparse.ArgumentParser(description="Estimate resource profile based on FASTQ file size.")
    parser.add_argument("fastq_path", type=Path, help="Path to the input FASTQ file.")
    parser.add_argument("--config", default="config/params.yaml", type=Path, help="Path to the params.yaml config file.")
    args = parser.parse_args()
    
    profile = get_resource_profile(args.fastq_path, args.config)
    print(profile)

if __name__ == "__main__":
    main()
