#!/usr/bin/env python3
"""
Intelligent caching system for RNASEQ-MINI pipeline.
Provides content-based hashing, cache management, and incremental processing capabilities.
"""

import hashlib
import json
from pathlib import Path
import shutil
import os

CACHE_DIR = Path(".cache")

def get_cache_dir() -> Path:
    """Returns the cache directory, creating it if it doesn't exist."""
    CACHE_DIR.mkdir(exist_ok=True)
    return CACHE_DIR

def compute_hash(input_files: list[Path], params: dict) -> str:
    """
    Computes a deterministic SHA256 hash for a stage based on its inputs and parameters.
    """
    hasher = hashlib.sha256()

    # 1. Hash based on the content of input files
    for file_path in sorted(input_files):
        if file_path.exists():
            with open(file_path, "rb") as f:
                while chunk := f.read(8192):
                    hasher.update(chunk)
    
    # 2. Hash based on the JSON representation of the parameters
    params_str = json.dumps(params, sort_keys=True)
    hasher.update(params_str.encode())
    
    return hasher.hexdigest()

def check_cache(cache_hash: str) -> bool:
    """Checks if a result exists in the cache for a given hash."""
    cache_path = get_cache_dir() / cache_hash
    return cache_path.exists()

def retrieve_from_cache(cache_hash: str, output_path: Path):
    """
    Retrieves a result from the cache by creating a symbolic link.
    Assumes the output path's parent directory exists.
    """
    cache_path = get_cache_dir() / cache_hash
    if output_path.exists() or output_path.is_symlink():
        output_path.unlink()
    
    # We create a symlink from the cache to the expected output path
    os.symlink(cache_path.resolve(), output_path)

def store_in_cache(cache_hash: str, source_path: Path):
    """
    Stores a new result in the cache by copying it.
    """
    cache_path = get_cache_dir() / cache_hash
    if not cache_path.exists():
        shutil.copy2(source_path, cache_path)
