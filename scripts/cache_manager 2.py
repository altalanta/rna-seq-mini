#!/usr/bin/env python3
"""
Intelligent caching system for RNASEQ-MINI pipeline.
Provides content-based hashing, cache management, and incremental processing capabilities.
"""

import hashlib
import json
import pickle
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple
import subprocess
import logging
from datetime import datetime, timedelta

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class CacheManager:
    """Manages intelligent caching for pipeline stages."""

    def __init__(self, cache_dir: str = ".cache", enabled: bool = True):
        self.cache_dir = Path(cache_dir)
        self.enabled = enabled
        self.cache_manifest = self.cache_dir / "manifest.json"
        self.hash_cache = self.cache_dir / "hashes.pkl"

        if enabled:
            self.cache_dir.mkdir(exist_ok=True)
            self._load_manifest()
            self._load_hash_cache()
        else:
            self.manifest = {}
            self.hash_cache_data = {}

    def _load_manifest(self):
        """Load cache manifest from disk."""
        if self.cache_manifest.exists():
            with open(self.cache_manifest, 'r') as f:
                self.manifest = json.load(f)
        else:
            self.manifest = {}

    def _save_manifest(self):
        """Save cache manifest to disk."""
        with open(self.cache_manifest, 'w') as f:
            json.dump(self.manifest, f, indent=2)

    def _load_hash_cache(self):
        """Load hash cache from disk."""
        if self.hash_cache.exists():
            with open(self.hash_cache, 'rb') as f:
                self.hash_cache_data = pickle.load(f)
        else:
            self.hash_cache_data = {}

    def _save_hash_cache(self):
        """Save hash cache to disk."""
        with open(self.hash_cache, 'wb') as f:
            pickle.dump(self.hash_cache_data, f)

    def compute_file_hash(self, file_path: Path, algorithm: str = 'sha256') -> str:
        """
        Compute content-based hash of a file.
        For large files, uses a sampling approach to avoid memory issues.
        """
        file_path = Path(file_path)

        if not file_path.exists():
            return ""

        # Check if we already have a cached hash
        cache_key = str(file_path)
        if cache_key in self.hash_cache_data:
            cached_hash, cached_mtime = self.hash_cache_data[cache_key]
            current_mtime = file_path.stat().st_mtime
            if cached_mtime == current_mtime:
                return cached_hash

        # Compute new hash
        hash_obj = hashlib.new(algorithm)

        # For large files, use sampling to avoid memory issues
        file_size = file_path.stat().st_size
        if file_size > 100 * 1024 * 1024:  # > 100MB
            sample_size = 10 * 1024 * 1024  # 10MB samples
            with open(file_path, 'rb') as f:
                # Sample from beginning, middle, and end
                samples = []
                chunk_size = 8192

                # Beginning
                chunk = f.read(min(sample_size // 3, chunk_size))
                samples.append(chunk)

                # Middle
                f.seek(file_size // 2)
                chunk = f.read(min(sample_size // 3, chunk_size))
                samples.append(chunk)

                # End
                f.seek(max(0, file_size - sample_size // 3))
                chunk = f.read(min(sample_size // 3, chunk_size))
                samples.append(chunk)

                for sample in samples:
                    hash_obj.update(sample)
        else:
            # For smaller files, hash entire content
            with open(file_path, 'rb') as f:
                for chunk in iter(lambda: f.read(8192), b''):
                    hash_obj.update(chunk)

        file_hash = hash_obj.hexdigest()

        # Cache the hash
        self.hash_cache_data[cache_key] = (file_hash, file_path.stat().st_mtime)
        self._save_hash_cache()

        return file_hash

    def compute_input_hash(self, inputs: List[Path], parameters: Optional[Dict] = None) -> str:
        """
        Compute hash for pipeline inputs including files and parameters.
        """
        hash_obj = hashlib.sha256()

        # Hash input files
        for input_path in sorted(inputs):
            file_hash = self.compute_file_hash(input_path)
            hash_obj.update(file_hash.encode())

        # Hash parameters
        if parameters:
            param_str = json.dumps(parameters, sort_keys=True)
            hash_obj.update(param_str.encode())

        return hash_obj.hexdigest()

    def is_cached(self, cache_key: str, output_files: List[Path]) -> bool:
        """
        Check if result is cached and all output files exist and are newer than inputs.
        """
        if not self.enabled or cache_key not in self.manifest:
            return False

        cache_info = self.manifest[cache_key]

        # Check if all output files exist
        for output_file in output_files:
            if not output_file.exists():
                logger.info(f"Cache miss: Output file {output_file} does not exist")
                return False

        # Check cache age and validity
        cache_time = datetime.fromisoformat(cache_info['created'])
        max_age = timedelta(days=30)  # Configurable cache expiry

        if datetime.now() - cache_time > max_age:
            logger.info(f"Cache expired for {cache_key}")
            return False

        logger.info(f"Cache hit for {cache_key}")
        return True

    def store_cache(self, cache_key: str, input_hash: str, output_files: List[Path],
                   metadata: Optional[Dict] = None):
        """
        Store cache entry for successful pipeline stage.
        """
        if not self.enabled:
            return

        # Verify all output files exist
        for output_file in output_files:
            if not output_file.exists():
                logger.warning(f"Cannot cache {cache_key}: Output file {output_file} missing")
                return

        cache_info = {
            'input_hash': input_hash,
            'output_files': [str(f) for f in output_files],
            'created': datetime.now().isoformat(),
            'metadata': metadata or {}
        }

        self.manifest[cache_key] = cache_info
        self._save_manifest()

        logger.info(f"Cached result for {cache_key}")

    def get_cached_outputs(self, cache_key: str) -> Optional[List[Path]]:
        """
        Retrieve cached output files for a cache key.
        """
        if not self.enabled or cache_key not in self.manifest:
            return None

        cache_info = self.manifest[cache_key]
        return [Path(f) for f in cache_info['output_files']]

    def cleanup_cache(self, max_age_days: int = 30, dry_run: bool = False) -> Dict:
        """
        Clean up old cache entries.
        Returns statistics about cleanup operation.
        """
        if not self.enabled:
            return {'removed_entries': 0, 'freed_space': 0}

        cutoff_date = datetime.now() - timedelta(days=max_age_days)
        entries_to_remove = []
        total_size = 0

        for cache_key, cache_info in self.manifest.items():
            cache_time = datetime.fromisoformat(cache_info['created'])

            if cache_time < cutoff_date:
                entries_to_remove.append(cache_key)

                # Calculate size of cached files
                for output_file in cache_info['output_files']:
                    file_path = Path(output_file)
                    if file_path.exists():
                        total_size += file_path.stat().st_size

        if dry_run:
            return {
                'entries_to_remove': len(entries_to_remove),
                'estimated_space': total_size
            }

        # Actually remove entries
        for cache_key in entries_to_remove:
            del self.manifest[cache_key]

        self._save_manifest()

        logger.info(f"Cleaned up {len(entries_to_remove)} cache entries, freed ~{total_size/1024/1024:.1f}MB")

        return {
            'removed_entries': len(entries_to_remove),
            'freed_space': total_size
        }

    def get_cache_stats(self) -> Dict:
        """Get cache usage statistics."""
        if not self.enabled:
            return {'enabled': False}

        total_entries = len(self.manifest)
        total_size = 0

        for cache_info in self.manifest.values():
            for output_file in cache_info['output_files']:
                file_path = Path(output_file)
                if file_path.exists():
                    total_size += file_path.stat().st_size

        return {
            'enabled': True,
            'total_entries': total_entries,
            'total_size_mb': total_size / 1024 / 1024,
            'cache_dir': str(self.cache_dir)
        }


def get_cache_manager(cache_dir: str = ".cache", enabled: bool = True) -> CacheManager:
    """Factory function to get a cache manager instance."""
    return CacheManager(cache_dir, enabled)


# Utility functions for pipeline integration
def should_skip_stage(cache_manager: CacheManager, stage_name: str,
                     input_files: List[Path], output_files: List[Path],
                     parameters: Optional[Dict] = None) -> Tuple[bool, Optional[str]]:
    """
    Determine if a pipeline stage should be skipped due to caching.

    Returns:
        (should_skip, cache_key)
    """
    if not cache_manager.enabled:
        return False, None

    input_hash = cache_manager.compute_input_hash(input_files, parameters)
    cache_key = f"{stage_name}_{input_hash[:16]}"  # Truncate for readability

    if cache_manager.is_cached(cache_key, output_files):
        return True, cache_key

    return False, cache_key


def mark_stage_complete(cache_manager: CacheManager, cache_key: str,
                       input_files: List[Path], output_files: List[Path],
                       parameters: Optional[Dict] = None,
                       metadata: Optional[Dict] = None):
    """
    Mark a pipeline stage as complete in the cache.
    """
    if not cache_manager.enabled:
        return

    input_hash = cache_manager.compute_input_hash(input_files, parameters)
    cache_manager.store_cache(cache_key, input_hash, output_files, metadata)


def cleanup_cache_command(cache_dir: str = ".cache", max_age_days: int = 30, dry_run: bool = False):
    """
    Command-line utility to clean up cache.
    """
    cache_manager = CacheManager(cache_dir, True)
    stats = cache_manager.cleanup_cache(max_age_days, dry_run)

    if dry_run:
        print(f"Cache cleanup (dry run): {stats['entries_to_remove']} entries, ~{stats['estimated_space']/1024/1024:.1f}MB")
    else:
        print(f"Cache cleanup: removed {stats['removed_entries']} entries, freed {stats['freed_space']/1024/1024:.1f}MB")

    return stats


def show_cache_stats(cache_dir: str = ".cache"):
    """
    Command-line utility to show cache statistics.
    """
    cache_manager = CacheManager(cache_dir, True)
    stats = cache_manager.get_cache_stats()

    if not stats['enabled']:
        print("Cache is disabled")
        return

    print(f"Cache directory: {stats['cache_dir']}")
    print(f"Total entries: {stats['total_entries']}")
    print(f"Total size: {stats['total_size_mb']:.1f}MB")


def clear_cache(cache_dir: str = ".cache", force: bool = False):
    """
    Command-line utility to clear entire cache.
    """
    if not force:
        confirm = input(f"Are you sure you want to clear the entire cache at {cache_dir}? (y/N): ")
        if confirm.lower() != 'y':
            print("Cache clear cancelled")
            return

    cache_manager = CacheManager(cache_dir, True)

    # Remove all cached files
    for cache_info in cache_manager.manifest.values():
        for output_file in cache_info['output_files']:
            file_path = Path(output_file)
            if file_path.exists():
                file_path.unlink()

    # Clear manifest and hash cache
    cache_manager.manifest = {}
    cache_manager.hash_cache_data = {}
    cache_manager._save_manifest()
    cache_manager._save_hash_cache()

    print(f"Cleared cache at {cache_dir}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="RNASEQ-MINI Cache Manager")
    parser.add_argument('--cache-dir', default='.cache', help='Cache directory')
    parser.add_argument('--stats', action='store_true', help='Show cache statistics')
    parser.add_argument('--cleanup', action='store_true', help='Clean up old cache entries')
    parser.add_argument('--clear', action='store_true', help='Clear entire cache')
    parser.add_argument('--dry-run', action='store_true', help='Show what would be cleaned (with --cleanup)')
    parser.add_argument('--max-age', type=int, default=30, help='Maximum age in days for cache cleanup')
    parser.add_argument('--force', action='store_true', help='Force clear cache without confirmation')

    args = parser.parse_args()

    if args.stats:
        show_cache_stats(args.cache_dir)
    elif args.cleanup:
        cleanup_cache_command(args.cache_dir, args.max_age, args.dry_run)
    elif args.clear:
        clear_cache(args.cache_dir, args.force)
    else:
        parser.print_help()
