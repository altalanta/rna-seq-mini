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


class CacheManager:
    """
    Cache manager class for managing pipeline stage caching.
    Provides methods for cache statistics, cleanup, and stage tracking.
    """
    
    def __init__(self, cache_dir: Path = None, enabled: bool = True):
        self.cache_dir = Path(cache_dir) if cache_dir else CACHE_DIR
        self.enabled = enabled
        self.cache_dir.mkdir(exist_ok=True)
        self._stage_status_file = self.cache_dir / ".stage_status.json"
        self._stage_status = self._load_stage_status()
    
    def _load_stage_status(self) -> dict:
        """Load stage completion status from file."""
        if self._stage_status_file.exists():
            try:
                with open(self._stage_status_file, 'r') as f:
                    return json.load(f)
            except (json.JSONDecodeError, IOError):
                return {}
        return {}
    
    def _save_stage_status(self):
        """Save stage completion status to file."""
        with open(self._stage_status_file, 'w') as f:
            json.dump(self._stage_status, f, indent=2)
    
    def get_cache_stats(self) -> dict:
        """Get statistics about the cache."""
        if not self.cache_dir.exists():
            return {"total_files": 0, "total_size_mb": 0, "enabled": self.enabled}
        
        files = list(self.cache_dir.glob("*"))
        # Exclude hidden files like .stage_status.json
        cache_files = [f for f in files if f.is_file() and not f.name.startswith('.')]
        total_size = sum(f.stat().st_size for f in cache_files)
        
        return {
            "total_files": len(cache_files),
            "total_size_mb": round(total_size / (1024 * 1024), 2),
            "cache_dir": str(self.cache_dir),
            "enabled": self.enabled,
            "stages_completed": len(self._stage_status)
        }
    
    def cleanup_cache(self, max_age_days: int = 30, dry_run: bool = False) -> dict:
        """Clean up old cache entries."""
        import time
        
        now = time.time()
        max_age_seconds = max_age_days * 24 * 60 * 60
        
        files_removed = 0
        bytes_freed = 0
        
        if self.cache_dir.exists():
            for cache_file in self.cache_dir.glob("*"):
                if cache_file.is_file() and not cache_file.name.startswith('.'):
                    file_age = now - cache_file.stat().st_mtime
                    if file_age > max_age_seconds:
                        file_size = cache_file.stat().st_size
                        if not dry_run:
                            cache_file.unlink()
                        files_removed += 1
                        bytes_freed += file_size
        
        return {
            "files_removed": files_removed,
            "bytes_freed": bytes_freed,
            "mb_freed": round(bytes_freed / (1024 * 1024), 2),
            "dry_run": dry_run
        }
    
    def should_skip_stage(self, stage_name: str, input_hash: str) -> bool:
        """Check if a stage should be skipped based on cache."""
        if not self.enabled:
            return False
        
        cached_hash = self._stage_status.get(stage_name, {}).get("hash")
        return cached_hash == input_hash
    
    def mark_stage_complete(self, stage_name: str, input_hash: str):
        """Mark a stage as complete in the cache."""
        self._stage_status[stage_name] = {
            "hash": input_hash,
            "completed_at": json.dumps({"timestamp": str(Path().stat().st_mtime)})
        }
        self._save_stage_status()
    
    def invalidate_stage(self, stage_name: str):
        """Invalidate a cached stage."""
        if stage_name in self._stage_status:
            del self._stage_status[stage_name]
            self._save_stage_status()
    
    def clear_all(self):
        """Clear all cache entries."""
        if self.cache_dir.exists():
            shutil.rmtree(self.cache_dir)
        self.cache_dir.mkdir(exist_ok=True)
        self._stage_status = {}
        self._save_stage_status()


# Global cache manager instance
_cache_manager_instance = None


def get_cache_manager(cache_dir: Path = None, enabled: bool = True) -> CacheManager:
    """
    Get or create a global CacheManager instance.
    
    Args:
        cache_dir: Optional custom cache directory
        enabled: Whether caching is enabled
    
    Returns:
        CacheManager instance
    """
    global _cache_manager_instance
    
    if _cache_manager_instance is None:
        _cache_manager_instance = CacheManager(cache_dir, enabled)
    
    return _cache_manager_instance


def should_skip_stage(stage_name: str, input_hash: str) -> bool:
    """Convenience function to check if a stage should be skipped."""
    return get_cache_manager().should_skip_stage(stage_name, input_hash)


def mark_stage_complete(stage_name: str, input_hash: str):
    """Convenience function to mark a stage as complete."""
    get_cache_manager().mark_stage_complete(stage_name, input_hash)
