#!/usr/bin/env python3
"""
RNASEQ-MINI Python SDK - Programmatic client for API interactions.
Provides a high-level interface for automation and integration.
"""

import json
import asyncio
import aiohttp
import logging
from typing import Dict, List, Optional, Any, Union
from pathlib import Path
from datetime import datetime
import time

logger = logging.getLogger(__name__)


class RNASEQMiniClient:
    """Python SDK client for RNASEQ-MINI API."""

    def __init__(self, base_url: str = "http://localhost:8001", api_key: str = None, timeout: int = 30):
        """
        Initialize RNASEQ-MINI API client.

        Args:
            base_url: Base URL of the API server
            api_key: API key for authentication (if required)
            timeout: Request timeout in seconds
        """
        self.base_url = base_url.rstrip('/')
        self.api_key = api_key
        self.timeout = timeout
        self.session = None

    async def __aenter__(self):
        """Async context manager entry."""
        await self._ensure_session()
        return self

    async def __aexit__(self, exc_type, exc_val, exc_tb):
        """Async context manager exit."""
        await self.close()

    async def _ensure_session(self):
        """Ensure aiohttp session is created."""
        if self.session is None:
            headers = {}
            if self.api_key:
                headers['Authorization'] = f'Bearer {self.api_key}'

            timeout = aiohttp.ClientTimeout(total=self.timeout)
            connector = aiohttp.TCPConnector(limit=100, limit_per_host=30)

            self.session = aiohttp.ClientSession(
                headers=headers,
                timeout=timeout,
                connector=connector
            )

    async def close(self):
        """Close the HTTP session."""
        if self.session:
            await self.session.close()
            self.session = None

    def _get_headers(self) -> Dict[str, str]:
        """Get headers for requests."""
        headers = {'Content-Type': 'application/json'}
        if self.api_key:
            headers['Authorization'] = f'Bearer {self.api_key}'
        return headers

    async def health_check(self) -> Dict[str, Any]:
        """Check API server health."""
        await self._ensure_session()
        async with self.session.get(f"{self.base_url}/health") as response:
            response.raise_for_status()
            return await response.json()

    async def get_api_status(self) -> Dict[str, Any]:
        """Get API server status and capabilities."""
        await self._ensure_session()
        async with self.session.get(f"{self.base_url}/api/v1/status") as response:
            response.raise_for_status()
            return await response.json()

    # Pipeline execution methods
    async def run_pipeline(self, config: Dict[str, Any], wait: bool = False) -> Dict[str, Any]:
        """
        Run RNA-seq pipeline analysis.

        Args:
            config: Pipeline configuration parameters
            wait: Whether to wait for completion

        Returns:
            Job information dictionary
        """
        await self._ensure_session()

        payload = {"config": config, "async": not wait}

        async with self.session.post(
            f"{self.base_url}/api/v1/pipeline/run",
            json=payload,
            headers=self._get_headers()
        ) as response:
            response.raise_for_status()
            job_info = await response.json()

        if wait:
            job_id = job_info["job_id"]
            return await self.wait_for_job_completion(job_id)

        return job_info

    async def get_job_status(self, job_id: str) -> Dict[str, Any]:
        """Get status of a specific job."""
        await self._ensure_session()
        async with self.session.get(f"{self.base_url}/api/v1/pipeline/jobs/{job_id}") as response:
            response.raise_for_status()
            return await response.json()

    async def list_jobs(self, limit: int = 50, offset: int = 0) -> Dict[str, Any]:
        """List recent pipeline jobs."""
        await self._ensure_session()
        params = {'limit': limit, 'offset': offset}
        async with self.session.get(f"{self.base_url}/api/v1/pipeline/jobs", params=params) as response:
            response.raise_for_status()
            return await response.json()

    async def wait_for_job_completion(self, job_id: str, check_interval: int = 5) -> Dict[str, Any]:
        """
        Wait for a job to complete and return results.

        Args:
            job_id: Job ID to wait for
            check_interval: Seconds between status checks

        Returns:
            Final job status and results
        """
        while True:
            job_status = await self.get_job_status(job_id)

            if job_status["status"] in ["completed", "failed"]:
                return job_status

            await asyncio.sleep(check_interval)

    # Configuration methods
    async def validate_config(self, config_file: str = "config/params.yaml") -> Dict[str, Any]:
        """Validate pipeline configuration file."""
        await self._ensure_session()
        params = {'config_file': config_file}
        async with self.session.get(f"{self.base_url}/api/v1/config/validate", params=params) as response:
            response.raise_for_status()
            return await response.json()

    async def optimize_config(self, fastq_files: List[str]) -> Dict[str, Any]:
        """
        Optimize configuration based on FASTQ file characteristics.

        Args:
            fastq_files: List of FASTQ file paths

        Returns:
            Optimized configuration parameters
        """
        await self._ensure_session()

        payload = {"fastq_files": fastq_files}

        async with self.session.post(
            f"{self.base_url}/api/v1/config/optimize",
            json=payload,
            headers=self._get_headers()
        ) as response:
            response.raise_for_status()
            return await response.json()

    # Quality assessment methods
    async def assess_quality(self, results_dir: str = "results") -> Dict[str, Any]:
        """
        Run comprehensive quality assessment.

        Args:
            results_dir: Directory containing analysis results

        Returns:
            Quality assessment report
        """
        await self._ensure_session()

        payload = {"results_dir": results_dir}

        async with self.session.post(
            f"{self.base_url}/api/v1/quality/assess",
            json=payload,
            headers=self._get_headers()
        ) as response:
            response.raise_for_status()
            return await response.json()

    # Multi-omics methods
    async def integrate_multiomics(self, datasets: Dict[str, Any]) -> Dict[str, Any]:
        """
        Integrate multiple omics datasets.

        Args:
            datasets: Dictionary of omics datasets

        Returns:
            Integration results
        """
        await self._ensure_session()

        payload = {"datasets": datasets}

        async with self.session.post(
            f"{self.base_url}/api/v1/multiomics/integrate",
            json=payload,
            headers=self._get_headers()
        ) as response:
            response.raise_for_status()
            return await response.json()

    # Cache management methods
    async def get_cache_stats(self) -> Dict[str, Any]:
        """Get cache statistics."""
        await self._ensure_session()
        async with self.session.get(f"{self.base_url}/api/v1/cache/stats") as response:
            response.raise_for_status()
            return await response.json()

    async def cleanup_cache(self, max_age_days: int = 30, dry_run: bool = False) -> Dict[str, Any]:
        """Clean up cache entries."""
        await self._ensure_session()

        payload = {"max_age_days": max_age_days, "dry_run": dry_run}

        async with self.session.post(
            f"{self.base_url}/api/v1/cache/cleanup",
            json=payload,
            headers=self._get_headers()
        ) as response:
            response.raise_for_status()
            return await response.json()

    # Export methods
    async def export_results(self, format: str = "json", results_dir: str = "results") -> bytes:
        """
        Export analysis results.

        Args:
            format: Export format (json, csv, zip)
            results_dir: Results directory

        Returns:
            File contents as bytes
        """
        await self._ensure_session()

        params = {'results_dir': results_dir}

        async with self.session.get(
            f"{self.base_url}/api/v1/export/results/{format}",
            params=params
        ) as response:
            response.raise_for_status()
            return await response.read()

    # Webhook methods
    async def trigger_webhook(self, event_type: str, data: Dict[str, Any]) -> Dict[str, Any]:
        """Trigger a webhook."""
        await self._ensure_session()

        payload = data

        async with self.session.post(
            f"{self.base_url}/api/v1/webhooks/{event_type}",
            json=payload,
            headers=self._get_headers()
        ) as response:
            response.raise_for_status()
            return await response.json()

    # Plugin methods
    async def list_plugins(self) -> Dict[str, Any]:
        """List available plugins."""
        await self._ensure_session()
        async with self.session.get(f"{self.base_url}/api/v1/plugins") as response:
            response.raise_for_status()
            return await response.json()

    async def execute_plugin(self, plugin_name: str, parameters: Dict[str, Any]) -> Dict[str, Any]:
        """Execute a plugin."""
        await self._ensure_session()

        payload = parameters

        async with self.session.post(
            f"{self.base_url}/api/v1/plugins/{plugin_name}/execute",
            json=payload,
            headers=self._get_headers()
        ) as response:
            response.raise_for_status()
            return await response.json()


# Synchronous wrapper for convenience
class RNASEQMiniClientSync:
    """Synchronous wrapper for RNASEQ-MINI API client."""

    def __init__(self, base_url: str = "http://localhost:8001", api_key: str = None, timeout: int = 30):
        self.client = RNASEQMiniClient(base_url, api_key, timeout)

    def __enter__(self):
        """Sync context manager entry."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Sync context manager exit."""
        # Note: In real implementation, you'd need to handle async cleanup
        pass

    def health_check(self) -> Dict[str, Any]:
        """Synchronous health check."""
        return asyncio.run(self.client.health_check())

    def run_pipeline(self, config: Dict[str, Any], wait: bool = False) -> Dict[str, Any]:
        """Synchronous pipeline execution."""
        return asyncio.run(self.client.run_pipeline(config, wait))

    def get_job_status(self, job_id: str) -> Dict[str, Any]:
        """Synchronous job status check."""
        return asyncio.run(self.client.get_job_status(job_id))

    def list_jobs(self, limit: int = 50, offset: int = 0) -> Dict[str, Any]:
        """Synchronous job listing."""
        return asyncio.run(self.client.list_jobs(limit, offset))

    def wait_for_job_completion(self, job_id: str, check_interval: int = 5) -> Dict[str, Any]:
        """Synchronous job completion waiting."""
        return asyncio.run(self.client.wait_for_job_completion(job_id, check_interval))

    def validate_config(self, config_file: str = "config/params.yaml") -> Dict[str, Any]:
        """Synchronous config validation."""
        return asyncio.run(self.client.validate_config(config_file))

    def optimize_config(self, fastq_files: List[str]) -> Dict[str, Any]:
        """Synchronous config optimization."""
        return asyncio.run(self.client.optimize_config(fastq_files))

    def assess_quality(self, results_dir: str = "results") -> Dict[str, Any]:
        """Synchronous quality assessment."""
        return asyncio.run(self.client.assess_quality(results_dir))

    def integrate_multiomics(self, datasets: Dict[str, Any]) -> Dict[str, Any]:
        """Synchronous multi-omics integration."""
        return asyncio.run(self.client.integrate_multiomics(datasets))

    def get_cache_stats(self) -> Dict[str, Any]:
        """Synchronous cache statistics."""
        return asyncio.run(self.client.get_cache_stats())

    def cleanup_cache(self, max_age_days: int = 30, dry_run: bool = False) -> Dict[str, Any]:
        """Synchronous cache cleanup."""
        return asyncio.run(self.client.cleanup_cache(max_age_days, dry_run))

    def export_results(self, format: str = "json", results_dir: str = "results") -> bytes:
        """Synchronous results export."""
        return asyncio.run(self.client.export_results(format, results_dir))

    def trigger_webhook(self, event_type: str, data: Dict[str, Any]) -> Dict[str, Any]:
        """Synchronous webhook trigger."""
        return asyncio.run(self.client.trigger_webhook(event_type, data))

    def list_plugins(self) -> Dict[str, Any]:
        """Synchronous plugin listing."""
        return asyncio.run(self.client.list_plugins())

    def execute_plugin(self, plugin_name: str, parameters: Dict[str, Any]) -> Dict[str, Any]:
        """Synchronous plugin execution."""
        return asyncio.run(self.client.execute_plugin(plugin_name, parameters))


# Utility functions for common operations
class RNASEQMiniSDK:
    """High-level SDK with convenience methods."""

    def __init__(self, base_url: str = "http://localhost:8001", api_key: str = None):
        self.client = RNASEQMiniClientSync(base_url, api_key)

    def run_complete_analysis(self, fastq_files: List[str], config: Dict[str, Any] = None) -> Dict[str, Any]:
        """
        Run complete RNA-seq analysis from FASTQ files to results.

        Args:
            fastq_files: List of FASTQ file paths
            config: Optional configuration overrides

        Returns:
            Complete analysis results
        """
        # Optimize configuration based on FASTQ files
        if not config:
            optimization = self.client.optimize_config(fastq_files)
            config = optimization["optimized_config"]

        # Run pipeline
        job_result = self.client.run_pipeline(config, wait=True)

        if job_result["status"] == "completed":
            # Run quality assessment
            quality_report = self.client.assess_quality()

            return {
                "job": job_result,
                "quality": quality_report,
                "config_used": config,
                "fastq_files": fastq_files
            }
        else:
            return {"job": job_result, "error": "Pipeline execution failed"}

    def monitor_analysis_progress(self, job_id: str) -> Dict[str, Any]:
        """Monitor progress of a running analysis."""
        return self.client.get_job_status(job_id)

    def export_analysis_results(self, job_id: str, format: str = "json") -> bytes:
        """Export results from a completed analysis."""
        return self.client.export_results(format)

    def setup_webhooks(self, webhook_config: Dict[str, str]) -> bool:
        """
        Setup webhooks for analysis events.

        Args:
            webhook_config: Dictionary mapping event types to URLs

        Returns:
            Success status
        """
        for event_type, url in webhook_config.items():
            try:
                self.client.trigger_webhook(event_type, {"url": url, "setup": True})
            except Exception as e:
                logger.error(f"Failed to setup webhook for {event_type}: {e}")
                return False

        return True

    def get_system_status(self) -> Dict[str, Any]:
        """Get comprehensive system status."""
        health = self.client.health_check()
        api_status = self.client.get_api_status()
        cache_stats = self.client.get_cache_stats()

        return {
            "health": health,
            "api": api_status,
            "cache": cache_stats,
            "timestamp": datetime.now().isoformat()
        }


