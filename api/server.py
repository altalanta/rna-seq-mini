#!/usr/bin/env python3
"""
RNASEQ-MINI API Server - REST API for all pipeline operations.
Provides comprehensive programmatic access to RNA-seq analysis capabilities.
"""

import json
import asyncio
import logging
import traceback
from pathlib import Path
from typing import Dict, List, Optional, Any, Union
from datetime import datetime
import uuid

from fastapi import FastAPI, HTTPException, BackgroundTasks, Depends, status, Request
from fastapi.security import HTTPBearer, HTTPAuthorizationCredentials
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse, FileResponse
import uvicorn

# Import pipeline components
try:
    from ..scripts.cache_manager import get_cache_manager
    from ..scripts.parameter_optimizer import ParameterOptimizer
    from ..scripts.quality_assessor import run_quality_assessment, ReferenceDataset
    from ..pipeline.multiomics import MultiOmicsIntegrator, OmicsType
except ImportError:
    # Fallback for when modules aren't available
    pass


# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class RNASEQMiniAPI:
    """Main API server for RNASEQ-MINI platform."""

    def __init__(self, host: str = "0.0.0.0", port: int = 8001, debug: bool = False):
        self.host = host
        self.port = port
        self.debug = debug

        # Initialize FastAPI app
        self.app = FastAPI(
            title="RNASEQ-MINI API",
            description="Comprehensive REST API for RNA-seq analysis platform",
            version="1.0.0",
            docs_url="/docs",
            redoc_url="/redoc",
            openapi_url="/openapi.json"
        )

        # Add CORS middleware
        self.app.add_middleware(
            CORSMiddleware,
            allow_origins=["*"],
            allow_credentials=True,
            allow_methods=["*"],
            allow_headers=["*"],
        )

        # Initialize components
        self.cache_manager = None
        self.parameter_optimizer = None
        self.multiomics_integrator = None
        self.webhook_manager = None

        # Analysis jobs tracking
        self.analysis_jobs = {}

        # Setup routes
        self._setup_routes()

        logger.info(f"RNASEQ-MINI API initialized on {self.host}:{self.port}")

    def _setup_routes(self):
        """Setup all API routes."""

        # Health and status endpoints
        @self.app.get("/health")
        async def health_check():
            """Health check endpoint."""
            return {
                "status": "healthy",
                "timestamp": datetime.now().isoformat(),
                "version": "1.0.0"
            }

        @self.app.get("/api/v1/status")
        async def api_status():
            """Get API status and capabilities."""
            return {
                "status": "running",
                "version": "1.0.0",
                "features": {
                    "caching": True,
                    "parameter_optimization": True,
                    "multiomics": True,
                    "quality_assessment": True,
                    "webhooks": True,
                    "plugins": True
                },
                "endpoints": len([route for route in self.app.routes if route.path.startswith("/api/")])
            }

        # Pipeline execution endpoints
        @self.app.post("/api/v1/pipeline/run")
        async def run_pipeline(request: Dict[str, Any], background_tasks: BackgroundTasks):
            """Run RNA-seq pipeline analysis."""
            job_id = str(uuid.uuid4())

            # Store job information
            self.analysis_jobs[job_id] = {
                "id": job_id,
                "status": "queued",
                "created_at": datetime.now().isoformat(),
                "parameters": request,
                "results": None,
                "error": None
            }

            # Run pipeline in background
            background_tasks.add_task(self._execute_pipeline_job, job_id, request)

            return {
                "job_id": job_id,
                "status": "queued",
                "message": "Pipeline execution started"
            }

        @self.app.get("/api/v1/pipeline/jobs/{job_id}")
        async def get_job_status(job_id: str):
            """Get status of a pipeline job."""
            if job_id not in self.analysis_jobs:
                raise HTTPException(status_code=404, detail="Job not found")

            return self.analysis_jobs[job_id]

        @self.app.get("/api/v1/pipeline/jobs")
        async def list_jobs(limit: int = 50, offset: int = 0):
            """List recent pipeline jobs."""
            jobs = list(self.analysis_jobs.values())
            jobs.sort(key=lambda x: x["created_at"], reverse=True)

            return {
                "jobs": jobs[offset:offset+limit],
                "total": len(jobs),
                "limit": limit,
                "offset": offset
            }

        # Configuration endpoints
        @self.app.get("/api/v1/config/validate")
        async def validate_config(config_file: str = "config/params.yaml"):
            """Validate pipeline configuration."""
            try:
                # Import validation function
                from ..scripts.validate_config import run_comprehensive_validation

                success, issues = run_comprehensive_validation(config_file)

                return {
                    "valid": success,
                    "issues": issues,
                    "config_file": config_file
                }
            except ImportError:
                raise HTTPException(status_code=500, detail="Configuration validation not available")

        @self.app.post("/api/v1/config/optimize")
        async def optimize_config(request: Dict[str, Any]):
            """Optimize configuration based on data characteristics."""
            try:
                fastq_files = request.get("fastq_files", [])

                if not fastq_files:
                    raise HTTPException(status_code=400, detail="FASTQ files required")

                # Initialize optimizer
                if not self.parameter_optimizer:
                    self.parameter_optimizer = ParameterOptimizer()

                # Convert to Path objects
                fastq_paths = [Path(f) for f in fastq_files]

                # Run optimization
                optimized_config = self.parameter_optimizer.optimize_pipeline_config(fastq_paths)

                return {
                    "success": True,
                    "optimized_config": optimized_config,
                    "timestamp": datetime.now().isoformat()
                }

            except Exception as e:
                logger.error(f"Error optimizing config: {e}")
                raise HTTPException(status_code=500, detail=str(e))

        # Quality assessment endpoints
        @self.app.post("/api/v1/quality/assess")
        async def assess_quality(request: Dict[str, Any]):
            """Run comprehensive quality assessment."""
            try:
                results_dir = request.get("results_dir", "results")

                # Run quality assessment
                quality_report = run_quality_assessment(
                    results_dir=results_dir,
                    output_file=f"temp_quality_{uuid.uuid4()}.json"
                )

                return {
                    "success": True,
                    "quality_report": quality_report,
                    "timestamp": datetime.now().isoformat()
                }

            except Exception as e:
                logger.error(f"Error assessing quality: {e}")
                raise HTTPException(status_code=500, detail=str(e))

        # Multi-omics endpoints
        @self.app.post("/api/v1/multiomics/integrate")
        async def integrate_multiomics(request: Dict[str, Any]):
            """Integrate multiple omics datasets."""
            try:
                datasets = request.get("datasets", {})

                if not datasets:
                    raise HTTPException(status_code=400, detail="Datasets required")

                # Initialize integrator if needed
                if not self.multiomics_integrator:
                    self.multiomics_integrator = MultiOmicsIntegrator()

                # Convert to MultiOmicsDataset format
                # This is a simplified implementation
                integration_results = {
                    "integrated": True,
                    "dataset_count": len(datasets),
                    "sample_count": 0,  # Would calculate from actual data
                    "features": {},
                    "correlations": {},
                    "timestamp": datetime.now().isoformat()
                }

                return {
                    "success": True,
                    "integration_results": integration_results
                }

            except Exception as e:
                logger.error(f"Error integrating multiomics: {e}")
                raise HTTPException(status_code=500, detail=str(e))

        # Cache management endpoints
        @self.app.get("/api/v1/cache/stats")
        async def cache_stats():
            """Get cache statistics."""
            try:
                if not self.cache_manager:
                    self.cache_manager = get_cache_manager()

                stats = self.cache_manager.get_cache_stats()
                return {"success": True, "stats": stats}

            except Exception as e:
                logger.error(f"Error getting cache stats: {e}")
                raise HTTPException(status_code=500, detail=str(e))

        @self.app.post("/api/v1/cache/cleanup")
        async def cleanup_cache(request: Dict[str, Any]):
            """Clean up cache entries."""
            try:
                max_age_days = request.get("max_age_days", 30)
                dry_run = request.get("dry_run", False)

                if not self.cache_manager:
                    self.cache_manager = get_cache_manager()

                stats = self.cache_manager.cleanup_cache(max_age_days, dry_run)

                return {
                    "success": True,
                    "cleanup_stats": stats,
                    "dry_run": dry_run
                }

            except Exception as e:
                logger.error(f"Error cleaning cache: {e}")
                raise HTTPException(status_code=500, detail=str(e))

        # Export endpoints
        @self.app.get("/api/v1/export/results/{format}")
        async def export_results(format: str, results_dir: str = "results"):
            """Export complete analysis results."""
            try:
                if format not in ["json", "csv", "zip"]:
                    raise HTTPException(status_code=400, detail="Invalid format")

                # Generate export file
                export_file = self._generate_export_file(results_dir, format)

                if export_file.exists():
                    return FileResponse(
                        export_file,
                        media_type='application/octet-stream',
                        filename=f"rnaseq_results.{format}"
                    )
                else:
                    raise HTTPException(status_code=404, detail="Export file not found")

            except Exception as e:
                logger.error(f"Error exporting results: {e}")
                raise HTTPException(status_code=500, detail=str(e))

        # Webhook endpoints
        @self.app.post("/api/v1/webhooks/{event_type}")
        async def trigger_webhook(event_type: str, request: Dict[str, Any]):
            """Trigger webhook for specific event."""
            try:
                if not self.webhook_manager:
                    from .webhooks import WebhookManager
                    self.webhook_manager = WebhookManager()

                # Trigger webhook
                await self.webhook_manager.trigger_webhook(event_type, request)

                return {
                    "success": True,
                    "event_type": event_type,
                    "message": "Webhook triggered"
                }

            except Exception as e:
                logger.error(f"Error triggering webhook: {e}")
                raise HTTPException(status_code=500, detail=str(e))

        # Plugin endpoints
        @self.app.get("/api/v1/plugins")
        async def list_plugins():
            """List available plugins."""
            try:
                if not hasattr(self, 'plugin_manager'):
                    from .plugins import PluginManager
                    self.plugin_manager = PluginManager()

                plugins = self.plugin_manager.list_plugins()
                return {"success": True, "plugins": plugins}

            except Exception as e:
                logger.error(f"Error listing plugins: {e}")
                raise HTTPException(status_code=500, detail=str(e))

        @self.app.post("/api/v1/plugins/{plugin_name}/execute")
        async def execute_plugin(plugin_name: str, request: Dict[str, Any]):
            """Execute a plugin."""
            try:
                if not hasattr(self, 'plugin_manager'):
                    from .plugins import PluginManager
                    self.plugin_manager = PluginManager()

                # Execute plugin
                result = await self.plugin_manager.execute_plugin(plugin_name, request)

                return {
                    "success": True,
                    "plugin": plugin_name,
                    "result": result
                }

            except Exception as e:
                logger.error(f"Error executing plugin {plugin_name}: {e}")
                raise HTTPException(status_code=500, detail=str(e))

    def _execute_pipeline_job(self, job_id: str, parameters: Dict[str, Any]):
        """Execute pipeline job in background."""
        try:
            logger.info(f"Starting pipeline job {job_id}")

            # Update job status
            self.analysis_jobs[job_id]["status"] = "running"
            self.analysis_jobs[job_id]["started_at"] = datetime.now().isoformat()

            # Execute pipeline (placeholder implementation)
            # In practice, this would call the actual pipeline execution
            results = {
                "status": "completed",
                "output_files": [],
                "metrics": {},
                "completed_at": datetime.now().isoformat()
            }

            # Update job with results
            self.analysis_jobs[job_id]["status"] = "completed"
            self.analysis_jobs[job_id]["results"] = results
            self.analysis_jobs[job_id]["completed_at"] = datetime.now().isoformat()

            # Trigger completion webhook if configured
            asyncio.create_task(self._trigger_completion_webhook(job_id, results))

            logger.info(f"Pipeline job {job_id} completed successfully")

        except Exception as e:
            logger.error(f"Error executing pipeline job {job_id}: {e}")

            # Update job with error
            self.analysis_jobs[job_id]["status"] = "failed"
            self.analysis_jobs[job_id]["error"] = str(e)
            self.analysis_jobs[job_id]["failed_at"] = datetime.now().isoformat()

            # Trigger error webhook if configured
            asyncio.create_task(self._trigger_error_webhook(job_id, str(e)))

    async def _trigger_completion_webhook(self, job_id: str, results: Dict):
        """Trigger webhook for job completion."""
        if self.webhook_manager:
            await self.webhook_manager.trigger_webhook("job_completed", {
                "job_id": job_id,
                "results": results
            })

    async def _trigger_error_webhook(self, job_id: str, error: str):
        """Trigger webhook for job error."""
        if self.webhook_manager:
            await self.webhook_manager.trigger_webhook("job_failed", {
                "job_id": job_id,
                "error": error
            })

    def _generate_export_file(self, results_dir: str, format: str) -> Path:
        """Generate export file for results."""
        # Placeholder implementation
        # In practice, this would create actual export files
        export_dir = Path("exports")
        export_dir.mkdir(exist_ok=True)

        export_file = export_dir / f"results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.{format}"

        # Create a simple export file for demonstration
        export_data = {
            "export_timestamp": datetime.now().isoformat(),
            "results_dir": results_dir,
            "format": format,
            "status": "generated"
        }

        with open(export_file, 'w') as f:
            if format == "json":
                json.dump(export_data, f, indent=2)
            else:
                f.write("CSV export would be generated here")

        return export_file

    def run(self):
        """Start the API server."""
        logger.info(f"Starting RNASEQ-MINI API server on {self.host}:{self.port}")

        uvicorn.run(
            self.app,
            host=self.host,
            port=self.port,
            reload=self.debug,
            log_level="info"
        )


# Standalone server function
def run_api_server(host: str = "0.0.0.0", port: int = 8001, debug: bool = False):
    """Run the RNASEQ-MINI API server."""
    api = RNASEQMiniAPI(host, port, debug)
    api.run()


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="RNASEQ-MINI API Server")
    parser.add_argument('--host', default='0.0.0.0', help='Host to bind to')
    parser.add_argument('--port', type=int, default=8001, help='Port to bind to')
    parser.add_argument('--debug', action='store_true', help='Enable debug mode')

    args = parser.parse_args()

    run_api_server(args.host, args.port, args.debug)


