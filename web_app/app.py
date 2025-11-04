#!/usr/bin/env python3
"""
Simple web server for RNASEQ-MINI dashboard and tutorials.
Serves static files and templates for the interactive web interface.
"""

import os
from pathlib import Path
from fastapi import FastAPI, Request
from fastapi.responses import HTMLResponse, FileResponse
from fastapi.staticfiles import StaticFiles
from fastapi.middleware.wsgi import WSGIMiddleware
from fastapi.templating import Jinja2Templates
import uvicorn
from api.server import RNASEQMiniAPI

# Import the dashboard creation function
from web_app.dashboard import create_dashboard

app = FastAPI(title="RNASEQ-MINI Web Interface")

# --- API Integration ---
# In a real app, you might run the API server in a separate process,
# but for simplicity, we'll instantiate it here to access job data.
api = RNASEQMiniAPI()
app.include_router(api.app.router, prefix="/api")


# --- Dashboard Mounting ---
# The main dashboard will show the latest default results
dash_app = create_dashboard(app)
app.mount("/interactive", WSGIMiddleware(dash_app.server))

# --- Static HTML and File Serving ---
static_dir = Path(__file__).parent / "static"
templates_dir = Path(__file__).parent / "templates"
app.mount("/static", StaticFiles(directory=str(static_dir)), name="static")
templates = Jinja2Templates(directory=str(templates_dir))


@app.get("/", response_class=HTMLResponse)
async def root(request: Request):
    """Serve the main landing page which lists available jobs."""
    jobs = list(api.analysis_jobs.values())
    return templates.TemplateResponse("index.html", {"request": request, "jobs": jobs})

@app.get("/dashboard/{job_id}", response_class=HTMLResponse)
async def job_dashboard(request: Request, job_id: str):
    """Serve an interactive dashboard for a specific job ID."""
    if job_id not in api.analysis_jobs:
        return HTMLResponse("Job not found.", status_code=404)
    
    results_path = api.analysis_jobs[job_id].get("results_dir", "results")
    
    # Create a new Dash app instance for this specific job
    dash_app_job = create_dashboard(app, results_dir_path=results_path)
    
    # Mount it on a unique path
    mount_path = f"/interactive/{job_id}"
    app.mount(mount_path, WSGIMiddleware(dash_app_job.server))

    # Redirect the user to the interactive dashboard
    # In a real application, you might render a template that includes the Dash app iframe
    from fastapi.responses import RedirectResponse
    return RedirectResponse(url=mount_path)


@app.get("/tutorial", response_class=HTMLResponse)
async def tutorial(request: Request):
    """Serve the interactive tutorial."""
    return templates.TemplateResponse("tutorial.html", {"request": request})

@app.get("/dashboard", response_class=HTMLResponse)
async def dashboard_redirect():
    """Redirect to main dashboard."""
    return await root(Request(scope={"type": "http", "path": "/", "query_string": b""}))

@app.get("/api/results")
async def api_results():
    """Mock API endpoint for results."""
    return {
        "status": "success",
        "message": "API endpoint for RNA-seq results",
        "data": {
            "samples": [],
            "contrasts": [],
            "results": {}
        }
    }

if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=8000)
