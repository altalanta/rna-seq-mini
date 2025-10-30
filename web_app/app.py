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
import uvicorn

# Import the dashboard creation function
from web_app.dashboard import create_dashboard

app = FastAPI(title="RNASEQ-MINI Web Interface")

# Create and mount the Dash dashboard
dash_app = create_dashboard(app)
app.mount("/interactive", WSGIMiddleware(dash_app.server))

# Mount static files
static_dir = Path(__file__).parent / "static"
app.mount("/static", StaticFiles(directory=str(static_dir)), name="static")

# Mount templates (for serving HTML files)
templates_dir = Path(__file__).parent / "templates"

@app.get("/", response_class=HTMLResponse)
async def dashboard():
    """Serve the main dashboard."""
    dashboard_file = templates_dir / "dashboard.html"
    if dashboard_file.exists():
        with open(dashboard_file, 'r') as f:
            content = f.read()
        return HTMLResponse(content)
    else:
        return HTMLResponse("<h1>Dashboard not found. Please run 'make serve-results'</h1>")

@app.get("/tutorial", response_class=HTMLResponse)
async def tutorial():
    """Serve the interactive tutorial."""
    tutorial_file = templates_dir / "tutorial.html"
    if tutorial_file.exists():
        with open(tutorial_file, 'r') as f:
            content = f.read()
        return HTMLResponse(content)
    else:
        return HTMLResponse("<h1>Tutorial not found</h1>")

@app.get("/dashboard", response_class=HTMLResponse)
async def dashboard_redirect():
    """Redirect to main dashboard."""
    return await dashboard()

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
