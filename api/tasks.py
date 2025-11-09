import os
import subprocess
import traceback
from pathlib import Path
import yaml
from datetime import datetime

# Add src to the Python path to allow importing modules
import sys
sys.path.append(str(Path(__file__).resolve().parent.parent / "src"))

from rnaseq_mini.logger import get_logger
from .celery_app import celery_app
from .db import SessionLocal, Job

log = get_logger("worker")

@celery_app.task(bind=True)
def execute_pipeline_task(self, job_id: str):
    """Celery task to execute a pipeline run."""
    db = SessionLocal()
    job = db.query(Job).filter(Job.id == job_id).first()
    if not job:
        log.error("job_not_found_in_db", job_id=job_id)
        return

    try:
        log.info("pipeline_task_started", job_id=job_id)
        job.status = "running"
        job.started_at = datetime.utcnow()
        db.commit()

        results_dir = Path("results") / f"job_{job_id}"
        results_dir.mkdir(parents=True, exist_ok=True)
        job.results_dir = str(results_dir)

        params_file = results_dir / "params.yaml"
        with open("config/params.yaml") as f:
            run_config = yaml.safe_load(f)
        
        # In a real app, job.parameters would be a JSON string
        import json
        request_params = json.loads(job.parameters)
        run_config.update(request_params.get("config", {}))
        run_config["results_dir"] = str(results_dir)
        
        with open(params_file, 'w') as f:
            yaml.dump(run_config, f)

        engine = run_config.get("engine", "snakemake")
        log_file = results_dir / "pipeline.log"
        job.log_file = str(log_file)

        cmd = []
        if engine == "snakemake":
            cmd = ["snakemake", "-s", "pipeline/snakemake/Snakefile", "--configfile", str(params_file), "--use-conda", "--cores", str(run_config.get("threads", 2))]
        elif engine == "nextflow":
            cmd = ["nextflow", "run", "pipeline/nextflow/main.nf", "-params-file", str(params_file), "-with-conda", "-profile", "local", "--results_dir", str(results_dir)]

        with open(log_file, "w") as log_f:
            process = subprocess.Popen(cmd, stdout=log_f, stderr=subprocess.STDOUT)
        
        process.wait()

        if process.returncode == 0:
            job.status = "completed"
            job.completed_at = datetime.utcnow()
            log.info("pipeline_task_succeeded", job_id=job_id)
        else:
            raise RuntimeError(f"Pipeline exited with non-zero status: {process.returncode}")

    except Exception as e:
        error_message = f"{e}\n\n{traceback.format_exc()}"
        log.error("pipeline_task_failed", job_id=job_id, error=str(e), traceback=traceback.format_exc())
        job.status = "failed"
        job.error = error_message
    finally:
        db.commit()
        db.close()
