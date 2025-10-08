#!/usr/bin/env python3
"""Minimal Slurm status checker used for Snakemake cluster mode.
This implementation is intentionally lightweight for smoke tests.
"""
from __future__ import annotations

import subprocess
import sys


def main() -> int:
    job_id = sys.argv[1] if len(sys.argv) > 1 else None
    if not job_id:
        return 0
    try:
        cmd = ["squeue", "-j", str(job_id), "-h", "-o", "%T"]
        status = subprocess.check_output(cmd, text=True).strip().lower()
    except Exception:
        status = "unknown"
    mapping = {
        "pending": "running",
        "configuring": "running",
        "running": "running",
        "completing": "running",
        "completed": "success",
        "failed": "failed",
        "timeout": "failed",
        "cancelled": "failed",
        "preempted": "failed",
        "unknown": "running",
    }
    print(mapping.get(status, "running"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
