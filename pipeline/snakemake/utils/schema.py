from __future__ import annotations

from pathlib import Path
from typing import Any, Dict


REQUIRED_KEYS = [
    "engine",
    "paths",
    "reference",
    "salmon",
    "r",
]


def validate_config(config: Dict[str, Any]) -> None:
    missing = [key for key in REQUIRED_KEYS if key not in config]
    if missing:
        raise ValueError(f"Missing required config keys: {', '.join(missing)}")
    for required_path in [
        config["paths"]["samples"],
        config["reference"]["transcripts_fa"],
        config["reference"]["annotation_gtf"],
    ]:
        path = Path(required_path)
        if not path.exists():
            raise FileNotFoundError(f"Expected file not found: {path}")
