from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict, List, Tuple


def read_samples(path: str) -> List[Dict[str, str]]:
    sample_path = Path(path)
    if not sample_path.exists():
        raise FileNotFoundError(f"Sample sheet not found: {sample_path}")
    with sample_path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = [row for row in reader]
    for row in rows:
        if "sample" not in row or "condition" not in row:
            raise ValueError("Sample sheet must contain 'sample' and 'condition' columns")
    return rows


def read_contrasts(path: str) -> List[Tuple[str, str]]:
    contrast_path = Path(path)
    if not contrast_path.exists():
        raise FileNotFoundError(f"Contrasts file not found: {contrast_path}")
    with contrast_path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "groupA" not in reader.fieldnames or "groupB" not in reader.fieldnames:
            raise ValueError("Contrast file must contain 'groupA' and 'groupB'")
        return [(row["groupA"], row["groupB"]) for row in reader]


def ensure_directory(path: str) -> None:
    Path(path).mkdir(parents=True, exist_ok=True)


def sample_is_single_end(config: Dict) -> bool:
    return bool(config.get("se", False))


def salmon_quant_dir(config: Dict) -> Path:
    return Path(config["paths"]["salmon"]) / "{sample}"


def results_dir(config: Dict, key: str) -> Path:
    return Path(config["paths"][key])
