from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict, List, Tuple
import yaml


def load_config(config_path: str) -> dict:
    """
    Loads a YAML configuration file that may contain 'include' directives.
    """
    config_path = Path(config_path)
    with config_path.open() as f:
        main_config = yaml.safe_load(f)

    if "include" not in main_config:
        return main_config

    base_dir = config_path.parent
    merged_config = {}
    for include_file in main_config["include"]:
        include_path = base_dir / include_file
        with include_path.open() as f:
            included_config = yaml.safe_load(f)
            merged_config.update(included_config)
    
    # Add any top-level keys from the main config file (that are not 'include')
    for key, value in main_config.items():
        if key != "include":
            merged_config[key] = value
            
    return merged_config


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


def get_yaml_loader():
    """Returns a YAML loader that supports the !include directive."""
    try:
        from yaml.loader import SafeLoader
        import yaml
    except ImportError:
        raise ImportError("PyYAML is not installed. Please install it to continue.")

    class YamlIncludeLoader(SafeLoader):
        def __init__(self, stream):
            self._root = Path(stream.name).parent
            super().__init__(stream)

        def include(self, node):
            filename = self._root / self.construct_scalar(node)
            with open(filename, 'r') as f:
                return yaml.load(f, YamlIncludeLoader)

    YamlIncludeLoader.add_constructor("!include", YamlIncludeLoader.include)
    return YamlIncludeLoader
