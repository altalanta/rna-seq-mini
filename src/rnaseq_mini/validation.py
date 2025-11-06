#!/usr/bin/env python3
"""
Config validation script for RNASEQ-MINI.
Validates config/params.yaml against a formal JSON schema.
"""

import yaml
import sys
import json
import argparse
from pathlib import Path
from jsonschema import validate, exceptions
from typing import List, Optional

from pydantic import BaseModel, Field, FilePath, DirectoryPath, field_validator
import pandas as pd

# --- Helper Validators ---

class SamplesTsv(BaseModel):
    sample: str
    fastq_1: FilePath
    fastq_2: Optional[FilePath] = None
    condition: str
    
    # Allow extra columns like 'batch'
    class Config:
        extra = 'allow'

class ContrastsTsv(BaseModel):
    contrast_a: str
    contrast_b: str

# --- Main Configuration Models (params.yaml) ---

class ReferenceConfig(BaseModel):
    transcripts_fa: FilePath
    annotation_gtf: FilePath
    decoy_fasta: Optional[FilePath] = None
    salmon_index: Path

class SalmonConfig(BaseModel):
    libtype: str = "A"
    extra: str = "--validateMappings --gcBias"
    threads: int = 4

class RConfig(BaseModel):
    design: str = "~ condition"
    contrast_variable: str = "condition"
    contrasts_file: FilePath

class PathsConfig(BaseModel):
    samples: FilePath
    outdir: Path = "results"
    logs: Path = "logs"

class SingleCellAnalysisConfig(BaseModel):
    min_genes_per_cell: int = 200
    min_cells_per_gene: int = 3
    max_mito_percent: float = 15.0
    n_pcs: int = 30
    n_neighbors: int = 15
    clustering_resolution: float = 0.5

class SingleCellPathsConfig(BaseModel):
    output_dir: Path = "results/singlecell/analysis"

class SingleCellConfig(BaseModel):
    enabled: bool = False
    input_dir: Optional[DirectoryPath] = None
    analysis: SingleCellAnalysisConfig = Field(default_factory=SingleCellAnalysisConfig)
    paths: SingleCellPathsConfig = Field(default_factory=SingleCellPathsConfig)

    @field_validator('input_dir')
    def check_input_dir_if_enabled(cls, v, values):
        data = values.data
        if data.get('enabled') and v is None:
            raise ValueError("`input_dir` must be specified when single-cell analysis is enabled.")
        return v

class FullConfig(BaseModel):
    engine: str = "snakemake"
    threads: int
    organism_name: str
    reference: ReferenceConfig
    salmon: SalmonConfig
    r: RConfig
    paths: PathsConfig
    singlecell: SingleCellConfig = Field(default_factory=SingleCellConfig)

    @field_validator('engine')
    def engine_must_be_valid(cls, v):
        if v not in ["snakemake", "nextflow"]:
            raise ValueError("Engine must be 'snakemake' or 'nextflow'")
        return v

# --- Main Validation Function ---

def load_and_validate_config(config_path: Path) -> FullConfig:
    """Loads, parses, and validates all configuration files."""
    import yaml
    
    # 1. Load and validate params.yaml
    with open(config_path) as f:
        config_dict = yaml.safe_load(f)
    config = FullConfig.parse_obj(config_dict)

    # 2. Load and validate samples.tsv
    samples_df = pd.read_csv(config.paths.samples, sep='	')
    samples_data = [SamplesTsv.parse_obj(row) for _, row in samples_df.iterrows()]

    # 3. Load and validate contrasts.tsv
    contrasts_df = pd.read_csv(config.r.contrasts_file, sep='	', header=None, names=['contrast_a', 'contrast_b'])
    contrasts_data = [ContrastsTsv.parse_obj(row) for _, row in contrasts_df.iterrows()]
    
    # 4. Perform cross-validation
    sample_conditions = {row.condition for row in samples_data}
    for contrast in contrasts_data:
        if contrast.contrast_a not in sample_conditions:
            raise ValueError(f"Contrast value '{contrast.contrast_a}' not found in samples sheet condition column.")
        if contrast.contrast_b not in sample_conditions:
            raise ValueError(f"Contrast value '{contrast.contrast_b}' not found in samples sheet condition column.")
            
    if config.r.contrast_variable not in samples_df.columns:
        raise ValueError(f"Contrast variable '{config.r.contrast_variable}' not found in samples sheet columns.")

    return config

