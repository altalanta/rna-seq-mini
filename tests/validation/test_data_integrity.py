import pytest
import pandas as pd
from pathlib import Path
import numpy as np

# Define the expected location of smoke test results
# In the CI workflow, we run the smoke test into a 'results-snakemake' dir
SMOKE_RESULTS_DIR = Path.cwd() / "results-snakemake"

@pytest.mark.integration
def test_smoke_results_exist():
    """Check that essential smoke test output files were created."""
    assert SMOKE_RESULTS_DIR.is_dir(), "Smoke test results directory not found."
    
    expected_files = [
        "report.html",
        "counts/counts.tsv",
        "counts/tpm.tsv",
        "de/de_summary.tsv",
        "qc/multiqc/multiqc_report.html"
    ]
    
    for f in expected_files:
        file_path = SMOKE_RESULTS_DIR / f
        assert file_path.is_file(), f"Expected output file not found: {file_path}"

@pytest.mark.integration
def test_counts_table_integrity():
    """Validate the structure and content of the gene counts table."""
    counts_file = SMOKE_RESULTS_DIR / "counts/counts.tsv"
    assert counts_file.is_file()
    
    df = pd.read_csv(counts_file, sep='\t')
    
    # Check structure
    assert not df.empty, "Counts table is empty."
    expected_cols = ['gene', 'A_1', 'B_1'] # Based on smoke test data
    for col in expected_cols:
        assert col in df.columns, f"Expected column '{col}' not in counts table."
        
    # Check content
    numeric_cols = df.select_dtypes(include=np.number).columns
    assert (df[numeric_cols] >= 0).all().all(), "Found negative values in counts table."
    assert df[numeric_cols].apply(lambda s: pd.to_numeric(s, errors='coerce').notnull().all()).all(), "Found non-integer values in counts table."


@pytest.mark.integration
def test_tpm_table_normalization():
    """Validate that TPM values are correctly normalized (sum to 1e6)."""
    tpm_file = SMOKE_RESULTS_DIR / "counts/tpm.tsv"
    assert tpm_file.is_file()
    
    df = pd.read_csv(tpm_file, sep='\t')
    
    sample_cols = [col for col in df.columns if col != 'gene']
    
    for col in sample_cols:
        # Check that the sum of TPM values is close to 1 million
        assert np.isclose(df[col].sum(), 1e6, rtol=1e-3), f"TPM column '{col}' does not sum to 1e6."
