# Repository Guidelines

## Project Overview (For Non-Technical Contributors)
- This repository automates RNA sequencing analysis end-to-end: raw FASTQ reads go through quality checks, expression quantification, statistical comparisons, and pathway summaries, finishing with an HTML report.
- Two workflow engines are supported—Snakemake and Nextflow—so the same analysis can run on a laptop or an HPC cluster without code changes.
- Key moving parts:
  - *Quality control*: FastQC and MultiQC inspect read quality.
  - *Quantification*: Salmon estimates transcript abundance and gene counts.
  - *Differential expression & pathways*: R scripts (DESeq2 + fgsea) compare conditions and highlight enriched pathways.
  - *Reporting*: A templated R Markdown document compiles plots, tables, and session info.
- To explore results, open `results/report.html` after any run; it links out to QC dashboards, DE tables, and pathway findings.
- New contributors can start by reviewing `config/params.yaml` (pipeline settings) and the sample fixtures under `tests/` to see minimal working inputs before tackling larger datasets.

## Project Structure & Module Organization
- **Config & presets:** `config/params.yaml`, organism presets in `config/genome.yaml`, Snakemake/Nextflow profiles under `config/profiles/`.
- **Workflows:** Snakemake rules live in `pipeline/snakemake/` (modular `rules/*.smk`, helpers in `utils/`); Nextflow DSL2 modules are under `pipeline/nextflow/modules/` with the entrypoint `main.nf`.
- **Analysis scripts:** R utilities (`tximport_deseq2.R`, `fgsea_pathways.R`, `render_report.R`) live in `scripts/`; shell helpers (e.g., `build_salmon_index.sh`) share the same directory.
- **Assets & tests:** Yeast fixtures in `references/yeast/`; synthetic FASTQs under `tests/data/fastq/`; smoke harness in `tests/run_smoke.sh`; report template in `report/`.

## Build, Test, and Development Commands
- `make setup` — create Conda envs and install pre-commit hooks.
- `make smoke` — run Snakemake then Nextflow on the tiny dataset, asserting core artifacts.
- `make run` — execute the engine selected in `config/params.yaml`.
- `make lint` — run ruff, yamllint, Snakemake lint, and Nextflow stub run.
- `tests/run_smoke.sh` — standalone smoke script (accepts `PIPELINE_ENGINE=snakemake|nextflow`).

## Coding Style & Naming Conventions
- **Python**: 4-space indent, ruff/black profiles defined in `pyproject.toml`; module names are snake_case.
- **R**: 2-space indent preferred; keep scripts functional and parameter-driven via `optparse`.
- **YAML**: 2-space indent; keep keys lowercase with hyphenated lists.
- Follow directory naming already in place (`results/<stage>/`).

## Testing Guidelines
- Primary check is the smoke workflow (`tests/run_smoke.sh`); ensure both engines pass before opening a PR.
- Add larger regression tests only if they fit within CI resource limits (<5 MB inputs).
- Name new fixtures descriptively (`tests/data/fastq/<sample>_<read>.fastq.gz`).

## Commit & Pull Request Guidelines
- Write commits in present-tense imperative (“Add Nextflow multiqc module”).
- Keep related changes together (config + rules + scripts).
- PRs should include: summary of pipeline stage touched, validation evidence (`make smoke` or targeted logs), and links to relevant issues.
- Attach sample output paths or hashes when expecting reviewers to verify determinism.

## Project Build Plan & Progress
- Use the checkboxes below to track delivery. When a task is finished, flip `[ ]` to `[x]` and mention the validating command or CI run in the same line.
- Revisit this list after every substantive PR so new work is reflected.

### Cross-Engine Validation
- [ ] Run `make smoke` locally (or in a controlled VM) and attach logs confirming both engines succeed.
- [ ] Add deterministic hash verification to CI (extend `.github/workflows/ci.yml` to run `make hash` after each engine).

### Documentation & Onboarding
- [ ] Produce a step-by-step Quickstart video or annotated screenshot tour of `results/report.html` for non-technical collaborators.
- [ ] Expand README with a short FAQ covering common errors (missing envs, failed Salmon index).

### Pipeline Hardening
- [ ] Test the pipelines on a modest real dataset (≥4 samples per condition) and note runtime/resource usage.
- [ ] Implement optional download automation for human references (script or Snakemake rule gated by config).
- [ ] Verify Slurm profiles on a staging cluster; document any required `--cluster-config` tweaks.

## Configuration & Operational Notes
- Always keep Snakemake and Nextflow consuming the same `config/params.yaml`; avoid duplicating settings in engine-specific files.
- Rebuild Salmon indexes when references change (`scripts/build_salmon_index.sh -t … -a … -d … -o …`).
- Use `config/profiles/slurm.*` profiles for HPC submission; never hard-code cluster options in rules or processes.
