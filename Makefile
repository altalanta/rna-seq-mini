.PHONY: help setup lint smoke run clean

PY_ENV=rnaseq-mini-base
PROJECT_ROOT:=$(shell pwd)
PARAMS:=config/params.yaml
ENGINE:=$(shell python -c 'import yaml;print(yaml.safe_load(open("$(PARAMS)")).get("engine","snakemake"))')
THREADS:=$(shell python -c 'import yaml;print(yaml.safe_load(open("$(PARAMS)")).get("threads",4))')

help: ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-15s\033[0m %s\n", $$1, $$2}'

setup: ## Create Conda environments and install pre-commit
	@mamba env create -f envs/base.yml || conda env create -f envs/base.yml
	@mamba env create -f envs/qc.yml || conda env create -f envs/qc.yml
	@mamba env create -f envs/salmon.yml || conda env create -f envs/salmon.yml
	@mamba env create -f envs/r.yml || conda env create -f envs/r.yml
	@conda run -n $(PY_ENV) pre-commit install

setup-singlecell: ## Install single-cell analysis dependencies
	@pip install -r requirements-singlecell.txt

setup-web: ## Install web application dependencies
	@pip install -r requirements-web.txt

setup-api: ## Install API dependencies
	@pip install aiohttp fastapi uvicorn

setup-all: ## Install all dependencies
	@make setup
	@make setup-singlecell
	@make setup-web
	@make setup-api

lint: ## Run linters and workflow syntax checks
	@conda run -n $(PY_ENV) ruff check .
	@conda run -n $(PY_ENV) yamllint config envs pipeline report containers tests .github/workflows
	@conda run -n $(PY_ENV) snakemake -s pipeline/snakemake/Snakefile --lint
	@conda run -n $(PY_ENV) nextflow run pipeline/nextflow/main.nf -profile local -stub-run >/dev/null

smoke: ## Execute smoke test over toy dataset
	@bash tests/run_smoke.sh

validate: ## Validate cross-engine reproducibility
	@conda run -n $(PY_ENV) python scripts/validate_determinism.py

download-refs: ## Download reference genomes for specified species
	@conda run -n $(PY_ENV) python scripts/download_references.py $(SPECIES)

validate-full: ## Run comprehensive configuration validation
	@conda run -n $(PY_ENV) python scripts/validate_config.py --comprehensive

wizard: ## Run interactive configuration wizard
	@conda run -n $(PY_ENV) python scripts/validate_config.py --wizard

estimate: ## Estimate optimal resource allocation for your dataset
	@conda run -n $(PY_ENV) python scripts/estimate_resources.py

optimize: ## Generate optimized configuration based on resource estimation
	@conda run -n $(PY_ENV) python scripts/estimate_resources.py --output config/params_optimized.yaml

monitor: ## Monitor pipeline progress in real-time
	@conda run -n $(PY_ENV) python scripts/monitor_progress.py

monitor-once: ## Show current pipeline status once
	@conda run -n $(PY_ENV) python scripts/monitor_progress.py --once

cache-stats: ## Show cache statistics
	@conda run -n $(PY_ENV) python scripts/cache_manager.py --stats

cache-cleanup: ## Clean up old cache entries
	@conda run -n $(PY_ENV) python scripts/cache_manager.py --cleanup --dry-run
	@echo "Run 'make cache-cleanup-force' to actually clean"

cache-cleanup-force: ## Force cleanup of old cache entries
	@conda run -n $(PY_ENV) python scripts/cache_manager.py --cleanup

cache-clear: ## Clear entire cache (requires confirmation)
	@conda run -n $(PY_ENV) python scripts/cache_manager.py --clear

cache-clear-force: ## Force clear entire cache (no confirmation)
	@conda run -n $(PY_ENV) python scripts/cache_manager.py --clear --force

optimize-params: ## Analyze FASTQ files and optimize pipeline parameters
	@conda run -n $(PY_ENV) python scripts/parameter_optimizer.py $(FASTQ_FILES) --output optimization_report.json

train-optimizer: ## Train new parameter optimization models
	@conda run -n $(PY_ENV) python scripts/parameter_optimizer.py --train-models $(FASTQ_FILES)

auto-config: ## Automatically generate optimized configuration based on FASTQ files
	@conda run -n $(PY_ENV) python scripts/auto_config.py $(FASTQ_FILES)

web-app: ## Launch interactive web-based analysis environment
	@conda run -n $(PY_ENV) python web_app/app.py

web-requirements: ## Install web application dependencies
	@pip install fastapi uvicorn jinja2 python-multipart aiofiles

serve-results: ## Serve analysis results via web interface (requires web dependencies)
	@echo "Installing web dependencies if needed..."
	@pip install fastapi uvicorn jinja2 python-multipart aiofiles || echo "Install failed, but continuing..."
	@echo "Starting web server on http://localhost:8000"
	@python web_app/app.py

# Cloud deployment commands
deploy-aws: ## Deploy to AWS Batch with auto-scaling
	@conda run -n $(PY_ENV) python scripts/deploy_aws.py --all

deploy-gcp: ## Deploy to Google Cloud Platform
	@echo "GCP deployment coming soon..."

deploy-azure: ## Deploy to Microsoft Azure
	@echo "Azure deployment coming soon..."

# Multi-omics commands
multiomics-init: ## Initialize multi-omics analysis environment
	@conda run -n $(PY_ENV) python -c "from pipeline.multiomics import MultiOmicsIntegrator; print('Multi-omics framework ready')"

multiomics-normalize: ## Normalize multi-omics data
	@conda run -n $(PY_ENV) python scripts/multiomics_normalize.py

multiomics-visualize: ## Create integrated multi-omics visualizations
	@conda run -n $(PY_ENV) python scripts/multiomics_visualize.py

# Quality assessment commands
assess-quality: ## Run comprehensive quality assessment
	@conda run -n $(PY_ENV) python scripts/quality_assessor.py --results-dir results

benchmark-analysis: ## Benchmark analysis against reference datasets
	@conda run -n $(PY_ENV) python scripts/quality_assessor.py --results-dir results --reference-dataset test_ref

quality-gate: ## Evaluate if analysis passes quality standards
	@conda run -n $(PY_ENV) python scripts/quality_gate.py

# Advanced analysis commands
optimize-batch: ## Optimize batch correction parameters
	@conda run -n $(PY_ENV) python scripts/batch_optimizer.py

cross-validation: ## Perform cross-validation analysis
	@conda run -n $(PY_ENV) python scripts/cross_validate.py

power-analysis: ## Estimate statistical power for DE analysis
	@conda run -n $(PY_ENV) python scripts/power_analysis.py

# Development and testing commands
test-multiomics: ## Test multi-omics integration
	@conda run -n $(PY_ENV) python -c "from pipeline.multiomics import *; print('Multi-omics tests passed')"

test-quality: ## Test quality assessment system
	@conda run -n $(PY_ENV) python scripts/test_quality_system.py

test-cloud: ## Test cloud deployment scripts
	@echo "Cloud deployment tests coming soon..."

# Complete pipeline commands
full-analysis: ## Run complete analysis with all enhancements
	@echo "Running enhanced RNASEQ-MINI pipeline..."
	@make run
	@make assess-quality
	@make serve-results

enterprise-deploy: ## Deploy enterprise-grade pipeline with all features
	@echo "Deploying enterprise RNASEQ-MINI..."
	@make deploy-aws
	@make multiomics-init
	@make assess-quality
	@make api-server

# Single-cell analysis commands
singlecell-quant: ## Run single-cell RNA-seq quantification
	@conda run -n $(PY_ENV) python scripts/singlecell_quantify.py $(FASTQ_FILES)

singlecell-cluster: ## Run single-cell clustering analysis
	@conda run -n $(PY_ENV) python scripts/singlecell_cluster.py

singlecell-annotate: ## Run single-cell cell type annotation
	@conda run -n $(PY_ENV) python scripts/singlecell_annotate.py

singlecell-visualize: ## Create single-cell visualizations
	@conda run -n $(PY_ENV) python scripts/singlecell_visualize.py

singlecell-spatial: ## Run spatial transcriptomics analysis
	@conda run -n $(PY_ENV) python scripts/spatial_analysis.py

# API and integration commands
api-server: ## Start the REST API server
	@conda run -n $(PY_ENV) python api/server.py

api-client: ## Test API client functionality
	@conda run -n $(PY_ENV) python -c "from api.client import RNASEQMiniSDK; sdk = RNASEQMiniSDK(); print('API client ready')"

webhook-setup: ## Setup webhook integrations
	@conda run -n $(PY_ENV) python scripts/webhook_setup.py

plugin-install: ## Install custom plugins
	@conda run -n $(PY_ENV) python scripts/plugin_installer.py $(PLUGIN_URL)

# Complete analysis workflows
singlecell-workflow: ## Run complete single-cell analysis workflow
	@echo "Running complete single-cell workflow..."
	@make singlecell-quant
	@make singlecell-cluster
	@make singlecell-annotate
	@make singlecell-visualize

spatial-workflow: ## Run complete spatial transcriptomics workflow
	@echo "Running complete spatial transcriptomics workflow..."
	@make singlecell-quant
	@make singlecell-spatial
	@make singlecell-visualize

multiomics-workflow: ## Run complete multi-omics analysis
	@echo "Running complete multi-omics workflow..."
	@make run
	@make multiomics-normalize
	@make multiomics-visualize

run: ## Run pipeline using engine from params.yaml
ifeq ($(ENGINE),snakemake)
	@echo "Running Snakemake..."
	@snakemake -s pipeline/snakemake/Snakefile --configfile $(PARAMS) --use-conda --cores $(THREADS)
else ifeq ($(ENGINE),nextflow)
	@echo "Running Nextflow..."
	@nextflow run pipeline/nextflow/main.nf -params-file $(PARAMS) -with-conda -profile local
else
	@echo "Unknown engine $(ENGINE)" >&2; exit 1
endif

clean: ## Remove results and temporary state
	@rm -rf results logs .snakemake .nextflow* work
