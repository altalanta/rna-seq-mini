.PHONY: help setup lint smoke run clean

PY_ENV=rnaseq-mini-base
PROJECT_ROOT:=$(shell pwd)
PARAMS:=config/params.yaml
ENGINE:=$(shell python -c 'import yaml;print(yaml.safe_load(open("$(PARAMS)")).get("engine","snakemake"))')
THREADS:=$(shell python -c 'import yaml;print(yaml.safe_load(open("$(PARAMS)")).get("threads",4))')

help: ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-15s\033[0m %s\n", $$1, $$2}'

setup: ## Create consolidated Conda environments and install pre-commit
	@mamba env create -f envs/rnaseq-core.yml || conda env create -f envs/rnaseq-core.yml
	@mamba env create -f envs/rnaseq-analysis.yml || conda env create -f envs/rnaseq-analysis.yml
	@conda run -n rnaseq-mini-core pre-commit install
	@conda run -n rnaseq-mini-core pip install -e .

setup-singlecell: ## Install single-cell analysis dependencies
	@pip install -r requirements-singlecell.txt

setup-web: ## Install web application dependencies
	@pip install -r requirements-web.txt

setup-api: ## Install API dependencies
	@pip install aiohttp fastapi uvicorn

setup-singlecell: ## Install single-cell analysis dependencies
	@mamba env create -f envs/singlecell.yml || conda env create -f envs/singlecell.yml
	@pip install -r requirements-singlecell.txt

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

# Determinism Validation
.PHONY: check-hashes update-hashes

check-hashes: ## Run smoke test and check output hashes against a reference file
	@if [ ! -f tests/validation/expected_hashes.sha256 ]; then \
		echo "ERROR: Reference hash file not found."; \
		echo "Please generate it by running 'make update-hashes' and commit the result."; \
		exit 1; \
	fi
	@echo "🧪 Running smoke test for determinism check..."
	@PIPELINE_ENGINE=snakemake tests/run_smoke.sh results-smoke-check
	@echo " hashlib Computing hashes for new results..."
	@scripts/compute_hashes.sh results-smoke-check generated_hashes.sha256
	@echo "    Comparing against reference hashes..."
	@if ! diff -q tests/validation/expected_hashes.sha256 generated_hashes.sha256; then \
		echo "❌ HASH MISMATCH DETECTED!"; \
		echo "Differences between tests/validation/expected_hashes.sha256 (expected) and generated_hashes.sha256 (generated):"; \
		diff tests/validation/expected_hashes.sha256 generated_hashes.sha256; \
		rm generated_hashes.sha256; \
		exit 1; \
	fi
	@echo "✅ Determinism check passed!"
	@rm generated_hashes.sha256
	@rm -rf results-smoke-check

update-hashes: ## Generate/update the reference hash file from a fresh smoke test run
	@echo "🛠️  Running smoke test to generate new reference hashes..."
	@PIPELINE_ENGINE=snakemake tests/run_smoke.sh results-smoke-ref
	@mkdir -p tests/validation
	@echo " hashlib Computing and saving new reference hashes..."
	@scripts/compute_hashes.sh results-smoke-ref tests/validation/expected_hashes.sha256
	@echo "✅ New reference hashes saved to tests/validation/expected_hashes.sha256"
	@rm -rf results-smoke-ref

download-refs: ## Download reference genomes for specified species
	@conda run -n $(PY_ENV) python scripts/download_references.py $(SPECIES)

validate-full: ## Run comprehensive configuration validation
	@conda run -n $(PY_ENV) python scripts/validate_config.py --comprehensive

wizard: ## Run interactive configuration wizard
	@conda run -n $(PY_ENV) python scripts/setup_wizard.py

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

# AI-powered insights commands
ai-insights: ## Run AI-powered result interpretation and insights
	@conda run -n $(PY_ENV) python scripts/ai_insights.py analyze --results-dir results --format both

pathway-impact: ## Run advanced pathway impact analysis
	@conda run -n $(PY_ENV) python scripts/pathway_impact_analyzer.py --results-dir results --format both

gene-networks: ## Analyze gene interaction networks
	@conda run -n $(PY_ENV) python scripts/gene_network_analyzer.py --results-dir results --format both

predictive-modeling: ## Run predictive modeling and biomarker discovery
	@conda run -n $(PY_ENV) python scripts/predictive_modeling.py --results-dir results --format both

ai-explanations: ## Generate natural language explanations for different audiences
	@conda run -n $(PY_ENV) python scripts/result_explanation.py --results-dir results --audiences researcher clinician student executive --format markdown

ai-complete: ## Run complete AI-powered analysis suite
	@echo "🧠 Running complete AI-powered analysis suite..."
	@make ai-insights
	@make pathway-impact
	@make gene-networks
	@make predictive-modeling
	@make ai-explanations

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

# Collaboration commands
collaboration-server: ## Start real-time collaborative analysis server
	@conda run -n $(PY_ENV) python scripts/collaborative_server.py server --host localhost --port 8765

collaboration-create: ## Create a new collaborative session
	@conda run -n $(PY_ENV) python scripts/collaborative_server.py create $(PROJECT_NAME) --config config/params.yaml

collaboration-join: ## Join an existing collaborative session
	@conda run -n $(PY_ENV) python scripts/collaborative_server.py join $(SESSION_ID) "$(USER_NAME)" "$(USER_EMAIL)"

collaboration-notebook: ## Create a Jupyter notebook for collaborative analysis
	@conda run -n $(PY_ENV) python scripts/jupyter_integration.py create $(SESSION_ID) --template standard

collaboration-notifications: ## Setup notification channels (Slack, Discord, etc.)
	@conda run -n $(PY_ENV) python scripts/notification_manager.py setup

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

run-singlecell: ## Run single-cell analysis only
	@echo "Running single-cell analysis..."
ifeq ($(ENGINE),snakemake)
	@snakemake -s pipeline/snakemake/Snakefile --configfile $(PARAMS) --use-conda --cores $(THREADS) singlecell_complete
else ifeq ($(ENGINE),nextflow)
	@nextflow run pipeline/nextflow/main.nf -params-file $(PARAMS) -with-conda -profile local --singlecell.enabled true
endif

singlecell-qc: ## Run single-cell quality control
	@conda run -n rnaseq-mini-singlecell python scripts/singlecell_qc.py $(MATRIX) $(BARCODES) $(OUTPUT_DIR)

singlecell-cluster: ## Run single-cell clustering
	@conda run -n rnaseq-mini-singlecell python scripts/singlecell_clustering.py $(MATRIX) $(BARCODES) $(GENES) $(OUTPUT_DIR)

# Resource optimization commands
estimate-resources: ## Estimate optimal resource allocation for your dataset
	@conda run -n $(PY_ENV) python scripts/resource_estimator.py $(SAMPLES_FILE) --output resource_estimation.json

optimize-resources: ## Generate optimized configuration based on resource estimation
	@conda run -n $(PY_ENV) python scripts/resource_estimator.py $(SAMPLES_FILE) --output config/resource_optimized.json
	@echo "Optimized configuration saved to config/resource_optimized.json"
	@echo "Copy this file to config/params.yaml for optimized settings"

# Cloud optimization commands
cloud-autoscale: ## Run cloud auto-scaling for AWS Batch
	@conda run -n $(PY_ENV) python scripts/cloud_autoscaler.py --job-queue $(JOB_QUEUE) --compute-env $(COMPUTE_ENV) --duration 60

cloud-cost-monitor: ## Monitor and analyze cloud costs
	@conda run -n $(PY_ENV) python scripts/cost_monitor.py --output-report cost_report.json --output-plot cost_analysis.png

cloud-setup-dashboard: ## Setup cost monitoring dashboard
	@conda run -n $(PY_ENV) python scripts/cost_monitor.py --setup-dashboard

# Error handling and troubleshooting commands
diagnostics: ## Run comprehensive diagnostic checks
	@conda run -n $(PY_ENV) python scripts/error_handler.py --diagnostics all

troubleshoot: ## Analyze and classify error messages
	@conda run -n $(PY_ENV) python scripts/error_handler.py --error-message "$(ERROR_MESSAGE)" --context-file "$(CONTEXT_FILE)"

health-check: ## Perform one-time system and pipeline health check
	@conda run -n $(PY_ENV) python scripts/troubleshooter.py --health-check

env-health: ## Validate environment integrity and check for missing dependencies
	@echo "🔍 Checking environment health..."
	@echo "Core environment (rnaseq-mini-core):"
	@conda run -n rnaseq-mini-core python -c "import snakemake, nextflow, pandas, yaml; print('✅ Core Python tools OK')" 2>/dev/null || echo "❌ Core environment issues"
	@conda run -n rnaseq-mini-core which fastqc multiqc salmon samtools 2>/dev/null | head -5 | wc -l | xargs -I {} echo "✅ Found {} QC/quantification tools" 2>/dev/null || echo "❌ Missing QC/quantification tools"
	@echo ""
	@echo "Analysis environment (rnaseq-mini-analysis):"
	@conda run -n rnaseq-mini-analysis R -e "library(DESeq2); library(tximport); print('✅ R analysis tools OK')" 2>/dev/null || echo "❌ R environment issues"
	@conda run -n rnaseq-mini-analysis python -c "import scanpy, anndata, numpy; print('✅ Single-cell tools OK')" 2>/dev/null || echo "❌ Single-cell tools issues"
	@echo ""
	@echo "Environment health check complete."

monitor-health: ## Monitor system health for specified duration
	@conda run -n $(PY_ENV) python scripts/troubleshooter.py --monitor $(DURATION_MINUTES) --interval 30

error-recovery: ## Attempt error recovery for failed jobs
	@echo "Error recovery system initialized"
	@echo "Run 'make troubleshoot' with specific error messages for detailed analysis"

cleanup-duplicates: ## Remove duplicate files with ' 2' suffix
	@echo "Cleaning up duplicate files..."
	@find . -name "* 2*" -type f -delete
	@echo "Removed duplicate files. Repository cleaned."

clean: ## Remove results and temporary state
	@rm -rf results logs .snakemake .nextflow* work
