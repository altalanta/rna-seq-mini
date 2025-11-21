# Makefile for the RNASEQ-MINI Platform and Core Pipeline
# ==========================================================
# This Makefile provides a unified command-line interface for both the
# bioinformatics pipeline and the surrounding platform components.
#
# Targets are namespaced to clarify their scope:
#   - `pipeline-*`: Commands that operate on the core pipeline (in pipeline-core/).
#   - `platform-*`: Commands related to the web app, API, and other platform services.
#

.PHONY: help \
	pipeline-setup pipeline-lint pipeline-smoke pipeline-run pipeline-validate pipeline-clean \
	platform-serve-webapp platform-serve-api

help: ## Show this help
	@echo "Makefile for the RNASEQ-MINI Platform"
	@echo ""
	@echo "Usage: make [target]"
	@echo ""
	@echo "Targets:"
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-25s\033[0m %s\n", $$1, $$2}'

# --- Pipeline Commands ---
PIPELINE_DIR := pipeline-core
PIPELINE_MAKE := $(MAKE) -C $(PIPELINE_DIR)

pipeline-setup: ## Setup Conda environments for the pipeline
	@echo "--> Setting up pipeline environments in $(PIPELINE_DIR)..."
	@$(PIPELINE_MAKE) setup

pipeline-lint: ## Lint the core pipeline
	@echo "--> Linting pipeline in $(PIPELINE_DIR)..."
	@$(PIPELINE_MAKE) lint

pipeline-smoke: ## Run the core pipeline smoke test
	@echo "--> Running smoke test in $(PIPELINE_DIR)..."
	@$(PIPELINE_MAKE) smoke

pipeline-run: ## Run the core analysis pipeline
	@echo "--> Executing pipeline in $(PIPELINE_DIR)..."
	@$(PIPELINE_MAKE) run

pipeline-validate: ## Validate the reproducibility of the pipeline
	@echo "--> Validating pipeline in $(PIPELINE_DIR)..."
	@$(PIPELINE_MAKE) validate

pipeline-clean: ## Clean all pipeline-generated files
	@echo "--> Cleaning pipeline in $(PIPELINE_DIR)..."
	@$(PIPELINE_MAKE) clean

# --- Platform Commands ---
PY_ENV=rnaseq-mini-base

platform-serve-webapp: ## Launch the web application
	@echo "--> Launching web application..."
	@conda run -n $(PY_ENV) python web_app/app.py

platform-serve-api: ## Launch the API server
	@echo "--> Launching API server..."
	@conda run -n $(PY_ENV) python api/server.py
