.PHONY: help pipeline-run pipeline-lint pipeline-smoke

help: ## Show this help
	@echo "Makefile for the RNASEQ-MINI Platform"
	@echo ""
	@echo "Usage: make [target]"
	@echo ""
	@echo "Targets:"
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'
	@echo ""
	@echo "Pipeline commands are delegated to the core pipeline."
	@echo "Run 'cd pipeline-core && make help' to see all available pipeline targets."

pipeline-run: ## Run the core analysis pipeline
	@echo "--> Executing pipeline in pipeline-core/..."
	@$(MAKE) -C pipeline-core run

pipeline-lint: ## Lint the core pipeline
	@echo "--> Linting pipeline in pipeline-core/..."
	@$(MAKE) -C pipeline-core lint

pipeline-smoke: ## Run the core pipeline smoke test
	@echo "--> Running smoke test in pipeline-core/..."
	@$(MAKE) -C pipeline-core smoke

web-app: ## Launch interactive web-based analysis environment
	@conda run -n rnaseq-mini-core python web_app/app.py

api-server: ## Start the REST API server
	@conda run -n rnaseq-mini-core python api/server.py
