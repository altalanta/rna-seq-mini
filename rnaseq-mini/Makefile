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

lint: ## Run linters and workflow syntax checks
	@conda run -n $(PY_ENV) ruff check .
	@conda run -n $(PY_ENV) yamllint config envs pipeline report containers tests .github/workflows
	@conda run -n $(PY_ENV) snakemake -s pipeline/snakemake/Snakefile --lint
	@conda run -n $(PY_ENV) nextflow run pipeline/nextflow/main.nf -profile local -stub-run >/dev/null

smoke: ## Execute smoke test over toy dataset
	@bash tests/run_smoke.sh

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
