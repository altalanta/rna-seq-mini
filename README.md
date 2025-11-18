# 🔬 RNA-seq Mini Platform

This repository contains the platform components for the RNASEQ-MINI project, including a web application, an API, and advanced deployment configurations.

The core bioinformatics pipeline is developed and maintained as a separate component in the `pipeline-core/` directory.

## 🚀 Quick Start

### 1. Setup the Environment

Ensure you have Conda/Mamba installed. Then, set up the necessary environments by running the setup command from the core pipeline:

```bash
cd pipeline-core
make setup
cd ..
```

### 2. Run the Analysis Pipeline

To run the bioinformatics pipeline, use the top-level `make` command, which delegates to the pipeline's internal `Makefile`:

```bash
make pipeline-run
```

All pipeline-related commands are prefixed with `pipeline-`. For a full list of pipeline commands, run:

```bash
cd pipeline-core && make help
```

### 3. Launch the Platform

The platform provides a web application and an API to interact with the pipeline results.

**Launch the Web App:**
```bash
make web-app
```

**Launch the API Server:**
```bash
make api-server
```

## 🏗️ Architecture

The repository is structured to separate the core scientific workflow from the surrounding platform infrastructure:

- **`pipeline-core/`**: Contains the standalone Snakemake pipeline, including all scripts, configuration, and testing required to run the analysis. This component is designed to be portable and independently executable.

- **`api/`**: A REST API for programmatically interacting with the pipeline.

- **`web_app/`**: An interactive web-based dashboard for exploring and visualizing analysis results.

- **`scripts/`**: Contains Python scripts that support the platform, such as advanced analytics, cloud deployment, and collaboration features.

- **`deployments/` & `helm/`**: Configuration for deploying the platform using Docker Compose and Kubernetes.

This modular design allows for independent development and testing of the core pipeline and the platform components.
