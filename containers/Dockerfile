# Multi-stage Dockerfile for RNASEQ-MINI
FROM mambaorg/micromamba:1.5.1 as builder

# Install system dependencies
USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
    procps \
    tini \
    git \
    wget \
    curl \
    unzip \
    && rm -rf /var/lib/apt/lists/*

# Create base environment with core tools
COPY envs/rnaseq-core.yml /tmp/core.yml
RUN micromamba create -y -n rnaseq-core -f /tmp/core.yml && \
    micromamba clean -a -y

# Create analysis environment with R and single-cell tools
COPY envs/rnaseq-analysis.yml /tmp/analysis.yml
RUN micromamba create -y -n rnaseq-analysis -f /tmp/analysis.yml && \
    micromamba clean -a -y

# Production stage - minimal runtime image
FROM mambaorg/micromamba:1.5.1

# Install minimal system dependencies
USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
    procps \
    tini \
    && rm -rf /var/lib/apt/lists/*

# Copy environments from builder stage
COPY --from=builder /opt/conda/envs/rnaseq-core /opt/conda/envs/rnaseq-core
COPY --from=builder /opt/conda/envs/rnaseq-analysis /opt/conda/envs/rnaseq-analysis

# Create symlink for default environment
RUN ln -sf /opt/conda/envs/rnaseq-core /opt/conda/envs/default

# Set environment variables
ENV PATH=/opt/conda/envs/rnaseq-core/bin:/opt/conda/envs/rnaseq-analysis/bin:$PATH
ENV MAMBA_ROOT_PREFIX=/opt/conda
WORKDIR /workspace

# Default command
CMD ["/bin/bash"]
