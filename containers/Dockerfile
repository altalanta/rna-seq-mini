# Multi-stage Dockerfile for RNASEQ-MINI
FROM mambaorg/micromamba:1.5.1 as builder

# Install system dependencies
USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
    procps \
    tini \
    gosu \
    git \
    wget \
    curl \
    unzip \
    && rm -rf /var/lib/apt/lists/*

# Create a single, unified environment
COPY --chown=micromamba:micromamba containers/environment.yml /tmp/env.yml
RUN micromamba create -y -n rnaseq-mini -f /tmp/env.yml && \
    micromamba clean -a -y

# Copy source code and install the package
USER micromamba
WORKDIR /tmp/install
COPY --chown=micromamba:micromamba src ./src
COPY --chown=micromamba:micromamba pyproject.toml .
RUN micromamba run -n rnaseq-mini pip install .


# Production stage - minimal runtime image
FROM mambaorg/micromamba:1.5.1

# Install minimal system dependencies
USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
    procps \
    tini \
    gosu \
    && rm -rf /var/lib/apt/lists/*

# Copy the unified environment from the builder stage
COPY --from=builder /opt/conda/envs/rnaseq-mini /opt/conda/envs/rnaseq-mini

# Copy the entire project directory into the container
COPY . /app
WORKDIR /app

# Copy the entrypoint script and make it executable
COPY containers/entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /usr/local/bin/entrypoint.sh

# Set environment variables to use the new environment by default
ENV PATH=/opt/conda/envs/rnaseq-mini/bin:$PATH
ENV MAMBA_DEFAULT_ENV=rnaseq-mini
ENV MAMBA_ROOT_PREFIX=/opt/conda
WORKDIR /workspace

# Use tini as the subreaper and run the entrypoint
ENTRYPOINT ["/usr/bin/tini", "--", "/usr/local/bin/entrypoint.sh"]

# Default command
CMD ["bash"]
