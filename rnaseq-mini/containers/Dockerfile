FROM mambaorg/micromamba:1.5.1

USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
    procps \
    tini \
    && rm -rf /var/lib/apt/lists/*

COPY containers/environment.lock.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean -a -y

ENV PATH=/opt/conda/bin:$PATH
WORKDIR /workspace
CMD ["/bin/bash"]
