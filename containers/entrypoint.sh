#!/usr/bin/env bash
set -euo pipefail

# Default to the micromamba user, but allow overrides
USER_ID=${LOCAL_USER_ID:-9001}
GROUP_ID=${LOCAL_GROUP_ID:-9001}

# Create a user and group with the specified IDs
echo "Starting with UID : $USER_ID, GID: $GROUP_ID"
if ! getent group micromamba > /dev/null; then
    groupadd -g "$GROUP_ID" micromamba
fi
if ! getent passwd micromamba > /dev/null; then
    useradd --shell /bin/bash -u "$USER_ID" -g micromamba micromamba
fi
export HOME=/home/micromamba

# If the first argument is 'snakemake' or 'nextflow', execute the pipeline
if [ "$1" = "snakemake" ] || [ "$1" = "nextflow" ]; then
    # Ensure the user has permissions for the mounted workspace
    chown -R micromamba:micromamba /workspace
    # Execute the command as the micromamba user
    exec /usr/sbin/gosu micromamba "$@"
# If the command is 'bash' or 'sh', start an interactive shell
elif [ "$1" = "bash" ] || [ "$1" = "sh" ]; then
    exec /usr/sbin/gosu micromamba "$@"
# Otherwise, show a help message
else
    echo "Usage: docker run --rm -it rnaseq-mini <command>"
    echo
    echo "Commands:"
    echo "  snakemake [args...]   - Run the Snakemake pipeline with specified arguments."
    echo "  nextflow [args...]    - Run the Nextflow pipeline with specified arguments."
    echo "  bash/sh               - Start an interactive shell inside the container."
    echo
    echo "Example:"
    echo "  docker run --rm -it -v \$(pwd)/config:/workspace/config -v \$(pwd)/results:/workspace/results rnaseq-mini snakemake --cores 2"
    exit 1
fi





