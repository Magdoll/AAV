#!/bin/bash
# Usage: run-workflow-local.sh
set -o pipefail
set -o errexit

SCRIPT_DIR=$(cd -P -- "$(dirname -- "$0")" && pwd -P)

export NXF_VER=23.10.1

# Run nextflow locally
nextflow run -profile local main.nf -params-file params-local-small.json
