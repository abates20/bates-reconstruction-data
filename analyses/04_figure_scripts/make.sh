#!/bin/bash

# Filepaths
SCRIPT_DIR=$(dirname "$(realpath "${BASH_SOURCE[0]}")")
PROJECT_DIR=$(cd "${SCRIPT_DIR}/../../" && pwd)
FIGURES="${PROJECT_DIR}/dist/figures"
BIN="${PROJECT_DIR}/bin"
DATA="${PROJECT_DIR}/data"
CACHE="${SCRIPT_DIR}/cache"


# Check provided argument and print help if not provided or invalid
help () {
    echo "Run script(s) to produce figures"
    echo ""
    echo "usage: bash make.sh FIGURE"

    exit 1
}

if [ -z "$1" ]; then
    help
elif [[ ! -e "${SCRIPT_DIR}/figure_${1}.py" && "$1" != "all" ]]; then
    echo "Invalid argument: '${1}'. Run script without arguments for help."
    exit 1
fi


# Make sure needed directories exist
mkdir -p "$FIGURES"
mkdir -p "$CACHE"

# Ignore generated files in git
echo "*" > "${CACHE}/.gitignore"

# Load .env variables (if .env exists)
ENV_FILE="${PROJECT_DIR}/.env"
if [ -f "$ENV_FILE" ]; then
    set -a
    source "$ENV_FILE"
    set +a
fi

# Set default python path if not set in .env
if [ -z "$PYTHON" ]; then
    PYTHON="${PROJECT_DIR}/.venv/bin/python"
fi


# Run specified figure scripts
if [[ "$1" == "2a" || "$1" == "all" ]]; then
    echo "Generating figure 2a..."
    "$PYTHON" "${SCRIPT_DIR}/figure_2a.py" \
        --model-info "${DATA}/model_info.csv" \
        --bigg-data "${DATA}/bigg_bacteria.csv" \
        --exp-data "${DATA}/substrate_utilization_data.csv" \
        --out "${FIGURES}/figure_2a.svg" \
        --cache "$CACHE"
fi
