#!/bin/bash

# Filepaths
SCRIPT_DIR=$(dirname "$(realpath "${BASH_SOURCE[0]}")")
PROJECT_DIR=$(cd "${SCRIPT_DIR}/../../" && pwd)
ANALYSES="${PROJECT_DIR}/analyses"
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
        "${ANALYSES}/01_reconstructor_vs_carveme/results/reaction_matrix.csv" \
        --set-info "${ANALYSES}/01_reconstructor_vs_carveme/datasets.csv" \
        --out "${FIGURES}/figure_2a.svg" \
        --cache "$CACHE"
fi

if [[ "$1" == "2b" || "$1" == "all" ]]; then
    echo "Generating figure 2b..."
    "$PYTHON" "${SCRIPT_DIR}/figure_2b.py" \
        "${ANALYSES}/01_reconstructor_vs_carveme/results/reaction_matrix.csv" \
        --set-info "${ANALYSES}/01_reconstructor_vs_carveme/datasets.csv" \
        --out "${FIGURES}/figure_2b.svg" \
        --cache "$CACHE"
fi

if [[ "$1" == "2c" || "$1" == "all" ]]; then
    echo "Generating figure 2c..."
    "$PYTHON" "${SCRIPT_DIR}/figure_2c.py" \
        "${ANALYSES}/01_reconstructor_vs_carveme/results/growth_data.csv" \
        --set-info "${ANALYSES}/01_reconstructor_vs_carveme/datasets.csv" \
        --out "${FIGURES}/figure_2c.svg" \
        --cache "$CACHE"
fi

if [[ "$1" == "3b" || "$1" == "all" ]]; then
    echo "Generating figure 3b..."
    "$PYTHON" "${SCRIPT_DIR}/figure_3b.py" \
        "${ANALYSES}/02_annotation_uncertainty/results/reaction_matrix.csv" \
        --set-info "${ANALYSES}/02_annotation_uncertainty/datasets.csv" \
        --model-info "${DATA}/model_info.csv" \
        --out "${FIGURES}/figure_3b.svg" \
        --cache "$CACHE"
fi

if [[ "$1" == "3c" || "$1" == "all" ]]; then
    echo "Generating figure 3c..."
    "$PYTHON" "${SCRIPT_DIR}/figure_3c.py" \
        "${ANALYSES}/02_annotation_uncertainty/results/growth_data.csv" \
        --set-info "${ANALYSES}/02_annotation_uncertainty/datasets.csv" \
        --model-info "${DATA}/model_info.csv" \
        --out "${FIGURES}/figure_3c.svg" \
        --cache "$CACHE"
fi

if [[ "$1" == "5c" || "$1" == "all" ]]; then
    echo "Generating figure 5c..."
    "$PYTHON" "${SCRIPT_DIR}/figure_5c.py" \
        "${ANALYSES}/03_taxonomic_representation/results/reaction_matrix.csv" \
        --set-info "${ANALYSES}/03_taxonomic_representation/results/dataset_metrics.csv" \
        --model-info "${DATA}/model_info.csv" \
        --out "${FIGURES}/figure_5c.svg" \
        --cache "$CACHE"
fi

if [[ "$1" == "5d" || "$1" == "all" ]]; then
    echo "Generating figure 5d..."
    "$PYTHON" "${SCRIPT_DIR}/figure_5d.py" \
        "${ANALYSES}/03_taxonomic_representation/results/reaction_matrix.csv" \
        --set-info "${ANALYSES}/03_taxonomic_representation/results/dataset_metrics.csv" \
        --model-info "${DATA}/model_info.csv" \
        --out "${FIGURES}/figure_5d.svg" \
        --cache "$CACHE"
fi

if [[ "$1" == "5e" || "$1" == "all" ]]; then
    echo "Generating figure 5e..."
    "$PYTHON" "${SCRIPT_DIR}/figure_5e.py" \
        "${ANALYSES}/03_taxonomic_representation/results/growth_data.csv" \
        --set-info "${ANALYSES}/03_taxonomic_representation/results/dataset_metrics.csv" \
        --model-info "${DATA}/model_info.csv" \
        --out "${FIGURES}/figure_5e.svg" \
        --cache "$CACHE"
fi

if [[ "$1" == "5f" || "$1" == "all" ]]; then
    echo "Generating figure 5f..."
    "$PYTHON" "${SCRIPT_DIR}/figure_5f.py" \
        "${ANALYSES}/03_taxonomic_representation/results/growth_data.csv" \
        --set-info "${ANALYSES}/03_taxonomic_representation/results/dataset_metrics.csv" \
        --model-info "${DATA}/model_info.csv" \
        --out "${FIGURES}/figure_5f.svg" \
        --cache "$CACHE"
fi
