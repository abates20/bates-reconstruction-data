#!/bin/bash

# Filepaths
SCRIPT_DIR=$(dirname "$(realpath "${BASH_SOURCE[0]}")")
PROJECT_DIR=$(cd "${SCRIPT_DIR}/../../" && pwd)
BIN="${PROJECT_DIR}/bin"
DATA="${PROJECT_DIR}/data"
NCBI_DMP="${SCRIPT_DIR}/ncbi_dmp"
RESULTS="${SCRIPT_DIR}/results"

# Make sure needed directories exist
mkdir -p "$NCBI_DMP"
mkdir -p "$RESULTS"

# Ignore generated files in git
echo "*" > "${NCBI_DMP}/.gitignore"
echo "*" > "${RESULTS}/.gitignore"

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


# Download NCBI taxonomic data
if [ ! -f "${NCBI_DMP}/readme.txt" ]; then # Use README as marker to see if already downloaded
    curl https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz | gunzip -c | tar xf - -C "$NCBI_DMP"
fi

# Download data from KEGG
kegg_counts="${RESULTS}/kegg_counts.csv"
kegg_functions="${RESULTS}/kegg_functions.csv"
if [[ ! -f "$kegg_counts" || ! -f "$kegg_functions" ]]; then
    echo "Downloading data from KEGG..."
    "$PYTHON" "${SCRIPT_DIR}/get_kegg_data.py" \
        --tax-nodes "${NCBI_DMP}/nodes.dmp" \
        --lineages "${NCBI_DMP}/lineages.dmp" \
        --counts-out "$kegg_counts" \
        --functions-out "$kegg_functions"
fi

# Get data from BiGG
bigg_counts="${RESULTS}/bigg_counts.csv"
bigg_functions="${RESULTS}/bigg_functions.csv"
if [[ ! -f "$bigg_data" || ! -f "$bigg_functions" ]]; then
    echo "Getting data from BiGG..."
    "$PYTHON" "${SCRIPT_DIR}/get_bigg_data.py" \
        --counts-out "$bigg_counts" \
        --functions-out "$bigg_functions"
fi
