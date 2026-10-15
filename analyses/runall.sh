#!/bin/bash

SCRIPT_DIR=$(dirname "$(realpath "${BASH_SOURCE[0]}")")

bash "${SCRIPT_DIR}/01_reconstructor_vs_carveme/run.sh"
bash "${SCRIPT_DIR}/02_database_comparison/run.sh"
bash "${SCRIPT_DIR}/03_misannotation/run.sh"
bash "${SCRIPT_DIR}/04_single_misannotations/run.sh"
bash "${SCRIPT_DIR}/05_taxonomic_representation/run.sh"
