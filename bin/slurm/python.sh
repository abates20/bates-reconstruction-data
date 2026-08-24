#!/bin/bash

#SBATCH --partition=standard


# -------------------------------------------------------------------------- #
# NOTE: when running this script, the path to the script to activate the     #
# virtual environment needs to be set using the sbatch argument              #
# --export=ENV="/path/to/venv/bin/activate"                                  #
# -------------------------------------------------------------------------- #


# Activate virtual environment and run python command
if source "$ENV"; then
    echo "python ${@}"
    python "$@"

else
    echo "Failed to run python command due to error activating venv."
    exit 1
fi