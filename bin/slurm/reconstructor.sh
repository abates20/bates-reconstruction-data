#!/bin/bash

#SBATCH --partition=standard
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

# -------------------------------------------------------------------------- #
# NOTE: when running this script, the path to the script to activate the     #
# virtual environment needs to be set using the sbatch argument              #
# --export=ENV="/path/to/venv/bin/activate"                                  #
# -------------------------------------------------------------------------- #


module purge
module load diamond/2.1.6 gurobi/13.0.0

# Activate virtual environment and run reconstructor
if source "$ENV"; then
    python -m bates_reconstruction_data.reconstructor "$@"
else
    echo "Failed to run reconstructor due to error activating venv."
    exit 1
fi
