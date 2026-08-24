#!/bin/bash

#SBATCH --partition=standard
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

# -------------------------------------------------------------------------- #
# NOTE: when running this script, the path to the python interpreter needs   #
# to be set using the sbatch argument                                        #
# --export=PYTHON="/path/to/venv/bin/python"                                 #
# -------------------------------------------------------------------------- #


module purge
module load diamond/2.1.6

# Run reconstructor
"$PYTHON" -m bates_reconstruction_data.reconstructor "$@"
