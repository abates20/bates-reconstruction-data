#!/bin/bash

#SBATCH --partition=standard
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G

# -------------------------------------------------------------------------- #
# NOTE: when running this script, the path to the script to activate the     #
# virtual environment needs to be set using the sbatch argument              #
# --export=ENV="/path/to/venv/bin/activate"                                  #
# -------------------------------------------------------------------------- #

module purge
module load diamond/2.1.6 gurobi/13.0.0

# Pull out output path from args
args=()
while [[ $# -gt 0 ]]; do
    case $1 in
        -o|--output)
            out="$2"
            args+=("$1")
            args+=("$2")
            shift 2
            ;;
        *)
            args+=("$1")
            shift
            ;;
    esac
done

# Activate virtual environment and run CarveMe
if source "$ENV"; then
    echo "Carving model and saving to ${out}..."
    echo "python -m bates_reconstruction_data.carveme.cli.carve ${args[@]}"
    python -m bates_reconstruction_data.carveme.cli.carve "${args[@]}"
    
    echo "Fixing model to remove warnings when loading with cobra..."
    sed -i 's/boundaryCondition="true"/boundaryCondition="false"/g' "$out"
    echo "Done."
    
else
    echo "Failed to run CarveMe due to error activating venv."
    exit 1
fi
