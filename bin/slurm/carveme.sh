#!/bin/bash

#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G

# -------------------------------------------------------------------------- #
# NOTE: when running this script, the path to the python interpreter needs   #
# to be set using the sbatch argument                                        #
# --export=PYTHON="/path/to/venv/bin/python"                                 #
# -------------------------------------------------------------------------- #

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

# Run CarveMe
echo "Carving model and saving to ${out}..."
echo "${PYTHON} -m bates_reconstruction_data.carveme.cli.carve ${args[@]}"
"$PYTHON" -m bates_reconstruction_data.carveme.cli.carve "${args[@]}"

echo "Fixing model to remove warnings when loading with cobra..."
sed -i 's/boundaryCondition="true"/boundaryCondition="false"/g' "$out"
echo "Done."
