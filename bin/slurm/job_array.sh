#!/bin/bash

#SBATCH --partition=standard


# Get the appropriate command from the supplied file
CMD_FILE="$1"
CMD=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$CMD_FILE")

# Safety check: make sure the line has a command
if [ -z "$CMD" ]; then
    echo "No command found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

echo "Running task ${SLURM_ARRAY_TASK_ID}: $CMD"

# Run the command
bash -c "$CMD"
