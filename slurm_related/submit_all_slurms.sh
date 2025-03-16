#!/bin/bash

# Iterate through all .slurm files in the current directory
for slurm_file in *.slurm; do
    if [[ -f "$slurm_file" ]]; then
        echo "Submitting $slurm_file..."
        sbatch "$slurm_file"
        sleep 1 
    fi
done

echo "All Slurm jobs have been submitted."

