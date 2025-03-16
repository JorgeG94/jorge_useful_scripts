#!/usr/bin/env python3

import os
import sys
import argparse

def print_help():
    print("Usage: python3 make_slurms.py --path <path_to_jsons> --max-nodes-per-slurm <max_nodes> --nnodes-per-calc <nodes_per_calc> --ntasks-per-node <tasks_per_node> --wall-time <wall_time>")
    print("Generates batched Slurm scripts based on JSON files in the specified directory.")
    print("\nArguments:")
    print("  --path                Path to the directory containing JSON files.")
    print("  --max-nodes-per-slurm Maximum number of nodes to use per Slurm script.")
    print("  --nnodes-per-calc     Number of nodes required per calculation.")
    print("  --ntasks-per-node     Number of tasks per node.")
    print("  --wall-time           Wall time for each Slurm job in HH:MM:SS format.")
    print("  --help                Show this help message and exit.")
    sys.exit(0)

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Generate Slurm scripts for JSON calculations.")
parser.add_argument("--path", required=True, help="Path to the directory containing JSON files.")
parser.add_argument("--max-nodes-per-slurm", type=int, required=True, help="Maximum number of nodes per Slurm script.")
parser.add_argument("--nnodes-per-calc", type=int, required=True, help="Number of nodes required per calculation.")
parser.add_argument("--ntasks-per-node", type=int, required=True, help="Number of tasks per node.")
parser.add_argument("--wall-time", required=False, default="00:30:00",help="Wall time for the job (HH:MM:SS format), default: 00:30:00)..")

args = parser.parse_args()

json_dir = args.path
max_nodes_per_batch = args.max_nodes_per_slurm
nnodes_per_calc = args.nnodes_per_calc
tasks_per_node = args.ntasks_per_node
wall_time = args.wall_time
gpus_per_node = 8  # Fixed for your system

# Check if the directory exists
if not os.path.isdir(json_dir):
    print(f"Error: Directory '{json_dir}' does not exist.")
    sys.exit(1)

# Count the number of JSON files
json_files = [f for f in os.listdir(json_dir) if f.endswith(".json")]
num_json_files = len(json_files)

if num_json_files == 0:
    print(f"Error: No JSON files found in '{json_dir}'.")
    sys.exit(1)

# Maximum calculations per batch
max_calcs_per_batch = max_nodes_per_batch // nnodes_per_calc

# Split JSON files into batches
batches = [json_files[i:i+max_calcs_per_batch] for i in range(0, num_json_files, max_calcs_per_batch)]

# Generate Slurm scripts for each batch
for batch_idx, batch in enumerate(batches):
    slurm_filename = f"run_exess_batch_{batch_idx+1}.slurm"
    num_calcs = len(batch)
    total_nodes = num_calcs * nnodes_per_calc  # Total nodes required for this batch

    with open(slurm_filename, "w") as slurm_file:
        slurm_file.write(f"""#!/bin/bash

#SBATCH -A CHM213
#SBATCH -J exess_batch_{batch_idx+1}
#SBATCH -o %x-%j.out
#SBATCH -t {wall_time}
#SBATCH -p batch
#SBATCH -N {total_nodes}

module load cpe/24.07
module load cmake/3.27.9
module load gcc-native/13.2
module load rocm/6.3.1
module load cmake/3.27.9
module load cray-hdf5/1.14.3.1
module load craype-accel-amd-gfx90a
module load cray-libsci/23.12.5
export MAGMA_ROOT=/lustre/orion/proj-shared/chm213/software/magma-2.8.0-wrocm-6.3.1
export LD_LIBRARY_PATH=${{CRAY_LD_LIBRARY_PATH}}:${{LD_LIBRARY_PATH}}
export MPI_ROOT=$MPICH_DIR
export MPICH_GPU_SUPPORT_ENABLED=1

""")
        
        for json_file in batch:
            slurm_file.write(f'srun -N {nnodes_per_calc} --ntasks={tasks_per_node*nnodes_per_calc} --ntasks-per-node={tasks_per_node} --gpus-per-node={gpus_per_node} ./exess "{json_dir}/{json_file}" &\n')
            slurm_file.write("sleep 1\n")

        slurm_file.write("\nwait\n")

    print(f"Slurm script '{slurm_filename}' generated successfully with {num_calcs} calculations using {total_nodes} nodes.")

print("All batch Slurm scripts generated successfully.")

