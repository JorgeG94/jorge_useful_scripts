import os
import sys
import json
import h5py
import argparse
import numpy as np
import subprocess
from input_output import collect_gradients
import shutil

# Function to call the shell script after JSON is generated
def run_shell_command(json_file):
    command = ['./run.sh', json_file, '2', '1']
    try:
        result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(f"Shell command executed successfully:\n{result.stdout.decode()}")
    except subprocess.CalledProcessError as e:
        print(f"Error executing shell command:\n{e.stderr.decode()}")


def correlate_gradients_with_xyz(hdf5_file, json_file):
    # Load the JSON to get the xyz order
    with open(json_file, 'r') as f:
        data = json.load(f)

    # Get the list of xyz files
    xyz_files = [topology["xyz"] for topology in data["topologies"]]

    # Collect the gradients from the HDF5 file
    gradients = collect_gradients(hdf5_file)

    # Correlate the gradients with the XYZ files
    print("\nCorrelating gradients with XYZ files:\n")
    gradient_pairs = []  # To store positive and negative gradient pairs
    unperturbed_gradient = None  # For storing the unperturbed gradient

    # Process the unperturbed state (assumed to be the first entry)
    unperturbed_file = xyz_files[0]
    unperturbed_topology = "topology_0"
    
    if unperturbed_topology in gradients:
        unperturbed_gradient = gradients[unperturbed_topology]
        print(f"Unperturbed state gradient from {unperturbed_file}")

    # Process the perturbed geometries (positive and negative pairs)
    for i in range(1, len(xyz_files), 2):  # Start at 1, skip the unperturbed state
        pos_file = xyz_files[i]
        neg_file = xyz_files[i+1]
        topology_pos = f"topology_{i}"
        topology_neg = f"topology_{i+1}"

        if topology_pos in gradients and topology_neg in gradients:
            #print(f"Correlating gradients for {pos_file} and {neg_file}")
            gradient_pairs.append((gradients[topology_pos], gradients[topology_neg]))
    return unperturbed_gradient, gradient_pairs
