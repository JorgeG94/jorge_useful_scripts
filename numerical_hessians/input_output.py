import os
import sys
import json
import h5py
import argparse
import numpy as np
import subprocess

import shutil

def cleanup_directory(directory):
    # Check if the directory exists
    if os.path.exists(directory):
        # Remove all files and directories within the directory
        for filename in os.listdir(directory):
            file_path = os.path.join(directory, filename)
            try:
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.unlink(file_path)  # Remove file or symbolic link
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)  # Remove directory
            except Exception as e:
                print(f"Failed to delete {file_path}. Reason: {e}")
    else:
        # If the directory does not exist, create it
        os.makedirs(directory)

# Function to pretty-print the Hessian matrix
def print_pretty_hessian(hessian):
    np.set_printoptions(precision=6, suppress=True)  # Set precision and suppress scientific notation
    for row in hessian:
        print("  ".join(f"{val:10.6f}" for val in row))


# Function to read an XYZ file
def read_xyz(xyz_file):
    with open(xyz_file, 'r') as f:
        lines = f.readlines()
    atom_count = int(lines[0].strip())
    atoms = []
    coordinates = []
    for i in range(2, 2 + atom_count):
        parts = lines[i].split()
        atoms.append(parts[0])
        coordinates.append([float(x) for x in parts[1:]])
    return atoms, np.array(coordinates)

# Function to write an XYZ file
def write_xyz(atoms, coordinates, filename):
    with open(filename, 'w') as f:
        f.write(f"{len(atoms)}\n\n")
        for atom, coord in zip(atoms, coordinates):
            f.write(f"{atom}  {coord[0]:.10f}  {coord[1]:.10f}  {coord[2]:.10f}\n")

# Function to generate the JSON file with the correct paths
# Function to generate the JSON file with the correct paths, now including unperturbed state
def generate_hessian_json(unperturbed_file, perturbed_geometries, output_file="hessian_input.json"):
    json_data = {
        "topologies": [{"xyz": os.path.relpath(unperturbed_file)}] + [{"xyz": os.path.relpath(filename)} for _, _, filename in perturbed_geometries],
        "model": {
            "method": "RestrictedHF",
            "basis": "STO-3G",
            "aux_basis": "cc-pVTZ-RIFIT",
            "standard_orientation": "None"
        },
        "keywords": {
            "scf": {
                "max_iters": 100,
                "fock_build_type": "HGP",
                "max_diis_history_length": 8,
                "convergence_threshold": 1e-10,
                "density_threshold": 1e-12,
                "convergence_metric": "DIIS"
            },
            "log": {"console": {"level": "Performance"}},
            "export": {
                "export_gradient": True
            }
        },
        "system": {
            "max_gpu_memory_mb": 16000
        },
        "driver": "Gradient"
    }

    # Write the JSON data to file
    with open(output_file, 'w') as f:
        json.dump(json_data, f, indent=2)
    print(f"JSON file '{output_file}' created successfully.")


# Function to collect gradients from the HDF5 file
def collect_gradients(hdf5_file):
    gradients = {}
    with h5py.File(hdf5_file, 'r') as f:
        # Iterate through all topology groups
        for topology in f.keys():
            group = f[topology]
            if "gradient" in group:
                # Read the gradient dataset
                gradient_data = group["gradient"][:]
                gradients[topology] = gradient_data

    return gradients