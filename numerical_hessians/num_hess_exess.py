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


bohr_radius = 0.52917721092
atomic_masses = {
    "H": 1.00784,
    "O": 15.99977,
    "C": 12.0096
    # Add masses for other elements as needed
}

# Function to get the atomic masses repeated for each coordinate (x, y, z)
def get_mass_weight_vector(atoms):
    mass_weight_vector = []
    for atom in atoms:
        mass_weight_vector.extend([atomic_masses[atom]] * 3)  # 3 coordinates (x, y, z) per atom
    return np.array(mass_weight_vector)

# Function to compute the mass-weighted Hessian
def mass_weight_hessian(hessian, atoms):
    # Get the mass vector for each coordinate
    mass_vector = get_mass_weight_vector(atoms)

    # Inverse square root of the mass matrix (diagonal matrix)
    mass_sqrt_inv = 1.0 / np.sqrt(mass_vector)

    # Create a diagonal matrix with the inverse sqrt of masses
    mass_sqrt_inv_matrix = np.diag(mass_sqrt_inv)

    # Apply mass-weighting: M^{-1/2} * H * M^{-1/2}
    mass_weighted_hessian = np.dot(mass_sqrt_inv_matrix, np.dot(hessian, mass_sqrt_inv_matrix))

    return mass_weighted_hessian
# Function to compute vibrational frequencies from the Hessian
def compute_vibrational_frequencies(hessian, atoms):
    # Mass-weight the Hessian
    mass_weighted_hessian = mass_weight_hessian(hessian, atoms)
    #print("\nMass weighted Hessian ")
    #print_pretty_hessian(mass_weighted_hessian)

    #print("\nMass weighted Hessian GAMESS ")
    #print_pretty_hessian(gms_mass_weight)

    #print("\n Diff between hess")
    #print_pretty_hessian(diff)

    # Diagonalize the mass-weighted Hessian to get eigenvalues
    np.set_printoptions(precision=12)
    #eigenvalues, _ = np.linalg.eigh(gms_mass_weight)
    eigenvalues, _ = np.linalg.eig(mass_weighted_hessian)

    # Remove small or negative eigenvalues (these correspond to translations/rotations)
    positive_eigenvalues = eigenvalues[eigenvalues > 0]

    # Convert eigenvalues to vibrational frequencies in cm^-1
    # 1 atomic unit of frequency = 2.1947 * 10^5 cm^-1
    conversion_factor = 2.642461e7
    frequencies_cm1 = np.sqrt(abs(eigenvalues) * conversion_factor)

    return frequencies_cm1

# Function to pretty-print the Hessian matrix
def print_pretty_hessian(hessian):
    np.set_printoptions(precision=6, suppress=True)  # Set precision and suppress scientific notation
    for row in hessian:
        print("  ".join(f"{val:10.6f}" for val in row))

def center_of_mass(atoms, positions):
    total_mass = 0
    weighted_position_sum = [0.0, 0.0, 0.0]

    for atom, position in zip(atoms, positions):
        mass = atomic_masses.get(atom, 0)  # Look up atomic mass, 0 if not found
        total_mass += mass
        for i in range(3):
            weighted_position_sum[i] += mass * position[i]

    if total_mass == 0:
        raise ValueError("Total mass cannot be zero!")

    # Calculate the center of mass
    center_of_mass = [coord / total_mass for coord in weighted_position_sum]
    return center_of_mass

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

# Function to generate perturbed geometries for finite differences
def generate_finite_difference_geometries(atoms, coordinates, delta, output_dir="num_hess"):
    perturbed_geometries = []
    num_atoms = len(atoms)
    num_coords = 3 * num_atoms

    # Ensure the output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    print(f"The delta is {delta}")
    # Save unperturbed state as the first geometry
    unperturbed_file = os.path.join(output_dir, "unperturbed.xyz")
    write_xyz(atoms, coordinates, unperturbed_file)

    for i in range(num_coords):
        # Generate positive displacement
        displaced_coords_pos = np.copy(coordinates)
        displaced_coords_pos[i // 3, i % 3] += delta
        filename_pos = os.path.join(output_dir, f"geom_pos_{i}.xyz")
        perturbed_geometries.append((atoms, displaced_coords_pos, filename_pos))

        # Generate negative displacement
        displaced_coords_neg = np.copy(coordinates)
        displaced_coords_neg[i // 3, i % 3] -= delta
        filename_neg = os.path.join(output_dir, f"geom_neg_{i}.xyz")
        perturbed_geometries.append((atoms, displaced_coords_neg, filename_neg))

    return unperturbed_file, perturbed_geometries

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

# Function to call the shell script after JSON is generated
def run_shell_command(json_file):
    command = ['./run.sh', json_file, '2', '1']
    try:
        result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(f"Shell command executed successfully:\n{result.stdout.decode()}")
    except subprocess.CalledProcessError as e:
        print(f"Error executing shell command:\n{e.stderr.decode()}")

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
def collect_gradients_fp64(hdf5_file):
    gradients = {}
    
    with h5py.File(hdf5_file, 'r') as f:
        # Iterate through all topology groups
        for topology in f.keys():
            group = f[topology]
            
            if "gradient" in group:
                # Read the gradient dataset
                gradient_dataset = group["gradient"]

                # Check the dtype of the gradient dataset
                if gradient_dataset.dtype == np.float64:
                    print(f"Gradient in {topology} is in double precision (float64).")

                    print(f"Warning: Gradient in {topology} is not in float64. It is stored as {gradient_dataset.dtype}.")
                    # Optionally, cast it to float64 if needed
                    gradient_data = np.array(gradient_dataset[:], dtype=np.float64)
                    print("Gradient data has been cast to float64.")
                else:
                    gradient_data = gradient_dataset[:]
                
                # Store the gradient data in the dictionary
                gradients[topology] = gradient_data

    return gradients

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

# Function to build the correct 9x9 Hessian matrix from gradient pairs
def build_hessian(grad_unp, gradient_pairs, delta):
    num_coords = len(gradient_pairs)  # Number of coordinates (half the size of gradient_pairs)
    num_atoms = int(num_coords / 3)  # Each atom has 3 coordinates (x, y, z)
    
    # Initialize a 9x9 Hessian matrix
    hessian = np.zeros((num_coords, num_coords))

    # Iterate over each displacement to fill the Hessian
    for i in range(num_coords):
        # Access the positive and negative gradient for this coordinate displacement
        grad_pos = np.array(gradient_pairs[i][0])  # Positive displacement
        grad_neg = np.array(gradient_pairs[i][1])  # Negative displacement

        # Compute the second derivative using finite differences
        second_derivative = (grad_pos - grad_neg ) / ( 2* delta)

        # Place this second derivative in the appropriate row of the Hessian
        # We want the i-th row of the Hessian, and the full row is determined by the flattened gradient
        hessian[i, :] = second_derivative.flatten()

    return hessian
def shift_to_center_of_mass(atoms, positions, center_of_mass):
    shifted_positions = []
    for position in positions:
        shifted_position = [position[i] - center_of_mass[i] for i in range(3)]
        shifted_positions.append(shifted_position)
    return shifted_positions
def calculate_inertia_tensor_programmatic(atoms, positions):
    # Convert positions to a NumPy array for vectorized operations
    positions = np.array(positions)
    masses = np.array([atomic_masses.get(atom, 0) for atom in atoms])

    # Extract x, y, z coordinates as separate arrays
    x, y, z = positions[:, 0], positions[:, 1], positions[:, 2]

    # Calculate diagonal elements (Ixx, Iyy, Izz)
    Ixx = np.sum(masses * (y**2 + z**2))
    Iyy = np.sum(masses * (x**2 + z**2))
    Izz = np.sum(masses * (x**2 + y**2))

    # Calculate off-diagonal elements (Ixy, Ixz, Iyz)
    Ixy = -np.sum(masses * x * y)
    Ixz = -np.sum(masses * x * z)
    Iyz = -np.sum(masses * y * z)

    # Construct the inertia tensor as a symmetric matrix
    inertia_tensor = np.array([[Ixx, Ixy, Ixz],
                               [Ixy, Iyy, Iyz],
                               [Ixz, Iyz, Izz]])
    
    return inertia_tensor

def calculate_eigenvalues_eigenvectors(inertia_tensor):
    # Use NumPy to calculate eigenvalues and eigenvectors
    eigenvalues, eigenvectors = np.linalg.eig(inertia_tensor)
    return eigenvalues, eigenvectors

def generate_D_vectors(atoms, positions):
    N = len(atoms)  # Number of atoms
    D1 = np.zeros(3 * N)
    D2 = np.zeros(3 * N)
    D3 = np.zeros(3 * N)

    for i, atom in enumerate(atoms):
        mass_sqrt = np.sqrt(atomic_masses.get(atom, 0))  # sqrt(mass)
        # Assign values for each axis (x, y, z)
        D1[3 * i] = mass_sqrt * positions[i][0]  # D1 corresponds to x axis
        D2[3 * i + 1] = mass_sqrt * positions[i][1]  # D2 corresponds to y axis
        D3[3 * i + 2] = mass_sqrt * positions[i][2]  # D3 corresponds to z axis

    return D1, D2, D3

def generate_rotational_D_vectors(atoms, positions, eigenvectors):
    N = len(atoms)  # Number of atoms
    D4 = np.zeros(3 * N)
    D5 = np.zeros(3 * N)
    D6 = np.zeros(3 * N)

    for i, atom in enumerate(atoms):
        mass_sqrt = np.sqrt(atomic_masses.get(atom, 0))  # sqrt(mass)
        Px, Py, Pz = positions[i]  # Atom's position with respect to center of mass

        # Extract the eigenvectors (X matrix)
        X1, X2, X3 = eigenvectors[:, 0], eigenvectors[:, 1], eigenvectors[:, 2]

        # Calculate D4, D5, and D6 for the i-th atom using the provided formulas
        D4[3 * i]     = (Py * X3[0] - Pz * X2[0]) * mass_sqrt  # x component
        D4[3 * i + 1] = (Py * X3[1] - Pz * X2[1]) * mass_sqrt  # y component
        D4[3 * i + 2] = (Py * X3[2] - Pz * X2[2]) * mass_sqrt  # z component

        D5[3 * i]     = (Pz * X1[0] - Px * X3[0]) * mass_sqrt  # x component
        D5[3 * i + 1] = (Pz * X1[1] - Px * X3[1]) * mass_sqrt  # y component
        D5[3 * i + 2] = (Pz * X1[2] - Px * X3[2]) * mass_sqrt  # z component

        D6[3 * i]     = (Px * X2[0] - Py * X1[0]) * mass_sqrt  # x component
        D6[3 * i + 1] = (Px * X2[1] - Py * X1[1]) * mass_sqrt  # y component
        D6[3 * i + 2] = (Px * X2[2] - Py * X1[2]) * mass_sqrt  # z component

    return D4, D5, D6

def normalize_vector(vector):
    norm = np.linalg.norm(vector)
    if norm != 0:
        return vector / norm
    else:
        return vector  # In case norm is zero, we return the original vector

def normalize_D_vectors(D1, D2, D3, D4, D5, D6):
    D1_normalized = normalize_vector(D1)
    D2_normalized = normalize_vector(D2)
    D3_normalized = normalize_vector(D3)
    D4_normalized = normalize_vector(D4)
    D5_normalized = normalize_vector(D5)
    D6_normalized = normalize_vector(D6)

    return D1_normalized, D2_normalized, D3_normalized, D4_normalized, D5_normalized, D6_normalized

def remove_spurious_vectors(D_vectors, threshold=1e-6):
    valid_D_vectors = []
    for D in D_vectors:
        dot_product = np.dot(D, D)  # Dot product of vector with itself
        if dot_product > threshold:  # Keep the vector if it has a non-zero norm
            valid_D_vectors.append(D)
    return valid_D_vectors

# Main function to prepare geometries for numerical Hessian
def prepare_hessian_geometries(xyz_file, delta, output_dir="num_hess", json_file="hessian_input.json", redo=True):
    atoms, coordinates = read_xyz(xyz_file)
    num_atoms = len(atoms)
    coordinates_in_bohr = coordinates / bohr_radius
    print(atoms)
    com = center_of_mass(atoms, coordinates_in_bohr)
    print("Center of Mass:", com)
    if redo:

        # Number of expected geometries
        expected_geometries = 6 * num_atoms 

        cleanup_directory(output_dir)

        # Generate perturbed geometries
        unperturbed_file, perturbed_geometries = generate_finite_difference_geometries(atoms, coordinates, delta, output_dir)

        # Write all perturbed geometries to files
        for atoms, coords, filename in perturbed_geometries:
            write_xyz(atoms, coords, filename)

        # Check if the number of generated geometries is correct
        if len(perturbed_geometries) == expected_geometries:
            print(f"Success: Generated {len(perturbed_geometries)} geometries, as expected.")
        else:
            print(f"Warning: Expected {expected_geometries} geometries, but generated {len(perturbed_geometries)}.")

        # Generate the JSON file
        generate_hessian_json(unperturbed_file, perturbed_geometries, json_file)

        # Call the shell script with the generated JSON file
        run_shell_command(json_file)

    # Correlate the gradients with XYZ files after the run
    hdf5_file = "hdf5_logs/hessian_input.hdf5"
    unperturbed_gradient, gradient_pairs = correlate_gradients_with_xyz(hdf5_file, json_file)

    # Build the Hessian matrix
    hessian_matrix = build_hessian(unperturbed_gradient, gradient_pairs, delta)
    #print("\nHessian matrix:")
    # Bohr radius in Ångströms

    hessian_corrected = hessian_matrix * bohr_radius
    #print_pretty_hessian(hessian_corrected)

    # Compute vibrational frequencies
    vibrational_frequencies = compute_vibrational_frequencies(hessian_corrected, atoms)
    print("\nVibrational Frequencies (in cm^-1):\n", vibrational_frequencies)
'''
    shifted_positions = shift_to_center_of_mass(atoms, coordinates_in_bohr, com)
    print("Shifted Positions:", shifted_positions)

    inertia_tensor = calculate_inertia_tensor_programmatic(atoms, shifted_positions)
    print("Inertia Tensor:\n", inertia_tensor)

    # Get eigenvalues and eigenvectors of the inertia tensor
    eigenvalues, eigenvectors = calculate_eigenvalues_eigenvectors(inertia_tensor)

    print("Eigenvalues (Principal Moments of Inertia):", eigenvalues)
    print("Eigenvectors (Principal Axes of Rotation):\n", eigenvectors)

    D1, D2, D3 = generate_D_vectors(atoms, shifted_positions)
    D4, D5, D6 = generate_rotational_D_vectors(atoms, shifted_positions, eigenvectors)

    print("D1:", D1)
    print("D2:", D2)
    print("D3:", D3)
    print("D4:", D4)
    print("D5:", D5)
    print("D6:", D6)

    D1_norm, D2_norm, D3_norm, D4_norm, D5_norm, D6_norm = normalize_D_vectors(D1, D2, D3, D4, D5, D6)

    print("Normalized D1:", D1_norm)
    print("Normalized D2:", D2_norm)
    print("Normalized D3:", D3_norm)
    print("Normalized D4:", D4_norm)
    print("Normalized D5:", D5_norm)
    print("Normalized D6:", D6_norm)

    D_vectors = [D1_norm, D2_norm, D3_norm, D4_norm, D5_norm, D6_norm]
    valid_D_vectors = remove_spurious_vectors(D_vectors)

    print("Valid D vectors:", valid_D_vectors, "number of vectors", len(valid_D_vectors))
'''
# Entry point for the script
if __name__ == "__main__":
    # Use argparse for parsing command line arguments
    parser = argparse.ArgumentParser(description="Prepare geometries for numerical Hessian.")
    parser.add_argument('--input', required=True, help='Path to the input XYZ file.')
    parser.add_argument('--redo', dest='redo', action='store_true', help='Redo geometry generation and calculation.')
    parser.add_argument('--no-redo', dest='redo', action='store_false', help='Skip geometry generation and use existing data.')

    parser.set_defaults(redo=True)  # Default is to redo the steps unless --no-redo is provided

    args = parser.parse_args()

    # Call the function with the parsed arguments
    delta = 0.005 * bohr_radius
    prepare_hessian_geometries(args.input, delta, redo=args.redo)
