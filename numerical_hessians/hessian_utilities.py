import os
import sys
import json
import h5py
import argparse
import numpy as np
import subprocess

import shutil

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
def calculate_inertia_tensor(atoms, positions):
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
    print(positions)
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
        mass_sqrt = np.sqrt(atomic_masses.get(atom, 0)) 
        # Calculate dot products: Px, Py, Pz
        Px = np.dot(positions[i], eigenvectors[:, 0])  # Dot product with the first eigenvector (X1)
        Py = np.dot(positions[i], eigenvectors[:, 1])  # Dot product with the second eigenvector (X2)
        Pz = np.dot(positions[i], eigenvectors[:, 2])  # Dot product with the third eigenvector (X3)

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

# Function to normalize vectors using the reciprocal square root of the scalar product
def normalize_vector(vector):
    norm_squared = np.dot(vector, vector)  # Scalar product of vector with itself
    if norm_squared > 0:  # Avoid division by zero
        norm_factor = 1 / np.sqrt(norm_squared)
        return vector * norm_factor
    else:
        return vector  # Return the vector unchanged if the norm is zero

# Normalize a list of D vectors
def normalize_D_vectors(D_vectors):
    normalized_D_vectors = [normalize_vector(D) for D in D_vectors]
    return normalized_D_vectors

def remove_spurious_vectors(D_vectors, threshold=1e-6):
    valid_D_vectors = []
    for D in D_vectors:
        dot_product = np.dot(D, D)  # Dot product of vector with itself
        if dot_product > threshold:  # Keep the vector if it has a non-zero norm
            valid_D_vectors.append(D)
    return valid_D_vectors

# Generate arbitrary displacement vectors for all 3N degrees of freedom
def generate_displacement_vectors(num_atoms):
    num_dof = 3 * num_atoms  # Total degrees of freedom (3 per atom)
    displacement_vectors = []
    for i in range(num_dof):
        vector = np.random.rand(num_dof)  # Create a random vector of length 3N
        displacement_vectors.append(vector)
    return displacement_vectors

# Apply Gram-Schmidt to get 3N - X orthogonalized vectors
def gram_schmidt_orthogonalization(vectors, reference_vectors):
    orthogonal_vectors = []
    
    for v in vectors:
        # Start with the current vector
        orthogonalized_v = np.copy(v)
        
        # Subtract projections onto reference vectors (translational and rotational)
        for ref in reference_vectors:
            projection = np.dot(orthogonalized_v, ref) / np.dot(ref, ref) * ref
            orthogonalized_v -= projection
        
        # Subtract projections onto previously orthogonalized vectors
        for u in orthogonal_vectors:
            projection = np.dot(orthogonalized_v, u) / np.dot(u, u) * u
            orthogonalized_v -= projection
        
        # Add the orthogonalized vector to the list if it is non-zero
        if np.linalg.norm(orthogonalized_v) > 1e-6:  # Avoid near-zero vectors
            orthogonal_vectors.append(orthogonalized_v)
    
    return orthogonal_vectors