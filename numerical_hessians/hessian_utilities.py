import os
import sys
import json
import h5py
import argparse
import numpy as np
from input_output import *
import subprocess
from scipy.constants import c, pi
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
    # Seminumerical hessians -> gradients and thee energy are analytic 
    # atom x y z -> 3N coordinates, N = number of atoms 
    # for each(coordinate):
    # coord += delta
    # coords -= delta
    # 6N 
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

def print_pretty_for_cpp(vector):
    """
    Function to print the elements of a vector in C++ format {value, value, value}.
    Handles both 1D and 2D numpy arrays.
    """
    # Flatten the vector if it's a multi-dimensional array
    if isinstance(vector, np.ndarray) and vector.ndim > 1:
        vector = vector.flatten()

    # Iterate over the elements of the vector and format each element
    formatted_vector = ', '.join(f'{val:.16e}' for val in vector)
    print(f'{{ {formatted_vector} }},')

# Function to build the correct 9x9 Hessian matrix from gradient pairs
def build_hessian(grad_unp, gradient_pairs, delta):
    num_coords = len(gradient_pairs)  # Number of coordinates 
    num_atoms = int(num_coords / 3)  # Each atom has 3 coordinates (x, y, z)
    
    # Initialize a 9x9 Hessian matrix
    hessian = np.zeros((num_coords, num_coords))

    # Iterate over each displacement to fill the Hessian
    print(f"delta is {delta}")
    for i in range(num_coords):
        grad_pos = np.array(gradient_pairs[i][0])  # Positive displacement
        grad_neg = np.array(gradient_pairs[i][1])  # Negative displacement
    
        # Compute the second derivative using finite differences
        second_derivative = (grad_pos - grad_neg) / (2 * delta)
    
        # Flatten the second derivative
        flat_second_derivative = second_derivative.flatten()
    
        # Assign values to ensure symmetry in the Hessian
        for j in range(len(flat_second_derivative)):
            hessian[i, j] = flat_second_derivative[j]
            hessian[j, i] = flat_second_derivative[j]  # Symmetric assignment
    
    return hessian


'''
    # Iterate over each displacement to fill the Hessian
    for i in range(num_coords):
        #print(gradient_pairs[i])
        # print_pretty_for_cpp(gradient_pairs[i][0])
        # print_pretty_for_cpp(gradient_pairs[i][1])
        # Access the positive and negative gradient for this coordinate displacement
        grad_pos = np.array(gradient_pairs[i][0])  # Positive displacement
        grad_neg = np.array(gradient_pairs[i][1])  # Negative displacement

        # Compute the second derivative using finite differences
        
        second_derivative = (grad_pos - grad_neg ) / ( 2* delta)

        # Place this second derivative in the appropriate row of the Hessian
        # We want the i-th row of the Hessian, and the full row is determined by the flattened gradient
        print("\n second der = ", second_derivative.flatten())
        hessian[i, :] = second_derivative.flatten()

    return hessian
'''
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

def calculate_eigenvalues_eigenvectors(matrix):
    # Use NumPy to calculate eigenvalues and eigenvectors
    eigenvalues, eigenvectors = np.linalg.eigh(matrix)
    return eigenvalues, eigenvectors

def generate_translational_rotational_D_vectors(atoms, positions):

    m = np.array([atomic_masses[atom] for atom in atoms])
    print(f"m shape: {m.shape}")
    z = np.zeros_like(m)
    i = np.ones_like(m)
    ux = np.ravel([i, z, z], order='F')
    uy = np.ravel([z, i, z], order='F')
    uz = np.ravel([z, z, i], order='F')

    geom = np.array(positions)
    xxx = np.repeat(geom[:, 0], 3)
    yyy = np.repeat(geom[:, 1], 3)
    zzz = np.repeat(geom[:, 2], 3)

    sqrtmmm = np.repeat(np.sqrt(m), 3)
    R4 = sqrtmmm * (yyy * uz - zzz * uy)
    R5 = sqrtmmm * (zzz * ux - xxx * uz)
    R6 = sqrtmmm * (xxx * uy - yyy * ux)

    
    T1 = sqrtmmm * ux
    T2 = sqrtmmm * uy
    T3 = sqrtmmm * uz

    print(f"T1: {T1}\n")
    print(f"T2: {T2}\n")
    print(f"T3: {T3}\n")

    return T1, T2, T3, R4, R5, R6

def orthogonalise_matrix(matrix, tol=1e-6):
    u, s, _ = np.linalg.svd(matrix, full_matrices=False)
    M, N = matrix.shape
    eps = np.finfo(float).eps
    if tol is None:
        tol = max(M, N) * np.amax(s) * eps
    num = np.sum(s > tol, dtype=int)
    Q = u[:, :num]
    return Q.T

def projection_matrix(matrix, N):
    P = np.identity(3 * N)
    for irt in matrix:
        P -= np.outer(irt, irt)
    return P

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

def generate_small_displacement_vectors(num_atoms, displacement_magnitude=0.005):
    num_dof = 3 * num_atoms  # Total degrees of freedom (3 per atom)
    displacement_vectors = []
    displacement_magnitude=0.005 * bohr_radius
    for i in range(num_dof):
        # Create a zero vector for each degree of freedom
        vector = np.zeros(num_dof)
        # Set the displacement in one degree of freedom (x, y, or z) to the specified magnitude
        vector[i] = displacement_magnitude
        displacement_vectors.append(vector)
    
    return displacement_vectors

def print_matrix_curly_braces(matrix):
    # Convert each row to the desired format
    formatted_rows = ["{" + ", ".join(f"{value:.19f}" for value in row) + "}" for row in matrix]
    # Join all rows with newline characters and wrap with curly braces for the entire matrix
    formatted_matrix = "{" + ",\n ".join(formatted_rows) + "}"
    print(formatted_matrix)

def calculate_wavenumbers(eigenvalues):
    conversion_factor = 2.642461e7 
    wavenumbers = np.sqrt((np.abs(eigenvalues)*conversion_factor)) 
    wavenumbers[eigenvalues < 0] *= -1  # Multiply by -1 for negative eigenvalues
    return wavenumbers

# Function to compute vibrational frequencies from the Hessian
def compute_vibrational_frequencies(hessian, atoms):
    # Mass-weight the Hessian
    mass_weighted_hessian = mass_weight_hessian(hessian, atoms)
    print("\nMass weighted Hessian ")
    print(np.array2string(mass_weighted_hessian, separator=', ', precision=6, suppress_small=False))
    #print(mass_weighted_hessian)

    # Diagonalize the mass-weighted Hessian to get eigenvalues
    np.set_printoptions(precision=12)
    #eigenvalues, _ = np.linalg.eigh(gms_mass_weight)
    eigenvalues, _ = np.linalg.eig(mass_weighted_hessian)
    real_eigenvalues = np.real(eigenvalues)
    print("\nEigenvalues")
    print_pretty(real_eigenvalues)

    # Remove small or negative eigenvalues (these correspond to translations/rotations)
    positive_eigenvalues = eigenvalues[real_eigenvalues > 0]

    # Convert eigenvalues to vibrational frequencies in cm^-1
    # 1 atomic unit of frequency = 2.1947 * 10^5 cm^-1
    conversion_factor = 2.642461e7

    frequencies_cm1 = np.sqrt(abs(real_eigenvalues) * conversion_factor)

    return mass_weighted_hessian, frequencies_cm1

def compute_vibrational_frequencies_internal(internal_hessian):
    eigenvalues, _ = np.linalg.eig(internal_hessian)
    real_eigenvalues = np.real(eigenvalues)
    print("\nEigenvalues from Internal Hessian")
    print_pretty(real_eigenvalues)

    # Convert eigenvalues to vibrational frequencies in cm^-1
    # 1 atomic unit of frequency = 2.1947 * 10^5 cm^-1
    conversion_factor = 2.642461e7

    frequencies_cm1 = np.sqrt(abs(real_eigenvalues) * conversion_factor)

    return frequencies_cm1