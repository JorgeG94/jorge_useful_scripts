import os
import sys
import json
import h5py
import argparse
from scipy.constants import c, pi
import numpy as np
import subprocess
import shutil
from input_output import read_xyz, print_pretty_hessian, print_pretty, cleanup_directory, write_xyz, generate_hessian_json
from hessian_utilities import center_of_mass, build_hessian, compute_vibrational_frequencies,shift_to_center_of_mass,\
calculate_inertia_tensor, calculate_eigenvalues_eigenvectors, generate_translational_rotational_D_vectors,\
normalize_D_vectors, remove_spurious_vectors, orthogonalise_matrix, projection_matrix, \
calculate_wavenumbers, generate_small_displacement_vectors, generate_finite_difference_geometries, \
compute_vibrational_frequencies_internal
from utilities import correlate_gradients_with_xyz, run_shell_command
bohr_radius = 0.52917721092
atomic_masses = {
    "H": 1.00784,
    "O": 15.99977,
    "C": 12.0096
    # Add masses for other elements as needed
}


# Main function to prepare geometries for numerical Hessian
def prepare_hessian_geometries(xyz_file, delta, output_dir="num_hess", redo=True, debug=True):
    print(xyz_file)
    atoms, coordinates = read_xyz(xyz_file)
    num_atoms = len(atoms)
    coordinates_in_bohr = coordinates / bohr_radius
    print(atoms)
    N = len(atoms)
    tmp_debug = True
    com = center_of_mass(atoms, coordinates_in_bohr)
    if debug:
        print("Center of Mass:", com)
    # Extract the base name of the xyz_file (without extension)
    base_name = os.path.splitext(os.path.basename(xyz_file))[0]
    
    # Set json_file and hdf5_file dynamically based on xyz_file name
    json_file = f"num_hess_{base_name}.json"
    hdf5_file = f"hdf5_logs/num_hess_{base_name}.hdf5"

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
    #hdf5_file = "hdf5_logs/hessian_input.hdf5"
    unperturbed_gradient, gradient_pairs = correlate_gradients_with_xyz(hdf5_file, json_file)

    # Build the Hessian matrix
    
    delta = delta / bohr_radius # this turns it to bohrs
    hessian_matrix = build_hessian(unperturbed_gradient, gradient_pairs, delta)
    print("\nHessian matrix:")
    # Bohr radius in Ångströms

    hessian_corrected = hessian_matrix #* bohr_radius #* bohr_radius

    print_pretty(hessian_corrected)

    # Compute vibrational frequencies
    mass_weighted_hessian, vibrational_frequencies = compute_vibrational_frequencies(hessian_corrected, atoms)
    print("\n Frequencies (in cm^-1):\n")
    print_pretty(vibrational_frequencies)

    if debug:
        print("\n Mass weighted hessian : ")
        print_pretty(mass_weighted_hessian)
    
    # 3N - X vibrational modes 
    # X = 0, atom
    # X = 5 if molecule is linear 
    # X = 6 for anything else 

    # Molecules -> rotation, translations, vibrations 
    # Vibrations 3N-X
    # rotations and translations 


    shifted_positions = shift_to_center_of_mass(atoms, coordinates_in_bohr, com)
    if debug:
        print("Shifted Positions:", shifted_positions)

    inertia_tensor = calculate_inertia_tensor(atoms, shifted_positions)
    if debug:
        print("Inertia Tensor:\n", inertia_tensor)

    # Get eigenvalues and eigenvectors of the inertia tensor
    eigenvalues, eigenvectors = calculate_eigenvalues_eigenvectors(inertia_tensor)
    if debug:
        print("Eigenvalues (Principal Moments of Inertia):", eigenvalues)
        print("Eigenvectors (Principal Axes of Rotation):\n", eigenvectors)

    D1, D2, D3, D4, D5, D6 = generate_translational_rotational_D_vectors(atoms, shifted_positions)
    
    if debug:
        print(" \n D vectors for transformations: ")
        print("D1:", D1)
        print("D2:", D2)
        print("D3:", D3)
        print("D4:", D4)
        print("D5:", D5)
        print("D6:", D6)

    D_vectors = [D1, D2, D3, D4, D5, D6]

    trans_rot_space = np.vstack(D_vectors)
    trans_rot_space_orth = orthogonalise_matrix(trans_rot_space.T)
    projection = projection_matrix(trans_rot_space_orth, N)


    mwhess_proj = np.dot(projection.T, mass_weighted_hessian).dot(projection)

    internal_hessian = mwhess_proj

    vibrational_frequencies_internal = compute_vibrational_frequencies_internal(internal_hessian)
    if debug:
        print("\n Frequencies (in cm^-1):\n")
        print_pretty(vibrational_frequencies_internal)

# Entry point for the script
if __name__ == "__main__":
    # Use argparse for parsing command line arguments
    parser = argparse.ArgumentParser(description="Prepare geometries for numerical Hessian.")
    parser.add_argument('--input', required=True, help='Path to the input XYZ file.')
    parser.add_argument('--redo', dest='redo', action='store_true', help='Redo geometry generation and calculation.')
    parser.add_argument('--no-redo', dest='redo', action='store_false', help='Skip geometry generation and use existing data.')
    parser.add_argument('--debug-print', dest='debug', action='store_true', help='Enable verbose printing' )

    parser.set_defaults(redo=True)  # Default is to redo the steps unless --no-redo is provided
    parser.set_defaults(debug=False)

    args = parser.parse_args()

    # Call the function with the parsed arguments
    # this is angstroms 
    delta = 0.005 * bohr_radius
    prepare_hessian_geometries(args.input, delta, redo=args.redo, debug=args.debug)
