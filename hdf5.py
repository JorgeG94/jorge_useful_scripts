import h5py
import numpy as np




def find_atom_labels(file_path):
    with h5py.File(file_path, 'r') as hdf_file:
        # Call the recursive function to find the "atom_labels" dataset
        atom_labels_dataset = find_dataset_recursive(hdf_file, "atom_labels")
        
        if atom_labels_dataset:
            # Print information about the "atom_labels" dataset
            print("Atom Labels Dataset:")
            print(atom_labels_dataset[:])
        else:
            print("'atom_labels' dataset not found in the HDF5 file.")

def find_dataset_recursive(group, target_dataset_name):
    # Check if the current group has the target dataset
    if target_dataset_name in group:
        return group[target_dataset_name]

    # If not, recursively search in subgroups
    for key, value in group.items():
        if isinstance(value, h5py.Group):
            target_dataset = find_dataset_recursive(value, target_dataset_name)
            if target_dataset:
                return target_dataset

    # Return None if the target dataset is not found
    return None


def extract_hessian_and_positions(file_path, iteration_group_name):
    with h5py.File(file_path, 'r') as hdf_file:
        # Call the recursive function to find the iteration group
        iteration_group = find_group_recursive(hdf_file, iteration_group_name)
        
        if iteration_group:

            if "gradient" in iteration_group:
                gradient_dataset = iteration_group["gradient"]
                # Extract and print information about the "hessian" dataset
                print("\nGradient Dataset:")
                print(f"Shape: {gradient_dataset.shape}, Dtype: {gradient_dataset.dtype}")
                # Optionally, you can print or visualize the hessian data
                print(gradient_dataset[:])
            # Extract "hessian" and "positions" datasets within "iteration_11"
            if "hessian" in iteration_group and "positions" in iteration_group:
                hessian_dataset = iteration_group["hessian"]
                positions_dataset = iteration_group["positions"]

                # Extract and print information about the "hessian" dataset
                print("\nHessian Dataset:")
                print(f"Shape: {hessian_dataset.shape}, Dtype: {hessian_dataset.dtype}")
                # Optionally, you can print or visualize the hessian data
                hessian_array = np.array(hessian_dataset)

                # Print the Hessian matrix as an array
                print("Hessian Matrix as NumPy Array:")
                np.set_printoptions(formatter={'float': '{: 0.6f}'.format}, linewidth=150)

                print(hessian_array)
                print(hessian_dataset[:])
                hessian_transposed = np.transpose(hessian_dataset)

                # Diagonalize the Hessian matrix
                eigenvalues, eigenvectors = np.linalg.eigh(hessian_dataset)

                # Print the eigenvalues and optionally eigenvectors
                print("\nEigenvalues:")
                print(219474.6305 * np.sqrt(eigenvalues))

                # Uncomment the following lines to print eigenvectors
                print("\nEigenvectors:")
                print(eigenvectors)

                # Extract and print information about the "positions" dataset
                print("\nPositions Dataset:")
                print(f"Shape: {positions_dataset.shape}, Dtype: {positions_dataset.dtype}")
                # Optionally, you can print or visualize the positions data
                print(positions_dataset[:])
            else:
                print(f"'hessian' or 'positions' dataset not found within '{iteration_group_name}'.")
        else:
            print(f"Group '{iteration_group_name}' not found in the HDF5 file.")

def find_group_recursive(group, target_group_name):
    # Check if the current group is the target group
    if target_group_name in group:
        return group[target_group_name]

    # If not, recursively search in subgroups
    for key, value in group.items():
        if isinstance(value, h5py.Group):
            target_group = find_group_recursive(value, target_group_name)
            if target_group:
                return target_group

    # Return None if the target group is not found
    return None

if __name__ == "__main__":
    hdf5_file_path = "hdf5_logs/validation_inputs_geo_opt_qn_w1_hessian.hdf5"
    iteration_group_name = "iteration_9"
    extract_hessian_and_positions(hdf5_file_path, iteration_group_name)

    find_atom_labels(hdf5_file_path)
