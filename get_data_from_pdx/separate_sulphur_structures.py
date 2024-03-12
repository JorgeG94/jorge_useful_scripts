import os
import shutil

def separate_sulphur_xyz(reactants_dir, no_sulphur_dir):
    # Create the no_sulphur directory if it doesn't exist
    os.makedirs(no_sulphur_dir, exist_ok=True)
    
    # Iterate over all files in the reactants directory
    for filename in os.listdir(reactants_dir):
        if filename.endswith(".xyz"):  # Ensure we're only processing .xyz files
            file_path = os.path.join(reactants_dir, filename)
            
            # Read the contents of the XYZ file
            with open(file_path, 'r', encoding='utf-8') as file:
                contents = file.readlines()
            
            # Check if "S" is present in the file
            if not any('S ' in line for line in contents):
                # Copy the file to the no_sulphur directory if "S" is not found
                shutil.copy(file_path, os.path.join(no_sulphur_dir, filename))

if __name__ == "__main__":
    reactants_dir = 'structures/reactants'  # Source directory
    no_sulphur_dir = 'structures/reactants/no_sulphur'  # Destination directory for files without sulfur
    separate_sulphur_xyz(reactants_dir, no_sulphur_dir)
    print("Files without sulfur atoms have been separated into", no_sulphur_dir)

