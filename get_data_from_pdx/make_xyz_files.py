import os

def process_file(input_file_path, output_directory):
    with open(input_file_path, 'r', encoding='utf-8') as file:
        # Read the entire file into a string and split into blocks based on double newlines
        blocks = file.read().split('\n\n')

    for block in blocks:
        # Split block into lines and find the identifier
        lines = block.split('\n')
        identifier = lines[0]

        # Find the start of the coordinates by searching for the line "Coordinates"
        try:
            start_index = lines.index('Coordinates') + 1
        except ValueError:
            continue  # Skip the block if "Coordinates" not found
        
        # Extract coordinates
        coordinates = lines[start_index:]
        
        # Prepare content for the XYZ file
        xyz_content = f"{len(coordinates)}\n\n" + "\n".join(coordinates)
        
        # Create filename and path
        filename = f"{identifier}.xyz"
        file_path = os.path.join(output_directory, filename)
        
        # Write to the XYZ file
        with open(file_path, 'w', encoding='utf-8') as xyz_file:
            xyz_file.write(xyz_content)

def create_xyz_files(reactant_file_path, ts_file_path, base_directory):
    reactants_dir = os.path.join(base_directory, 'reactants')
    ts_dir = os.path.join(base_directory, 'transition_states')
    # Create directories if they don't exist
    os.makedirs(reactants_dir, exist_ok=True)
    os.makedirs(ts_dir, exist_ok=True)

    # Process each file
    process_file(reactant_file_path, reactants_dir)
    process_file(ts_file_path, ts_dir)

if __name__ == "__main__":
    base_directory = 'structures'  # Base directory for reactants and TS
    reactant_file_path = 'reactant_structures.txt'  # Adjust to your file path
    ts_file_path = 'TS_structures.txt'  # Adjust to your file path
    create_xyz_files(reactant_file_path, ts_file_path, base_directory)
    print("XYZ files created in", base_directory)

