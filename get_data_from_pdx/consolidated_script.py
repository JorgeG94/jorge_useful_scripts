import os
import shutil
import pdfplumber

# Assuming functions for steps 2 to 5 are defined here (clean_text_file, add_separators_to_file, etc.)

def clean_text_file(input_file_path, output_file_path):
    with open(input_file_path, 'r', encoding='utf-8') as file:
        lines = file.readlines()
    
    cleaned_lines = [line for line in lines if ".md" not in line]

    with open(output_file_path, 'w', encoding='utf-8') as file:
        file.writelines(cleaned_lines)

def add_separators_to_file(input_file_path, output_file_path, separator='\n'):
    with open(input_file_path, 'r', encoding='utf-8') as file:
        lines = file.readlines()
    
    with open(output_file_path, 'w', encoding='utf-8') as file:
        for i, line in enumerate(lines):
            # Check if the line looks like the start of a new block
            if '_' in line and i != 0:  # Avoid adding separator before the first block
                file.write(separator)
            file.write(line)

def classify_and_write_data(input_file_path, ts_file_path, reactant_file_path, separator='\n'):
    with open(input_file_path, 'r', encoding='utf-8') as file:
        # Read the entire file into a single string
        data = file.read()

    # Split the data into blocks using the separator
    blocks = data.split(separator * 2)  # Assuming double newline as block separator
    
    with open(ts_file_path, 'w', encoding='utf-8') as ts_file, open(reactant_file_path, 'w', encoding='utf-8') as reactant_file:
        for block in blocks:
            if 'TS' in block.split('\n', 1)[0]:  # Checking if 'TS' is in the identifier
                ts_file.write(block + separator * 2)
            else:
                reactant_file.write(block + separator * 2)

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

def extract_text_from_pdf(pdf_path):
    text = ''
    with pdfplumber.open(pdf_path) as pdf:
        for page in pdf.pages:
            text += page.extract_text() + '\n'
    return text

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
def save_extracted_text_to_file(pdf_path, text_file_path):
    # Extract text from PDF
    text = ''
    with pdfplumber.open(pdf_path) as pdf:
        for page in pdf.pages:
            text += page.extract_text() + '\n'
    
    # Save extracted text to file
    with open(text_file_path, 'w', encoding='utf-8') as file:
        file.write(text)

def process_pdf_to_xyz(pdf_path, base_directory):
    # Step 1: Extract text from PDF
    extracted_text_file_path = os.path.join(base_directory, 'extracted_text.txt')
    
    # Step 1: Extract text from PDF and save it to a file
    save_extracted_text_to_file(pdf_path, extracted_text_file_path) 
    # Define paths for intermediate files
    cleaned_file_path = os.path.join(base_directory, 'cleaned_file.txt')
    separated_data_path = os.path.join(base_directory, 'separated_data.txt')
    reactant_file_path = os.path.join(base_directory, 'reactant_structures.txt')
    ts_file_path = os.path.join(base_directory, 'TS_structures.txt')
    
    # Step 2: Clean the extracted text
    # Assume clean_text_file function is defined
    clean_text_file(extracted_text_file_path, cleaned_file_path)
    
    # Step 3: Separate the cleaned text into blocks
    # Assume add_separators_to_file function is defined
    add_separators_to_file(cleaned_file_path, separated_data_path)
    
    # Step 4: Classify blocks into reactants and transition states
    # Assume classify_and_write_data function is defined
    classify_and_write_data(separated_data_path, reactant_file_path, ts_file_path)
    
    # Step 5: Generate XYZ files
    # Assume create_xyz_files function is defined
    create_xyz_files(reactant_file_path, ts_file_path, base_directory)
    
    # Step 6: Separate reactants without sulfur
    reactants_dir = os.path.join(base_directory, 'reactants')
    no_sulphur_dir = os.path.join(reactants_dir, 'no_sulphur')
    os.makedirs(no_sulphur_dir, exist_ok=True)
    for filename in os.listdir(reactants_dir):
        if filename.endswith(".xyz"):
            file_path = os.path.join(reactants_dir, filename)
            with open(file_path, 'r') as file:
                if not any(' S ' in line for line in file):
                    shutil.copy(file_path, os.path.join(no_sulphur_dir, filename))

if __name__ == "__main__":
    pdf_path = 'si.pdf'
    base_directory = '.'
    process_pdf_to_xyz(pdf_path, base_directory)

