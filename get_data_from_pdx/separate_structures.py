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

if __name__ == "__main__":
    input_file_path = 'separated_cleaned_data.txt'  # Adjust to your file
    ts_file_path = 'TS_structures.txt'  # File to write TS data
    reactant_file_path = 'reactant_structures.txt'  # File to write reactant data
    classify_and_write_data(input_file_path, ts_file_path, reactant_file_path)
    print("Data classified and saved to", ts_file_path, "and", reactant_file_path)

