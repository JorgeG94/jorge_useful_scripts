def add_separators_to_file(input_file_path, output_file_path, separator='\n'):
    with open(input_file_path, 'r', encoding='utf-8') as file:
        lines = file.readlines()
    
    with open(output_file_path, 'w', encoding='utf-8') as file:
        for i, line in enumerate(lines):
            # Check if the line looks like the start of a new block
            if '_' in line and i != 0:  # Avoid adding separator before the first block
                file.write(separator)
            file.write(line)

if __name__ == "__main__":
    input_file_path = 'cleaned_data.txt'  # Adjust to your file
    output_file_path = 'separated_cleaned_data.txt'  # Adjust to your preferred output file name
    add_separators_to_file(input_file_path, output_file_path)
    print("Data separated and saved to", output_file_path)

