def clean_text_file(input_file_path, output_file_path):
    with open(input_file_path, 'r', encoding='utf-8') as file:
        lines = file.readlines()
    
    cleaned_lines = [line for line in lines if ".md" not in line]

    with open(output_file_path, 'w', encoding='utf-8') as file:
        file.writelines(cleaned_lines)

if __name__ == "__main__":
    input_file_path = 'data.txt'  # The path to your input file
    output_file_path = 'cleaned_data.txt'  # The path for the cleaned output file
    clean_text_file(input_file_path, output_file_path)
    print("File cleaned and saved as", output_file_path)

