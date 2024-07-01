import os
import shutil

def copy_and_rename_files(input_dir):
    for filename in os.listdir(input_dir):
        if filename.endswith(".inp"):
            input_path = os.path.join(input_dir, filename)

            # Construct the output path by adding _fmo3 before the extension
            output_filename = os.path.splitext(filename)[0] + '_fmo3' + '.inp'
            output_path = os.path.join(input_dir, output_filename)

            # Copy the file with the new name
            shutil.copy(input_path, output_path)

if __name__ == "__main__":
    # Specify the directory containing the *.inp files
    input_directory = '.'

    # Call the function to copy and rename files
    copy_and_rename_files(input_directory)

