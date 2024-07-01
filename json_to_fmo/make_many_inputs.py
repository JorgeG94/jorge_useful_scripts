import os
import subprocess

def process_files(root_dir):
    for subdir, dirs, files in os.walk(root_dir):
        
        directory_to_skip = "tests"
        if directory_to_skip in dirs:
            dirs.remove(directory_to_skip)

        for file in files:
            # Check if the file has the extension *.qdxf
            if file.endswith('.qdxf'):
                input_path = os.path.join(subdir, file)
                output_path = os.path.splitext(input_path)[0] + '_fmo2_gamess.inp'

                # Construct the command with --from qdxf
                print("converting", input_path)
                command = [
                    'python3',
                    'convert_fragmented_system.py',
                    '--input', input_path,
                    '--output', output_path,
                    '--from', 'qdxf',
                    '--to', 'gamess',
                    '--frag_level', '2'
                ]

                # Execute the command
                subprocess.run(command)
            elif file.endswith('.json'):
                input_path = os.path.join(subdir, file)

                level = [1,2]
                print("converting", input_path)
                output_path = os.path.join(subdir, file.replace('.json', '_fmo2_gamess.inp'))

                # Construct the command with --from json
                command = [
                    'python3',
                    'convert_fragmented_system.py',
                    '--input', input_path,
                    '--output', output_path,
                    '--from', 'json',
                    '--to', 'gamess',
                    '--frag_level', '2'
                ]

                # Execute the command
                subprocess.run(command)

if __name__ == "__main__":
    # Specify the root directory to start the search
    root_directory = '.'

    # Call the function to process files
    process_files(root_directory)
