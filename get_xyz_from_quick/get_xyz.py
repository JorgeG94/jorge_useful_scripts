import os

def extract_coordinates(file_path):
    coordinates = []
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for line in lines[1:]:
            if line.strip():  # Ignore empty lines
                parts = line.split()
                element = parts[0]
                xyz = [float(coord) for coord in parts[1:4]]
                coordinates.append((element, xyz))
    return coordinates

def write_xyz_file(output_path, coordinates):
    with open(output_path, 'w') as output_file:
        output_file.write(str(len(coordinates)) + '\n\n')
        for element, xyz in coordinates:
            output_file.write(f'{element} {xyz[0]:.6f} {xyz[1]:.6f} {xyz[2]:.6f}\n')

def main():
    input_directory = '/home/jorgegv/work/scripts_jorge/get_xyz_from_quick'
    output_directory = '/home/jorgegv/work/scripts_jorge/get_xyz_from_quick/xyz_files'

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    for filename in os.listdir(input_directory):
        if filename.endswith('.in'):
            file_path = os.path.join(input_directory, filename)
            output_path = os.path.join(output_directory, filename.replace('.in', '.xyz'))

            coordinates = extract_coordinates(file_path)
            write_xyz_file(output_path, coordinates)

if __name__ == "__main__":
    main()

